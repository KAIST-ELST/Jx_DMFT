################## Fourier transformation of block G or V by FFTW subroutine

function FFTW_k2r_block(Klist,Rlist,atom12_list,F_k_iWn_block)
    F_R_iWn_block=[];
    for i=1:size(atom12_list)[1]
        F_R_iWn_atom12=[];
        for j=1:2
            F_R_iWn_block_atomj=zeros(ComplexF64,size(F_k_iWn_block[i][j]))
            FFTW_k2r_paralle(Klist,Rlist,F_k_iWn_block[i][j],F_R_iWn_block_atomj) # 1 & 2 atom FFT
            push!(F_R_iWn_atom12,F_R_iWn_block_atomj)
        end
        push!(F_R_iWn_block,F_R_iWn_atom12)
    end
    return F_R_iWn_block
end



function FFTW_r2k_block(Rlist,Klist,atom12_list,F_R_iWn_block)
    F_k_iWn_block=[];
    for i=1:size(atom12_list)[1]
        F_k_iWn_atom12=[];
        for j=1:2
            F_k_iWn_block_atomj=zeros(ComplexF64,size(F_R_iWn_block[i][j]))
            FFTW_r2k_paralle(Rlist,Klist,F_R_iWn_block[i][j],F_k_iWn_block_atomj) # 1 & 2 atom FFT
            push!(F_k_iWn_atom12,F_k_iWn_block_atomj)
        end
        push!(F_k_iWn_block,F_k_iWn_atom12)
    end
    return F_k_iWn_block
end



################## Fourier transform by FFTW subroutine

function FFTW_k2r_paralle_Hk(Klist,Rlist,H_k,H_R)
    for l1=1:size(H_k)[4]
        for l2=1:size(H_k)[5]
            for s=1:size(H_k)[6]
                H_R[:,:,:,l1,l2,s]=FFTW.AbstractFFTs.ifft(H_k[:,:,:,l1,l2,s])
            end
        end
    end
    return H_R
end



function FFTW_k2r_paralle(Klist,Rlist,VG_k_,VG_R_)
    for w = 1:size(VG_k_)[4]
        for l1=1:size(VG_k_)[5]
            for l2=1:size(VG_k_)[6]
                for s=1:size(VG_k_)[7]
                    VG_R_[:,:,:,w,l1,l2,s]=FFTW.AbstractFFTs.ifft(VG_k_[:,:,:,w,l1,l2,s])
                end
            end
        end
    end
    return VG_R_
end



function FFTW_r2k_paralle(Rlist,Klist,VG_R2_,VG_k2_)
    for w = 1:size(VG_R2_)[4]
        for l1=1:size(VG_R2_)[5]
            for l2=1:size(VG_R2_)[6]
                for s=1:size(VG_R2_)[7]
                    VG_k2_[:,:,:,w,l1,l2,s]=FFTW.AbstractFFTs.fft(VG_R2_[:,:,:,w,l1,l2,s])
                end
            end
        end
    end
    return VG_k2_
end



################## Fourier transformdirectly 

@everywhere function fourier_k2r_internal(args) # Kpoint paralle
    Rvect = args[1]
    Klist = args[2]
    GV_k = args[3]
    
    
    GV_R_tmp =zeros(ComplexF64,size(GV_k)[4:end])
    #H_r_tmp =  (H_k[j]*exp(-2*pi*im*sum(Rvect.*Kvect)))
    for k1 = 1:size(Klist)[1]
        for k2 = 1:size(Klist)[2]
            for k3 = 1:size(Klist)[3]
                
                
                if length(size(GV_k)[4:end]) == 4
                    GV_R_tmp[:,:,:,:] += (GV_k[k1,k2,k3,:,:,:,:]*exp(2*pi*im*sum(Rvect.*Klist[k1,k2,k3,:])));
                end
                
                if length(size(GV_k)[4:end]) == 3
                    GV_R_tmp[:,:,:] += (GV_k[k1,k2,k3,:,:,:]*exp(2*pi*im*sum(Rvect.*Klist[k1,k2,k3,:])));
                end                
                
            end
        end
    end
    GV_R_tmp /= ( size(Klist)[1]*size(Klist)[2]*size(Klist)[3] );
    return GV_R_tmp
end



function fourier_k2r_paralle(Klist,Rlist,GV_k,GV_R)
    arg_list=[]
    
    for R1 = 1:size(Rlist)[1]
        for R2 = 1:size(Rlist)[2]
            for R3 = 1:size(Rlist)[3]
                Rvect = Rlist[R1,R2,R3,:]
                push!(arg_list,( Rvect, Klist, GV_k ) )
            end
        end
    end
    GV_r_list_tmp = pmap(fourier_k2r_internal, arg_list)
    i=0
    for R1 = 1:size(Rlist)[1]
        for R2 = 1:size(Rlist)[2]
            for R3 = 1:size(Rlist)[3]
                i=i+1
                
                if length(size(GV_k)[4:end]) == 4
                    GV_R[R1,R2,R3,:,:,:,:] = GV_r_list_tmp[i];
                end
                
                if length(size(GV_k)[4:end]) == 3
                    GV_R[R1,R2,R3,:,:,:] = GV_r_list_tmp[i];
                end                                  
                
            end
        end
    end
    return GV_R
end




@everywhere function fourier_r2k_internal(args) # Kpoint paralle
    Kvect = args[1]
    Rlist = args[2]
    GV_R = args[3]
    
    GV_k_tmp =zeros(ComplexF64,size(GV_R)[4:end])
    #H_r_tmp =  (H_k[j]*exp(-2*pi*im*sum(Rvect.*Kvect)))
    for R1 = 1:size(Rlist)[1]
        for R2 = 1:size(Rlist)[2]
            for R3 = 1:size(Rlist)[3]
                if length(size(GV_R)[4:end]) == 4
                    GV_k_tmp[:,:,:,:] += (GV_R[R1,R2,R3,:,:,:,:]*exp(-2*pi*im*sum(Rlist[R1,R2,R3,:].*Kvect)));
                end
                
                if length(size(GV_R)[4:end]) == 3
                    GV_k_tmp[:,:,:] += (GV_R[R1,R2,R3,:,:,:]*exp(-2*pi*im*sum(Rlist[R1,R2,R3,:].*Kvect)));
                end
            end
        end
    end
    return GV_k_tmp
end



function fourier_r2k_paralle(Rlist,Klist,GV_R,GV_k)
    arg_list=[]
    
    for k1 = 1:size(Klist)[1]
        for k2 = 1:size(Klist)[2]
            for k3 = 1:size(Klist)[3]
                Kvect = Klist[k1,k2,k3,:]
                push!(arg_list,( Kvect, Rlist, GV_R ) )
            end
        end
    end
    GV_k_list_tmp = pmap(fourier_r2k_internal, arg_list)
    i=0
    for k1 = 1:size(Klist)[1]
        for k2 = 1:size(Klist)[2]
            for k3 = 1:size(Klist)[3]
                i=i+1
                if length(size(GV_R)[4:end]) == 4
                    GV_k[k1,k2,k3,:,:,:,:] = GV_k_list_tmp[i];
                end
                
                if length(size(GV_R)[4:end]) == 3
                    GV_k[k1,k2,k3,:,:,:] = GV_k_list_tmp[i];
                end                   
            end
        end
    end
    return GV_k
end

