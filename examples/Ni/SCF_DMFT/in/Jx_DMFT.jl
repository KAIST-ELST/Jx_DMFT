import DFTforge

import Plots
import FFTW
import JSON
import Dierckx
using JSON
using ImageFiltering
using Dierckx
using Distributed
using Printf
using LinearAlgebra

### EDMFTF related functions


function R2Csphericalharmonics(Trans_mat, Trans_mat_A, Corr_orbital_Ind)

    Trans_mat_C = []
    Trans_mat_A_C = Matrix{ComplexF64}(I,size(Trans_mat_A))
    for i=1:size(Corr_orbital_Ind)[1]
        Trans_mat_C_dum = zeros(ComplexF64,size(Trans_mat[i]))

        if (size(Trans_mat[i])[1]==5)
            Trans_mat_C_dum = zeros(ComplexF64,size(Trans_mat[i]))
            Trans_mat_C_dum = [1/sqrt(2)*im 0.0+0.0*im 0.0+0.0*im 0.0+0.0*im -1/sqrt(2)*im;
                               0.0+0.0*im 1/sqrt(2)*im 0.0+0.0*im 1/sqrt(2)*im  0.0+0.0*im;
                               0.0+0.0*im 0.0+0.0*im 1.0+0.0*im 0.0+0.0*im 0.0+0.0*im;
                               0.0+0.0*im 1/sqrt(2)+0.0*im 0.0+0.0*im -1/sqrt(2)+0.0*im 0.0+0.0*im;
                               1/sqrt(2)+0.0*im 0.0+0.0*im 0.0+0.0*im 0.0+0.0*im 1/sqrt(2)+0.0*im]
        else
            Trans_mat_C_dum += Trans_mat[i]
            println("Caution : Trans_mat_C is set to be same with Trans_mat")
        end
        
        push!(Trans_mat_C,Trans_mat_C_dum)
        Trans_mat_A_C[Corr_orbital_Ind[i],Corr_orbital_Ind[i]] = Trans_mat_C_dum
        
        
    end
    
    return Trans_mat_C, Trans_mat_A_C
end



function EDMFTF_write_trans(Trans_mat_C,Corr_atom_equiv,DMFT_spin_type)
    
    for i in unique(abs.(Corr_atom_equiv))
    
        open(string("imp_",i,"/Trans.dat"), "w") do io    
            @printf(io,"%3i %3i #  size of Sigind and CF",size(Trans_mat_C[1])[1]*DMFT_spin_type, size(Trans_mat_C[1])[1]*DMFT_spin_type)
            @printf(io,"\n")
            @printf(io,"#---- Sigind follows")
            @printf(io,"\n")
            
            ai = findall(x->x==i,abs.(Corr_atom_equiv))[1]
            
            for r = 1:size(Trans_mat_C[ai])[1]*DMFT_spin_type
                for c =1:size(Trans_mat_C[ai])[1]*DMFT_spin_type
                    
                    if (r <(size(imp_ind[i])[1] +1) && c <(size(imp_ind[i])[2] +1))
                        @printf(io,"%3i", imp_ind[i][r,c])
                    elseif (r <(size(imp_ind[i])[1] +1) && c >(size(imp_ind[i])[2] ))
                        @printf(io,"%3i", 0)
                    elseif (r >(size(imp_ind[i])[1]) && c <(size(imp_ind[i])[2] )+1)
                        @printf(io,"%3i", 0)
                    elseif (r >(size(imp_ind[i])[1]) && c >(size(imp_ind[i])[2] ))
                        if (imp_ind[i][r-size(imp_ind[i])[1], c-size(imp_ind[i])[2]] == 0)         
                            @printf(io,"%3i", 0)
                        else
                            @printf(io,"%3i", imp_ind[i][r-size(imp_ind[i])[1], c-size(imp_ind[i])[2]]+maximum(imp_ind[i]))
                        end
                    end
                end
                @printf(io,"\n")
            end
            
            @printf(io,"#---- CF follows")
            @printf(io,"\n")

            
            for r = 1:size(Trans_mat_C[ai])[1]*DMFT_spin_type
                for c =1:size(Trans_mat_C[ai])[1]*DMFT_spin_type
                    
                    if (r <(size(imp_ind[i])[1] +1) && c <(size(imp_ind[i])[2] +1))
                        @printf(io,"%13.9f %13.9f ", real(Trans_mat_C[ai][r,c]), imag(Trans_mat_C[ai][r,c]) )
                    elseif (r <(size(imp_ind[i])[1] +1) && c >(size(imp_ind[i])[2] ))
                        @printf(io,"%13.9f %13.9f ", 0.0, 0.0 )
                    elseif (r >(size(imp_ind[i])[1]) && c <(size(imp_ind[i])[2] )+1)
                        @printf(io,"%13.9f %13.9f ", 0.0, 0.0 )
                    elseif (r >(size(imp_ind[i])[1]) && c >(size(imp_ind[i])[2] ))
                        @printf(io,"%13.9f %13.9f ", real(Trans_mat_C[ai][r-size(imp_ind[i])[1],c-size(imp_ind[i])[2]]), imag(Trans_mat_C[ai][r-size(imp_ind[i])[1],c-size(imp_ind[i])[2]] ) )
                    end
                end
                @printf(io,"\n")
            end            

            
            
        end
        
        
        
    end
    
end



function regulate_Eimp_lev(H_loc, Corr_orbital_Ind, Corr_atom_equiv )
    regulated_Eimp = []
    regulated_ShiftE = []
    for i in unique(abs.(Corr_atom_equiv))
    
        orb_ind = Corr_orbital_Ind[findall(x->x==i,abs.(Corr_atom_equiv))[1]]
        push!(regulated_ShiftE, diag(real(H_loc[orb_ind,orb_ind,1]))[1])
        push!(regulated_Eimp, diag(real(H_loc[orb_ind,orb_ind,1])).-diag(real(H_loc[orb_ind,orb_ind,1]))[1])
    
    end
    return regulated_ShiftE, regulated_Eimp
end

function EDMFTF_write_PARAMS(Corr_atom_equiv,imp_U,imp_J,imp_dc_n,beta,imp_int_type,EDMFTF_MonteCar_step,EDMFTF_GlobalFlip,EDMFTF_warmup_step,mu,regulated_Eimp,regulated_ShiftE,EDMFTF_tsample,EDMFTF_nom,EDMFTF_PChangeOrder)

    for i in unique(abs.(Corr_atom_equiv))
        open(string("imp_",i,"/PARAMS"), "w") do io
            @printf(io,"OffDiagonal  real\n")
            @printf(io,"aom %2i\n", 1)
            @printf(io,"Sig  Sig.out\n")
            @printf(io,"Delta  Delta.inp\n")
            @printf(io,"cix  actqmc.cix\n")

            @printf(io,"Ncout %8i\n",1000000)
            @printf(io,"nom %5i\n",EDMFTF_nom)
            @printf(io,"svd_lmax %3i\n",25)
            @printf(io,"exe  ctqmc\n")
            @printf(io,"tsample %5i\n",EDMFTF_tsample)
            @printf(io,"warmup %14.2f\n", EDMFTF_warmup_step)
            @printf(io,"Ed  [")

            if DMFT_spin_type > 0
                for s = 1:DMFT_spin_type
                    for ii in unique(diag(imp_ind[i]))
                        if ((s != DMFT_spin_type) || (ii != maximum(diag(imp_ind[i])) ))
                            @printf(io,"%14.9f,", regulated_Eimp[i][findall(x->x==ii, diag(imp_ind[i]))[1]]+regulated_ShiftE[i]-mu )
                        else
                            @printf(io,"%14.9f]\n", regulated_Eimp[i][findall(x->x==ii, diag(imp_ind[i]))[1]]+regulated_ShiftE[i]-mu )
                        end
                    end
                end
            end
            
            if EDMFTF_PChangeOrder != 0.9
                @printf(io,"PChangeOrder  %f\n", EDMFTF_PChangeOrder)
            end
            @printf(io,"U  %f\n", imp_U[i])                
            @printf(io,"J  %f\n", imp_J[i])
            @printf(io,"M  %14.2f\n", EDMFTF_MonteCar_step)

            if (imp_int_type[i] == "full" || imp_int_type[i] == "Full" || imp_int_type[i] == "FULL")       
                @printf(io,"CoulombF  '%s'\n", "Full")
            elseif(imp_int_type[i] == "ising" || imp_int_type[i] == "Ising" || imp_int_type[i] == "ISING") 
                @printf(io,"CoulombF  '%s'\n", "Ising")            
            end
                
            @printf(io,"mu  %14.9f\n", -regulated_Eimp[i][findall(x->x==1, diag(imp_ind[i]))[1]]-regulated_ShiftE[i]+mu )
            @printf(io,"beta  %14.9f\n", beta)                
            @printf(io,"mode  SM\n")                
            @printf(io,"GlobalFlip  %i\n", EDMFTF_GlobalFlip)
            @printf(io, "nf0  %f\n", imp_dc_n[i])
            
        end
    end
    
end



#### etc. function

function Shift_H_k(H_k,lev_shift_mat)
    
    H_k_dum = deepcopy(H_k)
    
    if( size(H_k[1,1,1,:,:,1]) == size(lev_shift_mat) )
        for k1 = 1:size(H_k)[1]
            for k2 = 1:size(H_k)[2]
                for k3 = 1:size(H_k)[3]
                    for s =1:size(H_k)[6]
                        H_k_dum[k1,k2,k3,:,:,s] += lev_shift_mat
                    end
                end
            end
        end
        
    else
        error("Hamiltonian H_k and H_R size is not equivalent !!")
    end
    
    
    return H_k_dum
end


function imp_level_shift(imp_lev_shift, H_R0, Corr_atom_Ind, Corr_atom_equiv, imp_ind)
    lev_shift_mat=zeros(ComplexF64, size(H_R0))
    for j = 1:size(Corr_atom_Ind)[1]
        lev_shift_block = zeros(ComplexF64, size(imp_ind[abs(Corr_atom_equiv[j])]))
        lev_shift_block = imp_lev_shift[abs(Corr_atom_equiv[j])]

        
        lev_shift_mat[Corr_orbital_Ind[j],Corr_orbital_Ind[j]]= lev_shift_block;

    end
    return lev_shift_mat    
end


function cal_corr_ineq_orbital_N(imp_ind)
    Corr_ineq_orbital_N = []
    for i = 1:size(imp_ind)[1]
        push!(Corr_ineq_orbital_N,count(x->x>0, unique(imp_ind[i])));
    end
    
    return Corr_ineq_orbital_N
end



function imp_ind_trans_by_degeneracy(imp_ind,H_loc,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind)
    imp_ind_new = []
    for i =1:size(imp_ind)[1]
    
        orb_ind = Corr_orbital_Ind[findall(x->x==i,Corr_atom_equiv)[1]]
        dum_imp_lev = deepcopy(real(H_loc[orb_ind,orb_ind,1]))
        
        for l = 1:size(diag(dum_imp_lev))[1]
            sel_orb=findall(x-> (real(x)!=0) && abs.(real(x)-dum_imp_lev[l,l])<0.003, dum_imp_lev )
            dum_imp_lev[sel_orb].=real(dum_imp_lev[l,l])
        end
        
        dum_imp_block = deepcopy(imp_ind[i])
        
        for (ii,lev) in enumerate(unique(diag(dum_imp_lev)))
            sel_orb=findall(x->x==lev, dum_imp_lev)
            dum_imp_block[sel_orb].= ii
        end
        
        push!(imp_ind_new, dum_imp_block)
    end
    println("")
    println(" == You turned on the 'consider_degeneracy = true', so from now on degenerate states will be considered equivalently == ")
    println("")
    
    println("Old imp ind :")
    for i=1:size(imp_ind)[1]
        println("             imp-",i)
        for j =1:size(imp_ind[i])[1]
            println("                        ",imp_ind[i][j,:]) 
        end
    end
    println("")
    println("New imp ind :")
    for i=1:size(imp_ind)[1]
        println("             imp-",i)
        for j =1:size(imp_ind[i])[1]
            println("                        ",imp_ind_new[i][j,:]) 
        end
    end
    println("")
    println(" == You turn on the consider_degeneracy = true, so from now on degenerate states will be considered equivalently == ")
    println("")
   
    return imp_ind_new
end



function Tvec_print_from_reading(H0, transMat, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)
    Trans_mat = []
    Trans_mat_A = Matrix{Float64}(I,size(H_k)[4],size(H_k)[5])
    for i = 1:size(Corr_atom_Ind)[1]
        orb_ind = Corr_orbital_Ind[i]
        Trans_mat_A[orb_ind,orb_ind] = deepcopy(transMat[i])
    end
    return Trans_mat_A
end


function read_transforMat()
    re = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$";
    transMat=[]
    rowind =[]
    open("transform_Mat.dat", "r") do f 
        Nstring=readlines(f)
        #println(Nstring[1])
        #println(size(Nstring[1])[1])
        for i =1:size(Nstring)[1]
            if( (size(split(Nstring[i]))[1]>0) && occursin(re,split(Nstring[i])[1])) 
                push!(rowind, i)

            elseif(size(split(Nstring[i]))[1]<=0 && (size(rowind)[1]>0))
                #println(rowind)
                transMat_dum=zeros(Float64,size(rowind)[1], convert(Int64,(size(split(Nstring[rowind[end]]))[1])/2))
                    
                for j=1:size(rowind)[1]
                    for k=1:convert(Int64,(size(split(Nstring[rowind[end]]))[1])/2)
                        transMat_dum[j,k]=parse(Float64,split(Nstring[rowind[j]])[2*(k-1)+1])#+parse(Float64,split(Nstring[rowind[j]])[2*(k-1)+2])*im
                    end
                        
                end
                push!(transMat,transMat_dum)
                rowind =[]
            end

            
        end
        #println(split(Nstring[3]))
    end
    
    return transMat

end



function read_DFT_corrNele()
    dft_corrN_ = []
    open("DFT_CorrNele.dat", "r") do f
        Nstring=readlines(f)
        for i=1:size(Nstring)[1]
            push!(dft_corrN_, parse(Float64, split(Nstring[i])[4]))
        end
    end

    return dft_corrN_
end


function print_DFT_corrNele(Ineq_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Occup)
    dft_corrN = []
    open("DFT_CorrNele.dat", "w") do io
        for it in Ineq_atom_Ind
            orb_ind = Corr_orbital_Ind[findall(x->x==it, Corr_atom_equiv)][1]
            @printf(io, "latt. %2i :", it )
            @printf(io," %10.6f\n", sum(diag(real(Occup[orb_ind,orb_ind,1]))) * 2.0)
            push!(dft_corrN, sum(diag(real(Occup[orb_ind,orb_ind,1]))) * 2.0)
        end
    end 
    return dft_corrN
end


function beta_cal(Temp)
    BoltzmannC = 8.617333262145/100000;
    beta_=1/(BoltzmannC*Temp)
    return beta_
end


function FLL_dc(imp_dc_n, imp_U, imp_J)
    FLL_dc=[];
    for i=1:size(imp_dc_n)[1]
        dc=imp_U[i].*(imp_dc_n[i].-1/2).-imp_J[i]/2.0.*(imp_dc_n[i].-1);
        push!(FLL_dc, dc)
    end
    return FLL_dc
end


function dc_change(imp_totN, imp_dc_n, dft_corrN, imp_dc_type, imp_U, imp_J)
    if (imp_dc_type == "nominal" || imp_dc_type == "Norminal" || imp_dc_type == "NORMINAL")
        dc = FLL_dc(imp_dc_n, imp_U, imp_J)
    elseif (imp_dc_type == "fll" || imp_dc_type == "Fll" || imp_dc_type == "FLL")
        dc = FLL_dc(imp_totN, imp_U, imp_J)
    elseif (imp_dc_type == "FLL-DFT" || imp_dc_type == "fll-dft" || imp_dc_type == "Fll-dft")
        dc = FLL_dc(dft_corrN, imp_U, imp_J)
    end

    return dc
end


function go_bareSelfE(num,Ineq_atom_Ind,DMFT_solver)
    det_bare_selfE=[]
    for i in Ineq_atom_Ind
        if DMFT_solver == "ComCTQMC"
            file_n=string("imp_",i,"/params.obs.json")
        elseif DMFT_solver == "EDMFTF_ctqmc"
            file_n=string("imp_",i,"/Sig.out")
        end
        yn=maxtag_larger(file_n,num)
        push!(det_bare_selfE,yn)
    end
    return prod(det_bare_selfE)
end


function maxtag_larger(file_n,num)
    file_n_pre = file_n[1:findlast(".",file_n)[1]-1]
    file_n_suf = file_n[findlast(".",file_n)[1]:end]
    tag=[]
    for i=1:300
        f_str=string(file_n_pre,"_",i,file_n_suf)
        if isfile(f_str)
           push!(tag,i) 
        end
    end
    return (maximum(tag)>=num)
end


function readmu_fromscflog()
    mulog=[]
    for i =1:size(readlines("scf.log"))[1]
        if (size(split(readlines("scf.log")[i]))[1]>4)
            numchar = split(readlines("scf.log")[i])[4]
            if isa(tryparse(Float64,numchar), Number)
                mu = parse(Float64,numchar)
                push!(mulog,mu)
            end
        end
    end
    return mulog
    
end


function read_mu_NonIntH()
    open("mu_NonIntH.dat", "r") do f
        Nstring=read(f,String)
    end
end

function write_mu_NonIntH(mu)
    open("mu_NonIntH.dat", "w") do io
        @printf(io," %30.24f", mu)
    end
end


function read_Nele()
    open("Nele.dat", "r") do f
        Nstring=read(f,String)
    end
end

function write_Nele(targetN)
    open("Nele.dat", "w") do io
        @printf(io," %10.6f", targetN)
    end
end


function init_trans_mat(Corr_atom_Ind, Corr_orbital_Ind)
    Trans_mat = []
    for i = 1:size(Corr_atom_Ind)[1]
        orbsize = size( Corr_orbital_Ind[i] )[1]
        Tmat = Matrix{Float64}(I,orbsize,orbsize)
        push!(Trans_mat,Tmat)
        
    end
    return Trans_mat
end



function Trans_H_k(H_k,Trans_mat_A)
    dum_H_k = zeros(size(H_k))+ zeros(size(H_k))*im
    
    for k1 = 1:size(H_k)[1]
        for k2 = 1:size(H_k)[2]
            for k3 = 1:size(H_k)[3]
                for s= 1:size(H_k)[6]
                    dum_H_k[k1,k2,k3,:,:,s] = Trans_mat_A * H_k[k1,k2,k3,:,:,s] * Trans_mat_A'
    
                end
            end
        end
    end
    return dum_H_k
end




function Hloc_print(filename, H_loc ,Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)
    #H_loc_rot = []
    
    file_n = filename
    #file_n2 = string("Hloc_rot.dat")
    open(file_n, "w") do io
        for i = 1:size(Corr_atom_Ind)[1]
            @printf(io,"atom[%2i]\n\n",i)
            orb_ind = Corr_orbital_Ind[i]
            Hloc = H_loc[orb_ind,orb_ind,1]
            for j = 1:size(Hloc)[1]
                for jj = 1:size(Hloc)[2]
                    if jj != size(Hloc)[2]
                        @printf(io,"%14.8f %14.8f     ", real(Hloc[j,jj]),imag(Hloc[j,jj])  )
                    else
                        @printf(io,"%14.8f %14.8f\n", real(Hloc[j,jj]),imag(Hloc[j,jj])  )
                    end
                end
            end
            @printf(io,"\n\n"  )
        end
    end
    
    #open(file_n2, "w") do io
    #    for i = 1:size(Corr_atom_Ind)[1]
    #        @printf(io,"atom[%2i]\n\n",i)

    #        for j = 1:size(H_loc_rot[i])[1]
    #            for jj = 1:size(H_loc_rot[i])[2]
    #                if jj != size(H_loc_rot[i])[2]
    #                    @printf(io,"%14.8f %14.8f     ", real(H_loc_rot[i][j,jj]),imag(H_loc_rot[i][j,jj])  )
    #                else
    #                    @printf(io,"%14.8f %14.8f\n", real(H_loc_rot[i][j,jj]),imag(H_loc_rot[i][j,jj])  )
    #                end
    #            end
    #        end
    #        @printf(io,"\n\n"  )
    #    end
    #end    
    
    
    #return H_loc_rot
end



function Spec_atEf_print(Spec_atEf, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)
    
    file_n = string("Spec_atEf.dat")
    open(file_n, "w") do io
        for i = 1:size(Corr_atom_Ind)[1]
            @printf(io,"atom[%2i]\n\n",i)
            orb_ind = Corr_orbital_Ind[i]
            Spec_ = Spec_atEf[orb_ind,orb_ind,1]
            for j = 1:size(Spec_)[1]
                for jj = 1:size(Spec_)[2]
                    if jj != size(Spec_)[2]
                        @printf(io,"%14.8f %14.8f     ", real(Spec_[j,jj]),imag(Spec_[j,jj])  )
                    else
                        @printf(io,"%14.8f %14.8f\n", real(Spec_[j,jj]),imag(Spec_[j,jj])  )
                    end
                end
            end
            @printf(io,"\n\n"  )
        end
    end
    
end



function Tvec_print(Aw, nomal_implev, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind, basis_transform)
    println("===== basis transformation to diagonalizing sepc/hyb ... ======")

    Trans_mat = []
    Trans_mat_A = Matrix{Float64}(I,size(Aw)[1],size(Aw)[2])
    file_n = string("transform_Mat.dat")
    open(file_n, "w") do io
        for i = 1:size(Corr_atom_Ind)[1]
            @printf(io,"atom[%2i]\n\n",i)
            orb_ind = Corr_orbital_Ind[i]
            
            if (basis_transform == "spec")
                Tvec = eigvecs(real(Aw[orb_ind,orb_ind]))'
                Trans_mat_A[orb_ind,orb_ind] = deepcopy(Tvec)
                push!(Trans_mat, Tvec)
                for j = 1:size(Tvec)[1]
                    for jj = 1:size(Tvec)[2]
                        if jj != size(Tvec)[2]
                            @printf(io,"%20.16f %20.16f   ", Tvec[j,jj], 0.0  )
                        else
                            @printf(io,"%20.16f %20.16f\n", Tvec[j,jj], 0.0  )
                        end
                    end
                end            
                @printf(io,"\n\n"  )
            elseif (basis_transform == "hyb")
                Tvec = eigvecs(real(nomal_implev[i]))'
                Trans_mat_A[orb_ind,orb_ind] = deepcopy(Tvec)
                push!(Trans_mat, Tvec)
                for j = 1:size(Tvec)[1]
                    for jj = 1:size(Tvec)[2]
                        if jj != size(Tvec)[2]
                            @printf(io,"%20.16f %20.16f   ", Tvec[j,jj], 0.0  )
                        else
                            @printf(io,"%20.16f %20.16f\n", Tvec[j,jj], 0.0  )
                        end
                    end
                end            
                @printf(io,"\n\n"  )                
                
            end
            
        end
    end
    return Trans_mat, Trans_mat_A
end



function Tvec_print_I(H_k, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)
    Trans_mat = []
    Trans_mat_A = Matrix{Float64}(I,size(H_k)[4],size(H_k)[5])
    file_n = string("transform_Mat.dat")
    open(file_n, "w") do io
        for i = 1:size(Corr_atom_Ind)[1]
            @printf(io,"atom[%2i]\n\n",i)
            orb_ind = Corr_orbital_Ind[i]
            Tvec = deepcopy(Trans_mat_A[orb_ind,orb_ind])

            push!(Trans_mat, Tvec)
            for j = 1:size(Tvec)[1]
                for jj = 1:size(Tvec)[2]
                    if jj != size(Tvec)[2]
                        @printf(io,"%20.16f %20.16f   ", Tvec[j,jj], 0.0  )
                    else
                        @printf(io,"%20.16f %20.16f\n", Tvec[j,jj], 0.0  )
                    end
                end
            end
            @printf(io,"\n\n"  )
        end
    end
    return Trans_mat, Trans_mat_A
end




function input_handler(DMFT_Jx_option)

   if (DMFT_Jx_option.Calculation_mode == "DMFT")
        Calculation_mode = DMFT_Jx_option.Calculation_mode
        println(" Calculation mode :", Calculation_mode )
        ## Constant
        BoltzmannC = 8.617333262145/100000;
        delta =1e-10; # for grennf occupation count


        ### input parameter
        K_grid_num = DMFT_Jx_option.KgridNum;
        R_grid_num = DMFT_Jx_option.RgridNum;
        iW_grid_cut = DMFT_Jx_option.iWgridCut;
        Tem= DMFT_Jx_option.Temperature;
        beta = beta_cal(Tem);
        beta_for_occ = beta_cal(100.0);

        DMFT_spin_type = DMFT_Jx_option.Spin_type;
        DMFT_solver = DMFT_Jx_option.DMFT_solver;


        Corr_atom_Ind = DMFT_Jx_option.Corr_atom_Ind;
        Corr_orbital_Ind = DMFT_Jx_option.Corr_orbital_Ind;
        Corr_atom_equiv = DMFT_Jx_option.Corr_atom_equiv;

        Ineq_atom_Ind=unique(abs.(Corr_atom_equiv));

        if DMFT_spin_type == 2
            mag_order = ones(Int64, size(Corr_atom_equiv)[1]);
            mag_order[findall(x->x <0, Corr_atom_equiv)].=-1;
        end

        DMFT_loop_N = DMFT_Jx_option.DMFT_loop_N;


        Mix_selfE = DMFT_Jx_option.Mix_selfE;
        init_bias = DMFT_Jx_option.init_bias;
        smth_step = DMFT_Jx_option.smth_step;
        cal_susc = DMFT_Jx_option.cal_susc;

        imp_ind = DMFT_Jx_option.imp_block;
        imp_lev_shift = DMFT_Jx_option.imp_lev_shift;
        imp_U =  DMFT_Jx_option.imp_U;
        imp_J =  DMFT_Jx_option.imp_J;
        imp_dc_n =  DMFT_Jx_option.imp_dc;
        imp_dc = FLL_dc(imp_dc_n, imp_U, imp_J);
        #imp_U.*(imp_dc_n.-1/2).-imp_J/2.0.*(imp_dc_n.-1);
        imp_dc_type = DMFT_Jx_option.imp_dc_type;
        imp_int_type =  DMFT_Jx_option.imp_int_type;
        imp_int_parameterisation =  DMFT_Jx_option.imp_int_parameterisation;
        
        
        Corr_ineq_orbital_Num = []
        for i = 1:size(imp_ind)[1]
          push!(Corr_ineq_orbital_Num,count(x->x>0, unique(imp_ind[i])));
        end
    
        Solver_green_Cut = DMFT_Jx_option.Solver_green_Cut;


        F0 = imp_U;
        F2 = imp_J*14/2.6*1.6;
        F4 = imp_J*14/2.6*1.0;

        imp_Measure_time = DMFT_Jx_option.imp_Measure_time;
        imp_Thermal_time = DMFT_Jx_option.imp_Thermal_time;
        
        
        basis_transform = DMFT_Jx_option.basis_transform;
        consider_degeneracy =  DMFT_Jx_option.consider_degeneracy;
        
        green_basis = DMFT_Jx_option.green_basis;
        green_legendre_cutoff = DMFT_Jx_option.green_legendre_cutoff;
        
        EDMFTF_MonteCar_step = DMFT_Jx_option.EDMFTF_MonteCar_step
        EDMFTF_warmup_step = DMFT_Jx_option.EDMFTF_warmup_step
        EDMFTF_GlobalFlip = DMFT_Jx_option.EDMFTF_GlobalFlip
        EDMFTF_tsample = DMFT_Jx_option.EDMFTF_tsample
        EDMFTF_nom = DMFT_Jx_option.EDMFTF_nom
        EDMFTF_PChangeOrder = DMFT_Jx_option.EDMFTF_PChangeOrder
        
        compute_EF = DMFT_Jx_option.compute_EF
        
        return (Calculation_mode, BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut, Tem, beta, beta_for_occ,
                DMFT_spin_type, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind,
                DMFT_loop_N, Mix_selfE, init_bias,smth_step,cal_susc,compute_EF, DMFT_solver, imp_dc_type, imp_ind, imp_lev_shift, imp_U, imp_J, imp_dc_n, imp_dc, imp_int_type, imp_int_parameterisation
                , Corr_ineq_orbital_Num, Solver_green_Cut ,F0, F2, F4, imp_Measure_time, imp_Thermal_time, basis_transform,consider_degeneracy, green_basis, green_legendre_cutoff
                ,EDMFTF_MonteCar_step, EDMFTF_warmup_step, EDMFTF_GlobalFlip, EDMFTF_tsample, EDMFTF_nom, EDMFTF_PChangeOrder );
        
        
    elseif (DMFT_Jx_option.Calculation_mode == "Jx-DMFT")
        
        Calculation_mode = DMFT_Jx_option.Calculation_mode
        println(" Calculation mode :", Calculation_mode)
        ## Constant
        BoltzmannC = 8.617333262145/100000;
        delta =1e-10; # for grennf occupation count


        ### input parameter
        K_grid_num = DMFT_Jx_option.KgridNum;
        R_grid_num = DMFT_Jx_option.RgridNum;
        iW_grid_cut = DMFT_Jx_option.iWgridCut;
        Tem= DMFT_Jx_option.Temperature;
        beta = beta_cal(Tem);
        beta_for_occ = beta_cal(100.0);


        DMFT_spin_type = DMFT_Jx_option.Spin_type;
        DMFT_solver = DMFT_Jx_option.DMFT_solver


        Corr_atom_Ind = DMFT_Jx_option.Corr_atom_Ind;
        Corr_orbital_Ind = DMFT_Jx_option.Corr_orbital_Ind;
        Corr_atom_equiv = DMFT_Jx_option.Corr_atom_equiv;

        Ineq_atom_Ind=unique(abs.(Corr_atom_equiv));

        if DMFT_spin_type == 2
            mag_order = ones(Int64, size(atom12_list)[1]);
            for (a12ind, a12) in enumerate(atom12_list)
                a1ind = findall(x->x==a12[1],Corr_atom_Ind)[1]
                a2ind = findall(x->x==a12[2],Corr_atom_Ind)[1]
                
                if (Corr_atom_equiv[a1ind]*Corr_atom_equiv[a2ind])<0
                    mag_order[a12ind]=-1;
                end
            end
        end

        DMFT_loop_N = DMFT_Jx_option.DMFT_loop_N;


        Mix_selfE = DMFT_Jx_option.Mix_selfE;
        init_bias = DMFT_Jx_option.init_bias;
        smth_step = DMFT_Jx_option.smth_step;
        cal_susc = DMFT_Jx_option.cal_susc;

        imp_ind = DMFT_Jx_option.imp_block;
        imp_lev_shift = DMFT_Jx_option.imp_lev_shift;
        imp_U =  DMFT_Jx_option.imp_U;
        imp_J =  DMFT_Jx_option.imp_J;
        imp_dc_n =  DMFT_Jx_option.imp_dc;
        imp_dc_type = DMFT_Jx_option.imp_dc_type;
        imp_dc = FLL_dc(imp_dc_n, imp_U, imp_J);
        imp_int_type =  DMFT_Jx_option.imp_int_type;
        imp_int_parameterisation =  DMFT_Jx_option.imp_int_parameterisation;

        Corr_ineq_orbital_Num = []
        for i = 1:size(imp_ind)[1]
          push!(Corr_ineq_orbital_Num,count(x->x>0, unique(imp_ind[i])));
        end
    
        Solver_green_Cut = DMFT_Jx_option.Solver_green_Cut;


        F0 = imp_U;
        F2 = imp_J*14/2.6*1.6;
        F4 = imp_J*14/2.6*1.0;

        imp_Measure_time = DMFT_Jx_option.imp_Measure_time;
        imp_Thermal_time = DMFT_Jx_option.imp_Thermal_time;
        
        
        basis_transform = DMFT_Jx_option.basis_transform;
        consider_degeneracy =  DMFT_Jx_option.consider_degeneracy;
        
        green_basis = DMFT_Jx_option.green_basis;
        green_legendre_cutoff = DMFT_Jx_option.green_legendre_cutoff;
        
                
        EDMFTF_MonteCar_step = DMFT_Jx_option.EDMFTF_MonteCar_step
        EDMFTF_warmup_step = DMFT_Jx_option.EDMFTF_warmup_step
        EDMFTF_GlobalFlip = DMFT_Jx_option.EDMFTF_GlobalFlip
        EDMFTF_tsample = DMFT_Jx_option.EDMFTF_tsample
        EDMFTF_nom = DMFT_Jx_option.EDMFTF_nom
        EDMFTF_PChangeOrder = DMFT_Jx_option.EDMFTF_PChangeOrder
        
        compute_EF = DMFT_Jx_option.compute_EF
        
        
        return (Calculation_mode,BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut, Tem, beta, beta_for_occ,
                DMFT_spin_type, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind,
                DMFT_loop_N, Mix_selfE, init_bias,smth_step,cal_susc, compute_EF, DMFT_solver, imp_dc_type, imp_ind, imp_lev_shift, imp_U, imp_J, imp_dc_n, imp_dc, imp_int_type, imp_int_parameterisation
                , Corr_ineq_orbital_Num, Solver_green_Cut ,F0, F2, F4, imp_Measure_time, imp_Thermal_time, basis_transform,consider_degeneracy, green_basis, green_legendre_cutoff,mag_order
                ,EDMFTF_MonteCar_step, EDMFTF_warmup_step, EDMFTF_GlobalFlip, EDMFTF_tsample, EDMFTF_nom, EDMFTF_PChangeOrder );

        
    elseif (DMFT_Jx_option.Calculation_mode == "Jx0")
        Calculation_mode = DMFT_Jx_option.Calculation_mode
        println(" Calculation mode :", Calculation_mode)

        BoltzmannC = 8.617333262145/100000;
        delta =1e-10; # for grennf occupation count


        ### input parameter
        K_grid_num = DMFT_Jx_option.KgridNum;
        R_grid_num = DMFT_Jx_option.RgridNum;
        iW_grid_cut = DMFT_Jx_option.iWgridCut;
        
        Tem= DMFT_Jx_option.Temperature;
        beta = beta_cal(Tem);

        Corr_atom_Ind = DMFT_Jx_option.Corr_atom_Ind;
        Corr_orbital_Ind = DMFT_Jx_option.Corr_orbital_Ind;
        Corr_atom_equiv = DMFT_Jx_option.Corr_atom_equiv;

        Ineq_atom_Ind=unique(abs.(Corr_atom_equiv));

        imp_ind=[]
        Corr_ineq_orbital_Num = []
        for i = 1:size(Ineq_atom_Ind)[1]
            a =  Matrix{Int64}(I,size(Corr_orbital_Ind[i])[1],size(Corr_orbital_Ind[i])[1])
            push!(Corr_ineq_orbital_Num,size(Corr_orbital_Ind[i])[1]);
            for j = 1: size(Corr_orbital_Ind[i])[1]
               a[j,j]=j 
            end
            push!(imp_ind,a)
        end
        #imp_dc_n =  DMFT_Jx_opion.imp_dc;
        #imp_dc = imp_U.*(imp_dc_n.-1/2).-imp_J/2.0.*(imp_dc_n.-1);
    
        #mag_order = ones(Int64, size(Corr_atom_equiv)[1]);
        #mag_order[findall(x->x <0, Corr_atom_equiv)].=-1;        
        mag_order = ones(Int64, size(atom12_list)[1]);
        for (a12ind, a12) in enumerate(atom12_list)
            a1ind = findall(x->x==a12[1],Corr_atom_Ind)[1]
            a2ind = findall(x->x==a12[2],Corr_atom_Ind)[1]

            if (Corr_atom_equiv[a1ind]*Corr_atom_equiv[a2ind])<0
                mag_order[a12ind]=-1;
            end
        end        
            

        
        DMFT_spin_type = 0;
        
        

        
        return (Calculation_mode, BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut,Tem, beta,DMFT_spin_type,
                 Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind, Corr_ineq_orbital_Num,imp_ind,mag_order);
        

        
    elseif (DMFT_Jx_option.Calculation_mode == "Magnon")
        Calculation_mode = DMFT_Jx_option.Calculation_mode
        println(" Calculation mode :", Calculation_mode)

        qmode = DMFT_Jx_option.qmode;
        Neighbors_cut = DMFT_Jx_option.Neighbors_cut
        BoltzmannC = 8.617333262145/100000;
        delta =1e-10; # for grennf occupation count


        ### input parameter
        K_grid_num = DMFT_Jx_option.KgridNum;
        R_grid_num = DMFT_Jx_option.RgridNum;
        iW_grid_cut = DMFT_Jx_option.iWgridCut;
        
        Tem= DMFT_Jx_option.Temperature;
        beta = beta_cal(Tem);

        Corr_atom_Ind = DMFT_Jx_option.Corr_atom_Ind;
        Corr_orbital_Ind = DMFT_Jx_option.Corr_orbital_Ind;
        Corr_atom_equiv = DMFT_Jx_option.Corr_atom_equiv;

        Ineq_atom_Ind=unique(abs.(Corr_atom_equiv));

        imp_ind=[]
        Corr_ineq_orbital_Num = []
        for i = 1:size(Ineq_atom_Ind)[1]
            a =  Matrix{Int64}(I,size(Corr_orbital_Ind[i])[1],size(Corr_orbital_Ind[i])[1])
            push!(Corr_ineq_orbital_Num,size(Corr_orbital_Ind[i])[1]);
            for j = 1: size(Corr_orbital_Ind[i])[1]
               a[j,j]=j 
            end
            push!(imp_ind,a)
        end
        #imp_dc_n =  DMFT_Jx_opion.imp_dc;
        #imp_dc = imp_U.*(imp_dc_n.-1/2).-imp_J/2.0.*(imp_dc_n.-1);
    
        #mag_order = ones(Int64, size(Corr_atom_equiv)[1]);
        #mag_order[findall(x->x <0, Corr_atom_equiv)].=-1;        
        mag_order = ones(Int64, size(atom12_list)[1]);
        for (a12ind, a12) in enumerate(atom12_list)
            a1ind = findall(x->x==a12[1],Corr_atom_Ind)[1]
            a2ind = findall(x->x==a12[2],Corr_atom_Ind)[1]

            if (Corr_atom_equiv[a1ind]*Corr_atom_equiv[a2ind])<0
                mag_order[a12ind]=-1;
            end
        end        
            

        
        DMFT_spin_type = 0;
        
        

        
        return (Calculation_mode, qmode,Neighbors_cut, BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut,Tem, beta,DMFT_spin_type,
                 Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind, Corr_ineq_orbital_Num,imp_ind,mag_order);
        
        
    end
         
        
end




function Mix_SelfE_(SelfE_w,SelfE_w_new,Mix_selfE)
    SelfE_dum = deepcopy(SelfE_w)
    
    for i = 1:size(SelfE_w)[1]
        SelfE_dum[i] = (Mix_selfE) .* SelfE_w[i] .+ (1.0 - Mix_selfE) * SelfE_w_new[i]
    end
    
    return SelfE_dum
end




function mu_guess(DMFT_loop, mulog, DMFT_spin_type)
    
    if DMFT_spin_type <2
        if DMFT_loop < 2
            guess_mu_range = [-15.0, 15.0]
        else
            del = abs(mulog[end]-mulog[end-1])
            if del>0.2
                guess_mu_range = [mulog[end]-del*5.0, mulog[end]+del*5.0 ]
            else
                guess_mu_range = [mulog[end]-0.7, mulog[end]+0.7 ]
            end
        end
    
    else
        if DMFT_loop < 2
            guess_mu_range = [-15.0, 15.0]
        else
            del = abs(mulog[end]-mulog[end-1])
            if del>0.2
                guess_mu_range = [mulog[end]-del*5.0, mulog[end]+del*5.0 ]
            else
                guess_mu_range = [mulog[end]-0.7, mulog[end]+0.7 ]
            end
        end        
    end
    
    
    return guess_mu_range 
end


function change_mu_paramsjson(Ineq_atom_Ind,imp_dc,mu,regulated_ShiftE)
    for i in Ineq_atom_Ind
        str =string("./imp_",i,"/params.json")
        par=readJSON(str)

        par["mu"] = -regulated_ShiftE[i]+mu #imp_dc[i]+mu
        str =string("imp_",i)
        params_string = json(par, 4)
        open(string(str,"/","params.json"),"w") do f
            print(f, params_string)
        end           
    end
end



function write_occlog_imp(Ineq_atom_Ind, imp_ind, DMFT_spin_type, Corr_ineq_orbital_Num, imp_occ, DMFT_solver)
    open("occ.log", "a") do io   #### HEEEERRRRREE

        for i in Ineq_atom_Ind        
            occup_mat = zeros(Float64,size(imp_ind[i])[1],size(imp_ind[i])[1],DMFT_spin_type)
            if DMFT_spin_type == 1
                occup_mat[:,:,1] = deepcopy(imp_ind[i])
            elseif DMFT_spin_type == 2
                occup_mat[1:size(imp_ind[i])[1],1:size(imp_ind[i])[1],1] = deepcopy(imp_ind[i])
                occup_mat[1:size(imp_ind[i])[1],1:size(imp_ind[i])[1],2] = deepcopy(imp_ind[i])
            end
            occup_mat = deepcopy(convert(Array{Float64,3}, occup_mat))
            @printf(io, "imp. %3i :", i )
            for or = 1:Corr_ineq_orbital_Num[i]
                for s = 1:DMFT_spin_type
                    o_ind = findall(x->x==convert(Float64,or),occup_mat[:,:,s])
                    if DMFT_solver == "ComCTQMC"
                        occup_mat[o_ind,s].= imp_occ[i][or,s]
                    elseif DMFT_solver == "EDMFTF_ctqmc"
                        occup_mat[o_ind,s].= imp_occ[i][or,s]./(size(o_ind)[1])
                        #println("or:",or,", o_ind:",size(o_ind)[1])
                    end
                end
            end
            
            for or =1:size(imp_ind[i])[1]
                if DMFT_spin_type == 1
                    if DMFT_solver == "ComCTQMC"
                        if or != size(imp_ind[i])[1]
                            @printf(io," %10.6f,", occup_mat[or,or,1] * 2.0)
                        else
                            @printf(io," %10.6f,  n : %10.6f\n", occup_mat[or,or,1] * 2.0, sum( diag(occup_mat[:,:,1]) )*2.0  )
                        end
                    elseif DMFT_solver == "EDMFTF_ctqmc"
                        if or != size(imp_ind[i])[1]
                            @printf(io," %10.6f,", occup_mat[or,or,1] )
                        else
                            @printf(io," %10.6f,  n : %10.6f\n", occup_mat[or,or,1] , sum( diag(occup_mat[:,:,1]) ) )
                        end
                    end
                elseif DMFT_spin_type ==2
                    if or != size(imp_ind[i])[1]
                        @printf(io," %10.6f,", occup_mat[or,or,1] +occup_mat[or,or,2])
                    else
                        @printf(io," %10.6f,  [↑ : %5.3f, ↓ : %5.3f, m : %5.3f, n : %10.6f]\n", occup_mat[or,or,1] +occup_mat[or,or,2], sum( diag(occup_mat[:,:,1]) ), sum( diag(occup_mat[:,:,2]) ), sum( diag(occup_mat[:,:,1]) )-sum( diag(occup_mat[:,:,2]) ), sum( diag(occup_mat[:,:,1]) )+sum( diag(occup_mat[:,:,2]) )     )   
                    end
                end
            end
        end
        @printf(io, "-------------------------------------------------------------------------------------------\n\n\n",  )
   
        
    end
                    
    
end




function write_susclog_fromObs(DMFT_loop,Ineq_atom_Ind,DMFT_solver)
    open("susc.log", "a") do io
        @printf(io, "Iter.[%3i]|   ", DMFT_loop)
        
        if DMFT_solver == "ComCTQMC"
            for it in Ineq_atom_Ind
                str =string("./imp_",it,"/params.obs.json")
                obser=readJSON(str)
            
                if it != Ineq_atom_Ind[end]
                    @printf(io,"(imp %2i)   \"N\": %10.6f,   \"Sz\": %10.6f,     ",it, obser["partition"]["susceptibility"]["N"]["function"][1], obser["partition"]["susceptibility"]["Sz"]["function"][1])
                else
                    @printf(io,"(imp %2i)   \"N\": %10.6f,   \"Sz\": %10.6f |\n",it, obser["partition"]["susceptibility"]["N"]["function"][1], obser["partition"]["susceptibility"]["Sz"]["function"][1])
                end
            end
        elseif DMFT_solver == "EDMFTF_ctqmc"
            for it in Ineq_atom_Ind
                str =string("./imp_",it,"/ctqmc.log")
                  
               open(str, "r") do f
                    Nstring=readlines(f)
    
                    pat = pat = r"[+-]?\d+\.?\d*"
                    txt = Nstring[end]  
                    numbers = [parse(Float64, m.match) for m in eachmatch(pat, txt)]
 
                
                    if it != Ineq_atom_Ind[end]
                   
                        @printf(io,"(imp %2i)   \"N\": %10.6f,   \"Sz\": %10.6f,     ",it, numbers[4], numbers[3])
               
                    else
                    
                        @printf(io,"(imp %2i)   \"N\": %10.6f,   \"Sz\": %10.6f |\n",it, numbers[4], numbers[3])
                
                    end
               
                end
            end            

        end
    end
end



function load_suscw0(DMFT_loop,Ineq_atom_Ind,DMFT_solver)
    susc_N=[]
    susc_Sz=[]
    if DMFT_solver == "ComCTQMC"
        for it in Ineq_atom_Ind
            str =string("./imp_",it,"/params.obs.json")
            obser=readJSON(str)
            push!(susc_N, obser["partition"]["susceptibility"]["N"]["function"][1])
            push!(susc_Sz,obser["partition"]["susceptibility"]["Sz"]["function"][1])
        end
    elseif DMFT_solver == "EDMFTF_ctqmc"
        for it in Ineq_atom_Ind
            str =string("./imp_",it,"/ctqmc.log")
            open(str, "r") do f
                Nstring=readlines(f)
    
                pat = pat = r"[+-]?\d+\.?\d*"
                txt = Nstring[end]  
                numbers = [parse(Float64, m.match) for m in eachmatch(pat, txt)]

                #println(numbers)
                push!(susc_N, numbers[4])
                push!(susc_Sz,numbers[3])
            end
        end
    end
    return susc_N,susc_Sz
end







function loadOccup_from_obs(Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,DMFT_solver)
    imp_occ =[]
    
    if DMFT_solver == "ComCTQMC"
        for i in Ineq_atom_Ind
            str =string("./imp_",i,"/params.obs.json")
            obser=readJSON(str)
            imp_occ_block = zeros(Corr_ineq_orbital_Num[i],DMFT_spin_type)
            if DMFT_spin_type > 0
                for s = 1:DMFT_spin_type
                    for l=1:Corr_ineq_orbital_Num[i]
                        imp_occ_block[l,s] = obser["partition"]["occupation"][string(Corr_ineq_orbital_Num[i] * (s-1) + l)][1]
                    end
                end
            end
        push!(imp_occ, imp_occ_block)
    end
    elseif DMFT_solver == "EDMFTF_ctqmc"
        for i in Ineq_atom_Ind
            str =string("./imp_",i,"/Gf.out")
            imp_occ_block = zeros(Corr_ineq_orbital_Num[i],DMFT_spin_type)
            open(str, "r") do f
                Nstring=readlines(f)
            
                pat = pat = r"[+-]?\d+\.?\d*"
                txt = split(Nstring[1])[end]
                numbers = [parse(Float64, m.match) for m in eachmatch(pat, txt)]
            
                if DMFT_spin_type > 0
                    for s = 1:DMFT_spin_type
                        for l=1:Corr_ineq_orbital_Num[i]
                            imp_occ_block[l,s] = numbers[l+(s-1)*Corr_ineq_orbital_Num[i]] 
                        end
                    end
                end
            end
            push!(imp_occ, imp_occ_block)        
        end
    end
    return imp_occ
end


function loadOccup_from_obs_scalar(Ineq_atom_Ind, DMFT_solver)
    imp_occ =[]
    if DMFT_solver == "ComCTQMC"
        for i in Ineq_atom_Ind
            str =string("./imp_",i,"/params.obs.json")
            obser=readJSON(str)
            tN = obser["partition"]["scalar"]["N"][1]
            push!(imp_occ, tN)
        end
    elseif DMFT_solver == "EDMFTF_ctqmc"
        for i in Ineq_atom_Ind
            str =string("./imp_",i,"/Gf.out")
            open(str, "r") do f
                Nstring=readlines(f)
                tN=parse(Float64,split(split(Nstring[1])[2],'=')[2])  
                push!(imp_occ, tN)
            end

        end
    end
    return imp_occ
end



function run_DMFT(DMFT_solver, DMFT_spin_type, Ineq_atom_Ind, imp_J, imp_ind, imp_int_type, mpi_prefix, regulated_Eimp)
    
    if (DMFT_solver == "ComCTQMC")
        for i in Ineq_atom_Ind
            println("\n\n")
            print(" CTQMC run ... imp[",i,"]\n")
            mpi_prefix = string(mpi_prefix)
            excut_file ="CTQMC"
            input_file = string("./imp_",i,"/params")
        
            run_string = split(mpi_prefix)
            push!(run_string, excut_file)
            push!(run_string, input_file)
        
            run(`$run_string`)
            GC.gc()
    
            cp("hloc.json",string("./imp_",i,"/hloc.json"); force=true)
      
            println("\n\n")
            print(" EVALSIM run ... imp[",i,"]\n")
            excut_file ="EVALSIM"
        
            #run_string =[]
            run_string = split(mpi_prefix)
            push!(run_string, excut_file)
            push!(run_string, input_file)
       
            run(`$run_string`)
            GC.gc()
   
        end
        
        
    elseif (DMFT_solver == "EDMFTF_ctqmc")
         
        for i in Ineq_atom_Ind
            
            println("\n\n")
            print(" CTQMC run ... imp[",i,"]\n")
            str = string("imp_",i)
            cd(str)

            excut_file = "atom_d.py"

            run_string = [excut_file]
            push!(run_string, "qOCA=1")
            push!(run_string, "Eoca=1.0")
            push!(run_string, "l=2")
            push!(run_string, "cx=0.0")
            push!(run_string, "Ex=[0.5, 0.5, 0.5]")
            push!(run_string, "Ep=[3.0, 3.0, 3.0]")
            push!(run_string, "mOCA=0.001")
            push!(run_string, "OCA_G=False")
            push!(run_string, "qatom=0")
            if (imp_int_type[i] == "full" || imp_int_type[i] == "Full" || imp_int_type[i] == "FULL")
                push!(run_string, "CoulombF='Full'")
            elseif (imp_int_type[i] == "ising" || imp_int_type[i] == "Ising" || imp_int_type[i] == "ISING")
                push!(run_string, "CoulombF='Ising'")
            end
            push!(run_string, "HB2=True")




            
            Eimp_lev_tag = "Eimp=["
            if DMFT_spin_type > 0
                for s=1:DMFT_spin_type
                    for o in unique(diag(imp_ind[i]))
                        oi = findall(x->x == o, diag(imp_ind[i]))[1]
                        if ( (s != DMFT_spin_type) || ( diag(imp_ind[i])[oi] != maximum(diag(imp_ind[i])) ) )
                            Eimp_lev_tag = string(Eimp_lev_tag, regulated_Eimp[i][oi], ", ")
                        else
                            Eimp_lev_tag = string(Eimp_lev_tag, regulated_Eimp[i][oi], "]")
                        end
                    end
                end
            end

            J_tag = string("J=", imp_J[i])

            push!(run_string, Eimp_lev_tag)
            push!(run_string, J_tag)


            println(" Preparing the input file of EDFMTF : actqmc.cix ...")
            run(`$run_string`)
            GC.gc()

            mpi_prefix = string(mpi_prefix)
            excut_file ="ctqmc"
            input_file ="PARAMS"

            run_string = split(mpi_prefix)
            push!(run_string, excut_file)
            push!(run_string, input_file)


            run(`$run_string`)
            GC.gc()

            cd("../")

        end
    end

end








function write_paramsjson(cal_susc,Ineq_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv, imp_ind,imp_int_type, imp_dc,mu, H_loc, regulated_ShiftE, regulated_Eimp
    ,imp_Measure_time,imp_Thermal_time,beta, Solver_green_Cut, imp_U, imp_J, imp_int_parameterisation, F0, F2, F4, Trans_mat,green_basis, green_legendre_cutoff)
    
    print("Writing params.json ...\n")
    for i in Ineq_atom_Ind
        orb_ind = Corr_orbital_Ind[findall(x->x==i, Corr_atom_equiv)][1]
        sel_block_ind = [] 

        corr_ind = findall(x->x>0, imp_ind[i])
        for ori = 1:size(corr_ind)[1]
           push!(sel_block_ind, corr_ind[ori][1]) 
        end
        
        
        #transformation = Matrix{Float64}(I,  size(orb_ind)[1]*2,  size(sel_block_ind)[1]*2 )
        transformation = zeros(size(orb_ind)[1]*2, size(sel_block_ind)[1]*2)
    
        transformation[1:size(orb_ind)[1], 1:size(sel_block_ind)[1]] = Trans_mat[findall(x->x==i, Corr_atom_equiv)[1]][sel_block_ind,:]'
        transformation[size(orb_ind)[1]+1:size(orb_ind)[1]*2, size(sel_block_ind)[1]+1:size(sel_block_ind)[1]*2] = Trans_mat[findall(x->x==i, Corr_atom_equiv)[1]][sel_block_ind,:]'
        
        
        one_body = zeros(Float64, size(sel_block_ind)[1]*2, size(sel_block_ind)[1]*2)
        one_body[1:size(sel_block_ind)[1],1:size(sel_block_ind)[1]] = deepcopy(imp_ind[i][sel_block_ind,sel_block_ind])
        one_body[size(sel_block_ind)[1]+1:size(sel_block_ind)[1]*2,size(sel_block_ind)[1]+1:size(sel_block_ind)[1]*2] = deepcopy(imp_ind[i][sel_block_ind,sel_block_ind])
        for j in unique(imp_ind[i])
            if j>0
                mat_ind=findall(x->x==j,one_body)
                #orb_in_=findall(x->x==j,imp_ind[i])
                #one_body[mat_ind].= real(H_loc[orb_ind,orb_ind,1][orb_in_[1]])
                orb_in_=findall(x->x==j,diag(imp_ind[i]))
                one_body[mat_ind].= regulated_Eimp[i][orb_in_[1]] 
            end
        
        end
    
        # to Correct almost identical impurity levels to exactly the same
        #for ii= 1:size(sel_block_ind)[1]*2
        #    dia_ind=findall(x->abs.(x-one_body[ii,ii])<0.00003 && x != 0.0, one_body  )
        #    one_body[dia_ind].= one_body[ii,ii]
        #end
    
       
    
        imp_ind_block = zeros(Int64, size(sel_block_ind)[1]*2, size(sel_block_ind)[1]*2)
        imp_ind_block[1:size(sel_block_ind)[1],1:size(sel_block_ind)[1]] = deepcopy(imp_ind[i][sel_block_ind,sel_block_ind])
        if DMFT_spin_type == 1
            imp_ind_block[size(sel_block_ind)[1]+1:size(sel_block_ind)[1]*2,size(sel_block_ind)[1]+1:size(sel_block_ind)[1]*2] = deepcopy(imp_ind[i][sel_block_ind,sel_block_ind])
        elseif DMFT_spin_type ==2
           imp_ind_block[size(sel_block_ind)[1]+1:size(sel_block_ind)[1]*2,size(sel_block_ind)[1]+1:size(sel_block_ind)[1]*2] = maximum(imp_ind[i])*Matrix{Int64}(I,size(sel_block_ind)[1],size(sel_block_ind)[1])+ imp_ind[i][sel_block_ind,sel_block_ind]
        end
        imp_ind_str=string.(imp_ind_block)
        imp_ind_str[findall(x-> x =="0",imp_ind_str)].=""
        dir_str = string("imp_",i,"/")

    

        
        
        if imp_int_parameterisation[i] == "slater-condon"
            if (imp_int_type[i] == "ising" || imp_int_type[i] == "Ising" || imp_int_type[i] == "ISING")
                two_body = Dict("F0"=>F0[i], "F2"=>F2[i], "F4"=> F4[i], "approximation"=> "ising", "parametrisation"=> imp_int_parameterisation[i] )
                observables = Dict()
                probabilities = ["N", "energy", "Sz"]
                basis_type = "product"
            elseif (imp_int_type[i] == "full" || imp_int_type[i] == "Full" || imp_int_type[i] == "FULL")
                two_body = Dict("F0"=>F0[i], "F2"=>F2[i], "F4"=> F4[i], "approximation"=> "none", "parametrisation"=> imp_int_parameterisation[i] )
                observables = Dict("S2"=>Dict())
                probabilities = ["N", "energy","S2", "Sz"]
                basis_type = "product"
            end
        elseif imp_int_parameterisation[i] == "kanamori"
            if (imp_int_type[i] == "ising" || imp_int_type[i] == "Ising" || imp_int_type[i] == "ISING")
                two_body = Dict("U"=>imp_U[i], "J"=>imp_J[i], "approximation"=> "ising", "parametrisation"=> imp_int_parameterisation[i] )
                observables = Dict()
                probabilities = ["N", "energy", "Sz"]
                basis_type = ""
                transformation =[]
            elseif (imp_int_type[i] == "full" || imp_int_type[i] == "Full" || imp_int_type[i] == "FULL")
                two_body = Dict("U"=>imp_U[i], "J"=>imp_J[i], "approximation"=> "none", "parametrisation"=> imp_int_parameterisation[i] )
                observables = Dict("S2"=>Dict())
                probabilities = ["N", "energy","S2", "Sz"]
                basis_type = ""
                transformation =[]
            end            
        end    
        
        
            
        if cal_susc
                observables = Dict("S2"=>Dict())
                occupation_susceptibility_bulla = true
                occupation_susceptibility_direct = true
                probabilities = ["N", "energy","S2", "Sz"]
                quantum_number_susceptibility = true
                quantum_numbers = Dict("Sz"=>Dict())
                green_bulla = true

        elseif !cal_susc
                observables = Dict()
                occupation_susceptibility_bulla = false
                occupation_susceptibility_direct = false
                probabilities = ["N", "energy"]
                quantum_number_susceptibility = false
                quantum_numbers = Dict()
                green_bulla = false
        end
            
        
        
        
        if !( (green_basis == "legendre") || (green_basis == "Legendre") || (green_basis == "LEGENDRE") ) # matsubara case
            params =Dict(
                "basis"=> Dict("orbitals" => "d",   "transforamtion" => transformation,   "type" => basis_type)   
                ,"beta"=>beta
                ,"green basis"=>green_basis
                #,"green bulla"=> green_bulla
                #,"green matsubara cutoff"=> Solver_green_Cut
                #,"hloc" => Dict("one body"=>one_body, "quantum numbers" => Dict("N"=>Dict(), "Sz"=>Dict()), "two body"=> two_body  )
                ,"hloc" => Dict("one body"=>one_body, "two body"=> two_body  )
                ,"hybridisation" => Dict( "functions" => string(dir_str,"hyb.json"), "matrix" => imp_ind_str )
                ,"measurement time" => imp_Measure_time
                ,"mu" => -regulated_ShiftE[i]+mu #imp_dc[i]+mu
                #,"observables" => observables
                #,"occupation susceptibility bulla" => occupation_susceptibility_bulla
                ,"occupation susceptibility direct"=> occupation_susceptibility_direct
                #,"probabilities" => probabilities
                #,"quantum number susceptibility" => quantum_number_susceptibility
                #,"quantum numbers" => quantum_numbers
                #,"susceptibility cutoff" => 50
                #,"susceptibility tail" => 300
                ,"thermalisation time" => imp_Thermal_time
                ,"partition"=> Dict("green bulla"=> green_bulla, "density matrix precise"=> true, 
                    "green matsubara cutoff"=> Solver_green_Cut, "observables" => observables, 
                    "occupation susceptibility bulla" => occupation_susceptibility_bulla, 
                    "print density matrix" => true, "print eigenstates" => true, "probabilities" => probabilities,
                    "quantum number susceptibility" => quantum_number_susceptibility, "quantum numbers" => quantum_numbers,
                    "susceptibility cutoff" => 50, "susceptibility tail" => 300 )
    
            )            
            
        elseif ( (green_basis == "legendre") || (green_basis == "Legendre") || (green_basis == "LEGENDRE") )  # legendre case
            params =Dict(
                "basis"=> Dict("orbitals" => "d",   "transforamtion" => transformation,   "type" => "product")   
                ,"beta"=>beta
                ,"green basis"=>green_basis
                ,"green legendre cutoff"=>green_legendre_cutoff
                #,"green bulla"=> green_bulla
                #,"green matsubara cutoff"=> Solver_green_Cut
                #,"hloc" => Dict("one body"=>one_body, "quantum numbers" => Dict("N"=>Dict(), "Sz"=>Dict()), "two body"=> two_body  )
                ,"hloc" => Dict("one body"=>one_body, "two body"=> two_body  )
                ,"hybridisation" => Dict( "functions" => string(dir_str,"hyb.json"), "matrix" => imp_ind_str )
                ,"measurement time" => imp_Measure_time
                ,"mu" => -regulated_ShiftE[i]+mu #imp_dc[i]+mu
                #,"observables" => observables
                #,"occupation susceptibility bulla" => occupation_susceptibility_bulla
                ,"occupation susceptibility direct"=> occupation_susceptibility_direct
                #,"probabilities" => probabilities
                #,"quantum number susceptibility" => quantum_number_susceptibility
                #,"quantum numbers" => quantum_numbers
                #,"susceptibility cutoff" => 50
                #,"susceptibility tail" => 300
                ,"thermalisation time" => imp_Thermal_time
                ,"partition"=> Dict("green bulla"=> green_bulla, "density matrix precise"=> true, 
                    "green matsubara cutoff"=> Solver_green_Cut, "observables" => observables, 
                    "occupation susceptibility bulla" => occupation_susceptibility_bulla, 
                    "print density matrix" => true, "print eigenstates" => true, "probabilities" => probabilities,
                    "quantum number susceptibility" => quantum_number_susceptibility, "quantum numbers" => quantum_numbers,
                    "susceptibility cutoff" => 50, "susceptibility tail" => 300 )    
            )                  
            
            
        end
        
 
        #HERE paste

        if size(orb_ind) == 1
            params["basis"]["orbitals"] = "s"
        elseif size(orb_ind) == 3
            params["basis"]["orbitals"] = "p"
        elseif size(orb_ind) == 5
            params["basis"]["orbitals"] = "d"
        elseif size(orb_ind) == 7
            params["basis"]["orbitals"] = "f"
        end    
    
    
        if imp_int_parameterisation[i] == "kanamori"
            params["basis"]["orbitals"] = count(x->x>0, imp_ind[i] )
        end
        
    
        str =string("imp_",i)
        params_string = json(params, 4)
        open(string(str,"/","params.json"),"w") do f
            print(f, params_string)
        end
        
    end
    
    
end





function write_hybjson(Ineq_atom_Ind, hyb)
    for i in Ineq_atom_Ind
        str =string("imp_",i)
        json_string = json(hyb[string(i)],4)
        open(string(str,"/","hyb.json"),"w") do f
            print(f, json_string)
        end
    end
end


function make_imp_folder(Ineq_atom_Ind)
    for i in unique(abs.(Corr_atom_equiv))
        str =string("imp_",i)
        if !isdir(str)
            mkdir(str)
        end
    end
end



function hyb_array2dic(Ineq_atom_Ind, hyb_iWn, beta)
    
    print("Writing hyb.json ...\n")
    hyb= Dict{String,Dict{String,Dict{String,Any}}}()
    for i in Ineq_atom_Ind
        hyb[string(i)]=Dict{String,Dict{String,Any}}()
        for j =1:(size(hyb_iWn[i])[2])
            for s=1:size(hyb_iWn[i])[3]
                if s==1
                    hyb[string(i)][string(j)] = Dict("beta"=>beta, "real"=> real(hyb_iWn[i][:,j,s]), "imag"=> imag(hyb_iWn[i][:,j,s]) ) 
                elseif s==2
                    hyb[string(i)][string((size(hyb_iWn[i])[2]+j))] = Dict("beta"=>beta, "real"=> real(hyb_iWn[i][:,j,s]), "imag"=> imag(hyb_iWn[i][:,j,s]) ) 
                end
            end
        end
    end
    return hyb
end




function write_occlog_latt(DMFT_loop, Ineq_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, imp_dc, mu, Occup, DMFT_spin_type)
    open("occ.log", "a") do io
        @printf(io, "Iteration [%3i],    latt_mu : %10.6f,    imp_mu :", DMFT_loop, mu)
        for it in Ineq_atom_Ind
            @printf(io," %10.6f,", mu )#imp_dc[it]+mu  
        end
        
        if DMFT_spin_type == 1
            @printf(io, ",    Total_N : %10.6f\n", sum(real(diag(Occup[:,:,1])))*2 )
        elseif DMFT_spin_type ==2
            @printf(io, ",    Total_N : %10.6f\n", sum(real(diag(Occup[:,:,1]+Occup[:,:,2]))) )
        end
        
        
        @printf(io, "===========================================================================================\n",  )
        for it in Ineq_atom_Ind
            orb_ind = Corr_orbital_Ind[findall(x->x==it, Corr_atom_equiv)][1]
            @printf(io, "latt. %2i :", it )
            for or = 1:size(orb_ind)[1]
                if DMFT_spin_type == 1
                    if or != size(orb_ind)[1]
                        @printf(io," %10.6f,", diag(real(Occup[orb_ind,orb_ind,1]))[or] * 2.0)
                    else
                        @printf(io," %10.6f,  n : %10.6f\n", diag(real(Occup[orb_ind,orb_ind,1]))[or] *2.0, sum( diag(real(Occup[orb_ind,orb_ind,1])))*2.0  )
                    end
                elseif DMFT_spin_type == 2
                    if or != size(orb_ind)[1]
                        @printf(io," %10.6f,", diag(real(Occup[orb_ind,orb_ind,1]))[or]+ diag(real(Occup[orb_ind,orb_ind,2]))[or])
                    else
                        @printf(io," %10.6f,  [↑ : %5.3f, ↓ : %5.3f, m : %5.3f, n : %10.6f]\n", diag(real(Occup[orb_ind,orb_ind,1]))[or]+ diag(real(Occup[orb_ind,orb_ind,2]))[or],  sum( diag(real(Occup[orb_ind,orb_ind,1]))),  sum( diag(real(Occup[orb_ind,orb_ind,2]))),sum( diag(real(Occup[orb_ind,orb_ind,1])))- sum( diag(real(Occup[orb_ind,orb_ind,2]))) , sum( diag(real(Occup[orb_ind,orb_ind,1])))+ sum( diag(real(Occup[orb_ind,orb_ind,2])))  )   
                    end  
            
                end
            end

        end
        @printf(io, "-------------------------------------------------------------------------------------------\n",  )
    end
    
end



function write_scflog_latt(DMFT_loop, Ineq_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, imp_dc, mu, Occup, DMFT_spin_type,cal_susc,DMFT_solver)
    open("scf.log", "a") do io
        if DMFT_loop == 0
            @printf(io,"\n =======================")
            for i =1:maximum(Ineq_atom_Ind)*5
                if i!=maximum(Ineq_atom_Ind)*5
                    @printf(io, "============")
                else
                    @printf(io, "============\n")
                end
            end  
        
            @printf(io, "|Iter. [ n ]|    latt_mu|")
            
            
            for i in Ineq_atom_Ind
                @printf(io, "  latt%1i_occ|", i)
            end
            
            
            for i in Ineq_atom_Ind
                @printf(io, "    imp%1i_mu|", i)
            end
            
            for i in Ineq_atom_Ind
                @printf(io, "   imp%1i_occ|", i)
            end
        
            if ((cal_susc) || (DMFT_solver == "EDMFTF_ctqmc"))
                for i in Ineq_atom_Ind
                    @printf(io, "     imp%1i_N|", i)
                end
        
                for i in Ineq_atom_Ind
                    if i!=maximum(Ineq_atom_Ind)
                        @printf(io, "     imp%1i_S|", i)
                    else
                        @printf(io, "     imp%1i_S|\n", i)
                    end
                end
            else
                @printf(io, "\n")
            end
      
            @printf(io," =======================")
            for i =1:maximum(Ineq_atom_Ind)*5
                if i!=maximum(Ineq_atom_Ind)*5
                    @printf(io, "============")
                else
                    @printf(io, "============\n")
                end
            end        
        
        
   
        end
        @printf(io,"|Iter. [%3i]|%10.5f |",DMFT_loop, mu)
        
        for it in Ineq_atom_Ind
            orb_ind = Corr_orbital_Ind[findall(x->x==it, Corr_atom_equiv)][1]
            if DMFT_spin_type == 1
                @printf(io,"%10.5f |",  sum( diag(real(Occup[orb_ind,orb_ind,1])))*2.0  )
            elseif DMFT_spin_type == 2
                @printf(io,"%10.5f |", sum( diag(real(Occup[orb_ind,orb_ind,1])))+ sum( diag(real(Occup[orb_ind,orb_ind,2])))  )   
            end
        end
        
        
        
        for i in Ineq_atom_Ind
            @printf(io, "%10.5f |", mu ) # imp_dc[i]+mu)
        end
        
    
    end
    
end




function write_scflog_imp(Ineq_atom_Ind, imp_ind, DMFT_spin_type, Corr_ineq_orbital_Num, imp_occ, susc_N0, susc_S0, cal_susc, DMFT_solver)
    open("scf.log", "a") do io
        
        
        for i in Ineq_atom_Ind        
            occup_mat = zeros(Float64,size(imp_ind[i])[1],size(imp_ind[i])[1],DMFT_spin_type)
            
            if DMFT_spin_type == 1
                occup_mat[:,:,1] = deepcopy(imp_ind[i])
            elseif DMFT_spin_type == 2
                occup_mat[1:size(imp_ind[i])[1],1:size(imp_ind[i])[1],1] = deepcopy(imp_ind[i])
                occup_mat[1:size(imp_ind[i])[1],1:size(imp_ind[i])[1],2] = deepcopy(imp_ind[i])
            end
        
            occup_mat = deepcopy(convert(Array{Float64,3}, occup_mat))
        
            for or = 1:Corr_ineq_orbital_Num[i]
                for s = 1:DMFT_spin_type
                    o_ind = findall(x->x==convert(Float64,or),occup_mat[:,:,s])
                    if DMFT_solver == "ComCTQMC"
                        occup_mat[o_ind,s].= imp_occ[i][or,s]
                    elseif DMFT_solver == "EDMFTF_ctqmc"
                        occup_mat[o_ind,s].= imp_occ[i][or,s]./(size(o_ind)[1])
                    end                    
                end
            end
            
            if DMFT_spin_type == 1
                if DMFT_solver == "ComCTQMC"
                    @printf(io,"%10.5f |", sum( diag(occup_mat[:,:,1]) )*2.0  )
                elseif DMFT_solver == "EDMFTF_ctqmc"
                    @printf(io,"%10.5f |", sum( diag(occup_mat[:,:,1]) )  )
                end
            elseif DMFT_spin_type ==2
                @printf(io,"%10.5f |",sum( diag(occup_mat[:,:,1]) )+sum( diag(occup_mat[:,:,2]) )     )   
            end
        end
    
    
        if ( (cal_susc) || (DMFT_solver == "EDMFTF_ctqmc") )
            for i in Ineq_atom_Ind
                @printf(io, "%10.5f |", susc_N0[i])
            end
    
            for i in Ineq_atom_Ind
                if i!=Ineq_atom_Ind[end]
                    @printf(io, "%10.5f |", susc_S0[i])
                else
                    @printf(io, "%10.5f |\n", susc_S0[i])
                end
            end
        else
            @printf(io, "\n")
        end

    end
end









function Save_logfile(filename,Ineq_atom_Ind)
    if (filename == "params.obs" || filename == "hyb")
        for i in Ineq_atom_Ind
            save = true
            for j=1:500
                if !isfile(string("./imp_",i,"/",filename,"_",j,".json")) && save
                    save =false
                    folder=string("./imp_",i,"/")
                    file_1 = string(folder,filename,".json")
                    file_2 = string(folder,filename,"_",j-1,".json")
                    str=`cmp --quie $file_1 $file_2`
                    if !success(str)
                        cp(file_1,string("./imp_",i,"/",filename,"_",j,".json"))
                    end
                end
                
            end
        end    
    elseif (filename == "Sig" || filename == "Gf")
        for i in Ineq_atom_Ind
            save = true
            for j=1:500
                if !isfile(string("./imp_",i,"/",filename,"_",j,".out")) && save
                    save =false
                    folder=string("./imp_",i,"/")
                    file_1 = string(folder,filename,".out")
                    file_2 = string(folder,filename,"_",j-1,".out")
                    str=`cmp --quie $file_1 $file_2`
                    if !success(str)
                        cp(file_1,string("./imp_",i,"/",filename,"_",j,".out"))
                    end
                end
                
            end
        end   
    elseif (filename == "Delta")
        for i in Ineq_atom_Ind
            save = true
            for j=1:500
                if !isfile(string("./imp_",i,"/",filename,"_",j,".inp")) && save
                    save =false
                    folder=string("./imp_",i,"/")
                    file_1 = string(folder,filename,".inp")
                    file_2 = string(folder,filename,"_",j-1,".inp")
                    str=`cmp --quie $file_1 $file_2`
                    if !success(str)
                        cp(file_1,string("./imp_",i,"/",filename,"_",j,".inp"))
                    end
                end
                
            end
        end   
    elseif (filename == "ctqmc")
        for i in Ineq_atom_Ind
            save = true
            for j=1:500
                if !isfile(string("./imp_",i,"/",filename,"_",j,".log")) && save
                    save =false
                    folder=string("./imp_",i,"/")
                    file_1 = string(folder,filename,".log")
                    file_2 = string(folder,filename,"_",j-1,".log")
                    str=`cmp --quie $file_1 $file_2`
                    if !success(str)
                        cp(file_1,string("./imp_",i,"/",filename,"_",j,".log"))
                    end
                end
                
            end
        end   
    elseif (filename == "dc")
        for i in Ineq_atom_Ind
            save = true
            for j=1:500
                if !isfile(string("./imp_",i,"/",filename,"_",j,".dat")) && save
                    save =false
                    folder=string("./imp_",i,"/")
                    file_1 = string(folder,filename,".dat")
                    file_2 = string(folder,filename,"_",j-1,".dat")
                    str=`cmp --quie $file_1 $file_2`
                    if !success(str)
                        cp(file_1,string("./imp_",i,"/",filename,"_",j,".dat"))
                    end
                end
                
            end
        end            
    end       
    
end



function write_dc(Ineq_atom_Ind, imp_dc, imp_dc_type, imp_dc_n, imp_totN, dft_corrN)
    for i in Ineq_atom_Ind
        folder=string("./imp_",i,"/")
        filename = string("dc.dat")
        
        open(string(folder,filename), "w") do io
            if (imp_dc_type == "nominal" || imp_dc_type == "Norminal" || imp_dc_type == "NORMINAL")
                @printf(io,"occ : %10.5f\n", imp_dc_n[i]  )
                @printf(io,"dc (eV) : %10.5f", imp_dc[i]  )
            elseif (imp_dc_type == "fll" || imp_dc_type == "Fll" || imp_dc_type == "FLL")
                @printf(io,"occ : %10.5f\n", imp_totN[i]  )
                @printf(io,"dc (eV) : %10.5f", imp_dc[i]  )
            elseif (imp_dc_type == "FLL-DFT" || imp_dc_type == "fll-dft" || imp_dc_type == "Fll-dft")
                @printf(io,"occ : %10.5f\n", dft_corrN[i]  )
                @printf(io,"dc (eV) : %10.5f", imp_dc[i]  )
            end
        end
        
    end
end




function loadSelfE(Ineq_atom_Ind,DMFT_spin_type,SelfE_w,DMFT_solver)
    SelfE_w_dum = deepcopy(SelfE_w)
    
    if DMFT_solver == "ComCTQMC"
        for i in Ineq_atom_Ind
            print("Reading SelfE_w ... imp[",i,"]\n")
            str =string("./imp_",i,"/params.obs.json")
            obser=readJSON(str)
        
            ineq_orbital_Num = size(SelfE_w[i])[2]
            if DMFT_spin_type > 0
                for s = 1:DMFT_spin_type
                    for l=1:ineq_orbital_Num
                        SelfE_w_dum[i][:,l,s]=obser["partition"]["self-energy"][string( ineq_orbital_Num * (s-1) + l)]["function"]["real"].+obser["partition"]["self-energy"][string( ineq_orbital_Num * (s-1) + l)]["function"]["imag"]*im  # lattice selfE
                    end
                end
            end
        end
    elseif DMFT_solver == "EDMFTF_ctqmc"
        for i in Ineq_atom_Ind
            print("Reading SelfE_w ... imp[",i,"]\n")
            str = string("./imp_",i,"/Sig.out")
            open(str, "r") do f
                Nstring=readlines(f)
                linenum = size(Nstring)[1]
                
                ineq_orbital_Num = size(SelfE_w[i])[2]
                            
                if DMFT_spin_type > 0
                    for s = 1:DMFT_spin_type
                        for l=1:ineq_orbital_Num
                            for w = 1:(linenum-1)
                                SelfE_w_dum[i][w,l,s]= parse(Float64,split(Nstring[w+1])[2*(ineq_orbital_Num * (s-1) + l)]) + parse(Float64,split(Nstring[w+1])[2*(ineq_orbital_Num * (s-1) + l)+1])*im # lattice selfE
                            end
                        end
                    end
                end
                
            end
            
        end

    end
        
  
        
    return SelfE_w_dum
end




function loadSelfE_from_obs_cut(Ineq_atom_Ind,DMFT_spin_type,SelfE_w,SelfE_w_cut,iWnlist)
    SelfE_w_dum = deepcopy(SelfE_w)
    for i in Ineq_atom_Ind
        print("Reading SelfE_w ... imp[",i,"]\n")
        
        str =string("./imp_",i,"/params.obs.json")
        obser=readJSON(str)
        
        ineq_orbital_Num = size(SelfE_w[i])[2]
        if DMFT_spin_type > 0
           for s = 1:DMFT_spin_type
                for l=1:ineq_orbital_Num
                    SelfE_w_dum[i][:,l,s]=obser["partition"]["self-energy"][string( ineq_orbital_Num * (s-1) + l)]["function"]["real"].+obser["partition"]["self-energy"][string( ineq_orbital_Num * (s-1) + l)]["function"]["imag"]*im  # lattice selfE
                    SelfE_w_cut[i][1:size(iWnlist)[1],l,s]=SelfE_w_dum[i][1:size(iWnlist)[1],l,s]
                end
            end
        end

    end
    return SelfE_w_cut
end



function Check_restart_wSelfE(Ineq_atom_Ind, DMFT_solver)
    true_table=[]
    if (DMFT_solver == "ComCTQMC")
        for i in Ineq_atom_Ind
            if isfile(string("./imp_",i,"/params.obs.json"))
                push!(true_table,true)
            else
                push!(true_table,false)
            end
        end
    elseif (DMFT_solver == "EDMFTF_ctqmc")
        for i in Ineq_atom_Ind
            if isfile(string("./imp_",i,"/Sig.out"))
                push!(true_table,true)
            else
                push!(true_table,false)
            end
        end
    end
    
    if size(findall(true_table))[1]==size(Ineq_atom_Ind)[1]
        restart_wSelfE = true
    else
        restart_wSelfE = false
    end
    return restart_wSelfE
end



function Print_IntOccup(DMFT_spin_type,Occup)
    targetN =0.0
    if DMFT_spin_type == 1
        println("Calculated Occup : ",sum(diag(Occup[:,:,1]))+sum(diag(Occup[:,:,1])))
        targetN = round(real(sum(diag(Occup[:,:,1]))+sum(diag(Occup[:,:,1]))))  # This value is not intger, but if we change it as intger, it makes real(hyb (w->inf)) .ne. 0
        println("From now on system charge # : ",targetN )
        println("")
    else
        println("Calculated Occup : ",sum(diag(Occup[:,:,1]))+sum(diag(Occup[:,:,2])))
        targetN = round(real(sum(diag(Occup[:,:,1]))+sum(diag(Occup[:,:,2]))))  # This value is not intger, but if we change it as intger, it makes real(hyb (w->inf)) .ne. 0
        println("From now on system charge # : ",targetN )
        println("")
    end
    return targetN
end





function show_input(K_grid_num, R_grid_num, iW_grid_cut, Tem, DMFT_spin_type, Corr_atom_Ind, Corr_atom_equiv, Solver_green_Cut, basis_transform, consider_degeneracy
  ,green_basis,init_bias,Mix_selfE ,imp_dc_type ,imp_Measure_time, imp_Thermal_time, imp_U, imp_J, imp_dc_n, imp_dc, imp_int_type, imp_int_parameterisation, imp_ind, smth_step, cal_susc)
    
    println("")
    println("")
    println("")
    println("")
    println("")
    println("                 =============================================  ")
    println("                 ================    Jx DMFT  ================  ")
    println("                 =============================================  ")
    println("\n\n\n------------- input -------------")
    println("K-points grid : ", K_grid_num)
    println("R-points grid : ", R_grid_num)
    println("iW_grid_cut : ", iW_grid_cut)
    println("Temperature : ", Tem)
    println("DMFT spin type : ", DMFT_spin_type)
    println("Correlated atoms : ", Corr_atom_Ind)
    println("Correlated atoms equiv. : ", Corr_atom_equiv)
    println("")
    println("Green-function basis : ", green_basis)
    println("Solver green cut w: ", Solver_green_Cut)
    println("Solver Measure time: ", imp_Measure_time)
    println("Solver thermalization time: ", imp_Thermal_time)
    println("First smth. step for self energy: ", smth_step)
    println("Cal. susc. ?: ", cal_susc)
    println("initial bias of selfE for spin-polarized : +/-", init_bias)
    println("self-energy simple mixing ratio of old selfE : ", Mix_selfE)
    println("")
    println("transform basis ?: ", basis_transform)
    println("consider degeneracy ?: ", consider_degeneracy)   
    println("")
    println("double counting type: ", imp_dc_type)
    
    
    for i in unique(abs.(Corr_atom_equiv))
        println("")
        println("                imp [",i,"]")
        println("U     : ", imp_U[i])
        println("J     : ", imp_J[i])
        println("interaction type  :", imp_int_type[i])
        println("interaction parameterisation  :", imp_int_parameterisation[i])
        println("dc_n0 : ", imp_dc_n[i],",   [ nominal DC (eV) : ", imp_dc[i], " ]")
        println("imp-",i," index :")
        for j =1:size(imp_ind[i])[1]
            println("            ",imp_ind[i][j,:])
        end

        println("")

    end
    println("------------- input -------------")
    println("")

end


function init_selfE_dc(imp_dc,Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,iWnlist)

    SelfE_w=[]
    #print("===== SelfE = 0.0+0.0im block initialize ... ===== \n")
    if string(hamiltonian_info.dfttype) == "Wannier90"
        for i in Ineq_atom_Ind
            if DMFT_spin_type >0
                SelfE_block= zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],DMFT_spin_type)+zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],DMFT_spin_type)*im;
                SelfE_block .= imp_dc[i]
                push!(SelfE_w,SelfE_block)
            elseif(DMFT_spin_type==0)
                if string(hamiltonian_info.dfttype) == "OpenMX" && hamiltonian_info.scf_r.SpinP_switch == 1
                    SelfE_block= zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],hamiltonian_info.scf_r.SpinP_switch+1)+zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],hamiltonian_info.scf_r.SpinP_switch+1)*im;
                elseif string(hamiltonian_info.dfttype) == "Wannier90"
                    SelfE_block= zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],hamiltonian_info.scf_r.SpinP_switch)+zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],hamiltonian_info.scf_r.SpinP_switch)*im;                
                end
                push!(SelfE_w,SelfE_block)
            end
        end
    println("")
    elseif string(hamiltonian_info.dfttype) == "OpenMX"
        for i in Ineq_atom_Ind
            j = findall(x->x ==i, Corr_atom_equiv)[1]
            orb_size =size(Corr_orbital_Ind[j])[1]
            SelfE_block = zeros(size(iWnlist)[1],orb_size,2)+zeros(size(iWnlist)[1],orb_size,2)*im
            push!(SelfE_w,SelfE_block)
        end
    end
    
    
    return SelfE_w
end



function init_selfE_zero(Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,iWnlist)

    SelfE_w=[]
    #print("===== SelfE = 0.0+0.0im block initialize ... ===== \n")
    if string(hamiltonian_info.dfttype) == "Wannier90"
        for i in Ineq_atom_Ind
            if DMFT_spin_type >0
                SelfE_block= zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],DMFT_spin_type)+zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],DMFT_spin_type)*im;
                push!(SelfE_w,SelfE_block)
            elseif(DMFT_spin_type==0)
                if string(hamiltonian_info.dfttype) == "OpenMX" && hamiltonian_info.scf_r.SpinP_switch == 1
                    SelfE_block= zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],hamiltonian_info.scf_r.SpinP_switch+1)+zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],hamiltonian_info.scf_r.SpinP_switch+1)*im;
                elseif string(hamiltonian_info.dfttype) == "Wannier90"
                    SelfE_block= zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],hamiltonian_info.scf_r.SpinP_switch)+zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],hamiltonian_info.scf_r.SpinP_switch)*im;                
                end
                push!(SelfE_w,SelfE_block)
            end
        end
    println("")
    elseif string(hamiltonian_info.dfttype) == "OpenMX"
        for i in Ineq_atom_Ind
            j = findall(x->x ==i, Corr_atom_equiv)[1]
            orb_size =size(Corr_orbital_Ind[j])[1]
            SelfE_block = zeros(size(iWnlist)[1],orb_size,2)+zeros(size(iWnlist)[1],orb_size,2)*im
            push!(SelfE_w,SelfE_block)
        end
    end
    
    
    return SelfE_w
end



function init_selfE_dcbias(imp_dc,Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,iWnlist,init_bias)

    SelfE_w=[]
    print("===== SelfE = imp_dc block initialize ... ===== \n")
    for i in Ineq_atom_Ind
        if DMFT_spin_type >0
            print("Impurity [",i,"] : orbital # :", Corr_ineq_orbital_Num[i], ",    Spin type :", DMFT_spin_type,"\n" )
            SelfE_block= zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],DMFT_spin_type)+zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],DMFT_spin_type)*im;
            SelfE_block .= imp_dc[i]
            if DMFT_spin_type == 2
                SelfE_block[:,:,1] .= imp_dc[i]-init_bias
                SelfE_block[:,:,2] .= imp_dc[i]+init_bias  # initial splitting condition
            end
            push!(SelfE_w,SelfE_block)
    elseif(DMFT_spin_type==0)
            print("Impurity [",i,"] : orbital # :", Corr_ineq_orbital_Num[i], ",    Spin type :", hamiltonian_info.scf_r.SpinP_switch,"\n" )
            SelfE_block= zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],hamiltonian_info.scf_r.SpinP_switch)+zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],hamiltonian_info.scf_r.SpinP_switch)*im;
            push!(SelfE_w,SelfE_block)
        end
    end
    println("")


    return SelfE_w
end


function init_selfE_bias(Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,iWnlist,init_bias)

    SelfE_w=[]
    print("===== SelfE = 0.0+0.0im block initialize ... ===== \n")
    for i in Ineq_atom_Ind
        if DMFT_spin_type >0
            print("Impurity [",i,"] : orbital # :", Corr_ineq_orbital_Num[i], ",    Spin type :", DMFT_spin_type,"\n" )
            SelfE_block= zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],DMFT_spin_type)+zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],DMFT_spin_type)*im;
            if DMFT_spin_type == 2
                SelfE_block[:,:,1] .= -init_bias
                SelfE_block[:,:,2] .= +init_bias  # initial splitting condition
            end
            push!(SelfE_w,SelfE_block)
    elseif(DMFT_spin_type==0)
            print("Impurity [",i,"] : orbital # :", Corr_ineq_orbital_Num[i], ",    Spin type :", hamiltonian_info.scf_r.SpinP_switch,"\n" )
            SelfE_block= zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],hamiltonian_info.scf_r.SpinP_switch)+zeros(size(iWnlist)[1],Corr_ineq_orbital_Num[i],hamiltonian_info.scf_r.SpinP_switch)*im;
            push!(SelfE_w,SelfE_block)
        end
    end
    println("")


    return SelfE_w
end








function remove_spikes(Corr_atom_equiv,Tem,iWnlist,SelfE_w;spike_ratio_cut=0.03)
    esp = 1e-5
    SelfE_w_smooth = []
    tem_energy = Tem/11604
    remove_mask_arr = []

    for i in unique(abs.(Corr_atom_equiv))
        SelfE_w_smooth_block = zeros(ComplexF64, size(SelfE_w[i]))

        for l=1:size(SelfE_w_smooth_block)[2]
            for s =1:size(SelfE_w_smooth_block)[3]
                
                SelfE_w_smooth_block_dum = zeros(ComplexF64, size(SelfE_w_smooth_block)[1])
                
                wsec1=1:3
                ker = ImageFiltering.Kernel.gaussian((0.1,))
                SelfE_w_smooth_block_dum[wsec1] = imfilter(real(SelfE_w[i][wsec1,l,s]), ker).+imfilter(imag(SelfE_w[i][wsec1,l,s]), ker).*im
           
                wsec2=(wsec1[end]+1):size(iWnlist)[1]
                #wsec2=findall(x->x>= (tem_energy*30), imag(iWnlist[:]))
                ker = ImageFiltering.Kernel.gaussian((2.0,))
                SelfE_w_smooth_block_dum[wsec2] = imfilter(real(SelfE_w[i][wsec2,l,s]), ker).+imfilter(imag(SelfE_w[i][wsec2,l,s]), ker).*im


                remove_mask = zeros(Bool,size(SelfE_w_smooth_block_dum)[1])
                remove_mask = remove_mask .| (0. .< imag(SelfE_w[i][:,l,s])) # if Self E is positive
                
                spike_ratio = (abs.(imag(SelfE_w[i][:,l,s])) .-abs.(imag(SelfE_w_smooth_block_dum)) .+ esp) ./ (abs.(imag(SelfE_w_smooth_block_dum)) .+ esp);
                    
                remove_mask = remove_mask .| (spike_ratio_cut .< abs.(spike_ratio)) # if spike is to large
                push!(remove_mask_arr,sum(remove_mask))
                remove_mask[wsec1] .= false
                keep_mask = .! remove_mask;
                
                
                
                sig_real_fit = Spline1D(imag(iWnlist[keep_mask]), real(SelfE_w[i][keep_mask,l,s]); k=2,s=0.01)
                sig_imag_fit = Spline1D(imag(iWnlist[keep_mask]), imag(SelfE_w[i][keep_mask,l,s]); k=2,s=0.01)

                SelfE_w_smooth_block[:,l,s] = sig_real_fit(imag(iWnlist[:])) .+ sig_imag_fit(imag(iWnlist[:])).*im
                

                SelfE_w_smooth_block[wsec1,l,s] = SelfE_w[i][wsec1,l,s]                 
                
                dupoint = [wsec2[1]-3,wsec2[1]-2,wsec2[1]-1,wsec2[1]+1,wsec2[1]+2,wsec2[1]+3,wsec2[1]+4,wsec2[1]+5,wsec2[1]+6]
                sig_real_fit_2 = Spline1D(imag(iWnlist[dupoint]), real(SelfE_w_smooth_block[dupoint,l,s]); k=2, s=0.0)
                sig_imag_fit_2 = Spline1D(imag(iWnlist[dupoint]), imag(SelfE_w_smooth_block[dupoint,l,s]); k=2, s=0.0)
                
                SelfE_w_smooth_block[dupoint,l,s] = sig_real_fit_2(imag(iWnlist[dupoint])).+sig_imag_fit_2(imag(iWnlist[dupoint])).*im
                SelfE_w_smooth_block[wsec2[1],l,s] = sig_real_fit_2(imag(iWnlist[wsec2[1]]))+sig_imag_fit_2(imag(iWnlist[wsec2[1]]))*im

                SelfE_w_smooth_block[wsec1,l,s] = SelfE_w[i][wsec1,l,s] 



            end
        end
        push!(SelfE_w_smooth,SelfE_w_smooth_block)
        
    end
    return SelfE_w_smooth, remove_mask_arr
end



#function SelfE_smoothing(Corr_atom_equiv,Tem,iWnlist,SelfE_w)
#    SelfE_w_sm, remove_mask_arr = remove_spikes(Corr_atom_equiv,Tem,iWnlist,SelfE_w;spike_ratio_cut=0.5) # save smoothing version of selfE
#    for i = 1:10
#        if (sum(remove_mask_arr) > 0 )|| (i < 3)
#            SelfE_w_sm, remove_mask_arr = remove_spikes(Corr_atom_equiv,Tem,iWnlist,SelfE_w_sm;spike_ratio_cut=0.5)
#            println("Removed spike #:", convert(Array{Int64,1},remove_mask_arr[1:end]))
#        end
#    end
#    return SelfE_w_sm
#end


function SelfE_smoothing(Corr_atom_equiv,Tem,iWnlist,SelfE_w)
    SelfE_w_sm, remove_mask_arr = remove_spikes(Corr_atom_equiv,Tem,iWnlist,SelfE_w;spike_ratio_cut=0.5) # save smoothing version of selfE        
    println("Removed spike #:", convert(Array{Int64,1},remove_mask_arr[1:end]))
    #SelfE_w_sm, remove_mask_arr = remove_spikes(Corr_atom_equiv,Tem,iWnlist,SelfE_w_sm,spike_ratio_cut=0.5)
    #println("Removed spike #:", convert(Array{Int64,1},remove_mask_arr[1:end]))

    #for i = 1:2
    #    if (sum(remove_mask_arr) > 0 )|| (i < 2)
    #        #SelfE_w_sm, remove_mask_arr = remove_spikes(Corr_atom_equiv,Tem,iWnlist,SelfE_w_sm;spike_ratio_cut=0.5)
    #        SelfE_w_sm, remove_mask_arr = remove_spikes(Corr_atom_equiv,Tem,iWnlist,SelfE_w_sm,spike_ratio_cut=0.5)
    #        println("Removed spike #:", convert(Array{Int64,1},remove_mask_arr[1:end]))
    #    end
    #end

    #println("Removed spike #:", convert(Array{Int64,1},remove_mask_arr[1:end]))

    return SelfE_w_sm
end



function write_outfile(str,iWnlist,funct,Corr_atom_equiv,DMFT_spin_type,imp_ind)

    if str == "SelfE"
        filename="sig_bare"
    elseif str == "SelfE_sm"
        filename="sig_smth"
    elseif str == "SelfE_mix"
        filename="sig_mix"
    elseif str == "Delta"
        filename="delta"
    elseif str == "GreenF"
        filename="Gloc"
    end
    
    open(string(filename,".dat"), "w") do io
        strRe=""
        strIm=""
        iwn_str="Im{iWn}"
        @printf(io," %20s",iwn_str)
        if DMFT_spin_type>0
            for i in unique(abs.(Corr_atom_equiv))
                for s =1:DMFT_spin_type
                    for l in unique(imp_ind[i])
                        if l>0
                            strRe = string("Re{",str,"_",i,"}[",l,",",s,"]")
                            strIm = string("Im{",str,"_",i,"}[",l,",",s,"]")   
                            @printf(io, " %20s %20s", strRe, strIm)
                        end
                    end
                end
            end
        end
        
        @printf(io,"\n")

        for w = 1:size(iWnlist)[1]
            @printf(io, " %20.14f", imag(iWnlist[w]))
            if DMFT_spin_type>0
                for i in unique(abs.(Corr_atom_equiv))
                    for s =1:DMFT_spin_type
                        for l in unique(imp_ind[i])
                            if l>0
                                @printf(io, " %20.14f %20.14f", real(funct[i][w,l,s]), imag(funct[i][w,l,s]))
                            end
                        end
                    end
                end
            end
            @printf(io,"\n")
        
        end
    
        
    end;

        
    save = true   
    for j=1:500
        if !isfile(string(filename,"_",j,".dat")) && save
            save =false
            file_1 = string(filename,".dat")
            file_2 = string(filename,"_",j-1,".dat")
            str=`cmp --quie $file_1 $file_2`
            if !success(str)
                cp(file_1,string(filename,"_",j,".dat"))
            end
                
        end
    end

end



function write_outfile_mat(str,iWnlist,funct)

    if str == "Gfmat"
        filename="Gloc_mat"
    elseif str == "Dmat"
        filename="delta_mat"
    end
    
    for i = 1:size(funct)[1]
        
        open(string(filename,i,".dat"), "w") do io
            strRe=""
            strIm=""
            iwn_str="Im{iWn}"
            @printf(io," %20s",iwn_str)
            if DMFT_spin_type>0              
                for s =1:DMFT_spin_type
                    for l1 = 1:size(funct[i])[2]
                        for l2 = 1:size(funct[i])[3]
                            strRe = string("Re{",str,"_",l1,"_",l2,"}[",s,"]")
                            strIm = string("Im{",str,"_",l1,"_",l2,"}[",s,"]")
                            @printf(io, " %20s %20s", strRe, strIm)
                        end
                    end
                end
            end
            @printf(io,"\n")
            for w = 1:size(iWnlist)[1]
                @printf(io, " %20.14f", imag(iWnlist[w]))
                if DMFT_spin_type>0
                    for s =1:DMFT_spin_type
                        for l1 = 1:size(funct[i])[2]
                            for l2 = 1:size(funct[i])[3]
                                @printf(io, " %20.14f %20.14f", real(funct[i][w,l1,l2,s]), imag(funct[i][w,l1,l2,s]))

                            end
                        end
                    end
                end
                @printf(io,"\n")
            end
        end
    end

end




function EDMFTF_write_Deltainp(iWnlist,funct,Corr_atom_equiv,DMFT_spin_type,imp_ind)

    for i in unique(abs.(Corr_atom_equiv))
        open(string("imp_",i,"/Delta.inp"), "w") do io
            for w = 1:size(iWnlist)[1]
                @printf(io, " %20.14f", imag(iWnlist[w]))
                if DMFT_spin_type>0
                    for s =1:DMFT_spin_type
                        for l in unique(imp_ind[i])
                            if l>0
                                @printf(io, " %20.14f %20.14f", real(funct[i][w,l,s]), imag(funct[i][w,l,s]))
                            end
                        end
                    end
                end
                @printf(io,"\n")
            end
            @printf(io,"\n")
        end 
    end;


end





function readJSON(jsonfile)
    ob=Dict()
    open(jsonfile, "r") do f
        a=read(f,String)
        ob=JSON.parse(a)
    end
    
    return ob
end


function duplicatedfnt_atomlist(atom12_list)
    atomlist_primitive=[]
    for i=1:size(atom12_list)[1]
        for j=1:length(atom12_list[i])
            if !(atom12_list[i][j] in atomlist_primitive)
                push!(atomlist_primitive, atom12_list[i][j])
            end
        end
    end
    return atomlist_primitive
end


function atom_orb_list(hamiltonian_info)
    orbitalNums =  hamiltonian_info.scf_r.Total_NumOrbs
    orbitalStartIdx_list = cumsum(orbitalNums) .- orbitalNums[1]
    atom_orbitals=[]
    for i=1:hamiltonian_info.scf_r.atomnum
        push!(atom_orbitals,orbitalStartIdx_list[i] .+ (1:orbitalNums[i]))
    end
    return atom_orbitals
    
end

#### Grid section ####


function kPoint_gen(k_point_num)
    k_point_list=zeros(k_point_num[1],k_point_num[2],k_point_num[3],3)
    for k1 in ( 1 : k_point_num[1]  )
        for k2 in ( 1 : k_point_num[2]  )
            for k3 in ( 1 : k_point_num[3]  )
                k_point_list[k1,k2,k3,:]=[(k1-1)/k_point_num[1],(k2-1)/k_point_num[2],(k3-1)/k_point_num[3]]  
            end
        end
    end
    return k_point_list
end




function RPoint_gen(R_point_num)
    R_point_list=zeros(R_point_num[1],R_point_num[2],R_point_num[3],3)
    for R1 in ( 1 : R_point_num[1]  )
        for R2 in ( 1 : R_point_num[2]  )
            for R3 in ( 1 : R_point_num[3]  )
                R_point_list[R1,R2,R3,:]=[R1-1,R2-1,R3-1]
            end
        end
    end
    return R_point_list
end

function Matsubara_gen(iWnCut,beta)
    Matsubara_f_num=convert(Int64,floor(((iWnCut*beta/pi)+1)/2))
    iWnlist=zeros(Matsubara_f_num)*im
    for i =1:(Matsubara_f_num)
        iWnlist[i]=((2.0*(i-1)+1)*pi/beta)*im   
    end
    return iWnlist
end

#### real frequency grid calculation




function trapezoidalInt(wlist, Aw_ar, range)
    Aw = 0.0
    for i =1:size(wlist)[1]
        if (range[1]<wlist[i]<range[2])
            Aw += (wlist[i+1]-wlist[i])/2.0*(Aw_ar[i]+Aw_ar[i+1])
        end
    end
    return Aw
end




@everywhere function init_variables_grid_DMFT_realf(wlist,Klist)
    global g_wlist,g_Klist
   
    #g_SelfE_w = SelfE_w
    g_wlist = wlist;
    g_Klist = Klist;

end


function wlist_gen(wmax,wNum)
    
    if mod(wNum,2) == 0
    
        halfwNum = convert(Int64,wNum/2.0)
        wlist_h=zeros(wNum)
        r=(wmax)^(1/4)/(halfwNum)
        for i=1:halfwNum
                wlist_h[halfwNum+i]=(r*(i-2))^(4.0)
                wlist_h[halfwNum-i+1]=-(r*(i-2))^(4.0)
        end
        
        
    else
        halfwNum = convert(Int64,(wNum+1)/2.0)
        wlist_h=zeros(wNum)
        r=(wmax)^(1/4)/(halfwNum)
        for i=1:(halfwNum-1)
                wlist_h[halfwNum+i]=(r*(i))^(4.0)
                wlist_h[halfwNum-i]=-(r*(i))^(4.0)
        end        
        
        
    end
    
    return wlist_h
end




@everywhere function dcMat()
    dc_Mat = zeros(size(g_H_k[1,1,1,:,:,1]))+zeros(size(g_H_k[1,1,1,:,:,1]))*im
    for j = 1:size(g_Corr_atom_Ind)[1]
        
        block_dc = zeros(size(g_imp_ind[abs(g_Corr_atom_equiv[j])])[1],size(g_imp_ind[abs(g_Corr_atom_equiv[j])])[1])+zeros(size(g_imp_ind[abs(g_Corr_atom_equiv[j])])[1],size(g_imp_ind[abs(g_Corr_atom_equiv[j])])[1])*im
        for l =1:g_Corr_ineq_orbital_Num[abs(g_Corr_atom_equiv[j])]
            if g_Corr_atom_equiv[j] >0
                block_dc[findall(x->x==l,  g_imp_ind[abs(g_Corr_atom_equiv[j])] )].= g_imp_dc[abs(g_Corr_atom_equiv[j])]
            else
                block_dc[findall(x->x==l,  g_imp_ind[abs(g_Corr_atom_equiv[j])] )].= g_imp_dc[abs(g_Corr_atom_equiv[j])]
            end
        end
        
        dc_Mat[g_Corr_orbital_Ind[j],g_Corr_orbital_Ind[j]]= block_dc;
        
    end
    return dc_Mat
end




function dcMat_noglob(H_k,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,imp_ind)
    dc_Mat = zeros(size(H_k[1,1,1,:,:,1]))+zeros(size(H_k[1,1,1,:,:,1]))*im
    for j = 1:size(Corr_atom_Ind)[1]
        
        block_dc = zeros(size(imp_ind[abs(Corr_atom_equiv[j])])[1],size(imp_ind[abs(Corr_atom_equiv[j])])[1])+zeros(size(imp_ind[abs(Corr_atom_equiv[j])])[1],size(imp_ind[abs(Corr_atom_equiv[j])])[1])*im
        for l =1:Corr_ineq_orbital_Num[abs(Corr_atom_equiv[j])]
            if Corr_atom_equiv[j] >0
                block_dc[findall(x->x==l,  imp_ind[abs(Corr_atom_equiv[j])] )].= imp_dc[abs(Corr_atom_equiv[j])]
            else
                block_dc[findall(x->x==l,  imp_ind[abs(Corr_atom_equiv[j])] )].= imp_dc[abs(Corr_atom_equiv[j])]
            end
        end
        
        dc_Mat[Corr_orbital_Ind[j],Corr_orbital_Ind[j]]= block_dc;
        
    end
    return dc_Mat
end




function local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind)
       
    H0=zeros(ComplexF64,size(H_k)[4],size(H_k)[5])
        
    for k1 =1:size(H_k)[1]
        for k2 =1:size(H_k)[2]
            for k3=1:size(H_k)[3]
                H0 += H_k[k1,k2,k3,:,:,1]/(size(H_k)[1]*size(H_k)[2]*size(H_k)[3])
            end
        end
    end

    for i=1:size(imp_ind)[1]
        orb_ind = findall(x->x==i,Corr_atom_equiv)[1]
        println("----------------------------------------------------")
        println("Inequv. orbital index of imp. ",i," :", Corr_orbital_Ind[orb_ind])
        println("Diagonal part of Hloc : ",real(diag(H0[Corr_orbital_Ind[orb_ind],Corr_orbital_Ind[orb_ind]])))
        println("\n")
    end
    
    return H0
end




@everywhere function GreenF_at_agrid_realf(w,Kind,spin,SelfE_agrid,testmu,delta)
      
    Gmat = zeros(size(g_H_k[1,1,1,:,:,1]))+zeros(size(g_H_k[1,1,1,:,:,1]))*im 
    Gmat =inv( ( (w+testmu+delta*im)*Matrix{Float64}(I,size(g_H_k[1,1,1,:,:,1]) )) - g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin] - SelfE_agrid )

    return Gmat
end


@everywhere function G_k_sum_realf(args)
    w = args[1]
    spin = args[2]
    testmu = args[3]
    delta = args[4]
    
    wnum =convert(Int64, size(g_wlist)[1]/5)
    if (mod(w,wnum) == 0)
        per = convert(Int64,w/wnum)
        println("Calculating G(w) ...   ",per*20, "%  done")
    end

        
    G_w_s=zeros(g_Total_orb_Num,g_Total_orb_Num)+zeros(g_Total_orb_Num,g_Total_orb_Num)*im
    for k1 = 1:size(g_Klist)[1]
        for k2 = 1:size(g_Klist)[2]
            for k3 = 1:size(g_Klist)[3]
                G_w_s += GreenF_at_agrid_realf(g_wlist[w],[k1 k2 k3], spin, SelfE_at_agrid(w,spin),testmu,delta)
            end
        end
    end
    return G_w_s
end

function Gloc_gen_realf(H_k,wlist,SelfE_w,testmu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind,delta)
    
    
    println("===== Start to calculate wannier occupancy in real frequency ... ======")
    G_loc_w = zeros(size(wlist)[1],size(H_k)[4],size(H_k)[5],size(H_k)[6])+zeros(size(wlist)[1],size(H_k)[4],size(H_k)[5],size(H_k)[6])*im
    arg_list = [];
    for w = 1:size(wlist)[1]
        for s = 1:size(H_k)[6]
            push!(arg_list,( w,s, testmu, delta) )
        end
    end
    println("Passed all variables to calculate G(w) ... ")

    G_w_s_list = pmap(G_k_sum_realf, arg_list)
    i=0
    for w = 1:size(wlist)[1]
        for s = 1:size(H_k)[6]
            i+=1
            G_loc_w[w,:,:,s]=G_w_s_list[i]/( size(H_k)[1]*size(H_k)[2]*size(H_k)[3] )

        end
    end

    return  G_loc_w
end






function Aw_orb(Corr_orbital_Ind,Corr_atom_equiv,wlist, G_loc_w)

    
    dft_corrN = []
    Aw = -imag(deepcopy(G_loc_w))./pi
    println("")
    println("")
    println("For non-mag cal, total d-orbital occupation 5 !")
    println("===============================================")
    open("DFT_CorrNele.dat", "w") do io
        for i in unique(Corr_atom_equiv)
            j = findall(x->x==i, Corr_atom_equiv)[1]
            sumDOS = zeros(size(G_loc_w)[1])
            for o in Corr_orbital_Ind[j]
                println("orb-",o,", occup : ", trapezoidalInt(wlist,  Aw[:,o,o,1], [wlist[1],0.0]) )
                sumDOS += Aw[:,o,o,1]
            end
            wan_corr_occ =trapezoidalInt(wlist,  sumDOS, [wlist[1],0.0])*2.0
            println("--------------------------------------")
            println("latt-",i,", occup(up+dn) : ", wan_corr_occ )
            println("======================================")
            
            @printf(io, "latt. %2i :", i )
            @printf(io," %10.6f\n", wan_corr_occ)
            push!(dft_corrN, wan_corr_occ)
            
        end
    end
    sumDOS = zeros(size(G_loc_w)[1])
    for o = 1:size(Aw)[2]
        sumDOS += Aw[:,o,o,1]
    end
    
    wan_occ = trapezoidalInt(wlist,  sumDOS, [wlist[1],0.0])*2.0
    println("")    
    println("*------------------------------------*")    
    println("Wannier-Tot-occup(up+dn):" , wan_occ )
    println("*------------------------------------*")    


    

    return Aw, round(wan_occ), dft_corrN
end






@everywhere function Cal_InvWiess_fromGloc_0(SelfE_realw0,w0ind,G_loc_w0,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,DMFT_spin_type)
    InvWiess_w0=[]
    arg_list=[];
    for i =1:size(Corr_atom_equiv)[1]
        o_ind = Corr_orbital_Ind[i]
        Gloc_inv_block = inv(deepcopy(G_loc_w0[o_ind,o_ind]))
        InvWiess_block = Gloc_inv_block + SelfE_realw0[o_ind,o_ind]
        

        push!(InvWiess_w0,InvWiess_block)
    end
    
    return InvWiess_w0
end








function Cal_hyb_fromInvWiess_0(InvWiess_w0,imp_ind,w0,mu,H0,DMFT_spin_type,Corr_atom_equiv,Corr_orbital_Ind,delta)
    hyb_w0=[]
    for i =1:size(InvWiess_w0)[1]
        hyb_iWn_Block = zeros(ComplexF64,size(InvWiess_w0[i]) )
        mat_ind = Corr_orbital_Ind[i]
        IdenMat= Matrix{ComplexF64}(I,size(mat_ind)[1],size(mat_ind)[1])
        hyb_iWn_Block=(w0+mu+delta*im)*IdenMat-InvWiess_w0[i]-H0[mat_ind,mat_ind] #Diagonal(diag(H_loc[Corr_orbital_Ind[i],Corr_orbital_Ind[i],s]))

        push!(hyb_w0,hyb_iWn_Block)
        
    end
    return hyb_w0
end



function cal_norm_implev(hyb_w0,H0,Corr_orbital_Ind) 
    norm_implev = []
    
    for  i =1:size(hyb_w0)[1]
        norm_implev_block = zeros(ComplexF64,size(hyb_w0[i]))
        mat_ind = Corr_orbital_Ind[i]
    
        norm_implev_block = H0[mat_ind,mat_ind] + 0.5.* ( hyb_w0[i]+hyb_w0[i]' )
        push!(norm_implev, norm_implev_block)
    end
    
    
    return norm_implev
    
    
    
end




function SelfE_at_w0(SelfE_w,w0ind,Hstart,Corr_atom_Ind,imp_ind,Corr_atom_equiv)
    SelfE_w0 = zeros(size(Hstart))+zeros(size(Hstart))*im
    for j = 1:size(Corr_atom_Ind)[1]
        
        block_S = zeros(size(imp_ind[abs(Corr_atom_equiv[j])])[1],size(imp_ind[abs(Corr_atom_equiv[j])])[1])+zeros(size(imp_ind[abs(Corr_atom_equiv[j])])[1],size(imp_ind[abs(Corr_atom_equiv[j])])[1])*im
        for l =1:Corr_ineq_orbital_Num[abs(Corr_atom_equiv[j])]
                block_S[findall(x->x==l,  imp_ind[abs(Corr_atom_equiv[j])] )].= SelfE_w[abs(Corr_atom_equiv[j])][w0ind,l,1]
        end
        
        SelfE_w0[Corr_orbital_Ind[j],Corr_orbital_Ind[j]] = block_S;
        

    end
    return SelfE_w0
end




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




#### Non-interacting Greenfunction section ####

@everywhere function nonInt_H_R2k_internal(args) # Kpoint paralle
    Kvect = args[1]
    Rvect = args[2]
    H_R = args[3]
    
    H_k_tmp =zeros(ComplexF64,size(H_R[1]))
    #H_r_tmp =  (H_k[j]*exp(-2*pi*im*sum(Rvect.*Kvect)))
    for R = 1:size(Rvect)[1]
        H_k_tmp += (H_R[R]*exp(-2*pi*im*sum(Rvect[R,:].*Kvect)));
    end
    return H_k_tmp
end


@everywhere function nonInt_H_k_OpenMX(args)
    Ktuple = args
    
    totorb =sum(hamiltonian_info.scf_r.Total_NumOrbs)
    H_dum = zeros(ComplexF64, totorb, totorb,2)
    
    Hinfo = DFTforge.cal_colinear_eigenstate(Ktuple,hamiltonian_info,[1,2])
    
    H_dum[:,:,1] = Hinfo[1].Hamiltonian
    H_dum[:,:,2] = Hinfo[2].Hamiltonian
    return H_dum
end




function nonInt_H_k(Klist,hamiltonian_info,Calculation_mode,DMFT_spin_type)
    #initiaize H_k
    
    if string(hamiltonian_info.dfttype) == "Wannier90"
        
        if ((Calculation_mode == "Jx-DMFT") || (Calculation_mode == "Jx0") || (DMFT_spin_type==2) )
            H_k=zeros(size(Klist)[1],size(Klist)[2],size(Klist)[3],sum(hamiltonian_info.scf_r.Total_NumOrbs),sum(hamiltonian_info.scf_r.Total_NumOrbs),2)+zeros(size(Klist)[1],size(Klist)[2],size(Klist)[3],sum(hamiltonian_info.scf_r.Total_NumOrbs),sum(hamiltonian_info.scf_r.Total_NumOrbs),2)*im
        else
            H_k=zeros(size(Klist)[1],size(Klist)[2],size(Klist)[3],sum(hamiltonian_info.scf_r.Total_NumOrbs),sum(hamiltonian_info.scf_r.Total_NumOrbs),1)+zeros(size(Klist)[1],size(Klist)[2],size(Klist)[3],sum(hamiltonian_info.scf_r.Total_NumOrbs),sum(hamiltonian_info.scf_r.Total_NumOrbs),1)*im
        end
        arg_list=[]
        #load H_k
        println("===== loading H_k from tight-binding parameter ... ======") 
        for k1 = (1:size(Klist)[1])
            for k2 = (1:size(Klist)[2])
                for k3 = (1:size(Klist)[3])
                    for s=1:size(H_k)[6]
                        Kvect = Klist[k1,k2,k3,:]
                        push!(arg_list,( Kvect, hamiltonian_info.scf_r.R_vector_mat[s], hamiltonian_info.scf_r.Hks_R[s] ) )
                    
                        #Hinfo = DFTforge.cal_colinear_Hamiltonian(Tuple(Klist[k1,k2,k3,:]),hamiltonian_info.dfttype,hamiltonian_info,s)
                        #H_k[k1,k2,k3,:,:,s]=Hinfo;
                    end
                end
            end
        end
    
        H_k_list_tmp = pmap(nonInt_H_R2k_internal, arg_list)
        i=0
        for k1 = (1:size(Klist)[1])
            for k2 = (1:size(Klist)[2])
                for k3 = (1:size(Klist)[3])
                    for s=1:size(H_k)[6]
                        i=i+1
                        H_k[k1,k2,k3,:,:,s]=H_k_list_tmp[i];
                        #Hinfo = DFTforge.cal_colinear_Hamiltonian(Tuple(Klist[k1,k2,k3,:]),hamiltonian_info.dfttype,hamiltonian_info,s)
                        #H_k[k1,k2,k3,:,:,s]=Hinfo;
                    end
                end
            end
        end    
        if ((Calculation_mode == "Jx-DMFT") || (DMFT_spin_type == 2))
            H_k[:,:,:,:,:,2]= H_k[:,:,:,:,:,1]
        end
        
    elseif string(hamiltonian_info.dfttype) == "OpenMX"
        
        H_k=zeros(size(Klist)[1],size(Klist)[2],size(Klist)[3],sum(hamiltonian_info.scf_r.Total_NumOrbs),sum(hamiltonian_info.scf_r.Total_NumOrbs),hamiltonian_info.scf_r.SpinP_switch+1)+zeros(size(Klist)[1],size(Klist)[2],size(Klist)[3],sum(hamiltonian_info.scf_r.Total_NumOrbs),sum(hamiltonian_info.scf_r.Total_NumOrbs),hamiltonian_info.scf_r.SpinP_switch+1)*im

        arg_list=[]
        #load H_k
        println("===== loading H_k from tight-binding parameter ... ======") 
        for k1 = (1:size(Klist)[1])
            for k2 = (1:size(Klist)[2])
                for k3 = (1:size(Klist)[3])
                    Ktuple = tuple(Klist[k1,k2,k3,:]...)
                    push!(arg_list, Ktuple )
                    
                    #Hinfo = DFTforge.cal_colinear_Hamiltonian(Tuple(Klist[k1,k2,k3,:]),hamiltonian_info.dfttype,hamiltonian_info,s)
                    #H_k[k1,k2,k3,:,:,s]=Hinfo;
                end
            end
        end
        H_k_list_tmp = pmap(nonInt_H_k_OpenMX, arg_list)
        i=0
        for k1 = (1:size(Klist)[1])
            for k2 = (1:size(Klist)[2])
                for k3 = (1:size(Klist)[3])
                    i=i+1
                    H_k[k1,k2,k3,:,:,1]=H_k_list_tmp[i][:,:,1];
                    H_k[k1,k2,k3,:,:,2]=H_k_list_tmp[i][:,:,2];
                    #Hinfo = DFTforge.cal_colinear_Hamiltonian(Tuple(Klist[k1,k2,k3,:]),hamiltonian_info.dfttype,hamiltonian_info,s)
                    #H_k[k1,k2,k3,:,:,s]=Hinfo;
                end
            end
        end            
        
        
        
        #load H_k
        #for k1 = 1:size(Klist)[1]
        #    for k2 = 1:size(Klist)[2]
        #        for k3 = 1:size(Klist)[3]
        #            Hinfo = DFTforge.cal_colinear_eigenstate(tuple(Klist[k1,k2,k3,:]...),hamiltonian_info,[1,2])
        #            H_k[k1,k2,k3,:,:,1]=Hinfo[1].Hamiltonian
        #            H_k[k1,k2,k3,:,:,2]=Hinfo[2].Hamiltonian
        #        end
        #    end
        #end        #for TB hamiltonian

              
    
    end
    
    
    
    return H_k
end
        
function DC_shift(H_k,DCM)
    
    for k1 = 1:size(H_k)[1]
        for k2 = 1:size(H_k)[2]
            for k3 = 1:size(H_k)[3]
                for s =1:size(H_k)[6]
                
                    H_k[k1,k2,k3,:,:,s].-=DCM
                end
            end
        end
    end
    return H_k
end


@everywhere function init_variables_H_k(H_k)
    global g_H_k
   
    g_H_k = H_k;

end

@everywhere function init_variables_SelfE_w(SelfE_w)
    global g_SelfE_w
   
    #g_SelfE_w = SelfE_w
    g_SelfE_w = SelfE_w;

end


@everywhere function init_variables_grid_DMFT(iWnlist,Rlist,Klist)
    global g_iWnlist,g_Rlist,g_Klist
   
    #g_SelfE_w = SelfE_w
    g_iWnlist = iWnlist;
    g_Rlist = Rlist;
    g_Klist = Klist;

end


@everywhere function init_variables_grid(Rlist,Klist)
    global g_Rlist,g_Klist
   
    #g_SelfE_w = SelfE_w
    g_Rlist = Rlist;
    g_Klist = Klist;

end



@everywhere function init_variables_DMFT(Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv, imp_ind, Corr_ineq_orbital_Num, beta, delta)
    global g_Corr_atom_Ind, g_Corr_orbital_Ind, g_Corr_atom_equiv, g_imp_ind, g_Corr_ineq_orbital_Num, g_beta, g_delta
   
    g_Corr_atom_Ind = Corr_atom_Ind;
    g_Corr_orbital_Ind = Corr_orbital_Ind;
    g_Corr_atom_equiv = Corr_atom_equiv;
    g_imp_ind = imp_ind;
    g_Corr_ineq_orbital_Num = Corr_ineq_orbital_Num;
    g_beta = beta;
    g_delta = delta;
    
end

@everywhere function passing_dc(imp_dc)
    global g_imp_dc
   
    g_imp_dc = imp_dc
end

@everywhere function init_variables_Jx(Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind, beta, Corr_ineq_orbital_Num, delta)
    global g_Corr_atom_Ind, g_Corr_orbital_Ind, g_Corr_atom_equiv,g_imp_ind,g_beta, g_Corr_ineq_orbital_Num, g_delta
   
    g_Corr_atom_Ind = Corr_atom_Ind;
    g_Corr_orbital_Ind = Corr_orbital_Ind;
    g_Corr_atom_equiv = Corr_atom_equiv;
    g_imp_ind = imp_ind;
    g_Corr_ineq_orbital_Num = Corr_ineq_orbital_Num;
    g_beta = beta;
    g_delta = delta;
    
end


@everywhere function init_variables_Jx2(Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind, beta, Corr_ineq_orbital_Num, delta)
    global g_Corr_atom_Ind, g_Corr_orbital_Ind, g_Corr_atom_equiv,g_imp_ind,g_beta, g_Corr_ineq_orbital_Num, g_delta
   
    g_Corr_atom_Ind = Corr_atom_Ind;
    g_Corr_orbital_Ind = Corr_orbital_Ind;
    g_Corr_atom_equiv = Corr_atom_equiv;
    g_imp_ind = imp_ind;
    g_Corr_ineq_orbital_Num = Corr_ineq_orbital_Num;
    g_beta = beta;
    g_delta = delta;
    
end



@everywhere function init_variables_totorb(Total_orb_Num)
    global g_Total_orb_Num
   
    #g_SelfE_w = SelfE_w
    g_Total_orb_Num = Total_orb_Num;


end




#### +DMFT section ####


@everywhere function GreenF_at_agrid_dc(iWn,Kind,spin,SelfE_agrid,testmu)
    Gmat = zeros(size(g_H_k[1,1,1,:,:,1]))+zeros(size(g_H_k[1,1,1,:,:,1]))*im 
    Gmat = inv( ( (iWn+testmu)*Matrix{Float64}(I,size(g_H_k[1,1,1,:,:,1]) )) - g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin] - SelfE_agrid + dcMat() )
    return Gmat
end




@everywhere function GreenF_at_agrid(iWn,Kind,spin,SelfE_agrid,testmu)
    Gmat = zeros(size(g_H_k[1,1,1,:,:,1]))+zeros(size(g_H_k[1,1,1,:,:,1]))*im 
    Gmat = inv( ( (iWn+testmu)*Matrix{Float64}(I,size(g_H_k[1,1,1,:,:,1]) )) - g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin] - SelfE_agrid )
    return Gmat
end


@everywhere function GreenF_at_agrid_OpenMX(iWn,Kind,spin,testmu)
    Gmat = zeros(size(g_H_k[1,1,1,:,:,1]))+zeros(size(g_H_k[1,1,1,:,:,1]))*im 
    Gmat = inv( ( (iWn+testmu)*Matrix{Float64}(I,size(g_H_k[1,1,1,:,:,1]) )) - g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin]  )
    return Gmat
end


@everywhere function SelfE_at_agrid(Wind,spin)
    SelfE_agrid = zeros(size(g_H_k[1,1,1,:,:,1]))+zeros(size(g_H_k[1,1,1,:,:,1]))*im
    for j = 1:size(g_Corr_atom_Ind)[1]
        
        block_S = zeros(size(g_imp_ind[abs(g_Corr_atom_equiv[j])])[1],size(g_imp_ind[abs(g_Corr_atom_equiv[j])])[1])+zeros(size(g_imp_ind[abs(g_Corr_atom_equiv[j])])[1],size(g_imp_ind[abs(g_Corr_atom_equiv[j])])[1])*im
        for l =1:g_Corr_ineq_orbital_Num[abs(g_Corr_atom_equiv[j])]
            if g_Corr_atom_equiv[j] >0
                block_S[findall(x->x==l,  g_imp_ind[abs(g_Corr_atom_equiv[j])] )].= g_SelfE_w[abs(g_Corr_atom_equiv[j])][Wind,l,spin]
            else
                block_S[findall(x->x==l,  g_imp_ind[abs(g_Corr_atom_equiv[j])] )].= g_SelfE_w[abs(g_Corr_atom_equiv[j])][Wind,l,mod(spin,2)+1]
            end
        end
        
        SelfE_agrid[g_Corr_orbital_Ind[j],g_Corr_orbital_Ind[j]] = block_S;
        
        #for orb_n =1:size(g_Corr_orbital_Ind[j])[1]
        #    SelfE_agrid[g_Corr_orbital_Ind[j][orb_n],g_Corr_orbital_Ind[j][orb_n]] =  block_S[findall(x->x>0,  g_imp_ind[abs(g_Corr_atom_equiv[j])] )][orb_n]
        #end
        
        #if g_Corr_atom_equiv[j] >0
        #   SelfE_agrid[g_Corr_orbital_Ind[j],g_Corr_orbital_Ind[j]] += Diagonal(g_SelfE_w[abs(g_Corr_atom_equiv[j])][Wind,:,spin])
        #else
        #   SelfE_agrid[g_Corr_orbital_Ind[j],g_Corr_orbital_Ind[j]] += Diagonal(g_SelfE_w[abs(g_Corr_atom_equiv[j])][Wind,:,mod(spin,2)+1])
        #end
    end
    return SelfE_agrid
end





@everywhere function G_k_sum(args)
    w = args[1]
    spin = args[2]
    testmu = args[3]

    G_w_s=zeros(g_Total_orb_Num,g_Total_orb_Num)+zeros(g_Total_orb_Num,g_Total_orb_Num)*im
    for k1 = 1:size(g_Klist)[1]
        for k2 = 1:size(g_Klist)[2]
            for k3 = 1:size(g_Klist)[3]
                G_w_s += GreenF_at_agrid(g_iWnlist[w],[k1 k2 k3], spin, SelfE_at_agrid(w,spin),testmu)
            end
        end
    end
    return G_w_s
end


function Gloc_gen(H_k,iWnlist,SelfE_w,testmu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind)
    G_loc_iWn = zeros(size(iWnlist)[1],size(H_k)[4],size(H_k)[5],size(H_k)[6])+zeros(size(iWnlist)[1],size(H_k)[4],size(H_k)[5],size(H_k)[6])*im
    arg_list = [];
    for w = 1:size(iWnlist)[1]
        for s = 1:size(H_k)[6]
        push!(arg_list,( w,s, testmu,) )
        end
    end

    G_w_s_list = pmap(G_k_sum, arg_list)
    i=0
    for w = 1:size(iWnlist)[1]
        for s = 1:size(H_k)[6]
            i+=1
            G_loc_iWn[w,:,:,s]=G_w_s_list[i] /( size(H_k)[1]*size(H_k)[2]*size(H_k)[3] )
        end
    end

    return  G_loc_iWn
end




function reduce_G_loc(G_loc_iWn,Corr_atom_equiv,Corr_orbital_Ind,imp_ind,DMFT_spin_type)
    red_G_iWn =[]
    for i in unique(abs.(Corr_atom_equiv))
        orb_ind= Corr_orbital_Ind[findall(x->x==i,Corr_atom_equiv)][1]
        matdim=count(x->x>0, unique(imp_ind[i]) )
        red_G_iWn_block =zeros(ComplexF64,size(G_loc_iWn)[1],matdim,DMFT_spin_type)
        add_ind = findall(x->x==1,imp_ind[i])[1]
        for j in unique(imp_ind[i])
            if j>0
                for s =1:DMFT_spin_type
                    add_ind = findall(x->x==j,imp_ind[i])[1]
                    red_G_iWn_block[:,j,s]=G_loc_iWn[:,orb_ind[add_ind[1]],orb_ind[add_ind[2]],s]
                end
            end
        end

       push!(red_G_iWn,red_G_iWn_block)
    end
    return red_G_iWn
end


function reduce_G_loc_mat(G_loc_iWn,Corr_atom_Ind,Corr_orbital_Ind,DMFT_spin_type)
    red_G_iWn_mat =[]
    for i = 1:size(Corr_atom_Ind)[1]
        orb_ind= Corr_orbital_Ind[i]
        matdim=size(orb_ind)[1]
        red_G_iWn_block =zeros(ComplexF64,size(G_loc_iWn)[1],matdim,matdim,DMFT_spin_type)
        for s = 1:DMFT_spin_type
            red_G_iWn_block[:,:,:,s] = G_loc_iWn[:, orb_ind,orb_ind,s]
        end
        
        #add_ind = findall(x->x==1,imp_ind[i])[1]
        #for j in unique(imp_ind[i])
        #    if j>0
        #        for s =1:DMFT_spin_type
        #            add_ind = findall(x->x==j,imp_ind[i])[1]
        #            red_G_iWn_block[:,j,s]=G_loc_iWn[:,orb_ind[add_ind[1]],orb_ind[add_ind[2]],s]
        #        end
        #    end
        #end

       push!(red_G_iWn_mat,red_G_iWn_block)
    end
    return red_G_iWn_mat
end






@everywhere function Cal_Spectral_atEf_givenMu_givenK(args)
    Kind = args[1];
    spin = args[2];
    testmu = args[3];

    
    dummy = zeros(size(g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin]))+ zeros(size(g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin]))*im;
    Iden=Matrix{ComplexF64}(I,size(g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin]))

    MatG_cut_R = GreenF_at_agrid_dc(g_iWnlist[end],Kind,spin,SelfE_at_agrid(size(g_iWnlist)[1],spin),testmu)
    MatG_cut_L = transpose(conj(MatG_cut_R))
    #MatG_cut_L = GreenF_at_agrid(g_iWnlist[1],Kind,spin,SelfE_at_agrid(1,spin),testmu)
    gtilda2 = (g_iWnlist[end]*g_iWnlist[end]) * (MatG_cut_L+MatG_cut_R)/2.0
    gtilda3 = (g_iWnlist[end]*g_iWnlist[end]*g_iWnlist[end]) * ( (MatG_cut_R-MatG_cut_L)/2.0 - Iden/g_iWnlist[end] )
    for w =  1 : (size(g_iWnlist)[1])
        #println("w:",w,"end-w+1:",size(iWnlist)[1]-w+1)
        G_w = GreenF_at_agrid_dc(g_iWnlist[w],Kind,spin,SelfE_at_agrid(w,spin),testmu)
        G_w_t = transpose(conj(GreenF_at_agrid_dc(g_iWnlist[w],Kind,spin,SelfE_at_agrid(w,spin),testmu)))
        #G_rev_w = transpose(conj(GreenF_at_agrid(g_iWnlist[end-w+1],Kind,spin,SelfE_at_agrid(size(g_iWnlist)[1]-w+1,spin),testmu)))
        #dummy=dummy - (1/beta)*(G_k_iWn[k1,k2,k3,w,:,:,s] - Iden/iWnlist[w] - gtilda2/(iWnlist[w]*iWnlist[w]) - gtilda3/(iWnlist[w]*iWnlist[w]*iWnlist[w]) ) * exp(-iWnlist[w]*(beta-delta)) 
        #dummy=dummy - (1/beta)*(G_k_iWn[k1,k2,k3,end-w+1,:,:,s] + Iden/iWnlist[w] - gtilda2/(iWnlist[w]*iWnlist[w]) + gtilda3/(iWnlist[w]*iWnlist[w]*iWnlist[w]) ) * exp(iWnlist[w]*(beta-delta)) 
        dummy=dummy + (1/g_beta)*(G_w - Iden/g_iWnlist[w] - gtilda2/(g_iWnlist[w]*g_iWnlist[w]) - gtilda3/(g_iWnlist[w]*g_iWnlist[w]*g_iWnlist[w]) ) * exp(-g_iWnlist[w]*(g_beta-g_delta)/2.0) 
        dummy=dummy + (1/g_beta)*(G_w_t + Iden/g_iWnlist[w] - gtilda2/(g_iWnlist[w]*g_iWnlist[w]) + gtilda3/(g_iWnlist[w]*g_iWnlist[w]*g_iWnlist[w]) ) * exp(g_iWnlist[w]*(g_beta-g_delta)/2.0) 
    end
    asymtoticG= -Iden/2.0 + gtilda3*g_beta*g_beta/16.0
    dummy=dummy+asymtoticG
    return -(g_beta)/pi*dummy
    
end


function Cal_Spectral_atEf_givenMu(iWnlist,H_k,SelfE_w,testmu,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv)
    Spec_atEf=zeros(size(H_k)[4],size(H_k)[5],size(H_k)[6])+zeros(size(H_k)[4],size(H_k)[5],size(H_k)[6])*im
    Iden=Matrix{ComplexF64}(I,size(H_k)[4],size(H_k)[5])
    
    
    for s= 1 : size(H_k)[6]
        arg_list=[]
        for k1 =1 : size(H_k)[1]
            for k2 =1 : size(H_k)[2]
                for k3 =1 : size(H_k)[3]
                    #push!(arg_list,( iWnlist, H_k[k1,k2,k3,:,:,s], SelfE_w, testmu, beta, s,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv, delta) )
                    push!(arg_list,([k1 k2 k3],s,testmu) )

                end
            end
        end
        
        dummy_list_tmp = pmap(Cal_Spectral_atEf_givenMu_givenK, arg_list)
        
        Spec_atEf[:,:,s]=sum(dummy_list_tmp)/(size(H_k)[1]*size(H_k)[2]*size(H_k)[3])
        
    end
    return Spec_atEf
end








@everywhere function Cal_Occupation_givenMu_givenK(args)
    Kind = args[1];
    spin = args[2];
    testmu = args[3];

    
    dummy = zeros(size(g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin]))+ zeros(size(g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin]))*im;
    Iden=Matrix{ComplexF64}(I,size(g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin]))

    MatG_cut_R = GreenF_at_agrid_dc(g_iWnlist[end],Kind,spin,SelfE_at_agrid(size(g_iWnlist)[1],spin),testmu)
    MatG_cut_L = transpose(conj(MatG_cut_R))
    #MatG_cut_L = GreenF_at_agrid(g_iWnlist[1],Kind,spin,SelfE_at_agrid(1,spin),testmu)
    gtilda2 = (g_iWnlist[end]*g_iWnlist[end]) * (MatG_cut_L+MatG_cut_R)/2.0
    gtilda3 = (g_iWnlist[end]*g_iWnlist[end]*g_iWnlist[end]) * ( (MatG_cut_R-MatG_cut_L)/2.0 - Iden/g_iWnlist[end] )
    for w =  1 : (size(g_iWnlist)[1]-1)
        #println("w:",w,"end-w+1:",size(iWnlist)[1]-w+1)
        G_w = GreenF_at_agrid_dc(g_iWnlist[w],Kind,spin,SelfE_at_agrid(w,spin),testmu)
        G_w_t = transpose(conj(GreenF_at_agrid_dc(g_iWnlist[w],Kind,spin,SelfE_at_agrid(w,spin),testmu)))
        #G_rev_w = transpose(conj(GreenF_at_agrid(g_iWnlist[end-w+1],Kind,spin,SelfE_at_agrid(size(g_iWnlist)[1]-w+1,spin),testmu)))
        #dummy=dummy - (1/beta)*(G_k_iWn[k1,k2,k3,w,:,:,s] - Iden/iWnlist[w] - gtilda2/(iWnlist[w]*iWnlist[w]) - gtilda3/(iWnlist[w]*iWnlist[w]*iWnlist[w]) ) * exp(-iWnlist[w]*(beta-delta)) 
        #dummy=dummy - (1/beta)*(G_k_iWn[k1,k2,k3,end-w+1,:,:,s] + Iden/iWnlist[w] - gtilda2/(iWnlist[w]*iWnlist[w]) + gtilda3/(iWnlist[w]*iWnlist[w]*iWnlist[w]) ) * exp(iWnlist[w]*(beta-delta)) 
        dummy=dummy - (1/g_beta)*(G_w - Iden/g_iWnlist[w] - gtilda2/(g_iWnlist[w]*g_iWnlist[w]) - gtilda3/(g_iWnlist[w]*g_iWnlist[w]*g_iWnlist[w]) ) * exp(-g_iWnlist[w]*(g_beta-g_delta)) 
        dummy=dummy - (1/g_beta)*(G_w_t + Iden/g_iWnlist[w] - gtilda2/(g_iWnlist[w]*g_iWnlist[w]) + gtilda3/(g_iWnlist[w]*g_iWnlist[w]*g_iWnlist[w]) ) * exp(g_iWnlist[w]*(g_beta-g_delta)) 
    end
    asymtoticG= Iden/2.0 - gtilda2*g_beta*((g_beta-g_delta)/(2.0*g_beta)-1.0/4.0) + gtilda3*(g_beta*g_beta/4.0)*((g_beta-g_delta)*(g_beta-g_delta)/(g_beta*g_beta) - (g_beta-g_delta)/(g_beta))
    dummy=dummy+asymtoticG
    return dummy
    
end


function Cal_Occupation_givenMu(iWnlist,H_k,SelfE_w,testmu,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv)
    Occup=zeros(size(H_k)[4],size(H_k)[5],size(H_k)[6])+zeros(size(H_k)[4],size(H_k)[5],size(H_k)[6])*im
    H_loc=zeros(size(H_k)[4],size(H_k)[5],size(H_k)[6])+zeros(size(H_k)[4],size(H_k)[5],size(H_k)[6])*im
    Iden=Matrix{ComplexF64}(I,size(H_k)[4],size(H_k)[5])
    
    
    for s= 1 : size(H_k)[6]
        arg_list=[]
        for k1 =1 : size(H_k)[1]
            for k2 =1 : size(H_k)[2]
                for k3 =1 : size(H_k)[3]
                    #push!(arg_list,( iWnlist, H_k[k1,k2,k3,:,:,s], SelfE_w, testmu, beta, s,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv, delta) )
                    push!(arg_list,([k1 k2 k3],s,testmu) )

                    H_loc[:,:,s] += H_k[k1,k2,k3,:,:,s]/(size(H_k)[1]*size(H_k)[2]*size(H_k)[3])
                end
            end
        end
        
        dummy_list_tmp = pmap(Cal_Occupation_givenMu_givenK, arg_list)
        
        Occup[:,:,s]=sum(dummy_list_tmp)/(size(H_k)[1]*size(H_k)[2]*size(H_k)[3])
        
    end
    return Occup, H_loc
end







@everywhere function Cal_Occupation_givenMu_givenK_WOdc(args)
    Kind = args[1];
    spin = args[2];
    testmu = args[3];

    
    dummy = zeros(size(g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin]))+ zeros(size(g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin]))*im;
    Iden=Matrix{ComplexF64}(I,size(g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin]))

    MatG_cut_R = GreenF_at_agrid(g_iWnlist[end],Kind,spin,SelfE_at_agrid(size(g_iWnlist)[1],spin),testmu)
    MatG_cut_L = transpose(conj(MatG_cut_R))
    #MatG_cut_L = GreenF_at_agrid(g_iWnlist[1],Kind,spin,SelfE_at_agrid(1,spin),testmu)
    gtilda2 = (g_iWnlist[end]*g_iWnlist[end]) * (MatG_cut_L+MatG_cut_R)/2.0
    gtilda3 = (g_iWnlist[end]*g_iWnlist[end]*g_iWnlist[end]) * ( (MatG_cut_R-MatG_cut_L)/2.0 - Iden/g_iWnlist[end] )
    for w =  1 : (size(g_iWnlist)[1])
        #println("w:",w,"end-w+1:",size(iWnlist)[1]-w+1)
        G_w = GreenF_at_agrid(g_iWnlist[w],Kind,spin,SelfE_at_agrid(w,spin),testmu)
        G_w_t = transpose(conj(GreenF_at_agrid(g_iWnlist[w],Kind,spin,SelfE_at_agrid(w,spin),testmu)))
        #G_rev_w = transpose(conj(GreenF_at_agrid(g_iWnlist[end-w+1],Kind,spin,SelfE_at_agrid(size(g_iWnlist)[1]-w+1,spin),testmu)))
        #dummy=dummy - (1/beta)*(G_k_iWn[k1,k2,k3,w,:,:,s] - Iden/iWnlist[w] - gtilda2/(iWnlist[w]*iWnlist[w]) - gtilda3/(iWnlist[w]*iWnlist[w]*iWnlist[w]) ) * exp(-iWnlist[w]*(beta-delta)) 
        #dummy=dummy - (1/beta)*(G_k_iWn[k1,k2,k3,end-w+1,:,:,s] + Iden/iWnlist[w] - gtilda2/(iWnlist[w]*iWnlist[w]) + gtilda3/(iWnlist[w]*iWnlist[w]*iWnlist[w]) ) * exp(iWnlist[w]*(beta-delta)) 
        dummy=dummy - (1/g_beta)*(G_w - Iden/g_iWnlist[w] - gtilda2/(g_iWnlist[w]*g_iWnlist[w]) - gtilda3/(g_iWnlist[w]*g_iWnlist[w]*g_iWnlist[w]) ) * exp(-g_iWnlist[w]*(g_beta-g_delta)) 
        dummy=dummy - (1/g_beta)*(G_w_t + Iden/g_iWnlist[w] - gtilda2/(g_iWnlist[w]*g_iWnlist[w]) + gtilda3/(g_iWnlist[w]*g_iWnlist[w]*g_iWnlist[w]) ) * exp(g_iWnlist[w]*(g_beta-g_delta)) 
    end
    asymtoticG= Iden/2.0 - gtilda2*g_beta/4.0
    dummy=dummy+asymtoticG
    return dummy
    
end


function Cal_Occupation_givenMu_WOdc(iWnlist,H_k,SelfE_w,testmu,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv)
    Occup=zeros(size(H_k)[4],size(H_k)[5],size(H_k)[6])+zeros(size(H_k)[4],size(H_k)[5],size(H_k)[6])*im
    H_loc=zeros(size(H_k)[4],size(H_k)[5],size(H_k)[6])+zeros(size(H_k)[4],size(H_k)[5],size(H_k)[6])*im
    Iden=Matrix{ComplexF64}(I,size(H_k)[4],size(H_k)[5])
    
    
    for s= 1 : size(H_k)[6]
        arg_list=[]
        for k1 =1 : size(H_k)[1]
            for k2 =1 : size(H_k)[2]
                for k3 =1 : size(H_k)[3]
                    #push!(arg_list,( iWnlist, H_k[k1,k2,k3,:,:,s], SelfE_w, testmu, beta, s,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv, delta) )
                    push!(arg_list,([k1 k2 k3],s,testmu) )

                    H_loc[:,:,s] += H_k[k1,k2,k3,:,:,s]/(size(H_k)[1]*size(H_k)[2]*size(H_k)[3])
                end
            end
        end
        
        dummy_list_tmp = pmap(Cal_Occupation_givenMu_givenK_WOdc, arg_list)
        
        Occup[:,:,s]=sum(dummy_list_tmp)/(size(H_k)[1]*size(H_k)[2]*size(H_k)[3])
        
    end
    return Occup, H_loc
end





function Find_mu(iWnlist,H_k,SelfE_w,targetN,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv,DMFT_spin_type,init_mu)
    criterian = 0.000001 # occupation criterian
    
    #Emin= -15.0
    #Emax =15.0
    Emin=init_mu[1]
    Emax=init_mu[2]
        

    xgrid = range(Emin,stop=Emax,length=1000000)


    x = []
    y = []
    xd =[]
    yd =[]
    
    
    print("========= Chemical potential finding in range [",Emin,",",Emax, "] ==========\n")
    print("")

    testmu = Emin
    Occup, Hloc=Cal_Occupation_givenMu_WOdc(iWnlist,H_k,SelfE_w,testmu,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv);
    push!(x,testmu)

    if DMFT_spin_type == 1
        Tot_Occup =  sum(real(diag(Occup[:,:,1])))*2
    else
        Tot_Occup =  sum(real(diag(Occup[:,:,1]+Occup[:,:,2])))
    end    
    push!(y,  Tot_Occup )
    i = 1
    @printf " [%3.0f ] Test mu : %11.7f  ==>  Occup : %10.6f  |  target N : %4.1f \n" i testmu  Tot_Occup  targetN
    
    if ( abs.(targetN - Tot_Occup) <  criterian )
         print("-----------------------------------------------------------------------\n")
         print("          Chemical potential is founded :", testmu,"\n"                   )
         print("-----------------------------------------------------------------------\n")
         
        return testmu, Occup, Hloc
    end

    
    
    

    testmu = Emax
    Occup, Hloc=Cal_Occupation_givenMu_WOdc(iWnlist,H_k,SelfE_w,testmu,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv);
    push!(x,testmu)

    if DMFT_spin_type == 1
        Tot_Occup =  sum(real(diag(Occup[:,:,1])))*2
    else
        Tot_Occup =  sum(real(diag(Occup[:,:,1]+Occup[:,:,2])))
    end    
    push!(y,  Tot_Occup )
    i = 2
    @printf " [%3.0f ] Test mu : %11.7f  ==>  Occup : %10.6f  |  target N : %4.1f \n" i testmu  Tot_Occup  targetN
    
    
    if ( abs.(targetN - Tot_Occup) <  criterian )
         print("-----------------------------------------------------------------------\n")
         print("          Chemical potential is founded :", testmu,"\n"                   )
         print("-----------------------------------------------------------------------\n")
         
        return testmu, Occup, Hloc
    end    

    

    spl = Spline1D(x, y;  k=1)



    for i=3:100
        delta = abs(x[end]-x[end-1])*2.0
        rightend = minimum([testmu+delta, Emax])
        leftend = maximum([testmu-delta, Emin])
        xgrid = range(leftend,stop=rightend,length=1000000)
    
        testmu = xgrid[argmin(abs.(spl(xgrid).-targetN))]

        Occup, Hloc=Cal_Occupation_givenMu_WOdc(iWnlist,H_k,SelfE_w,testmu,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv);
        push!(x,testmu)

        if DMFT_spin_type == 1
            Tot_Occup =  sum(real(diag(Occup[:,:,1])))*2
        else
            Tot_Occup =  sum(real(diag(Occup[:,:,1]+Occup[:,:,2])))
        end    
        push!(y,  Tot_Occup )
        @printf " [%3.0f ] Test mu : %11.7f  ==>  Occup : %10.6f  |  target N : %4.1f \n" i testmu  Tot_Occup  targetN
    
        
        
        xd = deepcopy(x);
        yd = deepcopy(y);
    
        p = sortperm(xd)
        sort!(xd)
        yd=yd[p]
    
        if ( abs.(targetN - Tot_Occup) <  criterian )
            print("-----------------------------------------------------------------------\n")
            print("          Chemical potential is founded :", testmu,"\n"                   )
            print("-----------------------------------------------------------------------\n")
            
            return testmu, Occup, Hloc
        end

        if i > 5
            spl = Spline1D(xd, yd;  k=3)
        elseif i > 4
            spl = Spline1D(xd, yd;  k=2)
        else
            spl = Spline1D(xd, yd;  k=1)
        end

    end
    
    
    
end



function Cal_hyb_fromInvWiess(InvWiess_iWn,imp_ind,iWnlist,mu,H_loc,DMFT_spin_type,Corr_atom_equiv,Corr_orbital_Ind)
    hyb_iWn=[]
    for i =1:size(InvWiess_iWn)[1]
        hyb_iWn_Block = zeros(ComplexF64,size(InvWiess_iWn[i]) )
        mat_ind = Corr_orbital_Ind[findall(x->x==i, Corr_atom_equiv)[1]]
        for w = 1:size(iWnlist)[1]
            for s =1:DMFT_spin_type
                for j in unique(imp_ind[ i ])
                    if j != 0
                        orb_ind =  findall(x->x == j, imp_ind[ i ] )[1]
                        hyb_iWn_Block[w,j,s]=(iWnlist[w]+mu)-InvWiess_iWn[i][w,j,s]-H_loc[mat_ind,mat_ind,s][orb_ind] #Diagonal(diag(H_loc[Corr_orbital_Ind[i],Corr_orbital_Ind[i],s]))
                    end
                end
            end
        end
        push!(hyb_iWn,hyb_iWn_Block)
        
    end
    return hyb_iWn
end


function Cal_hybMat_fromInvWiessMat(InvWiess_iWn_mat,imp_ind,iWnlist,mu,H_loc,DMFT_spin_type,Corr_atom_equiv,Corr_orbital_Ind)
    hyb_iWn_mat=[]
    for i =1:size(InvWiess_iWn_mat)[1]
        hyb_mat_dum = zeros(ComplexF64,size(InvWiess_iWn_mat[i]) )
        mat_ind = Corr_orbital_Ind[i]
        for w = 1:size(iWnlist)[1]
            for s =1:DMFT_spin_type                    
                hyb_mat_dum[w,:,:,s]=(iWnlist[w]+mu)*Matrix{ComplexF64}(I,size(mat_ind)[1],size(mat_ind)[1])-InvWiess_iWn_mat[i][w,:,:,s]-H_loc[mat_ind,mat_ind,s] #Diagonal(diag(H_loc[Corr_orbital_Ind[i],Corr_orbital_Ind[i],s]))
            end
        end
        push!(hyb_iWn_mat,hyb_mat_dum)
        
    end
    return hyb_iWn_mat
end





@everywhere function Cal_InvWiess_fromGloc_at_w(args)
    i = args[1];
    w = args[2];
    G_loc_w = args[3];
    spin = args[4];    
    
            
    imp_block_ind = findall(x->x>0, g_imp_ind[i])
    oin=[]
    for bs =1:size(imp_block_ind)[1]
        push!(oin,imp_block_ind[bs][1])
    end
    
    matdim=count(x->x>0, unique(g_imp_ind[i]) )
    cormatdim = size(G_loc_w)[1]
    Gloc_inv_block_agrid = zeros(ComplexF64,cormatdim,cormatdim)
    Gloc_inv_block_agrid[oin,oin]=inv(G_loc_w[oin,oin])
    
    #InvWiess_iWn_block_agrid= inv(G_loc_w);
    InvWiess_ineq_comp = zeros(ComplexF64,matdim)
    for j =1:matdim #in unique(g_imp_ind[ i ])
        if j > 0
            orb_ind = findall(x->x == j, g_imp_ind[ i ] )[1]
            InvWiess_ineq_comp[j] = Gloc_inv_block_agrid[orb_ind] + g_SelfE_w[i][w,j,spin]
        end
    end
    
    return InvWiess_ineq_comp
end



function Cal_InvWiess_fromGloc(iWnlist,G_loc_iWn,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,DMFT_spin_type)
    InvWiess_iWn=[]
    arg_list=[];
    for i in unique(abs.(Corr_atom_equiv))
        

        o_ind = Corr_orbital_Ind[findall(x->x==i, Corr_atom_equiv)[1]]
        
        for w = 1:size(iWnlist)[1]
            for s = 1:DMFT_spin_type
                push!(arg_list,( i,w, G_loc_iWn[w,o_ind,o_ind,s],s ) )
                #push!(arg_list,( i,w, G_loc_iWn[w,o_ind,o_ind,s],s) )
            end
        end
    end
    
    InvWiess_w_list = pmap(Cal_InvWiess_fromGloc_at_w, arg_list)
    c = 0;
    for i in unique(abs.(Corr_atom_equiv))
        matdim=count(x->x>0, unique(imp_ind[i]) )
        InvWiess_iWn_block=zeros(ComplexF64,size(G_loc_iWn)[1],matdim,DMFT_spin_type)
        for w = 1:size(iWnlist)[1]
            for s = 1:DMFT_spin_type
                c=c+1
                InvWiess_iWn_block[w,:,s]=InvWiess_w_list[c]
            end
        end
        push!(InvWiess_iWn,InvWiess_iWn_block)
    end
    
    return InvWiess_iWn
end




@everywhere function Cal_InvWiessMat_fromGlocMat_at_w(args)
    o_ind = args[1];
    w = args[2];
    red_G_loc_iWn_mat = args[3];
    spin = args[4];


    InvWiessMat_dum = inv(red_G_loc_iWn_mat) + SelfE_at_agrid(w,spin)[o_ind,o_ind]

    return InvWiessMat_dum
end



@everywhere function Cal_InvWiessMat_fromGlocMat(iWnlist,red_G_loc_iWn_mat,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,DMFT_spin_type)
    InvWiess_iWn_mat=[]
    arg_list = []

    for i = 1:size(red_G_loc_iWn_mat)[1]
        o_ind = Corr_orbital_Ind[i]

        for w = 1:size(iWnlist)[1]
            for s = 1:DMFT_spin_type
                push!(arg_list,( o_ind, w, red_G_loc_iWn_mat[i][w,:,:,s], s ) )
                #InvWiessMat_dum[w,:,:,s] = inv(red_G_loc_iWn_mat[i][w,:,:,s]) + SelfE_at_agrid(w,s)[o_ind,o_ind]
                #push!(arg_list,( i,w, G_loc_iWn[w,o_ind,o_ind,s],s ) )
                #push!(arg_list,( i,w, G_loc_iWn[w,o_ind,o_ind,s],s) )
            end
        end
    end



    InvWiessMat_w_list = pmap(Cal_InvWiessMat_fromGlocMat_at_w, arg_list)
    c = 0;
    for i = 1:size(red_G_loc_iWn_mat)[1]
        InvWiessMat_dum = zeros(ComplexF64,size(red_G_loc_iWn_mat[i]))
        for w = 1:size(iWnlist)[1]
            for s = 1:DMFT_spin_type
                c=c+1
                InvWiessMat_dum[w,:,:,s] = InvWiessMat_w_list[c]
            end
        end
        push!(InvWiess_iWn_mat,InvWiessMat_dum)
    end



    return InvWiess_iWn_mat
end



function Cal_interacting_Greenf(G0_k_iWn,SelfE_diag_w)
    G_k_iWn=zeros(size(G0_k_iWn)[1],size(G0_k_iWn)[2],size(G0_k_iWn)[3],size(G0_k_iWn)[4],size(G0_k_iWn)[5],size(G0_k_iWn)[6],size(G0_k_iWn)[7])+zeros(size(G0_k_iWn)[1],size(G0_k_iWn)[2],size(G0_k_iWn)[3],size(G0_k_iWn)[4],size(G0_k_iWn)[5],size(G0_k_iWn)[6],size(G0_k_iWn)[7])*im
    for s = 1:size(G0_k_iWn)[7]
        for w = 1:size(G0_k_iWn)[4]
            for k1 = 1: size(G0_k_iWn)[1]
                for k2 = 1: size(G0_k_iWn)[2]
                    for k3 = 1: size(G0_k_iWn)[3]
                        G_k_iWn[k1,k2,k3,w,:,:,s]= G0_k_iWn[k1,k2,k3,w,:,:,s]+Diagonal(SelfE_diag_w[w,:,s])
                    end
                end
            end
        end 
    end
    return G_k_iWn
end




@everywhere function V_k_iWn_agrid(args)
   kind = args[1]
   wind = args[2]


   V_k_iWn_ =  ( g_H_k[kind[1],kind[2],kind[3],:,:,1] + SelfE_at_agrid(wind,1) ) - (g_H_k[kind[1],kind[2],kind[3],:,:,2] +  SelfE_at_agrid(wind,2) )

   return V_k_iWn_
end


function Cal_spinpotential_V(H_k,iWnlist,wind)

    V_k_iWn = zeros(size(H_k)[1],size(H_k)[2],size(H_k)[3],size(wind)[1],size(H_k)[4],size(H_k)[5])+zeros(size(H_k)[1],size(H_k)[2],size(H_k)[3],size(wind)[1],size(H_k)[4],size(H_k)[5])*im
    arg_list =[]
    for k1 = 1:size(V_k_iWn)[1]
        for k2 = 1:size(V_k_iWn)[2]
            for k3 = 1:size(V_k_iWn)[3]
                for wi in wind
                    push!(arg_list, ([k1 k2 k3], wi))
                end
            end
        end
    end
    V_k_iWn_list = pmap(V_k_iWn_agrid, arg_list)
    c=0

    for k1 = 1:size(V_k_iWn)[1]
        for k2 = 1:size(V_k_iWn)[2]
            for k3 = 1:size(V_k_iWn)[3]
                for n = 1:size(wind)[1]
                    c=c+1
                    V_k_iWn[k1,k2,k3,n,:,:]=V_k_iWn_list[c]
                end
            end
        end
    end

    return V_k_iWn
end


#### J cal section ####



function iWnlist_div(iWnlist,div)
    divleng = ceil(Int64,size(iWnlist)[1]/div)
    divind = []
    for i=1:div
        if i<div
            push!(divind, divleng*(i-1)+1:divleng*(i)) 
        else
            push!(divind, divleng*(i-1)+1:size(iWnlist)[1]-2) 
        end
    end
    
    
    
    return divind
end


function init_JR(R_grid_num, atom12_list, Corr_atom_Ind, Corr_orbital_Ind)
    JR = []
    JR_orb = []
    distance = []
    for i =1:size(atom12_list)[1]
        conv_orbind1 = findall(x->x == atom12_list[i][1], Corr_atom_Ind)[1]
        conv_orbind2 = findall(x->x == atom12_list[i][2], Corr_atom_Ind)[1]
    
        orb1_size = size(Corr_orbital_Ind[conv_orbind1])[1]
        orb2_size = size(Corr_orbital_Ind[conv_orbind2])[1]
    
        
        
        JR_=zeros(R_grid_num[1],R_grid_num[2],R_grid_num[3])+zeros(R_grid_num[1],R_grid_num[2],R_grid_num[3])*im
        JR_orb_=zeros(R_grid_num[1],R_grid_num[2],R_grid_num[3],orb1_size,orb2_size )+zeros(R_grid_num[1],R_grid_num[2],R_grid_num[3],orb1_size,orb2_size)*im
        distance_=zeros(R_grid_num[1],R_grid_num[2],R_grid_num[3])

        push!(JR,JR_)
        push!(JR_orb,JR_orb_)
        push!(distance,distance_)
        
    end
        
    return JR, JR_orb, distance
end



function J_print(atom12_list,sorted_cell_ind, sorted_R_grid_ind,sorted_dist, sorted_J, sorted_J_orb)
    for i = 1:size(atom12_list)[1]
        file_n = string("Jx_raw_",atom12_list[i][1],"_",atom12_list[i][2],".dat")
        open(file_n, "w") do io
            @printf(io,"  R1   R2   R3      I1   I2   I3      A1   A2           dist.               Jij\n")
            for j = 1:size(sorted_J[i])[1]
                #if sorted_dist[i][j] != 0.0
                    #@printf(io,"\n\n\n\n\nDistance: %10.6f,    Cell Ind. [%2i,%2i,%2i],    Atom Ind (%2i,%2i)\n\n",sorted_dist[i][j],sorted_cell_ind[i][j][1],sorted_cell_ind[i][j][2],sorted_cell_ind[i][j][3],atom12_list[i][1],atom12_list[i][2] )
                    @printf(io,"%4i %4i %4i    %4i %4i %4i    %4i %4i    %18.14f %18.14f\n",sorted_cell_ind[i][j][1],sorted_cell_ind[i][j][2],sorted_cell_ind[i][j][3],sorted_R_grid_ind[i][j][1],sorted_R_grid_ind[i][j][2],sorted_R_grid_ind[i][j][3],atom12_list[i][1],atom12_list[i][2], sorted_dist[i][j], real(sorted_J[i][j])  )
                #end
            end    
        end
    end    
    
    
    for i = 1:size(atom12_list)[1]
        file_n = string("Jx_",atom12_list[i][1],"_",atom12_list[i][2],".dat")
        open(file_n, "w") do io
            for j = 1:size(sorted_J[i])[1]
                if sorted_dist[i][j] != 0.0
                    @printf(io,"%10.6f %10.6f\n", sorted_dist[i][j], real(sorted_J[i][j])  )
                end
            end    
        end
    end
    
    open("zeroline.dat", "w") do io
        @printf(io,"%10.6f %10.6f\n", 0.0, 0.0)
        @printf(io,"%10.6f %10.6f\n", 1000.0, 0.0)
    end


    open("Jx.gnu", "w") do io
        @printf(io,"set xra [0:20]\n")
        @printf(io,"plot \"zeroline.dat\" u 1:2 w l lw 2 lt 0,")
        for i = 1:size(atom12_list)[1]
            file_n = string("Jx_",atom12_list[i][1],"_",atom12_list[i][2],".dat")
            if i !=size(atom12_list)[1]
                @printf(io,"\"%s\" u 1:2 w linespoints pt 7 ps 1.7 lw 1.3, ", file_n)
            else
                @printf(io,"\"%s\" u 1:2 w linespoints pt 7 ps 1.7 lw 1.3\n ", file_n)
            end
        end
        @printf(io,"pause -1")
    end
    

    for i = 1:size(atom12_list)[1]
        file_n = string("Jx_orb_",atom12_list[i][1],"_",atom12_list[i][2],".dat")
        open(file_n, "w") do io
            for j = 1:20
                if sorted_dist[i][j] != 0.0
                    @printf(io,"\n\n\n\n\nDistance: %10.6f,    Cell Ind. [%2i,%2i,%2i],    Atom Ind (%2i,%2i)\n\n",sorted_dist[i][j],sorted_cell_ind[i][j][1],sorted_cell_ind[i][j][2],sorted_cell_ind[i][j][3],atom12_list[i][1],atom12_list[i][2] )
                    @printf(io,"   J =          ")
                    for orb2 = 1:size(sorted_J_orb[1][1])[2]
                        @printf(io,"orb%1i    ",orb2)
                    end
                    @printf(io,"\n")
                    @printf(io,"        ")
                    for orb1 =1:size(sorted_J_orb[1][1])[1]
                        @printf(io,"orb%1i ",orb1)
                        for orb2 = 1:size(sorted_J_orb[1][1])[2]
                            @printf(io,"%8.4f",real(sorted_J_orb[i][j][orb1,orb2]))
                        end
                        @printf(io,"\n")
                        @printf(io,"        ")
                        
                        if orb1 == size(sorted_J_orb[1][1])[1]
                            @printf(io,"\n")
                            @printf(io,"          J_tot : %8.4f,    J_sum : %8.4f\n", real(sorted_J[i][j]), sum(real(sorted_J_orb[i][j][:,:])) ) 
                            
                        end
                        
                    end
                end
            end    
        end
    end
    
    
    
    
    
    
end


@everywhere function GreenF_gen_at_kw(args)
    iWn = args[1]
    kind = args[2]
    wind = args[3]
    spin = args[4]
    mu =args[5]

    G_k_iWn = GreenF_at_agrid(iWn,kind,spin,SelfE_at_agrid(wind,spin),mu)

    return G_k_iWn
end


@everywhere function GreenF_gen_wDC_at_kw(args)
    iWn = args[1]
    kind = args[2]
    wind = args[3]
    spin = args[4]
    mu =args[5]

    G_k_iWn = GreenF_at_agrid_dc(iWn,kind,spin,SelfE_at_agrid(wind,spin),mu)

    return G_k_iWn
end

function GreenF_gen(H_k,iWnlist,wind,SelfE_w,testmu,Corr_atom_Ind,Corr_orbital_Ind)
    G_k_iWn=zeros(size(H_k)[1],size(H_k)[2],size(H_k)[3],size(wind)[1],size(H_k)[4],size(H_k)[5],size(H_k)[6])+zeros(size(H_k)[1],size(H_k)[2],size(H_k)[3],size(wind)[1],size(H_k)[4],size(H_k)[5],size(H_k)[6])*im
    G_loc_iWn = zeros(size(wind)[1],size(H_k)[4],size(H_k)[5],size(H_k)[6])+zeros(size(wind)[1],size(H_k)[4],size(H_k)[5],size(H_k)[6])*im

    arg_list =[]
    for k1 = 1: size(G_k_iWn)[1]
        for k2 =1: size(G_k_iWn)[2]
            for k3 =1: size(G_k_iWn)[3]
                for wi in wind
                    for s = 1:size(G_k_iWn)[7]
                        push!(arg_list, (iWnlist[wi],[k1 k2 k3],wi,s,testmu ))  
                        

                    end
                end
            end
        end
    end
    
    G_list = pmap(GreenF_gen_at_kw, arg_list)

    c = 0
    for k1 = 1: size(G_k_iWn)[1]
        for k2 =1: size(G_k_iWn)[2]
            for k3 =1: size(G_k_iWn)[3]
                for (n,wi) in enumerate(wind)
                    for s = 1:size(G_k_iWn)[7]
                        c = c+1
                        G_k_iWn[k1,k2,k3,n,:,:,s]= G_list[c]
                        G_loc_iWn[n,:,:,s] += G_k_iWn[k1,k2,k3,n,:,:,s] /(size(H_k)[1]*size(H_k)[2]*size(H_k)[3])

                    end
                end
            end
        end
    end    
    
    
    return G_k_iWn,  G_loc_iWn
end



function GreenF_gen_OpenMX(H_k,iWnlist,testmu,Corr_atom_Ind,Corr_orbital_Ind)
    G_k_iWn=zeros(size(H_k)[1],size(H_k)[2],size(H_k)[3],size(iWnlist)[1],size(H_k)[4],size(H_k)[5],size(H_k)[6])+zeros(size(H_k)[1],size(H_k)[2],size(H_k)[3],size(iWnlist)[1],size(H_k)[4],size(H_k)[5],size(H_k)[6])*im
    G_loc_iWn = zeros(size(iWnlist)[1],size(H_k)[4],size(H_k)[5],size(H_k)[6])+zeros(size(iWnlist)[1],size(H_k)[4],size(H_k)[5],size(H_k)[6])*im


    for k1 = 1: size(G_k_iWn)[1]
        for k2 =1: size(G_k_iWn)[2]
            for k3 =1: size(G_k_iWn)[3]
                for (w,iWn) in enumerate(iWnlist)
                    for s = 1:size(G_k_iWn)[7]
                                         
                        G_k_iWn[k1,k2,k3,w,:,:,s]= GreenF_at_agrid_OpenMX(iWn,[k1 k2 k3],s,testmu)
                        G_loc_iWn[w,:,:,s] += G_k_iWn[k1,k2,k3,w,:,:,s] /(size(H_k)[1]*size(H_k)[2]*size(H_k)[3])

                    end
                end
            end
        end
    end
    return G_k_iWn,  G_loc_iWn
end





function Take_block_function_G(G_k_iWn,atomlist_primitive,Corr_orbital_Ind,atom12_list)
    G_k_iWn_block=[]
    for i=1:size(atom12_list)[1]
        G_k_iWn_atom12=[]
        conv_orbind1 = findall(x->x == atom12_list[i][1], Corr_atom_Ind)[1]
        conv_orbind2 = findall(x->x == atom12_list[i][2], Corr_atom_Ind)[1]

        block_Gij=zeros(size(G_k_iWn)[1],size(G_k_iWn)[2],size(G_k_iWn)[3],size(G_k_iWn)[4],size(Corr_orbital_Ind[conv_orbind1])[1] ,size(Corr_orbital_Ind[conv_orbind2])[1] ,size(G_k_iWn)[7])+zeros(size(G_k_iWn)[1],size(G_k_iWn)[2],size(G_k_iWn)[3],size(G_k_iWn)[4],size(Corr_orbital_Ind[conv_orbind1])[1] ,size(Corr_orbital_Ind[conv_orbind2])[1] ,size(G_k_iWn)[7])*im
        block_Gji=zeros(size(G_k_iWn)[1],size(G_k_iWn)[2],size(G_k_iWn)[3],size(G_k_iWn)[4],size(Corr_orbital_Ind[conv_orbind2])[1] ,size(Corr_orbital_Ind[conv_orbind1])[1] ,size(G_k_iWn)[7])+zeros(size(G_k_iWn)[1],size(G_k_iWn)[2],size(G_k_iWn)[3],size(G_k_iWn)[4],size(Corr_orbital_Ind[conv_orbind2])[1] ,size(Corr_orbital_Ind[conv_orbind1])[1] ,size(G_k_iWn)[7])*im
        block_Gij=deepcopy(G_k_iWn[:,:,:,:,Corr_orbital_Ind[conv_orbind1],Corr_orbital_Ind[conv_orbind2],:])  # Gij of atom12_list[x] = (i,j)
        block_Gji=deepcopy(G_k_iWn[:,:,:,:,Corr_orbital_Ind[conv_orbind2],Corr_orbital_Ind[conv_orbind1],:])  # Gji of atom12_list[x] = (i,j)
        push!(G_k_iWn_atom12,block_Gij)
        push!(G_k_iWn_atom12,block_Gji)

        push!(G_k_iWn_block,G_k_iWn_atom12)
    end
    return G_k_iWn_block
end



function Take_block_function_V(V_k_iWn,atomlist_primitive,Corr_orbital_Ind,atom12_list,spinmode)
    V_k_iWn_block=[]
    for i=1:size(atom12_list)[1]
        V_k_iWn_atom12=[]
        conv_orbind1 = findall(x->x == atom12_list[i][1], Corr_atom_Ind)[1]
        conv_orbind2 = findall(x->x == atom12_list[i][2], Corr_atom_Ind)[1]

        block_Vii=zeros(size(V_k_iWn)[1],size(V_k_iWn)[2],size(V_k_iWn)[3],size(V_k_iWn)[4],size(Corr_orbital_Ind[conv_orbind1])[1] ,size(Corr_orbital_Ind[conv_orbind1])[1],spinmode )+zeros(size(V_k_iWn)[1],size(V_k_iWn)[2],size(V_k_iWn)[3],size(V_k_iWn)[4],size(Corr_orbital_Ind[conv_orbind1])[1] ,size(Corr_orbital_Ind[conv_orbind1])[1] ,spinmode)*im
        block_Vjj=zeros(size(V_k_iWn)[1],size(V_k_iWn)[2],size(V_k_iWn)[3],size(V_k_iWn)[4],size(Corr_orbital_Ind[conv_orbind2])[1] ,size(Corr_orbital_Ind[conv_orbind2])[1],spinmode )+zeros(size(V_k_iWn)[1],size(V_k_iWn)[2],size(V_k_iWn)[3],size(V_k_iWn)[4],size(Corr_orbital_Ind[conv_orbind2])[1] ,size(Corr_orbital_Ind[conv_orbind2])[1] ,spinmode)*im
        block_Vii[:,:,:,:,:,:,1]=deepcopy(V_k_iWn[:,:,:,:,Corr_orbital_Ind[conv_orbind1],Corr_orbital_Ind[conv_orbind1] ])
        block_Vii[:,:,:,:,:,:,2]=deepcopy(-V_k_iWn[:,:,:,:,Corr_orbital_Ind[conv_orbind1],Corr_orbital_Ind[conv_orbind1] ])

        block_Vjj[:,:,:,:,:,:,1]=deepcopy(V_k_iWn[:,:,:,:,Corr_orbital_Ind[conv_orbind2],Corr_orbital_Ind[conv_orbind2] ])
        block_Vjj[:,:,:,:,:,:,2]=deepcopy(-V_k_iWn[:,:,:,:,Corr_orbital_Ind[conv_orbind2],Corr_orbital_Ind[conv_orbind2] ])
        push!(V_k_iWn_atom12,block_Vii)
        push!(V_k_iWn_atom12,block_Vjj)

        push!(V_k_iWn_block,V_k_iWn_atom12)
    end
    return V_k_iWn_block
end






function Cal_VG_product_in_k_iWn(V_k_iWn_block,G_k_iWn_block,atom12_list)
    VG_k_iWn_block=[]
    for i = 1:size(atom12_list)[1]
        VG_k_iWn_atom12=[]
        ViGij_k_iWn=zeros(size(V_k_iWn_block[i][1])[1],size(V_k_iWn_block[i][1])[2],size(V_k_iWn_block[i][1])[3],size(V_k_iWn_block[i][1])[4],size(V_k_iWn_block[i][1])[5] ,size(G_k_iWn_block[i][1])[6], size(G_k_iWn_block[i][1])[7])+zeros(size(V_k_iWn_block[i][1])[1],size(V_k_iWn_block[i][1])[2],size(V_k_iWn_block[i][1])[3],size(V_k_iWn_block[i][1])[4],size(V_k_iWn_block[i][1])[5] ,size(G_k_iWn_block[i][1])[6], size(G_k_iWn_block[i][1])[7])*im 
        VjGji_k_iWn=zeros(size(V_k_iWn_block[i][2])[1],size(V_k_iWn_block[i][2])[2],size(V_k_iWn_block[i][2])[3],size(V_k_iWn_block[i][2])[4],size(V_k_iWn_block[i][2])[5] ,size(G_k_iWn_block[i][2])[6], size(G_k_iWn_block[i][2])[7])+zeros(size(V_k_iWn_block[i][2])[1],size(V_k_iWn_block[i][2])[2],size(V_k_iWn_block[i][2])[3],size(V_k_iWn_block[i][2])[4],size(V_k_iWn_block[i][2])[5] ,size(G_k_iWn_block[i][2])[6], size(G_k_iWn_block[i][2])[7])*im 

        for k1 = 1: size(G_k_iWn_block[i][1])[1]
            for k2 = 1: size(G_k_iWn_block[i][1])[2]
                for k3 = 1: size(G_k_iWn_block[i][1])[3]
                    for w = 1: size(G_k_iWn_block[i][1])[4]
                        # V^{ud}G^{d}
                        ViGij_k_iWn[k1,k2,k3,w,:,:,1]= V_k_iWn_block[i][1][k1,k2,k3,w,:,:]*G_k_iWn_block[i][1][k1,k2,k3,w,:,:,2] # Vi^{ud}Gij^{d}
                        VjGji_k_iWn[k1,k2,k3,w,:,:,1]= V_k_iWn_block[i][2][k1,k2,k3,w,:,:]*G_k_iWn_block[i][2][k1,k2,k3,w,:,:,2] # Vj^{ud}Gji^{d}

                        # V^{du}G^{u}
                        ViGij_k_iWn[k1,k2,k3,w,:,:,2]= (-V_k_iWn_block[i][1][k1,k2,k3,w,:,:])*G_k_iWn_block[i][1][k1,k2,k3,w,:,:,1] # Vi^{du}Gij^{u}
                        VjGji_k_iWn[k1,k2,k3,w,:,:,2]= (-V_k_iWn_block[i][2][k1,k2,k3,w,:,:])*G_k_iWn_block[i][2][k1,k2,k3,w,:,:,1] # Vj^{du}Gji^{u}
                        

                    end
                end
            end
        end
        
        push!(VG_k_iWn_atom12,ViGij_k_iWn)
        push!(VG_k_iWn_atom12,VjGji_k_iWn)
                        
        push!(VG_k_iWn_block,VG_k_iWn_atom12)   
    end
    return VG_k_iWn_block
end



function Cal_J_R_print(V_R_iWn_block,G_R_iWn_block,Rv,iWnlist,beta,atom12_list,hamiltonian_info,mag_order)

    for atomI = 1:size(atom12_list)[1]
        J_ud=0.0+0.0*im
        J_du=0.0+0.0*im
        Jm_ud=0.0+0.0*im
        Jm_du=0.0+0.0*im
        atom12ind=atomI
        for w =1:size(iWnlist)[1]
            J_ud += mag_order[atomI]*1000/(4*beta)/2*tr(V_R_iWn_block[atom12ind][1][1,1,1,w,:,:,1]*G_R_iWn_block[atom12ind][1][1+Rv[1],1+Rv[2],1+Rv[3],w,:,:,1]*V_R_iWn_block[atom12ind][2][1,1,1,w,:,:,2]*G_mR_iWn_block[atom12ind][2][1+Rv[1],1+Rv[2],1+Rv[3],w,:,:,2])
          #  J_ud += mag_order[atomI]*1000/(4*beta)/2*tr(V_R_iWn_block[atom12ind][1][1,1,1,w,:,:,1]'*G_R_iWn_block[atom12ind][1][1+Rv[1],1+Rv[2],1+Rv[3],w,:,:,2]'*V_R_iWn_block[atom12ind][2][1,1,1,w,:,:,2]'*G_mR_iWn_block[atom12ind][2][1+Rv[1],1+Rv[2],1+Rv[3],w,:,:,1]')
            
            J_du += mag_order[atomI]*1000/(4*beta)/2*tr(V_R_iWn_block[atom12ind][1][1,1,1,w,:,:,2]*G_R_iWn_block[atom12ind][1][1+Rv[1],1+Rv[2],1+Rv[3],w,:,:,2]*V_R_iWn_block[atom12ind][2][1,1,1,w,:,:,1]*G_mR_iWn_block[atom12ind][2][1+Rv[1],1+Rv[2],1+Rv[3],w,:,:,1])
                    
            Jm_ud += mag_order[atomI]*1000/(4*beta)/2*tr(V_R_iWn_block[atom12ind][1][1,1,1,w,:,:,1]*G_mR_iWn_block[atom12ind][1][1+Rv[1],1+Rv[2],1+Rv[3],w,:,:,1]*V_R_iWn_block[atom12ind][2][1,1,1,w,:,:,2]*G_R_iWn_block[atom12ind][2][1+Rv[1],1+Rv[2],1+Rv[3],w,:,:,2])
            Jm_du += mag_order[atomI]*1000/(4*beta)/2*tr(V_R_iWn_block[atom12ind][1][1,1,1,w,:,:,2]*G_mR_iWn_block[atom12ind][1][1+Rv[1],1+Rv[2],1+Rv[3],w,:,:,2]*V_R_iWn_block[atom12ind][2][1,1,1,w,:,:,1]*G_R_iWn_block[atom12ind][2][1+Rv[1],1+Rv[2],1+Rv[3],w,:,:,1])
        end
        J=(J_ud+J_du)
        Jm=(Jm_ud+Jm_du)

        println("===============   Atom i,j index: ",atom12_list[atom12ind],"   ===============")
        println("  ----------      Rdiff vect:   ",Rv, "     ----------")
        pos=sqrt(sum((  (hamiltonian_info.scf_r.Gxyz.+Rv'*hamiltonian_info.scf_r.tv)[atom12_list[atom12ind][2],:]-(hamiltonian_info.scf_r.Gxyz.+[0,0,0]'*hamiltonian_info.scf_r.tv)[atom12_list[atom12ind][1],:]  ).^2))
        println("  |Distance:", pos,"|")
        println("")
        println("  J_ud :",J_ud )
        println("  J_du :",J_du )
        println("  J :",J )
        println("")
        println("  ----------      Rdiff vect:   ",-Rv, "     ----------")
        pos=sqrt(sum((  (hamiltonian_info.scf_r.Gxyz.-Rv'*hamiltonian_info.scf_r.tv)[atom12_list[atom12ind][2],:]-(hamiltonian_info.scf_r.Gxyz.+[0,0,0]'*hamiltonian_info.scf_r.tv)[atom12_list[atom12ind][1],:]  ).^2))
        println("  |Distance:", pos,"|")
        println("")
        println("  J_ud :",Jm_ud )
        println("  J_du :",Jm_du )
        println("  J :",Jm )
        println("===========================================================")
        println("")
        println("")

    end
    #return [Jm_ud,Jm_du]
    
end




function Cal_J_R(JR,JR_orb,V_R_iWn_block,G_R_iWn_block,G_mR_iWn_block,Rv,atom_list12_ind,R_grid_num,iWnlist,wind,beta,hamiltonian_info,mag_order)

    Rv_revised=deepcopy(Rv)
    pos=sqrt(sum((  (hamiltonian_info.scf_r.Gxyz.+Rv'*hamiltonian_info.scf_r.tv)[atom12_list[atom_list12_ind][2],:]-(hamiltonian_info.scf_r.Gxyz.+[0,0,0]'*hamiltonian_info.scf_r.tv)[atom12_list[atom_list12_ind][1],:]  ).^2))
    if Rv[1]<0
        Rv_revised[1]=R_grid_num[1]+Rv[1];
    end
    if Rv[2]<0
        Rv_revised[2]=R_grid_num[2]+Rv[2];
    end
    if Rv[3]<0
        Rv_revised[3]=R_grid_num[3]+Rv[3];
    end
    J_ud=0.0+0.0*im
    J_du=0.0+0.0*im
    
    or1_size = size(G_R_iWn_block[atom_list12_ind][1][1,1,1,1,:,1,2])[1]
    or2_size = size(G_R_iWn_block[atom_list12_ind][1][1,1,1,1,1,:,2])[1]
    J_ud_orb = zeros(or1_size,or2_size)+zeros(or1_size,or2_size)*im
    J_du_orb = zeros(or1_size,or2_size)+zeros(or1_size,or2_size)*im

    for w =1:size(wind)[1]
        J_ud += mag_order[atom_list12_ind]*1000/(4*beta)*tr(V_R_iWn_block[atom_list12_ind][1][1,1,1,w,:,:,1]*G_R_iWn_block[atom_list12_ind][1][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,:,1]*V_R_iWn_block[atom_list12_ind][2][1,1,1,w,:,:,2]*G_mR_iWn_block[atom_list12_ind][2][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,:,2])
        J_ud += mag_order[atom_list12_ind]*1000/(4*beta)*tr(V_R_iWn_block[atom_list12_ind][1][1,1,1,w,:,:,1]'*G_R_iWn_block[atom_list12_ind][1][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,:,1]'*V_R_iWn_block[atom_list12_ind][2][1,1,1,w,:,:,2]'*G_mR_iWn_block[atom_list12_ind][2][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,:,2]')  # Summation of minus matsubara frequency 

        J_du += mag_order[atom_list12_ind]*1000/(4*beta)*tr(V_R_iWn_block[atom_list12_ind][1][1,1,1,w,:,:,2]*G_R_iWn_block[atom_list12_ind][1][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,:,2]*V_R_iWn_block[atom_list12_ind][2][1,1,1,w,:,:,1]*G_mR_iWn_block[atom_list12_ind][2][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,:,1])
        J_du += mag_order[atom_list12_ind]*1000/(4*beta)*tr(V_R_iWn_block[atom_list12_ind][1][1,1,1,w,:,:,2]'*G_R_iWn_block[atom_list12_ind][1][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,:,2]'*V_R_iWn_block[atom_list12_ind][2][1,1,1,w,:,:,1]'*G_mR_iWn_block[atom_list12_ind][2][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,:,1]')  # Summation of minus matsubara frequency 

        for orb1=1:or1_size
            for orb2=1:or2_size
                J_ud_orb[orb1,orb2] += (mag_order[atom_list12_ind]*1000/(4*beta))*(transpose(V_R_iWn_block[atom_list12_ind][1][1,1,1,w,orb1,:,1])*(G_R_iWn_block[atom_list12_ind][1][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,orb2,1]))*(transpose(V_R_iWn_block[atom_list12_ind][2][1,1,1,w,orb2,:,2])*(G_mR_iWn_block[atom_list12_ind][2][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,orb1,2]))
                J_ud_orb[orb1,orb2] += (mag_order[atom_list12_ind]*1000/(4*beta))*(transpose(V_R_iWn_block[atom_list12_ind][1][1,1,1,w,:,:,1]'[orb1,:])*(G_R_iWn_block[atom_list12_ind][1][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,:,1]'[:,orb2]))*(transpose(V_R_iWn_block[atom_list12_ind][2][1,1,1,w,:,:,2]'[orb2,:])*(G_mR_iWn_block[atom_list12_ind][2][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,:,2]'[:,orb1]))  # Summation of minus matsubara frequency 

                J_du_orb[orb1,orb2] += (mag_order[atom_list12_ind]*1000/(4*beta))*(transpose(V_R_iWn_block[atom_list12_ind][1][1,1,1,w,orb1,:,2])*(G_R_iWn_block[atom_list12_ind][1][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,orb2,2]))*(transpose(V_R_iWn_block[atom_list12_ind][2][1,1,1,w,orb2,:,1])*(G_mR_iWn_block[atom_list12_ind][2][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,orb1,1]))
                J_du_orb[orb1,orb2] += (mag_order[atom_list12_ind]*1000/(4*beta))*(transpose(V_R_iWn_block[atom_list12_ind][1][1,1,1,w,:,:,2]'[orb1,:])*(G_R_iWn_block[atom_list12_ind][1][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,:,2]'[:,orb2]))*(transpose(V_R_iWn_block[atom_list12_ind][2][1,1,1,w,:,:,1]'[orb2,:])*(G_mR_iWn_block[atom_list12_ind][2][1+Rv_revised[1],1+Rv_revised[2],1+Rv_revised[3],w,:,:,1]'[:,orb1]))  # Summation of minus matsubara frequency 

            end
        end
    
    
    end
    JR[atom_list12_ind][Rv_revised[1]+1,Rv_revised[2]+1,Rv_revised[3]+1] += (J_ud+J_du)./2.0
    JR_orb[atom_list12_ind][Rv_revised[1]+1,Rv_revised[2]+1,Rv_revised[3]+1,:,:] += (J_ud_orb.+J_du_orb)./2.0

    
    return [[Rv_revised[1]+1,Rv_revised[2]+1,Rv_revised[3]+1],pos,JR,JR_orb]
    
end



function J_series(JR,JR_orb,distance,V_R_iWn_block,G_R_iWn_block,G_mR_iWn_block,R_grid_num,iWnlist,wind,beta,hamiltonian_info,mag_order)

    cell_ind=[]
    atom_ind=[]
    
    for atomI=1:size(atom12_list)[1]
        J_atomI=[]
        J_orb_atomI=[]
        distance_atomI=[]
        cell_ind_atomI=[]
        atom_ind_atomI=[]

        for Rx =convert(Int64,-(floor(R_grid_num[1]/2))):convert(Int64,ceil(R_grid_num[1]/2)-1)
            for Ry =convert(Int64,-(floor(R_grid_num[2]/2))):convert(Int64,ceil(R_grid_num[2]/2)-1)
                for Rz =convert(Int64,-(floor(R_grid_num[3]/2))):convert(Int64,ceil(R_grid_num[3]/2)-1)
                    Rv=[Rx,Ry,Rz]
                    posind,pos,JR,JR_orb= Cal_J_R(JR,JR_orb,V_R_iWn_block,G_R_iWn_block,G_mR_iWn_block,Rv,atomI,R_grid_num,iWnlist,wind,beta,hamiltonian_info,mag_order)
                    distance[atomI][posind[1],posind[2],posind[3]]=real(pos)

                end
            end
        end
        

    end
    return [JR,JR_orb,distance]
end



function Jx_calculation(JR, JR_orb,distance,hamiltonian_info,H_k,mag_order,R_grid_num,Rlist,Klist,iWnlist,wdivind,atom12_list,SelfE_w,testmu,Corr_atom_Ind,Corr_orbital_Ind)
    for i = 1:size(wdivind)[1]
            
        println("=========== Calculating J for w-segments [",i,"/",size(wdivind)[1],"] ===========")
        
        @time G_k_iWn, G_loc_iWn = GreenF_gen(H_k,iWnlist,wdivind[i],SelfE_w,testmu,Corr_atom_Ind,Corr_orbital_Ind)
        @time V_R_iWn = Cal_spinpotential_V(H_k,iWnlist,wdivind[i])
    
        GC.gc()

    
        #################  Cut the Gij, Vij with Corr_orbital_Ind block ##############
        println("Cut off the G & V for making block matrix ...")
        atomlist_primitive=duplicatedfnt_atomlist(atom12_list)
        @time G_k_iWn_block=Take_block_function_G(G_k_iWn,atomlist_primitive,Corr_orbital_Ind,atom12_list);
        @time V_R_iWn_block=Take_block_function_V(V_R_iWn,atomlist_primitive,Corr_orbital_Ind,atom12_list,2);

   
        #----- memory save -----#
        G_k_iWn=nothing
        V_k_iWn=nothing
        GC.gc()

    
        ################### fourier transformation of G_block :  k -> R  and k -> -R ###############
        println("Fourier transform of G : k -> R and k -> -R ...")
        @time G_R_iWn_block=FFTW_k2r_block(Klist,Rlist,atom12_list,G_k_iWn_block);
        @time G_mR_iWn_block_fac=FFTW_r2k_block(Rlist,Klist,atom12_list,G_k_iWn_block);
        G_mR_iWn_block=G_mR_iWn_block_fac/(size(Klist)[1]*size(Klist)[2]*size(Klist)[3]);
    
    
        ################### Calculate J  ###############
        @time JR,JR_orb,distance = J_series(JR,JR_orb,distance,V_R_iWn_block,G_R_iWn_block,G_mR_iWn_block,R_grid_num,iWnlist,wdivind[i],beta,hamiltonian_info,mag_order);
        GC.gc()
    end
    
    return [JR, JR_orb,distance]
end


function sort_dist(distance,JR,JR_orb, atom12_list, R_grid_num)
    sorted_cell_ind =[]
    sorted_R_grid_ind=[]
    sorted_JR =[]
    sorted_JR_orb =[]
    sorted_dist =[]
    for at=1:size(atom12_list)[1]
        dist_at=[]
        cell_ind_at=[]
        R_grid_ind = []
        JR_at =[]
        JR_orb_at=[]
        for r1=1:R_grid_num[1]
            for r2=1:R_grid_num[2]
                for r3=1:R_grid_num[3]
                    #if distance[at][r1,r2,r3]>0.0
                        push!(dist_at,distance[at][r1,r2,r3])
                        
                        rr1=copy(r1-1)
                        rr2=copy(r2-1)
                        rr3=copy(r3-1)
                        if (rr1>(ceil(R_grid_num[1]/2)-1))
                            rr1 -= R_grid_num[1]
                        end
                        if (rr2>(ceil(R_grid_num[2]/2)-1))  
                            rr2 -= R_grid_num[2]
                        end
                        if (rr3>(ceil(R_grid_num[3]/2)-1))  
                            rr3 -= R_grid_num[3]   
                        end
                                
                        push!(cell_ind_at,[rr1,rr2,rr3])
                        push!(R_grid_ind,[r1-1,r2-1,r3-1])
                        
                        push!(JR_at,JR[at][r1,r2,r3])
                        push!(JR_orb_at,JR_orb[at][r1,r2,r3,:,:])

                    #end
                end
            end
        end
        p=sortperm(dist_at)
        push!(sorted_dist, dist_at[p])
        push!(sorted_cell_ind, cell_ind_at[p])
        push!(sorted_R_grid_ind, R_grid_ind[p])
        push!(sorted_JR,JR_at[p])
        push!(sorted_JR_orb,JR_orb_at[p])

    end
    return sorted_dist,sorted_cell_ind,sorted_R_grid_ind,sorted_JR,sorted_JR_orb
end



X_VERSION = VersionNumber("0.9.3-pub+20190529");
println(" JX_VERSION: ",X_VERSION)
println(" Visit https://kaist-elst.github.io/DFTforge.jl/ for details & updates ")
println(" Tested with Julia v1.0 and v1.1 which the most recent version of Julia in 201906 https://julialang.org/")

@everywhere using LinearAlgebra
using Distributed
using DFTforge.DFTrefinery
using DFTforge.DFTcommon;
# Julia 1.0
using Statistics
#

@everywhere using Distributed
@everywhere import DFTforge
@everywhere using DFTforge.DFTrefinery
@everywhere using DFTforge.DFTcommon


##############################################################################
## 1. Read INPUT
## 1.1 Set Default values
## 1.2 Read input from argument & TOML file
## 1.3 Set values from intput (arg_input)
## 1.4 Set caluations type and ouput folder
##############################################################################

hdftmpdir = ""
## 1.1 Set Default values
#orbital_selection1 = Array{Int64,1}();
#orbital_selection2 = Array{Int64,1}();
orbital_selection1_list = Array{Array{Int}}(undef,0);
orbital_selection1_names = Array{AbstractString}(undef,0);
orbital_selection2_list = Array{Array{Int}}(undef,0);
orbital_selection2_names = Array{AbstractString}(undef,0);

orbital_selection3_list = Array{Array{Int}}(undef,0);
orbital_selection3_names = Array{AbstractString}(undef,0);
orbital_selection4_list = Array{Array{Int}}(undef,0);
orbital_selection4_names = Array{AbstractString}(undef,0);

orbital_selection_option = DFTcommon.nomask;
orbital_selection_on = false;

band_selection_on = false;
band_selection_upper = 0;
band_selection_lower = 0;

k_point_num = [3,3,3]
q_point_num = [3,3,3]
ChemP_delta_ev = 0.0
DFT_type = DFTcommon.OpenMX

## 1.2 Read input from argument & TOML file
arg_input = DFTcommon.Arg_Inputs();
arg_input = parse_input(ARGS,arg_input)
arg_input = parse_TOML(arg_input.TOMLinput,arg_input)

# let argument override
arg_input = parse_input(ARGS,arg_input)


## 1.3 Set values from intput (arg_input)
DFT_type = arg_input.DFT_type
Wannier90_type = arg_input.Wannier90_type
spin_type = arg_input.spin_type

result_file = arg_input.result_file
result_file_dict = arg_input.result_file_dict;

ChemP_delta_ev = arg_input.ChemP_delta_ev
 # k,q point num
k_point_num = arg_input.k_point_num
q_point_num = arg_input.q_point_num
 # atom 12
#atom1 = arg_input.atom1;
#atom2 = arg_input.atom2;
atom12_list = arg_input.atom12_list;
hdftmpdir = arg_input.hdftmpdir;

# orbital mask
orbital_selection_option = arg_input.orbital_selection_option;

orbital_selection1_list = arg_input.orbital_selection1_list;
orbital_selection1_names = arg_input.orbital_selection1_names;
orbital_selection2_list = arg_input.orbital_selection2_list;
orbital_selection2_names = arg_input.orbital_selection2_names;

orbital_selection3_list = arg_input.orbital_selection3_list;
orbital_selection3_names = arg_input.orbital_selection3_names;
orbital_selection4_list = arg_input.orbital_selection4_list;
orbital_selection4_names = arg_input.orbital_selection4_names;

println(orbital_selection1_list," ",orbital_selection1_names)
println(orbital_selection2_list," ",orbital_selection2_names)
println(orbital_selection3_list," ",orbital_selection3_names)
println(orbital_selection4_list," ",orbital_selection4_names)


@assert(length(orbital_selection1_list) == length(orbital_selection1_names));
@assert(length(orbital_selection2_list) == length(orbital_selection2_names));
@assert(length(orbital_selection3_list) == length(orbital_selection3_names));
@assert(length(orbital_selection4_list) == length(orbital_selection4_names));
# Band selection
if haskey(arg_input.Optional,"band_selection")
  band_selection_on =  arg_input.Optional["band_selection"]
  band_selection_lower =  arg_input.Optional["band_selection_boundary"][1]
  band_selection_upper =  arg_input.Optional["band_selection_boundary"][2]
end

if ((DFTcommon.unmask == orbital_selection_option) || (DFTcommon.mask == orbital_selection_option) )
 #orbital_selection_name = arg_input.orbital_selection_name
 orbital_selection_on = true
end

# Energy windows selection
energywindow_all_list = Array{Array{Float64,1},1}()
energywindow_1_list = Array{Array{Float64,1},1}()
energywindow_2_list = Array{Array{Float64,1},1}()
energywindow_3_list = Array{Array{Float64,1},1}()
energywindow_4_list = Array{Array{Float64,1},1}()
println(DFTcommon.bar_string) # print ====...====
println("energy windows")

if haskey(arg_input.Optional,"energywindow")
  energywindow_all_list = arg_input.Optional["energywindow"]["energywindow_all_list"]
  energywindow_1_list  = arg_input.Optional["energywindow"]["energywindow_1_list"]
  energywindow_2_list  = arg_input.Optional["energywindow"]["energywindow_2_list"]
  energywindow_3_list  = arg_input.Optional["energywindow"]["energywindow_3_list"]
  energywindow_4_list  = arg_input.Optional["energywindow"]["energywindow_4_list"]
  energywindow_name = arg_input.Optional["energywindow"]["energywindow_name"]
  println("energywindow_name: ",energywindow_name)
end
println("energywindow_all_list: ", energywindow_all_list)
println("energywindow_1_list: ", energywindow_1_list)
println("energywindow_2_list: ", energywindow_2_list)
println("energywindow_3_list: ", energywindow_3_list)
println("energywindow_4_list: ", energywindow_4_list)

# orbital orbital_reassign
basisTransform_rule = basisTransform_rule_type()
if haskey(arg_input.Optional,"basisTransform")
  basisTransform_rule = arg_input.Optional["basisTransform"]
end

println(DFTcommon.bar_string) # print ====...====
println("atom12_list: ",atom12_list)
println("q_point_num ",q_point_num, "\tk_point_num ",k_point_num)
println(string("DFT_type ",DFT_type))
println(string("orbital_selection_option ",orbital_selection_option))
println("mask1list ",orbital_selection1_list,"\tmask2list ",orbital_selection2_list)
println("basisTransform", basisTransform_rule)


cal_type = "jx.col.dmft.spin" # xq, ...
## 1.4 Set caluations type and ouput folder
if haskey(arg_input.Optional,"energywindow")
    energywindow_name = arg_input.Optional["energywindow"]["energywindow_name"]
    cal_type = string(cal_type,".Erange_",energywindow_name);
end

if (DFTcommon.Wannier90 == DFT_type)
  cal_type = string(cal_type,".wannier")
end
root_dir = dirname(result_file)
result_file_only = splitext(basename(result_file))
cal_name = result_file_only[1];
jq_output_dir =  joinpath(root_dir,string(cal_type,"_" ,ChemP_delta_ev))
if (!isdir(jq_output_dir))
  mkdir(jq_output_dir)
end
if ("" == hdftmpdir || !isdir(hdftmpdir) )
  hdftmpdir = jq_output_dir
end
hdf_cache_name = joinpath(hdftmpdir,string(cal_name,".hdf5"))
println(hdf_cache_name)
println(DFTcommon.bar_string) # print ====...====



#result_file=string(arg_input.result_file_dict["result_file"],"_hr.dat")

##############################################################################
## 2. Calculate & Store k,q points information
## 2.1 Set Input info
## 2.2 Generate k,q points
## 2.3 Calculate Eigenstate & Store Eigenstate into file in HDF5 format
## 2.4 Send Eigenstate info to child processes
##############################################################################

## 2.1 Set Input info
hamiltonian_info = [];
if (DFTcommon.OpenMX == DFT_type || DFTcommon.EcalJ == DFT_type)
  hamiltonian_info = set_current_dftdataset(result_file,result_file_dict, DFT_type, spin_type,basisTransform_rule)
elseif (DFTcommon.Wannier90 == DFT_type)
  atomnum = arg_input.Wannier_Optional_Info.atomnum
  atompos = arg_input.Wannier_Optional_Info.atompos
  atoms_orbitals_list = arg_input.Wannier_Optional_Info.atoms_orbitals_list

  #hamiltonian_info = DFTforge.read_dftresult(result_file,DFT_type,Wannier90_type,atoms_orbitals_list,atomnum,atompos,basisTransform_rule)
  hamiltonian_info = set_current_dftdataset(result_file,result_file_dict,
  DFT_type,Wannier90_type,spin_type,atoms_orbitals_list,atomnum,atompos,basisTransform_rule)
  #hamiltonian_info = set_current_dftdataset(scf_r, DFT_type, spin_type,basisTransform_rule)
end

DFTforge.pwork(set_current_dftdataset,(hamiltonian_info, 1));





if  ( arg_input.Optional["DMFT_Jx"].Calculation_mode == "DMFT" )

    (Calculation_mode, BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut, Tem, beta, beta_for_occ,
    DMFT_spin_type, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind,
    DMFT_loop_N, Mix_selfE,init_bias, smth_step, cal_susc,compute_EF, DMFT_solver, imp_dc_type, imp_ind, imp_lev_shift, imp_U, imp_J, imp_dc_n, imp_dc, imp_int_type, imp_int_parameterisation
    , Corr_ineq_orbital_Num, Solver_green_Cut ,F0, F2, F4, imp_Measure_time, imp_Thermal_time, basis_transform,consider_degeneracy, green_basis, green_legendre_cutoff
    , EDMFTF_MonteCar_step, EDMFTF_warmup_step, EDMFTF_GlobalFlip, EDMFTF_tsample, EDMFTF_nom, EDMFTF_PChangeOrder) = input_handler(arg_input.Optional["DMFT_Jx"]);
    

    
elseif (arg_input.Optional["DMFT_Jx"].Calculation_mode == "Jx-DMFT")

    (Calculation_mode, BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut, Tem, beta, beta_for_occ,
    DMFT_spin_type, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind,
    DMFT_loop_N, Mix_selfE, init_bias,smth_step,cal_susc,compute_EF, DMFT_solver, imp_dc_type, imp_ind, imp_lev_shift, imp_U, imp_J, imp_dc_n, imp_dc, imp_int_type, imp_int_parameterisation
    , Corr_ineq_orbital_Num, Solver_green_Cut ,F0, F2, F4, imp_Measure_time, imp_Thermal_time, basis_transform,consider_degeneracy, green_basis, green_legendre_cutoff,mag_order
    , EDMFTF_MonteCar_step, EDMFTF_warmup_step, EDMFTF_GlobalFlip, EDMFTF_tsample, EDMFTF_nom, EDMFTF_PChangeOrder) = input_handler(arg_input.Optional["DMFT_Jx"]);
    

    
elseif (arg_input.Optional["DMFT_Jx"].Calculation_mode == "Jx0")

    (Calculation_mode, BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut,Tem, beta, DMFT_spin_type,
                 Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind, Corr_ineq_orbital_Num,imp_ind,mag_order) = input_handler(arg_input.Optional["DMFT_Jx"]);

    ## DMFT part input
    @time DFTforge.pwork(init_variables_Jx2,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind, beta, Corr_ineq_orbital_Num, delta); # passing basic variable for DMFT ++ , here beta_for_occ means for initial occ

elseif (arg_input.Optional["DMFT_Jx"].Calculation_mode == "Magnon")

    (Calculation_mode, qmode, Neighbors_cut, BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut,Tem, beta, DMFT_spin_type,
                 Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind, Corr_ineq_orbital_Num,imp_ind,mag_order) = input_handler(arg_input.Optional["DMFT_Jx"]);

    ## DMFT part input
    @time DFTforge.pwork(init_variables_Jx2,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind, beta, Corr_ineq_orbital_Num, delta); # passing basic variable for DMFT ++ , here beta_for_occ means for initial occ
    
    
end
# SelfE_w = init_selfE_dcbias(imp_dc,Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,iWnlist,init_bias)


        
    #@time DFTforge.pwork(init_variables_DMFT,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind, Corr_ineq_orbital_Num, beta, delta); # passing basic variable for DMFT ++
    #@time DFTforge.pwork(passing_dc, imp_dc)    
    
    
    ################ GRID GENERATION PART ##################
    #atom_orbitals=atom_orb_list(hamiltonian_info);
    #Total_orb_Num=sum(hamiltonian_info.scf_r.Total_NumOrbs)
    #@time DFTforge.pwork(init_variables_totorb,Total_orb_Num); # passing total number of orbital to workerspace

    #Rlist=RPoint_gen(R_grid_num);
    #Klist=kPoint_gen(K_grid_num);
    #iWnlist=Matsubara_gen(iW_grid_cut,beta);
    
    ################ SHOW INPUT ################
    #show_input(K_grid_num, R_grid_num, iW_grid_cut, Tem, DMFT_spin_type, Corr_atom_Ind, Corr_atom_equiv, Solver_green_Cut, basis_transform, consider_degeneracy
    # ,green_basis,init_bias,Mix_selfE ,imp_dc_type,imp_Measure_time, imp_Thermal_time, imp_U, imp_J, imp_dc_n, imp_dc, imp_int_type, imp_int_parameterisation, imp_ind, smth_step, cal_susc)    

    
    
    
    
    ################  non-interacting hamiltonian construction ################ 
    #@time H_k=nonInt_H_k(Klist,hamiltonian_info,DMFT,Jx,DMFT_spin_type)
    #mu=hamiltonian_info.scf_r.ChemP
    #@time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
    #println("")
    

    #Trans_mat = read_transforMat()
    #Trans_mat_A =Tvec_print_from_reading(H_k, Trans_mat, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)
    #@time H_k = Trans_H_k(H_k,Trans_mat_A)



    #@time H_loc= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind)
    #imp_ind = imp_ind_trans_by_degeneracy(imp_ind,H0,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind)
    #Corr_ineq_orbital_Num = cal_corr_ineq_orbital_N(imp_ind)
    #DFTforge.pwork(init_variables_DMFT,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind, Corr_ineq_orbital_Num, beta, delta); # passing basic variable for DMFT ++
  
    #SelfE_w = init_selfE_dcbias(imp_dc,Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,iWnlist,0.0)            
    #GC.gc()



    #SelfE_w = loadSelfE_from_obs(Ineq_atom_Ind,DMFT_spin_type,SelfE_w)
    #SelfE_w_sm =SelfE_smoothing(Corr_atom_equiv,Tem,iWnlist,SelfE_w); # save smoothing version of SelfE
    #DFTforge.pwork(init_variables_SelfE_w,SelfE_w_sm); #passing the smoothing SelfE to worker


    #Plots.plot(imag(iWnlist),imag(SelfE_w_sm[1][:,1,1]))
    #Plots.plot!(imag(iWnlist),imag(SelfE_w_sm[1][:,2,1]))
    #Plots.plot!(imag(iWnlist),imag(SelfE_w_sm[1][:,3,1]))
    #Plots.plot!(imag(iWnlist),imag(SelfE_w_sm[2][:,1,1]))



    #DFTforge.pwork(init_variables_grid_DMFT,iWnlist,Rlist,Klist); # passing grid variable to worker
    #@time G_loc_iWn=Gloc_gen(H_k,iWnlist,SelfE_w_sm,mu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind)



###########################################################################
######################### Self-consistent DMFT part #######################
###########################################################################

if  ( Calculation_mode == "DMFT" )


    println("#####################################################")
    println("      Calculation type : self-consistent DMFT        ")
    println("#####################################################\n\n\n")
        
        
    @time DFTforge.pwork(init_variables_DMFT,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind, Corr_ineq_orbital_Num, beta, delta); # passing basic variable for DMFT ++
    @time DFTforge.pwork(passing_dc, imp_dc)    
    
    
    ################ GRID GENERATION PART ##################
    atom_orbitals=atom_orb_list(hamiltonian_info);
    Total_orb_Num=sum(hamiltonian_info.scf_r.Total_NumOrbs)
    @time DFTforge.pwork(init_variables_totorb,Total_orb_Num); # passing total number of orbital to workerspace

    Rlist=RPoint_gen(R_grid_num);
    Klist=kPoint_gen(K_grid_num);
    iWnlist=Matsubara_gen(iW_grid_cut,beta);
    
    ################ SHOW INPUT ################
    show_input(K_grid_num, R_grid_num, iW_grid_cut, Tem, DMFT_spin_type, Corr_atom_Ind, Corr_atom_equiv, Solver_green_Cut, basis_transform, consider_degeneracy
     ,green_basis,init_bias,Mix_selfE ,imp_dc_type,imp_Measure_time, imp_Thermal_time, imp_U, imp_J, imp_dc_n, imp_dc, imp_int_type, imp_int_parameterisation, imp_ind, smth_step, cal_susc)    

    
    
    
    
    ################  non-interacting hamiltonian construction ################ 
    @time H_k=nonInt_H_k(Klist,hamiltonian_info,Calculation_mode,DMFT_spin_type)
    mu=hamiltonian_info.scf_r.ChemP
    DCM =  dcMat_noglob(H_k,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,imp_ind)
    H0_k = DC_shift(H_k,DCM)
    @time DFTforge.pwork(init_variables_H_k,H0_k);  #passing variable H_k to worker
    println("")
    
    



    ############### start with DMFT_loop = 0 ###############
    DMFT_loop = 0
    

    ##  create real frequency grid & selfenergy in real frequency
    wnum = 3000
    wlist= wlist_gen(15,wnum)
    w0ind=convert(Int64,wnum/2)
    SelfE_realw = init_selfE_dc(imp_dc,Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,wlist)
    lattmu=0.0
    gaussian_broaden=0.025
    GC.gc()  
            
    @time DFTforge.pwork(init_variables_grid_DMFT_realf,wlist,Klist); # passing grid variable to worker for determining initial occ
    @time DFTforge.pwork(init_variables_SelfE_w,SelfE_realw); #passing variable SelfE to worker for determining initial occ
    GC.gc()
    
    #--- Check if the file "Nele.dat" exist or not. If there is no "Nele.dat", calculate Nele at low temperature (T=100K)
    if (  !isfile("Nele.dat") || !isfile("transform_Mat.dat") )
        G_loc_w= Gloc_gen_realf(H0_k,wlist,SelfE_realw,lattmu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind,gaussian_broaden)  
        GC.gc()

        Aw, targetN, dft_corrN = Aw_orb(Corr_orbital_Ind,Corr_atom_equiv,wlist, G_loc_w);
        write_Nele(targetN)
        @time H0= local_H0(H0_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind)

        if !isfile("Hloc_orig.dat")
            Hloc_print("Hloc_orig.dat", H0 ,Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)
        end          

        
        if ( basis_transform == "spec" ||  basis_transform == "hyb" )
            SelfE_w0 = SelfE_at_w0(SelfE_realw,w0ind,H0_k[1,1,1,:,:,1],Corr_atom_Ind,imp_ind,Corr_atom_equiv)
            @time InvWiess_w0 = Cal_InvWiess_fromGloc_0(SelfE_w0,w0ind,G_loc_w[w0ind,:,:,1],Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,DMFT_spin_type);
            @time hyb_w0=Cal_hyb_fromInvWiess_0(InvWiess_w0,imp_ind,wlist[w0ind],lattmu,H0,DMFT_spin_type,Corr_atom_equiv,Corr_orbital_Ind,gaussian_broaden);
            @time norm_implev = cal_norm_implev(hyb_w0,H0,Corr_orbital_Ind);
        
            println(" ==== Basis transformation :", basis_transform, "====\n")
            Trans_mat, Trans_mat_A = Tvec_print(Aw[w0ind,:,:,1],norm_implev ,Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind, basis_transform)
            
        else
            Trans_mat, Trans_mat_A = Tvec_print_I(H0_k, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)
        end        
        
    
    else # If "Nele.dat" exist, just read it 
        targetN=parse(Float64, read_Nele())[1]
        dft_corrN= read_DFT_corrNele()     
        Trans_mat = read_transforMat()
        Trans_mat_A =Tvec_print_from_reading(H0_k, Trans_mat, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)
    end   
    GC.gc()           

    # Obtaining transformation matrix in complex spherical harmonics
    Trans_mat_C, Trans_mat_A_C = R2Csphericalharmonics(Trans_mat, Trans_mat_A, Corr_orbital_Ind)
   
 
    # generating (totalorb * totalorb) dim impurity level shift matrix 
    lev_shift_mat=imp_level_shift(imp_lev_shift,hamiltonian_info.scf_r.Hks_R[1][1], Corr_atom_Ind, Corr_atom_equiv, imp_ind)
                    

    # passing new hamiltonian transformed by spectral-diagonalized basis
    @time H0_k = Trans_H_k(H0_k,Trans_mat_A)
    @time H0_k = Shift_H_k(H0_k,lev_shift_mat)

    DFTforge.pwork(init_variables_H_k,H0_k);  #passing variable H_k to worker
    GC.gc() 


    #SelfE_w_for_occ = nothing
    #iWnlist_for_occ = nothing  # remove dense iw and selfE variable used in determining initial occ 
    #GC.gc()                 
    
    
    
    SelfE_realw = nothing;
    wlist = nothing;
    G_loc_w = nothing;
    Aw = nothing;
    InvWiess_w0 = nothing;
    hyb_w0 = nothing;
    norm_implev = nothing;
    GC.gc() 
    
    #----- init self-energy as double counting without bias for a given temperature written in input file   
    SelfE_w = init_selfE_dcbias(imp_dc,Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,iWnlist,0.0)            
   
            
    #----- passing grid inform & selfE & input var -----#
    DFTforge.pwork(init_variables_grid_DMFT,iWnlist,Rlist,Klist); # passing grid variable to worker
    DFTforge.pwork(init_variables_SelfE_w,SelfE_w); #passing variable SelfE to worker
    DFTforge.pwork(init_variables_DMFT,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind, Corr_ineq_orbital_Num, beta, delta); # passing basic variable for DMFT ++
    GC.gc()                
            
     
    if compute_EF
        init_mu = [mu-1.0, mu+1.0 ]
        @time begin
            print("          ===== Finding mu ... =====\n")
            mu, Occup, H_loc =Find_mu(iWnlist,H0_k,SelfE_w,targetN,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv,DMFT_spin_type,init_mu)
            println("")
        end
    else
        @time begin
            mu = 0.0
            Occup, H_loc=Cal_Occupation_givenMu_WOdc(iWnlist,H0_k,SelfE_w,mu,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv);
        end
    end

    print("-------------------------------------\n")
    print("Initial latt. mu in transform basis :", mu,"\n" ) 
    print("-------------------------------------\n")
    Hloc_print("Hloc_trans.dat", H_loc ,Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)
            
    #write_mu_NonIntH(mu) # wirte "mu_NonIntH.dat" file
    GC.gc()    
                
 
      
    #if "consider_degeneracy" turned on, the generate states would be treated equivalently 
    if consider_degeneracy == true
        imp_ind = imp_ind_trans_by_degeneracy(imp_ind,H_loc,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind)
        Corr_ineq_orbital_Num = cal_corr_ineq_orbital_N(imp_ind)
        DFTforge.pwork(init_variables_DMFT,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind, Corr_ineq_orbital_Num, beta, delta); # passing basic variable for DMFT ++
    end
 

    # ----- change double counintg if imp_dc_type = "FLL-DFT"  ------ #
    imp_totN = deepcopy(imp_dc_n)
    imp_dc = dc_change(imp_totN, imp_dc_n, dft_corrN, imp_dc_type, imp_U, imp_J) # calculate double counting
    DFTforge.pwork(passing_dc, imp_dc) # passing double counting
        
    # ----- initial selfenergy construction with DC+bias ------ #
    SelfE_w = init_selfE_dcbias(imp_dc,Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,iWnlist,init_bias)    
    DFTforge.pwork(init_variables_SelfE_w,SelfE_w); #passing variable SelfE to worker


    GC.gc()    
    
        
    ################ Check if Restart with already existing selfEnergy or not #################
    restart_wSelfE = Check_restart_wSelfE(Ineq_atom_Ind, DMFT_solver)
    if  (restart_wSelfE) # if in all imp_x folder there is output file params.obs.json, then just read self energy 
        SelfE_w_new = loadSelfE(Ineq_atom_Ind,DMFT_spin_type,SelfE_w,DMFT_solver)
        
        imp_totN = loadOccup_from_obs_scalar(Ineq_atom_Ind,DMFT_solver);
        
        imp_dc = dc_change(imp_totN, imp_dc_n, dft_corrN, imp_dc_type, imp_U, imp_J) # calculate double counting
        DFTforge.pwork(passing_dc, imp_dc) # passing double counting
        if DMFT_solver == "ComCTQMC"
            Save_logfile("params.obs",Ineq_atom_Ind)
        elseif DMFT_solver == "EDMFTF_ctqmc"
            Save_logfile("ctqmc",Ineq_atom_Ind)
            Save_logfile("Sig",Ineq_atom_Ind)
            Save_logfile("Gf",Ineq_atom_Ind)
        end
        write_outfile("SelfE",iWnlist,SelfE_w_new,Corr_atom_equiv,DMFT_spin_type,imp_ind) #write initial selfenergy sig.dat

        if go_bareSelfE(smth_step,Ineq_atom_Ind,DMFT_solver)
            SelfE_w = deepcopy(SelfE_w_new)
        else
            SelfE_w =SelfE_smoothing(Corr_atom_equiv,Tem,iWnlist,SelfE_w_new); # save smoothing version of SelfE
            write_outfile("SelfE_sm",iWnlist,SelfE_w,Corr_atom_equiv,DMFT_spin_type,imp_ind) #write initial selfenergy sig_smth.dat
        end
        
        DFTforge.pwork(init_variables_SelfE_w,SelfE_w); #passing the smoothing SelfE to worker
    
    else
       write_outfile("SelfE",iWnlist,SelfE_w,Corr_atom_equiv,DMFT_spin_type,imp_ind) #write initial selfenergy sig_bare.dat
       #SelfE_w_sm = deepcopy(SelfE_w)
    end

    GC.gc()
    #write_outfile("SelfE_sm",iWnlist,SelfE_w_sm,Corr_atom_equiv,DMFT_spin_type,imp_ind) #write initial selfenergy sig_smth.dat


    
    
    #----- Find mu again with non-zero self-energy (non-zero means non-constant_double_counting self-energy) -----#
           #----- From now on, self energy is readed selfenergy from params.obs.json or self energy biased 
    
    if isfile("scf.log")
        mulog = readmu_fromscflog()
    else
        mulog =[]
    end

    if compute_EF
        init_mu = mu_guess(DMFT_loop, mulog, DMFT_spin_type)
        @time begin
            print("          ===== Finding mu ... =====\n")
            mu, Occup, H_loc =Find_mu(iWnlist,H0_k,SelfE_w,targetN,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv,DMFT_spin_type,init_mu)
            println("")
        end
    else
        @time begin
            mu = 0.0
            Occup, H_loc=Cal_Occupation_givenMu_WOdc(iWnlist,H0_k,SelfE_w,mu,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv);
        end
    end
    
    push!(mulog,mu)
    

    print("-------------------------\n")
    print("Initial latt. mu :", mu,"\n" ) 
    print("-------------------------\n")

    write_occlog_latt(DMFT_loop, Ineq_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, imp_dc, mu, Occup, DMFT_spin_type)
    write_scflog_latt(DMFT_loop, Ineq_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, imp_dc, mu, Occup, DMFT_spin_type, cal_susc,DMFT_solver)
    
    GC.gc()    
    


    #----- Calculating local greenfunction -----#
    print("===== Calculating local green function ... =====\n")
    @time G_loc_iWn=Gloc_gen(H0_k,iWnlist,SelfE_w,mu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind)
    println("\n")

    red_G_loc_iWn =reduce_G_loc(G_loc_iWn,Corr_atom_equiv,Corr_orbital_Ind,imp_ind,DMFT_spin_type)
    write_outfile("GreenF",iWnlist,red_G_loc_iWn,Corr_atom_equiv,DMFT_spin_type,imp_ind)
    GC.gc()
   
    
    #----- calculating local grenn function & hybridization matrix, and print them -----#
    G_loc_iWn_mat = reduce_G_loc_mat(G_loc_iWn,Corr_atom_equiv,Corr_orbital_Ind,DMFT_spin_type)
    InvWiess_iWn_mat = Cal_InvWiessMat_fromGlocMat(iWnlist,G_loc_iWn_mat,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,DMFT_spin_type)
    hyb_iWn_mat = Cal_hybMat_fromInvWiessMat(InvWiess_iWn_mat,imp_ind,iWnlist,mu,H_loc,DMFT_spin_type,Corr_atom_equiv,Corr_orbital_Ind)
    write_outfile_mat("Dmat",iWnlist,hyb_iWn_mat)
    write_outfile_mat("Gfmat",iWnlist,G_loc_iWn_mat)    
    GC.gc()

    #----- Calculating inverse weiss field -----#
    println("===== Calculating weiss field from G_loc_iWn ... =====")
    @time InvWiess_iWn=Cal_InvWiess_fromGloc(iWnlist,G_loc_iWn,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,DMFT_spin_type)
    println("")

    #----- Calculating hybrid function -----#
    println("===== Calculating hybrid function from weiss field ... =====")
    @time hyb_iWn=Cal_hyb_fromInvWiess(InvWiess_iWn,imp_ind,iWnlist,mu,H_loc,DMFT_spin_type,Corr_atom_equiv,Corr_orbital_Ind)
    println("")

    write_outfile("Delta",iWnlist,hyb_iWn,Corr_atom_equiv,DMFT_spin_type,imp_ind) #save hyb function, delta.dat
    GC.gc()

  
    #----- make hybrid function as a form of dictionary -----#
    print("Writing hyb.json ...\n")
    hyb =  hyb_array2dic(Ineq_atom_Ind, hyb_iWn, beta)
    make_imp_folder(Ineq_atom_Ind)

    EDMFTF_write_trans(Trans_mat_C,Corr_atom_equiv,DMFT_spin_type)
    regulated_ShiftE, regulated_Eimp = regulate_Eimp_lev(H_loc, Corr_orbital_Ind, Corr_atom_equiv )    
    EDMFTF_write_PARAMS(Corr_atom_equiv,imp_U,imp_J,imp_dc_n,beta,imp_int_type,EDMFTF_MonteCar_step,EDMFTF_GlobalFlip,EDMFTF_warmup_step,mu,regulated_Eimp,regulated_ShiftE,EDMFTF_tsample,EDMFTF_nom,EDMFTF_PChangeOrder)
    EDMFTF_write_Deltainp(iWnlist,hyb_iWn,Corr_atom_equiv,DMFT_spin_type,imp_ind) # write Delta.inp for EDMFTF CTQMC solver

    
    write_hybjson(Ineq_atom_Ind, hyb)

    write_dc(Ineq_atom_Ind, imp_dc, imp_dc_type, imp_dc_n, imp_totN,  dft_corrN)        
    Save_logfile("dc",Ineq_atom_Ind)            
    
    if DMFT_solver == "ComCTQMC"
        Save_logfile("hyb",Ineq_atom_Ind)
    elseif DMFT_solver == "EDMFTF_ctqmc"
        Save_logfile("Delta",Ineq_atom_Ind)
    end
    GC.gc()


    #----- write input file of CTQMC : params.json -----#
    write_paramsjson(cal_susc,Ineq_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv, imp_ind,imp_int_type, imp_dc,mu, H_loc, regulated_ShiftE, regulated_Eimp
                    ,imp_Measure_time,imp_Thermal_time,beta, Solver_green_Cut, imp_U, imp_J, imp_int_parameterisation, F0, F2, F4, Trans_mat,green_basis, green_legendre_cutoff)

    
    #----- run CTQMC -----#
    run_DMFT(DMFT_solver, DMFT_spin_type, Ineq_atom_Ind, imp_J, imp_ind, imp_int_type, arg_input.Optional["DMFT_Jx"].mpi_prefix, regulated_Eimp)
  
    GC.gc()
        
    if DMFT_solver == "ComCTQMC"
        Save_logfile("params.obs",Ineq_atom_Ind)
    elseif DMFT_solver == "EDMFTF_ctqmc"
        Save_logfile("ctqmc",Ineq_atom_Ind)
        Save_logfile("Sig",Ineq_atom_Ind)
        Save_logfile("Gf",Ineq_atom_Ind)
    end


    ################ START DMFT part [DMFT_loop = 1 ~] ################ 

    for loopN=1:DMFT_loop_N
        global SelfE_w, mulog
    
    
        DMFT_loop=loopN 
        print("\n")
        print("       ============================\n")
        print("          DMFT_loop [",DMFT_loop,"] \n")
        print("       ============================\n")
        print("\n")      
    
        # ----- selfenergy construction by reading obs.json file ------ #
        SelfE_w_new = loadSelfE(Ineq_atom_Ind,DMFT_spin_type,SelfE_w,DMFT_solver)
        write_outfile("SelfE",iWnlist,SelfE_w_new,Corr_atom_equiv,DMFT_spin_type,imp_ind) #write initial selfenergy sig.dat       
        
        if go_bareSelfE(smth_step,Ineq_atom_Ind,DMFT_solver)
            SelfE_w_new_sm = deepcopy(SelfE_w_new)
        else
            SelfE_w_new_sm =SelfE_smoothing(Corr_atom_equiv,Tem,iWnlist,SelfE_w_new); #save smoothing version of SelfE
            write_outfile("SelfE_sm",iWnlist,SelfE_w_new_sm,Corr_atom_equiv,DMFT_spin_type,imp_ind) #write initial selfenergy sig.dat
        end
        #write_outfile("SelfE_sm",iWnlist,SelfE_w_new_sm,Corr_atom_equiv,DMFT_spin_type,imp_ind) #write initial selfenergy sig.dat
        
        
        imp_occ= loadOccup_from_obs(Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,DMFT_solver)
        if ( (cal_susc) || (DMFT_solver == "EDMFTF_ctqmc") ) 
            susc_N0, susc_S0 = load_suscw0(DMFT_loop,Ineq_atom_Ind,DMFT_solver)
        else
            susc_N0 = 0.0
            susc_S0 = 0.0
        end
        write_occlog_imp(Ineq_atom_Ind, imp_ind, DMFT_spin_type, Corr_ineq_orbital_Num, imp_occ, DMFT_solver)
        write_scflog_imp(Ineq_atom_Ind, imp_ind, DMFT_spin_type, Corr_ineq_orbital_Num, imp_occ, susc_N0, susc_S0, cal_susc, DMFT_solver)
        if ( (cal_susc) || (DMFT_solver == "EDMFTF_ctqmc") ) 
            write_susclog_fromObs(DMFT_loop,Ineq_atom_Ind,DMFT_solver)
        end
        #Save_logfile("params.obs",Ineq_atom_Ind)
    
        SelfE_w = Mix_SelfE_(SelfE_w,SelfE_w_new_sm,Mix_selfE)  # Mixing selfE, 
        write_outfile("SelfE_mix",iWnlist,SelfE_w,Corr_atom_equiv,DMFT_spin_type,imp_ind) #write initial selfenergy sig.dat
    
        #if go_bareSelfE(smth_step,Ineq_atom_Ind)
        #    SelfE_w_sm = deepcopy(SelfE_w)
        #else
        #    SelfE_w_sm =SelfE_smoothing(Corr_atom_equiv,Tem,iWnlist,SelfE_w); #save smoothing version of SelfE
        #end
    
        DFTforge.pwork(init_variables_SelfE_w,SelfE_w); #passing variable SelfE to worker
        #write_outfile("SelfE_sm",iWnlist,SelfE_w_sm,Corr_atom_equiv,DMFT_spin_type,imp_ind) #write initial selfenergy sig.dat

        
        imp_totN = loadOccup_from_obs_scalar(Ineq_atom_Ind,DMFT_solver);
        imp_dc = dc_change(imp_totN, imp_dc_n, dft_corrN, imp_dc_type, imp_U, imp_J) # calculate double counting
        DFTforge.pwork(passing_dc, imp_dc) # passing double counting
        write_dc(Ineq_atom_Ind, imp_dc, imp_dc_type, imp_dc_n, imp_totN, dft_corrN)        
        Save_logfile("dc",Ineq_atom_Ind)                
        
    
        GC.gc()
    
    
        println("")
    
    
        #----- Find mu with new SelfE-----#
    
        if compute_EF
            init_mu = mu_guess(DMFT_loop, mulog, DMFT_spin_type)
            @time begin
                print("          ===== Finding mu ... =====\n")
                mu, Occup, H_loc =Find_mu(iWnlist,H0_k,SelfE_w,targetN,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv,DMFT_spin_type,init_mu)
                println("")
            end
        else
            @time begin
                mu = 0.0
                Occup, H_loc=Cal_Occupation_givenMu_WOdc(iWnlist,H0_k,SelfE_w,mu,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv);
            end
        end
        
    
        GC.gc()
    
        push!(mulog,mu)
    
        write_occlog_latt(DMFT_loop, Ineq_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, imp_dc, mu, Occup, DMFT_spin_type)
        write_scflog_latt(DMFT_loop, Ineq_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, imp_dc, mu, Occup, DMFT_spin_type,cal_susc,DMFT_solver)

    
        #----- Calculating local greenfunction -----#
        print("===== Calculating local green function ... =====\n")
        @time G_loc_iWn=Gloc_gen(H0_k,iWnlist,SelfE_w,mu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind)
        println("")
    
    
        red_G_loc_iWn =reduce_G_loc(G_loc_iWn,Corr_atom_equiv,Corr_orbital_Ind,imp_ind,DMFT_spin_type)
        write_outfile("GreenF",iWnlist,red_G_loc_iWn,Corr_atom_equiv,DMFT_spin_type,imp_ind)
        GC.gc()

   
        #----- calculating local grenn function & hybridization matrix, and print them -----#
        G_loc_iWn_mat = reduce_G_loc_mat(G_loc_iWn,Corr_atom_equiv,Corr_orbital_Ind,DMFT_spin_type)
        InvWiess_iWn_mat = Cal_InvWiessMat_fromGlocMat(iWnlist,G_loc_iWn_mat,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,DMFT_spin_type)
        hyb_iWn_mat = Cal_hybMat_fromInvWiessMat(InvWiess_iWn_mat,imp_ind,iWnlist,mu,H_loc,DMFT_spin_type,Corr_atom_equiv,Corr_orbital_Ind)
        write_outfile_mat("Dmat",iWnlist,hyb_iWn_mat)
        write_outfile_mat("Gfmat",iWnlist,G_loc_iWn_mat)    
        GC.gc()
        
        
    
        #----- Calculating inverse weiss field -----#
        println("===== Calculating weiss field from G_loc_iWn ... =====")
        @time InvWiess_iWn=Cal_InvWiess_fromGloc(iWnlist, G_loc_iWn,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,DMFT_spin_type)
        println("")

        #----- Calculating hybrid function -----#
        println("===== Calculating hybrid function from weiss field ... =====")
        @time hyb_iWn=Cal_hyb_fromInvWiess(InvWiess_iWn,imp_ind,iWnlist,mu,H_loc,DMFT_spin_type,Corr_atom_equiv,Corr_orbital_Ind)
        println("")    
        GC.gc()
    
        write_outfile("Delta",iWnlist,hyb_iWn,Corr_atom_equiv,DMFT_spin_type,imp_ind)
    
        #----- make hybrid function as a form of dictionary -----#
        print("Writing hyb.json ...\n")
        hyb =  hyb_array2dic(Ineq_atom_Ind, hyb_iWn, beta)

        make_imp_folder(Ineq_atom_Ind)
        write_hybjson(Ineq_atom_Ind, hyb)

   
        if DMFT_solver == "ComCTQMC"
            Save_logfile("hyb",Ineq_atom_Ind)
        elseif DMFT_solver == "EDMFTF_ctqmc"
            Save_logfile("Delta",Ineq_atom_Ind)
        end
        #----- write input file of CTQMC with changed mu : params.json -----#
        change_mu_paramsjson(Ineq_atom_Ind,imp_dc,mu,regulated_ShiftE)
        EDMFTF_write_PARAMS(Corr_atom_equiv,imp_U,imp_J,imp_dc_n,beta,imp_int_type,EDMFTF_MonteCar_step,EDMFTF_GlobalFlip,EDMFTF_warmup_step,mu,regulated_Eimp,regulated_ShiftE,EDMFTF_tsample,EDMFTF_nom,EDMFTF_PChangeOrder)
        EDMFTF_write_Deltainp(iWnlist,hyb_iWn,Corr_atom_equiv,DMFT_spin_type,imp_ind) # write Delta.inp for EDMFTF CTQMC solver
    
    
        #----- run CTQMC -----#
        run_DMFT(DMFT_solver, DMFT_spin_type, Ineq_atom_Ind, imp_J, imp_ind, imp_int_type, arg_input.Optional["DMFT_Jx"].mpi_prefix, regulated_Eimp)
        GC.gc()
        if DMFT_solver == "ComCTQMC"
            Save_logfile("params.obs",Ineq_atom_Ind)
        elseif DMFT_solver == "EDMFTF_ctqmc"
            Save_logfile("ctqmc",Ineq_atom_Ind)
            Save_logfile("Sig",Ineq_atom_Ind)
            Save_logfile("Gf",Ineq_atom_Ind)
        end
    
    
    
    end


end




################################################################################
######################### Jx calculation part w/ SelfE  ########################
################################################################################


if ( Calculation_mode == "Jx-DMFT" )

    println("#####################################################")
    println(" Calculation type : Jx calculation with self-energy ")
    println("#####################################################\n\n\n")

    
    ## DMFT part input
    @time DFTforge.pwork(init_variables_Jx,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind,beta, Corr_ineq_orbital_Num, delta); # passing basic variable for DMFT ++ , here beta_for_occ means for initial occ
    @time DFTforge.pwork(passing_dc, imp_dc)
    
    ################ GRID GENERATION PART ##################

    atom_orbitals=atom_orb_list(hamiltonian_info);
    Total_orb_Num=sum(hamiltonian_info.scf_r.Total_NumOrbs)

    Rlist=RPoint_gen(R_grid_num);
    Klist=kPoint_gen(K_grid_num);
    iWnlist=Matsubara_gen(iW_grid_cut,beta);
    
    if !isfile("Nele.dat")
        error("!!! No file : Nele.dat, please check if you conducted DMFT loop or not !!!")
    else
        @time DFTforge.pwork(init_variables_grid_DMFT,iWnlist,Rlist,Klist); # passing grid variable to worker
    end
    
    @time DFTforge.pwork(init_variables_totorb,Total_orb_Num); # passing total number of orbital to workerspace


    

    
    
    ################  non-interacting hamiltonian construction ################ 
    @time H_k=nonInt_H_k(Klist,hamiltonian_info,Calculation_mode,DMFT_spin_type)
    mu=hamiltonian_info.scf_r.ChemP
    DCM =  dcMat_noglob(H_k,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,imp_ind)
    H0_k = DC_shift(H_k,DCM)
    @time DFTforge.pwork(init_variables_H_k,H0_k);  #passing variable H_k to worker
    println("")



    ################ START DMFT part [DMFT_loop = 0] ################ 

    DMFT_loop= 0

    # ----- initial selfenergy construction Self =0.0+0.0im ------ #
    SelfE_w = init_selfE_dcbias(imp_dc,Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,iWnlist,0.0)
    

    targetN=parse(Float64, read_Nele())[1]


    #----- passing grid inform & selfE & input var -----#
    DFTforge.pwork(init_variables_grid_DMFT,iWnlist,Rlist,Klist); # passing grid variable to worker
    DFTforge.pwork(init_variables_SelfE_w,SelfE_w); #passing variable SelfE to worker
    DFTforge.pwork(init_variables_DMFT,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind, Corr_ineq_orbital_Num, beta, delta); # passing basic variable for DMFT ++


    GC.gc()

    
    
    dft_corrN = deepcopy(imp_dc_n)
    if (imp_dc_type == "FLL-DFT" || imp_dc_type == "fll-dft" || imp_dc_type == "Fll-dft")
        if !isfile("DFT_CorrNele.dat")
            error("!!! No file : DFT_CorrNele.dat, please check if you conducted FLL-DFT double counitng !!!")
        else
            dft_corrN= read_DFT_corrNele()
        end
    end    
    

    
    Trans_mat=init_trans_mat(Corr_atom_Ind, Corr_orbital_Ind)
    if !isfile("transform_Mat.dat")
        error("There is no file 'transform_Mat.dat', please check that you need basis_tranform or not")
    else
        Trans_mat = read_transforMat()
        Trans_mat_A =Tvec_print_from_reading(H0_k, Trans_mat, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)
    end
            
    
    
    @time H0_k = Trans_H_k(H0_k,Trans_mat_A)
    DFTforge.pwork(init_variables_H_k,H0_k);  #passing variable H_k to worker
 
        
    
    @time H_loc= local_H0(H0_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind)
    #if "consider_degeneracy" turned on, the generate states would be treated equivalently 
    if consider_degeneracy == true
        imp_ind = imp_ind_trans_by_degeneracy(imp_ind,H_loc,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind)
        Corr_ineq_orbital_Num = cal_corr_ineq_orbital_N(imp_ind)
        DFTforge.pwork(init_variables_DMFT,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind, Corr_ineq_orbital_Num, beta, delta); # passing basic variable for DMFT ++
    end
 

    # ----- change double counintg if imp_dc_type = "FLL-DFT"  ------ #
    imp_totN = deepcopy(imp_dc_n)
    imp_dc = dc_change(imp_totN, imp_dc_n, dft_corrN, imp_dc_type, imp_U, imp_J) # calculate double counting
    DFTforge.pwork(passing_dc, imp_dc) # passing double counting
        
    # ----- initial selfenergy construction with DC+bias ------ #
    SelfE_w = init_selfE_dcbias(imp_dc,Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,iWnlist,init_bias)    
    DFTforge.pwork(init_variables_SelfE_w,SelfE_w); #passing variable SelfE to worker

    println("--------------------------------------")
    println("imp_dc : ", imp_dc)
    println("--------------------------------------")


    GC.gc()    


    #----- check if do restart or not -----#
    imp_totN = deepcopy(imp_dc_n)
    restart_wSelfE = Check_restart_wSelfE(Ineq_atom_Ind, DMFT_solver)
    if  (restart_wSelfE) # if in all imp_x folder there is output file params.obs.json, then just read self energy 
        SelfE_w_new = loadSelfE(Ineq_atom_Ind,DMFT_spin_type,SelfE_w,DMFT_solver)

        imp_totN = loadOccup_from_obs_scalar(Ineq_atom_Ind,DMFT_solver);
        imp_dc = dc_change(imp_totN, imp_dc_n, dft_corrN, imp_dc_type, imp_U, imp_J) # calculate double counting
        DFTforge.pwork(passing_dc, imp_dc) # passing double counting        
        
        
        if go_bareSelfE(smth_step,Ineq_atom_Ind,DMFT_solver)
            SelfE_w_sm = deepcopy(SelfE_w_new)
        else
            SelfE_w_sm =SelfE_smoothing(Corr_atom_equiv,Tem,iWnlist,SelfE_w_new);
        end
        
        DFTforge.pwork(init_variables_SelfE_w,SelfE_w_sm); #passing the smoothing SelfE to worker
    else
        if DMFT_solver == "ComCTQMC"
            error("!!! No file : imp_*/params.obs.json, please check if you conducted DMFT loop or not !!!")
        elseif DMFT_solver == "EDMFTF_ctqmc"
            error("!!! No file : imp_*/Sig.out, please check if you conducted DMFT loop or not !!!")
        end
    end

    SelfE_w_new = nothing
    GC.gc()


    #----- Find mu -----#

    if isfile("scf.log")
        mulog = readmu_fromscflog()
    else
        error("!!! No file : scf.log, please check if you conducted DMFT loop or not !!!")
    end
    
    init_mu = mu_guess(10, mulog, DMFT_spin_type)
    @time begin
        print("          ===== Finding mu ... =====\n")
        mu, Occup, H_loc =Find_mu(iWnlist,H0_k,SelfE_w_sm,targetN,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv,DMFT_spin_type,init_mu)
        println("")
    end
    

    print("-------------------------\n")
    print("Initial latt. mu :", mu,"\n") 
    print("-------------------------\n")
    
    
    println("Up-spin occup :", diag(real(Occup[:,:,1])))
    println("")
    println("dn-spin occup :", diag(real(Occup[:,:,2])))
    println("")
    println("****************************************")
    println("Total charge #:", sum(diag(real(Occup[:,:,1]+Occup[:,:,2])) ) )
    println("****************************************")    

    
    


    wdivind = iWnlist_div(iWnlist,10)
    
    #------ initialize J_R  -------#
    JR,JR_orb,distance=init_JR(R_grid_num, atom12_list, Corr_atom_Ind, Corr_orbital_Ind)
    ################### J Calculation ##############
    JR,JR_orb,distance=Jx_calculation(JR, JR_orb,distance,hamiltonian_info,H0_k,mag_order,R_grid_num,Rlist,Klist,iWnlist,wdivind,atom12_list,SelfE_w_sm,mu,Corr_atom_Ind,Corr_orbital_Ind)
    
    
    sorted_dist,sorted_cell_ind,sorted_R_grid_ind,sorted_JR,sorted_JR_orb=sort_dist(distance,JR,JR_orb,atom12_list, R_grid_num)

    
    
    @printf "=================================== AFM order (-) | FM order (+) ===================================\n"
    for jj=1:size(atom12_list)[1]
        @printf "------------------------------------ Site ind. in cell : %s ------------------------------------\n" atom12_list[jj]
        for ii=1:10
    
            @printf "       ||  Distance: %9.4f (Ang)  |  J: %9.4f (meV)  |  Cell ind.: %14.12s  || \n" sorted_dist[jj][ii] real(sorted_JR[jj][ii]) sorted_cell_ind[jj][ii] 
        end
    end
    

    J_print(atom12_list,sorted_cell_ind,sorted_R_grid_ind, sorted_dist, sorted_JR, sorted_JR_orb)    
end    
    
    
    

#################################################################################
######################### Jx calculation part w/o SelfE  ########################
#################################################################################


if ( Calculation_mode == "Jx0" )

    println("#####################################################")
    println("  Calculation type : Jx calculation w/o self-energy  ")
    println("#####################################################\n\n\n")
    
    
    
    
    ################ GRID GENERATION PART ##################

    atom_orbitals=atom_orb_list(hamiltonian_info);
    Total_orb_Num=sum(hamiltonian_info.scf_r.Total_NumOrbs)

    Rlist=RPoint_gen(R_grid_num);
    Klist=kPoint_gen(K_grid_num);    
    iWnlist=Matsubara_gen(iW_grid_cut,beta);

    
    #@time DFTforge.pwork(init_test,1,2,H_k)
    @time DFTforge.pwork(init_variables_grid_DMFT,iWnlist,Rlist,Klist); # passing grid variable to worker for determining initial occ
    @time DFTforge.pwork(init_variables_totorb,Total_orb_Num); # passing total number of orbital to workerspace

    
    ################  non-interacting hamiltonian construction ################     
    @time H_k = nonInt_H_k(Klist,hamiltonian_info,Calculation_mode,DMFT_spin_type);
    DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
    
    
    
    #------ Fourier transform k to R -----#
    println("Fourier transform of H_k : k -> R ...")
    H_R=zeros(ComplexF64,size(H_k))
    @time FFTW_k2r_paralle_Hk(Klist,Rlist,H_k,H_R)
    #@time fourier_k2r_paralle(Klist,Rlist,H_k,H_R);
    println("")

    

    
    ################ Set SelfE = 0.0 ######################
    println("Setting the selfE = 0.0...")
    SelfE_w = init_selfE_zero(Ineq_atom_Ind,2,Corr_ineq_orbital_Num,hamiltonian_info,iWnlist)
    @time DFTforge.pwork(init_variables_SelfE_w,SelfE_w); #passing variable SelfE to worker for determining initial occ
    println("")

    GC.gc()


    ################  Calculating Occupation ##############
    testmu=hamiltonian_info.scf_r.ChemP
    println("Calculating initial Occupation of wannier functions ...")
    @time Occup, H_loc=Cal_Occupation_givenMu_WOdc(iWnlist,H_k,SelfE_w,testmu,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv);
    println("")    
    
    
    println("Up-spin occup :", diag(real(Occup[:,:,1])))
    println("")
    println("dn-spin occup :", diag(real(Occup[:,:,2])))
    println("")
    println("****************************************")
    println("Total charge #:", sum(diag(real(Occup[:,:,1]+Occup[:,:,2])) ) )
    println("****************************************")
    
    
    

    wdivind = iWnlist_div(iWnlist,20)
    
    #------ initialize J_R  -------#
    JR,JR_orb,distance=init_JR(R_grid_num, atom12_list, Corr_atom_Ind, Corr_orbital_Ind)
    ################### J Calculation ##############
    JR,JR_orb,distance=Jx_calculation(JR, JR_orb,distance,hamiltonian_info,H_k,mag_order,R_grid_num,Rlist,Klist,iWnlist,wdivind,atom12_list,SelfE_w,testmu,Corr_atom_Ind,Corr_orbital_Ind)
    
    
    sorted_dist,sorted_cell_ind,sorted_R_grid_ind,sorted_JR,sorted_JR_orb=sort_dist(distance,JR,JR_orb,atom12_list, R_grid_num)

    
    
    @printf "=================================== AFM order (-) | FM order (+) ===================================\n"
    for jj=1:size(atom12_list)[1]
        @printf "------------------------------------ Site ind. in cell : %s ------------------------------------\n" atom12_list[jj]
        for ii=1:10
    
            @printf "       ||  Distance: %9.4f (Ang)  |  J: %9.4f (meV)  |  Cell ind.: %14.12s  || \n" sorted_dist[jj][ii] real(sorted_JR[jj][ii]) sorted_cell_ind[jj][ii] 
        end
    end
    

    J_print(atom12_list,sorted_cell_ind,sorted_R_grid_ind, sorted_dist, sorted_JR, sorted_JR_orb)    
end

###################################################
### functions for magnon dispersion calculation ###
###################################################



function relative_position(frac_coord)
    anum = size(frac_coord)[1]
    rel_pos_frac = zeros(anum, anum, 3)
    
    for i=1:anum
        for j=1:anum
            rel_pos_frac[i,j,:]=frac_coord[i,:].-frac_coord[j,:]
        end
    end
    
    return rel_pos_frac
end


function check_isradata(atom12_list,Corr_atom_Ind)
    
    check_Box =[]
    nexist_file =[]
    readfilelist = []
    for i in Corr_atom_Ind
        for j in Corr_atom_Ind
                
            filen = string("Jx_raw_",i,"_",j,".dat")
            push!(readfilelist, filen)
            push!(check_Box, isfile(filen))
            if !isfile(filen)
                push!(nexist_file, filen) 
            end
        end
    end
    
    if !prod(check_Box)
        error("*** Error : You should check if the following files exist or not : ",nexist_file, "***")
    end
     
    return readfilelist
end


function Jmat_construct(readfilelist,Corr_atom_Ind,R_grid_num)
    rubyscan(pat, str) = (m.match for m in eachmatch(pat, str)) 
    pat = r"[+-]?\d+\.?\d*"
    rawnum = []
    
    for filename in readfilelist
        re = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$";
        open(filename, "r") do f 
            Nstring=readlines(f)
            datalinenum = 0
            for i =1:size(Nstring)[1]
                if( occursin(re,split(Nstring[i])[1]) ) 
                    datalinenum +=1    
                end
            end
            push!(rawnum,datalinenum)
        end
    end
    
    checkequal = rawnum./maximum(rawnum)
    if prod(checkequal.== 1)
        Jmat = zeros(R_grid_num[1],R_grid_num[2],R_grid_num[3], size(Corr_atom_Ind)[1], size(Corr_atom_Ind)[1])
        cellvec = zeros(R_grid_num[1],R_grid_num[2],R_grid_num[3],3)
        dist = zeros(R_grid_num[1],R_grid_num[2],R_grid_num[3], size(Corr_atom_Ind)[1], size(Corr_atom_Ind)[1])
    else
        error("Error : The number of raw of files is not equal !!!")
    end
    
    
    
    for filename in readfilelist
        re = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$";
        matind = Int.(parse.(Float64, rubyscan(pat, filename)))
        open(filename, "r") do f 
            Nstring=readlines(f)
            datalinenum = 0
            for i =1:size(Nstring)[1]
                if( occursin(re,split(Nstring[i])[1]) ) 
                    cellind =  parse.(Int64,split(Nstring[i])[1:3])
                    Rind = parse.(Int64,split(Nstring[i])[4:6])
                    matind2 = parse.(Int64,split(Nstring[i])[7:8])
                    #println(matind, ", ", matind2)
                    dist_dum = parse(Float64,split(Nstring[i])[9])
                    Jval = parse(Float64,split(Nstring[i])[10])
                    
                    if matind == matind2
                    #if true
                        aind1 = findall(x->x == matind[1], Corr_atom_Ind)
                        aind2 = findall(x->x == matind[2], Corr_atom_Ind)
                        
                        
                        dist[Rind[1]+1,Rind[2]+1,Rind[3]+1,aind1[1],aind2[1]] = dist_dum
                        #dist[Rind[1]+1,Rind[2]+1,Rind[3]+1,aind2[1],aind1[1]] = dist_dum
                        
                        Jmat[Rind[1]+1,Rind[2]+1,Rind[3]+1,aind1[1],aind2[1]] = Jval
                        #Jmat[Rind[1]+1,Rind[2]+1,Rind[3]+1,aind2[1],aind1[1]] = Jval
                        cellvec[Rind[1]+1,Rind[2]+1,Rind[3]+1,:] = cellind
            
                    else
                        error("Error : Jmat index is mismatched !!!")
                    end
                end
            end
            push!(rawnum,datalinenum)
        end
    end    
    
    ordered_dist = []
    for i =1:size(Corr_atom_Ind)[1]
        push!(ordered_dist, sort(unique(trunc.(dist[:,:,:,i,:], digits = 5))))
    end
    
    
    return ordered_dist, dist, cellvec, Jmat
end



function read_kpath(filename)
    kpath=[]
    kind=[]
    open(filename, "r") do f 
        Nstring=readlines(f)
        for i =1:size(Nstring)[1]
            println(split(Nstring[i]))
            push!(kind,split(Nstring[i])[1])
            push!(kpath,[parse(Float64,split(Nstring[i])[2]) parse(Float64,split(Nstring[i])[3]) parse(Float64,split(Nstring[i])[4]) ])
            push!(kind,split(Nstring[i])[5])
            push!(kpath,[parse(Float64,split(Nstring[i])[6]) parse(Float64,split(Nstring[i])[7]) parse(Float64,split(Nstring[i])[8]) ])
        end
    end
    return kind, kpath
end
    



function band_kpts2(kpath,totgrid)
    
    kpt = zeros(totgrid,3)
    xleng =[]
    lattind = []
    push!(lattind, 1)
    push!(xleng,0.0)
    totlength=0
    
    kind = 1
    kpt[1,:]=kpath[1]
    dumleng = 0
    
    seg_num = Int(size(kpath)[1]/2)
    
    for i =1:seg_num
        length_ = ((kpath[i*2][1]-kpath[i*2-1][1])^2+(kpath[i*2][2]-kpath[i*2-1][2])^2+(kpath[i*2][3]-kpath[i*2-1][3])^2)^(1/2)
        dumleng += length_
        push!(xleng,dumleng)
        totlength += length_
    end
    xpt = range(0, totlength, length=totgrid)
    xpt_group = zeros(Int64,size(xpt)[1])
    for i = 1:size(xpt)[1]
        for j =1:seg_num
            if(xleng[j]<xpt[i]<=xleng[j+1])
                xpt_group[i]=j
            end
        end
    end
    
    for i=1:seg_num
        ind_list = findall(x->x==i, xpt_group)
     
        grid = size(ind_list)[1]
        k1=range(kpath[i*2-1][1], stop=kpath[i*2][1], length =grid+1)
        k2=range(kpath[i*2-1][2], stop=kpath[i*2][2], length =grid+1)
        k3=range(kpath[i*2-1][3], stop=kpath[i*2][3], length =grid+1)
        
        
        for j=2:grid+1
        
            kind=kind+1
            kpt[kind,:]=[k1[j],k2[j],k3[j]]
            if (j==grid+1)
                push!(lattind, kind) 
            end
        end
    end
    
    return kpt, xpt,lattind
end




function Cal_MagMat(qmode, Lattconst, qpt_frac, cell_vec, r_vec, rel_pos_frac, R_grid_num, cellvec_fac, Corr_atom_Ind, Corr_atom_equiv, Moment, Neighbors_cut, ordered_dist, dist, Jmat)
    Orbind_num = [size(Jmat)[4],  size(Jmat)[5]]
    Magmat_ = zeros(ComplexF64, size(qpt_frac)[1], Orbind_num[1], Orbind_num[2])
    considered_Neigh = zeros(Int64, size(Corr_atom_Ind)[1], Neighbors_cut)
  
    sign = Corr_atom_equiv./(abs.(Corr_atom_equiv))

    Magnon_disp = zeros(ComplexF64, size(qpt_frac)[1], Orbind_num[1])
    
    for qind =1:size(qpt_frac)[1]
        qvec = qpt_frac[qind,1]*r_vec[1,:] + qpt_frac[qind,2]*r_vec[2,:] + qpt_frac[qind,3]*r_vec[3,:]
        
        for (o1,o1val) in enumerate(Corr_atom_Ind)
            for (o2,o2val) in enumerate(Corr_atom_Ind)
                for R1 = 1:R_grid_num[1]
                    for R2 = 1:R_grid_num[2]
                        for R3 = 1:R_grid_num[3]
                            if Neighbors_cut == 300
                                Neighbors_cut = size(ordered_dist[o1])[1]-1
                            end
 
                            if ( dist[R1,R2,R3,o1,o2] <= ordered_dist[o1][Neighbors_cut+1]+0.00001 )
                                for n = 1:Neighbors_cut
                                    if ( abs(dist[R1,R2,R3,o1,o2]-ordered_dist[o1][n+1]) <0.00001 && (qind ==1) )
                                        considered_Neigh[o1,n]+=1
                                    end
                                end
                                cvec = cellvec_fac[R1,R2,R3,1]*cell_vec[1,:] + cellvec_fac[R1,R2,R3,2]*cell_vec[2,:] + cellvec_fac[R1,R2,R3,3]*cell_vec[3,:]
                                
                                Magmat_[qind, o1,o1] += (4.0/Moment[o1])*(Jmat[R1,R2,R3,o1,o2])*sign[o2]  
                                #Magmat_[qind, o1,o2] -= (4.0/Moment[o1])*(Jmat[R1,R2,R3,o1,o2])*sign[o1]*exp( 2*pi*dot(qpoint[qind,:], cellvec[R1,R2,R3,:])*im )
                                if qmode == 1
                                    Magmat_[qind, o1,o2] -= (4.0/Moment[o1])*(Jmat[R1,R2,R3,o1,o2])*sign[o1]*exp( dot(qvec, cvec)*im )
                                elseif qmode == 2 
                                    Magmat_[qind, o1,o2] -= (4.0/Moment[o1])*(Jmat[R1,R2,R3,o1,o2])*sign[o1]*exp( 2*pi/Lattconst*dot(qpt_frac[qind,:], cvec)*im )
                                end
                                
                            end

                            
                        end
                    end
                end
            end
        end
        
        Magnon_disp[qind,:] = eigvals(Magmat_[qind,:,:])
    end
            
    
    return considered_Neigh, Magmat_, Magnon_disp
end

function reciprocal_vec(cell_vec)
    r_vec = deepcopy(cell_vec)
    for i=1:size(r_vec)[1]
        cvec1=cross(cell_vec[2,:], cell_vec[3,:])
        r_vec[1,:] =2*pi*cvec1/dot(cell_vec[1,:],cvec1)
        
        cvec2=cross(cell_vec[3,:], cell_vec[1,:])
        r_vec[2,:] =2*pi*cvec2/dot(cell_vec[2,:],cvec2)
        
        cvec3=cross(cell_vec[1,:], cell_vec[2,:])
        r_vec[3,:] =2*pi*cvec3/dot(cell_vec[3,:],cvec3)
    end
    
    return r_vec 
end

function print_c_r_vec(cell_vec, r_vec)
   
    println("Unit cell vectors : ")
    for i =1:size(cell_vec)[1]
        println("  ", cell_vec[i,:])
    end

    println("")
    println("Reciprocal vectors : ")
    for i =1:size(cell_vec)[1]
        println("  ", r_vec[i,:])
    end
    println("")
   
end


function loadmomentsize(filename, Corr_atom_equiv)
        
    lattmom = []
    impmom =[]
    Ineq_atom_ind = unique(abs.(Corr_atom_equiv))
    
    Iterline_ind = []
    complete_sec_Iterline_ind =[]
    
    latt_mom =zeros(size(Corr_atom_equiv)[1])
    imp_mom =zeros(size(Corr_atom_equiv)[1])
    
    
    open(filename, "r") do f
        Nstring=readlines(f)
        for i=1:size(Nstring)[1]
            if !(split(Nstring[i]) ==[])
                if (split(Nstring[i])[1] == "Iteration" && i+(2*size(Ineq_atom_ind)[1])+2 <= size(Nstring)[1] )
                    push!(Iterline_ind, i)
                end
            end
        end
        
        for i in Iterline_ind
            lattnum = 0
            impnum = 0
            for j = i:i+(2*size(Ineq_atom_ind)[1])+2
                if split(Nstring[j])[1] == "latt."
                    lattnum += 1
                end
                
                if split(Nstring[j])[1] == "imp."
                    impnum += 1
                end
                
                if (impnum == size(Ineq_atom_ind)[1] && impnum == size(Ineq_atom_ind)[1])
                    push!(complete_sec_Iterline_ind, i)
                end
                
            end
        end
        
        for i = Iterline_ind[end]:Iterline_ind[end]+(2*size(Ineq_atom_ind)[1])+2
                
            if split(Nstring[i])[1] == "latt."
                mind = findall(x->x=="m",split(Nstring[i]))[1]
                push!(lattmom, parse(Float64, split(split(Nstring[i])[mind+2], ",")[1]) )
            end
            
                
            if split(Nstring[i])[1] == "imp."
                mind = findall(x->x=="m",split(Nstring[i]))[1]
                push!(impmom, parse(Float64, split(split(Nstring[i])[mind+2], ",")[1]) )
            end
            
        end
        
    end
    
    for (i,v) in enumerate(Ineq_atom_ind)
        ind = findall(x->x == v, abs.(Corr_atom_equiv))
        latt_mom[ind] .= lattmom[i]
        imp_mom[ind] .= impmom[i]
    end
        
    mom_size =0.5*(latt_mom.+imp_mom)
    
    
    return mom_size
end  



function print_magnon(xpt, Magnon_disp, lattind, kindex)
    max_ene_arr =[]
    open("magnon.dat", "w") do io
        @printf(io,"    %14s", "kpt")
        for i=1:size(Magnon_disp)[2]
            push!(max_ene_arr, maximum(abs.(real(Magnon_disp[:,i]))) )
            tag = string("disp_",i,)
            if i !=size(Magnon_disp)[2]
                @printf(io,"    %14s", tag)
            else
                @printf(io,"    %14s\n", tag)
            end            
        end
        
        for i=1:size(Magnon_disp)[1]
            @printf(io,"    %14.8f", xpt[i])
            for j=1:size(Magnon_disp)[2]
                tag = string("disp_",j,)
                if j !=size(Magnon_disp)[2]
                    @printf(io,"    %14.8f", real(Magnon_disp[i,j]))
                else
                    @printf(io,"    %14.8f\n", real(Magnon_disp[i,j]))
                end            
            end           
        end
    end

    
    max_ene = maximum(max_ene_arr)
    open("vline.dat", "w") do io
        for i =2:(size(lattind)[1]-1)
            @printf(io,"%8.4f  %8.4f\n",xpt[lattind[i]],0.0)
            @printf(io,"%8.4f  %8.4f\n\n",xpt[lattind[i]],max_ene*1.3)
        end
    end
        

    

    open("magnon.gnu", "w") do io
        @printf(io,"set xra [%8.4f: %8.4f]\n",xpt[1],xpt[end])
        @printf(io,"set yra [%8.4f: %8.4f]\n",0.0,max_ene*1.3)
        @printf(io, "set ylabel \"Energy (meV)\" \n")
        @printf(io,"set xtics (")
        for (i,v) in enumerate(lattind)
            if i == 1
                str = kindex[i]
            elseif i == size(lattind)[1]
                str = kindex[2*(i-1)]
            else
                if kindex[2*(i-1)] == kindex[2*(i-1)+1]
                    str = kindex[2*(i-1)]
                else
                    str = string(kindex[2*(i-1)],"|",kindex[2*(i-1)+1])
                end
            end
            
            @printf(io,"\"%s\" ",str)
            if i != size(lattind)[1]
                @printf(io," %8.4f,",xpt[v])
            else
                @printf(io," %8.4f)\n",xpt[v])
            end
   
            
        end
        
        for i=1:size(Magnon_disp)[2]
            if i==1
                @printf(io,"plot \"magnon.dat\" u 1:%i w l \n",i+1)
            else
                @printf(io,"replot \"magnon.dat\" u 1:%i w l \n",i+1)
            end            
        end
        @printf(io,"replot \"vline.dat\" u 1:2 w l lc 0 \n")
        
        @printf(io,"pause -1\n")        

    end
end

######################################################################
######################### Magnono dispersion  ########################
######################################################################


if ( Calculation_mode == "Magnon" )
    println("    #############################################")
    println("       Calculation type : Spin-wave calculation  ")
    println("    #############################################\n\n\n")  
    
    println("               DFT type : ", hamiltonian_info.dfttype) 
    println("")
    println("!!! If DFT type == Wannier90, Unit cell vector  !!!")
    println("!!! and atom position are read from *.win and   !!!")
    println("!!! *.toml respectively                         !!!")
    println("")
    
    println("qmode : ", qmode)
    println("")

    if (qmode != 1)
        println("Please input the Lattice constant (L), to make grids 2pi/L*qpoint_factors ...")
        Lattconst_str = readline()
        Lattconst = parse(Float64, Lattconst_str)
    else
        Lattconst = 1
    end

    cell_vec = deepcopy(hamiltonian_info.scf_r.tv)
    r_vec = reciprocal_vec(cell_vec)
   
    print_c_r_vec(cell_vec, r_vec)
    
    Moment_size = loadmomentsize("occ.log", Corr_atom_equiv)  #magnetic moment size, unit :Bohr magneton
    println("Magnetic moment size :", Moment_size, " (unit : Bohr magneton)")
    println("")
    frac_coord = hamiltonian_info.scf_r.Gxyz*inv(hamiltonian_info.scf_r.tv)
    rel_pos_frac = relative_position(frac_coord)
    readfilelist = check_isradata(atom12_list,Corr_atom_Ind)
    println("Loading Jmat from following file list ...")
    println("    :", readfilelist)
    println("")
    ordered_dist, dist, cellvec, Jmat = Jmat_construct(readfilelist,Corr_atom_Ind,R_grid_num)
    
    tot_grid = 1000   # grid # between two high symmetry points
    println("Magnon k-path w/ grid num :",tot_grid)
    kindex, kpath= read_kpath("kpath")
    kpt_frac, xpt,lattind= band_kpts2(kpath,tot_grid)
    
    
    considered_Neigh, Magmat, Magnon_disp =Cal_MagMat(qmode, Lattconst, kpt_frac, cell_vec, r_vec, rel_pos_frac, R_grid_num, cellvec, Corr_atom_Ind, Corr_atom_equiv, Moment_size, Neighbors_cut, ordered_dist, dist, Jmat)
    println()
    
    println("Considered Neighbors : ")
    for i=1:size(considered_Neigh)[1]
        ind = findall(x->x>0, considered_Neigh[i,:])
        println("  # of 1st NN, 2nd NN ... :", considered_Neigh[i,ind])
    end
    println()
    
    print_magnon(xpt, Magnon_disp, lattind, kindex)
    
end

#Plots.plot(xpt,real(Magnon_disp[:,1]), width =3, yrange =[0])
#Plots.plot!([xpt[lattind]], seriestype="vline")
#Plots.plot(xpt,real(Magnon_disp[:,2]), width =3, color = "red", yrange = [0,125])


