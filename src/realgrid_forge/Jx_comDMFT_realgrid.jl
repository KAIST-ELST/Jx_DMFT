import DFTforge

import Plots
import ColorSchemes
import FFTW
import JSON
import Dierckx
using JSON
using ImageFiltering
using Dierckx
using Distributed
using Printf
using LinearAlgebra
using Plots
using ColorSchemes

function Findmu(tot_Aw, targ_N, findmu_range, wlist)

    criteria = 0.00001
    
    println("---- we will find the mu making total occ. = ",targ_N," ----")
    println("")
    
    left_lim = findmu_range[1]
    right_lim = findmu_range[2]
    mid_lim = (findmu_range[1] + findmu_range[2] )/2.0
    
    left_N =  trapezoidalInt(wlist,tot_Aw,[wlist[1],left_lim])
    right_N =  trapezoidalInt(wlist,tot_Aw,[wlist[1],right_lim])
    mid_N = trapezoidalInt(wlist,tot_Aw,[wlist[1],mid_lim])
    
    i=1
    println("[",i,"] trial mu : ", mid_lim," , occ : ", mid_N ) 
    
    if (abs(mid_N-targ_N)<criteria)
        println("chemical potential is founded !!  mu = ",mid_lim)
        return mid_lim
    end
    
    
    for i = 2:100
        if (mid_N >= targ_N)
            right_lim = mid_lim
            mid_lim = (left_lim + right_lim)/2.0
        else
            left_lim = mid_lim
            mid_lim = (left_lim + right_lim)/2.0
        end

        left_N =  trapezoidalInt(wlist,tot_Aw,[wlist[1],left_lim])
        right_N =  trapezoidalInt(wlist,tot_Aw,[wlist[1],right_lim])
        mid_N = trapezoidalInt(wlist,tot_Aw,[wlist[1],mid_lim])
    
        println("[",i,"] trial mu : ", mid_lim," , occ : ", mid_N ) 
    
   
        if (abs(mid_N-targ_N)<criteria)
            println("chemical potential is founded !!  mu = ",mid_lim)
            return mid_lim
        end    
        
    end
    
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



function Awsum_for_FS(Akw)
    Awsum = zeros(size(Akw)[1],1,size(Akw)[3],size(Akw)[4])
    for k=1:size(Akw)[1]
        for l =1:size(Akw)[3]
            for s=1:size(Akw)[4]
                Awsum[k,1,l,s]=Akw[k,1,l,s]+Akw[k,2,l,s]
            end
        end
    end
    return Awsum
end

function FS_kgrid(K_grid_num)

    k1list=range(-0.5, stop=0.5, length =K_grid_num[1])
    k2list=range(-0.5, stop=0.5, length =K_grid_num[2])
    k3list=range(-0.5, stop=0.5, length =K_grid_num[3])

    klist=[]


    for k3 =1:(K_grid_num[3])
        klist_dum =zeros(K_grid_num[1]*K_grid_num[2],3)
        ind =0
        for k1 = 1:(K_grid_num[1])
            for k2 = 1:(K_grid_num[2])
                ind += 1
                klist_dum[ind,:] = [k1list[k1],k2list[k2],k3list[k3]]
            end
        end
        push!(klist,klist_dum)
    end


    return klist
end


function zeroselfE(selfE,wlist)

    if (wlist[findmin(abs.(wlist))[2]]>0)
        zind = findmin(abs.(wlist))[2]
    else
        zind = findmin(abs.(wlist))[2]+1
    end

    zwlist=[(wlist[zind]+wlist[zind])/2.0, wlist[zind]]

    zselfE = []
    for i=1:size(selfE)[1]
        zselfE_dum = zeros(ComplexF64,2,size(selfE[i])[2], size(selfE[i])[3])
        for s=1:size(selfE[i])[3]
            zselfE_dum[1,:,s] =(selfE[i][zind,:,s]+selfE[i][zind-1,:,s])/2.0
            zselfE_dum[2,:,s] =selfE[i][zind,:,s]
        end
        push!(zselfE, zselfE_dum)
    end

    return zwlist, zselfE
end






function percent_lowerbound(criper,spectralDat)
    delstart = 0.0
    pct = 0.0
    for i = 1:5000
        a=findall(x->x>delstart , spectralDat/maximum(spectralDat)*10)
        pct = size(a)[1]/size(Matdat)[1]/size(Matdat)[2]*100
        if pct > criper
            delstart+=0.002
        else
            break;
        end
    end

    println("\n Lower bound of the upper ",pct,"% of the normalized spectral weigh :", delstart)
    return pct, delstart
end


function read_Akw(filename)
    kptN = 0
    wptN = 0
    orbN = 0
    spin = 0

    xlist =[]
    wlist =[]

    open(filename, "r") do f
        Nstring=readlines(f)
        first_kpt = parse(Float64, split(Nstring[2])[1])
        spin = parse(Int64, string(split(Nstring[1])[end][end]))
        orbN =convert(Int64, (size(split(Nstring[1]))[1]-2)/spin)
        for i =2:size(Nstring)[1]
            if first_kpt<parse(Float64, split(Nstring[i])[1])
                kptN = convert(Int64, (size(Nstring)[1]-1)/(i-2) )
                wptN = (i-2)
                println("# of kpt: ", kptN)
                println("# of wpt: ", wptN)
                println("# of orb: ", orbN)
                println("spin ind: ", spin)
                break;
            end
        end

    end


    Akw=zeros(kptN,wptN,orbN,spin)

    open(filename, "r") do f
        Nstring=readlines(f)

        lind = 2
        for k = 1:kptN
            push!(xlist,   parse(Float64, split(Nstring[lind])[1]))

            for w = 1:wptN
                if k ==1
                    push!(wlist,   parse(Float64, split(Nstring[lind])[2]))
                end

                for s=1:size(Akw)[4]
                    for l = 1:size(Akw)[3]
                        Akw[k,w,l,s] = parse(Float64, split(Nstring[lind])[l+((s-1)*l)+2])
                    end
                end
                lind+=1
            end
        end

    end



    return xlist, wlist, Akw
end

function print_Aw(Aw,wlist)
    open("Dos_Aw.dat", "w") do io

        @printf(io,"      %8s      ","wpt")
        for s =1:size(Aw)[3]
            for l=1:size(Aw)[2]
                tag_string=string("ORB_",l,"_SPIN_",s)
                @printf(io,"  %16s  ",tag_string)
                if l==size(Aw)[2] && s==size(Aw)[3]
                    @printf(io, "\n")
                end
            end
        end


        for w = 1:size(Aw)[1]
            @printf(io,"  %16.10f  ",wlist[w])
            for s=1:size(Aw)[3]
                for l = 1:size(Aw)[2]
                    @printf(io,"  %16.10f  ",Aw[w,l,s])
                    if l==size(Aw)[2] && s==size(Aw)[3]
                        @printf(io, "\n")

                    end
                end

            end
        end
    end

end




function print_Akw(Akw,kpt,xpt,wlist,calculation_mode)

    if ((calculation_mode=="band") || (calculation_mode=="integrated_band"))

        println("")
        println("Writing 'Band_Akw.dat' file... ")
        println("")

        open("Band_Akw.dat", "w") do io
            @printf(io,"      %8s      ","xpt")
            @printf(io,"      %8s      ","wpt")

            for s =1:size(Akw)[4]
                for l=1:size(Akw)[3]
                    tag_string=string("ORB_",l,"_SPIN_",s)
                    @printf(io,"  %16s  ",tag_string)
                    if l==size(Akw)[3] && s==size(Akw)[4]
                        @printf(io, "\n")
                    end
                end
            end


            for k = 1:size(Akw)[1]
                for w = 1:size(Akw)[2]
                    @printf(io,"  %16.10f  ",xpt[k])
                    @printf(io,"  %16.10f  ",wlist[w])
                    for s=1:size(Akw)[4]
                        for l = 1:size(Akw)[3]
                            @printf(io,"  %16.10f  ",Akw[k,w,l,s])

                            if l==size(Akw)[3] && s==size(Akw)[4]
                                @printf(io, "\n")
                            end
                        end
                    end
                end
            end
        end

    elseif (calculation_mode=="FS")

        println("")
        println("Writing 'FS_Ak.dat' file... ")
        println("")
        open("FS_Ak.dat", "w") do io
            @printf(io,"      %8s      ","k1pt")
            @printf(io,"      %8s      ","k2pt")

            for s =1:size(Akw)[4]
                for l=1:size(Akw)[3]
                    tag_string=string("ORB_",l,"_SPIN_",s)
                    @printf(io,"  %16s  ",tag_string)
                    if l==size(Akw)[3] && s==size(Akw)[4]
                        @printf(io, "\n")
                    end
                end
            end


            k=0
            for k1 = 1:size(Akw)[1]
                for k2 = 1:size(Akw)[2]
                    k+=1
                    @printf(io,"  %16.10f  ",kpt[k,1])
                    @printf(io,"  %16.10f  ",kpt[k,2])
                    for s=1:size(Akw)[4]
                        for l = 1:size(Akw)[3]
                            @printf(io,"  %16.10f  ",Akw[k1,k2,l,s])
                            if l==size(Akw)[3] && s==size(Akw)[4]
                                @printf(io, "\n")
                            end
                        end
                    end
                end
            end

        end

    end

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
            #InvWiess_ineq_comp[j] = Gloc_inv_block_agrid[orb_ind] + (g_SelfE_w[i][w,j,spin] - g_imp_dc[i])
            InvWiess_ineq_comp[j] = Gloc_inv_block_agrid[orb_ind] + (g_SelfE_w[i][w,j,spin])

        end
    end
    
    return InvWiess_ineq_comp
end



function Cal_InvWiess_fromGloc(iWnlist, G_loc_iWn,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,DMFT_spin_type)
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



function Cal_hyb_fromInvWiess(InvWiess_iWn,imp_ind,iWnlist,mu,H_loc,DMFT_spin_type,Corr_atom_equiv,Corr_orbital_Ind,delta)
    hyb_iWn=[]
    for i =1:size(InvWiess_iWn)[1]
        hyb_iWn_Block = zeros(ComplexF64,size(InvWiess_iWn[i]) )
        mat_ind = Corr_orbital_Ind[findall(x->x==i, Corr_atom_equiv)[1]]
        for w = 1:size(iWnlist)[1]
            for s =1:DMFT_spin_type
                for j in unique(imp_ind[ i ])
                    if j != 0
                        orb_ind =  findall(x->x == j, imp_ind[ i ] )[1]
                        hyb_iWn_Block[w,j,s]=(iWnlist[w]+mu+delta*im)-InvWiess_iWn[i][w,j,s]-H0[mat_ind,mat_ind][orb_ind] #Diagonal(diag(H_loc[Corr_orbital_Ind[i],Corr_orbital_Ind[i],s]))

                    end
                end
            end
        end
        push!(hyb_iWn,hyb_iWn_Block)
        
    end
    return hyb_iWn
end



function print_hyb(hyb_w,wlist,Corr_atom_equiv,DMFT_spin_type,imp_ind)
    open("hyb_w.dat", "w") do io
        strRe=""
        strIm=""
        iwn_str="wpts"
        @printf(io," %20s",iwn_str)
        if DMFT_spin_type>0
            for i in unique(abs.(Corr_atom_equiv))
                for s =1:DMFT_spin_type
                    for l in unique(imp_ind[i])
                        if l>0
                            strRe = string("Re{hyb_w_",i,"}[",l,",",s,"]")
                            strIm = string("Im{hyb_w_",i,"}[",l,",",s,"]")   
                            @printf(io, " %20s %20s", strRe, strIm)
                        end
                    end
                end
            end
        end
        
        @printf(io,"\n")

        for w = 1:size(wlist)[1]
            @printf(io, " %20.14f",wlist[w])
            if DMFT_spin_type>0
                for i in unique(abs.(Corr_atom_equiv))
                    for s =1:DMFT_spin_type
                        for l in unique(imp_ind[i])
                            if l>0
                                @printf(io, " %20.14f %20.14f", real(hyb_w[i][w,l,s]), imag(hyb_w[i][w,l,s]))
                            end
                        end
                    end
                end
            end
            @printf(io,"\n")
        
        end
    
        
    end;
end






@everywhere function passing_dc(imp_dc)
    global g_imp_dc
   
    g_imp_dc = imp_dc
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
    

function band_kpts(kpath,grid)
    
    kpt = zeros((size(kpath)[1]-1)*(grid-1)+1,3)
    xpt = []
    lattind = []
    push!(lattind, 1)
    push!(xpt,0)
    totlength=0
    
    kind = 1
    kpt[1,:]=kpath[1]
    for i=2:size(kpath)[1]
     
        k1=range(kpath[i-1][1], stop=kpath[i][1], length =grid)
        k2=range(kpath[i-1][2], stop=kpath[i][2], length =grid)
        k3=range(kpath[i-1][3], stop=kpath[i][3], length =grid)
        
        length_ = ((kpath[i][1]-kpath[i-1][1])^2+(kpath[i][2]-kpath[i-1][2])^2+(kpath[i][3]-kpath[i-1][3])^2)^(1/2)
        totlength += length_
        dum_leng=range(totlength-length_,stop=totlength, length=grid)
        
        for j=2:grid
        
            kind=kind+1
            kpt[kind,:]=[k1[j],k2[j],k3[j]]
            push!(xpt,dum_leng[j])
            if (j==grid)
                push!(lattind, kind) 
            end
        end
    end
    
    return kpt, xpt,lattind
end



function integrated_band_kpts(kpath,grid,kz_grid_for_integrated)
    println("\n Caution !!! : For integrated band mode, your kpath file should contain only the in-plane points (kz=0) !!! \n  ")
    int_kpt =[]
    kpt = zeros((size(kpath)[1]-1)*(grid-1)+1,3)
    xpt = []
    lattind = []
    push!(lattind, 1)
    push!(xpt,0)
    totlength=0
    
    kind = 1
    kpt[1,:]=kpath[1]
    for i=2:size(kpath)[1]
     
        k1=range(kpath[i-1][1], stop=kpath[i][1], length =grid)
        k2=range(kpath[i-1][2], stop=kpath[i][2], length =grid)
        k3=range(kpath[i-1][3], stop=kpath[i][3], length =grid)
        
        length_ = ((kpath[i][1]-kpath[i-1][1])^2+(kpath[i][2]-kpath[i-1][2])^2+(kpath[i][3]-kpath[i-1][3])^2)^(1/2)
        totlength += length_
        dum_leng=range(totlength-length_,stop=totlength, length=grid)
        
        for j=2:grid
        
            kind=kind+1
            kpt[kind,:]=[k1[j],k2[j],k3[j]]
            push!(xpt,dum_leng[j])
            if (j==grid)
                push!(lattind, kind) 
            end
        end
    end
    
    kz_grid = range(0.0, stop=0.5, length =kz_grid_for_integrated)
    for kz in kz_grid
        kpt_dum = deepcopy(kpt)
        kpt_dum[:,3].= kz
        push!(int_kpt,kpt_dum)
    end
    
    
    return int_kpt, xpt,lattind
end




function beta_cal(Temp)
    BoltzmannC = 8.617333262/100000;
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




function input_handler(DMFT_Jx_opion)

   if (DMFT_Jx_opion.go_DMFT && !DMFT_Jx_opion.go_Jx)
        println(" go_DMFT :", DMFT_Jx_opion.go_DMFT, ", go_Jx :",DMFT_Jx_opion.go_Jx)
        ## Constant
        BoltzmannC = 8.617333262/100000;
        delta =1e-7; # for grennf occupation count


        ### input parameter
        K_grid_num = DMFT_Jx_opion.KgridNum;
        R_grid_num = DMFT_Jx_opion.RgridNum;
        iW_grid_cut = DMFT_Jx_opion.iWgridCut;
        Tem= DMFT_Jx_opion.Temperature;
        beta = beta_cal(Tem);
        beta_for_occ = beta_cal(100.0);

        DMFT_spin_type = DMFT_Jx_opion.Spin_type;


        Corr_atom_Ind = DMFT_Jx_opion.Corr_atom_Ind;
        Corr_orbital_Ind = DMFT_Jx_opion.Corr_orbital_Ind;
        Corr_atom_equiv = DMFT_Jx_opion.Corr_atom_equiv;

        Ineq_atom_Ind=unique(abs.(Corr_atom_equiv));

        if DMFT_spin_type == 2
            mag_order = ones(Int64, size(Corr_atom_equiv)[1]);
            mag_order[findall(x->x <0, Corr_atom_equiv)].=-1;
        end

        DMFT_loop_N = DMFT_Jx_opion.DMFT_loop_N;


        Mix_selfE = DMFT_Jx_opion.Mix_selfE;
        init_bias = DMFT_Jx_opion.init_bias;
        smth_step = DMFT_Jx_opion.smth_step;
        cal_susc = DMFT_Jx_opion.cal_susc;

        imp_ind = DMFT_Jx_opion.imp_block;
        imp_lev_shift = DMFT_Jx_opion.imp_lev_shift;
        imp_U =  DMFT_Jx_opion.imp_U;
        imp_J =  DMFT_Jx_opion.imp_J;
        imp_dc_n =  DMFT_Jx_opion.imp_dc;
        imp_dc = FLL_dc(imp_dc_n, imp_U, imp_J);
        #imp_U.*(imp_dc_n.-1/2).-imp_J/2.0.*(imp_dc_n.-1);
        imp_dc_type = DMFT_Jx_opion.imp_dc_type;
        imp_int_type =  DMFT_Jx_opion.imp_int_type;
        imp_int_parameterisation =  DMFT_Jx_opion.imp_int_parameterisation;
        
        
        Corr_ineq_orbital_Num = []
        for i = 1:size(imp_ind)[1]
          push!(Corr_ineq_orbital_Num,count(x->x>0, unique(imp_ind[i])));
        end
    
        Solver_green_Cut = DMFT_Jx_opion.Solver_green_Cut;


        F0 = imp_U;
        F2 = imp_J*14/2.6*1.6;
        F4 = imp_J*14/2.6*1.0;

        imp_Measure_time = DMFT_Jx_opion.imp_Measure_time;
        imp_Thermal_time = DMFT_Jx_opion.imp_Thermal_time;
        
        
        basis_transform = DMFT_Jx_opion.basis_transform;
        
        green_basis = DMFT_Jx_opion.green_basis;
        green_legendre_cutoff = DMFT_Jx_opion.green_legendre_cutoff;
        return (BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut, Tem, beta, beta_for_occ,
                DMFT_spin_type, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind,
                DMFT_loop_N, Mix_selfE, init_bias,smth_step,cal_susc,imp_dc_type, imp_ind, imp_lev_shift, imp_U, imp_J, imp_dc_n, imp_dc, imp_int_type, imp_int_parameterisation
                , Corr_ineq_orbital_Num, Solver_green_Cut ,F0, F2, F4, imp_Measure_time, imp_Thermal_time, basis_transform, green_basis, green_legendre_cutoff);
        
        
    elseif (DMFT_Jx_opion.go_DMFT && DMFT_Jx_opion.go_Jx)
          println(" go_DMFT :", DMFT_Jx_opion.go_DMFT, ", go_Jx :",DMFT_Jx_opion.go_Jx)
        ## Constant
        BoltzmannC = 8.617333262/100000;
        delta =1e-7; # for grennf occupation count


        ### input parameter
        K_grid_num = DMFT_Jx_opion.KgridNum;
        R_grid_num = DMFT_Jx_opion.RgridNum;
        iW_grid_cut = DMFT_Jx_opion.iWgridCut;
        Tem= DMFT_Jx_opion.Temperature;
        beta = beta_cal(Tem);
        beta_for_occ = beta_cal(100.0);


        DMFT_spin_type = DMFT_Jx_opion.Spin_type;


        Corr_atom_Ind = DMFT_Jx_opion.Corr_atom_Ind;
        Corr_orbital_Ind = DMFT_Jx_opion.Corr_orbital_Ind;
        Corr_atom_equiv = DMFT_Jx_opion.Corr_atom_equiv;

        Ineq_atom_Ind=unique(abs.(Corr_atom_equiv));

        if DMFT_spin_type == 2
            mag_order = ones(Int64, size(Corr_atom_equiv)[1]);
            mag_order[findall(x->x <0, Corr_atom_equiv)].=-1;
        end

        DMFT_loop_N = DMFT_Jx_opion.DMFT_loop_N;


        Mix_selfE = DMFT_Jx_opion.Mix_selfE;
        init_bias = DMFT_Jx_opion.init_bias;
        smth_step = DMFT_Jx_opion.smth_step;
        cal_susc = DMFT_Jx_opion.cal_susc;

        imp_ind = DMFT_Jx_opion.imp_block;
        imp_U =  DMFT_Jx_opion.imp_U;
        imp_J =  DMFT_Jx_opion.imp_J;
        imp_dc_n =  DMFT_Jx_opion.imp_dc;
        imp_dc = FLL_dc(imp_dc_n, imp_U, imp_J);
        imp_int_type =  DMFT_Jx_opion.imp_int_type;
        imp_int_parameterisation =  DMFT_Jx_opion.imp_int_parameterisation;

        Corr_ineq_orbital_Num = []
        for i = 1:size(imp_ind)[1]
          push!(Corr_ineq_orbital_Num,count(x->x>0, unique(imp_ind[i])));
        end
    
        Solver_green_Cut = DMFT_Jx_opion.Solver_green_Cut;


        F0 = imp_U;
        F2 = imp_J*14/2.6*1.6;
        F4 = imp_J*14/2.6*1.0;

        imp_Measure_time = DMFT_Jx_opion.imp_Measure_time;
        imp_Thermal_time = DMFT_Jx_opion.imp_Thermal_time;
        
        
        basis_transform = DMFT_Jx_opion.basis_transform;
        
        green_basis = DMFT_Jx_opion.green_basis;
        green_legendre_cutoff = DMFT_Jx_opion.green_legendre_cutoff;
        return (BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut, Tem, beta, beta_for_occ,
                DMFT_spin_type, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind,
                DMFT_loop_N, Mix_selfE, init_bias,smth_step,cal_susc,imp_ind, imp_U, imp_J, imp_dc_n, imp_dc, imp_int_type,imp_int_parameterisation
                , Corr_ineq_orbital_Num, Solver_green_Cut ,F0, F2, F4, imp_Measure_time, imp_Thermal_time, basis_transform, green_basis, green_legendre_cutoff,mag_order);
        
        
    elseif (!DMFT_Jx_opion.go_DMFT && DMFT_Jx_opion.go_Jx)
        println(" go_DMFT :", DMFT_Jx_opion.go_DMFT, ", go_Jx :",DMFT_Jx_opion.go_Jx)

        BoltzmannC = 8.617333262/100000;
        delta =1e-7; # for grennf occupation count


        ### input parameter
        K_grid_num = DMFT_Jx_opion.KgridNum;
        R_grid_num = DMFT_Jx_opion.RgridNum;
        iW_grid_cut = DMFT_Jx_opion.iWgridCut;
        
        Tem= DMFT_Jx_opion.Temperature;
        beta = beta_cal(Tem);

        Corr_atom_Ind = DMFT_Jx_opion.Corr_atom_Ind;
        Corr_orbital_Ind = DMFT_Jx_opion.Corr_orbital_Ind;
        Corr_atom_equiv = DMFT_Jx_opion.Corr_atom_equiv;

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
    
        mag_order = ones(Int64, size(Corr_atom_equiv)[1]);
        mag_order[findall(x->x <0, Corr_atom_equiv)].=-1;        

        DMFT_spin_type = 0;
        
        

        
        return (BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut,Tem, beta,DMFT_spin_type,
                 Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind, Corr_ineq_orbital_Num,imp_ind,mag_order);
        
        
    end
end


@everywhere function init_variables_H_k(H_k)
    global g_H_k
   
    #g_SelfE_w = SelfE_w
    g_H_k = H_k;

end



@everywhere function init_variables_DMFT(Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv, imp_ind, Corr_ineq_orbital_Num, delta)
    global g_Corr_atom_Ind, g_Corr_orbital_Ind, g_Corr_atom_equiv, g_imp_ind, g_Corr_ineq_orbital_Num, g_delta
   
    g_Corr_atom_Ind = Corr_atom_Ind;
    g_Corr_orbital_Ind = Corr_orbital_Ind;
    g_Corr_atom_equiv = Corr_atom_equiv;
    g_imp_ind = imp_ind;
    g_Corr_ineq_orbital_Num = Corr_ineq_orbital_Num;
    g_delta = delta;
    
end

@everywhere function init_variables_SelfE_w(SelfE_w)
    global g_SelfE_w
   
    #g_SelfE_w = SelfE_w
    g_SelfE_w = SelfE_w;

end




@everywhere function init_variables_grid_DMFT_realf(wlist,Klist)
    global g_wlist,g_Klist
   
    #g_SelfE_w = SelfE_w
    g_wlist = wlist;
    g_Klist = Klist;

end



@everywhere function init_variables_totorb(Total_orb_Num)
    global g_Total_orb_Num
   
    #g_SelfE_w = SelfE_w
    g_Total_orb_Num = Total_orb_Num;


end

@everywhere function init_variables_calmode(calculation_mode)
    global g_calculation_mode
   
    #g_SelfE_w = SelfE_w
    g_calculation_mode = calculation_mode;


end


function Tvec_print(H0, transMat, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)
    Trans_mat = []
    Trans_mat_A = Matrix{Float64}(I,size(H0)[1],size(H0)[2])
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


function isfile_tranformMat()
   
    if isfile("transform_Mat.dat")
       
        println("\n         =========================================================================")
        println("There is 'transform_Mat.dat' file, wannier hamiltonian will be transformed by using the transform matrix")
        println("         =========================================================================\n")
    
        tof=true
    else
        println("\n =========================================")
        println("--- There is no 'transform_Mat.dat' file ---")
        println(" =========================================\n")

        tof=false
    end
    
    return tof
end



function local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
    if ((calculation_mode =="dos") || (calculation_mode =="DFT-dos"))
        H0=zeros(ComplexF64,size(H_k)[4],size(H_k)[5])
        for k1 =1:size(H_k)[1]
            for k2 =1:size(H_k)[2]
                for k3=1:size(H_k)[3]
                    H0 += H_k[k1,k2,k3,:,:,1]/(size(H_k)[1]*size(H_k)[2]*size(H_k)[3])
                end
            end
        end
    elseif ((calculation_mode =="band") || (calculation_mode =="integrated_band") || (calculation_mode =="FS") )
        H0=zeros(ComplexF64,size(H_k)[2],size(H_k)[3])  
        for k =1:size(H_k)[1]
            H0 += H_k[k,:,:,1]/(size(H_k)[1])
        end    
    end

    H0_lev =[]
    for i=1:size(imp_ind)[1]
        orb_ind = findall(x->x==i,Corr_atom_equiv)[1]
        println("----------------------------------------------------")
        println("Inequv. orbital index of imp. ",i," :", Corr_orbital_Ind[orb_ind])
        println("Diagonal part of Hloc : ",real(diag(H0[Corr_orbital_Ind[orb_ind],Corr_orbital_Ind[orb_ind]])))
        println("\n")
    end
    
    return H0
end


function Trans_H_k(H_k,Trans_mat_A, calculation_mode)
    dum_H_k = zeros(size(H_k))+ zeros(size(H_k))*im

    if ((calculation_mode == "dos")||(calculation_mode == "DFT-dos"))
        for k1 = 1:size(H_k)[1]
            for k2 = 1:size(H_k)[2]
                for k3 = 1:size(H_k)[3]
                    for s= 1:size(H_k)[6]
                        dum_H_k[k1,k2,k3,:,:,s] = Trans_mat_A * H_k[k1,k2,k3,:,:,s] * Trans_mat_A'
                    end
                end
            end
        end
    elseif ((calculation_mode == "band") || (calculation_mode == "integrated_band") || (calculation_mode == "FS") )
        for k = 1:size(H_k)[1]
            for s=1:size(H_k)[4]
            
                dum_H_k[k,:,:,s] = Trans_mat_A * H_k[k,:,:,s] * Trans_mat_A'
            end
        end        
    end
    
    return dum_H_k
end

    



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








function nonInt_H_k(Klist,hamiltonian_info,go_DMFT,go_Jx,DMFT_spin_type)
    #initiaize H_k
    
    if string(hamiltonian_info.dfttype) == "Wannier90"
        
        if (go_Jx || (DMFT_spin_type==2) )
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
                    for s=1:hamiltonian_info.scf_r.SpinP_switch
                        Kvect = Klist[k1,k2,k3,:]
                        push!(arg_list,( Kvect, hamiltonian_info.scf_r.R_vector_mat[s], hamiltonian_info.scf_r.Hks_R[s] ) )
                    

                    end
                end
            end
        end
    
        H_k_list_tmp = pmap(nonInt_H_R2k_internal, arg_list)
        i=0
        for k1 = (1:size(Klist)[1])
            for k2 = (1:size(Klist)[2])
                for k3 = (1:size(Klist)[3])
                    for s=1:hamiltonian_info.scf_r.SpinP_switch
                        i=i+1
                        H_k[k1,k2,k3,:,:,s]=H_k_list_tmp[i];

                    end
                end
            end
        end    
        if ((go_DMFT && go_Jx) || (DMFT_spin_type == 2))
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
        
        

              
    
    end
    
    
    
    return H_k
end
        



function nonInt_H_k_band(Klist,hamiltonian_info,go_DMFT,go_Jx,DMFT_spin_type)
    #initiaize H_k
    
    if string(hamiltonian_info.dfttype) == "Wannier90"
        
        if (go_Jx || (DMFT_spin_type==2) )
            H_k=zeros(size(Klist)[1],sum(hamiltonian_info.scf_r.Total_NumOrbs),sum(hamiltonian_info.scf_r.Total_NumOrbs),2)+zeros(size(Klist)[1],sum(hamiltonian_info.scf_r.Total_NumOrbs),sum(hamiltonian_info.scf_r.Total_NumOrbs),2)*im
        else
            H_k=zeros(size(Klist)[1],sum(hamiltonian_info.scf_r.Total_NumOrbs),sum(hamiltonian_info.scf_r.Total_NumOrbs),1)+zeros(size(Klist)[1],sum(hamiltonian_info.scf_r.Total_NumOrbs),sum(hamiltonian_info.scf_r.Total_NumOrbs),1)*im
        end
        arg_list=[]
        #load H_k
        println("===== loading H_k from tight-binding parameter ... ======") 
        for k = (1:size(Klist)[1])
            for s=1:hamiltonian_info.scf_r.SpinP_switch
                Kvect = Klist[k,:]
                push!(arg_list,( Kvect, hamiltonian_info.scf_r.R_vector_mat[s], hamiltonian_info.scf_r.Hks_R[s] ) )

            end
        end
    
        H_k_list_tmp = pmap(nonInt_H_R2k_internal, arg_list)
        i=0
        for k = (1:size(Klist)[1])
            for s=1:hamiltonian_info.scf_r.SpinP_switch
                i=i+1
                H_k[k,:,:,s]=H_k_list_tmp[i];
            end
        end    
        if ((go_DMFT && go_Jx) || (DMFT_spin_type == 2))
            H_k[:,:,:,2]= H_k[:,:,:,1]
        end
        
    elseif string(hamiltonian_info.dfttype) == "OpenMX"
        
        H_k=zeros(size(Klist)[1],sum(hamiltonian_info.scf_r.Total_NumOrbs),sum(hamiltonian_info.scf_r.Total_NumOrbs),hamiltonian_info.scf_r.SpinP_switch+1)+zeros(size(Klist)[1],sum(hamiltonian_info.scf_r.Total_NumOrbs),sum(hamiltonian_info.scf_r.Total_NumOrbs),hamiltonian_info.scf_r.SpinP_switch+1)*im

        arg_list=[]
        #load H_k
        println("===== loading H_k from tight-binding parameter ... ======") 
        for k = (1:size(Klist)[1])
            Ktuple = tuple(Klist[k,:]...)
            push!(arg_list, Ktuple )

        end
        H_k_list_tmp = pmap(nonInt_H_k_OpenMX, arg_list)
        i=0
        for k1 = (1:size(Klist)[1])
            i=i+1
            H_k[k,:,:,1]=H_k_list_tmp[i][:,:,1];
            H_k[k,:,:,2]=H_k_list_tmp[i][:,:,2];
  
        end            
        
    end
    
    
    
    return H_k
end



function w_Num_orb_Num(filename,imp_ind,DMFT_spin_type)

    orb_size = 0
    Num_freq = 0
    re = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$";
    stind = 0    
    open(filename, "r") do f 
        Nstring=readlines(f)
        
        for i=1:size(Nstring)[1]
            if occursin(re,split(Nstring[i])[1])
                stind=stind+1
                break;
            else
               stind = stind+1 
            end
        end
        
        orb_size = convert(Int64, size(split(Nstring[3])[2:end])[1]/2) 
        Num_freq =  size(Nstring)[1]-stind+1

    end

    println("# of orbital :  ", orb_size)
    println("# of freq    :  ", Num_freq)
    
    real_freq = zeros(Num_freq)

        
        
    open(filename, "r") do f 
        Nstring=readlines(f)
        for i=stind:size(Nstring)[1]
            real_freq[i+1-stind]=parse(Float64,split(Nstring[i])[1])
        end
    end    
    
    return real_freq, orb_size, Num_freq, stind
end




function load_real_selfE(filename, DMFT_spin_type, imp_ind, orb_size, Num_freq, stind)

    tot_orb_num =0
    selfE_on_realf =[]
    orb_N=[]
    for i=1:size(imp_ind)[1]
        orb_N_imp = size(findall(x->x>0,unique(imp_ind[i])))[1]
        tot_orb_num+=orb_N_imp
        selfE_on_realf_imp = zeros(ComplexF64, Num_freq,orb_N_imp,DMFT_spin_type )
        push!(selfE_on_realf,selfE_on_realf_imp)
        push!(orb_N,orb_N_imp)
    end
    


    open(filename, "r") do f 
        Nstring=readlines(f)
        for i=stind:size(Nstring)[1]
            orb_N_sum = 0;
            for im=1:size(imp_ind)[1]
                for s=1:DMFT_spin_type
                    for j=1:orb_N[im]
                        selfE_on_realf[im][i+1-stind,j,s]= complex(parse(Float64,split(Nstring[i])[((orb_N_sum+j)*2)+(orb_N[im]*(s-1)*2)]),parse(Float64,split(Nstring[i])[((orb_N_sum+j)*2+1)+(orb_N[im]*(s-1)*2)]))
                    end
                end
                orb_N_sum += orb_N[im]
            end
        end
    end    


   



    return selfE_on_realf
end




function init_selfE_zero_realf(Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,w_num)

    SelfE_w=[]
    for i in Ineq_atom_Ind
   
        SelfE_block= zeros(w_num,Corr_ineq_orbital_Num[i],DMFT_spin_type)+zeros(w_num,Corr_ineq_orbital_Num[i],DMFT_spin_type)*im;
   
        push!(SelfE_w,SelfE_block)
    end
    
    return SelfE_w
end




@everywhere function SelfE_at_agrid(Wind,spin)

    if ((g_calculation_mode == "dos") || (g_calculation_mode == "DFT-dos"))
        SelfE_agrid=zeros(size(g_H_k[1,1,1,:,:,1]))+zeros(size(g_H_k[1,1,1,:,:,1]))*im
    elseif ((g_calculation_mode == "band") || (g_calculation_mode == "integrated_band") || (calculation_mode == "FS") )
        SelfE_agrid=zeros(size(g_H_k[1,:,:,1]))+zeros(size(g_H_k[1,:,:,1]))*im
    end
       
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

    end
    return SelfE_agrid
end


@everywhere function gaussian_selfE(delta,SelfE_agrid)

    #broadening method1
    #SelfE=deepcopy(SelfE_agrid)
    #oi1 = findall(x->(imag(x)!=0.0 && -imag(x)<delta), SelfE_agrid)
    #SelfE[oi1].=real(SelfE[oi1])
    #oi2 = findall(x->(x==1.0+0.0im || (imag(x)!=0.0 && -imag(x)<delta)) , SelfE_agrid+Matrix{Float64}(I,size(g_H_k[1,1,1,:,:,1]) ))
    #SelfE[oi2].-= delta*im
    
    #broadening method2
    SelfE=deepcopy(SelfE_agrid)
    SelfE = SelfE-Matrix{Float64}(I,size(g_H_k[1,1,1,:,:,1]) )*delta*im
    
    
    return SelfE
end




@everywhere function GreenF_at_agrid_realf(w,Kind,spin,SelfE_agrid,testmu,delta)
    #deltaM = deltaMat(delta,SelfE_agrid)
    
    #SelfE = gaussian_selfE(delta,SelfE_agrid)
    
    #deltaM = zeros(ComplexF64,size(g_H_k[1,1,1,:,:,1]))
    #Gmat=inv( ( (w+testmu)*Matrix{Float64}(I,size(g_H_k[1,1,1,:,:,1]) )) + deltaM - g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin] - SelfE_agrid )

    if ((g_calculation_mode =="dos") || (g_calculation_mode == "DFT-dos") )
        Gmat = zeros(size(g_H_k[1,1,1,:,:,1]))+zeros(size(g_H_k[1,1,1,:,:,1]))*im 
        Gmat =inv( ( (w+testmu+delta*im)*Matrix{Float64}(I,size(g_H_k[1,1,1,:,:,1]) )) - g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin] - SelfE_agrid )
    elseif ((g_calculation_mode =="band") || (g_calculation_mode =="integrated_band") || (calculation_mode == "FS") )
        Gmat = zeros(size(g_H_k[1,:,:,1]))+zeros(size(g_H_k[1,:,:,1]))*im 
        Gmat =inv( ( (w+testmu+delta*im)*Matrix{Float64}(I,size(g_H_k[1,:,:,1]) )) - g_H_k[Kind,:,:,spin] - SelfE_agrid )       
    end
        #Gmat -= Gmat'

    #Gmat-=inv( ( (w+testmu-delta*im)*Matrix{Float64}(I,size(g_H_k[1,1,1,:,:,1]) )) - g_H_k[Kind[1],Kind[2],Kind[3],:,:,spin] - SelfE_agrid )

    #Gmat./=2.0
    
    return Gmat
end



@everywhere function G_k_sum_realf(args)
    w = args[1]
    spin = args[2]
    testmu = args[3]
    delta = args[4]
    println("Calculating G(w) at w =",w, "... ")

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


function Gloc_gen_realf(H_k,iWnlist,SelfE_w,testmu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind,delta)
    G_loc_iWn = zeros(size(iWnlist)[1],size(H_k)[4],size(H_k)[5],size(H_k)[6])+zeros(size(iWnlist)[1],size(H_k)[4],size(H_k)[5],size(H_k)[6])*im
    arg_list = [];
    for w = 1:size(iWnlist)[1]
        for s = 1:size(H_k)[6]
            push!(arg_list,( w,s, testmu, delta) )
        end
    end
    println("passed all variables to calculate G(w) ... ")

    G_w_s_list = pmap(G_k_sum_realf, arg_list)
    i=0
    for w = 1:size(iWnlist)[1]
        for s = 1:size(H_k)[6]
            i+=1
            G_loc_iWn[w,:,:,s]=G_w_s_list[i]/( size(H_k)[1]*size(H_k)[2]*size(H_k)[3] )

        end
    end

    return  G_loc_iWn
end






@everywhere function G_w_k_realf(args)
    w = args[1]
    k = args[2]
    spin = args[3]
    testmu = args[4]
    delta = args[5]
    if k==1
        println("Calculating G(w,k) at w =",w,"... ")
    end

    G_w_k_s=zeros(g_Total_orb_Num,g_Total_orb_Num)+zeros(g_Total_orb_Num,g_Total_orb_Num)*im
    G_w_k_s= GreenF_at_agrid_realf(g_wlist[w],k, spin, SelfE_at_agrid(w,spin),testmu,delta)


    return G_w_k_s
end


function G_gen_realf(H_k,iWnlist,SelfE_w,testmu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind,delta)
    G_w_k_s = zeros(size(iWnlist)[1],size(H_k)[1],size(H_k)[2],size(H_k)[3],size(H_k)[4])+zeros(size(iWnlist)[1],size(H_k)[1],size(H_k)[2],size(H_k)[3],size(H_k)[4])*im
    arg_list = [];
    for w = 1:size(iWnlist)[1]
        for k =1:size(H_k)[1]
            for s = 1:size(H_k)[4]
                push!(arg_list,( w,k,s, testmu, delta) )
            end
        end
    end
    println("passed all variables to calculate G(w) ... ")

    G_w_k_s_list = pmap(G_w_k_realf, arg_list)
    i=0
    for w = 1:size(iWnlist)[1]
        for k =1:size(H_k)[1]
            for s = 1:size(H_k)[4]
                i+=1
                G_w_k_s[w,k,:,:,s]=G_w_k_s_list[i]
            end
        end
    end

    return  G_w_k_s
end










function Aw_orb(orb_selection,wlist, G_loc_w)

    Aw = zeros(size(G_loc_w)[1],size(orb_selection)[1],size(G_loc_w)[4])
    println("")
    println("")
    println("For non-mag cal, total d-orbital occupation 5 !")
    for i = 1:size(orb_selection)[1]
        for s = 1:size(G_loc_w)[4]
            for j in orb_selection[i][1]
                Aw[:,i,s] += -imag(G_loc_w[:,j,j,s])./pi
            end
            
            if size(G_loc_w)[4] == 1
                println("Orb",i,", occup : ", trapezoidalInt(wlist,  Aw[:,i,s], [wlist[1],0.0]) )
            else size(G_loc_w)[4] == 2
                if s == 1
                    println("Orb",i,", up-occup : ", trapezoidalInt(wlist,  Aw[:,i,s], [wlist[1],0.0]) )
                    
                else
                    println("Orb",i,", dn-occup : ", trapezoidalInt(wlist,  Aw[:,i,s], [wlist[1],0.0]) )
                    
                end
            end
            
        end
    end
    
    return Aw
end



function Akw_orb(orb_selection,wlist, G_w_k_s)

    Akw = zeros(size(G_w_k_s)[2],size(G_w_k_s)[1],size(orb_selection)[1],size(G_w_k_s)[5])
    println("")
    println("")
    println("For non-mag cal, total d-orbital occupation 5 !")
    for i = 1:size(orb_selection)[1]
        for k = 1:size(G_w_k_s)[2]
            for w = 1:size(G_w_k_s)[1]
                for j in orb_selection[i][1]
                    for s = 1:size(G_w_k_s)[5]
                        Akw[k,w,i,s] += -imag(G_w_k_s[w,k,j,j,s])./pi
                    end
                end
            end
        end
    end
            
    
    return Akw
end



function Akw_orb_FS(orb_selection,wlist, G_w_k_s, K_grid_num)

    Akw = zeros(K_grid_num[1],K_grid_num[2],size(orb_selection)[1],size(G_w_k_s)[5])
    println("")
    println("")
    println("For non-mag cal, total d-orbital occupation 5 !")
    for i = 1:size(orb_selection)[1]
        k=0
        for k1 = 1:K_grid_num[1]
            for k2 = 1:K_grid_num[2]
                k+=1
                for j in orb_selection[i][1]
                    for s = 1:size(G_w_k_s)[5]
                        Akw[k1,k2,i,s] += -(imag(G_w_k_s[1,k,j,j,s])+imag(G_w_k_s[2,k,j,j,s]))./pi
                    end
                end
            end
        end
    end
            
    
    return Akw
end



function trapezoidalInt(wlist, Aw_ar, range)
    Aw = 0.0
    for i =1:size(wlist)[1]
        if (range[1]<wlist[i]<range[2])
            Aw += (wlist[i+1]-wlist[i])/2.0*(Aw_ar[i]+Aw_ar[i+1])
        end
    end
    return Aw
end

###psudo opt cal

function Avg_Zfactor(selfE,wlist)
    Zfact=zeros(size(selfE)[1])
    for i=1:size(selfE)[1]
        for l = 1:size(selfE[i])[2]
            for s = 1:size(selfE[i])[3]
                mind = maximum(findall(x->x<0, wlist))
                Zfact[i]+=1.0/(1.0-(real(selfE[i][mind,l,s])-real(selfE[i][mind+1,l,s]))/(wlist[mind]-wlist[mind+1]) )
            
                #println(1.0/(1.0-(real(selfE[i][mind,l,s])-real(selfE[i][mind+1,l,s]))/(wlist[mind]-wlist[mind+1]) ))
            end
        end
        Zfact[i]/=size(selfE[i])[2]*size(selfE[i])[3]
    end
    
    return Zfact
end

function wcut_for_optic(wlist, opt_range)
    wlist_opt = wlist[findall(x->x>0 && x<opt_range, wlist)] 
    return wlist_opt
end


function wcalist(wlist, wlist_opt)
    wcal_pack = []
    for i=1:size(wlist_opt)[1]
        w_dum =[]
        for ii in (findall(x->x<0 && x>-wlist_opt[i],wlist))
            nearestV, nearestind = findmin( abs.(wlist.-(wlist[ii]+wlist_opt[i])) );
            if ((wlist[nearestind] >0 )&& (nearestV/wlist_opt[i])<0.01)
               push!(w_dum, [ ii , nearestind]) 
            end
            
        end
        push!(wcal_pack, w_dum)
    end
    return wcal_pack
end


function cal_opt(Aw1, Aw2, Zfact1, Zfact2, wcal_pack )
    opt_intensity=zeros(size(wcal_pack)[1])
    for i=1:size(wcal_pack)[1]
        if size(wcal_pack[i])[1]>1
            for j=1:size(wcal_pack[i])[1]
                opt_intensity[i]+=Aw1[wcal_pack[i][j][1]]*Aw2[wcal_pack[i][j][2]]*Zfact1*Zfact2
            end
        end
        
    end
    return opt_intensity
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



arg_input.TOMLinput = "/home/tj-kim/Dropbox/Jx.jl/src/realgrid_forge/ndnio2_J_wannier.toml"
arg_input.TOMLinput = "/home/tj-kim/Dropbox/Jx.jl/src/realgrid_forge/fgt_J_wannier.toml"
arg_input.TOMLinput = "/home/tj-kim/Dropbox/Jx.jl/examples/FGT/DMFT_U5_J1.0_dc6.0_6.0_1000K_Full/fgt_J_wannier.toml"


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

#result_file= "/home/tjkim/Desktop/software_installation/Jx.jl/src/DMFT_code/wannier_hr.dat"
#result_file= "/home/tj-kim/바탕화면/Jx/Jx.jl/examples/FGT/wannier_hr.dat"
result_file= "/home/tj-kim/바탕화면/Jx/Jx.jl/examples/NiO_G-AFM.OpenMx/nio.HWR"
#result_file= "/home/tj-kim/바탕화면/Jx/Jx.jl/examples/NiO_G_PM.OpenMx/wannier_hr.dat"
#result_file= "/home/tjkim/Desktop/software_installation/Jx.jl/examples/NiO_G-PM.comsuite/wannier_hr.dat"
#result_file = "/home/tj-kim/바탕화면/Jx/Jx.jl/examples/LaNiO2_AFM_G/AFM_G.GGA.openmx.a3.92_c3.31/lanio.scfout"
#result_file = "/home/tjkim/Desktop/software_installation/Jx.jl/examples/LaNiO2_AFM_G/AFM_G.GGA.openmx.a3.92_c3.31/lanio.scfout"


#result_file = "/home/tj-kim/바탕화면/Jx/Jx.jl/examples/LaNiO2_AFM_G/AFM_G.GGA.openmx.a3.92_c3.31/lanio.scfout"

result_file ="/home/tj-kim/바탕화면/Jx/Jx.jl/examples/LaNiO2_AFM_G/Fe6.0-s2p2d2f1-0.0_0.0-2.784/fe.scfout"
result_file ="/home/tj-kim/바탕화면/Jx/Jx.jl/examples/Fe_primitive/fe.scfout"
#result_file ="/home/tj-kim/바탕화면/Jx/Jx.jl/examples/LaNiO2_AFM_G/AFM_G.GGA.openmx.a3.92_c3.31/lanio.scfout"
result_file= "/home/tj-kim/바탕화면/Jx/Jx.jl/examples/NiO_G-AFM.OpenMx/nio.HWR"
result_file= "/home/tj-kim/Dropbox/Jx.jl/examples/FGT/wannier_hr.dat"
result_file= "/home/tj-kim/Dropbox/Jx.jl/src/realgrid_forge/wannier_hr.dat"
result_file= "/home/tj-kim/Dropbox/Jx.jl/examples/FGT/DMFT_U5_J1.0_dc6.0_6.0_1000K_Full/wannier_hr.dat"
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


K_point_groups = arg_input.Optional["bandplot"]
K_point_groups[1].K_point_list

K_point_list_all = Array{k_point_Tuple}(undef,0);

for (i,v) in enumerate(K_point_groups)
  map(x -> push!(K_point_list_all,x) ,v.K_point_list)
end


K_point_list_all_mat =  zeros(length(K_point_list_all),3)
K_point_name_list = Array{Tuple{Int64,String} }(undef,0)


################## removed by TJ ################
#Cnt = 1
#for (i,v) in enumerate(K_point_groups)
    
#    for (i2,K_point) in enumerate(v.K_point_list)
#        if(1 == i2)
#            push!(K_point_name_list,(Cnt,v.K_start_point_name))
#        end
#        K_point_cat = hamiltonian_info.scf_r.rv * [K_point[1],K_point[2],K_point[3]];
#        K_point_cat = [K_point[1],K_point[2],K_point[3]]' * hamiltonian_info.scf_r.rv;
#        K_point_list_all_mat[Cnt,:] = K_point_cat;
        
#        if(i2 == length(v.K_point_list) && i == length(K_point_groups) )
#            push!(K_point_name_list,(Cnt,v.K_end_point_name))
#        end
#        Cnt += 1
#    end
#end
################## removed by TJ ################


# Convert K_point distances
K_point_dist_list = zeros(length(K_point_list_all))

for i in 2:length(K_point_list_all)
    dist = sum( (K_point_list_all_mat[i-1,:] -     K_point_list_all_mat[i,:]).^2)
    K_point_dist_list[i] = sqrt(dist);
end
#
K_point_dist_list = cumsum(K_point_dist_list);

#K_point_name_list


K_point_tick_pos = map(x -> K_point_dist_list[x[1]],K_point_name_list)
K_point_tick_name = map(x -> x[2],K_point_name_list)
print(K_point_tick_pos," ",K_point_tick_name)

if  ( arg_input.Optional["DMFT_Jx"].go_DMFT && !arg_input.Optional["DMFT_Jx"].go_Jx)
    go_DMFT = arg_input.Optional["DMFT_Jx"].go_DMFT
    go_Jx = arg_input.Optional["DMFT_Jx"].go_Jx    
    (BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut, Tem, beta, beta_for_occ,
    DMFT_spin_type, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind,
    DMFT_loop_N, Mix_selfE,init_bias, smth_step, cal_susc,imp_dc_type, imp_ind, imp_lev_shift, imp_U, imp_J, imp_dc_n, imp_dc, imp_int_type, imp_int_parameterisation
    , Corr_ineq_orbital_Num, Solver_green_Cut ,F0, F2, F4, imp_Measure_time, imp_Thermal_time, basis_transform, green_basis, green_legendre_cutoff) = input_handler(arg_input.Optional["DMFT_Jx"]);
    
    @time DFTforge.pwork(init_variables_DMFT,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind, Corr_ineq_orbital_Num, delta); # passing basic variable for DMFT ++ , here beta_for_occ means for initial occ

    
elseif (arg_input.Optional["DMFT_Jx"].go_DMFT && arg_input.Optional["DMFT_Jx"].go_Jx)
    go_DMFT = arg_input.Optional["DMFT_Jx"].go_DMFT
    go_Jx = arg_input.Optional["DMFT_Jx"].go_Jx
    (BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut, Tem, beta, beta_for_occ,
    DMFT_spin_type, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind,
    DMFT_loop_N, Mix_selfE,init_bias, smth_step, cal_susc,imp_ind, imp_U, imp_J, imp_dc_n, imp_dc, imp_int_type, imp_int_parameterisation
    , Corr_ineq_orbital_Num, Solver_green_Cut ,F0, F2, F4, imp_Measure_time, imp_Thermal_time, basis_transform, green_basis, green_legendre_cutoff,mag_order) = input_handler(arg_input.Optional["DMFT_Jx"]);
 
elseif (!arg_input.Optional["DMFT_Jx"].go_DMFT && arg_input.Optional["DMFT_Jx"].go_Jx)
    go_DMFT = arg_input.Optional["DMFT_Jx"].go_DMFT
    go_Jx = arg_input.Optional["DMFT_Jx"].go_Jx
    (BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut,Tem, beta, DMFT_spin_type,
                 Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind, Corr_ineq_orbital_Num,imp_ind,mag_order) = input_handler(arg_input.Optional["DMFT_Jx"]);


end
# SelfE_w = init_selfE_dcbias(imp_dc,Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,iWnlist,init_bias)



#======================= input tag ===========================#

calculation_mode = "DFT-dos"  # dos, band, integrated_band, FS (kz-integrated), band_plot,    are possible

## common
filename = "Sig.out"    # real_frequency self energy file
K_grid_num = [9, 9, 3]
orb_selection = [[1:108],[5:9],[14:18],[23:27],[32:36],[41:45],[50:54],[5],[6],[7],[8]
    ,[9],[14],[15],[16],[17],[18],[23],[24],[25],[26],[27],[32],[33],[34],[35],[36]
    ,[41],[42],[43],[44],[45],[50],[51],[52],[53],[54]] # orbital indexes to calculate spectral weight
lattmu=0.0  # converged chemical potential of lattice green function
gaussian_broaden = 0.025  # gaussian broadening

targ_N = 35.8581290865675;
Num_wgrid = 2401  # if we need making w grid for DFT-dos
wmax = 15 # if we need making w grid for DFT-dos

## for pseudo optic
opt_range = 3 # for pseudo optical cond.


## for bans calculation
kpath_file = "kpath_2"  # for band or band_plot, we need kpath file
seg_grid = 30   # grid # between two high symmetry points
band_drawing_orb = 1
## for integrated band cal
kz_grid_for_integrated = 5 # grid # for integrated band along kz line


## for fermi surface

#======================= input tag ===========================#



####################### format of kpath file #############################

# G   0.00000000   0.00000000   0.00000000       X   0.00000000   0.00000000   0.50000000
# M   0.00000000   0.50000000   0.50000000       G   0.00000000   0.00000000   0.00000000
# Z   0.50000000   0.00000000   0.00000000       R   0.50000000   0.00000000   0.50000000
# A   0.50000000   0.50000000   0.50000000       Z   0.50000000   0.00000000   0.00000000
# X   0.00000000   0.00000000   0.50000000       R   0.50000000   0.00000000   0.50000000
# M   0.00000000   0.50000000   0.50000000       A   0.50000000   0.50000000   0.50000000

##########################################################################

################ MAIN CALCULATION PART ##################

if calculation_mode == "dos"
    @time DFTforge.pwork(init_variables_calmode,calculation_mode); # passing calculation mode to workerspace

    
    
    # construct non-int. hamiltonian
    Klist=kPoint_gen(K_grid_num);
    @time H_k=nonInt_H_k(Klist,hamiltonian_info,go_DMFT,go_Jx,DMFT_spin_type)
    @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker

    Total_orb_Num=sum(hamiltonian_info.scf_r.Total_NumOrbs)
    @time DFTforge.pwork(init_variables_totorb,Total_orb_Num); # passing total number of orbital to workerspace    
    @time DFTforge.pwork(passing_dc, imp_dc)
    GC.gc()
    
    
    # constructure real-frequency grid, dimension of self-energy on real-frequency
    wlist, orb_size, Num_freq, stind = w_Num_orb_Num(filename,imp_ind,DMFT_spin_type)
    Num_freq = 1500   # if we need user freq grid without Sig.out file
    wlist= wlist_gen(15,Num_freq)  # if we need user freq grid without Sig.out file
    @time DFTforge.pwork(init_variables_grid_DMFT_realf,wlist,Klist); # passing grid variable to worker for determining initial occ
    GC.gc()
    
    
    # basis transform of non-int. hamiltonian 
    if isfile_tranformMat()
        transMat=read_transforMat()
        println("Old hamiltonian")
        @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
        Trans_mat_A= Tvec_print(H0, transMat, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)


        @time H_k = Trans_H_k(H_k,Trans_mat_A,calculation_mode)
        @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
        println("New hamiltonian")
        @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
        
    else
        @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
    end
    GC.gc()
    
    
    # load self-energy on real-freq.
    selfE=init_selfE_zero_realf(Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,Num_freq)
    selfE =load_real_selfE(filename,DMFT_spin_type,imp_ind,orb_size, Num_freq, stind);
    @time DFTforge.pwork(init_variables_SelfE_w,selfE); #passing variable SelfE to worker for determining initial occ
    GC.gc()
    
    
    # calculate local green function
    G_loc_w=Gloc_gen_realf(H_k,wlist,selfE,lattmu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind,gaussian_broaden)  
    GC.gc()
    
    # calculate local spectral weight
    Aw = Aw_orb(orb_selection,wlist, G_loc_w);
    GC.gc()
    
    
    # print the spectral weight
    print_Aw(Aw,wlist)
    Aw=nothing
    GC.gc()    

    
    # calculating hybridization function
    @time InvWiess_w= Cal_InvWiess_fromGloc(wlist,G_loc_w,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,DMFT_spin_type)
    @time hyb_w=Cal_hyb_fromInvWiess(InvWiess_w,imp_ind,wlist,lattmu,H0,DMFT_spin_type,Corr_atom_equiv,Corr_orbital_Ind,delta)
    print_hyb(hyb_w,wlist,Corr_atom_equiv,DMFT_spin_type,imp_ind)
    
elseif calculation_mode == "band"
    
    
    @time DFTforge.pwork(init_variables_calmode,calculation_mode); # passing calculation mode to workerspace
    GC.gc() 
    
    
    # make band kpath grid
    kind, kpath= read_kpath(kpath_file)
    kpt, xpt,lattind= band_kpts(kpath,seg_grid)

    Total_orb_Num=sum(hamiltonian_info.scf_r.Total_NumOrbs)
    @time DFTforge.pwork(init_variables_totorb,Total_orb_Num); # passing total number of orbital to workerspace    
        
    
    # construct non-int. hamiltonian along kpath
    @time H_k=nonInt_H_k_band(kpt,hamiltonian_info,go_DMFT,go_Jx,DMFT_spin_type)
    @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
    

    
    # constructure real-frequency grid, dimension of self-energy on real-frequency
    wlist, orb_size, Num_freq, stind = w_Num_orb_Num(filename,imp_ind,DMFT_spin_type)
    @time DFTforge.pwork(init_variables_grid_DMFT_realf,wlist,kpt); # passing grid variable to worker for determining initial occ
    GC.gc()
    
    
    # basis transform of non-int. hamiltonian 
    if isfile_tranformMat()
        transMat=read_transforMat()
        println("Old hamiltonian")
        @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
        Trans_mat_A= Tvec_print(H0, transMat, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)


        @time H_k = Trans_H_k(H_k,Trans_mat_A,calculation_mode)

        @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
        println("New hamiltonian")
        @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
        
    else
        @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
    end
    GC.gc()
        
    # load self-energy on real-freq.
    selfE=init_selfE_zero_realf(Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,Num_freq)
    selfE =load_real_selfE(filename,DMFT_spin_type,imp_ind,orb_size, Num_freq, stind);
    @time DFTforge.pwork(init_variables_SelfE_w,selfE); #passing variable SelfE to worker for determining initial occ
    GC.gc()    
    
    # calculate green function for all k path and w
    G_w_k_s = G_gen_realf(H_k,wlist,selfE,lattmu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind,delta)
    GC.gc() 
    
    # calculate spectral weight
    Akw=Akw_orb(orb_selection,wlist, G_w_k_s)
    G_w_k_s = nothing
    GC.gc() 
    
    # print the spectral weight
    print_Akw(Akw,kpt,xpt,wlist,calculation_mode);
    Akw=nothing
    GC.gc()
    

elseif calculation_mode == "integrated_band"
    
    
     
    @time DFTforge.pwork(init_variables_calmode,calculation_mode); # passing calculation mode to workerspace
    GC.gc() 
    

    kind, kpath= read_kpath(kpath_file)
    int_kpt, xpt,lattind= integrated_band_kpts(kpath,seg_grid,kz_grid_for_integrated)
    kz_grid = range(0.0, stop=0.5, length =kz_grid_for_integrated)
    
    Total_orb_Num=sum(hamiltonian_info.scf_r.Total_NumOrbs)
    @time DFTforge.pwork(init_variables_totorb,Total_orb_Num); # passing total number of orbital to workerspace    
        

    # construct non-int. hamiltonian along kpath
    @time H_k=nonInt_H_k_band(int_kpt[1],hamiltonian_info,go_DMFT,go_Jx,DMFT_spin_type)
    @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
    

    # constructure real-frequency grid, dimension of self-energy on real-frequency
    wlist, orb_size, Num_freq, stind = w_Num_orb_Num(filename,imp_ind,DMFT_spin_type)


    # passing grid variable to worker
    @time DFTforge.pwork(init_variables_grid_DMFT_realf,wlist,int_kpt[1]); 
    GC.gc()



    # basis transform of non-int. hamiltonian 
    if isfile_tranformMat()
        transMat=read_transforMat()
        println("Old hamiltonian")
        @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
        Trans_mat_A= Tvec_print(H0, transMat, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)


        @time H_k = Trans_H_k(H_k,Trans_mat_A,calculation_mode)

        @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
        println("New hamiltonian")
        @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
        
    else
        @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
    end
    GC.gc()


    # load self-energy on real-freq.
    selfE=init_selfE_zero_realf(Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,Num_freq)
    selfE =load_real_selfE(filename,DMFT_spin_type,imp_ind,orb_size, Num_freq, stind);
    @time DFTforge.pwork(init_variables_SelfE_w,selfE); #passing variable SelfE to worker for determining initial occ
    GC.gc()    

    println("             \n ======= Calculating kz = 0.0  plane =======\n")
    # calculate green function for all k path and w
    G_w_k_s = G_gen_realf(H_k,wlist,selfE,lattmu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind,delta)
    GC.gc() 
    
    # calculate spectral weight
    Akw=Akw_orb(orb_selection,wlist, G_w_k_s)
    G_w_k_s = nothing
    GC.gc() 



    for kz=2:size(int_kpt)[1]
        global Akw
        # construct non-int. hamiltonian along kpath
        @time H_k=nonInt_H_k_band(int_kpt[kz],hamiltonian_info,go_DMFT,go_Jx,DMFT_spin_type)
        @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
    
    
   
        # passing grid variable to worker
        @time DFTforge.pwork(init_variables_grid_DMFT_realf,wlist,int_kpt[kz]); 
        GC.gc()

    
        # basis transform of non-int. hamiltonian 
   
        if isfile_tranformMat()
            transMat=read_transforMat()
            println("Old hamiltonian")
            @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
            Trans_mat_A= Tvec_print(H0, transMat, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)

      
            @time H_k = Trans_H_k(H_k,Trans_mat_A,calculation_mode)
            @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
            println("New hamiltonian")
            @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
        else
            @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
        end
   
        GC.gc()
    
    
        println("             \n ======= Calculating kz =",kz_grid[kz]," plane =======\n")
        # calculate green function for all k path and w
        G_w_k_s = G_gen_realf(H_k,wlist,selfE,lattmu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind,delta)
        GC.gc() 

        Akw2=Akw_orb(orb_selection,wlist, G_w_k_s)
        G_w_k_s = nothing
    
        Akw+=Akw2
        Akw2=nothing
        GC.gc()

    end
    
    # print the spectral weight
    print_Akw(Akw,kpt,xpt,wlist,calculation_mode);
    Akw=nothing
    GC.gc()
    
    
elseif calculation_mode == "FS"
    
         

    @time DFTforge.pwork(init_variables_calmode,calculation_mode); # passing calculation mode to workerspace
    GC.gc() 
    

    klist =FS_kgrid(K_grid_num)

    Total_orb_Num=sum(hamiltonian_info.scf_r.Total_NumOrbs)
    @time DFTforge.pwork(init_variables_totorb,Total_orb_Num); # passing total number of orbital to workerspace    
        

    # construct non-int. hamiltonian along kpath
    @time H_k=nonInt_H_k_band(klist[1],hamiltonian_info,go_DMFT,go_Jx,DMFT_spin_type)
    @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
    
    # constructure real-frequency grid, dimension of self-energy on real-frequency
    wlist, orb_size, Num_freq, stind = w_Num_orb_Num(filename,imp_ind,DMFT_spin_type)

    # passing grid variable to worker
    @time DFTforge.pwork(init_variables_grid_DMFT_realf,wlist,klist[1]); 
    GC.gc()




    # basis transform of non-int. hamiltonian 

    if isfile_tranformMat()
    
        transMat=read_transforMat()
    
        println("Old hamiltonian")
        @time H0 = local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
        Trans_mat_A= Tvec_print(H0, transMat, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)
        @time H_k = Trans_H_k(H_k,Trans_mat_A,calculation_mode)

        @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
        println("New hamiltonian")
        @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)

    else
        @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
    end
    GC.gc()



    # load self-energy on real-freq.
    selfE=init_selfE_zero_realf(Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,Num_freq)
    selfE =load_real_selfE(filename,DMFT_spin_type,imp_ind,orb_size, Num_freq, stind);
    zwlist, zselfE = zeroselfE(selfE,wlist)

    @time DFTforge.pwork(init_variables_grid_DMFT_realf,zwlist,klist[1]); 

    @time DFTforge.pwork(init_variables_SelfE_w,zselfE); #passing variable SelfE to worker for determining initial occ
    GC.gc()    
    GC.gc()    

    println("             \n ======= Calculating kz =",klist[1][1,3]," plane =======\n")
    # calculate green function for all k path and w
    G_w_k_s = G_gen_realf(H_k,zwlist,zselfE,lattmu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind,delta)
    GC.gc() 
    
    # calculate spectral weight
    Ak=Akw_orb_FS(orb_selection,wlist, G_w_k_s, K_grid_num)
    G_w_k_s = nothing
    GC.gc() 



    for kz=2:size(klist)[1]
        global Ak, zwlist
        # construct non-int. hamiltonian along kpath
        @time H_k=nonInt_H_k_band(klist[kz],hamiltonian_info,go_DMFT,go_Jx,DMFT_spin_type)
        GC.gc()
        @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
        GC.gc()
    
   
        # passing grid variable to worker
        @time DFTforge.pwork(init_variables_grid_DMFT_realf,zwlist,klist[kz]); 
        GC.gc()

    
        # basis transform of non-int. hamiltonian 
   
        if isfile_tranformMat()
            transMat=read_transforMat()
            println("Old hamiltonian")
            @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
            Trans_mat_A= Tvec_print(H0, transMat, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)

      
            @time H_k = Trans_H_k(H_k,Trans_mat_A,calculation_mode)
            @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
            println("New hamiltonian")
            @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
        else
            @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
        end
   
        GC.gc()
    
    
        println("             \n ======= Calculating kz =",klist[kz][1,3]," plane =======\n")
        # calculate green function for all k path and w
        G_w_k_s = G_gen_realf(H_k,zwlist,zselfE,lattmu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind,delta)
        GC.gc() 

        Ak2=Akw_orb_FS(orb_selection,wlist, G_w_k_s, K_grid_num)
        G_w_k_s = nothing
    
        Ak+=Ak2
        Ak2=nothing
        GC.gc()

    end
    
    # print the spectral weight
    xpt=0
    print_Akw(Ak,klist[1],xpt,zwlist,calculation_mode);
    Ak=nothing
    GC.gc()
        
    
    
elseif calculation_mode == "band_plot"
    
    xlist, wlist, Akw= read_Akw("Band_Akw.dat");
    Matdat = Akw[:,:,band_drawing_orb,1]'
    
    pct, delstart = percent_lowerbound(40.0,Matdat)
    
  
    kind, kpath= read_kpath(kpath_file)
    int_kpt, xpt,lattind= integrated_band_kpts(kpath,seg_grid,kz_grid_for_integrated)

    
    Plots.heatmap(xlist, wlist, Matdat./(maximum(Matdat))*10,c=cgrad([:white, :blue], [delstart*0.3]), clims=(0, 10), yrange=(-6,1), xrange=(xlist[1],xlist[end]))
    Plots.plot!(xticks=(xpt[lattind],kind),xtickfontsize=15,ytickfontsize=15)
    Plots.plot!([xpt[lattind]], seriestype="vline",color=[:black])
    




elseif calculation_mode == "DFT-dos"
    @time DFTforge.pwork(init_variables_calmode,calculation_mode); # passing calculation mode to workerspace
    
    
    # construct non-int. hamiltonian
    Klist=kPoint_gen(K_grid_num);
    @time H_k=nonInt_H_k(Klist,hamiltonian_info,go_DMFT,go_Jx,DMFT_spin_type)
    @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker

    Total_orb_Num=sum(hamiltonian_info.scf_r.Total_NumOrbs)
    @time DFTforge.pwork(init_variables_totorb,Total_orb_Num); # passing total number of orbital to workerspace    
    @time DFTforge.pwork(passing_dc, imp_dc)
    GC.gc()
    
    
    # constructure real-frequency grid, dimension of self-energy on real-frequency
    wlist, orb_size, Num_freq, stind = w_Num_orb_Num(filename,imp_ind,DMFT_spin_type)
    Num_freq = Num_wgrid   # if we need user freq grid without Sig.out file
    wlist= wlist_gen(wmax,Num_freq)  # if we need user freq grid without Sig.out file
    @time DFTforge.pwork(init_variables_grid_DMFT_realf,wlist,Klist); # passing grid variable to worker for determining initial occ
    GC.gc()
    
    lev_shift_mat = imp_level_shift(imp_lev_shift, hamiltonian_info.scf_r.Hks_R[1][1], Corr_atom_Ind, Corr_atom_equiv, imp_ind)
    
    # basis transform of non-int. hamiltonian 
    if isfile_tranformMat()
        transMat=read_transforMat()
        println("Old hamiltonian")
        @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)

        Trans_mat_A= Tvec_print(H0, transMat, Corr_atom_Ind, Corr_atom_equiv, Corr_orbital_Ind)


        @time H_k = Trans_H_k(H_k,Trans_mat_A,calculation_mode)
        @time H_k = Shift_H_k(H_k,lev_shift_mat)
        @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
        println("New hamiltonian")
        @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
        
    else
        @time H_k = Shift_H_k(H_k,lev_shift_mat)
        @time H0= local_H0(H_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,calculation_mode)
        @time DFTforge.pwork(init_variables_H_k,H_k);  #passing variable H_k to worker
    end
    GC.gc()
    
    
    # load self-energy on real-freq.
    selfE=init_selfE_zero_realf(Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,Num_freq)
    @time DFTforge.pwork(init_variables_SelfE_w,selfE); #passing variable SelfE to worker for determining initial occ
    GC.gc()
    
    
    # calculate local green function
    G_loc_w=Gloc_gen_realf(H_k,wlist,selfE,lattmu,Ineq_atom_Ind,Corr_atom_Ind,Corr_orbital_Ind,gaussian_broaden)  
    GC.gc()
    
    # calculate local spectral weight
    Aw = Aw_orb(orb_selection,wlist, G_loc_w);
    GC.gc()
    
    
    # print the spectral weight
    print_Aw(Aw,wlist)
    
    if targ_N != 0.0
        Findmu(Aw[:,1,1], targ_N, [-0.3,0.5], wlist)
    end
    
    Aw=nothing
    GC.gc()    

    
    # calculating hybridization function
    @time InvWiess_w= Cal_InvWiess_fromGloc(wlist,G_loc_w,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,DMFT_spin_type)
    @time hyb_w=Cal_hyb_fromInvWiess(InvWiess_w,imp_ind,wlist,lattmu,H0,DMFT_spin_type,Corr_atom_equiv,Corr_orbital_Ind,delta)
    print_hyb(hyb_w,wlist,Corr_atom_equiv,DMFT_spin_type,imp_ind)

end



