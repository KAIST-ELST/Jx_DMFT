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



function loadSelfE(Ineq_atom_Ind,DMFT_spin_type,SelfE_w,DMFT_solver,SelfE_file; imp_ind=false )
    println( "Loading self-energy :: 'loadSelfE()'..." )
    SelfE_w_dum = deepcopy(SelfE_w)
    
    for i in Ineq_atom_Ind
        ineq_orbital_Num = size(SelfE_w[i])[2]
        if DMFT_spin_type > 0
            for s = 1:DMFT_spin_type
                for l=1:ineq_orbital_Num
                    @show ( ineq_orbital_Num , (s-1) , l)
                end
            end
        end
        if imp_ind != false
            println( "imp_ind ::" )
            println( imp_ind )
        end
        orbital_indices    = imp_ind[i]
        ineq_orbital_Num_red = maximum(orbital_indices)
        for l=1:ineq_orbital_Num
            @show ineq_orbital_Num_red, l, orbital_indices[l,l], ineq_orbital_Num_red + orbital_indices[l,l]
        end
    end

    if SelfE_file != ""
            
        for i in Ineq_atom_Ind
            
            str = string(SelfE_file,"_",i,".dat")
            println("Reading SelfE_w ..., [",str,"]")            
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

    else
        if DMFT_solver == "ComCTQMC"
            for i in Ineq_atom_Ind

                print("Reading SelfE_w (ComCTQMC) ... imp[",i,"]\n")
                str =string("./imp_",i,"/params.obs.json")
                obser=readJSON(str)
           
                ineq_orbital_Num = size(SelfE_w[i])[2]
                if imp_ind == false 
                    orbital_indices    = collect(1:ineq_orbital_Num)
                else
                    orbital_indices    = imp_ind[i]
                end
                ineq_orbital_Num_red = maximum(orbital_indices)
                if DMFT_spin_type > 0
                    for s = 1:DMFT_spin_type
                        for l=1:ineq_orbital_Num
                            SelfE_w_dum[i][:,l,s]=obser["partition"]["self-energy"][string( (s-1)*ineq_orbital_Num_red+orbital_indices[l,l] )]["function"]["real"].+obser["partition"]["self-energy"][string( (s-1)*ineq_orbital_Num_red+orbital_indices[l,l] )]["function"]["imag"]*im  # lattice selfE
                        end
                    end
                end
            end
        elseif DMFT_solver == "EDMFTF_ctqmc"
            for i in Ineq_atom_Ind
                print("Reading SelfE_w (EDMFTF_ctqmc) ... imp[",i,"]\n")
                str = string("./imp_",i,"/Sig.out")
                open(str, "r") do f
                    Nstring=readlines(f)
                    linenum = size(Nstring)[1]

                    ineq_orbital_Num = size(SelfE_w[i])[2]

                    @show (DMFT_spin_type, ineq_orbital_Num, linenum-1)
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



function Check_restart_wSelfE(Ineq_atom_Ind, DMFT_solver, SelfE_file)
    
    true_table=[]    
    if (SelfE_file != "")
        for i in Ineq_atom_Ind
            str = string(SelfE_file,"_",i,".dat")
            
            if isfile(str)
                push!(true_table,true)
            else
                push!(true_table,false)
            end
        end                
     
                
    else
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




function dcMat_noglob(H_k,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,imp_ind,imp_dc)
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



function check_print_version()
    X_VERSION = VersionNumber("0.9.3-pub+20190529");
    println(" JX_VERSION: ",X_VERSION)
    println(" Visit https://kaist-elst.github.io/DFTforge.jl/ for details & updates ")
    println(" Tested with Julia v1.0 and v1.1 which the most recent version of Julia in 201906 https://julialang.org/")
end

function check_write_occupancy( occu_raw, corr_orbital_ind, corr_atom_ind; foutname="occu_mag_tblatt.dat" ) 
    Occu_orb_up     =  diag( real( occu_raw[:,:,1] ) )
    Occu_orb_dn     =  diag( real( occu_raw[:,:,2] ) ) 
    Occu_up     =  sum(Occu_orb_up)
    Occu_dn     =  sum(Occu_orb_dn)
    Occu_tot    =  Occu_up + Occu_dn
    
    println("Up-spin occup :", Occu_orb_up)
    println("")
    println("dn-spin occup :", Occu_orb_dn)
    println("")
    println("****************************************")
    println("Total charge #:", Occu_tot )
    println("****************************************")

    
    println("Corr orb ind for atoms: ", corr_orbital_ind )

    Occu_orb_up_corr    =  Occu_orb_up[ [(corr_orbital_ind...)...] ]
    Occu_orb_dn_corr    =  Occu_orb_dn[ [(corr_orbital_ind...)...] ]
    Occu_orb_tot_corr   =  Occu_orb_up_corr + Occu_orb_dn_corr
    Occu_orb_diff_corr  =  Occu_orb_up_corr - Occu_orb_dn_corr
    Occu_up_corr    = sum(Occu_orb_up_corr  )
    Occu_dn_corr    = sum(Occu_orb_dn_corr  )
    Occu_tot_corr   = sum(Occu_orb_tot_corr )
    Occu_diff_corr  = sum(Occu_orb_diff_corr)
    println("Occ corr orb up : ", Occu_orb_up_corr)
    println("Occ corr orb dn : ", Occu_orb_dn_corr)
    println("Occ corr orb tot : ", Occu_orb_tot_corr)
    println("Occ corr orb up-dn : ", Occu_orb_up_corr - Occu_orb_dn_corr)

    Occu_atom_tot_corr   = []
    Occu_atom_diff_corr  = []
    for iatom in corr_atom_ind
        push!( Occu_atom_tot_corr  , sum(Occu_orb_up[ corr_orbital_ind[iatom] ] +  Occu_orb_dn[ corr_orbital_ind[iatom] ]) )
        push!( Occu_atom_diff_corr , sum(Occu_orb_up[ corr_orbital_ind[iatom] ] -  Occu_orb_dn[ corr_orbital_ind[iatom] ]) )
    end
    println("Occu atoms  :", Occu_atom_tot_corr )
    println("Up-dn atoms :", Occu_atom_diff_corr)
    writedlm( foutname , [ Occu_atom_tot_corr , Occu_atom_diff_corr ] )

    println("****************************************")
    println("Charge for correlated orbitals :", Occu_tot_corr)
    println("Spin-polrization for correlated orbitals :", 0.5*Occu_diff_corr)
    println("****************************************")
    println("")
end
