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


