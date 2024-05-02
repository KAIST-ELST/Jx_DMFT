#### EDMFTF related functions

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

