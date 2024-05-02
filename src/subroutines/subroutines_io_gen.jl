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
        SelfE_file = DMFT_Jx_option.SelfE_file


        return (Calculation_mode, BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut, Tem, beta, beta_for_occ,
                DMFT_spin_type, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind,
                DMFT_loop_N, Mix_selfE, init_bias,smth_step,cal_susc,compute_EF, DMFT_solver, imp_dc_type, imp_ind, imp_lev_shift, imp_U, imp_J, imp_dc_n, imp_dc, imp_int_type, imp_int_parameterisation
                , Corr_ineq_orbital_Num, Solver_green_Cut ,F0, F2, F4, imp_Measure_time, imp_Thermal_time, basis_transform,consider_degeneracy, green_basis, green_legendre_cutoff
                ,EDMFTF_MonteCar_step, EDMFTF_warmup_step, EDMFTF_GlobalFlip, EDMFTF_tsample, EDMFTF_nom, EDMFTF_PChangeOrder, SelfE_file );
        
        
    elseif ((DMFT_Jx_option.Calculation_mode == "Jx-DMFT") || (DMFT_Jx_option.Calculation_mode == "Jx-DMFT-private"))

        
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
        SelfE_file = DMFT_Jx_option.SelfE_file
        
        return (Calculation_mode,BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut, Tem, beta, beta_for_occ,
                DMFT_spin_type, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind,
                DMFT_loop_N, Mix_selfE, init_bias,smth_step,cal_susc, compute_EF, DMFT_solver, imp_dc_type, imp_ind, imp_lev_shift, imp_U, imp_J, imp_dc_n, imp_dc, imp_int_type, imp_int_parameterisation
                , Corr_ineq_orbital_Num, Solver_green_Cut ,F0, F2, F4, imp_Measure_time, imp_Thermal_time, basis_transform,consider_degeneracy, green_basis, green_legendre_cutoff,mag_order
                ,EDMFTF_MonteCar_step, EDMFTF_warmup_step, EDMFTF_GlobalFlip, EDMFTF_tsample, EDMFTF_nom, EDMFTF_PChangeOrder,SelfE_file );

        
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


function readJSON(jsonfile)
    ob=Dict()
    open(jsonfile, "r") do f
        a=read(f,String)
        ob=JSON.parse(a)
    end
    
    return ob
end

