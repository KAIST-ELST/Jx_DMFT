import DFTforge

import Plots
import FFTW
import JSON
import Dierckx
using JSON
using ImageFiltering
using Dierckx
using Distributed
using DelimitedFiles
using Printf
using LinearAlgebra

include("subroutines/subroutines_edmftf.jl")        #### EDMFTF related functions
include("subroutines/subroutines_etc.jl")           #### etc. function
include("subroutines/subroutines_io_gen.jl")        #### processing input/ouput
include("subroutines/subroutines_dmft.jl")          #### dmft
include("subroutines/subroutines_grid.jl")          
include("subroutines/subroutines_fourier.jl")          
include("subroutines/subroutines_jx.jl")          

@everywhere using LinearAlgebra
using Distributed
using DFTforge.DFTrefinery
using DFTforge.DFTcommon;

# Julia 1.0
using Statistics
#

check_print_version()

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
hdftmpdir = ""
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
    , EDMFTF_MonteCar_step, EDMFTF_warmup_step, EDMFTF_GlobalFlip, EDMFTF_tsample, EDMFTF_nom, EDMFTF_PChangeOrder,SelfE_file) = input_handler(arg_input.Optional["DMFT_Jx"]);
    

    
elseif ((arg_input.Optional["DMFT_Jx"].Calculation_mode == "Jx-DMFT") || (arg_input.Optional["DMFT_Jx"].Calculation_mode == "Jx-DMFT-private"))

    (Calculation_mode, BoltzmannC, delta, K_grid_num, R_grid_num, iW_grid_cut, Tem, beta, beta_for_occ,
    DMFT_spin_type, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Ineq_atom_Ind,
    DMFT_loop_N, Mix_selfE, init_bias,smth_step,cal_susc,compute_EF, DMFT_solver, imp_dc_type, imp_ind, imp_lev_shift, imp_U, imp_J, imp_dc_n, imp_dc, imp_int_type, imp_int_parameterisation
    , Corr_ineq_orbital_Num, Solver_green_Cut ,F0, F2, F4, imp_Measure_time, imp_Thermal_time, basis_transform,consider_degeneracy, green_basis, green_legendre_cutoff,mag_order
    , EDMFTF_MonteCar_step, EDMFTF_warmup_step, EDMFTF_GlobalFlip, EDMFTF_tsample, EDMFTF_nom, EDMFTF_PChangeOrder, SelfE_file) = input_handler(arg_input.Optional["DMFT_Jx"]);
    

    
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
    DCM =  dcMat_noglob(H_k,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,imp_ind, imp_dc)
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
    restart_wSelfE = Check_restart_wSelfE(Ineq_atom_Ind, DMFT_solver,SelfE_file)
    if  (restart_wSelfE) # if in all imp_x folder there is output file params.obs.json, then just read self energy 
        SelfE_w_new = loadSelfE(Ineq_atom_Ind,DMFT_spin_type,SelfE_w,DMFT_solver,SelfE_file)
        
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
        SelfE_w_new = loadSelfE(Ineq_atom_Ind,DMFT_spin_type,SelfE_w,DMFT_solver,SelfE_file)
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
    println( " mu=hamiltonian_info.scf_r.ChemP ", mu ) 
    DCM =  dcMat_noglob(H_k,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,imp_ind,imp_dc)
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
        for i in Ineq_atom_Ind
            @show imp_ind[i]
        end
        DFTforge.pwork(init_variables_DMFT,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind, Corr_ineq_orbital_Num, beta, delta); # passing basic variable for DMFT ++
    end
 

    GC.gc()    


    #----- check if do restart or not -----#
    imp_totN = deepcopy(imp_dc_n)
    restart_wSelfE = Check_restart_wSelfE(Ineq_atom_Ind, DMFT_solver, SelfE_file)
    if  (restart_wSelfE) # if in all imp_x folder there is output file params.obs.json, then just read self energy 
        @show Ineq_atom_Ind
        SelfE_w_new = loadSelfE(Ineq_atom_Ind,DMFT_spin_type,SelfE_w,DMFT_solver,SelfE_file; imp_ind=imp_ind)

        imp_totN = loadOccup_from_obs_scalar(Ineq_atom_Ind,DMFT_solver);
        imp_dc = dc_change(imp_totN, imp_dc_n, dft_corrN, imp_dc_type, imp_U, imp_J) # calculate double counting
        DFTforge.pwork(passing_dc, imp_dc) # passing double counting        
    
        println("--------------------------------------")
        println("imp_dc : ", imp_dc)
        println("--------------------------------------")
       
        
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
    
    check_write_occupancy( Occup, Corr_orbital_Ind, Corr_atom_Ind )
    


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
    
    
    


if ( Calculation_mode == "Jx-DMFT-private" )

    println("#####################################################")
    println(" Calculation type : Jx calculation with self-energy  ")
    println("#####################################################\n\n\n")

    
    Rlist=RPoint_gen(R_grid_num);
    Klist=kPoint_gen(K_grid_num);  
    iWnlist=Matsubara_gen(iW_grid_cut,beta);
    @time DFTforge.pwork(init_variables_grid_DMFT,iWnlist,Rlist,Klist); # passing grid variable to worker
    
    ################  non-interacting hamiltonian construction ################ 
    @time H_k=nonInt_H_k(Klist,hamiltonian_info,Calculation_mode,DMFT_spin_type)
    mu=hamiltonian_info.scf_r.ChemP
    println( " mu=hamiltonian_info.scf_r.ChemP ", mu ) 
    DCM =  dcMat_noglob(H_k,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind,imp_ind,imp_dc)
    H0_k = DC_shift(H_k,DCM)
    @time DFTforge.pwork(init_variables_H_k,H0_k);  #passing variable H_k to worker
    println("")
    
    
    ## DMFT part input
    @time DFTforge.pwork(init_variables_Jx,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind,beta, Corr_ineq_orbital_Num, delta); # passing basic variable for DMFT ++ , here beta_for_occ means for initial occ
    @time DFTforge.pwork(passing_dc, imp_dc)
    
    ################ GRID GENERATION PART ##################

    atom_orbitals=atom_orb_list(hamiltonian_info);
    Total_orb_Num=sum(hamiltonian_info.scf_r.Total_NumOrbs)
    @time DFTforge.pwork(init_variables_totorb,Total_orb_Num); # passing total number of orbital to workerspace


    
    
    ################ Determination of Nele, DFTcorr_ele ##################
    ##  create real frequency grid & selfenergy in real frequency
    wnum = 3000
    wlist= wlist_gen(15,wnum)
    w0ind=convert(Int64,wnum/2)
    SelfE_realw = init_selfE_dc(imp_dc,Ineq_atom_Ind,DMFT_spin_type,Corr_ineq_orbital_Num,hamiltonian_info,wlist)
    lattmu=mu  # 0.0
    gaussian_broaden=0.025
    GC.gc()  
            
    @time DFTforge.pwork(init_variables_grid_DMFT_realf,wlist,Klist); # passing grid variable to worker for determining initial occ
    @time DFTforge.pwork(init_variables_SelfE_w,SelfE_realw); #passing variable SelfE to worker for determining initial occ
    GC.gc()
    
    #--- Check if the file "Nele.dat" exist or not. If there is no "Nele.dat", calculate Nele at low temperature (T=100K)
    if (  !isfile("Nele.dat") || !isfile("transform_Mat.dat") || !isfile("DFT_CorrNele.dat"))
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
    
    SelfE_realw = nothing;
    wlist = nothing;
    G_loc_w = nothing;
    Aw = nothing;
    InvWiess_w0 = nothing;
    hyb_w0 = nothing;
    norm_implev = nothing;
    GC.gc() 
    
    
    
   
    



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

    
    

    
    @time H_loc= local_H0(H0_k,imp_ind,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind)
    #if "consider_degeneracy" turned on, the generate states would be treated equivalently 
    if consider_degeneracy == true
        imp_ind = imp_ind_trans_by_degeneracy(imp_ind,H_loc,Corr_atom_Ind,Corr_atom_equiv,Corr_orbital_Ind)
        Corr_ineq_orbital_Num = cal_corr_ineq_orbital_N(imp_ind)
        DFTforge.pwork(init_variables_DMFT,Corr_atom_Ind,Corr_orbital_Ind, Corr_atom_equiv,imp_ind, Corr_ineq_orbital_Num, beta, delta); # passing basic variable for DMFT ++
    end
 


    GC.gc()    


    #----- check if do restart or not -----#
    imp_totN = deepcopy(imp_dc_n)
    restart_wSelfE = Check_restart_wSelfE(Ineq_atom_Ind, DMFT_solver,SelfE_file)
    if  (restart_wSelfE) # if in all imp_x folder there is output file params.obs.json, then just read self energy 
        SelfE_w_new = loadSelfE(Ineq_atom_Ind,DMFT_spin_type,SelfE_w,DMFT_solver,SelfE_file)

        imp_dc = dc_change(imp_totN, imp_dc_n, dft_corrN, imp_dc_type, imp_U, imp_J) # calculate double counting
        DFTforge.pwork(passing_dc, imp_dc) # passing double counting        
     
        println("--------------------------------------")
        println("imp_dc : ", imp_dc)
        println("--------------------------------------")       
       
	            
	SelfE_w_sm = deepcopy(SelfE_w_new)
        
	#if go_bareSelfE(smth_step,Ineq_atom_Ind,DMFT_solver)
        #    SelfE_w_sm = deepcopy(SelfE_w_new)
        #else
        #    SelfE_w_sm =SelfE_smoothing(Corr_atom_equiv,Tem,iWnlist,SelfE_w_new);
        #end
        
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
    mulog =[] 
    init_mu = mu_guess(1, mulog, DMFT_spin_type)
    @time begin
        print("          ===== Finding mu ... =====\n")
        mu, Occup, H_loc =Find_mu(iWnlist,H0_k,SelfE_w_sm,targetN,beta,delta,Corr_atom_Ind,Corr_orbital_Ind,Corr_atom_equiv,DMFT_spin_type,init_mu)
        println("")
    end
    

    print("-------------------------\n")
    print("Initial latt. mu :", mu,"\n") 
    print("-------------------------\n")
    
    check_write_occupancy( Occup, Corr_orbital_Ind, Corr_atom_Ind )
    
    


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
    
    
    check_write_occupancy( Occup, Corr_orbital_Ind, Corr_atom_Ind )
    
    

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
    for (i,j) in atom12_list
        filen = string("Jx_raw_",i,"_",j,".dat")
        push!(readfilelist, filen)
        push!(check_Box, isfile(filen))
        if !isfile(filen)
            push!(nexist_file, filen) 
        end
    end
    #for i in Corr_atom_Ind
    #    for j in Corr_atom_Ind
    #            
    #        filen = string("Jx_raw_",i,"_",j,".dat")
    #        push!(readfilelist, filen)
    #        push!(check_Box, isfile(filen))
    #        if !isfile(filen)
    #            push!(nexist_file, filen) 
    #        end
    #    end
    #end
    
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

    n_Corr_atom_Ind = size(Corr_atom_Ind)[1]
    println( "n_Corr_atom_Ind = ", size(Corr_atom_Ind)[1] )
    
    
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
                    
                    if abs(Jval)<1e-5
                        Jval = 0.
                    end
                    if matind == matind2
                    #if true
                        aind1 = findall(x->x == matind[1], Corr_atom_Ind)
                        aind2 = findall(x->x == matind[2], Corr_atom_Ind)
                        
                        
                        dist[Rind[1]+1,Rind[2]+1,Rind[3]+1,aind1[1],aind2[1]] = dist_dum
                        #dist[Rind[1]+1,Rind[2]+1,Rind[3]+1,aind2[1],aind1[1]] = dist_dum
                        
                        Jmat[Rind[1]+1,Rind[2]+1,Rind[3]+1,aind1[1],aind2[1]] = Jval
                        #Jmat[Rind[1]+1,Rind[2]+1,Rind[3]+1,aind2[1],aind1[1]] = Jval
                        if ( aind1[1] == aind2[1] )
                            for iatom = 1:n_Corr_atom_Ind
                                Jmat[Rind[1]+1,Rind[2]+1,Rind[3]+1,iatom,iatom] = Jval
                            end
                        end

                        cellvec[Rind[1]+1,Rind[2]+1,Rind[3]+1,:] = cellind

                        rvec = cellind[1]*cell_vec[1,:] + cellind[2]*cell_vec[2,:] + cellind[3]*cell_vec[3,:]
                        println( filename, ": ", cellind, " ", Rind, " ", aind1, " ", aind2, " ", dist_dum, " ", Jval, " ; ", rvec  )

            
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
    
    println( "\ndist 11 :: " )
    println( dist[:,:,:,1,1] )
    if ( n_Corr_atom_Ind>1 ) 
        println( "\ndist 12 :: " )
        println( dist[:,:,:,1,2] )
    end

    #println( "\nJmat 11 :: " )
    #println( Jmat[:,:,:,1,1] )
    if ( n_Corr_atom_Ind>1 ) 
        #println( "\nJmat 22 :: " )
        #println( Jmat[:,:,:,2,2] )
        println( "\nJmat 12 :: " )
        println( Jmat[:,:,:,1,2] )
        println( "\nJmat 21 :: " )
        println( Jmat[:,:,:,2,1] )
        println( "" )
    
        println( "\nJmat 12 11n :: " )
        println( Jmat[1,1,:,1,2] )
    end
    
    println( "\nJmat 11 12n :: " )
    println( Jmat[1,2,:,1,1] )
    println( dist[1,2,:,1,1] )
    
    return ordered_dist, dist, cellvec, Jmat
end

function get_cellvec_minus(cellvec, R_grid_num, Jmat)
    cellvec_minus_ind = zeros(Int,R_grid_num[1],R_grid_num[2],R_grid_num[3],3)
    for ra1 = 1:R_grid_num[1]
        for ra2 = 1:R_grid_num[2]
            for ra3 = 1:R_grid_num[3]
                for rb1 = 1:R_grid_num[1]
                    for rb2 = 1:R_grid_num[2]
                        for rb3 = 1:R_grid_num[3]
                            rab = cellvec[ra1, ra2, ra3, :] + cellvec[rb1, rb2, rb3, :] 
                            if norm(rab)< 1e-5
                                cellvec_minus_ind[ra1,ra2,ra3,:] = [rb1,rb2,rb3]
                            end
                        end
                    end
                end
            end
        end
    end

    #for ra = Iterators.product(1:R_grid_num[1],1:R_grid_num[2],1:R_grid_num[3])
    #    ra_minus = cellvec_minus_ind[ra[1],ra[2],ra[3],:]
    #    #println( ra )
    #    #println( ra_minus )
    #    if sum(ra_minus) >= 3 
    #        println( ra, ra_minus, " r=", cellvec[ra[1],ra[2],ra[3],:],  ", -r=", cellvec[ra_minus[1],ra_minus[2],ra_minus[3],:],
    #                         ", J11(r)=",Jmat[ra[1],ra[2],ra[3],1,1], ", J12(r)=",Jmat[ra[1],ra[2],ra[3],1,2], 
    #                         ", J11(-r)=",Jmat[ra_minus[1],ra_minus[2],ra_minus[3],1,1], ", J12(-r)=",Jmat[ra_minus[1],ra_minus[2],ra_minus[3],1,2] )
    #    else
    #        println( ra, ra_minus, " r=", cellvec[ra[1],ra[2],ra[3],:],
    #                         ", J11(r)=",Jmat[ra[1],ra[2],ra[3],1,1], ", J12(r)=",Jmat[ra[1],ra[2],ra[3],1,2] )
    #    end
    #end
    return cellvec_minus_ind
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

function print_kpath_cart( kpath, r_vec )
    println("kpath in cart :")
    kpath_r = kpath_in_reciprocal(kpath,r_vec)
    for i=1:size(kpath_r)[1]
        println( kpath_r[i] ) 
    end
    println("")
end
    
function kpath_in_reciprocal(kpath,r_vec)
    kpath_r =[]
    for i=1:size(kpath)[1]
        kpath_dum = kpath[i][1].*r_vec[1,:] + kpath[i][2].*r_vec[2,:] + kpath[i][3].*r_vec[3,:]
    
        push!(kpath_r,kpath_dum)
    end
    return kpath_r
end


function band_kpts2(kpath,totgrid, r_vec)
    
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
    
    #qvec = qpt_frac[qind,1]*r_vec[1,:] + qpt_frac[qind,2]*r_vec[2,:] + qpt_frac[qind,3]*r_vec[3,:]
    kpath_r = kpath_in_reciprocal(kpath,r_vec)
    for i =1:seg_num
        length_ = ((kpath_r[i*2][1]-kpath_r[i*2-1][1])^2+(kpath_r[i*2][2]-kpath_r[i*2-1][2])^2+(kpath_r[i*2][3]-kpath_r[i*2-1][3])^2)^(1/2)
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


function Cal_MagMat(qmode, Lattconst, qpt_frac, cell_vec, r_vec, rel_pos_frac, R_grid_num, cellvec_fac, Corr_atom_Ind, Corr_atom_equiv, Moment, Neighbors_cut, ordered_dist, dist, Jmat, cellvec_minus_ind)
    Orbind_num = [size(Jmat)[4],  size(Jmat)[5]]
    Magmat_ = zeros(ComplexF64, size(qpt_frac)[1], Orbind_num[1], Orbind_num[2])
    considered_Neigh = zeros(Int64, size(Corr_atom_Ind)[1], Neighbors_cut)
  
    sign = Corr_atom_equiv./(abs.(Corr_atom_equiv))
    println( "Moment sign from 'Corr_atom_equiv' : ", sign )
    println( "Neighbors_cut : ", Neighbors_cut ) 
    println( "odist : ", ordered_dist[1] )
    println( "Neighbors_cut : ", Neighbors_cut,  ", ", size(ordered_dist[1]) )

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
 
                            if ( 0.0 < dist[R1,R2,R3,o1,o2] <= 20. )
                                #println( Neighbors_cut, "dist odist : ", dist[R1,R2,R3,o1,o2], ", ", ordered_dist[o1][Neighbors_cut+1]+0.00001 )
                                for n = 1:Neighbors_cut
                                    if ( abs(dist[R1,R2,R3,o1,o2]-ordered_dist[o1][n+1]) <0.00001 && (qind ==1) )
                                        considered_Neigh[o1,n]+=1
                                    end
                                end
                                #cvec = cellvec_fac[R1,R2,R3,1]*cell_vec[1,:] + cellvec_fac[R1,R2,R3,2]*cell_vec[2,:] + cellvec_fac[R1,R2,R3,3]*cell_vec[3,:]
                                cvec = (cellvec_fac[R1,R2,R3,1]+rel_pos_frac[o1,o2,1])*cell_vec[1,:] + (cellvec_fac[R1,R2,R3,2]+rel_pos_frac[o1,o2,2])*cell_vec[2,:] + (cellvec_fac[R1,R2,R3,3]+rel_pos_frac[o1,o2,3])*cell_vec[3,:]
                                #println( "R:", [R1,R2,R3], ", J:", Jmat[R1,R2,R3,o1,o2] )
                                #error("It is using wrong indices.")
    
                                Magmat_[qind, o1,o1] += (2.0/Moment[o1])*(Jmat[R1,R2,R3,o1,o2])*sign[o1]
                                Magmat_[qind, o2,o2] += (2.0/Moment[o2])*(Jmat[R1,R2,R3,o1,o2])*sign[o2]
                                # """Magmat_[qind, o1,o2] += (1.0/Moment[o1])*(Jmat[R1,R2,R3,o1,o2])*sign[o1]*exp( 2*pi*dot(qpoint[qind,:], cellvec[R1,R2,R3,:])*im )"""
                                if qmode == 1
                                    Magmat_[qind, o1,o2] -= (2.0/Moment[o2])*(Jmat[R1,R2,R3,o1,o2])*sign[o2]*exp( dot(qvec, -cvec)*im )
                                    Magmat_[qind, o2,o1] -= (2.0/Moment[o1])*(Jmat[R1,R2,R3,o1,o2])*sign[o1]*exp( dot(qvec,  cvec)*im )
                                elseif qmode == 2 
                                    Magmat_[qind, o1,o2] -= (2.0/Moment[o2])*(Jmat[R1,R2,R3,o1,o2])*sign[o2]*exp( 2*pi/Lattconst*dot(qpt_frac[qind,:], -cvec)*im )
                                    Magmat_[qind, o2,o1] -= (2.0/Moment[o1])*(Jmat[R1,R2,R3,o1,o2])*sign[o1]*exp( 2*pi/Lattconst*dot(qpt_frac[qind,:],  cvec)*im )
                                end
                            end

                            
                        end
                    end
                end
            end
        end
        
        Magnon_disp[qind,:] = eigvals(Magmat_[qind,:,:])
    end

    println( Magmat_[50,:,:] )
            
    
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

function unique_magnon(Magnon_disp)
 

    uniq_magnon_dum = []
    push!(uniq_magnon_dum, abs.(Magnon_disp[:,1]) )
    for i =1:size(Magnon_disp)[2]
        for j =1:size(uniq_magnon_dum)[1]
            iden_criteria = abs.(abs.(uniq_magnon_dum[j]) .- abs.(Magnon_disp[:,i]))
            ind = findall(x->x<0.0001, iden_criteria)
            if (size(ind)[1] != size(Magnon_disp)[1] )
                push!(uniq_magnon_dum, abs.(Magnon_disp[:,i]))
            end
        end
    end
    
    uniq_magnon = zeros(size(Magnon_disp)[1],size(uniq_magnon_dum)[1])
    for i =1:size(uniq_magnon)[2]
        uniq_magnon[:,i] = uniq_magnon_dum[i][:]
    end
    
    return uniq_magnon
end


function print_magnon(xpt, Magnon_disp, lattind, kindex)
    max_ene_arr =[]
    #uniq_magnon = unique_magnon(Magnon_disp)
    uniq_magnon = Magnon_disp
    
    open("magnon.dat", "w") do io
        @printf(io,"    %14s", "kpt")
        for i=1:size(uniq_magnon)[2]
            push!(max_ene_arr, maximum(abs.(real(uniq_magnon[:,i]))) )
            tag = string("disp_",i,)
            if i !=size(uniq_magnon)[2]
                @printf(io,"    %14s", tag)
            else
                @printf(io,"    %14s\n", tag)
            end            
        end
        
        for i=1:size(uniq_magnon)[1]
            @printf(io,"    %14.8f", xpt[i])
            for j=1:size(uniq_magnon)[2]
                tag = string("disp_",j,)
                if j !=size(uniq_magnon)[2]
                    @printf(io,"    %14.8f", real(uniq_magnon[i,j]))
                else
                    @printf(io,"    %14.8f\n", real(uniq_magnon[i,j]))
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
        
        for i=1:size(uniq_magnon)[2]
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
   
    """
        cell_vec 

    Lattice vector in Cartesian

        r_vec

    Reciprocal vector in Cartesian
    """
    print_c_r_vec(cell_vec, r_vec)
    
    if( isfile("occ.log") )
        Moment_size = loadmomentsize("occ.log", Corr_atom_equiv)  #magnetic moment size, unit :Bohr magneton
    else
        Moment_size = abs.( readdlm("occu_mag_tblatt.dat")[2,:] )
    end
    println("Magnetic moment size :", Moment_size, " (unit : Bohr magneton)")
    println("")

    frac_coord = hamiltonian_info.scf_r.Gxyz*inv(hamiltonian_info.scf_r.tv)
    rel_pos_frac = relative_position(frac_coord)
    println("Frac coordinate of atoms :")
    println( frac_coord )
    println("")
    println("Relative frac coordinate between (correlated) atoms :")
    for (o1,o1val) in enumerate(Corr_atom_Ind)
        for (o2,o2val) in enumerate(Corr_atom_Ind)
            println( o1val, o2val, " :: ", rel_pos_frac[o1val,o2val,:] )
        end
    end
    println("")

    readfilelist = check_isradata(atom12_list,Corr_atom_Ind)
    println("Loading Jmat from following file list ...")
    println("    :", readfilelist)
    println("")

    """
        cellvec

    Coefficients of the lattice vectors (or cell vector integers) in a Jmat data.
    """
    ordered_dist, dist, cellvec, Jmat = Jmat_construct(readfilelist,Corr_atom_Ind,R_grid_num)
    
    tot_grid = 1000   # grid # between two high symmetry points
    println("Magnon k-path w/ grid num :",tot_grid)
    kindex, kpath= read_kpath("kpath")
    print_kpath_cart(kpath,r_vec)
    kpt_frac, xpt,lattind= band_kpts2(kpath,tot_grid, r_vec)
    
    
    println("Getting minus cell vectors..")
    cellvec_minus_ind = get_cellvec_minus( cellvec, R_grid_num , Jmat)

    #for dum_obj in (qmode, Lattconst, kpt_frac, cell_vec, r_vec, rel_pos_frac, R_grid_num, cellvec, Corr_atom_Ind, Corr_atom_equiv, Moment_size, Neighbors_cut, ordered_dist, dist, Jmat)
    #    println( "########################" )
    #    println( "TYPE : ", typeof(dum_obj) )
    #    println( "SIZE : ", size(dum_obj) )
    #    println( "DATA : \n", dum_obj )
    #end
    #exit(1)
    considered_Neigh, Magmat, Magnon_disp =Cal_MagMat(qmode, Lattconst, kpt_frac, cell_vec, r_vec, rel_pos_frac, R_grid_num, cellvec, Corr_atom_Ind, Corr_atom_equiv, Moment_size, Neighbors_cut, ordered_dist, dist, Jmat, cellvec_minus_ind)
    println("Magmat size : ", size(Magmat) )
    println("Magnon_disp size : ", size(Magnon_disp) )
    println()
    
    println("Considered Neighbors : ")
    for i=1:size(considered_Neigh)[1]
        ind = findall(x->x>0, considered_Neigh[i,:])
        println("  # of 1st NN, 2nd NN ... :", considered_Neigh[i,ind])
    end
    println()
    
    print_magnon(xpt, Magnon_disp, lattind, kindex)
    
end

#Plots.plot(xpt,real(Magnon_disp[:,1]./2), width =3, yrange =[0])
#Plots.plot!([xpt[lattind]], seriestype="vline")
#Plots.plot!(xpt,real(Magnon_disp[:,2]./2), width =3, color = "red")


