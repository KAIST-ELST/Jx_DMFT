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

