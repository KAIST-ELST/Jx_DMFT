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


