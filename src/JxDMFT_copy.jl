SCF_DMFT_dir = map(x->string(x), ARGS)[1]
println("SCF-DMFT dir : ", SCF_DMFT_dir)
println("")

cpfdir = string(pwd(),"/", SCF_DMFT_dir)




for i = 1:10
    impdir = string(cpfdir,"/imp_",i)
    if isdir(impdir)

        pwdimpdir = string( pwd(), "/imp_",i)
        if !isdir(pwdimpdir)
            mkdir(pwdimpdir)
        end
 
        fname = "params.obs.json"
        if isfile( string(impdir,"/",fname) )
            cp(string(impdir,"/",fname), string(pwdimpdir,"/",fname),force=true)
            println("cp ", SCF_DMFT_dir,"/imp_",i,"/",fname,"  ./imp_",i,"/",fname)

            fnlist = filter!(x-> occursin(r"params.obs",x), readdir(impdir))
            
            for fn in fnlist
                cp(string(impdir,"/",fn), string(pwdimpdir,"/",fn),force=true)
            end         
            println("cp ", SCF_DMFT_dir,"/imp_",i,"/params.obs*  ./imp_",i,"/params.obs*")

        end

        fname = "Sig.out"
        if isfile( string(impdir,"/",fname) )
            cp(string(impdir,"/",fname), string(pwdimpdir,"/",fname),force=true)
            println("cp ", SCF_DMFT_dir,"/imp_",i,"/",fname,"  ./imp_",i,"/",fname)


            fnlist = filter!(x-> occursin(r"Sig",x), readdir(impdir))                
            for fn in fnlist
                cp(string(impdir,"/",fn), string(pwdimpdir,"/",fn),force=true)
            end
            println("cp ", SCF_DMFT_dir,"/imp_",i,"/Sig*  ./imp_",i,"/Sig*")

        end

        fname = "Gf.out"
        if isfile( string(impdir,"/",fname) )
            cp(string(impdir,"/",fname), string(pwdimpdir,"/",fname),force=true)
            println("cp ", SCF_DMFT_dir,"/imp_",i,"/",fname,"  ./imp_",i,"/",fname)


            fnlist = filter!(x-> occursin(r"Gf",x), readdir(impdir))                 
            for fn in fnlist
                cp(string(impdir,"/",fn), string(pwdimpdir,"/",fn),force=true)
            end
            println("cp ", SCF_DMFT_dir,"/imp_",i,"/Gf*  ./imp_",i,"/Gf*")

        end

    end
end


fname = "Jx_DMFT.jl"
cp(string(cpfdir,"/",fname), string(pwd(),"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)




fname = filter!(x-> occursin(r"toml",x), readdir(SCF_DMFT_dir))[1]           
cp(string(cpfdir,"/",fname), string(pwd(),"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)


fname = "wannier.win"
cp(string(cpfdir,"/",fname), string(pwd(),"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)


fname = "wannier_hr.dat"
cp(string(cpfdir,"/",fname), string(pwd(),"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)


fname = "Nele.dat"
cp(string(cpfdir,"/",fname), string(pwd(),"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)

fname = "DFT_CorrNele.dat"
cp(string(cpfdir,"/",fname), string(pwd(),"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)

fname = "transform_Mat.dat"
cp(string(cpfdir,"/",fname), string(pwd(),"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)


fname = "scf.log"
cp(string(cpfdir,"/",fname), string(pwd(),"/",fname),force=true)
println("cp ", SCF_DMFT_dir,"/",fname,"  ./",fname)







