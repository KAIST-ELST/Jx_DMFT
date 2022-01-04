Jx_DMFT_dir = map(x->string(x), ARGS)[1]
println("")
println("Jx-DMFT dir : ", Jx_DMFT_dir)
println("")

cpfdir = string(pwd(),"/", Jx_DMFT_dir)






fname = "Jx_DMFT.jl"
cp(string(cpfdir,"/",fname), string(pwd(),"/",fname),force=true)
println("cp ", Jx_DMFT_dir,"/",fname,"  ./",fname)


fname = filter!(x-> occursin(r"toml",x), readdir(Jx_DMFT_dir))[1]           
cp(string(cpfdir,"/",fname), string(pwd(),"/",fname),force=true)
println("cp ", Jx_DMFT_dir,"/",fname,"  ./",fname)


fname = "wannier.win"
cp(string(cpfdir,"/",fname), string(pwd(),"/",fname),force=true)
println("cp ", Jx_DMFT_dir,"/",fname,"  ./",fname)


fname = "wannier_hr.dat"
cp(string(cpfdir,"/",fname), string(pwd(),"/",fname),force=true)
println("cp ", Jx_DMFT_dir,"/",fname,"  ./",fname)

fname = "occ.log"
cp(string(cpfdir,"/",fname), string(pwd(),"/",fname),force=true)
println("cp ", Jx_DMFT_dir,"/",fname,"  ./",fname)


fname = filter!(x-> occursin(r"raw",x), readdir(Jx_DMFT_dir))
#println(fname)
for i =1:size(fname)[1]
    println("cp ", Jx_DMFT_dir,"/",fname[i],"  ./",fname[i])
    cp(string(cpfdir,"/",fname[i]), string(pwd(),"/",fname[i]),force=true)
end

println("")
println("-------------------------------------------------------------------------------")
println("To calculate magnon dispersion, change *.toml file 'Calculation_mode = \"Magnon\"")
println("and write a file \"kpath\" in the following format :")
println("")
println("                 G 0.00  0.00  0.00    N 0.00  0.00  0.50")
println("                 N 0.00  0.00  0.50    P 0.25  0.25  0.25")
println("                 P 0.25  0.25  0.25    G 0.00  0.00  0.00")
println("                 G 0.00  0.00  0.00    H 0.50 -0.50  0.50")
println("                 H 0.50 -0.50  0.50    N 0.00  0.00  0.50")
println("-------------------------------------------------------------------------------")
println("")
println("")
