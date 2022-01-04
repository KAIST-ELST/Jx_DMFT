set xra [0:20]
plot "zeroline.dat" u 1:2 w l lw 2 lt 0,"Jx_1_1.dat" u 1:2 w linespoints pt 7 ps 1.7 lw 1.3, "Jx_1_2.dat" u 1:2 w linespoints pt 7 ps 1.7 lw 1.3
 pause -1