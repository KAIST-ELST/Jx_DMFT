set xra [  0.0000:   2.9392]
set yra [  0.0000: 321.7375]
set ylabel "Energy (meV)" 
set xtics ("G"    0.0000,"N"    0.4972,"P"    0.9326,"G"    1.3651,"H"    2.2301,"N"    2.9392)
plot "magnon.dat" u 1:2 w l 
replot "vline.dat" u 1:2 w l lc 0 
pause -1
