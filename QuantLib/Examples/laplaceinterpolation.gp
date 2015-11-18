set hidden3d
splot 'laplaceinterpolation.dat' u 1:2:3 w l
splot 'laplaceinterpolation.dat' u 1:2:4 w l
splot 'laplaceinterpolation.dat' u 1:2:($3*1/(1.0-$5)) w p
