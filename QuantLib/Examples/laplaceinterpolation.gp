cd '/home/peter/quantlib/QuantLib/Examples'
! N=25 DELPERC=0.80 ./laplaceinterpolation >laplaceinterpolation.dat
set hidden3d
unset key
splot 'laplaceinterpolation.dat' u 1:2:3:3 w p pt 7 ps 0.5 palette # original
splot 'laplaceinterpolation.dat' u 1:2:($3*1/(1.0-$5)):3 w p pt 7 ps 0.5 palette # destructed
splot 'laplaceinterpolation.dat' u 1:2:4:4 w p pt 7 ps 0.5 palette # reconstructued
splot 'laplaceinterpolation.dat' u 1:2:($4-$3) w pm3d # error

splot 'laplaceinterpolation.dat' u 1:2:3 w pm3d # original
splot 'laplaceinterpolation.dat' u 1:2:4 w pm3d # reconstructued
