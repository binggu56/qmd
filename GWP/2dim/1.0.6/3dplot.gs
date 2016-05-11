#set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 500, 350 
# set output 'surface1.13.png'
set pm3d map  
#set dgrid3d 30,30 
#set samples 51, 51
set isosamples 100, 100
#set title "3D gnuplot demo" 
#set xlabel "X axis" 
#set xlabel  offset character -3, -2, 0 font "" textcolor lt -1 norotate
#set xrange [ -1.00000 : 1.00000 ] noreverse nowriteback
#set ylabel "Y axis" 
#set ylabel  offset character 3, -2, 0 font "" textcolor lt -1 rotate by -270
#set yrange [ -1.00000 : 1.00000 ] noreverse nowriteback
#set zlabel "Z axis" 
#set zlabel  offset character -5, 0, 0 font "" textcolor lt -1 norotate
#set zrange [ 0.00000 : 1.00000 ] noreverse nowriteback
#set hidden3d
splot 'wf0.dat' u 1:2:3
pause -1 
