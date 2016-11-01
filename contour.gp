reset
set grid
set title "Contours of z=0 and Newton-Raphson step points"
set xrange [0:0.0005]
set yrange [0.006:0.009]
set xlabel "X axis"
set ylabel "Y axis"
set isosamples 70,70
set contour both
set cntrparam levels discrete 0
set view map
#unset surface
#set hidden3d
splot (0.028*((0.5-44.6*x-64.2*y)/0.15)**6-x) nosurface, (0.54*((0.5-44.6*x-64.2*y)/0.15)**6-y) nosurface, "points.txt" with points nocontour
pause -1
