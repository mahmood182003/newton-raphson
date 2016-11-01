reset
set title "results with ν_1=6.6 and ν_2=6.2" font "Times Italic, 12"
set grid lw 1
set grid ztics
set xlabel "c_{eq1}" font "Times Italic, 12"
set ylabel "c_{eq2}" font "Times Italic, 12"
set zlabel "q" offset -1,1 font "Times Italic, 12" rotate by 90
set ytics offset -1,-1
set xyplane at 0
set view 54,288
set xrange [0:0.001] 
set yrange [0:0.0001]
set zrange [0:0.006]
#set dgrid3d 40,40 qnorm 2
set hidden3d
splot "q1.txt" with points, "q2.txt" with points linecolor rgb "green"
pause -1
