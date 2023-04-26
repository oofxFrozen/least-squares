set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'results.png'

set title "Least Squares Approximation"
set xlabel "x"
set ylabel "y"

set xrange [-20:20]
set yrange [-20:20]

plot 'results.dat' using 1:2 with points title "Points on Curve", \
     'data.dat' using 1:2 with points title "Data Points"
