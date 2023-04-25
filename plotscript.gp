# Set the output file format and name
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'results.png'

# Set the title and labels for the plot
set title "Least Squares Approximation"
set xlabel "x"
set ylabel "y"

# Set the range for the x and y axes
set xrange [-10:10]
set yrange [-10:10]

# Plot the data points and the points on the curve
plot 'results.dat' using 1:2 with points title "Points on Curve", \
     'data.dat' using 1:2 with points title "Data Points"
