set title 'Fluid Velocity'
set xlabel 'cell # along x-dimension'
set ylabel 'cell # along y-dimension'
set size ratio -1
set autoscale fix

set terminal png
set output './png/iter.png'

#set terminal postscript eps size 12,10 enhanced color font 'Verdana,30' linewidth 2
#set output 'final_state.eps'

plot './png/iter.dat' using 1:2:3 with image
