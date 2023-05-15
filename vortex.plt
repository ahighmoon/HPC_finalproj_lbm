set title 'Fluid Vorticity'
set xlabel 'cell # along x-dimension'
set ylabel 'cell # along y-dimension'
set size ratio -1
set autoscale fix

set terminal png
set output './png/iter.png'

zstep=0.01
set cbrange [-0.04:0.04 ]
set cbtics zstep
set palette defined (-0.03 "red",\
                     -0.002 "orange",\
                     0 "yellow",\
                     0.002 "light-green",\
                     0.03 "green")        


#set terminal postscript eps size 12,10 enhanced color font 'Verdana,30' linewidth 2
#set output 'final_state.eps'

plot './png/iter.dat' using 1:2:4 with image
