set title 'Fluid Vorticity'
set xlabel 'cell # along x-dimension'
set ylabel 'cell # along y-dimension'
set size ratio -1
set autoscale fix

set terminal png
set output './png/iter.png'

zstep=0.01
set cbrange [-0.015:0.015 ]
#set cbrange [-0.03:0.03]
set cbtics zstep
set palette defined (-0.04 "red",\
                     -0.03 "orange",\
                     -0.02 "goldenrod",\
                     -0.006 "yellow",\
                     -0.004 "khaki1",\
                     -0.002 "lemonchiffon",\
                     0 "gray100",\
                     0.002 "light-cyan",\
                     0.004 "light-turquoise",\
                     0.008 "cyan",\
                     0.01 "dark-cyan",\
                     0.02 "web-blue",\
                     0.03 "royalblue",\
                     0.04 "blue")  

plot './png/iter.dat' using 1:2:4 with image
