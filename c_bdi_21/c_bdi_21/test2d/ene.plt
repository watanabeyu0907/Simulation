# set size ratio -1

set multiplot layout 1,3

file=sprintf('test2d_001s_stat.csv')

set logscale y
p \
file u 1:2 w l lw 2,\
file u 1:3 w l lw 2,\

unset logscale y
p \
file u 1:15 w l lw 2,\


p \
file u 1:4 w l lw 2,\
file u 1:5 w l lw 2,\
file u 1:7 w l lw 2,\
file u 1:($4+$5+$7) w l lw 2,\

