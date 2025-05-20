set terminal epslatex standalone header \
" \
\\usepackage{arev}\n \
\\usepackage[]{arevmath}\n \
\\usepackage{amssymb, amsmath, bm, color, xcolor}\n \
\\usepackage{textcomp}\n \
\\usepackage{siunitx}[=v2]\n \
\\DeclareMathSymbol{I}{\\mathalpha}{extraitalic}{30}\n \
\\DeclareMathSymbol{a}{\\mathalpha}{extraitalic}{50}\n \
\\DeclareMathSymbol{f}{\\mathalpha}{extraitalic}{55}\n \
\\DeclareMathSymbol{i}{\\mathalpha}{extraitalic}{58}\n \
\\DeclareMathSymbol{l}{\\mathalpha}{extraitalic}{61}\n \
\\DeclareMathSymbol{u}{\\mathalpha}{extraitalic}{70}\n \
\\DeclareMathSymbol{v}{\\mathalpha}{extraitalic}{71}\n \
\\DeclareMathSymbol{w}{\\mathalpha}{extraitalic}{72}\n \
\\DeclareMathSymbol{x}{\\mathalpha}{extraitalic}{73}\n \
\\sisetup{math-micro={\\usefont{T1}{phv}{m}{n}\\text{\\textmu}}}\n \
"
# " \
# \\usepackage{amssymb, amsmath, bm}\n\
# \\usepackage{siunitx}\n \
# "


set terminal epslatex size 4, 2 standalone color solid 9
 
set encoding utf8
set datafile separator ','

#set view map
 
# set pm3d at b
# set pm3d map
# set pm3d interpolate 5, 5
 
set tics front 
set border front
set border lw 2
set size ratio 1

 
# unset border
# unset key
 
# nn = 900
if (!exists("file01")) file01='test'
oname = sprintf('%s.tex', file00)
file01 = sprintf('%s.csv', file01)
set output oname



# set output out01
# set multiplot layout 1, 2 title'' # row,col
set multiplot

set style line 100 lt 1 lw 4 ps 2 lc rgb '#7f7f7f' # gray
set style line 101 lt 1 lw 4 ps 2 lc rgb '#d62728' # red
set style line 102 lt 1 lw 4 ps 2 lc rgb '#1f77b4' # blue
set style line 103 lt 1 lw 4 ps 2 lc rgb '#2ca02c' # green
set style line 104 lt 1 lw 4 ps 2 lc rgb '#ff7f0e' # orange
set style line 105 lt 1 lw 4 ps 2 lc rgb '#17becf' # light-blue
set style line 106 lt 1 lw 4 ps 2 lc rgb '#9467bd' # purple
set style line 107 lt 1 lw 4 ps 2 lc rgb '#bcbd22' # yellow
set style line 108 lt 1 lw 4 ps 2 lc rgb '#e377c2' # orchid
set style line 109 lt 1 lw 4 ps 2 lc rgb '#8c564b' # brown
 
set style line 110 lt 1 lw 4 ps 2 lc rgb '#c7c7c7' # gray
set style line 111 lt 1 lw 4 ps 2 lc rgb '#ff9896' # red
set style line 112 lt 1 lw 4 ps 2 lc rgb '#aec7e8' # blue
set style line 113 lt 1 lw 4 ps 2 lc rgb '#98df8a' # green
set style line 114 lt 1 lw 4 ps 2 lc rgb '#ffbb78' # orange
set style line 115 lt 1 lw 4 ps 2 lc rgb '#9edae5' # light-blue
set style line 116 lt 1 lw 4 ps 2 lc rgb '#c5b0d5' # purple
set style line 117 lt 1 lw 4 ps 2 lc rgb '#dbdb8d' # yellow
set style line 118 lt 1 lw 4 ps 2 lc rgb '#f7b6d2' # orchid
set style line 119 lt 1 lw 4 ps 2 lc rgb '#c49c94' # brown

set style line 120 lt 1 lw 4 ps 2 lc rgb '#4d4d4d' #gray
set style line 121 lt 1 lw 4 ps 2 lc rgb '#991c1c' #red
set style line 122 lt 1 lw 4 ps 2 lc rgb '#165380' #blue
set style line 123 lt 1 lw 4 ps 2 lc rgb '#154d15' #green
set style line 124 lt 1 lw 4 ps 2 lc rgb '#cc650a' #orange
set style line 125 lt 1 lw 4 ps 2 lc rgb '#118b99' #light-blue
set style line 126 lt 1 lw 4 ps 2 lc rgb '#4f3766' #purple
set style line 127 lt 1 lw 4 ps 2 lc rgb '#808017' #yellow
set style line 128 lt 1 lw 4 ps 2 lc rgb '#995083' #orchid
set style line 129 lt 1 lw 4 ps 2 lc rgb '#4d2f29' #brown


set format '$%g$'

tht = 15.0*pi/180.0
cs  = cos(tht)
sn  = sin(tht)
 
aj = 1000

####################################################################################################
####################################################################################################
set size ratio 0.5



# set label 99 left at screen 0.005,0.975 '\LARGE (a)'

# set logscale y

set format x '$%g$'
set xlabel '{\Large $t$ [\SI{}{s}]}' offset 0,0.3
set xrange [0:1]
set xtics 0.5
set mxtics 5

set format y '$%g$'
set ylabel '{\Large Ene. bud., $\dot{E}$ [\SI{}{kgm/s^3}]}' offset 0.0,0.0
# set yrange [-0.03:0.06]
# set ytics 0.03 offset 0,0
# set mytics 3


set key width -7 left top box opaque spacing 1.2 samplen 1 Left reverse
# set key at graph 0.015,0.9775 
# unset key
set key maxrow 3

####################################################################################################
####################################################################################################

set lmargin screen 0.19
set rmargin screen 0.99
set bmargin screen 0.18
set tmargin screen 0.95

####################################################################################################
####################################################################################################

set palette defined (0 '#ffffff',\
                     1 '#00008b',\
					 2 '#2ca9e1',\
					 3 '#008000',\
					 4 '#c8c800',\
					 5 '#ff0000',\
					 6 '#1a1a1a') positive

set samples 400
e(x) = 1.0/(0.5-0.57722-log(x/8.0))

####################################################################################################

ri  = 2.5e-3
ell = 25e-3
rho = 998
omg = 3.14

ao = 1.0
ai = 0.5
		
####################################################################################################
####################################################################################################

plot \
0 ls 100 lw 2 dt(2,2) t '',\
file01 u 1:4  w l ls 101 title '\small $\dot{E}_{-\text{d}K/\text{d}t}$',\
file01 u 1:5  w l ls 102 title '\small $\dot{E}_\text{diss}$',\
file01 u 1:($7+$8)  w l ls 103 title '\small $\dot{E}_\text{ext}$',\
file01 u 1:9  w l ls 104 title '\small $\dot{E}_\text{surf}$',\
file01 u 1:10 w l ls 105 title '\small $\dot{E}_\text{grav}$',\
file01 u 1:11 w l ls 100 lc rgb 'black' title '\small $\dot{E}_\text{all}$',\

####################################################################################################
####################################################################################################

unset multiplot
reset
