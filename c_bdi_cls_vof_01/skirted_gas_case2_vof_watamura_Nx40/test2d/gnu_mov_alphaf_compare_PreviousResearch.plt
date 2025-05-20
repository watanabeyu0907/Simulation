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


set terminal epslatex size 4.2, 3.0 standalone color solid 9
 
set encoding utf8
# set datafile separator ','

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
if (!exists("nn")) nn=0
oname = sprintf('%s%06d.tex', file01, nn)
ifile = sprintf('%s%06d.csv', file01, nn)
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



####################################################################################################
####################################################################################################

set contour base
set cntrparam levels discrete 0.5
unset surface

set table 'cont.dat'
sp \
ifile u 1:2:($12)
unset table
unset contour
	
####################################################################################################
####################################################################################################
	
set format '$%g$'
set samples 1000

tht = 15.0*pi/180.0
cs  = cos(tht)
sn  = sin(tht)

aa = 3.0e-1
bb = aa*0.33

fps= 10.0

aj = 1.0

ro = 20
ri = 10
ar =  5

t = nn/fps
f1(x) = abs(x) < ro ? -(ro**2-x**2)**0.5 : 0
f2(x) = abs(x) < ro ? +(ro**2-x**2)**0.5 : 0

####################################################################################################
####################################################################################################
	
set size ratio -1

set lmargin screen 0.17
set rmargin screen 0.88
set bmargin screen 0.17
set tmargin screen 0.91

####################################################################################################
####################################################################################################

set xlabel '{\Large $x$ [\SI{}{-}]}' offset 0,0.5
set xrange [0:1]
set xtics 0.5 offset 0,0.1
set mxtics 2

set ylabel '{\Large $y$ [\SI{}{-}]}' offset 0.8,0.0
set yrange [0.5:1.5]
set ytics 0.5 offset 0.3,0
set mytics 2

set cblabel '{\Large $\phi_\text{gas}$ [\SI{}{-}]}' offset 0.2,0.0
set cbrange [0:1]
set cbtics 0.5 offset -0.3,0
set mcbtics 5
set palette defined (\
	0 'skyblue',\
	0.5 'white',\
	1 'light-pink')
	
set colorbox vert user origin graph 1.1,0 size graph 0.05,1

set key
set key box opaque spacing 1.1 samplen 2 width 2 Left reverse 
set key at screen 0.44, 0.90    # → 画面右下に近い位置


####################################################################################################
####################################################################################################

aaa = sprintf('{\Large $t =  %0.2f $ [\SI{}{-}]}',t)
set label 1 at graph 0.5,1.06 aaa front

mouth = 180*nn/fps -360*floor((180*nn/fps)/360)
# set obj 91 circle at ar*cos(3.14*nn/fps),ar*sin(3.14*nn/fps) size ri fs solid fc rgb 'khaki' lw 2 front
# set obj 92 circle at ar*cos(3.14*nn/fps),ar*sin(3.14*nn/fps) size ri fc rgb 'black' lw 2 front
# set obj 93 circle at (ar+ar*0.8)*cos(3.14*nn/fps)-ar*0.8*sin(3.14*nn/fps),(ar+ar*0.8)*sin(3.14*nn/fps)+ar*0.8*cos(3.14*nn/fps) size 0.5 fs solid noborder fc rgb 'black' front
# set obj 94 circle at (ar-ar*0.8)*cos(3.14*nn/fps)-ar*0.8*sin(3.14*nn/fps),(ar-ar*0.8)*sin(3.14*nn/fps)+ar*0.8*cos(3.14*nn/fps) size 0.5 fs solid noborder fc rgb 'black' front
# set obj 95 circle at ar*cos(3.14*nn/fps),ar*sin(3.14*nn/fps)size ri*0.6 arc [180+mouth:mouth] nowedge fc rgb 'black' lw 2 front
# set obj 96 circle at 0,0 size ro fc rgb 'black' lw 2 front
plot \
ifile u ($1*aj+0.5):($2*aj+1):12 w image t '',\
'cont.dat' u ($1*aj+0.5):($2*aj+1) w l lw 5 lc rgb '#1a5c8a' t 'Present',\
"case2_MooNMD.csv" using 1:2 with line lw 5 lc rgb '#ff7f0e' pt 7  title 'MooNMD',\
# "case2_FreeLIFE_left.csv" using 1:2 with line lw 5 lc rgb '#ff7f0e' pt 7  notitle ,\
# "case2_FreeLIFE_right.csv" using 1:2 with line lw 5 lc rgb '#ff7f0e' pt 7  notitle 

# ifile using ($1*aj+0.5):($2*aj+1):($4*aa*aj):($5*aa*aj) ev 8:8:1:1 with vector lw 2 lc rgb 'forest-green' nohead title "",\
# ifile using (($1+$4*aa)*aj+0.5):(($2+$5*aa)*aj+1):((-cs*$4+sn*$5)*aj*bb):(-( sn*$4+cs*$5)*aj*bb) ev 8:8:1:1 with vector lw 2 lc rgb 'forest-green' nohead title "",\
# ifile using (($1+$4*aa)*aj+0.5):(($2+$5*aa)*aj+1):((-cs*$4-sn*$5)*aj*bb):(-(-sn*$4+cs*$5)*aj*bb) ev 8:8:1:1 with vector lw 2 lc rgb 'forest-green' nohead title "",\



# plot \
# ifile u ($1*aj):($2*aj):12 w image t '',\
# ifile using ($1*aj):($2*aj):($4*aa*aj):($5*aa*aj) ev 8:8:1:1 with vector lw 2 lc rgb 'forest-green' nohead title "",\
# ifile using (($1+$4*aa)*aj):(($2+$5*aa)*aj):((-cs*$4+sn*$5)*aj*bb):(-( sn*$4+cs*$5)*aj*bb) ev 8:8:1:1 with vector lw 2 lc rgb 'forest-green' nohead title "",\
# ifile using (($1+$4*aa)*aj):(($2+$5*aa)*aj):((-cs*$4-sn*$5)*aj*bb):(-(-sn*$4+cs*$5)*aj*bb) ev 8:8:1:1 with vector lw 2 lc rgb 'forest-green' nohead title "",\
# 'cont.dat' u ($1*aj):($2*aj) w l lw 1 lc rgb 'white' t '',\
# '+' using 1:(f1($1)):(-25) with filledcurves closed lc rgb 'white' fs trans solid 0.5 t '',\
# '+' using 1:(f2($1)):( 25) with filledcurves closed lc rgb 'white' fs trans solid 0.5 t '',\

####################################################################################################
####################################################################################################

unset multiplot
reset
