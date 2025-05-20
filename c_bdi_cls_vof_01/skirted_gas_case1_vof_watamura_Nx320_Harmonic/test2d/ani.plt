reset

set size ratio -1
adj1 = 3e-2
adj2 = 2e-4

set palette defined (\
	0 'skyblue',\
	0.5 'black',\
	1 'light-pink')

unset key

set samples 1000
rr = 0.02
f1(x) = abs(x) < rr ? -(rr**2-x**2)**0.5 : 0
f2(x) = abs(x) < rr ? +(rr**2-x**2)**0.5 : 0

do for [i=0:20]{
	file=sprintf('test2d_001s_%06d.csv',i)
	
	set contour base
	set cntrparam levels discrete 0.5
	unset surface

	set table 'cont.dat'
	sp \
	file u 1:2:($12)
	unset table
	unset contour
	
	set xr [-0.021:0.021]
	set yr [-0.021:0.021]
	set cbr [0:1]

	set obj 1 circ at 0,0 size 0.02 fc rgb 'gray20' lw 2 front
	set obj 2 circ at 0.005*cos((i+1e-99)/20.0*3.14159265),0.005*sin((i+1e-99)/20.0*3.14159265) size 0.01 fc rgb 'khaki' fs trans solid 0.5 noborder front
	set obj 3 circ at 0.005*cos((i+1e-99)/20.0*3.14159265),0.005*sin((i+1e-99)/20.0*3.14159265) size 0.01 fc rgb 'gray20' lw 2 front
	
	p \
	file u 1:2:($12) w image,\
	'cont.dat' u 1:2 w l lw 4 lc rgb 'black',\
	'+' u 1:(0):(0) w l lw 2 lc rgb 'red' t '',\
	'+' u (0):1:(0) w l lw 2 lc rgb 'red' t '',\
	file u 1:2:(adj1*$4):(adj1*$5) ev 4:4 w vec lw 2 lc rgb 'dark-violet' ,\
	'+' using 1:(f1($1)):(-0.025) with filledcurves closed lc rgb 'white' fs trans solid 0.5 t '',\
	'+' using 1:(f2($1)):( 0.025) with filledcurves closed lc rgb 'white' fs trans solid 0.5 t '',\
	# f1(x) w l lw 4 lc rgb 'black' t '',\
	# f2(x) w l lw 4 lc rgb 'black' t '',\
	# # 0 lw 4 lc rgb 'red' t '',\
	# # # file u 1:2:(adj2*$19):(adj2*$20) ev 1:1 w vec lw 2 lc rgb 'dark-violet'
	pause 0.1
	
	unset obj 1
}
