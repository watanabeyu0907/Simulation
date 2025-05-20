reset

set size ratio -1
adj1 = 3e-2
adj2 = 5e-4

# set palette defined (\
	# 0 'khaki',\
	# 0.2 'gold',\
	# 0.5 'black',\
	# 0.8 'skyblue',\
	# 1 'light-cyan')
	
# set cbr [-4:4]

do for [i=0:20]{

	set obj 1 circ at 0,0 size 0.02 fc rgb 'red' lw 4 front
	set obj 2 circ at 0.005*cos((i+1e-99)/20.0*3.14159265*1e-0),0.005*sin((i+1e-99)/20.0*3.14159265*1e-0) size 0.01 fc rgb 'red' lw 4 front
	
	file=sprintf('test2d_001s_%06d.csv',i)
	p \
	file u 1:2:($10+$11) w image,\
	0 lw 4 lc rgb 'red' t '',\
	file u 1:2:(adj1*$4):(adj1*$5) ev 2:2 w vec lw 2 lc rgb 'forest-green' ,\
	# file u 1:2:($22) w image,\
	# file u 1:2:(1/(2.*cosh($22)**2)) w image,\
	# file u 1:2:(adj2*$19):(adj2*$20) ev 1:1 w vec lw 2 lc rgb 'dark-violet'
	# file u 1:2:(1-$10-$11-$12) w image,\
	# file u 1:2:((1-$10-$11-$12)*($12)) w image,\
	# file u 1:2:($14+((($1-0.005)**2+$2**2)**0.5-10e-3)/(40.0e-3/(64-8.0))) w image,\
	pause 0.1
	
	unset obj 1
}
	# file u 1:2:($12*(1-$10-$11-$12)+$10*(1-$10)+$11*(1-$11)) w image,\