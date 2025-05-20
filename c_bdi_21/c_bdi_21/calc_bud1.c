# include "watasys.h"
# include "watafnc.h"

double calc_bud1(coord *uu2, coord *uu3, double **phir, double **phim)
{
	int ii, jj;	
	int ex, ey;
	double tmp=0.0;
	double ans=0.0;
	double pw;
	coord0 vel2, vel3;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, vel2, vel3, tmp, pw) \
	reduction(+:ans) \
	shared(ex, ey, para, uu2, uu3, phir, phim)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		tmp = 0.0;
		
		vel2.xx = uu2->xx[jj  ][ii-1];
		vel3.xx = uu3->xx[jj  ][ii-1];
		
		pw = (  +phir [jj  ][ii-1] + phim [jj  ][ii-1]
				+phir [jj  ][ii  ] + phim [jj  ][ii  ])*0.5;
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		
		tmp += (+0.5*vel3.xx*vel3.xx
			    -0.5*vel2.xx*vel2.xx)* (1.0 - pw);
		
		
		vel2.xx = uu2->xx[jj  ][ii  ];
		vel3.xx = uu3->xx[jj  ][ii  ];
		
		pw = (  +phir [jj  ][ii  ] + phim [jj  ][ii  ]
				+phir [jj  ][ii+1] + phim [jj  ][ii+1])*0.5;
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		
		tmp += (+0.5*vel3.xx*vel3.xx
			    -0.5*vel2.xx*vel2.xx)* (1.0 - pw);
			   
			   
		vel2.yy = uu2->yy[jj-1][ii  ];
		vel3.yy = uu3->yy[jj-1][ii  ];
		
		pw = (  +phir [jj-1][ii  ] + phim [jj-1][ii  ]
				+phir [jj  ][ii  ] + phim [jj  ][ii  ])*0.5;
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		
		tmp += (+0.5*vel3.yy*vel3.yy
			    -0.5*vel2.yy*vel2.yy)* (1.0 - pw);
		
		
		vel2.yy = uu2->yy[jj  ][ii  ];
		vel3.yy = uu3->yy[jj  ][ii  ];
		pw = (  +phir [jj  ][ii  ] + phim [jj  ][ii  ]
				+phir [jj+1][ii  ] + phim [jj+1][ii  ])*0.5;
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
			
		tmp += (+0.5*vel3.yy*vel3.yy
			    -0.5*vel2.yy*vel2.yy)* (1.0 - pw);
		
		
		tmp *= 0.5*para.ro*para.idT*(para.dx*para.dy);
		ans += tmp;
	}
	}
	
	return(ans);

}

