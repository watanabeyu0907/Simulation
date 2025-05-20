# include "watasys.h"
# include "watafnc.h"

void calc_NSeq2_xx(coord *uu2, coord *uu3, coord *adv, double **phir, double **phim, coord *vvv, int rem, double *sum)
{
	int ii, jj;
	int ex, ey;
	int sx[8] = {0, 1, 1, 0};
	int sy[8] = {0, 1, 0, 1};

	double vis2, vis3;
	double nu, rst, pw;
	
	double epsi;
	double tmp = 0.0;
	
	
	ex = para.nx+1;
	ey = para.ny+2;
	
	
	// ( 1-a ) x vel
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, epsi, nu, vis2, vis3, rst, pw) \
	shared(para, sx, sy, ex, ey, rem) \
	shared(uu2, uu3, adv, phir, phim, vvv) \
	reduction(+: tmp)
	#endif
	for(jj=2+sx[rem]; jj<ey; jj+=2){
	for(ii=2+sy[rem]; ii<ex; ii+=2){
		
		epsi = uu3->xx[jj][ii];
		if(epsi==0.0) epsi=1.0;
		
		nu = para.mu/para.ro;
		
		vis2 = +1.0*(para.idx*para.idx)*uu2->xx[jj  ][ii-1]*2.0
			   -2.0*(para.idx*para.idx)*uu2->xx[jj  ][ii  ]*2.0
			   +1.0*(para.idx*para.idx)*uu2->xx[jj  ][ii+1]*2.0
			   
			   +1.0*(para.idy*para.idy)*uu2->xx[jj-1][ii  ]
			   -2.0*(para.idy*para.idy)*uu2->xx[jj  ][ii  ]
			   +1.0*(para.idy*para.idy)*uu2->xx[jj+1][ii  ]
			   
			   +1.0*(para.idx*para.idy)*uu2->yy[jj-1][ii  ]
			   -1.0*(para.idx*para.idy)*uu2->yy[jj-1][ii+1]
			   -1.0*(para.idx*para.idy)*uu2->yy[jj  ][ii  ]
			   +1.0*(para.idx*para.idy)*uu2->yy[jj  ][ii+1]
			   ;
		
		vis3 = +1.0*(para.idx*para.idx)*uu3->xx[jj  ][ii-1]*2.0
			   +1.0*(para.idx*para.idx)*uu3->xx[jj  ][ii+1]*2.0
			   
			   +1.0*(para.idy*para.idy)*uu3->xx[jj-1][ii  ]
			   +1.0*(para.idy*para.idy)*uu3->xx[jj+1][ii  ]
			   
			   +1.0*(para.idx*para.idy)*uu3->yy[jj-1][ii  ]
			   -1.0*(para.idx*para.idy)*uu3->yy[jj-1][ii+1]
			   -1.0*(para.idx*para.idy)*uu3->yy[jj  ][ii  ]
			   +1.0*(para.idx*para.idy)*uu3->yy[jj  ][ii+1]
			   ;
				
		rst =  -2.0*(para.idx*para.idx)*2.0
			   -2.0*(para.idy*para.idy);
			   
		pw = +0.5*(phir[jj][ii]+phir[jj][ii+1])
			 +0.5*(phim[jj][ii]+phim[jj][ii+1]);
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
			   
		
		// x vel
		uu3->xx[jj][ii]
			=+0.5*(phir[jj][ii]+phir[jj][ii+1]) * 0.0
			 +0.5*(phim[jj][ii]+phim[jj][ii+1]) * vvv->xx[jj][ii]
			+(1.0-pw)*(uu2->xx[jj][ii]
			 +para.dT*(
				+adv->xx[jj][ii]
				+0.5*nu*vis2
				+0.5*nu*vis3
			 )
			)
			/(
				+1.0
				-0.5*nu*rst*para.dT
			);
			
		epsi = 1.0 - uu3->xx[jj][ii]/epsi;
		
		tmp+=epsi*epsi;
	}
	}
	
	*sum = tmp;
	
}

