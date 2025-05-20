# include "watasys.h"
# include "watafnc.h"

void calc_NSeq2_yy(coord *uu2, coord *uu3, coord *adv, double **phir, double **phim, coord *vvv, int rem, double *sum)
{
	int ii, jj;
	int ex, ey;
	int sx[8] = {0, 1, 1, 0};
	int sy[8] = {0, 1, 0, 1};

	double vis2, vis3;
	double nu, rst, pw;
	
	double epsi;
	double tmp = 0.0;
	
	
	ex = para.nx+2;
	ey = para.ny+1;
	
	
	// ( 1-b ) y vel
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
		
		epsi = uu3->yy[jj][ii];
		if(epsi==0.0) epsi=1.0;
		
		nu = para.mu/para.ro;
		
		vis2 = +1.0*(para.idx*para.idx)*uu2->yy[jj  ][ii-1]
			   -2.0*(para.idx*para.idx)*uu2->yy[jj  ][ii  ]
			   +1.0*(para.idx*para.idx)*uu2->yy[jj  ][ii+1]
			   
			   +1.0*(para.idy*para.idy)*uu2->yy[jj-1][ii  ]*2.0
			   -2.0*(para.idy*para.idy)*uu2->yy[jj  ][ii  ]*2.0
			   +1.0*(para.idy*para.idy)*uu2->yy[jj+1][ii  ]*2.0
			   
			   +1.0*(para.idx*para.idy)*uu2->xx[jj  ][ii-1]
			   -1.0*(para.idx*para.idy)*uu2->xx[jj+1][ii-1]
			   -1.0*(para.idx*para.idy)*uu2->xx[jj  ][ii  ]
			   +1.0*(para.idx*para.idy)*uu2->xx[jj+1][ii  ]
			   ;
		
		vis3 = +1.0*(para.idx*para.idx)*uu3->yy[jj  ][ii-1]
			   +1.0*(para.idx*para.idx)*uu3->yy[jj  ][ii+1]
			   
			   +1.0*(para.idy*para.idy)*uu3->yy[jj-1][ii  ]*2.0
			   +1.0*(para.idy*para.idy)*uu3->yy[jj+1][ii  ]*2.0
			   
			   +1.0*(para.idx*para.idy)*uu3->xx[jj  ][ii-1]
			   -1.0*(para.idx*para.idy)*uu3->xx[jj+1][ii-1]
			   -1.0*(para.idx*para.idy)*uu3->xx[jj  ][ii  ]
			   +1.0*(para.idx*para.idy)*uu3->xx[jj+1][ii  ]
			   ;
				
		rst =  -2.0*(para.idx*para.idx)
			   -2.0*(para.idy*para.idy)*2.0;
			   
		pw = +0.5*(phir[jj][ii]+phir[jj+1][ii])
			 +0.5*(phim[jj][ii]+phim[jj+1][ii]);
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		
		
		// y vel
		uu3->yy[jj][ii]
			=+0.5*(phir[jj][ii]+phir[jj+1][ii]) * 0.0
			 +0.5*(phim[jj][ii]+phim[jj+1][ii]) * vvv->yy[jj][ii]
			+(1.0-pw)*(uu2->yy[jj][ii]
			 +para.dT*(
				+adv->yy[jj][ii]
				+0.5*nu*vis2
				+0.5*nu*vis3
			 )
			)
			/(
				+1.0
				-0.5*nu*rst*para.dT
			);			
			
		epsi = 1.0 - uu3->yy[jj][ii]/epsi;
		
		tmp+=epsi*epsi;
	}
	}
	
	*sum = tmp;

}

