void calc_NSeq2_xx(visc_fnc *dd, coord *uu2, coord *uu3, coord *adv, double **phir, double **phim, coord *vvv, int rem, double *sum)
{
	int ii, jj;
	int ex, ey;
	int sx[8] = {0, 1, 1, 0};
	int sy[8] = {0, 1, 0, 1};

	double vis2, vis3;
	double rst, ph;
	
	double epsi;
	double tmp = 0.0;
	
	
	ex = para.nx+1;
	ey = para.ny+2;
	
	
	// ( 1-a ) x vel
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, epsi, vis2, vis3, rst, ph) \
	shared(para, sx, sy, ex, ey, rem) \
	shared(dd, uu2, uu3, adv, phir, phim, vvv) \
	reduction(+: tmp)
	#endif
	for(jj=2+sx[rem]; jj<ey; jj+=2){
	for(ii=2+sy[rem]; ii<ex; ii+=2){
		
		epsi = uu3->xx[jj][ii];
		if(epsi==0.0) epsi=1.0;
		
		vis2 = +dd->xmxm[jj][ii]*uu2->xx[jj  ][ii-1]
			   +dd->ccxx[jj][ii]*uu2->xx[jj  ][ii  ]
			   +dd->xpxp[jj][ii]*uu2->xx[jj  ][ii+1]
			   
			   +dd->ymym[jj][ii]*uu2->xx[jj-1][ii  ]
			   +dd->ccyy[jj][ii]*uu2->xx[jj  ][ii  ]
			   +dd->ypyp[jj][ii]*uu2->xx[jj+1][ii  ]
			   
			   +dd->ymxm[jj][ii]*uu2->yy[jj-1][ii  ]
			   +dd->ymxp[jj][ii]*uu2->yy[jj-1][ii+1]
			   +dd->ypxm[jj][ii]*uu2->yy[jj  ][ii  ]
			   +dd->ypxp[jj][ii]*uu2->yy[jj  ][ii+1]
			   ;
		
		vis3 = +dd->xmxm[jj][ii]*uu3->xx[jj  ][ii-1]
			   +dd->xpxp[jj][ii]*uu3->xx[jj  ][ii+1]
			   
			   +dd->ymym[jj][ii]*uu3->xx[jj-1][ii  ]
			   +dd->ypyp[jj][ii]*uu3->xx[jj+1][ii  ]
			   
			   +dd->ymxm[jj][ii]*uu3->yy[jj-1][ii  ]
			   +dd->ymxp[jj][ii]*uu3->yy[jj-1][ii+1]
			   +dd->ypxm[jj][ii]*uu3->yy[jj  ][ii  ]
			   +dd->ypxp[jj][ii]*uu3->yy[jj  ][ii+1]
			   ;
				
		rst =  +dd->ccxx[jj][ii]
			   +dd->ccyy[jj][ii];
			   
		ph = +0.5*(phir[jj][ii]+phir[jj][ii+1])
			 +0.5*(phim[jj][ii]+phim[jj][ii+1]);
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
			   
		
		// // x vel
		uu3->xx[jj][ii]
			=+0.5*(phir[jj][ii]+phir[jj][ii+1]) * 0.0
			 +0.5*(phim[jj][ii]+phim[jj][ii+1]) * vvv->xx[jj][ii]
			+(1.0-ph)*(uu2->xx[jj][ii]
			 +para.dT*(
				+adv->xx[jj][ii]
				+0.5*vis2
				+0.5*vis3
			 )
			)
			/(
				+1.0
				-0.5*rst*para.dT
			);
			
		epsi = 1.0 - uu3->xx[jj][ii]/epsi;
		
		tmp+=epsi*epsi;
	}
	}
	
	*sum = tmp;
	
}

