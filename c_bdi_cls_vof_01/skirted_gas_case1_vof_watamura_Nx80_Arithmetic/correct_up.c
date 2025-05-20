void correct_up(coord *uu, double **pp, double **dp, double **phir, double **phim, double **phit2, double **phit3)
{
	int ii, jj;
	double ph, ro;
		
	
	// collect uu
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, ph, ro) \
	shared(para, uu, dp, phir, phim, phit2, phit3)
	#endif
	for(jj=2; jj<para.ny+2; jj++){
	for(ii=2; ii<para.nx+1; ii++){
		
		ph = +(phir[jj][ii] + phir[jj][ii+1])*0.5
			 +(phim[jj][ii] + phim[jj][ii+1])*0.5;
 		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
		
		ro = (+phit2[jj][ii] + phit2[jj][ii+1]
			  +phit3[jj][ii] + phit3[jj][ii+1])*0.25;
		ro = 1.0 + ro*(para.ro2/para.ro1-1.0);

		uu->xx[jj][ii] -= (-dp[jj  ][ii  ]
						   +dp[jj  ][ii+1])
						*(1.0-ph)/ro*para.idx*para.dT;
	}
	}
	
	// collect vv
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, ph, ro) \
	shared(para, uu, dp, phir, phim, phit2, phit3)
	#endif
	for(jj=2; jj<para.ny+1; jj++){
	for(ii=2; ii<para.nx+2; ii++){
		
		ph = +(phir[jj][ii] + phir[jj+1][ii])*0.5
			 +(phim[jj][ii] + phim[jj+1][ii])*0.5;
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
		
		ro = (+phit2[jj][ii] + phit2[jj+1][ii]
			  +phit3[jj][ii] + phit3[jj+1][ii])*0.25;
		ro = 1.0 + ro*(para.ro2/para.ro1-1.0);
			 
		uu->yy[jj][ii] -= (-dp[jj  ][ii  ]
						   +dp[jj+1][ii  ])
						*(1.0-ph)/ro*para.idy*para.dT;
	}
	}

	// correct pressure
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(para, dp, pp, phir)
	#endif
	for(jj=1; jj<para.ny+3; jj++){
	for(ii=1; ii<para.nx+3; ii++){
		
		pp[jj][ii] += dp[jj][ii];
		
	}
	}
	
}

