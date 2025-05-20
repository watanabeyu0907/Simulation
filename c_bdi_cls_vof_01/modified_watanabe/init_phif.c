void init_phif(double **phir, double **phim, double **phif)
{
	int ii, jj;
	int jlevl;
	double ycent, ylevl;
	
	
	ycent = para.ly*para.ny/(double)(para.ny-8)*0.5;
	ylevl = ycent -0.01;
	jlevl = (int)(ylevl/para.dy);
	

	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(para, jlevl, ylevl, phir, phim, phif)
	#endif
	for(jj=0; jj<para.ny+4; jj++){
		
		if(jj - 1.5 > jlevl + para.bet){
			for(ii=0; ii<para.nx+4; ii++){
				phif[jj][ii] = 1.0 - phir[jj][ii] - phim[jj][ii];
			}
		}
		else if(fabs(jj - 1.5 - jlevl) <= para.bet)
			for(ii=0; ii<para.nx+4; ii++){
				phif[jj][ii] = -0.5*(tanh(-para.bet*(double)(jj-1.5-jlevl))-1.0)
								   *(1.0 - phir[jj][ii] - phim[jj][ii]);
			}
		
	}
	
	// int ii, jj;
	// int ilevl;
	// double xcent, xlevl;
	
	
	// xcent = para.lx*para.nx/(double)(para.nx-8)*0.5;
	// xlevl = xcent + 0.0;
	// ilevl = (int)(xlevl/para.dx);
	

	// #ifdef _OPENMP
	// #pragma omp parallel for \
	// default(none) \
	// private(ii, jj) \
	// shared(para, ilevl, xlevl, phir, phim, phif)
	// #endif
	// for(ii=0; ii<para.nx+4; ii++){
		// if(ii >= ilevl + 3){
			// for(jj=0; jj<para.ny+4; jj++){
				// phif[jj][ii] = 1.0 - phir[jj][ii] - phim[jj][ii];
			// }
		// }
		// else if(ii == ilevl + 2){
			// for(jj=0; jj<para.ny+4; jj++){
				// phif[jj][ii] = (1.0 - phir[jj][ii] - phim[jj][ii])
							  // *((double)(ilevl+1)-xlevl*para.idy);
		// }
			// }
	// }
	
}

