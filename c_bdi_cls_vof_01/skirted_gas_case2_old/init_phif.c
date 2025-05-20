void init_phif(double **phir, double **phim, double **phif)
{
	int ii, jj;
	// int jlevl;
	// double ycent, ylevl;
	
	int ti,tj;
	int ni,nj;
	double radi;
	double dis_xm, dis_xp;
	double dis_ym, dis_yp;
	double rr_mm, rr_mp, rr_pm, rr_pp;
	double rr_ti, rr_tj, rr_ij;
	double radi_in, radi_out;
	double vol, dvol;
	
	radi = para.lx*0.5*para.Ri;
	radi_in = radi*(para.nx-8-para.bet*2)/(para.nx-8);
	radi_out = radi*(para.nx-8+para.bet*2)/(para.nx-8);
	ni = nj = 16;

	// xcent = para.lx*para.nx/(double)(para.nx-8)*0.5;
	// ycent = para.ly*para.ny/(double)(para.ny-8)*0.5;
	// ylevl = ycent + 0.0;
	// jlevl = (int)(ylevl/para.dy);
	
	

	

	// #ifdef _OPENMP
	// #pragma omp parallel for \
	// default(none) \
	// private(ii, jj) \
	// shared(para, jlevl, ylevl, phir, phim, phif)
	// #endif
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj,ti,tj) \
	private(dis_xm, dis_xp, dis_ym, dis_yp) \
	private(rr_mm, rr_mp, rr_pm, rr_pp, rr_ti, rr_tj, rr_ij, vol, dvol) \
	shared(para, phir, phim, phif, radi, radi_in, radi_out, ni, nj)
	#endif
	for(jj=0; jj<para.ny+4; jj++){
	for(ii=0; ii<para.nx+4; ii++){
		dis_xm = ((double)ii - 2.0)*para.dx - 0.5*para.lx*para.nx/(double)(para.nx-8);
		dis_xp = ((double)ii - 1.0)*para.dx - 0.5*para.lx*para.nx/(double)(para.nx-8);
		dis_ym = ((double)jj - 2.0)*para.dy - 0.25*para.ly*para.ny/(double)(para.ny-8);
		dis_yp = ((double)jj - 1.0)*para.dy - 0.25*para.ly*para.ny/(double)(para.ny-8);
		
		rr_mm = sqrt(pow(dis_ym, 2.0) + pow(dis_xm, 2.0));
		rr_mp = sqrt(pow(dis_ym, 2.0) + pow(dis_xp, 2.0));
		rr_pm = sqrt(pow(dis_yp, 2.0) + pow(dis_xm, 2.0));
		rr_pp = sqrt(pow(dis_yp, 2.0) + pow(dis_xp, 2.0));
		
			if(rr_mm < radi_in
			&& rr_mp < radi_in
			&& rr_pm < radi_in
			&& rr_pp < radi_in){	
				vol = 1.0 - phir[jj][ii] - phim[jj][ii];			
			}
			else if(   rr_mm > radi_out
					&& rr_mp > radi_out
					&& rr_pm > radi_out
					&& rr_pp > radi_out){
					vol = 0.0;
			}else {		
			vol =0.0;			
				for(tj = 0; tj < nj; tj++){
				for(ti = 0; ti < ni; ti++){
						rr_ti = para.dx/(double)ni*((double)ti+0.5) + dis_xm;
						rr_tj = para.dy/(double)nj*((double)tj+0.5) + dis_ym;
						rr_ij = sqrt(pow(rr_ti, 2.0) + pow(rr_tj, 2.0));
						dvol =-0.5*(tanh(-para.bet*(double)((radi-rr_ij)*para.idc))-1.0)
									   *(1.0 - phir[jj][ii] - phim[jj][ii])/(double)ni/(double)nj;		
						vol += 	dvol;
				}
				}
			}
			phif[jj][ii] =vol;		
	}
	}
	
	// for(jj=0; jj<para.ny+4; jj++){
		
		// if(jj - 1.5 > jlevl + para.bet){
			// for(ii=0; ii<para.nx+4; ii++){
				// phif[jj][ii] = 1.0 - phir[jj][ii] - phim[jj][ii];
			// }
		// }
		// else if(fabs(jj - 1.5 - jlevl) <= para.bet)
			// for(ii=0; ii<para.nx+4; ii++){
				// phif[jj][ii] = -0.5*(tanh(-para.bet*(double)(jj-1.5-jlevl))-1.0)
								   // *(1.0 - phir[jj][ii] - phim[jj][ii]);
			// }
		
	// }
	
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

