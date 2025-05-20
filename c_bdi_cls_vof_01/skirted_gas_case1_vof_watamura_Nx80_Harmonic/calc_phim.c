void calc_phim(double *hev, double **phir, double **phim, double **phif, double **phit, double **dlt, double **lvs, int flg)
{
	int ii, jj;
	int ex, ey;
	double tmp;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double sum3 = 0.0;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex, ey, phif) \
	reduction(+:sum1)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		sum1 += phif[jj][ii];
		
	}
	}
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, tmp) \
	shared(ex, ey, hev, phir, phim, phif, phit, dlt, lvs, flg) \
	reduction(+:sum2, sum3)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		// @ cell centre
		tmp = phim[jj][ii];
		dlt[jj][ii] = get_smd(lvs[jj][ii]);
		phim[jj][ii] = get_hev(lvs[jj][ii], hev);
		phim[jj][ii] = MIN(MAX(phim[jj][ii],1.0e-9),1.0-1.0e-9);
		if(phit[jj][ii] > 0.5 && flg != 0){
			phif[jj][ii] += tmp-phim[jj][ii];
		}
		sum2 += phif[jj][ii];
		sum3 += pow(phif[jj][ii],2.0)*pow(1.0-phif[jj][ii]-phir[jj][ii]-phim[jj][ii],2.0);
		
	}
	}
	
	sum1 = (sum1-sum2)/(sum3+1e-19);
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex, ey, phir, phim, phif, sum1)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		phif[jj][ii]+= sum1*pow(phif[jj][ii],2.0)*pow(1.0-phif[jj][ii]-phir[jj][ii]-phim[jj][ii],2.0);
		phif[jj][ii] = MAX(MIN(phif[jj][ii], 1.0), 0.0);
		
	}
	}
	
	boundary_pp(phif);


}

