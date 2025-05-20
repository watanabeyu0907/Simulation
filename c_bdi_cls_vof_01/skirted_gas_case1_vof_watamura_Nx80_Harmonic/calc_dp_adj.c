void calc_dp_adj(double **dpp, double **phir, double **phim)
{
	int ii, jj;
	double sum=0.0,frc=0.0;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(para, dpp, phir, phim) \
	reduction(+:sum,frc)
	#endif
	for(jj=2; jj<para.ny+2; jj++){
	for(ii=2; ii<para.nx+2; ii++){
		
		sum+=(1.0-phir[jj][ii]-phim[jj][ii])*dpp[jj][ii];
		frc+=(1.0-phir[jj][ii]-phim[jj][ii]);
		
	}
	}
	
	sum/=frc;
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(para, dpp, sum)
	#endif
	for(jj=2; jj<para.ny+2; jj++){
	for(ii=2; ii<para.nx+2; ii++){
		
		dpp[jj][ii] -= sum;
		
	}
	}


}

