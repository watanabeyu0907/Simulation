# include "watasys.h"
# include "watafnc.h"

double calc_div_rms(double **div)
{
	int ii, jj;
	double sum=0.0;
	double ans;
	int ex, ey;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(para, ii, jj) \
	shared(ex, ey, div) \
	reduction(+:sum)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		sum+=div[jj][ii]*div[jj][ii];
		
	}
	}
	
	ans = sqrt(sum/(double)(para.nx*para.ny));
	
	return(ans);
}

