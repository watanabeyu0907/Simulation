# include "watasys.h"
# include "watafnc.h"

void calc_phim(double *hev, double **phir, double **phim, double **dlt, double **lvs, int flg)
{
	int ii, jj;
	int ex, ey;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex, ey, hev, phir, phim, dlt, lvs, flg) 
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		// @ cell centre
		dlt [jj][ii] = get_smd(lvs[jj][ii]);
		phim[jj][ii] = get_hev(lvs[jj][ii], hev);
		phim[jj][ii] = MIN(MAX(phim[jj][ii],1.0e-9),1.0-1.0e-9);
	}
	}


}

