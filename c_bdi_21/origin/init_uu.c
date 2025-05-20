# include "watasys.h"
# include "watafnc.h"

void init_uu(coord *uu)
{
	int ii, jj;
	int ex, ey;
	
	
	ex = para.nx+4;
	ey = para.ny+4;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex, ey, uu)
	#endif
	for(jj=0; jj<ey; jj++){
	for(ii=0; ii<ex; ii++){
		uu->xx[jj][ii] = 0.0;
		uu->yy[jj][ii] = 0.0;
	}
	}
	
}
