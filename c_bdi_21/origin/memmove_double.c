# include "watasys.h"
# include "watafnc.h"

void memmove_double(double **d1, double **d2)
{
	int jj;
	
	#ifdef _OPENMP
	#pragma omp parallel for 
	#endif
	for(jj=0; jj<para.ny+4; jj++){
		memmove(d1[jj], d2[jj], (para.nx+4)*sizeof(double));
	}
	
}
