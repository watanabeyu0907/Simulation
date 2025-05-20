# include "watasys.h"
# include "watafnc.h"

void memmove_coord(coord *u1, coord *u2)
{
	int jj;
	
	#ifdef _OPENMP
	#pragma omp parallel for 
	#endif
	for(jj=0; jj<para.ny+4; jj++){
		memmove(u1->xx[jj], u2->xx[jj], (para.nx+4)*sizeof(double));
		memmove(u1->yy[jj], u2->yy[jj], (para.nx+4)*sizeof(double));
	}
	
}
