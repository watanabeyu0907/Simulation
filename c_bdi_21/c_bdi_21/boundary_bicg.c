# include "watasys.h"
# include "watafnc.h"

void boundary_bicg(bicg_fnc *bicg)
{
	int ii, jj;
	int ex, ey;
	
	
	ex = para.nx+2;
	ey = para.ny+2;


	// B.C. @ x direction
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(jj) \
	shared(ey, para, bicg)
	#endif
	for(jj=2; jj<ey; jj++){
		
		bicg->r0[jj][        0] = bicg->r0[jj][        3];
		bicg->rr[jj][        0] = bicg->rr[jj][        3];
		bicg->pp[jj][        0] = bicg->pp[jj][        3];
		bicg->ee[jj][        0] = bicg->ee[jj][        3];
		bicg->ap[jj][        0] = bicg->ap[jj][        3];
		bicg->ae[jj][        0] = bicg->ae[jj][        3];
		
		bicg->r0[jj][        1] = bicg->r0[jj][        2];
		bicg->rr[jj][        1] = bicg->rr[jj][        2];
		bicg->pp[jj][        1] = bicg->pp[jj][        2];
		bicg->ee[jj][        1] = bicg->ee[jj][        2];
		bicg->ap[jj][        1] = bicg->ap[jj][        2];
		bicg->ae[jj][        1] = bicg->ae[jj][        2];
		
		bicg->r0[jj][para.nx+2] = bicg->r0[jj][para.nx+1];
		bicg->rr[jj][para.nx+2] = bicg->rr[jj][para.nx+1];
		bicg->pp[jj][para.nx+2] = bicg->pp[jj][para.nx+1];
		bicg->ee[jj][para.nx+2] = bicg->ee[jj][para.nx+1];
		bicg->ap[jj][para.nx+2] = bicg->ap[jj][para.nx+1];
		bicg->ae[jj][para.nx+2] = bicg->ae[jj][para.nx+1];
		
		bicg->r0[jj][para.nx+3] = bicg->r0[jj][para.nx+0];
		bicg->rr[jj][para.nx+3] = bicg->rr[jj][para.nx+0];
		bicg->pp[jj][para.nx+3] = bicg->pp[jj][para.nx+0];
		bicg->ee[jj][para.nx+3] = bicg->ee[jj][para.nx+0];
		bicg->ap[jj][para.nx+3] = bicg->ap[jj][para.nx+0];
		bicg->ae[jj][para.nx+3] = bicg->ae[jj][para.nx+0];
		
	}
	
	
	// B.B. @ y direction
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii) \
	shared(ex, para, bicg)
	#endif
	for(ii=2; ii<ex; ii++){
		
		bicg->r0[         0][ii] = bicg->r0[       3][ii];
		bicg->rr[         0][ii] = bicg->rr[       3][ii];
		bicg->pp[         0][ii] = bicg->pp[       3][ii];
		bicg->ee[         0][ii] = bicg->ee[       3][ii];
		bicg->ap[         0][ii] = bicg->ap[       3][ii];
		bicg->ae[         0][ii] = bicg->ae[       3][ii];
		
		bicg->r0[         1][ii] = bicg->r0[       2][ii];
		bicg->rr[         1][ii] = bicg->rr[       2][ii];
		bicg->pp[         1][ii] = bicg->pp[       2][ii];
		bicg->ee[         1][ii] = bicg->ee[       2][ii];
		bicg->ap[         1][ii] = bicg->ap[       2][ii];
		bicg->ae[         1][ii] = bicg->ae[       2][ii];
		
		bicg->r0[para.ny+2][ii] = bicg->r0[para.ny+1][ii];
		bicg->rr[para.ny+2][ii] = bicg->rr[para.ny+1][ii];
		bicg->pp[para.ny+2][ii] = bicg->pp[para.ny+1][ii];
		bicg->ee[para.ny+2][ii] = bicg->ee[para.ny+1][ii];
		bicg->ap[para.ny+2][ii] = bicg->ap[para.ny+1][ii];
		bicg->ae[para.ny+2][ii] = bicg->ae[para.ny+1][ii];
		
		bicg->r0[para.ny+3][ii] = bicg->r0[para.ny+0][ii];
		bicg->rr[para.ny+3][ii] = bicg->rr[para.ny+0][ii];
		bicg->pp[para.ny+3][ii] = bicg->pp[para.ny+0][ii];
		bicg->ee[para.ny+3][ii] = bicg->ee[para.ny+0][ii];
		bicg->ap[para.ny+3][ii] = bicg->ap[para.ny+0][ii];
		bicg->ae[para.ny+3][ii] = bicg->ae[para.ny+0][ii];
		
	}
	
}

