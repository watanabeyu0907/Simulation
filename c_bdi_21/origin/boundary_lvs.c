# include "watasys.h"
# include "watafnc.h"

void boundary_lvs(double **lvs)
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
	shared(para, ey, lvs)
	#endif
	for(jj=2; jj<ey; jj++){
		lvs[jj][        0] = 3.0*lvs[jj][        2] - 2.0*lvs[jj][        3];
		lvs[jj][        1] = 2.0*lvs[jj][        2] - 1.0*lvs[jj][        3];
		lvs[jj][para.nx+2] = 2.0*lvs[jj][para.nx+1] - 1.0*lvs[jj][para.nx  ];
		lvs[jj][para.nx+3] = 3.0*lvs[jj][para.nx+1] - 2.0*lvs[jj][para.nx  ];
	}
	
// B.C. @ y direction
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii) \
	shared(para, ex, lvs)
	#endif
	for(ii=2; ii<ex; ii++){
		lvs[        0][ii] = 3.0*lvs[        2][ii] - 2.0*lvs[        3][ii];
		lvs[        1][ii] = 2.0*lvs[        2][ii] - 1.0*lvs[        3][ii];
		lvs[para.ny+2][ii] = 2.0*lvs[para.ny+1][ii] - 1.0*lvs[para.ny  ][ii];
		lvs[para.ny+3][ii] = 3.0*lvs[para.ny+1][ii] - 2.0*lvs[para.ny  ][ii];
	}
	
	lvs[        0][        0] = 3.0*lvs[        2][        2] - 2.0*lvs[        3][        3];
	lvs[        1][        1] = 2.0*lvs[        2][        2] - 1.0*lvs[        3][        3];
	lvs[        0][        1] = 0.5*lvs[        0][        0] + 0.5*lvs[        0][        2];
	lvs[        1][        0] = 0.5*lvs[        0][        0] + 0.5*lvs[        2][        0];
	
	lvs[para.ny+3][        0] = 3.0*lvs[para.ny+1][        2] - 2.0*lvs[para.ny  ][        3];
	lvs[para.ny+2][        1] = 2.0*lvs[para.ny+1][        2] - 1.0*lvs[para.ny  ][        3];
	lvs[para.ny+3][        1] = 0.5*lvs[para.ny+3][        0] + 0.5*lvs[para.ny+3][        2];
	lvs[para.ny+2][        0] = 0.5*lvs[para.ny+1][        0] + 0.5*lvs[para.ny+3][        0];
	
	lvs[        1][para.nx+2] = 2.0*lvs[        2][para.nx+1] - 1.0*lvs[        3][para.nx  ];
	lvs[        0][para.nx+3] = 3.0*lvs[        2][para.nx+1] - 2.0*lvs[        3][para.nx  ];
	lvs[        1][para.nx+3] = 0.5*lvs[        0][para.nx+3] + 0.5*lvs[        2][para.nx+3];
	lvs[        0][para.nx+2] = 0.5*lvs[        0][para.nx+1] + 0.5*lvs[        0][para.nx+3];
	
	lvs[para.ny+2][para.nx+2] = 2.0*lvs[para.ny+1][para.nx+1] - 1.0*lvs[para.ny  ][para.nx  ];
	lvs[para.ny+3][para.nx+3] = 3.0*lvs[para.ny+1][para.nx+1] - 2.0*lvs[para.ny  ][para.nx  ];
	lvs[para.ny+3][para.nx+2] = 0.5*lvs[para.ny+3][para.nx+1] + 0.5*lvs[para.ny+3][para.nx+3];
	lvs[para.ny+2][para.nx+3] = 0.5*lvs[para.ny+1][para.nx+3] + 0.5*lvs[para.ny+3][para.nx+3];
	
}

