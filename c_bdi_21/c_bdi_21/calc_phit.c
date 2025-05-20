# include "watasys.h"
# include "watafnc.h"

void calc_phit(double **phir, double **phim, double **phif, double **phit)
{
	int ii, jj, ll, mm;
	int tx[4] = {0, 1, 1, 0};
	int ty[4] = {0, 1, 0, 1};
	int sx, sy;
	int ex, ey;
	
	double pxm, pxp, pym, pyp, pcc;
	double tmp;
	double sum = 1.0;
	double err = 0.0;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex, ey, phir, phim, phif) \
	reduction(+:err)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		err += +pow(phir[jj][ii] + phim[jj][ii], 2.0)
			   +pow(phif[jj][ii], 2.0);
	}
	}
	
	err /= (double)(para.nx*para.ny);
	
	pxm = para.idx*para.idx;
	pxp = para.idx*para.idx;
	pym = para.idy*para.idy;
	pyp = para.idy*para.idy;
	
	pcc = pxm+pxp+pym+pyp;
	
	pxm /= pcc;
	pxp /= pcc;
	pym /= pcc;
	pyp /= pcc;
	
	
	mm=0;
	while (sum>para.Ep && mm < 500){
		sum=0.0;
		
		for(ll=0; ll<4; ll++){
			
			sx = 2+tx[ll];
			sy = 2+ty[ll];
			
			#ifdef _OPENMP
			#pragma omp parallel for \
			default(none) \
			private(ii, jj, tmp) \
			shared(para, sx, sy, ex, ey) \
			shared(pxm, pxp, pym, pyp) \
			shared(phir, phim, phif, phit) \
			reduction(+:sum)
			#endif
			for(jj=sx; jj<ey; jj+=2){
			for(ii=sy; ii<ex; ii+=2){
				
				tmp =  (+phit[jj  ][ii-1]*pxm
						+phit[jj  ][ii+1]*pxp
						+phit[jj-1][ii  ]*pym
						+phit[jj+1][ii  ]*pyp)
					   *(phir[jj  ][ii  ]+phim[jj][ii])
					   + phif[jj  ][ii  ]
					    -phit[jj  ][ii  ];
				
				sum += tmp*tmp;
				phit[jj][ii] += para.Cr*tmp;
				phit[jj][ii]  = MAX(MIN(phit[jj][ii],1.0),0.0);
				
			}
			}
		}
			
		boundary_pp(phit);
			
		sum /= err;
		
		if(mm<=10) sum=1.0;
		mm++;
	}
	// printf("pt=%d  sum=%le\n", mm, sum);
	
}

