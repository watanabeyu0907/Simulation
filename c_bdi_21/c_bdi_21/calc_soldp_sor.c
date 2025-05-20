# include "watasys.h"
# include "watafnc.h"

void calc_soldp_sor(pois_fnc *pois, double **dp, double **div, int ite, int *ans)
{
	int ii, jj, ll, mm;
	int tx[4] = {0, 1, 1, 0};
	int ty[4] = {0, 1, 0, 1};
	int sx, sy;
	int ex, ey;
	
	double tmp1, tmp2;
	double epsi;
	double sum = 1.0;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	mm=0;
	while (sum>para.Ep && mm < ite){
		sum=0.0;				
		
		for(ll=0; ll<4; ll++){
			
			sx = 2+tx[ll];
			sy = 2+ty[ll];
			
			#ifdef _OPENMP
			#pragma omp parallel for \
			default(none) \
			private(ii, jj, epsi, tmp1, tmp2) \
			shared(para, sx, sy, ex, ey, dp, div, pois) \
			reduction(+:sum)
			#endif
			for(jj=sx; jj<ey; jj+=2){
			for(ii=sy; ii<ex; ii+=2){
				
				epsi =   dp[jj  ][ii  ];
				tmp1 =  +dp[jj  ][ii-1]*pois->xm[jj][ii]
						+dp[jj  ][ii+1]*pois->xp[jj][ii]
						+dp[jj-1][ii  ]*pois->ym[jj][ii]
						+dp[jj+1][ii  ]*pois->yp[jj][ii];
				
				tmp2 = (-1.0*div[jj][ii]+tmp1)*pois->ic[jj][ii];
						
				dp[jj][ii] += para.Cr*(tmp2-dp[jj][ii]);
					
				epsi = 1.0 - dp[jj][ii]/epsi;		
				
				sum+=epsi*epsi;
			}
			}
		}
			
		boundary_pp(dp);
			
		sum=sqrt(sum/(double)(para.nx*para.ny));
		
		if(mm==0) sum=1.0;
		mm++;
	}
	// printf("iidp1 = %d \t (mgsor %d)\n", mm, num);

	*ans = mm;
	
}


