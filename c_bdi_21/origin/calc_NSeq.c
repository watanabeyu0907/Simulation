# include "watasys.h"
# include "watafnc.h"

void calc_NSeq(coord *uu1, coord *uu2, coord *uu3, double **ppp, coord *adv, double **phir, double **phim, coord *vvv)
{
	int ii, jj;
	double sum;
	
	ii=0;
	sum=1;
	
	calc_NSeq1_xx(uu1, uu2, ppp, adv, phir, phim);
	calc_NSeq1_yy(uu1, uu2, ppp, adv, phir, phim);
	
	while (sum>para.Ep && ii < 100){
		sum=0;
		
		// 4-color SOR
		for(jj=0; jj<4; jj++){
			calc_NSeq2_xx(uu2, uu3, adv, phir, phim, vvv, jj, &sum);
			calc_NSeq2_yy(uu2, uu3, adv, phir, phim, vvv, jj, &sum);
			boundary_uu (uu3);
		}
		
		sum=sqrt(sum/(double)(para.nx*para.ny));
		
		if(ii==0) sum=1;
		ii++;
	}
	// printf("iiNS = %d\n", ii);
	
}

