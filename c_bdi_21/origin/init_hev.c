# include "watasys.h"
# include "watafnc.h"

void init_hev(double *hev)
{
	int	ii;
	double aa;
	double *smd;
	
	array1D(smd, hevnum*4+1, double);

	for(ii=0; ii<hevnum*4+1; ii++){
		
		smd[ii] = get_smd(ii/(double)(hevnum)-2.0);
		
	}
	
	aa = 0.0;
	for(ii=0; ii<hevnum*4; ii++){
		
		aa += 0.5*(smd[ii]+smd[ii+1])/(double)hevnum;
		hev[ii+1] = aa;
		
	}
	
	free1D(smd, hevnum*4+1);
	
}

