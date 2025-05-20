# include "watasys.h"
# include "watafnc.h"
		   
void monitr_bud(
	int tt, double **data, coord  *uu2, coord *uu3, double **ppp,
	double **phir, double **phim, double **dlt, double **lvs, coord *vvv, coord *euu)
{
	double tmp[6];
	
	data[tt][2] = calc_bud1(uu2, uu3, phir, phim);
	data[tt][3] =+calc_bud2(uu2, phir, phim)
				 +calc_bud2(uu3, phir, phim);
	data[tt][4] =+calc_bud3(uu2)
				 +calc_bud3(uu3);
				  calc_bud4(vvv, euu, ppp, dlt, lvs, tmp);
	data[tt][5] = tmp[0] + tmp[3];
	
	data[tt][ 6] = tmp[0];
	data[tt][ 7] = tmp[1];
	data[tt][ 8] = tmp[2];
	data[tt][ 9] = tmp[3];
	data[tt][10] = tmp[4];
	data[tt][11] = tmp[5];
	
}

