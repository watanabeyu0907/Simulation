# include "watasys.h"
# include "watafnc.h"

void calc_sol_dp(double **div, double **dpp, pois_fnc *pois, bicg_fnc *bicg, double **phir, double **phim)
{
	int ii, jj;
	
	double err[4];
	
	
	bicg_00(pois, phir, phim);
	
	err[0] = calc_div_rms(div);
	
	
	/* -------------------- 4-color SOR --------------------*/
	calc_soldp_sor(pois, dpp, div, 50, &ii);
	// printf("iidp1 = %d\n", ii);
	/* -------------------- 4-color SOR --------------------*/
	
	
	/* ----------------------- BiCG ----------------------- */
	calc_soldp_bicg(pois, bicg, dpp, div, err[0], 500, &ii);
	// printf("iidp2 = %d\n", ii);
	/* ----------------------- BiCG ----------------------- */
	
	
	/* -------------------- 4-color SOR --------------------*/
	if(ii>1e9)			jj = 2000;
	else if(ii>=1000)	jj = 1000;
	else				jj = 500;
	
	calc_soldp_sor(pois, dpp, div, jj, &ii);
	// printf("iidp3 = %d\n", ii);
	/* -------------------- 4-color SOR --------------------*/
	
	calc_dp_adj(dpp, phir, phim);
	
}
