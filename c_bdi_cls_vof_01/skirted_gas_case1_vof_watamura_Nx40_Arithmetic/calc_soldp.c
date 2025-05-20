extern void calc_soldp_sor(pois_fnc *, double **, double **, int, int *);
extern void calc_dp_adj(double **, double **, double **);
extern double calc_div_rms(double **);

extern void calc_soldp_bicg(pois_fnc *, bicg_fnc *,	double **, double **, double, int, int *);
extern void bicg_00(pois_fnc *, double **, double **, double **, double **);

void calc_sol_dp(double **div, double **dpp, pois_fnc *pois, bicg_fnc *bicg, double **phir, double **phim, double **phit2, double **phit3)
{
	int ii, jj;
	
	double err[4];
	
	
	bicg_00(pois, phir, phim, phit2, phit3);
	
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
