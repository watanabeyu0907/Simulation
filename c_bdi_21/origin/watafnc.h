#ifndef FUNC_H
#define FUNC_H

void read_setting(char *);

void boundary_bicg(bicg_fnc *);
void boundary_lvs(double **);
void boundary_pp(double **);
void boundary_uu(coord *);

void write_stat(char   *, char    *, double **);
void write_data(char   *,  char   *, double   , double **,
	coord   *, coord   *, double **, double **,
	double **, double **, double **, double **, 
	coord   *, coord   *);
void write_core(char   *, char    *, coord   *, double **, 
	double **, double **, double **, double **, coord *);

void monitr_div(int, double **, double **);
void monitr_phi(int, double **, double **);
void monitr_bud(int, double **, coord   *, coord   *, double **,
	double **, double **, double **, double **, coord   *, coord   *);

double calc_wei(double, double, double);
double calc_dgesv(double *, double *, double *, double *, int *, int, double, int);
void calc_basis(int, int, double *, double *);
void init_mlsmt(double *, double *);
void calc_mlsmt(int, int, int, int, double, double *, double *, double *, double *, coord *, int);
void calc_ext(coord *, coord *, double **, double **);

void calc_lvsm(double **, double **, double);
void calc_phim(double *, double **, double **, double **, double **, int);

void calc_NSeq(coord *, coord *, coord *, double **, coord *, double **, double **, coord *);
void calc_NSeq1_xx(coord *, coord *, double **, coord *, double **, double **);
void calc_NSeq1_yy(coord *, coord *, double **, coord *, double **, double **);
void calc_NSeq2_xx(coord *, coord *, coord *, double **, double **, coord *, int, double *);
void calc_NSeq2_yy(coord *, coord *, coord *, double **, double **, coord *, int, double *);

void calc_sol_dp(double **, double **, pois_fnc *, bicg_fnc *, double **, double **);
void calc_soldp_sor (pois_fnc *, double **, double **, int, int *);
void calc_soldp_bicg(pois_fnc *, bicg_fnc *, double **, double **, double, int, int *);
void bicg_00(pois_fnc *, double **, double **);
void bicg_01(pois_fnc *, double **, double **, double *);
void bicg_02(double **, double **, double, double);
void bicg_03(pois_fnc *, bicg_fnc *, double **, double **, double *);
void bicg_04(pois_fnc *, bicg_fnc *, double *);
void bicg_05(bicg_fnc *, double);
void bicg_06(pois_fnc *, bicg_fnc *, double *);
void bicg_07(bicg_fnc *, double **, double, double, double *);
void bicg_08(pois_fnc *, double **, double **, double *);
void bicg_09(bicg_fnc *, double, double);


double calc_bud1(coord *, coord   *, double **, double **);
double calc_bud2(coord *, double **, double **);
double calc_bud3(coord *);
void   calc_bud4(coord *, coord   *, double **, double **, double **, double *);

void   calc_dp_adj(double **, double **, double **);
void   calc_div(coord *, double **);
double calc_div_rms(double **);

void correct_up(coord *, double **, double **, double **, double **);


double smooth_sgn(double);
double minmod(double, double);
double get_hev(double, double *);
double get_smd(double);


void init_pp(double **);
void init_uu(coord *);
void init_vv(coord *);
void init_hev(double *);
void init_lvss(double **, double **, double **);
void init_phim(double **);
void init_phir(double **);

double calc_bud1(coord *, coord   *, double **, double **);
double calc_bud2(coord *, double **, double **);
double calc_bud3(coord *);
void   calc_bud4(coord *, coord *, double **, double **, double **, double *);

void memmove_double(double **, double **);
void memmove_coord (coord   *, coord   *);

void SMAC(char   *, char *, double  **,
	coord    *, coord    *, coord    *, coord    *,
	double  **, double  **, double  **, double  **,
	pois_fnc *, bicg_fnc *, 
	double   *,  double **, double  **, double  **, double **, 
	double  **, double  **, coord    *, coord    *);

#endif