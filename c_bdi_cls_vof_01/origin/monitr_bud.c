extern double calc_bud1(coord *, coord *, double **, double **, double **, double **);
extern double calc_bud2(coord *, double **, double **, double **);
extern double calc_bud3(coord *);
extern void   calc_bud4(coord *, coord *, double **, double **, double **, double **, double *);
extern double calc_bud5(coord *, double **, double **, double **);
extern double calc_bud6(coord *, double **, double **, double **, coord *);
extern double calc_bud7(coord *, double **, double **, double **, double **);
		   
void monitr_bud(int tt, double **data,
				coord  *uu2, coord *uu3, double **ppp,
				double **phir, double **phim, double **phit2, double **phit3,
				double **dlt, double **lvs, 
				coord *vvv, coord *euu, coord *surf, 
				double tht)
{
	double tmp[6];
	
	data[tt][2] = calc_bud1(uu2, uu3, phir, phim, phit2, phit3);
	data[tt][3] =+calc_bud2(uu2, phir, phim, phit2)
				 +calc_bud2(uu3, phir, phim, phit3);
	data[tt][4] =+calc_bud3(uu2)
				 +calc_bud3(uu3);
				  calc_bud4(vvv, euu, ppp, phit3, dlt, lvs, tmp);
	data[tt][5] = tmp[0] + tmp[3];
	data[tt][6] = calc_bud5(uu3, phir, phim, phit3);
	data[tt][7] = calc_bud6(uu3, phir, phim, phit3, surf);
	data[tt][8] = calc_bud7(uu3, ppp, phir, phim, phit3);
				  
	data[tt][ 9] = tmp[1];
	data[tt][10] = tmp[2];
	data[tt][11] = tmp[4];
	data[tt][12] = tmp[5];
	data[tt][13] = tmp[0];
	data[tt][14] = tmp[3];
	
}

