extern void calc_NSeq0_xx(visc_fnc *, double **, double **, double **, double **);
extern void calc_NSeq0_yy(visc_fnc *, double **, double **, double **, double **);
extern void calc_NSeq1_xx(coord *, coord *, double **, coord *, double **, double **, double **, double **, coord *);
extern void calc_NSeq1_yy(coord *, coord *, double **, coord *, double **, double **, double **, double **, coord *);
extern void calc_NSeq2_xx(visc_fnc *, coord *, coord *, coord *, double **, double **, coord *, int, double *);
extern void calc_NSeq2_yy(visc_fnc *, coord *, coord *, coord *, double **, double **, coord *, int, double *);
extern void boundary_uu(coord *);

void calc_NSeq(coord *uu1, coord *uu2, coord *uu3, double **ppp, coord *adv, double **phir, double **phim, double **phit2, double **phit3, coord *vvv, visc_fnc *vfxx, visc_fnc *vfyy, coord *sft)
{
	int ii, jj;
	double sum;
	
	ii=0;
	sum=1;
	
	calc_NSeq0_xx(vfxx, phir, phim, phit2, phit3);
	calc_NSeq0_yy(vfyy, phir, phim, phit2, phit3);
	
	calc_NSeq1_xx(uu1, uu2, ppp, adv, phir, phim, phit2, phit3, sft);
	calc_NSeq1_yy(uu1, uu2, ppp, adv, phir, phim, phit2, phit3, sft);
	
	while (sum>para.Ep && ii < 100){
		sum=0;
		
		// 4-color SOR
		for(jj=0; jj<4; jj++){
			calc_NSeq2_xx(vfxx, uu2, uu3, adv, phir, phim, vvv, jj, &sum);
			calc_NSeq2_yy(vfyy, uu2, uu3, adv, phir, phim, vvv, jj, &sum);
			boundary_uu (uu3);
		}
		
		sum=sqrt(sum/(double)(para.nx*para.ny));
		
		if(ii==0) sum=1;
		ii++;
	}
	// printf("iiNS = %d\n", ii);
	
}

