extern void boundary_uu(coord *);
extern void boundary_pp(double **);
extern void init_uu(coord *);
extern void init_vv(coord *);
extern void init_pp(double **);
extern void memmove_double(double **, double **);
extern void memmove_coord(coord *, coord *);
extern void correct_up(coord *, double **, double **, double **, double **, double **, double **);
extern void init_hev(double *);
extern void init_phir(double **);
extern void init_phim(double **);
extern void init_phif(double **, double **, double **);
extern void init_lvss(double **, double **, double **);
extern void init_lvsf(double **, double **);
extern void calc_sol_dp(double **, double **, pois_fnc *, bicg_fnc *, double **, double **, double **, double **);
extern void calc_NSeq(coord *, coord *, coord *, double **, coord *, double **, double **,double **, double **, coord *, visc_fnc *, visc_fnc *, coord *);
extern void calc_phim(double *, double **, double **, double **, double **, double **, double **, int);
extern void calc_phif(coord *, double **, double **, double **, double **, double **, coord *, coord *, coord *, double **, double **, double **);
extern void calc_phit(double **, double **, double **, double **);
extern void calc_surf(double **, double **, double **, double **, double **, coord *, coord *, coord *, coord *, coord *, coord *, double **, coord *);
extern void calc_lvsm(double **, double **, double);
extern void calc_lvsf(double **, double **, double **, double **, double**, double**, double**);


void SMAC(  char   *File, char  *sfile, double **data,
	coord   *uu1, coord   *uu2, coord   *uu3, coord   *adv,
	double **ppp, double **dpp, double **div, double **tmp,
	visc_fnc *vfxx, visc_fnc *vfyy, pois_fnc *pois,	bicg_fnc *bicg, 
	double *hev,  double **dlt, double **lvr0, double **lvm0, double **lvm3, 
	double **phir, double **phim, double **phif, double **phil, double **phit2, double **phit3,
	double **sfd1, double **sfd2, double **sfd3, double **lvf3, 
	coord *rnf, coord *rns, coord *rnfc, coord *rnfa, coord *rnsc, coord *rnsa, coord *flx, 
	coord  *vvv, coord  *euu, coord *sft)
{	
	int cn;
	char fname[256], sys[256];
	int ii;
	double TT;
	int numd=100;
	
	
	
	
	cn=strrchr(File,'.')-File;
	strncpy(fname, File, cn);
	fname[cn]='\0';
	
	sprintf(sys,"mkdir %s\n", fname);
	system(sys);
	// sprintf(sys,"copy %s %s\\.\n",sfile, fname);
	sprintf(sys,"cp %s %s/.\n",sfile, fname);
	system(sys);
	
	
	init_uu(uu1);
	init_uu(uu2);
	init_uu(uu3);
	init_vv(vvv);
	init_pp(ppp);
	init_hev(hev);
	init_phir(phir);
	init_lvss(phir, tmp, lvr0);
	init_phim(phim);
	init_lvss(phim, tmp, lvm3);
	memmove_double(lvm0, lvm3);
	
	calc_lvsm(lvm0, lvm3, 0);
	calc_phim(hev, phir, phim, phif, phit3, dlt, lvm3, 0);
	init_phif(phir, phim, phif);
	memmove_double(phit3, phif);
	calc_phif(uu3, phir, phim, phif, phil, phit3, rnf, rns, flx, sfd1, sfd2, sfd3);
	calc_phit(phir, phim, phif, phit3);
	calc_lvsf(phir, phim, phil, phit3, sfd3, tmp, lvf3);
	
	write_data(File, fname,
			0, data, 
			uu2, uu3, ppp, div,
			phir, phim, phif, phit2, phit3, dlt, lvm3, 
			vvv, euu, sft, sfd3, lvf3);
	
	puts("");
	puts("  calculating...");
	TT=0.0;
	while (TT<para.cT){
		TT+=para.dT;
		
		// printf("TT = %le\n", TT);
		for(ii=0; ii<=TT/para.cT*numd; ii++){
			if(TT/para.cT*numd-ii<para.dT/para.cT*numd){
				printf("  %3.1f sec\t%3.1f%%\n",TT, TT/para.cT*100);
			}
		}
		
		// ( 1 ) calculation of N.-S. eq.
		calc_lvsm(lvm0, lvm3, para.og*TT);
		calc_phim(hev, phir, phim, phif, phit3, dlt, lvm3, 1);
		calc_phif(uu3, phir, phim, phif, phil, phit3, rnf, rns, flx, sfd1, sfd2, sfd3);
		calc_phit(phir, phim, phif, phit3);
		calc_lvsf(phir, phim, phil, phit3, sfd3, tmp, lvf3);
		calc_surf(phir, phim, lvr0, lvm3, lvf3, rnf, rns, rnfc, rnfa, rnsc, rnsa, phit3, sft);
		calc_NSeq(uu1, uu2, uu3, ppp, adv, phir, phim, phit2, phit3, vvv, vfxx, vfyy, sft);
		
		// // ( 2 ) calculation of continuity eq.
		calc_div   (uu3, div);
		calc_sol_dp(div, dpp, pois, bicg, phir, phim, phit2, phit3);
		correct_up (uu3, ppp, dpp, phir, phim, phit2, phit3);
		boundary_uu(uu3);
		boundary_pp(ppp);
		
		write_data(File, fname,
				TT, data, 
				uu2, uu3, ppp, div,
				phir, phim, phif, phit2, phit3, dlt, lvm3,
				vvv, euu, sft, sfd3, lvf3);
		
		memmove_coord(uu1, uu2);
		memmove_coord(uu2, uu3);
		memmove_double(phit2, phit3);
		
	}
	
	printf("  %3.1f sec\t%3.1f%%\n",TT, TT/para.cT*100);

}

