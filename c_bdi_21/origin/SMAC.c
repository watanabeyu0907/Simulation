# include "watasys.h"
# include "watafnc.h"

void SMAC(  char   *File, char  *sfile, double **data,
	coord   *uu1, coord   *uu2, coord   *uu3, coord   *adv,
	double **ppp, double **dpp, double **div, double **tmp,
	pois_fnc *pois, bicg_fnc *bicg, 
	double  *hev,  double **dlt, double **lvr0, double **lvm0, double **lvm3, 
	double **phir, double **phim, coord  *vvv, coord  *euu)
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
	calc_phim(hev, phir, phim, dlt, lvm3, 0);
	
	write_data(File, fname, 0, data, 
			uu2, uu3, ppp, div,
			phir, phim, dlt, lvm3, 
			vvv, euu);
	
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
		calc_phim(hev, phir, phim, dlt, lvm3, 1);
		calc_NSeq(uu1, uu2, uu3, ppp, adv, phir, phim, vvv);
		
		// // ( 2 ) calculation of continuity eq.
		calc_div   (uu3, div);
		calc_sol_dp(div, dpp, pois, bicg, phir, phim);
		correct_up (uu3, ppp, dpp, phir, phim);
		boundary_uu(uu3);
		boundary_pp(ppp);
		
		write_data(File, fname, TT, data, 
				uu2, uu3, ppp, div,
				phir, phim, dlt, lvm3,
				vvv, euu);
		
		memmove_coord(uu1, uu2);
		memmove_coord(uu2, uu3);
		
	}
	
	printf("  %3.1f sec\t%3.1f%%\n",TT, TT/para.cT*100);

}

