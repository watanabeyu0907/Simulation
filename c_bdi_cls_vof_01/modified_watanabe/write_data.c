extern void write_core(char *, char *, coord *, double **, double **, double **, double **, double **, double **, double **, coord *, coord *, double **, double **);
extern void monitr_phi(int, double **, double **, int);
extern void monitr_div(int, double **, double **);
extern void monitr_bud(int, double **, coord *, coord *, double **, double **, double **, double **, double **, double **, double **, coord *, coord *, coord *, double);
extern void write_stat(char *, char *, double **);
extern void calc_div(coord *, double **);
extern void calc_lvm3(double **, double **, double);
extern void calc_ext(coord *, coord *, double **, double **);

void write_data(char   File[], char   fname[],
				double TT,     double **data,
				coord   *uu2,  coord   *uu3, 
				double **ppp,  double **div,
				double **phir, double **phim, double **phif, double **phit2, double **phit3, 
				double **dlt, double **lvm3,
				coord   *vvv,  coord   *euu, coord *sft, double **sfd, double **lvf3)
{
	char out[256];
	char sys[256];
	int ii, num, flg;
	
	// fps
	num=para.cT*fpsnum;
	
	flg=0;
	for(ii=0; ii<=TT/para.cT*num; ii++){
		if(TT/para.cT*num-ii<para.dT/para.cT*num
		|| TT==para.cT){
			
			flg=1;
			calc_ext(uu3, euu, phim, lvm3);
			
			sprintf(out, "%s_%03ds_%06d.csv", fname, (int)para.cT, ii);
			write_core(out, fname, uu3, ppp, phir, phim, phif, phit3, dlt, lvm3, euu, sft, sfd, lvf3);
			break;
		}
	}
	
	if(((int)(TT/para.dT+0.5))%st_num == 0){
		
		calc_div(uu3, div);
		
		if(flg==0){
			calc_ext(uu3, euu, phim, lvm3);
		}		
		monitr_phi((int)(TT/para.dT+0.5)/st_num, data, phim, 0);
		monitr_phi((int)(TT/para.dT+0.5)/st_num, data, phif, 1);
		monitr_div((int)(TT/para.dT+0.5)/st_num, data, div);
		monitr_bud((int)(TT/para.dT+0.5)/st_num, data, uu2, uu3, ppp, phir, phim, phit2, phit3, dlt, lvm3, vvv, euu, sft, para.og*TT);
		
		if(((int)(TT/para.dT+0.5))%(st_num*5) == 0
		|| TT == para.dT){
			
			sprintf(out, "%s_%03ds_stat.csv", fname, (int)para.cT);
			write_stat(out, fname, data);
			// sprintf(sys,"move %s %s\n",out, fname);
			sprintf(sys,"mv %s %s\n",out, fname);
			system(sys);
		}
	}
	
}

