// --------------------------------------------------
// 
// SMAC simulation
// 
// 2021/02/01 Tomoaki Watamura @ Osaka UNIV.
//
// --------------------------------------------------

# include "mkl.h"
# include "watasys.h"
# include "const.c"

// input
# include "read_set.c"

// // output
# include "write_data.c"
# include "write_core.c"
# include "write_stat.c"
# include "monitr_div.c"
# include "monitr_bud.c"
# include "monitr_phi.c"

// // smac core
# include "SMAC.c"

# include "memmove_double.c"
# include "memmove_coord.c"
# include "get_smd.c"
# include "get_hev.c"

# include "calc_NSeq.c"
# include "calc_NSeq0_xx.c"
# include "calc_NSeq0_yy.c"
# include "calc_NSeq1_xx.c"
# include "calc_NSeq1_yy.c"
# include "calc_NSeq2_xx.c"
# include "calc_NSeq2_yy.c"

# include "calc_ext.c"
# include "calc_lvsm.c"
# include "calc_lvsf.c"
# include "calc_phim.c"
# include "calc_phif.c"
# include "calc_phit.c"
# include "calc_surf.c"
# include "calc_nrml.c"

# include "calc_soldp.c"
# include "calc_soldp_sor.c"
# include "calc_soldp_bicg.c"
# include "calc_dp_adj.c"
# include "calc_div.c"
# include "calc_div_rms.c"
# include "correct_up.c"

# include "init_uu.c"
# include "init_vv.c"
# include "init_pp.c"
# include "init_hev.c"
# include "init_phim.c"
# include "init_phir.c"
# include "init_phif.c"
# include "init_lvss.c"
# include "init_lvsf.c"

# include "boundary_uu.c"
# include "boundary_pp.c"
# include "boundary_lvs.c"
# include "boundary_bicg.c"

// calc energy budget
# include "calc_bud1.c"
# include "calc_bud2.c"
# include "calc_bud3.c"
# include "calc_bud4.c"
# include "calc_bud5.c"
# include "calc_bud6.c"
# include "calc_bud7.c"


int main(int argc, char *argv[])
{
	char out[256]="test2d.csv";
	char set[256]="set_cavity.csv";
	char ch, *tch;

	int count=0;
	
	coord  uu1, uu2, uu3; // at cell surface
	coord  adv;
	coord  vvv;
	coord  euu;
	coord  sft;

	double **ppp, **dpp, **div, **tmp; // at cell centre
	double **phim, **phir, **phif, **phil;
	double **phit2, **phit3;
	double **sfd1, **sfd2, **sfd3;
	double **dlt;
	double **lvm0, **lvm3;
	double **lvs0, **lvf3;
	
	// normal vector
	coord rnf, rns;		// at cell surface
	coord rnfc, rnsc;	// at cell centre
	coord rnfa, rnsa;	// at cell apex
	coord flx;
	
	pois_fnc pois;
	bicg_fnc bicg;
	visc_fnc vfxx, vfyy;
	
	double **data;
	double *hev;

	
	time_t time1, time2;
	double time3, time4;
	
	#ifdef _OPENMP
	time3 = omp_get_wtime();
	#else
	time1 = clock();
	#endif


	
	if (argc > 1){
		count=0;
		do {
			count++;
			tch = argv[count];
			if (*argv[count] == '-'){
				ch = *(++tch);
				switch (ch){
					case 's':
						count++;
						if(argv[count] != NULL){
							strcpy(set, argv[count]);
							printf("  -s  read parameter\t: %s\n", set);
						}
						break;
						
					case 'o':
						count++;
						if(argv[count] != NULL){
							strcpy(out, argv[count]);
							printf("  -o  read parameter\t: %s\n", out);
						}
						break;
						
					case 'f':
						count++;
						if(argv[count] != NULL){
							fpsnum = atoi(argv[count]);
							printf("  -f  read parameter\t: %d\n", fpsnum);
						}
						break;
						
					case 'n':
						count++;
						if(argv[count] != NULL){
							st_num = atoi(argv[count]);
							printf("  -n  read parameter\t: %d\n", st_num);
						}
						break;
						
					default:
						printf("error on command");
						exit(0);
				}
			}
		}while (count < argc - 1);
	}
	else{
		printf("  cavity.exe -s [setting] -o [output]\n");
		printf("  -s : setting file (default: %s)\n",set);
		printf("  -o : output file (default: %s)\n", out);
		printf("  -f : output frame/sec. (default: %d)\n", fpsnum);
		printf("  -n : output statistic interval (default: %d)\n", st_num);
		exit(-1);
	}
	
	puts("");
	#ifdef _OPENMP
	printf("  thread num in this machine: %d\n", omp_get_max_threads());
	puts("");
	#endif

	
	// read setting file
	read_setting(set);
	// cpypara();
	// set_gpu();
	
	
	printf("  lx = %8.3le\n",para.lx);
	printf("  ly = %8.3le\n",para.ly);
	printf("  lz = %8.3le\n",para.lz);
	printf("  wu = %8.3le\n",para.wu);
	printf("  nu1= %8.3le\n",para.mu1/para.ro1);
	printf("  nu2= %8.3le\n",para.mu2/para.ro2);
	printf("  gg = %8.3le\n",para.gg);
	printf("  si = %8.3le\n",para.si);
	printf("  ca = %8.3le\n",para.si);
	printf("  Re = %8.2le\n",para.Re);
	printf("  Cr = %8.2le\n",para.Cr);
	printf("  Ep = %8.2le\n",para.Ep);
	printf("  nx = %8d\n",para.nx);
	printf("  ny = %8d\n",para.ny);
	printf("  nz = %8d\n",para.nz);
	printf("  dx = %5.2le\n",para.dx);
	printf("  dy = %5.2le\n",para.dy);
	printf("  dz = %5.2le\n",para.dz);
	printf("  dt = %5.2le\n",para.dT);
	printf("  ct = %8.3le\n",para.cT);
	printf("  og = %8.3le\n",para.og);
	printf("  ar = %8.3le\n",para.ar);
	printf("  Ri = %8.3le\n",para.Ri);
	printf("  Ro = %8.3le\n",para.Ro);
	puts("");
	

	// making array
	array1D(hev, hevnum*4+1, double);
	
	array2D(uu1.xx, para.ny+4, para.nx+4, double);
	array2D(uu1.yy, para.ny+4, para.nx+4, double);
	array2D(uu2.xx, para.ny+4, para.nx+4, double);
	array2D(uu2.yy, para.ny+4, para.nx+4, double);
	array2D(uu3.xx, para.ny+4, para.nx+4, double);
	array2D(uu3.yy, para.ny+4, para.nx+4, double);
	array2D(adv.xx, para.ny+4, para.nx+4, double);
	array2D(adv.yy, para.ny+4, para.nx+4, double);
	array2D(vvv.xx, para.ny+4, para.nx+4, double);
	array2D(vvv.yy, para.ny+4, para.nx+4, double);
	array2D(euu.xx, para.ny+4, para.nx+4, double);
	array2D(euu.yy, para.ny+4, para.nx+4, double);
	array2D(sft.xx, para.ny+4, para.nx+4, double);
	array2D(sft.yy, para.ny+4, para.nx+4, double);
	array2D(rnf.xx, para.ny+4, para.nx+4, double);
	array2D(rnf.yy, para.ny+4, para.nx+4, double);
	array2D(rns.xx, para.ny+4, para.nx+4, double);
	array2D(rns.yy, para.ny+4, para.nx+4, double);
	array2D(rnfc.xx, para.ny+4, para.nx+4, double);
	array2D(rnfc.yy, para.ny+4, para.nx+4, double);
	array2D(rnfa.xx, para.ny+4, para.nx+4, double);
	array2D(rnfa.yy, para.ny+4, para.nx+4, double);
	array2D(rnsc.xx, para.ny+4, para.nx+4, double);
	array2D(rnsc.yy, para.ny+4, para.nx+4, double);
	array2D(rnsa.xx, para.ny+4, para.nx+4, double);
	array2D(rnsa.yy, para.ny+4, para.nx+4, double);
	array2D(flx.xx, para.ny+4, para.nx+4, double);
	array2D(flx.yy, para.ny+4, para.nx+4, double);
	array2D(ppp, para.ny+4, para.nx+4, double);
	array2D(tmp, para.ny+4, para.nx+4, double);
	array2D(dpp, para.ny+4, para.nx+4, double);
	array2D(div, para.ny+4, para.nx+4, double);
	array2D(dlt, para.ny+4, para.nx+4, double);
	array2D(lvs0, para.ny+4, para.nx+4, double);
	array2D(lvm0, para.ny+4, para.nx+4, double);
	array2D(lvm3, para.ny+4, para.nx+4, double);
	array2D(lvf3, para.ny+4, para.nx+4, double);
	array2D(sfd1, para.ny+4, para.nx+4, double);
	array2D(sfd2, para.ny+4, para.nx+4, double);
	array2D(sfd3, para.ny+4, para.nx+4, double);
	array2D(phim, para.ny+4, para.nx+4, double);
	array2D(phir, para.ny+4, para.nx+4, double);
	array2D(phif, para.ny+4, para.nx+4, double);
	array2D(phil, para.ny+4, para.nx+4, double);
	array2D(phit2, para.ny+4, para.nx+4, double);
	array2D(phit3, para.ny+4, para.nx+4, double);
	
	array2D(pois.xm, para.ny+4, para.nx+4, double);
	array2D(pois.xp, para.ny+4, para.nx+4, double);
	array2D(pois.ym, para.ny+4, para.nx+4, double);
	array2D(pois.yp, para.ny+4, para.nx+4, double);
	array2D(pois.cc, para.ny+4, para.nx+4, double);
	array2D(pois.ic, para.ny+4, para.nx+4, double);
	
	array2D(vfxx.xmxm, para.ny+4, para.nx+4, double);
	array2D(vfxx.xpxp, para.ny+4, para.nx+4, double);
	array2D(vfxx.ymym, para.ny+4, para.nx+4, double);
	array2D(vfxx.ypyp, para.ny+4, para.nx+4, double);
	array2D(vfxx.ymxm, para.ny+4, para.nx+4, double);
	array2D(vfxx.ymxp, para.ny+4, para.nx+4, double);
	array2D(vfxx.ypxm, para.ny+4, para.nx+4, double);
	array2D(vfxx.ypxp, para.ny+4, para.nx+4, double);
	array2D(vfxx.ccxx, para.ny+4, para.nx+4, double);
	array2D(vfxx.ccyy, para.ny+4, para.nx+4, double);
	
	array2D(vfyy.xmxm, para.ny+4, para.nx+4, double);
	array2D(vfyy.xpxp, para.ny+4, para.nx+4, double);
	array2D(vfyy.ymym, para.ny+4, para.nx+4, double);
	array2D(vfyy.ypyp, para.ny+4, para.nx+4, double);
	array2D(vfyy.ymxm, para.ny+4, para.nx+4, double);
	array2D(vfyy.ymxp, para.ny+4, para.nx+4, double);
	array2D(vfyy.ypxm, para.ny+4, para.nx+4, double);
	array2D(vfyy.ypxp, para.ny+4, para.nx+4, double);
	array2D(vfyy.ccxx, para.ny+4, para.nx+4, double);
	array2D(vfyy.ccyy, para.ny+4, para.nx+4, double);
	
	array2D(bicg.r0, para.ny+4, para.nx+4, double);
	array2D(bicg.rr, para.ny+4, para.nx+4, double);
	array2D(bicg.pp, para.ny+4, para.nx+4, double);
	array2D(bicg.ee, para.ny+4, para.nx+4, double);
	array2D(bicg.ap, para.ny+4, para.nx+4, double);
	array2D(bicg.ae, para.ny+4, para.nx+4, double);
	
	array2D(data, para.cT/para.dT/st_num+1, 20, double);
	
	
	// // SMAC method
	SMAC(out, set, data, 
		&uu1, &uu2,  &uu3, &adv,
		ppp, dpp, div, tmp, 
		&vfxx, &vfyy, &pois, &bicg, 
		hev, dlt, lvs0, lvm0, lvm3, 
		phir, phim, phif, phil, phit2, phit3,
		sfd1,  sfd2, sfd3, lvf3, 
		&rnf, &rns, &rnfc, &rnfa, &rnsc, &rnsa, &flx, 
		&vvv,  &euu, &sft);
	
	
	free1D(hev, hevnum*4+1);
	
	free2D(uu1.xx, para.ny+4, para.nx+4);
	free2D(uu1.yy, para.ny+4, para.nx+4);
	free2D(uu2.xx, para.ny+4, para.nx+4);
	free2D(uu2.yy, para.ny+4, para.nx+4);
	free2D(uu3.xx, para.ny+4, para.nx+4);
	free2D(uu3.yy, para.ny+4, para.nx+4);
	free2D(adv.xx, para.ny+4, para.nx+4);
	free2D(adv.yy, para.ny+4, para.nx+4);
	free2D(vvv.xx, para.ny+4, para.nx+4);
	free2D(vvv.yy, para.ny+4, para.nx+4);
	free2D(euu.xx, para.ny+4, para.nx+4);
	free2D(euu.yy, para.ny+4, para.nx+4);
	free2D(sft.xx, para.ny+4, para.nx+4);
	free2D(sft.yy, para.ny+4, para.nx+4);
	free2D(ppp, para.ny+4, para.nx+4);
	free2D(tmp, para.ny+4, para.nx+4);
	free2D(dpp, para.ny+4, para.nx+4);
	free2D(div, para.ny+4, para.nx+4);
	free2D(dlt, para.ny+4, para.nx+4);
	free2D(lvs0, para.ny+4, para.nx+4);
	free2D(lvm0, para.ny+4, para.nx+4);
	free2D(lvm3, para.ny+4, para.nx+4);
	free2D(lvf3, para.ny+4, para.nx+4);
	free2D(sfd1, para.ny+4, para.nx+4);
	free2D(sfd2, para.ny+4, para.nx+4);
	free2D(sfd3, para.ny+4, para.nx+4);
	free2D(rnf.xx, para.ny+4, para.nx+4);
	free2D(rnf.yy, para.ny+4, para.nx+4);
	free2D(rns.xx, para.ny+4, para.nx+4);
	free2D(rns.yy, para.ny+4, para.nx+4);
	free2D(rnfc.xx, para.ny+4, para.nx+4);
	free2D(rnfc.yy, para.ny+4, para.nx+4);
	free2D(rnfa.xx, para.ny+4, para.nx+4);
	free2D(rnfa.yy, para.ny+4, para.nx+4);
	free2D(rnsc.xx, para.ny+4, para.nx+4);
	free2D(rnsc.yy, para.ny+4, para.nx+4);
	free2D(rnsa.xx, para.ny+4, para.nx+4);
	free2D(rnsa.yy, para.ny+4, para.nx+4);
	free2D(flx.xx, para.ny+4, para.nx+4);
	free2D(flx.yy, para.ny+4, para.nx+4);
	free2D(phim, para.ny+4, para.nx+4);
	free2D(phir, para.ny+4, para.nx+4);
	free2D(phif, para.ny+4, para.nx+4);
	free2D(phil, para.ny+4, para.nx+4);
	free2D(phit2, para.ny+4, para.nx+4);
	free2D(phit3, para.ny+4, para.nx+4);
	
	free2D(pois.xm, para.ny+4, para.nx+4);
	free2D(pois.xp, para.ny+4, para.nx+4);
	free2D(pois.ym, para.ny+4, para.nx+4);
	free2D(pois.yp, para.ny+4, para.nx+4);
	free2D(pois.cc, para.ny+4, para.nx+4);
	free2D(pois.ic, para.ny+4, para.nx+4);
	
	free2D(vfxx.xmxm, para.ny+4, para.nx+4);
	free2D(vfxx.xpxp, para.ny+4, para.nx+4);
	free2D(vfxx.ymym, para.ny+4, para.nx+4);
	free2D(vfxx.ypyp, para.ny+4, para.nx+4);
	free2D(vfxx.ymxm, para.ny+4, para.nx+4);
	free2D(vfxx.ypxm, para.ny+4, para.nx+4);
	free2D(vfxx.ymxp, para.ny+4, para.nx+4);
	free2D(vfxx.ypxp, para.ny+4, para.nx+4);
	free2D(vfxx.ccxx, para.ny+4, para.nx+4);
	free2D(vfxx.ccyy, para.ny+4, para.nx+4);
	
	free2D(vfyy.xmxm, para.ny+4, para.nx+4);
	free2D(vfyy.xpxp, para.ny+4, para.nx+4);
	free2D(vfyy.ymym, para.ny+4, para.nx+4);
	free2D(vfyy.ypyp, para.ny+4, para.nx+4);
	free2D(vfyy.ymxm, para.ny+4, para.nx+4);
	free2D(vfyy.ypxm, para.ny+4, para.nx+4);
	free2D(vfyy.ymxp, para.ny+4, para.nx+4);
	free2D(vfyy.ypxp, para.ny+4, para.nx+4);
	free2D(vfyy.ccxx, para.ny+4, para.nx+4);
	free2D(vfyy.ccyy, para.ny+4, para.nx+4);
	
	free2D(bicg.r0, para.ny+4, para.nx+4);
	free2D(bicg.rr, para.ny+4, para.nx+4);
	free2D(bicg.pp, para.ny+4, para.nx+4);
	free2D(bicg.ee, para.ny+4, para.nx+4);
	free2D(bicg.ap, para.ny+4, para.nx+4);
	free2D(bicg.ae, para.ny+4, para.nx+4);
	
	free2D(data, para.cT/para.dT/st_num+1, 20);

	
	puts("");
	#ifdef _OPENMP
	time4 = omp_get_wtime();
	printf("  processing time: %f [sec]\n", (double)(time4-time3));
	#else
	time2 = clock();
	printf("  processing time: %f [sec]\n", (double)(time2-time1)/CLOCKS_PER_SEC);
	#endif
	puts("  done !!");
	
	return (0);
}



