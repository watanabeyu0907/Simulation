// --------------------------------------------------
// 
// SMAC simulation
// 
// 2021/02/01 Tomoaki Watamura @ Osaka UNIV.
//
// --------------------------------------------------

# include "watasys.h"
# include "watafnc.h"


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

	double **ppp, **dpp, **div, **tmp; // at cell centre
	double **phim, **phir;
	double **dlt;
	double **lvs0, **lvm0, **lvm3;
	
	
	pois_fnc pois;
	bicg_fnc bicg;
	
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
	printf("  ro = %8.3le\n",para.ro);
	printf("  nu = %8.3le\n",para.mu);
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
	array2D(ppp, para.ny+4, para.nx+4, double);
	array2D(tmp, para.ny+4, para.nx+4, double);
	array2D(dpp, para.ny+4, para.nx+4, double);
	array2D(div, para.ny+4, para.nx+4, double);
	array2D(dlt, para.ny+4, para.nx+4, double);
	array2D(lvs0, para.ny+4, para.nx+4, double);
	array2D(lvm0, para.ny+4, para.nx+4, double);
	array2D(lvm3, para.ny+4, para.nx+4, double);
	array2D(phim, para.ny+4, para.nx+4, double);
	array2D(phir, para.ny+4, para.nx+4, double);
	
	array2D(pois.xm, para.ny+4, para.nx+4, double);
	array2D(pois.xp, para.ny+4, para.nx+4, double);
	array2D(pois.ym, para.ny+4, para.nx+4, double);
	array2D(pois.yp, para.ny+4, para.nx+4, double);
	array2D(pois.cc, para.ny+4, para.nx+4, double);
	array2D(pois.ic, para.ny+4, para.nx+4, double);
	
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
		&pois, &bicg, 
		hev, dlt, lvs0, lvm0, lvm3, 
		phir, phim, &vvv,  &euu);
	
	
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
	free2D(ppp, para.ny+4, para.nx+4);
	free2D(tmp, para.ny+4, para.nx+4);
	free2D(dpp, para.ny+4, para.nx+4);
	free2D(div, para.ny+4, para.nx+4);
	free2D(dlt, para.ny+4, para.nx+4);
	free2D(lvs0, para.ny+4, para.nx+4);
	free2D(lvm0, para.ny+4, para.nx+4);
	free2D(lvm3, para.ny+4, para.nx+4);
	free2D(phim, para.ny+4, para.nx+4);
	free2D(phir, para.ny+4, para.nx+4);
	
	free2D(pois.xm, para.ny+4, para.nx+4);
	free2D(pois.xp, para.ny+4, para.nx+4);
	free2D(pois.ym, para.ny+4, para.nx+4);
	free2D(pois.yp, para.ny+4, para.nx+4);
	free2D(pois.cc, para.ny+4, para.nx+4);
	free2D(pois.ic, para.ny+4, para.nx+4);
	
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



