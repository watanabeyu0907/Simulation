void write_core(char File[], char Dir[], coord *uu, double **pp, double **hr, double **hm, double **hf, double **ht, double **dl, double **lm, coord *eu, coord *st, double **sd, double **lf)
{
	FILE *fp;
	char sys[256];
	int ii, jj;
	
	double posx, posy;
	double velx, vely;
	double extx, exty;
	double pres;
	double divu, dive;
	double rot;
	double phir, phim, phif, phit;
	double dlt, lvm, lvf;
	double sftx, sfty;
	double sfd;

	
	
	if((fp=fopen(File, "w")) == NULL) {
		printf("ERROR; can not open %s\n", File);
		exit(1);
	}
	printf("  writng: %s\n",File);
	
	fprintf(fp,"#dataformat, 5\n" );
	fprintf(fp,"#,%d, %d\n", para.nx, para.ny);
	fprintf(fp,"#        xx[m],          yy[m], zz[m],");
	fprintf(fp,"        uu[m/s],        vv[m/s], ww[m/s],");
	fprintf(fp,"      pp[m/s^2],      divu[1/s],       omg[1/s],");
	fprintf(fp,"        phir[-],        phim[-],        phif[-],        phit[-],");
	fprintf(fp,"         lvm[-],         dlt[-],");
	fprintf(fp,"        eu[m/s],        ev[m/s],      dive[m/s],");
	fprintf(fp,"      sx[m/s^2],      sy[m/s^2],");
	fprintf(fp,"         sfd[-],         lvf[-],");
	fprintf(fp,"\n");
	
	for (jj=0; jj<para.ny+3; jj++){
		for (ii=0; ii<para.nx+3; ii++){
			
			posx = para.dx*(ii-1)-0.5*para.lx*para.nx/(double)(para.nx-8);
			posy = para.dy*(jj-1)-0.5*para.ly*para.ny/(double)(para.ny-8);
			
			velx = +0.50*uu->xx[jj  ][ii  ]
				   +0.50*uu->xx[jj+1][ii  ];
			vely = +0.50*uu->yy[jj  ][ii  ]
				   +0.50*uu->yy[jj  ][ii+1];
				   
			extx = +0.50*eu->xx[jj  ][ii  ]
				   +0.50*eu->xx[jj+1][ii  ];
			exty = +0.50*eu->yy[jj  ][ii  ]
				   +0.50*eu->yy[jj  ][ii+1];
				   
			sftx = +0.50*st->xx[jj  ][ii  ]
				   +0.50*st->xx[jj+1][ii  ];
			sfty = +0.50*st->yy[jj  ][ii  ]
				   +0.50*st->yy[jj  ][ii+1];
				   
			pres = +0.25*pp[jj  ][ii  ]
				   +0.25*pp[jj  ][ii+1]
				   +0.25*pp[jj+1][ii  ]
				   +0.25*pp[jj+1][ii+1];
				   
			phir = +0.25*hr[jj  ][ii  ]
				   +0.25*hr[jj  ][ii+1]
				   +0.25*hr[jj+1][ii  ]
				   +0.25*hr[jj+1][ii+1];
				   
			phim = +0.25*hm[jj  ][ii  ]
				   +0.25*hm[jj  ][ii+1]
				   +0.25*hm[jj+1][ii  ]
				   +0.25*hm[jj+1][ii+1];
				   
			phif = +0.25*hf[jj  ][ii  ]
				   +0.25*hf[jj  ][ii+1]
				   +0.25*hf[jj+1][ii  ]
				   +0.25*hf[jj+1][ii+1];
				   
			phit = +0.25*ht[jj  ][ii  ]
				   +0.25*ht[jj  ][ii+1]
				   +0.25*ht[jj+1][ii  ]
				   +0.25*ht[jj+1][ii+1];
				   
			lvm  = +0.25*lm[jj  ][ii  ]
				   +0.25*lm[jj  ][ii+1]
				   +0.25*lm[jj+1][ii  ]
				   +0.25*lm[jj+1][ii+1];
				   
			lvf  = +0.25*lf[jj  ][ii  ]
				   +0.25*lf[jj  ][ii+1]
				   +0.25*lf[jj+1][ii  ]
				   +0.25*lf[jj+1][ii+1];
				   
			dlt  = +0.25*dl[jj  ][ii  ]
				   +0.25*dl[jj  ][ii+1]
				   +0.25*dl[jj+1][ii  ]
				   +0.25*dl[jj+1][ii+1];
				   
			sfd  = +0.25*sd[jj  ][ii  ]
				   +0.25*sd[jj  ][ii+1]
				   +0.25*sd[jj+1][ii  ]
				   +0.25*sd[jj+1][ii+1];
				   
			rot =-(-uu->xx[jj  ][ii  ]
				   +uu->xx[jj+1][ii  ])*para.idy
				 +(-uu->yy[jj  ][ii  ]
				   +uu->yy[jj  ][ii+1])*para.idx;
			
			if(ii==0 || ii==para.nx+2 || jj==0 || jj==para.ny+2){
				divu = 0.0;
				dive = 0.0;
			}
			else{
				divu=+0.25*(-uu->xx[jj+1][ii-1]
							-uu->xx[jj  ][ii-1]
							+uu->xx[jj+1][ii+1]
							+uu->xx[jj  ][ii+1])*para.idx
					 +0.25*(-uu->yy[jj-1][ii  ]
							-uu->yy[jj-1][ii+1]
							+uu->yy[jj+1][ii  ]
							+uu->yy[jj+1][ii+1])*para.idy;
				dive=+0.25*(-eu->xx[jj+1][ii-1]
							-eu->xx[jj  ][ii-1]
							+eu->xx[jj+1][ii+1]
							+eu->xx[jj  ][ii+1])*para.idx
					 +0.25*(-eu->yy[jj-1][ii  ]
							-eu->yy[jj-1][ii+1]
							+eu->yy[jj+1][ii  ]
							+eu->yy[jj+1][ii+1])*para.idy;
			}
									
			
			fprintf(fp, "% 8.7le, ", posx);	//  1
			fprintf(fp, "% 8.7le, ", posy);	//  2
			fprintf(fp, "    0, ");			//  3
			fprintf(fp, "% 8.7le, ", velx);	//  4
			fprintf(fp, "% 8.7le, ", vely);	//  5
			fprintf(fp, "      0, ");		//  6
			fprintf(fp, "% 8.7le, ", pres);	//  7
			fprintf(fp, "% 8.7le, ", divu);	//  8
			fprintf(fp, "% 8.7le, ", rot);	//  9
			fprintf(fp, "% 8.7le, ", phir);	// 10
			fprintf(fp, "% 8.7le, ", phim);	// 11
			fprintf(fp, "% 8.7le, ", phif);	// 12
			fprintf(fp, "% 8.7le, ", phit);	// 13
			fprintf(fp, "% 8.7le, ", lvm);	// 14
			fprintf(fp, "% 8.7le, ", dlt);	// 15
			fprintf(fp, "% 8.7le, ", extx);	// 16
			fprintf(fp, "% 8.7le, ", exty);	// 17
			fprintf(fp, "% 8.7le, ", dive);	// 18
			fprintf(fp, "% 8.7le, ", sftx);	// 19
			fprintf(fp, "% 8.7le, ", sfty);	// 20
			fprintf(fp, "% 8.7le, ", sfd);	// 21
			fprintf(fp, "% 8.7le, ", lvf);	// 22
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	

	fclose(fp);
	
	// sprintf(sys,"move %s %s\n",File, Dir);
	sprintf(sys,"mv %s %s\n",File, Dir);
	system(sys);
}

// void write_core(char File[], char Dir[], coord *uu, double **pp, double **ho, double **hi, double ***lv, coord *eu, double **ep)
// {
	// FILE *fp;
	// char sys[256];
	// int ii, jj;
	
	// double posx, posy;
	// double velx, vely;
	// double extx, exty;
	// double pres;
	// double extp;
	// double div;
	// double rot;
	// double phir, phim;
	// double lv0, lv2, lv3, lv4;

	
	
	// if((fp=fopen(File, "w")) == NULL) {
		// printf("ERROR; can not open %s\n", File);
		// exit(1);
	// }
	// printf("  writng: %s\n",File);
	
	// fprintf(fp,"#dataformat, 5\n" );
	// fprintf(fp,"#,%d, %d\n", para.nx, para.ny);
	// fprintf(fp,"#xx[m], yy[m], zz[m], uu[m/s], vv[m/s], ww[m/s], ");
	// fprintf(fp,"pp[m/s^2], phim[-], phir[-], ");
	// fprintf(fp,"lv0[-], lv1[-], lv2[-], lv3[-], lv4, [-]");
	// fprintf(fp,"\n");
	
	// for (jj=1; jj<para.ny+3; jj++){
		// for (ii=1; ii<para.nx+3; ii++){
			
			// posx = para.dx*(ii-1.5)-0.5*para.ly;
			// posy = para.dy*(jj-1.5)-0.5*para.ly;
			
			// velx = +0.50*uu->xx[jj  ][ii-1]
				   // +0.50*uu->xx[jj  ][ii  ];
			// vely = +0.50*uu->yy[jj-1][ii  ]
				   // +0.50*uu->yy[jj  ][ii  ];			
			
			// fprintf(fp, "% 8.7lf, ", posx);
			// fprintf(fp, "% 8.7lf, ", posy);
			// fprintf(fp, " 0, ");
			// fprintf(fp, "% 8.7lf, ", velx);
			// fprintf(fp, "% 8.7lf, ", vely);
			// fprintf(fp, " 0, ");
			// fprintf(fp, "% 8.7lf, ", pp[jj][ii]);
			// fprintf(fp, "% 8.7le, ", hi[jj][ii]);
			// fprintf(fp, "% 8.7le, ", ho[jj][ii]);
			// fprintf(fp, "% 8.7le, ", lv[0][jj][ii]);
			// fprintf(fp, "% 8.7le, ", lv[1][jj][ii]);
			// fprintf(fp, "% 8.7le, ", lv[2][jj][ii]);
			// fprintf(fp, "% 8.7le, ", lv[3][jj][ii]);
			// fprintf(fp, "% 8.7le, ", lv[4][jj][ii]);
			// fprintf(fp, "0, ");
			// fprintf(fp, "\n");
		// }
		// fprintf(fp, "\n");
	// }
	

	// fclose(fp);
	
	// sprintf(sys,"move %s %s\n",File, Dir);
	// system(sys);
// }

