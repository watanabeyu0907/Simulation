void write_stat(char File[], char Dir[], double **data)
{
	FILE *fp;
	int jj;

	
	if((fp=fopen(File, "w")) == NULL) {
		printf("ERROR; can not open %s\n", File);
		exit(1);
	}
	printf("  writng: %s\n",File);
	
	fprintf(fp,"#datanum, %d\n", (int)(para.cT/para.dT/st_num+1));
	fprintf(fp,"#        tt[s],    div_L2[1/s],  div_Linf[1/s],");
	fprintf(fp,"Edkdt[kg/m/s^3], Edis[kg/m/s^3], ESUw[kg/m/s^3], Eext[kg/m/s^3],");
	fprintf(fp," Egrv[kg/m/s^3], Esft[kg/m/s^3], Estp[kg/m/s^3],");
	fprintf(fp," Eall[kg/m/s^3],");
	fprintf(fp,"    fxp[kg/s^2],    fyp[kg/s^2],");
	fprintf(fp,"    fxv[kg/s^2],    fyv[kg/s^2],");
	fprintf(fp,"  Eep[kg/m/s^3],  Eev[kg/m/s^3],");
	fprintf(fp,"        phim[-],        phif[-],");
	fprintf(fp,"\n");
	
	for (jj=0; jj<para.cT*para.idT/st_num+1; jj++){
		
		fprintf(fp, "% 8.7le, ", jj*st_num*para.dT);//  1
		fprintf(fp, "% 8.7le, ", data[jj][0]);		//  2
		fprintf(fp, "% 8.7le, ", data[jj][1]);		//  3
		
		fprintf(fp, "% 8.7le, ",-data[jj][2]);		//  4
		fprintf(fp, "% 8.7le, ", data[jj][3]);		//  5
		fprintf(fp, "% 8.7le, ", data[jj][4]);		//  6
		fprintf(fp, "% 8.7le, ", data[jj][5]);		//  7
		fprintf(fp, "% 8.7le, ", data[jj][6]);		//  8
		fprintf(fp, "% 8.7le, ", data[jj][7]);		//  9
		fprintf(fp, "% 8.7le, ", data[jj][8]);		//  10
		
		fprintf(fp, "% 8.7le, ",-data[jj][2]
								+data[jj][3]
								+data[jj][4]
								+data[jj][5]
								+data[jj][6]
								+data[jj][7]
								+data[jj][8]);		// 11
		
		fprintf(fp, "% 8.7le, ", data[jj][ 9]);		// 12
		fprintf(fp, "% 8.7le, ", data[jj][10]);		// 13
		
		fprintf(fp, "% 8.7le, ", data[jj][11]);		// 14
		fprintf(fp, "% 8.7le, ", data[jj][12]);		// 15
		
		fprintf(fp, "% 8.7le, ", data[jj][13]);		// 16
		fprintf(fp, "% 8.7le, ", data[jj][14]);		// 17
		
		fprintf(fp, "% 8.7le, ", data[jj][15]);		// 18
		fprintf(fp, "% 8.7le, ", data[jj][16]);		// 19
		
		fprintf(fp, "\n");
	}
	

	fclose(fp);
	
}

