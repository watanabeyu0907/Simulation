void read_setting(char *File)
{
	FILE *fp;
	int ii;
	char line[1024];
	
	
	if((fp = fopen(File, "r")) == NULL) {
		printf("ERROR; can not open %s\n", File);
		exit(1);
	}
	
	for(ii=0; ii<22; ii++){
		fgets(line, sizeof(line), fp);
		strtok(line,",");
		switch (ii){
			case 0:
				para.lx = atof(strtok(NULL,","));
				break;
			case 1:
				para.ly = atof(strtok(NULL,","));
				break;
			case 2:
				para.lz = atof(strtok(NULL,","));
				break;
			case 3:
				para.wu = atof(strtok(NULL,","));
				break;
			case 4:
				para.ro1= atof(strtok(NULL,","));
				break;
			case 5:
				para.ro2= atof(strtok(NULL,","));
				break;
			case 6:
				para.mu1= atof(strtok(NULL,","));
				break;
			case 7:
				para.mu2= atof(strtok(NULL,","));
				break;
			case 8:
				para.gg = atof(strtok(NULL,","));
				break;
			case 9:
				para.si = atof(strtok(NULL,","));
				break;
			case 10:
				para.ca = atof(strtok(NULL,","));
				break;
			case 11:
				para.Cr = atof(strtok(NULL,","));
				break;
			case 12:
				para.Ep = atof(strtok(NULL,","));
				break;
			case 13:
				para.nx = atoi(strtok(NULL,","));
				break;
			case 14:
				para.ny = atoi(strtok(NULL,","));
				break;
			case 15:
				para.nz = atoi(strtok(NULL,","));
				break;
			case 16:
				para.cT = atof(strtok(NULL,","));
				break;
			case 17:
				para.dT = atof(strtok(NULL,","));
				break;
			case 18:
				para.og = atof(strtok(NULL,","));
				break;
			case 19:
				para.ar = atof(strtok(NULL,","));
				break;
			case 20:
				para.Ri = atof(strtok(NULL,","));
				break;
			case 21:
				para.Ro = atof(strtok(NULL,","));
				break;
		}
	}

	fclose(fp);

	para.Re = para.wu*para.lx/(para.mu1/para.ro1);
	para.dx = para.lx/(double)(para.nx-0);	
	para.dy = para.ly/(double)(para.ny-0);
	para.dz = para.lz/(double)(para.nz-0);
	
	para.idx = 1.0/para.dx;	
	para.idy = 1.0/para.dy;
	para.idz = 1.0/para.dz;
	para.idT = 1.0/para.dT;
	
	para.bet = 2.0;
	para.gp1 = -0.5/sqrt(3.0);
	para.gp2 =  0.5/sqrt(3.0);
	para.idc =  1.0/sqrt((pow(para.dx,2.0)+pow(para.dy,2.0))*0.5);
	para.eps_mtc = para.bet/(2.0*pow(cosh(-para.bet*1.5),2.0));
	para.eps_trn = 0.5*(tanh(-para.bet*2.5)+1.0);
}

