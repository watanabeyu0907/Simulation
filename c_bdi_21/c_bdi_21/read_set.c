# include "watasys.h"
# include "watafnc.h"

void read_setting(char *File)
{
	FILE *fp;
	int ii;
	char line[1024];
	
	
	if((fp = fopen(File, "r")) == NULL) {
		printf("ERROR; can not open %s\n", File);
		exit(1);
	}
	
	for(ii=0; ii<17; ii++){
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
				para.ro = atof(strtok(NULL,","));
				break;
			case 5:
				para.mu = atof(strtok(NULL,","));
				break;
			case 6:
				para.Cr = atof(strtok(NULL,","));
				break;
			case 7:
				para.Ep = atof(strtok(NULL,","));
				break;
			case 8:
				para.nx = atoi(strtok(NULL,","));
				break;
			case 9:
				para.ny = atoi(strtok(NULL,","));
				break;
			case 10:
				para.nz = atoi(strtok(NULL,","));
				break;
			case 11:
				para.cT = atof(strtok(NULL,","));
				break;
			case 12:
				para.dT = atof(strtok(NULL,","));
				break;
			case 13:
				para.og = atof(strtok(NULL,","));
				break;
			case 14:
				para.ar = atof(strtok(NULL,","));
				break;
			case 15:
				para.Ri = atof(strtok(NULL,","));
				break;
			case 16:
				para.Ro = atof(strtok(NULL,","));
				break;
		}
	}

	fclose(fp);

	para.Re = para.wu*para.lx/(para.mu/para.ro);
	para.dx = para.lx/(double)(para.nx-8);
	para.dy = para.ly/(double)(para.ny-8);
	para.dz = para.lz/(double)(para.nz-8);
	
	para.idx = 1.0/para.dx;
	para.idy = 1.0/para.dy;
	para.idz = 1.0/para.dz;
	para.idT = 1.0/para.dT;
}

