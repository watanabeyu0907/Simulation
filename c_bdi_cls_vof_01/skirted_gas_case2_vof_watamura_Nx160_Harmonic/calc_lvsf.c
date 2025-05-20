extern void boundary_lvs(double **);
extern double smooth_sgn(double );
extern double minmod(double, double);

void calc_lvsf(double **phir, double **phim, double **phil, double **phit, double **sfd, double **tmp, double **lvs3)
{
	int ii, jj, kk;
	int cnt;
	int ex2, ey2;
	int ex3, ey3;
	double eps, sum;
	double sgn;
	double gdn,      gdp;
	double ddn, ddc, ddp;
	double ggx, ggy;
	double qqq;
	
	
	ex2 = para.nx+2;
	ey2 = para.ny+2;
	ex3 = para.nx+3;
	ey3 = para.ny+3;
	
	
	init_lvsf(phit, tmp);
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, qqq) \
	shared(ex3, ey3, phir, phim, phil) 
	#endif
	for(jj=1; jj<ey3; jj++){
	for(ii=1; ii<ex3; ii++){
			
		qqq = phir[jj-1][ii-1]+phir[jj-1][ii  ]+phir[jj-1][ii+1]
			 +phir[jj  ][ii-1]+phir[jj  ][ii  ]+phir[jj  ][ii+1]
			 +phir[jj+1][ii-1]+phir[jj+1][ii  ]+phir[jj+1][ii+1]
			 +phim[jj-1][ii-1]+phim[jj-1][ii  ]+phim[jj-1][ii+1]
			 +phim[jj  ][ii-1]+phim[jj  ][ii  ]+phim[jj  ][ii+1]
			 +phim[jj+1][ii-1]+phim[jj+1][ii  ]+phim[jj+1][ii+1];
		
		phil[jj][ii] = qqq;
		
	}
	}

	
	
	kk = 0;
	sum= 1;
	// narrow band levelset in 2nd order upwind ENO scheme
	while (sum>para.Ep && kk < 200){
		
		sum = 0.0;
		cnt = 0;
		
		#ifdef _OPENMP
		#pragma omp parallel for \
		default(none) \
		private(ii, jj, kk, eps, sgn) \
		private(ddn, ddc, ddp, gdn, gdp, ggx, ggy, qqq) \
		shared(phil, phit, sfd, tmp, lvs3) \
		shared(para, ex2, ey2) \
		reduction(+:sum, cnt)
		#endif
		for(jj=2; jj<ey2; jj++){
		for(ii=2; ii<ex2; ii++){
			
			sgn = smooth_sgn(lvs3[jj][ii]);
			
			ddn = tmp[jj][ii-2] - 2.0*tmp[jj][ii-1] + tmp[jj][ii  ];
			ddc = tmp[jj][ii-1] - 2.0*tmp[jj][ii  ] + tmp[jj][ii+1];
			ddp = tmp[jj][ii  ] - 2.0*tmp[jj][ii+1] + tmp[jj][ii+2];
			
			gdn = - tmp[jj][ii-1] + tmp[jj][ii  ] + 0.5*minmod(ddc,ddn);
			gdp = - tmp[jj][ii  ] + tmp[jj][ii+1] + 0.5*minmod(ddc,ddp);
			
			ggx = 0.0;
			if     (gdp*sgn < 0 && gdn*sgn < -gdp*sgn)	ggx = gdp;
			else if(gdn*sgn > 0 && gdp*sgn > -gdn*sgn)	ggx = gdn;
			else if(gdn*sgn < 0 && gdp*sgn > 0)			ggx = 0.0;
			
			
			ddn = tmp[jj-2][ii] - 2.0*tmp[jj-1][ii] + tmp[jj  ][ii];
			ddc = tmp[jj-1][ii] - 2.0*tmp[jj  ][ii] + tmp[jj+1][ii];
			ddp = tmp[jj  ][ii] - 2.0*tmp[jj+1][ii] + tmp[jj+2][ii];
			
			gdn = - tmp[jj-1][ii] + tmp[jj  ][ii] + 0.5*minmod(ddc,ddn);
			gdp = - tmp[jj  ][ii] + tmp[jj+1][ii] + 0.5*minmod(ddc,ddp);
			
			ggy = 0.0;
			if     (gdp*sgn < 0 && gdn*sgn < -gdp*sgn)	ggy = gdp;
			else if(gdn*sgn > 0 && gdp*sgn > -gdn*sgn)	ggy = gdn;
			else if(gdn*sgn < 0 && gdp*sgn > 0)			ggy = 0.0;
			
			lvs3[jj][ii] = tmp[jj][ii] - 0.25*sgn*(sqrt(pow(ggx,2.0) + pow(ggy,2.0))-1.0);
			
			
			if(fabs(phit[jj][ii]-0.5) < 0.3
			&&      phil[jj][ii]      < 0.01)
				lvs3[jj][ii] = sfd[jj][ii];
			
			if(fabs(lvs3[jj][ii]) < 20.0
			&&      phil[jj][ii]  > 0.01){
				eps = 1.0 - (lvs3[jj][ii]+1.0e-99)/(tmp[jj][ii]+1.0e-99);
				sum += eps*eps;
				cnt ++;
			}
			
		}
		}
		sum = sqrt(sum/(double)cnt);
		
		
		boundary_lvs(lvs3);
		memmove_double(tmp, lvs3);
		
		
		if(kk<50) sum=1;
		kk++;
		
	}
	// printf("kk=%d  cnt=%d  sum=%le\n", kk, cnt, sum);
	
}

