# include "watasys.h"
# include "watafnc.h"

double smooth_sgn(double val)
{
	return (val)/sqrt(pow(val,2.0)+1.0);
}

double minmod(double aa, double bb)
{
	double ans;
	double sgn;
	
	if(aa<0)		sgn = -1.0;
	else if(aa>0)	sgn =  1.0;
	else			sgn =  0.0;
	
	if(aa*bb > 0){
		if (fabs(aa) <= fabs(bb))	ans = sgn*aa;
		else						ans = sgn*bb;
	}
	else{
		ans = 0.0;
	}
	
	return ans;
}

void init_lvss(double **phi, double **tmp, double **lvs)
{
	int ii, jj, kk;
	int cnt;
	int ex2, ey2;
	int ex4, ey4;
	double eps, sum;
	double sgn;
	double gdn,      gdp;
	double ddn, ddc, ddp;
	double ggx, ggy;
	
	
	ex2 = para.nx+2;
	ey2 = para.ny+2;
	ex4 = para.nx+4;
	ey4 = para.ny+4;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex4, ey4, phi, tmp)
	#endif
	for(jj=0; jj<ey4; jj++){
	for(ii=0; ii<ex4; ii++){
		
		tmp[jj][ii] = smooth_sgn(phi[jj][ii]-0.5);
		
	}
	}
	
	
	kk = 0;
	sum= 1;
	// narrow band levelset in 2nd order upwind ENO scheme
	while (sum>para.Ep && kk < 300){
		
		sum = 0.0;
		cnt = 0;
		
		#ifdef _OPENMP
		#pragma omp parallel for \
		default(none) \
		private(ii, jj, eps, sgn, ddn, ddc, ddp, gdn, gdp, ggx, ggy) \
		shared(kk, ex2, ey2, ex4, ey4, phi, tmp, lvs) \
		reduction(+:sum, cnt)
		#endif
		for(jj=2; jj<ey2; jj++){
		for(ii=2; ii<ex2; ii++){	
			
			sgn = smooth_sgn(phi[jj][ii]-0.5);
			
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
			
			lvs[jj][ii] = tmp[jj][ii] - 0.25*sgn*(sqrt(pow(ggx,2.0) + pow(ggy,2.0))-1.0);
			
			if(pow(sgn,2.0) < pow(0.75*0.5/sqrt(1.25),2.0))	lvs[jj][ii] = 2.0*phi[jj][ii]-1.0;
			
			if(fabs(lvs[jj][ii]) < 20.0){
				eps = 1.0 - (lvs[jj][ii]+1.0e-99)/(tmp[jj][ii]+1.0e-99);
				sum += eps*eps;
				cnt ++;
			}
			
		}
		}
		sum = sqrt(sum/(double)cnt);
		
		
		boundary_lvs(lvs);		
		memmove_double(tmp, lvs);
		
		if(kk<5) sum=1;
		kk++;
		
	}
	// printf("kk=%d  cnt=%d  sum=%le\n", kk, cnt, sum);
	
}

