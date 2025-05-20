extern void boundary_bicg(bicg_fnc *);
extern void init_pp(double **);

void bicg_00(pois_fnc *pois, double **phir, double **phim, double **phit2, double **phit3)
{
	int ii, jj;
	int ex, ey;
	double hhxm, hhxp;
	double hhym, hhyp;
	double roxm, roxp;
	double roym, royp;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	private(hhxm, hhxp, hhym, hhyp) \
	private(roxm, roxp, roym, royp) \
	shared(para, ex, ey, pois, phir, phim, phit2, phit3)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		hhxm = +(phir[jj  ][ii-1] + phir[jj  ][ii  ])*0.5
			   +(phim[jj  ][ii-1] + phim[jj  ][ii  ])*0.5;
		hhxp = +(phir[jj  ][ii  ] + phir[jj  ][ii+1])*0.5
			   +(phim[jj  ][ii  ] + phim[jj  ][ii+1])*0.5;
		hhym = +(phir[jj-1][ii  ] + phir[jj  ][ii  ])*0.5
			   +(phim[jj-1][ii  ] + phim[jj  ][ii  ])*0.5;
		hhyp = +(phir[jj  ][ii  ] + phir[jj+1][ii  ])*0.5
			   +(phim[jj  ][ii  ] + phim[jj+1][ii  ])*0.5;
		
		hhxm = MAX(MIN(hhxm,1.0-1.0e-9),1.0e-9);
		hhxp = MAX(MIN(hhxp,1.0-1.0e-9),1.0e-9);
		hhym = MAX(MIN(hhym,1.0-1.0e-9),1.0e-9);
		hhyp = MAX(MIN(hhyp,1.0-1.0e-9),1.0e-9);
		
		roxm = (+phit2[jj  ][ii-1]+phit2[jj  ][ii  ]
				+phit3[jj  ][ii-1]+phit3[jj  ][ii  ])*0.25;
		roxp = (+phit2[jj  ][ii  ]+phit2[jj  ][ii+1]
				+phit3[jj  ][ii  ]+phit3[jj  ][ii+1])*0.25;
		roym = (+phit2[jj-1][ii  ]+phit2[jj  ][ii  ]
				+phit3[jj-1][ii  ]+phit3[jj  ][ii  ])*0.25;
		royp = (+phit2[jj  ][ii  ]+phit2[jj+1][ii  ]
				+phit3[jj  ][ii  ]+phit3[jj+1][ii  ])*0.25;
		
		roxm = 1.0 + roxm*(para.ro2/para.ro1-1.0);
		roxp = 1.0 + roxp*(para.ro2/para.ro1-1.0);
		roym = 1.0 + roym*(para.ro2/para.ro1-1.0);
		royp = 1.0 + royp*(para.ro2/para.ro1-1.0);
		
		pois->xm[jj][ii] = para.idx*para.idx*(1.0-hhxm)/roxm;
		pois->xp[jj][ii] = para.idx*para.idx*(1.0-hhxp)/roxp;
		pois->ym[jj][ii] = para.idy*para.idy*(1.0-hhym)/roym;
		pois->yp[jj][ii] = para.idy*para.idy*(1.0-hhyp)/royp;
		
		pois->cc[jj][ii] =+pois->xm[jj][ii]
						  +pois->xp[jj][ii]
						  +pois->ym[jj][ii]
						  +pois->yp[jj][ii];
						
		
		if(ii == 2){
			pois->cc[jj][ii] = pois->cc[jj][ii] - pois->xm[jj][ii];
			pois->xm[jj][ii] = 0.0;
		}
		if(ii == para.nx+1){			
			pois->cc[jj][ii] = pois->cc[jj][ii] - pois->xp[jj][ii];
			pois->xp[jj][ii] = 0.0;
		}
		if(jj == 2){
			pois->cc[jj][ii] = pois->cc[jj][ii] - pois->ym[jj][ii];
			pois->ym[jj][ii] = 0.0;
		}
		if(jj == para.ny+1){
			pois->cc[jj][ii] = pois->cc[jj][ii] - pois->yp[jj][ii];
			pois->yp[jj][ii] = 0.0;
		}
		
		pois->ic[jj][ii] = 1.0/(pois->cc[jj][ii]+1.0e-99);
		
	}
	}

}

void bicg_01(pois_fnc *pois, double ** dp, double **div, double *sum)
{
	int ii, jj;
	double lap_d, lap_p;
	double sum0=0.0, sum1=0.0, sum2=0.0, sum3=0.0, sum4=0.0, sum5=0.0;
	int ex, ey;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, lap_d, lap_p) \
	shared(para, ex, ey, dp, div, pois) \
	reduction(+:sum0, sum1, sum2, sum3, sum4, sum5)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		lap_d = +div[jj  ][ii-1]*pois->xm[jj][ii]
				+div[jj-1][ii  ]*pois->ym[jj][ii]
				-div[jj  ][ii  ]*pois->cc[jj][ii]
				+div[jj  ][ii+1]*pois->xp[jj][ii]
				+div[jj+1][ii  ]*pois->yp[jj][ii];
				
		lap_p = +dp[jj  ][ii-1]*pois->xm[jj][ii]
				+dp[jj-1][ii  ]*pois->ym[jj][ii]
				-dp[jj  ][ii  ]*pois->cc[jj][ii]
				+dp[jj  ][ii+1]*pois->xp[jj][ii]
				+dp[jj+1][ii  ]*pois->yp[jj][ii];
		
		sum0 += lap_d*lap_d;
		sum1 += lap_p*lap_p;
		sum2 += div[jj][ii]*div[jj][ii];
		sum3 += lap_d*lap_p;
		sum4 += lap_d*div[jj][ii];
		sum5 += lap_p*div[jj][ii];
	
	}
	}
	
	sum[0] = sum0;
	sum[1] = sum1;
	sum[2] = sum2;
	sum[3] = sum3;
	sum[4] = sum4;
	sum[5] = sum5;
	
}

void bicg_02(double **dpp, double **div, double cf1, double cf2)
{
	int ii, jj;
	int ex, ey;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex, ey, para, dpp, div, cf1, cf2)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		dpp[jj][ii] = + cf1*div[jj][ii]
					  + cf2*dpp[jj][ii];
					 

	}
	}
	
}

void bicg_03(pois_fnc *pois, bicg_fnc *bicg, double **dpp, double **div, double *sum)
{
	int ii, jj;
	double ddp;
	double tmp = 0.0;
	int ex, ey;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, ddp) \
	shared(ex, ey, pois, bicg, dpp, div) \
	reduction(+:tmp)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
	
		ddp = +dpp[jj  ][ii-1]*pois->xm[jj][ii]
			  +dpp[jj-1][ii  ]*pois->ym[jj][ii]
			  -dpp[jj  ][ii  ]*pois->cc[jj][ii]
			  +dpp[jj  ][ii+1]*pois->xp[jj][ii]
			  +dpp[jj+1][ii  ]*pois->yp[jj][ii]
			  -div[jj  ][ii  ];
		
		bicg->r0[jj][ii] = ddp*pois->ic[jj][ii];
		bicg->pp[jj][ii] = bicg->r0[jj][ii];
		bicg->rr[jj][ii] = bicg->r0[jj][ii];
		
		tmp += bicg->r0[jj][ii]*bicg->r0[jj][ii];
		
	}
	}
	
	*sum = tmp;
	
}

void bicg_04(pois_fnc *pois, bicg_fnc *bicg, double *sum)
{
	int ii, jj;
	double tmp = 0.0;
	int ex, ey;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex, ey, pois, bicg) \
	reduction(+:tmp)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
	
		bicg->ap[jj][ii] =(-bicg->pp[jj  ][ii-1]*pois->xm[jj][ii]
						   -bicg->pp[jj-1][ii  ]*pois->ym[jj][ii]
						   +bicg->pp[jj  ][ii  ]*pois->cc[jj][ii]
						   -bicg->pp[jj  ][ii+1]*pois->xp[jj][ii]
						   -bicg->pp[jj+1][ii  ]*pois->yp[jj][ii])
						   *pois->ic[jj  ][ii  ];
		
		tmp += bicg->r0[jj][ii]*bicg->ap[jj][ii];
		
	}
	}
	
	*sum = tmp;
	
}

void bicg_05(bicg_fnc *bicg, double aa)
{
	int ii, jj;
	int ex, ey;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex, ey, bicg, aa)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
	
		bicg->ee[jj][ii] = +bicg->rr[jj][ii]
						   -bicg->ap[jj][ii]*aa;
		
	}
	}
	
}

void bicg_06(pois_fnc *pois, bicg_fnc *bicg, double *sum)
{
	int ii, jj;
	double sum0=0.0, sum1=0.0;
	int ex, ey;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex, ey, pois, bicg) \
	reduction(+:sum0, sum1)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		bicg->ae[jj][ii] =(-bicg->ee[jj  ][ii-1]*pois->xm[jj][ii]
						   -bicg->ee[jj-1][ii  ]*pois->ym[jj][ii]
						   +bicg->ee[jj  ][ii  ]*pois->cc[jj][ii]
						   -bicg->ee[jj  ][ii+1]*pois->xp[jj][ii]
						   -bicg->ee[jj+1][ii  ]*pois->yp[jj][ii])
						                        *pois->ic[jj][ii];
		
		sum0 += bicg->ae[jj][ii]*bicg->ee[jj][ii];
		sum1 += bicg->ae[jj][ii]*bicg->ae[jj][ii];
		
	
	}
	}
	
	sum[0] = sum0;
	sum[1] = sum1;
	
}

void bicg_07(bicg_fnc *bicg, double **dpp, double aa, double c3, double *sum)
{
	int ii, jj;
	double tmp = 0.0;
	int ex, ey;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex, ey, dpp, bicg, aa, c3) \
	reduction(+:tmp)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		dpp[jj][ii] += +c3*bicg->ee[jj][ii]
					   +aa*bicg->pp[jj][ii];
		
		bicg->rr[jj][ii] =     bicg->ee[jj][ii]
						 - c3*bicg->ae[jj][ii];
		
		tmp += bicg->r0[jj][ii]*bicg->rr[jj][ii];
		
	}
	}
	
	*sum = tmp;	
	
}

void bicg_08(pois_fnc *pois, double **dpp, double **div, double *sum)
{
	int ii, jj;
	double ddp;
	double tmp = 0.0;
	int ex, ey;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, ddp) \
	shared(ex, ey, dpp, div, pois) \
	reduction(+:tmp)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		ddp = +dpp[jj  ][ii-1]*pois->xm[jj][ii]
			  +dpp[jj-1][ii  ]*pois->ym[jj][ii]
			  -dpp[jj  ][ii  ]*pois->cc[jj][ii]
			  +dpp[jj  ][ii+1]*pois->xp[jj][ii]
			  +dpp[jj+1][ii  ]*pois->yp[jj][ii]
			  -div[jj][ii];
		
		tmp += ddp*ddp;
	
	}
	}	
	
	*sum = sqrt(tmp/(double)(para.nx*para.ny));;	
	
}

void bicg_09(bicg_fnc *bicg, double bb, double c3)
{
	int ii, jj;
	int ex, ey;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex, ey, bicg, bb, c3)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		bicg->pp[jj][ii] = +bicg->rr[jj][ii]
						   +bb*(bicg->pp[jj][ii] - c3*bicg->ap[jj][ii]);
		
	}
	}
	
	
}


void calc_soldp_bicg(
	pois_fnc *pois, bicg_fnc *bicg, 
	double **dpp, double ** div, 
	double err, int ite, int *ans)
{
	int ii, jj;
	double sum;
	double err0, err1, err2, err3;
	
	double sma[6];
	double det, cf1, cf2;
	double c1, c2, c3;
	double alpha, beta;
	double d1, d2;
	
	
	err0 = err;
	err3 = err;
	bicg_01(pois, dpp, div, sma);
	

	det =   sma[0]*sma[1]-sma[3]*sma[3];
	cf1 = (+sma[1]*sma[4]-sma[3]*sma[5])*det/(det*det+1.0e-99);
	cf2 = (-sma[3]*sma[4]+sma[0]*sma[5])*det/(det*det+1.0e-99);

	bicg_02(dpp, div, cf1, cf2);	
	bicg_03(pois, bicg, dpp, div, &sum);

	c1 = sum;

	boundary_bicg(bicg);
	
	ii = 0;
	jj = 0;
	sum = 1.0;
	while (sum>para.Ep && ii < ite){
		
		bicg_04(pois, bicg, &sum);

		c2 = sum;
		alpha = c1*c2/(c2*c2+1.0e-99);
		
		bicg_05(bicg, alpha);
		boundary_bicg(bicg);
		
		bicg_06(pois, bicg, sma);		
		
		d1 = sma[0];
		d2 = sma[1];
		c3 = d1*d2/(d2*d2+1.0e-99);
		
		bicg_07(bicg, dpp, alpha, c3, &sum);
		
		c1 = sum;
		beta = c1*(c2*c3)/((c2*c3)*(c2*c3)+1.0e-99);
		
		boundary_pp(dpp);
		boundary_bicg(bicg);
		
		bicg_08(pois, dpp, div, &sum);
		
		err1 = sum;
		sum  = err1/err0;
		
		if(ii==0) err2 = err1 + 1.0e-99;
		if(ii>=1 && err1>=err3)	jj++;
		if(jj>=0 && err1/err2>1000){
			
			init_pp (dpp);
			boundary_pp(dpp);
			ii = 1e9;
		}
		else{
			
			err3 = err1;
			bicg_09(bicg, beta, c3);
			boundary_bicg(bicg);
			
		}
		
		if(ii<5) sum=1.0;
		ii++;
	}
	
	*ans = ii;
	
	
	
	
}


