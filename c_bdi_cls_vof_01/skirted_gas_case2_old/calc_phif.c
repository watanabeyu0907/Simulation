extern void calc_nrml(double **, double **, double **, double **, double **, coord *, coord *);
	
void mthinc_00(double **phir, double **phim, double **phif, double **phil)
{
	int ii, jj;
	int ex4, ey4;
	
	
	ex4 = para.nx+4;
	ey4 = para.ny+4;
	

	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex4, ey4, phir, phim, phif, phil)
	#endif
	for(jj=0; jj<ey4; jj++){
	for(ii=0; ii<ex4; ii++){
		
		phil[jj][ii] = MAX(MIN(phir[jj][ii] + phim[jj][ii] + phif[jj][ii], 1.0), 0.0);
		
	}
	}
	
}

void mthinc_01(double **phi, coord *rnf, double **sfd)
{
	int ii, jj, ll;
	int ex, ey;
	double coef, rnx, rny, q00;
	double x11, x12, x21, x22;
	double veldt_dl, dd, ff, df;
	
		
	ex = para.nx+3;
	ey = para.ny+3;
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, ll) \
	private(coef, rnx, rny, q00) \
	private(x11, x12, x21, x22) \
	private(veldt_dl, dd, ff, df)\
	shared(ex, ey, para) \
	shared(phi, rnf, sfd)
	#endif
	for(jj=1; jj<ey; jj++){
	for(ii=1; ii<ex; ii++){
		
		rnx = ( +rnf->xx[jj-1][ii-1]
				+rnf->xx[jj-1][ii  ]
				+rnf->xx[jj  ][ii-1]
				+rnf->xx[jj  ][ii  ])*0.25;
		rny = ( +rnf->yy[jj-1][ii-1]
				+rnf->yy[jj-1][ii  ]
				+rnf->yy[jj  ][ii-1]
				+rnf->yy[jj  ][ii  ])*0.25;
		q00 = sqrt(pow(rnx,2.0)+pow(rny,2.0));
		
		if(q00 <= para.eps_mtc*para.idc){
			
			if(      phi[jj][ii] <= para.eps_trn) sfd[jj][ii] = -3.0;
			if(1.0 - phi[jj][ii] <= para.eps_trn) sfd[jj][ii] =  3.0;
			
		}
		else{
			
			coef = 1.0/(sqrt(pow(rnx,2.0) + pow(rny,2.0)+1.0e-99));
			rnx *= coef;
			rny *= coef;
			
			x11 = rny*para.gp1 + rnx*para.gp1;
			x12 = rny*para.gp1 + rnx*para.gp2;
			x21 = rny*para.gp2 + rnx*para.gp1;
			x22 = rny*para.gp2 + rnx*para.gp2;
			
			dd = 0.0;
		
			for(ll=0; ll<5; ll++){
				ff = -phi[jj][ii]
					+0.25/(1.0+exp(-2.0*para.bet*(x11+dd)))
					+0.25/(1.0+exp(-2.0*para.bet*(x12+dd)))
					+0.25/(1.0+exp(-2.0*para.bet*(x21+dd)))
					+0.25/(1.0+exp(-2.0*para.bet*(x22+dd)));
					
				df = 0.5*para.bet*(
					+0.25/pow(cosh(para.bet*(x11+dd)),2.0)
					+0.25/pow(cosh(para.bet*(x12+dd)),2.0)
					+0.25/pow(cosh(para.bet*(x21+dd)),2.0)
					+0.25/pow(cosh(para.bet*(x22+dd)),2.0));
					
				dd = dd-ff*df/(df*df+1.0e-99);
			}
			
			sfd[jj][ii] = dd;
			
		}
		
	}
	}
	
}

void mthinc_02(coord *uu, double **phi, coord *rnf, coord *flx, double **sfd)
{
	int ii, jj;
	int ex, ey;
	double coef, rnx, rny, q00;
	double x11, x12, x21, x22;
	double veldt_dl, dd;
	double xgp1, xgp2, ygp1, ygp2;
	double xst, xen, yst, yen;
	
	
	init_uu(flx);
	
		
	ex = para.nx+3;
	ey = para.ny+3;
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	private(coef, rnx, rny, q00) \
	private(x11, x12, x21, x22) \
	private(veldt_dl, dd)\
	private(xgp1, xgp2, ygp1, ygp2)\
	private(xst, xen, yst, yen)\
	shared(ex, ey, para) \
	shared(uu, phi, rnf, flx, sfd)
	#endif
	for(jj=1; jj<ey; jj++){
	for(ii=1; ii<ex; ii++){
		
		rnx = ( +rnf->xx[jj-1][ii-1]
				+rnf->xx[jj-1][ii  ]
				+rnf->xx[jj  ][ii-1]
				+rnf->xx[jj  ][ii  ])*0.25;
		rny = ( +rnf->yy[jj-1][ii-1]
				+rnf->yy[jj-1][ii  ]
				+rnf->yy[jj  ][ii-1]
				+rnf->yy[jj  ][ii  ])*0.25;
		q00 = sqrt(pow(rnx,2.0)+pow(rny,2.0));
		
		if(q00 <= para.eps_mtc*para.idc){
			
			if(      phi[jj][ii] <= para.eps_trn) phi[jj][ii] =  0.0;
			if(1.0 - phi[jj][ii] <= para.eps_trn) phi[jj][ii] =  1.0;
						
			if(uu->xx[jj  ][ii-1] < 0.0)
				flx->xx[jj  ][ii-1] = uu->xx[jj  ][ii-1]*phi[jj][ii]*para.dT*para.idx;
			if(uu->xx[jj  ][ii  ] > 0.0)
				flx->xx[jj  ][ii  ] = uu->xx[jj  ][ii  ]*phi[jj][ii]*para.dT*para.idx;
			if(uu->yy[jj-1][ii  ] < 0.0)
				flx->yy[jj-1][ii  ] = uu->yy[jj-1][ii  ]*phi[jj][ii]*para.dT*para.idy;
			if(uu->yy[jj  ][ii  ] > 0.0)
				flx->yy[jj  ][ii  ] = uu->yy[jj  ][ii  ]*phi[jj][ii]*para.dT*para.idy;
			
		}
		else{
			
			coef = 1.0/(sqrt(pow(rnx,2.0) + pow(rny,2.0)+1.0e-99));
			rnx *= coef;
			rny *= coef;
			
			x11 = rny*para.gp1 + rnx*para.gp1;
			x12 = rny*para.gp1 + rnx*para.gp2;
			x21 = rny*para.gp2 + rnx*para.gp1;
			x22 = rny*para.gp2 + rnx*para.gp2;
			
			dd = sfd[jj][ii];
			
			if(uu->xx[jj  ][ii-1] < 0.0){
				veldt_dl = uu->xx[jj][ii-1]*para.dT*para.idx;
				xst = -0.5;
				xen = -0.5-veldt_dl;
				xgp1= (xst+xen)*0.5+(-xst+xen)*para.gp1;
				xgp2= (xst+xen)*0.5+(-xst+xen)*para.gp2;
				ygp1= para.gp1;
				ygp2= para.gp2;
				x11 = rny*ygp1 + rnx*xgp1;
				x12 = rny*ygp1 + rnx*xgp2;
				x21 = rny*ygp2 + rnx*xgp1;
				x22 = rny*ygp2 + rnx*xgp2;
				flx->xx[jj][ii-1] = 0.25*veldt_dl*(
					+1.0/(1.0+exp(-2.0*para.bet*(x11+dd)))
					+1.0/(1.0+exp(-2.0*para.bet*(x12+dd)))
					+1.0/(1.0+exp(-2.0*para.bet*(x21+dd)))
					+1.0/(1.0+exp(-2.0*para.bet*(x22+dd))));				
			}
			
			if(uu->xx[jj  ][ii  ] > 0.0){
				veldt_dl = uu->xx[jj][ii]*para.dT*para.idx;
				xst = +0.5-veldt_dl;
				xen = +0.5;
				xgp1= (xst+xen)*0.5+(-xst+xen)*para.gp1;
				xgp2= (xst+xen)*0.5+(-xst+xen)*para.gp2;
				ygp1= para.gp1;
				ygp2= para.gp2;
				x11 = rny*ygp1 + rnx*xgp1;
				x12 = rny*ygp1 + rnx*xgp2;
				x21 = rny*ygp2 + rnx*xgp1;
				x22 = rny*ygp2 + rnx*xgp2;
				flx->xx[jj][ii] = 0.25*veldt_dl*(
					+1.0/(1.0+exp(-2.0*para.bet*(x11+dd)))
					+1.0/(1.0+exp(-2.0*para.bet*(x12+dd)))
					+1.0/(1.0+exp(-2.0*para.bet*(x21+dd)))
					+1.0/(1.0+exp(-2.0*para.bet*(x22+dd))));				
			}
			
			if(uu->yy[jj-1][ii  ] < 0.0){
				veldt_dl = uu->yy[jj-1][ii]*para.dT*para.idx;
				yst = -0.5;
				yen = -0.5-veldt_dl;
				xgp1= para.gp1;
				xgp2= para.gp2;
				ygp1= (yst+yen)*0.5+(-yst+yen)*para.gp1;
				ygp2= (yst+yen)*0.5+(-yst+yen)*para.gp2;
				x11 = rny*ygp1 + rnx*xgp1;
				x12 = rny*ygp1 + rnx*xgp2;
				x21 = rny*ygp2 + rnx*xgp1;
				x22 = rny*ygp2 + rnx*xgp2;
				flx->yy[jj-1][ii] = 0.25*veldt_dl*(
					+1.0/(1.0+exp(-2.0*para.bet*(x11+dd)))
					+1.0/(1.0+exp(-2.0*para.bet*(x12+dd)))
					+1.0/(1.0+exp(-2.0*para.bet*(x21+dd)))
					+1.0/(1.0+exp(-2.0*para.bet*(x22+dd))));				
			}
			
			if(uu->yy[jj  ][ii  ] > 0.0){
				veldt_dl = uu->yy[jj][ii]*para.dT*para.idx;
				yst = +0.5-veldt_dl;
				yen = +0.5;
				xgp1= para.gp1;
				xgp2= para.gp2;
				ygp1= (yst+yen)*0.5+(-yst+yen)*para.gp1;
				ygp2= (yst+yen)*0.5+(-yst+yen)*para.gp2;
				x11 = rny*ygp1 + rnx*xgp1;
				x12 = rny*ygp1 + rnx*xgp2;
				x21 = rny*ygp2 + rnx*xgp1;
				x22 = rny*ygp2 + rnx*xgp2;
				flx->yy[jj][ii] = 0.25*veldt_dl*(
					+1.0/(1.0+exp(-2.0*para.bet*(x11+dd)))
					+1.0/(1.0+exp(-2.0*para.bet*(x12+dd)))
					+1.0/(1.0+exp(-2.0*para.bet*(x21+dd)))
					+1.0/(1.0+exp(-2.0*para.bet*(x22+dd))));				
			}
			
		}
		
	}
	}
	
}

void calc_phif(coord *uu,
	double **phir, double **phim, double **phif, double **phil, double **phit,
	coord *rnf, coord *rns, coord *flx, double **sfd1, double **sfd2, double **sfd3
)
{
	int ii, jj;
	int ex2, ey2;
	int ex4, ey4;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double sum3 = 0.0;
	
	
	ex2 = para.nx+2;
	ey2 = para.ny+2;
	ex4 = para.nx+4;
	ey4 = para.ny+4;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex2, ey2, phif) \
	reduction(+:sum1)
	#endif
	for(jj=2; jj<ey2; jj++){
	for(ii=2; ii<ex2; ii++){
		
		sum1 += phif[jj][ii];
		
	}
	}
	
	
	mthinc_00(phir, phim, phif, phil);
	calc_nrml(phir, phim, phif, phil, phit, rnf, rns);
	mthinc_01(phif, rnf, sfd1);
	mthinc_02(uu, phif, rnf, flx, sfd1);
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex2, ey2, para) \
	shared(phir, phim, phif, phit, flx, uu)
	#endif
	for(jj=2; jj<ey2; jj++){
	for(ii=2; ii<ex2; ii++){
	
		if(phit[jj][ii] < 0.5){
			phif[jj][ii] = +phif[jj][ii]
						  +(+flx->xx[jj  ][ii-1]-flx->xx[jj][ii]
							+flx->yy[jj-1][ii  ]-flx->yy[jj][ii])
							+para.dT*phif[jj][ii]*
							(+(-uu->xx[jj  ][ii-1]+uu->xx[jj][ii])*para.idx
							 +(-uu->yy[jj-1][ii  ]+uu->yy[jj][ii])*para.idy);
			phif[jj][ii] = MIN(MAX(MIN(phif[jj][ii], 1.0), 0.0),1.0-phir[jj][ii]-phim[jj][ii]);
		}
	}
	}
	
	
	mthinc_01(phil, rnf, sfd2);
	mthinc_02(uu, phil, rnf, flx, sfd2);
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex2, ey2, para) \
	shared(phir, phim, phif, phil, phit, flx, uu)
	#endif
	for(jj=2; jj<ey2; jj++){
	for(ii=2; ii<ex2; ii++){
	
		if(phit[jj][ii] >= 0.5){
			phif[jj][ii] = -phir[jj][ii]-phim[jj][ii]+phil[jj][ii]
						  +(+flx->xx[jj  ][ii-1]-flx->xx[jj][ii]
							+flx->yy[jj-1][ii  ]-flx->yy[jj][ii])
							+para.dT*phil[jj][ii]*
							(+(-uu->xx[jj  ][ii-1]+uu->xx[jj][ii])*para.idx
							 +(-uu->yy[jj-1][ii  ]+uu->yy[jj][ii])*para.idy);
			phif[jj][ii] = MIN(MAX(MIN(phif[jj][ii], 1.0), 0.0),1.0-phir[jj][ii]-phim[jj][ii]);
		}
	}
	}
	
	
	boundary_pp(phif);

	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex2, ey2, phir, phim, phif) \
	reduction(+:sum2, sum3)
	#endif
	for(jj=2; jj<ey2; jj++){
	for(ii=2; ii<ex2; ii++){
		
		sum2 += phif[jj][ii];
		sum3 += pow(phif[jj][ii],2.0)*pow(1.0-phif[jj][ii]-phir[jj][ii]-phim[jj][ii],2.0);
		
	}
	}
	
	sum1 = (sum1-sum2)/(sum3+1e-19);
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex2, ey2, phir, phim, phif, sum1)
	#endif
	for(jj=2; jj<ey2; jj++){
	for(ii=2; ii<ex2; ii++){
		
		phif[jj][ii]+= sum1*pow(phif[jj][ii],2.0)*pow(1.0-phif[jj][ii]-phir[jj][ii]-phim[jj][ii],2.0);
		phif[jj][ii] = MAX(MIN(phif[jj][ii], 1.0), 0.0);
		
	}
	}
	
	boundary_pp(phif);
	
		
	mthinc_00(phir, phim, phif, phil);
	calc_nrml(phir, phim, phif, phil, phit, rnf, rns);
	mthinc_01(phif, rnf, sfd1);
	mthinc_01(phil, rnf, sfd2);
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex4, ey4, sfd1, sfd2, sfd3, phit, phir, phim) 
	#endif
	for(jj=0; jj<ey4; jj++){
	for(ii=0; ii<ex4; ii++){
		
		sfd3[jj][ii] = 0.5*(sfd1[jj][ii] + sfd2[jj][ii]);
		
	}
	}
	
}

