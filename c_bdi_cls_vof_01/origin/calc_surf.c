void calc_surf(double **phir, double **phim, double **lvs0, double **lvm3, double **lvf3, coord *rnf, coord *rns, coord *rnfc, coord *rnfa, coord *rnsc, coord *rnsa, double **phit3, coord *sft)
{
	int ii, jj;
	int ex1, ey1;
	int ex2, ey2;
	int ex3, ey3;
	int ex4, ey4;
	double q00, q01, q10, q11, qqq;
	double cs, sn, kp;
	double rnx, rny, pnx, pny, tmp;
	double qf, qs;
	double p11, p12, p21, p22;
	
	
	sn  = sin(para.ca*pi/180.0);
	cs  = cos(para.ca*pi/180.0);
	
	
	ex1 = para.nx+1;
	ey1 = para.ny+1;
	ex2 = para.nx+2;
	ey2 = para.ny+2;
	ex3 = para.nx+3;
	ey3 = para.ny+3;
	ex4 = para.nx+4;
	ey4 = para.ny+4;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, tmp) \
	private(q00, q01, q10, q11, qqq) \
	shared(ex3, ey3, para) \
	shared(para, lvs0, lvm3, lvf3, rnf, rns)
	#endif
	for(jj=0; jj<ey3; jj++){
	for(ii=0; ii<ex3; ii++){
		
		q00 = lvf3[jj  ][ii  ];
		q01 = lvf3[jj  ][ii+1];
		q10 = lvf3[jj+1][ii  ];
		q11 = lvf3[jj+1][ii+1];
		qqq = 0.25*(q00 + q01 + q10 + q11);
		
		rnf->xx[jj][ii] = (-q00+q01-q10+q11)*para.idx*0.5;
		rnf->yy[jj][ii] = (-q00-q01+q10+q11)*para.idy*0.5;
	
		q00 = lvm3[jj  ][ii  ];
		q01 = lvm3[jj  ][ii+1];
		q10 = lvm3[jj+1][ii  ];
		q11 = lvm3[jj+1][ii+1];
		qqq = 0.25*(q00 + q01 + q10 + q11);
		
		if(fabs(qqq) < 2.5){
			rns->xx[jj][ii] = (-q00+q01-q10+q11)*para.idx*0.5;
			rns->yy[jj][ii] = (-q00-q01+q10+q11)*para.idy*0.5;
		}
	
		q00 = lvs0[jj  ][ii  ];
		q01 = lvs0[jj  ][ii+1];
		q10 = lvs0[jj+1][ii  ];
		q11 = lvs0[jj+1][ii+1];
		qqq = 0.25*(q00 + q01 + q10 + q11);
		
		if(fabs(qqq) < 2.5){
			rns->xx[jj][ii] = (-q00+q01-q10+q11)*para.idx*0.5;
			rns->yy[jj][ii] = (-q00-q01+q10+q11)*para.idy*0.5;
		}
		
	}
	}
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	private(rnx, rny, tmp) \
	shared(ex4, ey4, rnf, rns, rnfc, rnsc) 
	#endif
	for(jj=1; jj<ey4; jj++){
	for(ii=1; ii<ex4; ii++){
		
		rnx = +(+rnf->xx[jj-1][ii-1]
				+rnf->xx[jj-1][ii  ]
				+rnf->xx[jj  ][ii-1]
				+rnf->xx[jj  ][ii  ])*0.25;
		rny = +(+rnf->yy[jj-1][ii-1]
				+rnf->yy[jj-1][ii  ]
				+rnf->yy[jj  ][ii-1]
				+rnf->yy[jj  ][ii  ])*0.25;
		tmp = 1.0/sqrt(pow(rnx,2.0)+pow(rny,2.0)+1.0e-99);
		rnfc->xx[jj][ii] = rnx*tmp;
		rnfc->yy[jj][ii] = rny*tmp;
		
		
		rnx = +(+rns->xx[jj-1][ii-1]
				+rns->xx[jj-1][ii  ]
				+rns->xx[jj  ][ii-1]
				+rns->xx[jj  ][ii  ])*0.25;
		rny = +(+rns->yy[jj-1][ii-1]
				+rns->yy[jj-1][ii  ]
				+rns->yy[jj  ][ii-1]
				+rns->yy[jj  ][ii  ])*0.25;
		tmp = 1.0/sqrt(pow(rnx,2.0)+pow(rny,2.0)+1.0e-99);
		rnsc->xx[jj][ii] = rnx*tmp;
		rnsc->yy[jj][ii] = rny*tmp;
		
	}
	}
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	private(rnx, rny, tmp) \
	shared(ex4, ey4, rnf, rns, rnfa, rnsa) 
	#endif
	for(jj=0; jj<ey4; jj++){
	for(ii=0; ii<ex4; ii++){
		
		rnx = rnf->xx[jj][ii];
		rny = rnf->yy[jj][ii];
		tmp = 1.0/sqrt(pow(rnx,2.0)+pow(rny,2.0)+1.0e-99);
		rnfa->xx[jj][ii] = rnx*tmp;
		rnfa->yy[jj][ii] = rny*tmp;
		
		rnx = rns->xx[jj][ii];
		rny = rns->yy[jj][ii];
		tmp = 1.0/sqrt(pow(rnx,2.0)+pow(rny,2.0)+1.0e-99);
		rnsa->xx[jj][ii] = rnx*tmp;
		rnsa->yy[jj][ii] = rny*tmp;
		
	}
	}
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	private(rnx, rny, pnx, pny, tmp, qf, qs) \
	shared(p11, p12, p21, p22, sn, cs)  \
	shared(para, ex4, ey4, rnf, rns, rnfc, rnsc, lvs0, lvf3, lvm3) 
	#endif
	for(jj=1; jj<ey4; jj++){
	for(ii=1; ii<ex4; ii++){
		
		qf = fabs(lvf3[jj][ii]);
		qs = MIN(fabs(lvs0[jj][ii]), fabs(lvm3[jj][ii]));
		
		if(qf < 2.5
		&& qs < 1.5)
		{
			p11 = 1.0   - rnsc->xx[jj][ii]*rnsc->xx[jj][ii];
			p22 = 1.0   - rnsc->yy[jj][ii]*rnsc->yy[jj][ii];
			p12 = p21 = - rnsc->xx[jj][ii]*rnsc->yy[jj][ii];
			rnx = rnfc->xx[jj][ii];
			rny = rnfc->yy[jj][ii];
			pnx = p11*rnx + p12*rny;
			pny = p21*rnx + p22*rny;
			tmp = 1.0/sqrt(pow(pnx,2.0)+pow(pny,2.0)+1.0e-99);
			rnfc->xx[jj][ii] = -cs*rnsc->xx[jj][ii] + fabs(sn)*tmp*pnx;
			rnfc->yy[jj][ii] = -cs*rnsc->yy[jj][ii] + fabs(sn)*tmp*pny;
		}
		
	}
	}
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	private(rnx, rny, pnx, pny, tmp, qf, qs) \
	shared(p11, p12, p21, p22, sn, cs)  \
	shared(para, ex3, ey3, rnf, rns, rnfa, rnsa, lvs0, lvm3, lvf3) 
	#endif
	for(jj=0; jj<ey3; jj++){
	for(ii=0; ii<ex3; ii++){
		
		qf =     fabs(0.25*(+lvf3[jj  ][ii  ]
							+lvf3[jj  ][ii+1]
							+lvf3[jj+1][ii  ]
							+lvf3[jj+1][ii+1]));
		qs = MIN(fabs((lvs0[jj  ][ii  ]
					  +lvs0[jj  ][ii+1]
					  +lvs0[jj+1][ii  ]
					  +lvs0[jj+1][ii+1])*0.25),
				 fabs((lvm3[jj  ][ii  ]
					  +lvm3[jj  ][ii+1]
					  +lvm3[jj+1][ii  ]
					  +lvm3[jj+1][ii+1])*0.25));
		
		if(qf < 2.5
		&& qs < 1.5)
		{
			p11 = 1.0   - rnsa->xx[jj][ii]*rnsa->xx[jj][ii];
			p22 = 1.0   - rnsa->yy[jj][ii]*rnsa->yy[jj][ii];
			p12 = p21 = - rnsa->xx[jj][ii]*rnsa->yy[jj][ii];
			rnx = rnfa->xx[jj][ii];
			rny = rnfa->yy[jj][ii];
			pnx = p11*rnx + p12*rny;
			pny = p21*rnx + p22*rny;
			tmp = 1.0/sqrt(pow(pnx,2.0)+pow(pny,2.0)+1.0e-99);
			rnfa->xx[jj][ii] = -cs*rnsa->xx[jj][ii] + fabs(sn)*tmp*pnx;
			rnfa->yy[jj][ii] = -cs*rnsa->yy[jj][ii] + fabs(sn)*tmp*pny;
		}
		
	}
	}
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	private(kp) \
	shared(ex1, ey2, para, rnfc, rnfa, phit3, sft, phir, phim, lvf3) 
	#endif
	for(jj=2; jj<ey2; jj++){
	for(ii=2; ii<ex1; ii++){
		kp = +(-rnfc->xx[jj  ][ii  ]+rnfc->xx[jj  ][ii+1])*para.idx
			 +(-rnfa->yy[jj-1][ii  ]+rnfa->yy[jj  ][ii  ])*para.idy;
		sft->xx[jj][ii] = -0.5*kp*(rnfc->xx[jj  ][ii  ]+rnfc->xx[jj  ][ii+1])*para.si/para.ro1*para.idx;
		sft->xx[jj][ii]*= para.bet/(2.0*pow(cosh(para.bet*(lvf3[jj][ii]+lvf3[jj][ii+1])),2.0));
		
		if(fabs(0.5*(lvf3[jj][ii]+lvf3[jj][ii+1])) > 3.5 ) sft->xx[jj][ii] = 0.0;
		if(fabs(0.5*(phir[jj][ii]+phir[jj][ii+1])
			   +0.5*(phim[jj][ii]+phim[jj][ii+1])) > 0.99) sft->xx[jj][ii] = 0.0;
	}
	}
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	private(kp) \
	shared(ex2, ey1, para, rnfc, rnfa, phit3, sft, phir, phim, lvf3) 
	#endif
	for(jj=2; jj<ey1; jj++){
	for(ii=2; ii<ex2; ii++){
		kp = +(-rnfa->xx[jj  ][ii-1]+rnfa->xx[jj  ][ii  ])*para.idx
			 +(-rnfc->yy[jj  ][ii  ]+rnfc->yy[jj+1][ii  ])*para.idy;
		sft->yy[jj][ii] = -0.5*kp*(rnfc->yy[jj  ][ii  ]+rnfc->yy[jj+1][ii  ])*para.si/para.ro1*para.idy;
		sft->yy[jj][ii]*= para.bet/(2.0*pow(cosh(para.bet*(lvf3[jj][ii]+lvf3[jj+1][ii])),2.0));
		
		if(fabs(0.5*(lvf3[jj][ii]+lvf3[jj+1][ii])) > 3.5 ) sft->yy[jj][ii] = 0.0;
		if(fabs(0.5*(phir[jj][ii]+phir[jj+1][ii])
			   +0.5*(phim[jj][ii]+phim[jj+1][ii])) > 0.99) sft->yy[jj][ii] = 0.0;
	}
	}
	
}

