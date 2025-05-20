void calc_NSeq0_xx(visc_fnc *dd, double **phir, double **phim, double **phit2, double **phit3)
{
	int ii, jj;
	int ex, ey;
	double rnxm, rnxp, rnym, rnyp;
	double nu1, nu2, ph, pw, ro, cf;
	
	
	ex = para.nx+1;
	ey = para.ny+2;
	
	
	// ( 1-a ) x vel
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	private(rnxm, rnxp, rnym, rnyp) \
	private(nu1, nu2, ph, pw, ro, cf) \
	shared(para, ex, ey, phir, phim, phit2, phit3, dd)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		nu1 = para.mu1/para.ro1;
		nu2 = para.mu2/para.ro2;
		
		rnxm = +phit3[jj  ][ii  ];
		rnxp = +phit3[jj  ][ii+1];
		rnym =(+phit3[jj-1][ii  ]
			   +phit3[jj-1][ii+1]
			   +phit3[jj  ][ii  ]
			   +phit3[jj  ][ii+1])*0.25;
		rnyp =(+phit3[jj  ][ii  ]
			   +phit3[jj  ][ii+1]
			   +phit3[jj+1][ii  ]
			   +phit3[jj+1][ii+1])*0.25;
			   
		// rnxm = nu1 / (1.0 + rnxm*(nu1/nu2-1.0));
		// rnxp = nu1 / (1.0 + rnxp*(nu1/nu2-1.0));
		// rnym = nu1 / (1.0 + rnym*(nu1/nu2-1.0));
		// rnyp = nu1 / (1.0 + rnyp*(nu1/nu2-1.0));
		
		rnxm = nu1 / (1.0 + rnxm*(para.mu1/para.mu2-1.0));
		rnxp = nu1 / (1.0 + rnxp*(para.mu1/para.mu2-1.0));
		rnym = nu1 / (1.0 + rnym*(para.mu1/para.mu2-1.0));
		rnyp = nu1 / (1.0 + rnyp*(para.mu1/para.mu2-1.0));
		
		ph = (+phit2[jj  ][ii  ] + phit2[jj  ][ii+1]
			  +phit3[jj  ][ii  ] + phit3[jj  ][ii+1])*0.25;
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
			  
		ro = 1.0 + ph*(para.ro2/para.ro1-1.0);
		
		pw = +(+phir[jj  ][ii  ] + phir[jj  ][ii+1])*0.5
			 +(+phim[jj  ][ii  ] + phim[jj  ][ii+1])*0.5;
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		
		cf = (1.0-pw)/ro;
		
		dd->xmxm[jj][ii] = +cf*para.idx*para.idx*rnxm*2.0;
		dd->xpxp[jj][ii] = +cf*para.idx*para.idx*rnxp*2.0;
		dd->ymym[jj][ii] = +cf*para.idy*para.idy*rnym;
		dd->ypyp[jj][ii] = +cf*para.idy*para.idy*rnyp;
		dd->ymxm[jj][ii] = +cf*para.idx*para.idy*rnym;
		dd->ymxp[jj][ii] = -cf*para.idx*para.idy*rnym;
		dd->ypxm[jj][ii] = -cf*para.idx*para.idy*rnyp;
		dd->ypxp[jj][ii] = +cf*para.idx*para.idy*rnyp;
		dd->ccxx[jj][ii] = -cf*para.idx*para.idx*(rnxm+rnxp)*2.0;
		dd->ccyy[jj][ii] = -cf*para.idy*para.idy*(rnym+rnyp);
		
	}
	}
	
}

