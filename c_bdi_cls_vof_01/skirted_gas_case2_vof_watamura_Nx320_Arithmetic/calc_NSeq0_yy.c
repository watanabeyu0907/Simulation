void calc_NSeq0_yy(visc_fnc *dd, double **phir, double **phim, double **phit2, double **phit3)
{
	int ii, jj;
	int ex, ey;
	double rnxm, rnxp, rnym, rnyp;
	double nu1, nu2, ph, pw, ro, cf;
	
	
	ex = para.nx+2;
	ey = para.ny+1;
	
	
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
		
		rnxm =(+phit3[jj  ][ii-1]
			   +phit3[jj  ][ii  ]
			   +phit3[jj+1][ii-1]
			   +phit3[jj+1][ii  ])*0.25;
		rnxp =(+phit3[jj  ][ii  ]
			   +phit3[jj  ][ii+1]
			   +phit3[jj+1][ii  ]
			   +phit3[jj+1][ii+1])*0.25;
		rnym = +phit3[jj  ][ii  ];
		rnyp = +phit3[jj+1][ii  ];

		// rnxm = nu1 * (1.0 + rnxm*(nu2/nu1-1.0));
		// rnxp = nu1 * (1.0 + rnxp*(nu2/nu1-1.0));
		// rnym = nu1 * (1.0 + rnym*(nu2/nu1-1.0));
		// rnyp = nu1 * (1.0 + rnyp*(nu2/nu1-1.0));
			   

		rnxm = nu1 * (1.0 + rnxm*(para.mu2/para.mu1-1.0));
		rnxp = nu1 * (1.0 + rnxp*(para.mu2/para.mu1-1.0));
		rnym = nu1 * (1.0 + rnym*(para.mu2/para.mu1-1.0));
		rnyp = nu1 * (1.0 + rnyp*(para.mu2/para.mu1-1.0));

		ph = (+phit2[jj  ][ii  ] + phit2[jj+1][ii  ]
			  +phit3[jj  ][ii  ] + phit3[jj+1][ii  ])*0.25;
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
			  
		ro = 1.0 + ph*(para.ro2/para.ro1-1.0);
		
		pw = +(+phir[jj  ][ii  ] + phir[jj+1][ii  ])*0.5
			 +(+phim[jj  ][ii  ] + phim[jj+1][ii  ])*0.5;
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		
		cf = (1.0-pw)/ro;
		
		dd->xmxm[jj][ii] = +cf*para.idx*para.idx*rnxm;
		dd->xpxp[jj][ii] = +cf*para.idx*para.idx*rnxp;
		dd->ymym[jj][ii] = +cf*para.idy*para.idy*rnym*2.0;
		dd->ypyp[jj][ii] = +cf*para.idy*para.idy*rnyp*2.0;
		dd->ymxm[jj][ii] = +cf*para.idx*para.idy*rnxm;
		dd->ymxp[jj][ii] = -cf*para.idx*para.idy*rnxm;
		dd->ypxm[jj][ii] = -cf*para.idx*para.idy*rnxp;
		dd->ypxp[jj][ii] = +cf*para.idx*para.idy*rnxp;
		dd->ccxx[jj][ii] = -cf*para.idx*para.idx*(rnxm+rnxp);
		dd->ccyy[jj][ii] = -cf*para.idy*para.idy*(rnym+rnyp)*2.0;
		
	}
	}
	
}

