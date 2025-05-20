double calc_bud5(coord *uu, double **phir, double **phim, double **phit)
{
	int ii, jj;	
	int ex, ey;
	double tmp;
	double ans=0.0;
	double ph, pw, ro;
	
	
	ex = para.nx+2;
	ey = para.ny+2;

	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	private(ph, pw, ro, tmp) \
	shared(ex, ey, para, uu, phir, phim, phit) \
	reduction(+:ans)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		tmp = 0.0;
		
		// 1/2
		ph = 0.5*(	+phit[jj-1][ii  ]
					+phit[jj  ][ii  ]);
		pw = 0.5*(	+phir[jj-1][ii  ]+phim[jj-1][ii  ]
					+phir[jj  ][ii  ]+phim[jj  ][ii  ]);
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		ro = ph*para.ro1;
		ro = para.ro1 + ph*(para.ro2-para.ro1);

		tmp += 0.5 * ro * (1.0-pw) * para.gg * uu->yy[jj-1][ii  ];
		
		// 2/2
		ph = 0.5*(	+phit[jj  ][ii  ]
					+phit[jj+1][ii  ]);
		pw = 0.5*(	+phir[jj  ][ii  ]+phim[jj  ][ii  ]
					+phir[jj+1][ii  ]+phim[jj+1][ii  ]);
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		ro = ph*para.ro1;
		ro = para.ro1 + ph*(para.ro2-para.ro1);

		tmp += 0.5 * ro * (1.0-pw) * para.gg * uu->yy[jj  ][ii  ];
		
		
		ans -= tmp*para.dx*para.dy;
		
	}
	}
	
	return(ans);

}

