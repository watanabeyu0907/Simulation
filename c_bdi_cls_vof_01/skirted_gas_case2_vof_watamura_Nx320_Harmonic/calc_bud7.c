double calc_bud7(coord *uu, double **pp, double **phir, double **phim, double **phit)
{
	int ii, jj;	
	int ex, ey;
	double tmp;
	double ans=0.0;
	double ph, pw, ro, pd;
	
	
	ex = para.nx+2;
	ey = para.ny+2;

	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	private(ph, pw, ro, pd, tmp) \
	shared(ex, ey, para, uu, pp, phir, phim, phit) \
	reduction(+:ans)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		tmp = 0.0;
		
		// xx
		ph = 0.5*(	+phit[jj  ][ii-1]
					+phit[jj  ][ii  ]);
		pw = 0.5*(	+phir[jj  ][ii-1]+phim[jj  ][ii-1]
					+phir[jj  ][ii  ]+phim[jj  ][ii  ]);
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		
		pd = (-pp[jj  ][ii-1]+pp[jj  ][ii  ])*para.ro1*para.idx;
		tmp += 0.5 * ph * (1.0-pw) * pd * uu->xx[jj  ][ii-1];
		
		// xx
		ph = 0.5*(	+phit[jj  ][ii  ]
					+phit[jj  ][ii+1]);
		pw = 0.5*(	+phir[jj  ][ii  ]+phim[jj  ][ii  ]
					+phir[jj  ][ii+1]+phim[jj  ][ii+1]);
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);

		pd = (-pp[jj  ][ii  ]+pp[jj  ][ii+1])*para.ro1*para.idx;
		tmp += 0.5 * ph * (1.0-pw) * pd * uu->xx[jj  ][ii  ];
		
		// yy
		ph = 0.5*(	+phit[jj-1][ii  ]
					+phit[jj  ][ii  ]);
		pw = 0.5*(	+phir[jj-1][ii  ]+phim[jj-1][ii  ]
					+phir[jj  ][ii  ]+phim[jj  ][ii  ]);
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		ro = ph*para.ro1;
		ro = para.ro1 + ph*(para.ro2-para.ro1);

		
		pd = (-pp[jj-1][ii  ]+pp[jj  ][ii  ])*para.ro1*para.idy;
		tmp += 0.5 * ph * (1.0-pw) * pd * uu->yy[jj-1][ii  ];
		tmp += 0.5 * ph * ro * (1.0-pw) * para.gg * uu->yy[jj-1][ii  ];
		
		// yy
		ph = 0.5*(	+phit[jj  ][ii  ]
					+phit[jj+1][ii  ]);
		pw = 0.5*(	+phir[jj  ][ii  ]+phim[jj  ][ii  ]
					+phir[jj+1][ii  ]+phim[jj+1][ii  ]);
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		ro = ph*para.ro1;
		ro = para.ro1 + ph*(para.ro2-para.ro1);

		pd = (-pp[jj  ][ii  ]+pp[jj+1][ii  ])*para.ro1*para.idy;
		tmp += 0.5 * ph * (1.0-pw) * pd * uu->yy[jj  ][ii  ];
		tmp += 0.5 * ph * ro * (1.0-pw) * para.gg * uu->yy[jj  ][ii  ];
		
		
		
		ans += tmp*para.dx*para.dy;
		
	}
	}
	
	return(ans);

}

