double calc_bud6(coord *uu, double **phir, double **phim, double **phit, coord *surf)
{
	int ii, jj;	
	int ex, ey;
	double tmp;
	double ans=0.0;
	double pw;
	
	
	ex = para.nx+2;
	ey = para.ny+2;

	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	private(pw, tmp) \
	shared(ex, ey, para, uu, phir, phim, surf) \
	reduction(+:ans)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		tmp = 0.0;
		
		// xx
		pw   = 0.5*(+phir[jj  ][ii-1]+phim[jj  ][ii-1]
					+phir[jj  ][ii  ]+phim[jj  ][ii  ]);
		pw   = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		tmp += 0.5  * (1.0-pw)
					* surf->xx[jj  ][ii-1]
					*   uu->xx[jj  ][ii-1];
		
		pw   = 0.5*(+phir[jj  ][ii  ]+phim[jj  ][ii  ]
					+phir[jj  ][ii+1]+phim[jj  ][ii+1]);
		pw   = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		tmp += 0.5  * (1.0-pw)
					* surf->xx[jj  ][ii  ]
					*   uu->xx[jj  ][ii  ];
		
		
		// yy
		pw   = 0.5*(+phir[jj-1][ii  ]+phim[jj-1][ii  ]
					+phir[jj  ][ii  ]+phim[jj  ][ii  ]);
		pw   = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		tmp += 0.5  * (1.0-pw)
					* surf->yy[jj-1][ii  ]
					*   uu->yy[jj-1][ii  ];
					
		pw   = 0.5*(+phir[jj  ][ii  ]+phim[jj  ][ii  ]
					+phir[jj+1][ii  ]+phim[jj+1][ii  ]);
		pw   = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		tmp += 0.5  * (1.0-pw)
					* surf->yy[jj  ][ii  ]
					*   uu->yy[jj  ][ii  ];
		
		
		ans += tmp*para.ro1*para.dx*para.dy;
		
	}
	}
	
	return(ans);

}

