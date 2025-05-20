double calc_bud2(coord *uu, double **phir, double **phim, double **phit)
{
	int ii, jj;	
	int ex, ey;
	double ans=0.0;
	double ph, pw, mu, mu1, mu2;
	
	double dudx, dudy;
	double dvdx, dvdy;
	
	double ss11, ss12;
	double ss21, ss22;
	
	double dis;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	mu1 = para.mu1;
	mu2 = para.mu2;

	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, dudx, dudy, dvdx, dvdy, ss11, ss12, ss21, ss22, dis) \
	private(ph, pw, mu) \
	shared(ex, ey, para, uu, phir, phim, phit, mu1, mu2) \
	reduction(+:ans)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		dis = 0.0;
		
		// uu
		// dudx + dvdy
		dudx = -uu->xx[jj  ][ii-1]
			   +uu->xx[jj  ][ii  ];
			   
		dvdy = -uu->yy[jj-1][ii  ]
			   +uu->yy[jj  ][ii  ];
		
		dudx *= para.idx;
		dvdy *= para.idy;
		
		ss11 = dudx;
		ss22 = dvdy;
		
		ph = phit[jj][ii];
		pw = phir[jj][ii]+phim[jj][ii];
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		mu = mu1 / (1.0 + ph*(mu1/mu2-1.0));
		
		dis += +0.5*mu*(1.0-pw)*ss11*ss11
			   +0.5*mu*(1.0-pw)*ss22*ss22;
		
		// dudy + dvdx
		// 1/4
		dudy = -uu->xx[jj-1][ii-1]
			   +uu->xx[jj  ][ii-1];
			   
		dvdx = -uu->yy[jj-1][ii-1]
			   +uu->yy[jj-1][ii  ];
		
		dudy *= para.idy;
		dvdx *= para.idx;	   
			   
		ss12 = 0.5*dudy+0.5*dvdx;
		ss21 = 0.5*dudy+0.5*dvdx;
		
		ph = 0.25*(	+phit[jj-1][ii-1]
					+phit[jj-1][ii  ]
					+phit[jj  ][ii-1]
					+phit[jj  ][ii  ]);
		pw = 0.25*(	+phir[jj-1][ii-1]+phim[jj-1][ii-1]
					+phir[jj-1][ii  ]+phim[jj-1][ii  ]
					+phir[jj  ][ii-1]+phim[jj  ][ii-1]
					+phir[jj  ][ii  ]+phim[jj  ][ii  ]);
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		mu = mu1 / (1.0 + ph*(mu1/mu2-1.0));
		
		dis += +0.125*mu*(1.0-pw)*ss12*ss12
			   +0.125*mu*(1.0-pw)*ss21*ss21;
		
		// 2/4
		dudy = -uu->xx[jj-1][ii  ]
			   +uu->xx[jj  ][ii  ];
			   
		dvdx = -uu->yy[jj-1][ii  ]
			   +uu->yy[jj-1][ii+1];
		
		dudy *= para.idy;
		dvdx *= para.idx;	   
			   
		ss12 = 0.5*dudy+0.5*dvdx;
		ss21 = 0.5*dudy+0.5*dvdx;
		
		ph = 0.25*(	+phit[jj-1][ii  ]
					+phit[jj-1][ii+1]
					+phit[jj  ][ii  ]
					+phit[jj  ][ii+1]);
		pw = 0.25*(	+phir[jj-1][ii  ]+phim[jj-1][ii  ]
					+phir[jj-1][ii+1]+phim[jj-1][ii+1]
					+phir[jj  ][ii  ]+phim[jj  ][ii  ]
					+phir[jj  ][ii+1]+phim[jj  ][ii+1]);
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		mu = mu1 / (1.0 + ph*(mu1/mu2-1.0));
		
		dis += +0.125*mu*(1.0-pw)*ss12*ss12
			   +0.125*mu*(1.0-pw)*ss21*ss21;
		
		// 3/4
		dudy = -uu->xx[jj  ][ii-1]
			   +uu->xx[jj+1][ii-1];
			   
		dvdx = -uu->yy[jj  ][ii-1]
			   +uu->yy[jj  ][ii  ];
		
		dudy *= para.idy;
		dvdx *= para.idx;	   
			   
		ss12 = 0.5*dudy+0.5*dvdx;
		ss21 = 0.5*dudy+0.5*dvdx;
		
		ph = 0.25*(	+phit[jj  ][ii-1]
					+phit[jj  ][ii  ]
					+phit[jj+1][ii-1]
					+phit[jj+1][ii  ]);
		pw = 0.25*(	+phir[jj  ][ii-1]+phim[jj  ][ii-1]
					+phir[jj  ][ii  ]+phim[jj  ][ii  ]
					+phir[jj+1][ii-1]+phim[jj+1][ii-1]
					+phir[jj+1][ii  ]+phim[jj+1][ii  ]);
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		mu = mu1 / (1.0 + ph*(mu1/mu2-1.0));
		
		dis += +0.125*mu*(1.0-pw)*ss12*ss12
			   +0.125*mu*(1.0-pw)*ss21*ss21;
		
		// 4/4
		dudy = -uu->xx[jj  ][ii  ]
			   +uu->xx[jj+1][ii  ];
			   
		dvdx = -uu->yy[jj  ][ii  ]
			   +uu->yy[jj  ][ii+1];
		
		dudy *= para.idy;
		dvdx *= para.idx;	   
			   
		ss12 = 0.5*dudy+0.5*dvdx;
		ss21 = 0.5*dudy+0.5*dvdx;
		
		ph = 0.25*(	+phit[jj  ][ii  ]
					+phit[jj  ][ii+1]
					+phit[jj+1][ii  ]
					+phit[jj+1][ii+1]);
		pw = 0.25*(	+phir[jj  ][ii  ]+phim[jj  ][ii  ]
					+phir[jj  ][ii+1]+phim[jj  ][ii+1]
					+phir[jj+1][ii  ]+phim[jj+1][ii  ]
					+phir[jj+1][ii+1]+phim[jj+1][ii+1]);
		ph = MAX(MIN(ph,1.0-1.0e-9),1.0e-9);
		pw = MAX(MIN(pw,1.0-1.0e-9),1.0e-9);
		mu = mu1 / (1.0 + ph*(mu1/mu2-1.0));
		
		dis += +0.125*mu*(1.0-pw)*ss12*ss12
			   +0.125*mu*(1.0-pw)*ss21*ss21;
			   
		
		
		// dis *= (1.0-phir[jj][ii]-phim[jj][ii]);
		
		ans += -dis*2.0*para.dx*para.dy;
		
	}
	}
	
	return(ans);

}

