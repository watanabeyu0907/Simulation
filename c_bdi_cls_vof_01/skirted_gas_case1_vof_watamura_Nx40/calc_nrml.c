void calc_nrml(double **phir, double **phim, double **phif, double **phil, double **phit,
	coord *rnf, coord *rns
)
{
	int ii, jj;
	int ex3, ey3;
	double q00, q01, q10, q11;
	double tmp;
	
	
	ex3 = para.nx+3;
	ey3 = para.ny+3;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, tmp) \
	private(q00, q01, q10, q11) \
	shared(ex3, ey3, para) \
	shared(phir, phim, phif, phil, phit, rnf, rns)
	#endif
	for(jj=0; jj<ey3; jj++){
	for(ii=0; ii<ex3; ii++){
		
		tmp = 0.25*(+phit[jj  ][ii  ]
					+phit[jj  ][ii+1]
					+phit[jj+1][ii  ]
					+phit[jj+1][ii+1]);
		
		if(tmp < 0.5){
			q00 = phif[jj  ][ii  ];
			q01 = phif[jj  ][ii+1];
			q10 = phif[jj+1][ii  ];
			q11 = phif[jj+1][ii+1];
		}
		else{
			q00 = phil[jj  ][ii  ];
			q01 = phil[jj  ][ii+1];
			q10 = phil[jj+1][ii  ];
			q11 = phil[jj+1][ii+1];
		}
		
		rnf->xx[jj][ii] = (-q00+q01-q10+q11)*para.idx*0.5;
		rnf->yy[jj][ii] = (-q00-q01+q10+q11)*para.idy*0.5;
	
	
	
		q00 = phir[jj  ][ii  ] + phim[jj  ][ii  ];
		q01 = phir[jj  ][ii+1] + phim[jj  ][ii+1];
		q10 = phir[jj+1][ii  ] + phim[jj+1][ii  ];
		q11 = phir[jj+1][ii+1] + phim[jj+1][ii+1];
		
		rns->xx[jj][ii] = (-q00+q01-q10+q11)*para.idx*0.5;
		rns->yy[jj][ii] = (-q00-q01+q10+q11)*para.idy*0.5;
		
	}
	}
	
}

