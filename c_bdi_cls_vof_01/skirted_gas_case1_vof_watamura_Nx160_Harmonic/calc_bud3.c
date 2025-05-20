double calc_bud3(coord *uu)
{
	int ii, jj;	
	double tmp=0.0;
	double ans=0.0;
	
	double dudx, dudy;
	double dvdx, dvdy;
	
	double ss11, ss12;
	double ss21, ss22;
	double velx, vely;
	

	// boundary @ left
	#ifdef _OPENMP
	#pragma omp parallel for \
	private(dudx, dudy, dvdx, dvdy, ss11, ss12, ss21, ss22, velx, vely) \
	reduction(+:tmp)
	#endif
	for(jj=2; jj<para.ny+2; jj++){
		
		// dudy + dvdx
		// 1/2
		vely = +uu->yy[jj-1][    1]
			   +uu->yy[jj-1][    2];
		vely*= 0.5;		
		dudy = -uu->xx[jj-1][    1]
			   +uu->xx[jj  ][    1];
		dvdx = -uu->yy[jj-1][    1]
			   +uu->yy[jj-1][    2];
		ss21 = 0.5*dudy+0.5*dvdx;
		
		tmp += -0.5*ss21*vely*para.dy;
		
		// 2/2
		vely = +uu->yy[jj  ][    1]
			   +uu->yy[jj  ][    2];
		vely*= 0.5;		
		dudy = -uu->xx[jj  ][    1]
			   +uu->xx[jj+1][    1];
		dvdx = -uu->yy[jj  ][    1]
			   +uu->yy[jj  ][    2];
		ss21 = 0.5*dudy+0.5*dvdx;
		
		tmp += -0.5*ss21*vely*para.dy;
		
		
		// dvdy
		velx = +uu->yy[jj  ][    1];
		dudx = -uu->xx[jj  ][    0]
			   +uu->xx[jj-1][    2];
		dudx *= para.idx*0.5;
		ss11 = dudx;
		
		tmp += -ss11*vely*para.dy;
	}
	
	// boundary @ right
	#ifdef _OPENMP
	#pragma omp parallel for \
	private(dudx, dudy, dvdx, dvdy, ss11, ss12, ss21, ss22, velx, vely) \
	reduction(+:tmp)
	#endif
	for(jj=2; jj<para.ny+2; jj++){
		
		// dudy + dvdx
		// 1/2
		vely = +uu->yy[jj-1][para.nx+1]
			   +uu->yy[jj-1][para.nx+2];
		vely*= 0.5;		
		dudy = -uu->xx[jj-1][para.nx+1]
			   +uu->xx[jj  ][para.nx+1];
		dvdx = -uu->yy[jj-1][para.nx+1]
			   +uu->yy[jj-1][para.nx+2];
		ss21 = 0.5*dudy+0.5*dvdx;
		
		tmp += +0.5*ss21*vely*para.dy;
		
		// 2/2
		vely = +uu->yy[jj  ][para.nx+1]
			   +uu->yy[jj  ][para.nx+2];
		vely*= 0.5;		
		dudy = -uu->xx[jj  ][para.nx+1]
			   +uu->xx[jj+1][para.nx+1];
		dvdx = -uu->yy[jj  ][para.nx+1]
			   +uu->yy[jj  ][para.nx+2];
		ss21 = 0.5*dudy+0.5*dvdx;
		
		tmp += +0.5*ss21*vely*para.dy;
		
		
		// dvdy
		velx = +uu->yy[jj  ][para.nx+1];
		dudx = -uu->xx[jj  ][para.nx  ]
			   +uu->xx[jj  ][para.nx+2];
		dudx *= para.idx*0.5;
		ss11 = dudx;
		
		tmp += +ss11*vely*para.dy;
	}

	// boundary @ low
	#ifdef _OPENMP
	#pragma omp parallel for \
	private(dudx, dudy, dvdx, dvdy, ss11, ss12, ss21, ss22, velx, vely) \
	reduction(+:tmp)
	#endif
	for(ii=2; ii<para.nx+2; ii++){
		
		// dudy + dvdx
		// 1/2
		velx = +uu->xx[   1][ii-1]
			   +uu->xx[   2][ii-1];
		velx*= 0.5;		
		dudy = -uu->xx[   1][ii-1]
			   +uu->xx[   2][ii-1];
		dvdx = -uu->yy[   1][ii-1]
			   +uu->yy[   1][ii  ];
		dudy *= para.idy;
		dvdx *= para.idx;
		ss12 = 0.5*dudy+0.5*dvdx;
		
		tmp += -0.5*ss12*velx*para.dx;
		
		// 2/2
		velx = +uu->xx[   1][ii  ]
			   +uu->xx[   2][ii  ];
		velx*= 0.5;		
		dudy = -uu->xx[   1][ii  ]
			   +uu->xx[   2][ii  ];
		dvdx = -uu->yy[   1][ii  ]
			   +uu->yy[   1][ii+1];
		dudy *= para.idy;
		dvdx *= para.idx;
		ss12 = 0.5*dudy+0.5*dvdx;
		
		tmp += -0.5*ss12*velx*para.dx;
		
		
		// dvdy
		vely = +uu->yy[   1][ii  ];
		dvdy = -uu->yy[   0][ii  ]
			   +uu->yy[   2][ii  ];
		dvdy *= para.idy*0.5;
		ss22 = dvdy;
		
		tmp += -ss22*vely*para.dx;
	}
	
	// boundary @ top
	#ifdef _OPENMP
	#pragma omp parallel for \
	private(dudx, dudy, dvdx, dvdy, ss11, ss12, ss21, ss22, velx, vely) \
	reduction(+:tmp)
	#endif
	for(ii=2; ii<para.nx+2; ii++){
		
		// dudy + dvdx
		// 1/2
		velx = +uu->xx[para.ny+1][ii-1]
			   +uu->xx[para.ny+2][ii-1];
		velx*= 0.5;
		dudy = -uu->xx[para.ny+1][ii-1]
			   +uu->xx[para.ny+2][ii-1];
		dvdx = -uu->yy[para.ny+1][ii-1]
			   +uu->yy[para.ny+1][ii  ];
		dudy *= para.idy;
		dvdx *= para.idx;
		ss12 = 0.5*dudy+0.5*dvdx;
		
		tmp += +0.5*ss12*velx*para.dx;
		
		// 2/2
		velx = +uu->xx[para.ny+1][ii  ]
			   +uu->xx[para.ny+2][ii  ];
		velx*= 0.5;
		dudy = -uu->xx[para.ny+1][ii  ]
			   +uu->xx[para.ny+2][ii  ];
		dvdx = -uu->yy[para.ny+1][ii  ]
			   +uu->yy[para.ny+1][ii+1];
		dudy *= para.idy;
		dvdx *= para.idx;
		ss12 = 0.5*dudy+0.5*dvdx;
		
		tmp += +0.5*ss12*velx*para.dx;
		
		
		// dvdy
		vely = +uu->yy[para.ny+1][ii  ];
		dvdy = -uu->yy[para.ny  ][ii  ]
			   +uu->yy[para.ny+2][ii  ];
		dvdy *= para.idy*0.5;
		ss22 = dvdy;
		
		tmp += +ss22*vely*para.dx;
	}

	ans = tmp*para.mu1;
	
	return(ans);

}

