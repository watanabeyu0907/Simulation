void boundary_uu(coord *uu)
{
	int ii, jj;


	// B.C. @ x direction
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(jj) \
	shared(para, uu)
	#endif
	for(jj=2; jj<para.ny+2; jj++){
		uu->xx[jj][        0] =  uu->xx[jj][        2];
		uu->xx[jj][        1] =  0.0;
		uu->xx[jj][para.nx+1] =  0.0;
		uu->xx[jj][para.nx+2] =  uu->xx[jj][para.nx  ];
		uu->xx[jj][para.nx+3] =  uu->xx[jj][para.nx-1];

		uu->yy[jj][        0] = +uu->yy[jj][        3];
		uu->yy[jj][        1] = +uu->yy[jj][        2];
		uu->yy[jj][para.nx+2] = +uu->yy[jj][para.nx+1];
		uu->yy[jj][para.nx+3] = +uu->yy[jj][para.nx  ];
		
		// uu->yy[jj][        0] = -uu->yy[jj][        3];
		// uu->yy[jj][        1] = -uu->yy[jj][        2];
		// uu->yy[jj][para.nx+2] = -uu->yy[jj][para.nx+1];
		// uu->yy[jj][para.nx+3] = -uu->yy[jj][para.nx  ];
	}
	
	// B.B. @ y direction	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii) \
	shared(para, uu)
	#endif
	for(ii=1; ii<para.nx+2; ii++){
		uu->xx[        0][ii] = -uu->xx[        3][ii];
		uu->xx[        1][ii] = -uu->xx[        2][ii];
		uu->xx[para.ny+2][ii] = -uu->xx[para.ny+1][ii] + 2.0*para.wu;
		uu->xx[para.ny+3][ii] = -uu->xx[para.ny  ][ii] + 2.0*para.wu;
							   
		uu->yy[        0][ii] =  uu->yy[        2][ii];
		uu->yy[        1][ii] =  0.0;
		uu->yy[para.ny+1][ii] =  0.0;
		uu->yy[para.ny+2][ii] =  uu->yy[para.ny  ][ii];
		uu->yy[para.ny+3][ii] =  uu->yy[para.ny-1][ii];
	}	

	
	// B.C. @ corner
	uu->xx[        0][        0] = uu->xx[        3][        2];
	uu->xx[        1][        0] = uu->xx[        2][        2];
	
	uu->yy[        0][        0] = uu->yy[        2][        3];
	uu->yy[        0][        1] = uu->yy[        2][        2];
	

	uu->xx[        0][para.nx+2] = uu->xx[        3][para.nx  ];
	uu->xx[        0][para.nx+3] = uu->xx[        3][para.nx-1];
	uu->xx[        1][para.nx+2] = uu->xx[        2][para.nx  ];
	uu->xx[        1][para.nx+3] = uu->xx[        2][para.nx-1];
	
	uu->yy[        0][para.nx+2] = uu->yy[        2][para.nx+1];
	uu->yy[        0][para.nx+3] = uu->yy[        2][para.nx  ];
	
	
	uu->xx[para.ny+2][        0] = 2.0*para.wu-uu->xx[para.ny+1][        2];
	uu->xx[para.ny+3][        0] = 2.0*para.wu-uu->xx[para.ny  ][        2];
	
	// uu->yy[para.ny+2][        0] = -uu->yy[para.ny  ][        3];
	// uu->yy[para.ny+2][        1] = -uu->yy[para.ny  ][        2];
	// uu->yy[para.ny+3][        0] = -uu->yy[para.ny-1][        3];
	// uu->yy[para.ny+3][        1] = -uu->yy[para.ny-1][        2];
	
	uu->yy[para.ny+2][        0] = +uu->yy[para.ny  ][        3];
	uu->yy[para.ny+2][        1] = +uu->yy[para.ny  ][        2];
	uu->yy[para.ny+3][        0] = +uu->yy[para.ny-1][        3];
	uu->yy[para.ny+3][        1] = +uu->yy[para.ny-1][        2];
	
	
	uu->xx[para.ny+2][para.nx+2] = 2.0*para.wu-uu->xx[para.ny+1][para.nx  ];
	uu->xx[para.ny+2][para.nx+3] = 2.0*para.wu-uu->xx[para.ny+1][para.nx-1];
	uu->xx[para.ny+3][para.nx+2] = 2.0*para.wu-uu->xx[para.ny  ][para.nx  ];
	uu->xx[para.ny+3][para.nx+3] = 2.0*para.wu-uu->xx[para.ny  ][para.nx-1];
	
	// uu->yy[para.ny+2][para.nx+2] = -uu->yy[para.ny  ][para.nx+1];
	// uu->yy[para.ny+2][para.nx+3] = -uu->yy[para.ny  ][para.nx  ];
	// uu->yy[para.ny+3][para.nx+2] = -uu->yy[para.ny-1][para.nx+1];
	// uu->yy[para.ny+3][para.nx+3] = -uu->yy[para.ny-1][para.nx  ];
	
	uu->yy[para.ny+2][para.nx+2] = +uu->yy[para.ny  ][para.nx+1];
	uu->yy[para.ny+2][para.nx+3] = +uu->yy[para.ny  ][para.nx  ];
	uu->yy[para.ny+3][para.nx+2] = +uu->yy[para.ny-1][para.nx+1];
	uu->yy[para.ny+3][para.nx+3] = +uu->yy[para.ny-1][para.nx  ];
	
}

