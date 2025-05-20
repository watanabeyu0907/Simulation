void boundary_pp(double **pp)
{
	int ii, jj;
	int ex, ey;
	
	ex = para.nx+2;
	ey = para.ny+2;


	// B.C. @ x direction
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(jj) \
	shared(para, ey, pp)
	#endif
	for(jj=2; jj<ey; jj++){
		pp[jj][        0] = pp[jj][        3];
		pp[jj][        1] = pp[jj][        2];
		pp[jj][para.nx+2] = pp[jj][para.nx+1];
		pp[jj][para.nx+3] = pp[jj][para.nx  ];
	}
	
	// B.C. @ y direction
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii) \
	shared(para, ex, pp)
	#endif
	for(ii=2; ii<ex; ii++){
		pp[        0][ii] = pp[        3][ii];
		pp[        1][ii] = pp[        2][ii];
		pp[para.ny+2][ii] = pp[para.ny+1][ii];
		pp[para.ny+3][ii] = pp[para.ny  ][ii];		
	}

	// B.C. @ corner
	pp[        0][        0] = pp[        3][        3];
	pp[        0][        1] = pp[        3][        2];
	pp[        1][        0] = pp[        2][        3];
	pp[        1][        1] = pp[        2][        2];
	
	pp[        0][para.nx+2] = pp[        3][para.nx+1];
	pp[        0][para.nx+3] = pp[        3][para.nx  ];
	pp[        1][para.nx+2] = pp[        2][para.nx+1];
	pp[        1][para.nx+3] = pp[        2][para.nx  ];
	
	pp[para.ny+2][        0] = pp[para.ny+1][        3];
	pp[para.ny+2][        1] = pp[para.ny+1][        2];
	pp[para.ny+3][        0] = pp[para.ny  ][        3];
	pp[para.ny+3][        1] = pp[para.ny  ][        2];
	
	pp[para.ny+2][para.nx+2] = pp[para.ny+1][para.nx+1];
	pp[para.ny+2][para.nx+3] = pp[para.ny+1][para.nx  ];
	pp[para.ny+3][para.nx+2] = pp[para.ny  ][para.nx+1];
	pp[para.ny+3][para.nx+3] = pp[para.ny  ][para.nx  ];
	
	
}

