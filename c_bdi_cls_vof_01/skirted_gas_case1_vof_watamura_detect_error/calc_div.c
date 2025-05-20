void calc_div(coord *uu, double **div)
{
	int ii, jj;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(para, uu, div)
	#endif
	for(jj=2; jj<para.ny+2; jj++){
	for(ii=2; ii<para.nx+2; ii++){

		// calculation of divergence
		div[jj][ii] = +(-uu->xx[jj  ][ii-1]
				        +uu->xx[jj  ][ii  ])*para.idx
				      +(-uu->yy[jj-1][ii  ]
				        +uu->yy[jj  ][ii  ])*para.idy;
		
		div[jj][ii] *= para.idT; 
					 
	}
	}

}

