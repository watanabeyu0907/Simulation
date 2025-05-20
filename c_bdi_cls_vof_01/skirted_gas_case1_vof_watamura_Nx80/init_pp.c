void init_pp(double **pp)
{
	int ii, jj;
	int ex, ey;
	
	
	ex = para.nx+4;
	ey = para.ny+4;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex, ey, pp)
	#endif
	for(jj=0; jj<ey; jj++){
	for(ii=0; ii<ex; ii++){
		pp[jj][ii] = 0.0;
	}
	}
	
}

