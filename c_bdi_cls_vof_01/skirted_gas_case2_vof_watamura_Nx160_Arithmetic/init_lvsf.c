void init_lvsf(double **phit, double **lvf3)
{
	int ii, jj;
	int ex4, ey4;
	
	
	ex4 = para.nx+4;
	ey4 = para.ny+4;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(ex4, ey4, phit, lvf3)
	#endif
	for(jj=0; jj<ey4; jj++){
	for(ii=0; ii<ex4; ii++){
		
		lvf3[jj][ii] = phit[jj][ii]-0.5;
		
	}
	}
		
}

