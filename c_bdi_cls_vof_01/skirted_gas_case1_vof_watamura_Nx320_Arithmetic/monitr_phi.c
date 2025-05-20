void monitr_phi(int tt, double **data, double **phi, int num)
{
	int ii, jj;
	double sum=0.0;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	shared(para, phi) \
	reduction(+:sum)
	#endif
	for(jj=2; jj<para.ny+2; jj++){
	for(ii=2; ii<para.nx+2; ii++){
		sum += phi[jj][ii];
	}
	}
	
	data[tt][15+num] = sum*(para.dx*para.dy);
	
}

