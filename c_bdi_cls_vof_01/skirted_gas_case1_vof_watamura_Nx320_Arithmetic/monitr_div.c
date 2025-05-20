void monitr_div(int tt, double **data, double **div)
{
	int ii, jj;
	
	double sum = 0.0;
	double max = 0.0;
	
	
	#pragma omp parallel \
	private(ii) \
	reduction(+: sum)
	{
		double tmp = -1.0;
		#pragma omp for
		for(jj=2; jj<para.ny+2; jj++){
		for(ii=2; ii<para.nx+2; ii++){
			sum += div[jj][ii]*div[jj][ii];
			if(div[jj][ii]*div[jj][ii] > tmp){
				tmp = div[jj][ii]*div[jj][ii];
			}
		}
		}
		#pragma omp critical
		{
			if(tmp > max) max = tmp;
		}
	}
	sum /= (double)(para.nx*para.ny);
	sum = sqrt(sum)*para.dT;
	max = sqrt(max)*para.dT;
	// printf("div_L2  = %le\n", sum);
	// printf("div_inf = %le\n", max);
	
	data[tt][0] = sum;
	data[tt][1] = max;
}

