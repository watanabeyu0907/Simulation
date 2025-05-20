double get_hev(double xx, double *hev)
{
	int ii;
	double aa, dd, rr;
	
	
	xx = MIN(MAX(xx,-2.0),2.0);
	
	dd = (xx+2)*hevnum;
	ii = (int)dd;
	rr = dd - (double)ii;
	
	aa = +(1.0 - rr)*hev[ii  ]
		 +(      rr)*hev[ii+1];

	aa = MIN(MAX(aa,1e-99),1.0-1e-99);

	return(aa);
}

