void calc_lvsm(double **lvs0, double **lvs3, double tht)
{
	int ii, jj;
	int ti, tj;
	int si, sj;
	int ex, ey;
	double cs, sn;
	double ax[4], ay[4];
	double ans;
	coord0 p3, p0;
	coord0 dd;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	cs  = cos(tht);
	sn  = sin(tht);
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj) \
	private(p0, p3, ti, tj, dd) \
	private(si, sj, ax, ay, ans) \
	shared(para, ex, ey, sn, cs, lvs0, lvs3)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		p3.xx = ((double)ii-1.5)*para.dx - 0.5*para.lx*para.nx/(double)(para.nx-8);
		p3.yy = ((double)jj-1.5)*para.dy - 0.5*para.ly*para.ny/(double)(para.ny-8);
		
		p0.xx =  p3.xx*cs + p3.yy*sn + 0.5*para.lx*para.nx/(double)(para.nx-8);
		p0.yy = -p3.xx*sn + p3.yy*cs + 0.5*para.ly*para.ny/(double)(para.ny-8);
		
		p0.xx = p0.xx*para.idx + 1.5;
		p0.yy = p0.yy*para.idy + 1.5;
		
		ti = (int)(p0.xx);
		tj = (int)(p0.yy);
		
		if(ti >= 1 && ti < para.nx-1
		&& tj >= 1 && tj < para.ny-1){
			
			dd.xx = p0.xx - (double)ti + 1;
			dd.yy = p0.yy - (double)tj + 1;
			
			ax[0] =             (dd.xx-1.0)*(dd.xx-2.0)*(dd.xx-3.0)
				  /(            (0.0  -1.0)*(0.0  -2.0)*(0.0  -3.0));
			ax[1] = (dd.xx-0.0)            *(dd.xx-2.0)*(dd.xx-3.0)
				  /((1.0  -0.0)            *(1.0  -2.0)*(1.0  -3.0));
			ax[2] = (dd.xx-0.0)*(dd.xx-1.0)            *(dd.xx-3.0)
				  /((2.0  -0.0)*(2.0  -1.0)            *(2.0  -3.0));
			ax[3] = (dd.xx-0.0)*(dd.xx-1.0)*(dd.xx-2.0)
				  /((3.0  -0.0)*(3.0  -1.0)*(3.0  -2.0)            );
			
			ay[0] =             (dd.yy-1.0)*(dd.yy-2.0)*(dd.yy-3.0)
				  /(            (0.0  -1.0)*(0.0  -2.0)*(0.0  -3.0));
			ay[1] = (dd.yy-0.0)            *(dd.yy-2.0)*(dd.yy-3.0)
				  /((1.0  -0.0)            *(1.0  -2.0)*(1.0  -3.0));
			ay[2] = (dd.yy-0.0)*(dd.yy-1.0)            *(dd.yy-3.0)
				  /((2.0  -0.0)*(2.0  -1.0)            *(2.0  -3.0));
			ay[3] = (dd.yy-0.0)*(dd.yy-1.0)*(dd.yy-2.0)
				  /((3.0  -0.0)*(3.0  -1.0)*(3.0  -2.0)            );
				  
			
			ans = 0.0;
			for(si=0; si<4; si++){
			for(sj=0; sj<4; sj++){
				ans += ay[sj]*ax[si]*lvs0[tj+sj-1][ti+si-1];
			}
			}
			lvs3[jj][ii] = ans;
			
		}
		
	}
	}
	
}

