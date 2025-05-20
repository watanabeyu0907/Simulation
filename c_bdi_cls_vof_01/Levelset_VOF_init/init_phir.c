void init_phir(double **phi)
{
	int ii, jj;
	int ni, nj;
	int ti, tj;
	int ex, ey;
	double dis_xm, dis_xp;
	double dis_ym, dis_yp;
	double rr_mm, rr_mp, rr_pm, rr_pp;
	double rr_ti, rr_tj;
	double vol, dvol;
	double radi_x, radi_y;
	
	
	ex = para.nx+4;
	ey = para.ny+4;
	
	
	ni = nj = 16;
	dvol = 1.0/(double)(ni*nj);
	radi_x = para.lx*0.5*para.Ro;
	radi_y = para.ly*0.5*para.Ro;

	// radi = para.lx*0.5*para.Ri;
	

	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, ti, tj) \
	private(dis_xm, dis_xp, dis_ym, dis_yp) \
	private(rr_mm, rr_mp, rr_pm, rr_pp, rr_ti, rr_tj) \
	private(vol) \
	shared(ex, ey, para, radi_x, radi_y, dvol, ni, nj, phi)
	#endif
	for(jj=0; jj<ey; jj++){
	for(ii=0; ii<ex; ii++){
		
		phi[jj][ii] = 1.0e-9;
		
		
		// @ cell centre
		dis_xm = ((double)ii - 2.0)*para.dx - 0.5*para.lx*para.nx/(double)(para.nx-8);
		dis_xp = ((double)ii - 1.0)*para.dx - 0.5*para.lx*para.nx/(double)(para.nx-8);
		dis_ym = ((double)jj - 2.0)*para.dy - 0.5*para.ly*para.ny/(double)(para.ny-8);
		dis_yp = ((double)jj - 1.0)*para.dy - 0.5*para.ly*para.ny/(double)(para.ny-8);
		
		// rr_mm = sqrt(pow(dis_ym, 2.0) + pow(dis_xm, 2.0));
		// rr_mp = sqrt(pow(dis_ym, 2.0) + pow(dis_xp, 2.0));
		// rr_pm = sqrt(pow(dis_yp, 2.0) + pow(dis_xm, 2.0));
		// rr_pp = sqrt(pow(dis_yp, 2.0) + pow(dis_xp, 2.0));
		
		// rr_mm = MAX(fabs(dis_ym),fabs(dis_xm));
		// rr_mp = MAX(fabs(dis_ym),fabs(dis_xp));
		// rr_pm = MAX(fabs(dis_yp),fabs(dis_xm));
		// rr_pp = MAX(fabs(dis_yp),fabs(dis_xp));
		

		
		vol = 0.0;
		if(fabs(dis_ym) < radi_y
		&& fabs(dis_yp) < radi_y){						
			vol = 1.0;
		}
		else if(fabs(dis_ym) < radi_y
			 || fabs(dis_yp) < radi_y){
			
			for(tj = 0; tj < nj; tj++){
			for(ti = 0; ti < ni; ti++){
				rr_ti = para.dx/(double)ni*((double)ti+0.5) + dis_xm;
				rr_tj = para.dy/(double)nj*((double)tj+0.5) + dis_ym;
				if(fabs(rr_tj) < radi_y){
						vol += dvol;
				}
			}
			}
			
		}
		
		// if(rr_mm < radi
		// && rr_mp < radi
		// && rr_pm < radi
		// && rr_pp < radi){						
			// vol = 1.0;
		// }
		// else if(
		   // rr_mm < radi
		// || rr_mp < radi
		// || rr_pm < radi
		// || rr_pp < radi){
			
			// for(tj = 0; tj < nj; tj++){
			// for(ti = 0; ti < ni; ti++){
				// rr_ti = para.dx/(double)ni*((double)ti+0.5) + dis_xm;
				// rr_tj = para.dy/(double)nj*((double)tj+0.5) + dis_ym;
				// if(sqrt(pow(rr_tj, 2.0) + pow(rr_ti, 2.0)) < radi){
						// vol += dvol;
				// }
			// }
			// }
			
		// }
		
		phi[jj][ii] = MAX(MIN(1.0-vol,1.0-1.0e-9),1.0e-9);
		// phi[jj][ii] = MAX(MIN(    vol,1.0-1.0e-9),1.0e-9);
		
	}
	}
	
}

