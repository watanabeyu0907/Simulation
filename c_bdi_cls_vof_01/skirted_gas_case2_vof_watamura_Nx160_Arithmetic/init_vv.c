void init_vv(coord *vv)
{
	int ii, jj;
	int ex, ey;
	double tx, ty;
	
	
	ex = para.nx+4;
	ey = para.ny+4;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, tx, ty) \
	shared(ex, ey, para, vv)
	#endif
	for(jj=0; jj<ey; jj++){
	for(ii=0; ii<ex; ii++){
		
		tx = ((double)ii - 1.5)*para.dx - 0.5*para.lx*para.nx/(double)(para.nx-0);
		ty = ((double)jj - 1.5)*para.dy - 0.5*para.ly*para.ny/(double)(para.ny-0);
		                             
		if(sqrt(tx*tx+ty*ty) < para.lx*0.5*para.nx/(double)(para.nx-0)){
			// vv->xx[jj][ii] = -ty*para.og;
			// vv->yy[jj][ii] =  tx*para.og;
			vv->xx[jj][ii] = 0;
			vv->yy[jj][ii] = 0;
		}
		else{
			vv->xx[jj][ii] = 0.0;
			vv->yy[jj][ii] = 0.0;
		}
		
	}
	}
	
}

