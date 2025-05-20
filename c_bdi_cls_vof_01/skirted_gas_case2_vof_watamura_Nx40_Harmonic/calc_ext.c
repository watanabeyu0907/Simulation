extern double calc_wei(double mm, double aa, double bb)
{
	double ans;
	
	ans = 0.5*(1.0+cos(pi*MIN(MAX(aa,-mm),mm)/mm))/mm;
	ans = MIN(MAX(aa,1.0e-9),1.0-1.0e-9);
	if(bb > 1.0e-2) ans=0.0;

	return (ans);
}

extern double calc_dgesv(double *a13, double *b13, double *a7, double *b7, int *ipv, int cnt, double val, int pp)
{
	int ni, nj;
	int n1=1, n3=3, n6=6, n7=7, n12=12, n13=13;
	double ans;
	int info;
	
	if(cnt>=7){
		dgesv(&n13, &n1, a13, &n13, ipv, b13, &n13, &info);
		if(pp==0)	ans = b13[0];
		if(pp==1)	ans = b13[6];
	}
	else if(cnt>=4){
		for(nj=0; nj<n3; nj++){
		for(ni=0; ni<n3; ni++){
			a7[n7*    nj  +     ni ] = a13[n13*    nj  +     ni ];
			a7[n7*    nj  + (n3+ni)] = a13[n13*    nj  + (n6+ni)];
			a7[n7*(n3+nj) +     ni ] = a13[n13*(n6+nj) +     ni ];
			a7[n7*(n3+nj) + (n3+ni)] = a13[n13*(n6+nj) + (n6+ni)];
		}
		}
		for(nj=0; nj<n3; nj++){
			a7[n7*    nj  +     n6 ] = a13[n13*    nj  +    n12 ];
			a7[n7*(n3+nj) +     n6 ] = a13[n13*(n6+nj) +    n12 ];
		}
		for(ni=0; ni<n3; ni++){
			a7[n7*    n6  +     ni ] = a13[n13*   n12  +     ni ];
			a7[n7*    n6  + (n3+ni)] = a13[n13*   n12  + (n6+ni)];
		}
		for(nj=0; nj<n3; nj++){
			b7[   nj] = b13[   nj];
			b7[n3+nj] = b13[n3+nj];
		}
		
		dgesv(&n7, &n1, a7, &n7, ipv, b7, &n7, &info);
		if(pp==0)	ans = b7[0];
		if(pp==1)	ans = b7[3];
	}
	else{
		ans = val;
	}
	
	
	return (ans);
}

void calc_basis(int kk, int ll, double *ff, double *gg)
{
					
	ff[ 0] = 1.0;
	ff[ 1] = (double)kk;
	ff[ 2] = (double)ll;
	ff[ 3] = (double)kk*kk;
	ff[ 4] = (double)ll*ll;
	ff[ 5] = (double)kk*ll;
	
	gg[ 0] = 0.0;
	gg[ 1] = 1.0;
	gg[ 2] = 0.0; 
	gg[ 3] = 2.0*(double)kk;
	gg[ 4] = 0.0;
	gg[ 5] = 1.0*(double)ll;
	gg[ 6] = 0.0;
	gg[ 7] = 0.0;
	gg[ 8] = 1.0;
	gg[ 9] = 0.0;
	gg[10] = 2.0*(double)ll;
	gg[11] = 1.0*(double)kk;
}

void init_mlsmt(double *a13, double *b13)
{
	int ii;
	int n13=13;

	for(ii=0; ii<n13*n13; ii++)	a13[ii] = 0.0;
	for(ii=0; ii<n13    ; ii++)	b13[ii] = 0.0;
}

void calc_mlsmt(int ii, int jj, int kk, int ll, double wei,
				double *ff, double *gg, double *a13, double *b13, 
				coord *uuu, int pp)
{
	int ni, nj;
	int n6=6, n12=12, n13=13;
	
	for(nj=0; nj<n6; nj++){
	for(ni=0; ni<n6; ni++){
		a13[n13*    nj +    ni ] += ff[nj]*ff[ni]*wei;
		a13[n13*(n6+nj)+(n6+ni)] += ff[nj]*ff[ni]*wei;
	}
	}
	for(nj=0; nj<n6; nj++){
		if(pp==0){
			b13[   nj] +=        uuu->xx[jj+ll  ][ii+kk  ] *ff[nj]*wei;
			b13[n6+nj] += 0.25*(+uuu->yy[jj+ll-1][ii+kk  ]
								+uuu->yy[jj+ll-1][ii+kk+1]
								+uuu->yy[jj+ll  ][ii+kk  ]
								+uuu->yy[jj+ll  ][ii+kk+1])*ff[nj]*wei;
		}
		if(pp==1){
			b13[   nj] += 0.25*(+uuu->xx[jj+ll  ][ii+kk-1]
								+uuu->xx[jj+ll+1][ii+kk-1]
								+uuu->xx[jj+ll  ][ii+kk  ]
								+uuu->xx[jj+ll+1][ii+kk  ])*ff[nj]*wei;
			b13[n6+nj] +=        uuu->yy[jj+ll  ][ii+kk  ] *ff[nj]*wei;
		}
	}
	for(nj=0; nj<n12; nj++){
		a13[n13*nj+12] += 0.5*gg[nj]*wei;
	}
	for(ni=0; ni<n12; ni++){
		a13[n13*12+ni] += 0.5*gg[ni]*wei;
	}
	
}

void calc_ext(coord *uuu, coord *euu, double **phi, double **lvs)
{
	int ii, jj;
	int kk, ll;
	int	mm=6, nn=3;
	int n7=7, n13=13;
	int cnt;
	
	double wei;
	double ff[6], gg[12];
	double a13[13*13], b13[13], a7[7*7], b7[7];
	int ipv[13];
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, kk, ll, cnt, wei)  \
	firstprivate(ff, gg, a13, b13, a7, b7, ipv) \
	shared(para, phi, lvs, uuu, euu, mm, nn, n7, n13) \
	schedule(dynamic, 4)
	#endif
	for(jj=0; jj<para.ny+4; jj++){
	for(ii=0; ii<para.nx+3; ii++){
		
		if(fabs(0.5*(lvs[jj][ii] + lvs[jj][ii+1])) > nn + 0.5
		||      0.5*(phi[jj][ii] + phi[jj][ii+1])  < 1.0e-2){
			euu->xx[jj][ii] = uuu->xx[jj][ii];
		}
		else{
			
			cnt = 0;
			init_mlsmt(a13, b13);
			
			for(ll=-mm; ll<=mm; ll++){
			for(kk=-mm; kk<=mm; kk++){
				if(ii+kk>=2 && ii+kk<para.nx+2
				&& jj+ll>=2 && jj+ll<para.ny+2){
					
					wei = calc_wei(
						(double)mm,
						sqrt(ll*ll+kk*kk), 
						0.5*(phi[jj+ll][ii+kk] + phi[jj+ll][ii+kk+1]));
					
					if(wei > 1.0e-2){
						cnt++;
						calc_basis(kk, ll, ff, gg);
						calc_mlsmt(ii, jj, kk, ll, wei, 
							ff, gg, a13, b13, uuu, 0);
					}
				}
			}
			}
			
			euu->xx[jj][ii] = calc_dgesv(a13, b13, a7, b7, ipv, cnt, uuu->xx[jj][ii], 0);
		}
		
	}
	}
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, kk, ll, cnt, wei)  \
	firstprivate(ff, gg, a13, b13, a7, b7, ipv) \
	shared(para, phi, lvs, uuu, euu, mm, nn, n7, n13) \
	schedule(dynamic, 4)
	#endif
	for(jj=0; jj<para.ny+3; jj++){
	for(ii=0; ii<para.nx+4; ii++){
		
		if(fabs(0.5*(lvs[jj][ii] + lvs[jj+1][ii])) > nn + 0.5
		||      0.5*(phi[jj][ii] + phi[jj+1][ii])  < 1.0e-2){
			euu->yy[jj][ii] = uuu->yy[jj][ii];
		}
		else{
			
			cnt = 0;
			init_mlsmt(a13, b13);
			
			for(ll=-mm; ll<=mm; ll++){
			for(kk=-mm; kk<=mm; kk++){
				if(ii+kk>=2 && ii+kk<para.nx+2
				&& jj+ll>=2 && jj+ll<para.ny+2){
					
					wei = calc_wei(
						(double)mm,
						sqrt(ll*ll+kk*kk), 
						0.5*(phi[jj+ll][ii+kk] + phi[jj+ll+1][ii+kk]));
					
					if(wei > 1.0e-2){
						cnt++;
						calc_basis(kk, ll, ff, gg);
						calc_mlsmt(ii, jj, kk, ll, wei, 
							ff, gg, a13, b13, uuu, 1);
					}
				}
			}
			}
			
			euu->yy[jj][ii] = calc_dgesv(a13, b13, a7, b7, ipv, cnt, uuu->yy[jj][ii], 1);
		}
		
	}
	}
	
}

