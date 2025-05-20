# include "watasys.h"
# include "watafnc.h"

void calc_bud4(coord *vv, coord *uu, double **pp, double **dl, double **lv, double *ans)
{
	int ii, jj;	
	int ex, ey;
	double tq_prs, fx_prs, fy_prs;
	double tq_vis, fx_vis, fy_vis;
	double atqp=0.0, afxp=0.0, afyp=0.0;
	double atqv=0.0, afxv=0.0, afyv=0.0;
	
	double dhdx, dhdy;
	double dudx, dudy;
	double dvdx, dvdy;
	double velx, vely;
	
	
	ex = para.nx+2;
	ey = para.ny+2;
	
	
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, dhdx, dhdy, dudx, dudy, dvdx, dvdy, velx, vely) \
	private(tq_prs, fx_prs, fy_prs, tq_vis, fx_vis, fy_vis) \
	shared(ex, ey, para, vv, uu, pp, dl, lv) \
	reduction(+:atqp, afxp, afyp, atqv, afxv, afyv) 
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
		
		tq_prs = 0.0;	fx_prs = 0.0;	fy_prs = 0.0;
		tq_vis = 0.0;	fx_vis = 0.0;	fy_vis = 0.0;
		
		// pp
		// 1
		velx   =  vv->xx[jj  ][ii-1];
		dhdx   = (-lv[jj  ][ii-1] + lv[jj  ][ii  ])*para.idx;
		tq_prs+= - velx * dhdx * pp[jj][ii] * para.ro * 0.5;
		fx_prs+= -        dhdx * pp[jj][ii] * para.ro * 0.5;
		// 2
		velx   =  vv->xx[jj  ][ii  ];
		dhdx   = (-lv[jj  ][ii  ] + lv[jj  ][ii+1])*para.idx; 
		tq_prs+= - velx * dhdx * pp[jj][ii] * para.ro * 0.5;
		fx_prs+= -        dhdx * pp[jj][ii] * para.ro * 0.5;
		
		// 1
		vely   =  vv->yy[jj-1][ii  ];
		dhdy   = (-lv[jj-1][ii  ] + lv[jj  ][ii  ])*para.idy;
		tq_prs+= - vely * dhdy * pp[jj][ii] * para.ro * 0.5;
		fy_prs+= -        dhdy * pp[jj][ii] * para.ro * 0.5;
		// 2
		vely   =  vv->yy[jj  ][ii  ];
		dhdy   = (-lv[jj  ][ii  ] + lv[jj+1][ii  ])*para.idy; 
		tq_prs+= - vely * dhdy * pp[jj][ii] * para.ro * 0.5;
		fy_prs+= -        dhdy * pp[jj][ii] * para.ro * 0.5;
		
		
		// ss11		
		// 1		
		velx   =  vv->xx[jj  ][ii-1];
		dudx   = (uu->xx[jj  ][ii  ] - uu->xx[jj  ][ii-1])*para.idx;
		dhdx   = (-lv[jj  ][ii-1] +  lv[jj  ][ii  ])*para.idx; 
		tq_vis+= velx * dhdx * dudx * para.mu;
		fx_vis+=        dhdx * dudx * para.mu;		
		// 2
		velx   =  vv->xx[jj  ][ii  ];
		dudx   = (uu->xx[jj  ][ii  ] - uu->xx[jj  ][ii-1])*para.idx;
		dhdx   = (-lv[jj  ][ii  ] +  lv[jj  ][ii+1])*para.idx;
		tq_vis+= velx * dhdx * dudx * para.mu;
		fx_vis+=        dhdx * dudx * para.mu;
		
		// ss22
		// 1
		vely   =  vv->yy[jj-1][ii  ];
		dvdy   = (uu->yy[jj  ][ii  ] - uu->yy[jj-1][ii  ])*para.idy;
		dhdy   = (-lv[jj-1][ii  ] +  lv[jj  ][ii  ])*para.idy;
		tq_vis+= vely * dhdy * dvdy * para.mu;		
		fy_vis+=        dhdy * dvdy * para.mu;		
		// 2
		vely   =  vv->yy[jj  ][ii  ];
		dvdy   = (uu->yy[jj  ][ii  ] - uu->yy[jj-1][ii  ])*para.idy;
		dhdy   = (-lv[jj  ][ii  ] +  lv[jj+1][ii  ])*para.idy; 
		tq_vis+= vely * dhdy * dvdy * para.mu;
		fy_vis+=        dhdy * dvdy * para.mu;		
		
		// ss12
		// 1		
		velx   = ( vv->xx[jj-1][ii-1] + vv->xx[jj  ][ii-1])*0.5;
		dudy   = (-uu->xx[jj-1][ii-1] + uu->xx[jj  ][ii-1])*para.idy;
		dvdx   = (-uu->yy[jj-1][ii-1] + uu->yy[jj-1][ii  ])*para.idx;
		dhdy   = ( -lv[jj-1][ii  ] +  lv[jj  ][ii  ])*para.idy;
		tq_vis+= velx * dhdy * (dudy + dvdx) * para.mu * 0.25;
		fx_vis+=        dhdy * (dudy + dvdx) * para.mu * 0.25;
		// 2
		velx   = ( vv->xx[jj-1][ii  ] + vv->xx[jj  ][ii  ])*0.5;
		dudy   = (-uu->xx[jj-1][ii  ] + uu->xx[jj  ][ii  ])*para.idy;
		dvdx   = (-uu->yy[jj-1][ii  ] + uu->yy[jj-1][ii+1])*para.idx;
		dhdy   = ( -lv[jj-1][ii  ] +  lv[jj  ][ii  ])*para.idy;
		tq_vis+= velx * dhdy * (dudy + dvdx) * para.mu * 0.25;
		fx_vis+=        dhdy * (dudy + dvdx) * para.mu * 0.25;
		// 3
		velx   = ( vv->xx[jj  ][ii-1] + vv->xx[jj+1][ii-1])*0.5;
		dudy   = (-uu->xx[jj  ][ii-1] + uu->xx[jj+1][ii-1])*para.idy;
		dvdx   = (-uu->yy[jj  ][ii-1] + uu->yy[jj  ][ii  ])*para.idx;
		dhdy   = ( -lv[jj  ][ii  ] +  lv[jj+1][ii  ])*para.idy;
		tq_vis+= velx * dhdy * (dudy + dvdx) * para.mu * 0.25;
		fx_vis+=        dhdy * (dudy + dvdx) * para.mu * 0.25;
		// 4
		velx   = ( vv->xx[jj  ][ii  ] + vv->xx[jj+1][ii  ])*0.5;
		dudy   = (-uu->xx[jj  ][ii  ] + uu->xx[jj+1][ii  ])*para.idy;
		dvdx   = (-uu->yy[jj  ][ii  ] + uu->yy[jj  ][ii+1])*para.idx;
		dhdy   = ( -lv[jj  ][ii  ] +  lv[jj+1][ii  ])*para.idy;
		tq_vis+= velx * dhdy * (dudy + dvdx) * para.mu * 0.25;
		fx_vis+=        dhdy * (dudy + dvdx) * para.mu * 0.25;
		
		// ss21
		// 1
		vely   = ( vv->yy[jj-1][ii-1] + vv->yy[jj-1][ii  ])*0.5;
		dudy   = (-uu->xx[jj-1][ii-1] + uu->xx[jj  ][ii-1])*para.idy;
		dvdx   = (-uu->yy[jj-1][ii-1] + uu->yy[jj-1][ii  ])*para.idx;
		dhdx   = ( -lv[jj  ][ii-1] +  lv[jj  ][ii  ])*para.idx;
		tq_vis+= vely * dhdx * (dudy + dvdx) * para.mu * 0.25;
		fy_vis+=        dhdx * (dudy + dvdx) * para.mu * 0.25;
		// 2
		vely   = ( vv->yy[jj-1][ii  ] + vv->yy[jj-1][ii+1])*0.5;
		dudy   = (-uu->xx[jj-1][ii  ] + uu->xx[jj  ][ii  ])*para.idy;
		dvdx   = (-uu->yy[jj-1][ii  ] + uu->yy[jj-1][ii+1])*para.idx;
		dhdx   = ( -lv[jj  ][ii  ] +  lv[jj  ][ii+1])*para.idx;
		tq_vis+= vely * dhdx * (dudy + dvdx) * para.mu * 0.25;
		fy_vis+=        dhdx * (dudy + dvdx) * para.mu * 0.25;
		// 3
		vely   = ( vv->yy[jj  ][ii-1] + vv->yy[jj  ][ii  ])*0.5;
		dudy   = (-uu->xx[jj  ][ii-1] + uu->xx[jj+1][ii-1])*para.idy;
		dvdx   = (-uu->yy[jj  ][ii-1] + uu->yy[jj  ][ii  ])*para.idx;
		dhdx   = ( -lv[jj  ][ii-1] +  lv[jj  ][ii  ])*para.idx;
		tq_vis+= vely * dhdx * (dudy + dvdx) * para.mu * 0.25;
		fy_vis+=        dhdx * (dudy + dvdx) * para.mu * 0.25;
		// 4
		vely   = ( vv->yy[jj  ][ii  ] + vv->yy[jj  ][ii+1])*0.5;
		dudy   = (-uu->xx[jj  ][ii  ] + uu->xx[jj+1][ii  ])*para.idy;
		dvdx   = (-uu->yy[jj  ][ii  ] + uu->yy[jj  ][ii+1])*para.idx;
		dhdx   = ( -lv[jj  ][ii  ] +  lv[jj  ][ii+1])*para.idx;
		tq_vis+= vely * dhdx * (dudy + dvdx) * para.mu * 0.25;
		fy_vis+=        dhdx * (dudy + dvdx) * para.mu * 0.25;
		
			   
		atqp += tq_prs*(para.dx*para.dy) * dl[jj][ii];
		afxp += fx_prs*(para.dx*para.dy) * dl[jj][ii];
		afyp += fy_prs*(para.dx*para.dy) * dl[jj][ii];
		
		atqv += tq_vis*(para.dx*para.dy) * dl[jj][ii];
		afxv += fx_vis*(para.dx*para.dy) * dl[jj][ii];
		afyv += fy_vis*(para.dx*para.dy) * dl[jj][ii];
		
	}
	}
	
	ans[0] = atqp;
	ans[1] = afxp;
	ans[2] = afyp;
	ans[3] = atqv;
	ans[4] = afxv;
	ans[5] = afyv;

}

