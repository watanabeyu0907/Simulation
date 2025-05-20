void calc_NSeq1_yy(coord *uu1, coord *uu2, double **ppp, coord *adv, double **phir, double **phim, double **phit2, double **phit3, coord *sft)
{
	int ii, jj;
	int ex, ey;
	coord0 adv1, adv2;
	double pd, ph, ro, gr;
	
	
	ex = para.nx+2;
	ey = para.ny+1;
	
	
	// ( 1-b ) y vel
	#ifdef _OPENMP
	#pragma omp parallel for \
	default(none) \
	private(ii, jj, adv1, adv2, pd, ph, ro, gr) \
	shared(para, ex, ey, uu1, uu2, ppp, adv, phit2, phit3, sft)
	#endif
	for(jj=2; jj<ey; jj++){
	for(ii=2; ii<ex; ii++){
				  
		pd = -1.0*para.idy*
					(ppp[jj+1][ii  ]
					-ppp[jj  ][ii  ]);
					
		ph = (+phit2[jj  ][ii  ]+phit2[jj+1][ii  ]
			  +phit3[jj  ][ii  ]+phit3[jj+1][ii  ])*0.25;
			  
		ro = 1.0 + ph*(para.ro2/para.ro1-1.0);
		gr = para.gg*ph;
		gr = para.gg;
		
		adv1.xx = +(uu1->xx[jj  ][ii-1]+uu1->xx[jj+1][ii-1])
				  *(uu1->yy[jj  ][ii-1]+uu1->yy[jj  ][ii  ])
				  *0.25*para.idx
				  -(uu1->xx[jj  ][ii  ]+uu1->xx[jj+1][ii  ])
				  *(uu1->yy[jj  ][ii  ]+uu1->yy[jj  ][ii+1])
				  *0.25*para.idx;
				  
		adv1.yy = +(uu1->yy[jj-1][ii  ]+uu1->yy[jj  ][ii  ])
				  *(uu1->yy[jj-1][ii  ]+uu1->yy[jj  ][ii  ])
				  *0.25*para.idy
				  -(uu1->yy[jj  ][ii  ]+uu1->yy[jj+1][ii  ])
				  *(uu1->yy[jj  ][ii  ]+uu1->yy[jj+1][ii  ])
				  *0.25*para.idy;		
		
		adv2.xx = +(uu2->xx[jj  ][ii-1]+uu2->xx[jj+1][ii-1])
				  *(uu2->yy[jj  ][ii-1]+uu2->yy[jj  ][ii  ])
				  *0.25*para.idx
				  -(uu2->xx[jj  ][ii  ]+uu2->xx[jj+1][ii  ])
				  *(uu2->yy[jj  ][ii  ]+uu2->yy[jj  ][ii+1])
				  *0.25*para.idx;
				  
		adv2.yy = +(uu2->yy[jj-1][ii  ]+uu2->yy[jj  ][ii  ])
				  *(uu2->yy[jj-1][ii  ]+uu2->yy[jj  ][ii  ])
				  *0.25*para.idy
				  -(uu2->yy[jj  ][ii  ]+uu2->yy[jj+1][ii  ])
				  *(uu2->yy[jj  ][ii  ]+uu2->yy[jj+1][ii  ])
				  *0.25*para.idy;				  
				  
		
		adv->yy[jj][ii]= -0.5*(adv1.xx+adv1.yy)
						 +1.5*(adv2.xx+adv2.yy)
						 +(+pd+sft->yy[jj][ii])/ro-gr;
						
	}
	}

}

