double get_smd(double xx)
{
	double aa;
	double cc[20];

	cc[00] = +17.0/48.0+sqrt(3.0)*pi/108.0;
	cc[01] = +1.0/4.0;
	cc[02] = -1.0/4.0;
	cc[03] = +1.0/16.0;
	cc[04] = -sqrt(3.0)/12.0;
	cc[05] = +sqrt(3.0)/2.0;
	
	cc[10] = +55.0/48.0-sqrt(3.0)*pi/108.0;
	cc[11] = -13.0/12.0;
	cc[12] = +1.0/4.0;
	cc[13] = +1.0/48.0;
	cc[14] = +sqrt(3.0)/36.0;
	cc[15] = +sqrt(3.0)/2.0;
	
	
	if(fabs(xx) <= 1.0){
		aa = +cc[00]
			 +cc[01]*fabs(xx)
			 +cc[02]*xx*xx
			 +cc[03]*(1.0-2.0*fabs(xx))
			 *sqrt(-12.0*xx*xx+12.0*fabs(xx)+1.0)
			 +cc[04]
			 *asin(cc[05]*(2.0*fabs(xx)-1.0));
	}
	else if(fabs(xx) < 2.0){
		aa = +cc[10]
			 +cc[11]*fabs(xx)
			 +cc[12]*xx*xx
			 +cc[13]*(2.0*fabs(xx)-3.0)
			 *sqrt(-12.0*xx*xx+36.0*fabs(xx)-23.0)
			 +cc[14]
			 *asin(cc[15]*(2.0*fabs(xx)-3.0));
	}
	else{
		aa = 0.0;
	}
	
	// if(fabs(xx) <= 1.0){
		// aa = 0.5*(1.0+cos(pi*xx));
	// }
	// else{
		// aa = 0.0;
	// }
	
	// if(fabs(xx) <= 0.5){
		// aa = 1.0;
	// }
	// else{
		// aa = 0.0;
	// }

	return(aa);
}

