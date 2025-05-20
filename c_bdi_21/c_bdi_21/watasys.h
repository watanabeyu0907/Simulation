// --------------------------------------------------
// 2021/02/01 Tomoaki Watamura @ Osaka UNIV.
// --------------------------------------------------

#ifndef WATASYS_H
#define WATASYS_H

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <math.h>

# include <omp.h>
# include <mkl.h>

#endif

# define MAX(A, B) ((A)>(B) ? (A):(B))
# define MIN(A, B) ((A)<(B) ? (A):(B))

# define TRUE	1
# define FALSE	!TRUE


# define array1D(array, xx, type){				\
	array = (type *)calloc(xx, sizeof(type));	\
}

# define free1D(array, xx){		\
	free(array);				\
}

# define array2D(array, yy, xx, type){					\
	int jj;												\
	array = (type **)calloc(yy, sizeof(type *));		\
	for(jj=0; jj<yy; jj++){								\
		array[jj] = (type *)calloc(xx, sizeof(type));	\
	}													\
}

# define free2D(array, yy, xx){		\
	int jj;							\
	for(jj=0; jj<yy; jj++){			\
		free(array[jj]);			\
	}								\
	free(array);					\
}

# define array3D(array, zz, yy, xx, type){						\
	int jj, kk;													\
	array = (type ***)calloc(zz, sizeof(type **));				\
	for(kk=0; kk<zz; kk++){										\
		array[kk] = (type **)calloc(yy, sizeof(type *));		\
		for(jj=0; jj<yy; jj++){									\
			array[kk][jj] = (type *)calloc(xx, sizeof(type));	\
		}														\
	}															\
}

# define free3D(array, zz, yy, xx){	\
	int jj, kk;						\
	for(kk=0; kk<zz; kk++){			\
		for(jj=0; jj<yy; jj++){		\
			free(array[kk][jj]);	\
		}							\
		free(array[kk]);			\
	}								\
	free(array);					\
}


typedef struct{
	double lx;
	double ly;
	double lz;
	double wu;
	double ro;
	double mu;
	
	double Cr;
	double Ep;
	double Re;
	int nx;
	int ny;
	int nz;
	double dx;
	double dy;
	double dz;
	double cT;
	double dT;
	
	double og;
	double ar;
	double Ri;
	double Ro;
	
	double idx;
	double idy;
	double idz;
	double idT;
	
}spec;

typedef struct{
	double **xx;
	double **yy;
}coord;

typedef struct{
	double xx;
	double yy;
}coord0;

typedef struct{
	double **xm;
	double **xp;
	double **ym;
	double **yp;
	double **cc;
	double **ic;
}pois_fnc;

typedef struct{
	double **r0;
	double **rr;
	double **pp;
	double **ee;
	double **ap;
	double **ae;
}bicg_fnc;


extern spec para;

extern int st_num;
extern int fpsnum;
extern int hevnum;

extern double pi;