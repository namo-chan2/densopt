/*
this is not wave function but density matrix

 dy[i]dt=sum(j=0,j=element)coef[i][j]yin[j]

 
 */ 
#define NRANSI
//#include "nrutil.h"
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>
#define NR_END 1
#define FREE_ARG char*
double complex *vector(long nl,long nh);
void free_vector(double complex *v,long nl, long nh);
void rk4(double complex y[], double complex dydx[], int element,int nstate, double x, double h, double complex yout[],double complex coef[nstate][nstate],
	void (*derivs)(double t, double complex yin[], double complex dydt[],int element,int nstate, double complex coef[nstate][nstate]))
{
	int i;
	double xh,hh,h6;
	double complex *dym,*dyt,*yt;
	
	dym=vector(1,element);
	dyt=vector(1,element);
	yt=vector(1,element);
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=0;i<=element-1;i++){ 
		yt[i]=y[i]+hh*dydx[i];
		(*derivs)(xh,yt,dyt,element,nstate,coef);
	}	
	for (i=0;i<=element-1;i++){ 	
		yt[i]=y[i]+hh*dyt[i];
		(*derivs)(xh,yt,dym,element,nstate,coef);
	}	
	for (i=0;i<=element-1;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt,element,nstate,coef);
	for (i=0;i<=element-1;i++){
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	}	

	free_vector(yt,1,element);
	free_vector(dyt,1,element);
	free_vector(dym,1,element);
}
#undef NRANSI

void nrerror(char error_text[]){
	fprintf(stderr,"Numericalrecipes run-time erroe..\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system\n");
	exit(1);
}	

double complex *vector(long nl,long nh){
	double complex *v;
	v=(double complex *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double complex)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}	

void free_vector(double complex *v,long nl, long nh){
	free((FREE_ARG) (v+nl-NR_END));
}	

