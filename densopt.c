#include <stdio.h>
#define NRANSI
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <complex.h>
#include <time.h>
#include "/home/nishikata/physics.h"
double complex *vector(long nl, long ndt);
void free_vector(double complex *v, long nl, long ndt);
void derivs(double t,double complex yin[],double complex dydt[],int element,int nstate,double complex coef[nstate][nstate]){
	int j,k,l;
	double complex matyin[nstate][nstate],matdydt[nstate][nstate];
	for(j=0;j<=element-1;j++){
		dydt[j]=0;
	}	
	for(j=0;j<=nstate-1;j++){
		for(k=0;k<=nstate-1;k++){
			matyin[j][k]=yin[nstate*j+k];
		}
	}			
	for(j=0;j<=nstate-1;j++){
		for(k=0;k<=nstate-1;k++){
			matdydt[j][k]=0;
			for(l=0;l<=nstate-1;l++){
				matdydt[j][k] += coef[j][l]*matyin[l][k]-matyin[j][l]*coef[l][k];
			}
		}
	}
	for(j=0;j<=nstate-1;j++){
		for(k=0;k<=nstate-1;k++){
			dydt[nstate*j+k]=matdydt[j][k];
			yin[nstate*j+k]=matyin[j][k];
		}
	}			
}
void rk4(double complex yin[], double complex dydt[], int element,int nstate, double t, double dt, double complex yout[],
	double complex coef[nstate][nstate], 
	void (*derivs)(double t,double complex yin[],double complex dydt[],int element,int nstate,double complex coef[nstate][nstate]));
void spinmatrix(int nspin,int nstate, double complex sperator[nspin][3][nstate][nstate],double complex spin[nspin+1][3]);
void structure(int nspin, double exconst[nspin][nspin]);
void coefham(int j,int nspin,int element,int nstate,double complex coef[nstate][nstate],double complex coef1[nstate][nstate]
	,double complex sperator[nspin][3][nstate][nstate],double exconst[nspin][nspin],double mRyd,double Bmf[3]);

int main(void){
//initial condition 
//[scan parameter
	int i,j,k,l,m,n,o;
	int nspin,mgrid,element,nstate;
	double fintime,Apena;
	FILE *fp1;
	if((fp1=fopen("para_quantum.dat","r"))==NULL){
		fprintf(stderr,"ERROR1\n");
		return EXIT_FAILURE;
	}
	
	fscanf(fp1,"%*s %d",&nspin);
	fscanf(fp1,"%*s %d",&mgrid);
	fscanf(fp1,"%*s %le",&fintime);
	fscanf(fp1,"%*s %le",&Apena);
	printf("n:%d\n",nspin);
	printf("m:%d\n",mgrid);
	printf("t:%1.15e\n",fintime);
	printf("A:%1.15e\n",Apena);
	
	nstate=pow(2,nspin);
	printf("nstate:%d\n",nstate);

	double omega[nstate],exconst[nspin][nspin];
	for(j=0;j<=nspin-1;j++){
		for(k=0;k<=nspin-1;k++){
			exconst[j][k]=0;
		}
	}		

//scan parameter]
//	element=nstate;
	element=nstate*nstate;
	printf("element=%d\n",element);

	double complex efkeisan,denmat[nstate][nstate][mgrid+1],ximat[nstate][nstate][mgrid+1];
	double complex *yin,*dydt,*yout,coef[nstate][nstate],spin[nspin+1][3], consv,norm,coef1[nstate][nstate];
	double complex sperator[nspin][3][nstate][nstate];//S1x[e*e] S2x[e*e] S1y[e*e] S2y S1z S2z
	double t,dt,mf0[mgrid],mf[mgrid+1],functional,before_functional,penalty,efintegral,efd,calclationtime,expected,Apenap[mgrid+1];
	double  mfst, Bmf[3];
	time_t time0,time1;
	FILE *fp2,*fp3,*fp4,*fp5,*fp6,*fp7,*fp8,*fp9;
	char filepath[256];
	struct tm *local;
	dt=fintime/mgrid;
	yin=vector(1,element);
	dydt=vector(1,element);
	yout=vector(1,element);
	time0=time(NULL);
	local=localtime(&time0);

//	if((fp2=fopen("ef.dat","r"))==NULL){
//		fprintf(stderr,"ERROR2\n");
//		return EXIT_FAILURE;
//	}	
//	if((fp3=fopen("ef_debug.dat","w"))==NULL){
//		fprintf(stderr,"ERROR3\n");
//		return EXIT_FAILURE;
//	}	
	sprintf(filepath,"pop%d_1st.dat",nspin);
	if((fp4=fopen(filepath,"w"))==NULL){
		fprintf(stderr,"ERROR4\n");
		return EXIT_FAILURE;
	}	
	sprintf(filepath,"sp%d_1st.dat",nspin);
	if((fp3=fopen(filepath,"w"))==NULL){
		fprintf(stderr,"ERROR3\n");
		return EXIT_FAILURE;
	}	
	sprintf(filepath,"%02d%02dconvergence_check.dat",local->tm_hour,local->tm_min);
	if((fp7=fopen(filepath,"w"))==NULL){
		fprintf(stderr,"ERROR7\n");
		return EXIT_FAILURE;
	}	

	printf("////\nStart Time : %02d%02d\n////\n",local->tm_hour,local->tm_min);
//[reed magnetic field 
/*
	for(j=0;j<=mgrid-1;j++){
			fscanf(fp2,"%le %le",&t,&mf0[j]);
			if(j==1){
				printf("dt%1.15e\n",t-dt);
			}	
			mf0[j]=mf0[j];
			fprintf(fp3,"%1.15e %1.15e\n",t,mf0[j]); 
	}	 	

	fclose(fp2);
	fclose(fp3);
*/
//reed magnetic field]
//initial setting	
	for(j=0;j<=nstate-1;j++){
		for(k=0;k<=nstate-1;k++){
			denmat[j][k][0]=0;
		}
	}
	denmat[nstate-1][nstate-1][0]=1;
	
	for(j=0;j<=nstate-1;j++){
		for(k=0;k<=nstate-1;k++){
			yin[nstate*j+k]=denmat[j][k][0];
		}
	}	
	norm=0;
	for(j=0;j<=nstate-1;j++){
		norm += denmat[j][j][0];
	}	
	printf("%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e \n",t,creal(norm),creal(denmat[0][0][0]),creal(denmat[0][1][0]),creal(denmat[1][0][0]),creal(denmat[1][1][0]));
	
	printf("norm=%1.15e\n",creal(norm));

	t=0;
	fprintf(fp4,"%1.15e %1.15e ",t,1-creal(norm));
	for(j=0;j<=nstate-1;j++){
		fprintf(fp4,"%1.15e ",creal(denmat[j][j][0]));
	}	
	fprintf(fp4,"\n");

	double mfomega=2*M_PI*3.8e+12;	
	double mRyd=1.e-3*Rydberg;
	double sigma=1./sqrt(2*M_PI)*4.e-13;  //2ps
	double t0=7.e-13;
	efintegral=0;
	for(j=0;j<=mgrid;j++){
		t=j*dt;
		mf0[j]=10*exp(-(t-t0)*(t-t0)/(2*sigma*sigma))*cos(mfomega*t); //Bx
		efintegral += (mf0[j]*mf0[j]+mf0[j+1]*mf0[j+1])/2*dt;
		if(j<mgrid/72){
			Apenap[j] = Apena*sin(mfomega*j*dt);
		}else if(j>mgrid-mgrid/72){
			Apenap[j] = Apena*sin(mfomega*(fintime-j*dt));
		}else{
			Apenap[j] = Apena;
		}
	
	}	

	spinmatrix(nspin,nstate,sperator,spin); //spin matrix calclation
	structure(nspin,exconst); //13spin
	
	for(k=0;k<=nstate-1;k++){
		for(l=0;l<=nstate-1;l++){
			coef1[k][l]=0;
		}
	}		

	for(k=0;k<=nstate-1;k++){
		for(l=0;l<=nstate-1;l++){
			for(m=0;m<=3-1;m++){
				for(n=0;n<=nstate-1;n++){
					for(i=0;i<=nspin-1;i++){
						for(o=0;o<=nspin-1;o++){
							coef1[k][l] += -2*2*exconst[i][o]*mRyd*sperator[i][m][k][n]*sperator[o][m][n][l];//io=01
						}
					}	
				}
			}	
		}	
	}

	for(k=0;k<=nspin-1;k++){
		for(l=0;l<=3-1;l++){
			spin[k][l]=0;
			for(m=0;m<=nstate-1;m++){
				for(n=0;n<=nstate-1;n++){
					spin[k][l] += denmat[m][n][0]*sperator[k][l][n][m];
				}
			}
		}			
	}

	printf("%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\n"
			,t,creal(spin[0][0]),creal(spin[0][1]),creal(spin[0][2]),creal(spin[1][0]),creal(spin[1][1]),creal(spin[1][2]));
//[iteration 0
//dinamics calculation j:time	
	for(j=0;j<=mgrid-1;j++){
		t=dt*(j+1);

		Bmf[0]=mf0[j];
		Bmf[1]=0;
		//Bmf[0]=0;
		//Bmf[1]=mf0[j];
		Bmf[2]=mfomega*plank/(2*M_PI*1.65*bohr_mag);

		coefham(j,nspin,element,nstate,coef,coef1,sperator, exconst,mRyd,Bmf);
//[runge=kutta
		derivs(t,yin,dydt,element,nstate,coef);
		rk4(yin,dydt,element,nstate,t,dt,yout,coef,derivs);
		for(k=0;k<=nstate-1;k++){
			for(l=0;l<=nstate-1;l++){
				denmat[k][l][j+1]=yout[nstate*k+l];
				yin[nstate*k+l]=yout[nstate*k+l];
			}
		}
//runge=kutta]
		//if(j%500==0){	
			norm=0;
			for(k=0;k<=nstate-1;k++){
				norm += denmat[k][k][j+1];
			}
			fprintf(fp4,"%1.15e %1.15e ",dt*(j+1),1-creal(norm));
			for(k=0;k<=nstate-1;k++){
				fprintf(fp4,"%1.15e ",creal(denmat[k][k][j+1]));
			}	
			fprintf(fp4,"%1.15e\n",Bmf[1]);

			for(k=0;k<=nspin-1;k++){
				for(l=0;l<=3-1;l++){
					spin[k][l] = 0;
					for(m=0;m<=nstate-1;m++){
						for(n=0;n<=nstate-1;n++){
							spin[k][l] += denmat[m][n][j+1]*sperator[k][l][n][m]; // S1x S1y S1z S2x ...
						}	
					}
				}	
			}
			consv=0;
			for(k=0;k<=nspin-1;k++){
				for(l=0;l<=3-1;l++){
					consv += spin[k][l]*spin[k][l];
				}
			}
			fprintf(fp3,"%1.15e %1.15e ",t,1./4*nspin-consv);
			for(k=0;k<=nspin-1;k++){
				fprintf(fp3,"%1.15e %1.15e %1.15e "
				,creal(spin[k][0]),creal(spin[k][1]),creal(spin[k][2]));
			}	
			fprintf(fp3,"\n");
	
		//}
//spin expected	
	}	

	fclose(fp4);
	fclose(fp3);
	penalty=efintegral*2*M_PI/(plank*Apena);
	functional=0;
	for(j=0;j<=nstate-1;j++){
		for(k=0;k<=nstate-1;k++){
			for(l=0;l<=nspin-1;l++){
				expected += denmat[j][k][mgrid]*(sperator[l][2][k][j]);
			}
		}
	}		
	expected =expected/nspin; 
	expected += 0.5;
	functional=expected-penalty;
	before_functional=functional;
//iteration 0]	
	fprintf(fp7,"%d %1.15e %1.15e %1.15e %1.15e\n"
		,0,expected,penalty,functional,functional-before_functional);
	
	printf("Bz=%1.15e\n",Bmf[2]);

	printf("coefBz=%1.15e\n",1.65*bohr_mag/2*Bmf[2]);
	printf("coefBx=%1.15e\n",1.65*bohr_mag/2*10);

	printf("F=%1.15e  ",functional);
	printf("expected=%1.15e  ",expected);
	printf("Pn=%1.15e  ",penalty);
	printf("dF=%1.15e\n",functional-before_functional);
		
	double complex transmoment[nstate][nstate];

	for(j=0;j<=nstate-1;j++){
		for(k=0;k<=nstate-1;k++){
			transmoment[j][k]=0;
		}
	}		
	for(j=0;j<=nstate-1;j++){
		for(k=0;k<=nstate-1;k++){
			for(l=0;l<=nspin-1;l++){
				transmoment[j][k] += sperator[l][1][j][k]*1.65*bohr_mag;
			}
		}
	}		
//[iteration J	
// j:iteration

	for(j=0;;j++){
	//	t=dt*mgrid
		printf("j=%d\n",j+1);
		for(k=0;k<=nstate-1;k++){
			for(l=0;l<=nstate-1;l++){
				ximat[k][l][mgrid]=0;
			}
		}	

		for(k=0;k<=nstate-1;k++){
			for(l=0;l<=nstate-1;l++){
				for(m=0;m<=nstate-1;m++){
					for(n=0;n<=nspin-1;n++){
						//ximat[k][l][mgrid] += sperator[n][2][k][m]*denmat[m][l][mgrid];
						ximat[k][l][mgrid] += denmat[k][m][mgrid]*sperator[n][2][m][l];
					}
				}	
				ximat[k][l][mgrid] = ximat[k][l][mgrid]/nspin;
			}
		}		

		for(k=0;k<=nstate-1;k++){
			for(l=0;l<=nstate-1;l++){
				ximat[k][l][mgrid] += (0.5)*denmat[k][l][mgrid];
			}
		}		

		for(k=0;k<=nstate-1;k++){
			for(l=0;l<=nstate-1;l++){
				yin[nstate*k+l]=ximat[k][l][mgrid];
			}
		}	
		
		efkeisan=0;		
		for(k=0;k<=nstate-1;k++){
			for(l=0;l<=nstate-1;l++){
				for(m=0;m<=nstate-1;m++){
					efkeisan += ximat[k][l][mgrid]*transmoment[l][m]*denmat[m][k][mgrid];
				}
			}
		}
		efd=-Apenap[mgrid]*cimag(efkeisan);
// k:time
//[time reverce
		for(k=0;k<=mgrid-1;k++){
			t=dt*(mgrid-k-1);
			Bmf[0]=0;
			Bmf[1]=efd;
			Bmf[2]=mfomega*plank/(2*M_PI*1.65*bohr_mag);

			coefham(k,nspin,element,nstate,coef,coef1,sperator,exconst,mRyd,Bmf);
	
//[runge=kutta
			derivs(t,yin,dydt,element,nstate,coef);
			rk4(yin,dydt,element,nstate,t,-dt,yout,coef,derivs);
			norm=0;
			efkeisan=0;		

			for(l=0;l<=nstate-1;l++){
				for(m=0;m<=nstate-1;m++){
					ximat[l][m][mgrid-k-1]=yout[nstate*l+m];
					yin[nstate*l+m]=yout[nstate*l+m];
				}
			}
			for(l=0;l<=nstate-1;l++){
				for(m=0;m<=nstate-1;m++){
					for(n=0;n<=nstate-1;n++){
						efkeisan += ximat[l][m][mgrid-1-k]*transmoment[m][n]*denmat[n][l][mgrid-1-k];///csqrt(ximat[l][l][mgrid-1-k]*denmat[m][m][mgrid-1-k]);
					}
				}
			}
			
			efd = -Apenap[mgrid-1-k]*cimag(efkeisan);
		}
//time reverce]
		t=dt*0;
		for(k=0;k<=nstate-1;k++){
			for(l=0;l<=nstate-1;l++){
				yin[nstate*k+l]=denmat[k][l][0];
			}
		}	
		norm=0;
		efkeisan=0;		
		for(k=0;k<=nstate-1;k++){
			for(l=0;l<=nstate-1;l++){
				for(m=0;m<=nstate-1;m++){
					efkeisan += ximat[k][l][0]*transmoment[l][m]*denmat[m][k][0];///csqrt(ximat[k][k][0]*denmat[l][l][0]);
				}
			}	
		}
		mf[0]=-Apenap[0]*cimag(efkeisan);
		efintegral=0;
	
//[time forward
		for(k=0;k<=mgrid-1;k++){
			t=dt*(k+1);
			Bmf[0]=0;
			Bmf[1]=mf[k];
			Bmf[2]=mfomega*plank/(2*M_PI*1.65*bohr_mag);

			coefham(k,nspin,element,nstate,coef,coef1,sperator, exconst,mRyd,Bmf);
//[runge=kutta
			derivs(t,yin,dydt,element,nstate,coef);
			rk4(yin,dydt,element,nstate,t,dt,yout,coef,derivs);
			
			for(l=0;l<=nstate-1;l++){
				for(m=0;m<=nstate-1;m++){
					denmat[l][m][k+1]=yout[nstate*l+m];
					yin[nstate*l+m]=yout[nstate*l+m];
				}
			}
			efkeisan=0;		
			for(l=0;l<=nstate-1;l++){
				for(m=0;m<=nstate-1;m++){
					for(n=0;n<=nstate-1;n++){
						efkeisan += ximat[l][m][k+1]*transmoment[m][n]*denmat[n][l][k+1];///csqrt(ximat[l][l][k+1]*denmat[m][m][k+1]);
						//efkeisan += -ximat[l][m][k+1]*transmoment[m][n]*denmat[n][l][k+1];///csqrt(ximat[l][l][k+1]*denmat[m][m][k+1]);
					}
				}
			}

			mf[k+1]=-Apenap[k+1]*cimag(efkeisan);
			efintegral += (mf[k]*mf[k]+mf[k+1]*mf[k+1])/2*dt;
		}
//time forward]	

		if(j%50==0){
		//	sprintf(filepath,"%02d%02dpop%d.dat",local->tm_hour,local->tm_min,nspin);
			sprintf(filepath,"dpop%d.dat",nspin);
			if((fp5=fopen(filepath,"w"))==NULL){
				fprintf(stderr,"ERROR6\n");
				return EXIT_FAILURE;
			}	
			sprintf(filepath,"dsp%d.dat",nspin);
			if((fp8=fopen(filepath,"w"))==NULL){
				fprintf(stderr,"ERROR6\n");
				return EXIT_FAILURE;
			}	
			fprintf(fp5,"%d\n",j);
			fprintf(fp8,"%d\n",j);
			for(k=0;k<=mgrid;k++){
				//if(k%500==0){
					t=dt*k;
					norm=0;
					double complex overlap=0;
					for(l=0;l<=nstate-1;l++){
						norm += denmat[l][l][k];
						//for(m=0;m<=nstate-1;m++){
						//	overlap += denmat[l][m][k]*ximat[m][l][k];
						//}
					}
						
					fprintf(fp5,"%1.15e %1.15e ",t,1-creal(norm));
					//fprintf(fp5,"%1.15e ",csqrt(overlap));
					for(l=0;l<=nstate-1;l++){
						fprintf(fp5,"%1.15e ",creal(denmat[l][l][k]));
					}
					fprintf(fp5,"%1.15e %1.15e\n",Apenap[k],mf[k]);

					for(l=0;l<=nspin-1;l++){
						for(m=0;m<=3-1;m++){
							spin[l][m] = 0;
							for(n=0;n<=nstate-1;n++){
								for(o=0;o<=nstate-1;o++){
									spin[l][m] += denmat[o][n][k]*sperator[l][m][n][o]; // S1x S1y S1z S2x ...
								}	
							}
						}	
					}
					consv=0;
					for(l=0;l<=nspin-1;l++){
						for(m=0;m<=3-1;m++){
							consv += spin[l][m]*spin[l][m];
						}
					}
					fprintf(fp8,"%1.15e %1.15e ",t,1./4*nspin-consv);
					for(l=0;l<=nspin-1;l++){
						fprintf(fp8,"%1.15e %1.15e %1.15e "
						,creal(spin[l][0]),creal(spin[l][1]),creal(spin[l][2]));
					}	
					fprintf(fp8,"\n");
				//}
			}	

		fclose(fp5);
		fclose(fp8);
		
		}
		if(j%50==1){
		//	sprintf(filepath,"%02d%02dpop%d_sub.dat",local->tm_hour,local->tm_min,nspin);
			sprintf(filepath,"dpop%d_sub.dat",nspin);
			if((fp6=fopen(filepath,"w"))==NULL){
				fprintf(stderr,"ERROR8\n");
				return EXIT_FAILURE;
			}	
			sprintf(filepath,"dsp%d_sub.dat",nspin);
			if((fp9=fopen(filepath,"w"))==NULL){
				fprintf(stderr,"ERROR8\n");
				return EXIT_FAILURE;
			}	
			fprintf(fp6,"%d\n",j);
			fprintf(fp9,"%d\n",j);
			for(k=0;k<=mgrid;k++){
				if(k%500==0){
					t=dt*k;
					norm=0;
					for(l=0;l<=nstate-1;l++){
						norm += denmat[l][l][k];
					}	
					fprintf(fp6,"%1.15e %1.15e ",t,1.-creal(norm));
					for(l=0;l<=nstate-1;l++){
						fprintf(fp6,"%1.15e ",creal(denmat[l][l][k]));
					}
					fprintf(fp6,"%1.15e\n",mf[k]);

					for(l=0;l<=nspin-1;l++){
						for(m=0;m<=3-1;m++){
							spin[l][m] = 0;
							for(n=0;n<=nstate-1;n++){
								for(o=0;o<=nstate-1;o++){
									spin[l][m] += denmat[o][n][k]*sperator[l][m][o][n]; // S1x S1y S1z S2x ...
								}	
							}
						}	
					}
					consv=0;
					for(l=0;l<=nspin-1;l++){
						for(m=0;m<=3-1;m++){
							consv += spin[l][m]*spin[l][m];
						}
					}
					fprintf(fp8,"%1.15e %1.15e ",t,1./4*nspin-consv);
					for(l=0;l<=nspin-1;l++){
						fprintf(fp8,"%1.15e %1.15e %1.15e "
						,creal(spin[l][0]),creal(spin[l][1]),creal(spin[l][2]));
					}	
					fprintf(fp8,"\n");
				}	
			}	
		fclose(fp6);
		fclose(fp9);
		}
		printf("%1.15e\n",efintegral);
		penalty=efintegral*2*M_PI/(plank*Apena);
		functional=0;
		expected=0;
		for(k=0;k<=nstate-1;k++){
			for(l=0;l<=nstate-1;l++){
				for(m=0;m<=nspin-1;m++){
					expected += creal(denmat[k][l][mgrid]*(sperator[m][2][l][k]));
				}
			}
		}		
		expected=expected/nspin;
		expected += 0.5;
		functional= expected-penalty;
		printf("F=%1.15e  ",functional);
		printf("expected=%1.15e  ",expected);
		printf("Pn=%1.15e  ",penalty);
		printf("dF=%1.15e ",functional-before_functional);
		printf("Sx=%1.15e Sy=%1.15e Sz=%1.15e\n",creal(spin[0][0]),creal(spin[0][1]),creal(spin[0][2]));
		fprintf(fp7,"%d %1.15e %1.15e %1.15e %1.15e\n"
			,j+1,expected,penalty,functional,functional-before_functional);
		if((functional-before_functional)<1.0e-11){
			break;
		}	
		before_functional=functional;
	}
	fclose(fp7);

	free_vector(yout,1,element);
	free_vector(dydt,1,element);
	free_vector(yin,1,element);
	time1=time(NULL);
	calclationtime=difftime(time1,time0);
	printf("calclationtime-%lf\n",calclationtime);

	return 0;
}

#undef NRANSI

void spinmatrix(int nspin,int nstate, double complex sperator[nspin][3][nstate][nstate],double complex spin[nspin+1][3]){
	int j,k,l,m,gyou,retu;
	double shoukou[nspin][2][nstate][nstate]; // 0+ 1-
	double complex denmat[nstate][nstate];
	
	//
	for(j=0;j<=nspin-1;j++){ 
		for(k=0;k<=nstate-1;k++){ 
			for(l=0;l<=nstate-1;l++){ 
				sperator[j][0][k][l]=0;
				sperator[j][1][k][l]=0;
				sperator[j][2][k][l]=0;
				shoukou[j][0][k][l]=0;
				shoukou[j][1][k][l]=0;
			}
		}
	}

	//Sz
	for(j=0;j<=nspin-1;j++){  //2spin 0 1
		for(k=0;k<=pow(2,nspin-j-1)-1;k++){ 
			for(l=0;l<=pow(2,j)-1;l++){ 
				gyou=pow(2,j+1)*k+l;
				sperator[j][2][gyou][gyou]=1./2;
				retu=pow(2,j+1)*k+l+pow(2,j);
				sperator[j][2][retu][retu]=-1./2;
			}
		}
	}	
	
	//S+-
	for(j=0;j<=nspin-1;j++){  //2spin 0 1
		for(k=0;k<=pow(2,nspin-j-1)-1;k++){ 
			for(l=0;l<=pow(2,j)-1;l++){ 
				gyou=pow(2,j+1)*k+l;
				retu=pow(2,j+1)*k+l+pow(2,j);
				shoukou[j][0][gyou][retu]=1;
				shoukou[j][1][retu][gyou]=1;
			}
		}
	}	
	//Sxy
	for(j=0;j<=nspin-1;j++){  //2spin 0 1
		for(k=0;k<=nstate-1;k++){ 
			for(l=0;l<=nstate-1;l++){ 
				sperator[j][0][k][l]=(shoukou[j][0][k][l]+shoukou[j][1][k][l])/2;
				sperator[j][1][k][l]=-I*(shoukou[j][0][k][l]-shoukou[j][1][k][l])/2;
			}
		}
	}		

	for(j=0;j<=nspin-1;j++){
		for(k=0;k<=3-1;k++){
			spin[j][k] = 0;
			for(l=0;l<=nspin-1;l++){
				for(m=0;m<=nspin-1;m++){
					spin[j][k]+=denmat[l][m]*sperator[j][k][m][l];
				}	
			}
		}	
	}

/*	
	printf("\n//S0z//\n");
	for(j=0;j<nstate;j++){
		for(k=0;k<nstate;k++){
			printf("%2e ",creal(sperator[0][2][j][k]));
		}
		printf("\n");
	}	
	printf("\n//S2z//\n");
	for(j=0;j<nstate;j++){
		for(k=0;k<nstate;k++){
			printf("%2e ",cimag(sperator[1][1][j][k]));
		}
		printf("\n");
	}	
*/
}

void structure(int nspin, double exconst[nspin][nspin]){
	exconst[0][0]=0;
	exconst[0][1]=1.116398;
	exconst[0][2]=0;
	exconst[1][0]=1.116398;
	exconst[1][1]=0;
	exconst[1][2]=1.116398;
	exconst[2][0]=0;
	exconst[2][1]=1.116398;
	exconst[2][2]=0;
/*
	int j,k;
	for(j=1;j<=13-1;j++){                                                                                                                               
		exconst[0][j]=1.116398;
	}
	exconst[1][2]=1.116398;
	exconst[1][4]=1.116398;
	exconst[1][5]=1.116398;
	exconst[1][8]=1.116398;

	exconst[2][3]=1.116398;
	exconst[2][5]=1.116398;
	exconst[2][6]=1.116398;

	exconst[3][4]=1.116398;
	exconst[3][6]=1.116398;
	exconst[3][7]=1.116398;

	exconst[4][7]=1.116398;
	exconst[4][8]=1.116398;

	exconst[5][9]=1.116398;
	exconst[5][10]=1.116398;

	exconst[6][10]=1.116398;
	exconst[6][11]=1.116398;

	exconst[7][11]=1.116398;
	exconst[7][12]=1.116398;

	exconst[8][9]=1.116398;
	exconst[8][12]=1.116398;

	exconst[9][10]=1.116398;
	exconst[9][12]=1.116398;

	exconst[10][11]=1.116398;
	exconst[11][12]=1.116398;

	for(j=0;j<=13-1;j++){
		for(k=0;k<=13-1;k++){
			exconst[k][j]=exconst[j][k];
		}
	}
*/
/*
	int j,k,l,x,y,z,xf,yf,zf,xs,ys,zs,check,nshell,sitef,sites;
	nshell=1;
	int site[2*nshell+1][2*nshell+1][2+nshell+1];
	printf("hoge\n");
	for(j=0;j<=2*nshell+1;j++){
		for(k=0;k<=2*nshell+1;k++){
			for(l=0;l<=2*nshell+1;l++){
				site[j][k][l]=0;
			}
		}
	}			
	sitef=-1;
	for(xf=-nshell;xf<=nshell;xf++){
		for(yf=-nshell;yf<=nshell;yf++){
			for(zf=-nshell;zf<=nshell;zf++){
				if((xf+yf+zf)%2==0){
					check=xf*xf+yf*yf+zf*zf;
					if(check<=2*nshell){
						x=xf+nshell;			
						y=yf+nshell;			
						z=zf+nshell;			
						sitef++;
						site[x][y][z]=sitef;
						printf("%d: %d %d %d\n",sitef,xf,yf,zf);
					}
				}	
			}
		}
	}
						

	sitef=-1;
	printf("hoge\n");
	for(xf=-nshell;xf<=nshell;xf++){
		for(yf=-nshell;yf<=nshell;yf++){
			for(zf=-nshell;zf<=nshell;zf++){
				if((xf+yf+zf)%2==0){
					check=xf*xf+yf*yf+zf*zf;
					if(check<=2*nshell){
						for(j=0;j<=3-1;j++){
					
							for(k=-1;k<=1;k++){
								for(l=-1;l<=1;l++){
									if(j==0){
										xs=xf;
										ys=yf+k;
										zs=zf+l;
									}else if(j==1){
										xs=xf+l;
										ys=yf;
										zs=zf+k;
									}else{
										xs=xf+k;
										ys=yf+l;
										zs=zf;
									}
						
									if((xs+ys+zs)%2==0){
										check=xs*xs+ys*ys+zs*zs;
										if(check<=2*nshell){
											x=xs+nshell;			
											y=ys+nshell;			
											z=zs+nshell;			
											sites=site[x][y][z];
											printf("%d: %d %d %d\n",sites,xs,ys,zs);
											exconst[sitef][sites]=1.116398;
										}	

									}
								}	
							}
						}
					}
				}
			}
		}			
	}

*/
}

void coefham(int j,int nspin,int element,int nstate,double complex coef[nstate][nstate],double complex coef1[nstate][nstate],
double complex sperator[nspin][3][nstate][nstate],double exconst[nspin][nspin],double mRyd,double Bmf[3]){

	int i,o,k,l,m,n;

	for(k=0;k<=nstate-1;k++){
		for(l=0;l<=nstate-1;l++){
			coef[k][l]=0;
		}
	}	
		
//unaxial anisotropy

/*	for(k=0;k<=nstate-1;k++){
		for(l=0;l<=nstate-1;l++){
			for(m=0;m<=nspin-1;m++){
				coef[k][l] += -2*3*1.e-6*-1.6e-19*pow(sperator[m][0][k][l],2);
			}	
		}	
	}
*/
//unaxial anisotropy
//[calc hamil mf
	for(k=0;k<=nstate-1;k++){
		for(l=0;l<=nstate-1;l++){
			for(m=0;m<=nspin-1;m++){
				for(n=0;n<=3-1;n++){
					coef[k][l] += -1.65*bohr_mag*(Bmf[n]*sperator[m][n][k][l]);
				}
			}
				
		}	
	}
//calc hamil mf]

	for(k=0;k<=nstate-1;k++){
		for(l=0;l<=nstate-1;l++){
			coef[k][l] = (coef1[k][l]+coef[k][l])*(-I*2*M_PI/plank);
		}
	}		

}
