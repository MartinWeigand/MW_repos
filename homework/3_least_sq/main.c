#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include"functions.h"
#define FMT "%7.3f"

int main(){
	double x[]  = {1,  2,  3, 4, 6, 9,   10,  13,  15};
	double y[]  = {117,100,88,72,53,29.5,25.2,15.2,11.1};
	int n = sizeof(x)/sizeof(x[0]);
	double dy[n];
	for(int i=0; i<n; i++){
		dy[i] = 0.05*y[i];
	}
	for(int i=0; i<n; i++){
		dy[i] /= y[i];
		y[i] = log(y[i]);
	}
//	printf("%10g %10g %10g %10g %10g %10g",dy[0],dy[1],dy[2],dy[3],dy[4],dy[5]);
	gsl_vector* vx = gsl_vector_alloc(n);
	gsl_vector* vy = gsl_vector_alloc(n);
	gsl_vector* vdy = gsl_vector_alloc(n);
	for(int i=0; i<n; i++){
		gsl_vector_set(vx,i,x[i]);
		gsl_vector_set(vy,i,y[i]);
		gsl_vector_set(vdy,i,dy[i]);
	}

	int m=2;
	double funs(int i, double t){
   		switch(i){
   			case 0: return 1; break;
   			case 1: return t; break;
   			default: return NAN;
   		}
	}

	gsl_vector* c = gsl_vector_alloc(m);
	gsl_matrix* S = gsl_matrix_alloc(m,m);
	lsfit(m,funs,vx,vy,vdy,c,S);
//	gsl_vector_fprintf(stdout,c,FMT);
//	gsl_vector* dc = gsl_vector_alloc(m);
//	for(int k=0;k<m;k++){
//		double skk=gsl_matrix_get(S,k,k);
//		gsl_vector_set(dc,k,sqrt(skk));
//		}
//
//	double fit(double x){
//		double s=0;
//		for(int k=0;k<m;k++)s+=gsl_vector_get(c,k)*funs(k,x);
//		return s;
//		}
//
//	double fit_plus(int i, double x){
//		return fit(x)+gsl_vector_get(dc,i)*funs(i,x);
//		}
//
//	double fit_minus(int i, double x){
//		return fit(x)-gsl_vector_get(dc,i)*funs(i,x);
//		}
//
//	double c1 =gsl_vector_get(c, 1);
//	double dc1=gsl_vector_get(dc,1);
//	double T=-1/c1*log(2.0);
//	double dT=dc1/c1/c1;
//	printf("# half-life = %.3g +- %.2g days\n",T,dT);
//
//	printf("# time log(activity) delta(log(activity))\n");
//	for(int i=0;i<n;i++)printf("%g %g %g\n",x[i],y[i],dy[i]);
//	printf("\n\n");
//
//	for(int i=0;i<m;i++){
//		printf("# time fit fit_plus fit_minus; k=%i\n",i);
//		for(double z=x[0],dz=(x[n-1]-x[0])/64;z<=x[n-1];z+=dz)
//		printf("%g %g %g %g\n",z,fit(z),fit_plus(i,z),fit_minus(i,z));
//	printf("\n\n");
//	}
//
return 0;
}
