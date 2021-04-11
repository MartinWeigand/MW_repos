#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include"functions.h"

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

	gsl_vector* delta_c = gsl_vector_alloc(m);
	for(int k=0; k<m; k++){
		double S_kk = gsl_matrix_get(S,k,k);
		gsl_vector_set(delta_c,k,sqrt(S_kk));
	}
	double c1 = gsl_vector_get(c,1);
	double delta_c1 = gsl_vector_get(delta_c,1);
	double half_life = -log(2)/c1;
	double delta_half_life = delta_c1/c1/c1; //from propogation of errors
	printf("# The fit gives ln(y)=%g(+-%g)%g(+-%g)*t. This gives a half-life of ln(2)/%g=%g+-%g days. According to wikipedia the modern value is 3,6319 days. Hence the half-life doesn't agree with the moderne value within the estimated uncertainty.\n",gsl_vector_get(c,0),gsl_vector_get(delta_c,0),c1,delta_c1,-c1,half_life,delta_half_life);
	printf("\n");
	printf("# time ln(activity) uncertainty\n");
	for(int i=0; i<n; i++){
		printf("%10g %10g %10g\n",x[i],y[i],dy[i]);
	}
	printf("\n\n");
	double fit(double x){
		double s=0;
		for(int k=0; k<m; k++){
			s+=gsl_vector_get(c,k)*funs(k,x);
		}
		return s;
	}
	printf("# time fit\n");
	for(double z=x[0],dz=(x[n-1]-x[0])/64;z<=x[n-1];z+=dz){
		printf("%g %g\n",z,fit(z));
	}

	double fit_plus(int i, double x){
		return fit(x)+gsl_vector_get(delta_c,i)*funs(i,x);
	}

	double fit_minus(int i, double x){
		return fit(x)-gsl_vector_get(delta_c,i)*funs(i,x);
	}
	printf("\n\n");
	for(int i=0;i<m;i++){
		printf("# time fit_plus fit_minus; k=%i\n",i);
		for(double z=x[0], dz=(x[n-1]-x[0])/64; z<=x[n-1]; z+=dz){
			printf("%g %g %g\n",z,fit_plus(i,z),fit_minus(i,z));
		}
	printf("\n\n");
	}

return 0;
}
