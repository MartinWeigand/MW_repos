#include <stdio.h>
#include"func.h"
#include <math.h>
#define R 0.8

double f(int dim, double* p){
	double x = p[0];
	double y = p[1];
	printf("%g %g\n",x,y);
	double r = x*x+y*y<R*R?1:0;
	return r;
	}

int main(){
	double a[] = {0,0};
	double b[] = {1,1};
	double acc = 1e-3;
	double eps = 0;
	int dim = sizeof(a)/sizeof(a[0]);
	double integ = stratified_sam(dim,f,a,b,acc,eps,0,0);
	double exact = M_PI*R*R/4;
	FILE* stream = fopen("results_taskC.txt","w");
	fprintf(stream,"Task C:\n\n");
	fprintf(stream,"Stratified integration of x*x+y*y<%g*%g (1 if true, otherwise 0): \n",R,R);
	fprintf(stream,"Integral = %g\n",integ);
	fprintf(stream,"Estimated error = %g\n",acc+fabs(integ)*eps);
	fprintf(stream,"Exact = %g\n",exact);
	fprintf(stream,"Actual error = %g\n\n",fabs(integ-exact));
	fprintf(stream,"This stratified sampling is illustrated in stra_sam.plot.png");
}
