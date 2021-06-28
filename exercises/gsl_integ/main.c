#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double f (double x, void* params){
	double f = log(x)/sqrt(x);
	return f;
}

double integration(){
	gsl_function F;
	F.function=&f;
	int limit=999;
	gsl_integration_workspace* w;
	w = gsl_integration_workspace_alloc (limit);
	double lower=0,upper=1,acc=1e-6,eps=1e-6,result,error;
	gsl_integration_qags(&F,lower,upper,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}
int main(){
	printf("The integral of ln(x)/sqrt(x) from 0 to 1 is %g\n",integration());
return 0;
}
