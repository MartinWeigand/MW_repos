#include <math.h>
#include <assert.h>
#include <stdio.h>
#include"integrator_func.h"
#include"cc_func.h"
#include"infint_func.h"
#include <gsl/gsl_integration.h>

int calls;

double f(double x){
	calls++;
	return sqrt(x);
}

double fn(double x){
	calls++;
	return 4*sqrt(1-x*x);
}

double fn_gsl(double x, void* p){
	return fn(x);
}

double fm(double x){
	calls++;
	return exp(sqrt(x));
}

double fh(double x){
	calls++;
	return 1/sqrt(x);
}

double fg(double x){
	calls++;
	return log(x)/sqrt(x);
}

double fj(double x){
	calls++;
	return exp(-x);
}

double fj_gsl(double x, void* p){
	return fj(x);
}
int main(){
	//Test for sqrt(x)
	double a=0,b=1,acc=0.001,eps=0.001;
	calls = 0;
	double Q = adapt(f,a,b,acc,eps);
	double exact = 0.66666666667;
	printf("The following integrations is performed with absolute and relative accuracy goal on 0.001.\n\n");
	printf("--------------------------------------------------\n");
	printf("Task A:\n\n");
	printf("Integration from 0->1 of sqrt(x) gives:\n");
	printf("Calculated: Q = %g\n",Q);
	printf("Exact: Q = %g\n",exact);
	printf("Number of calls: %d\n",calls);
	printf("Estimated error = %g\n",acc+fabs(Q)*eps);
	printf("Actual error = %g\n\n",fabs(Q-exact));

	//test for 4*sqrt(1-x*x)
	a=0,b=1,acc=0.001,eps=0.001;
	calls = 0;
	Q = adapt(fn,a,b,acc,eps);
	exact = 3.14159265359;
	printf("Integration from 0->1 of 4*sqrt(1-x*x) gives:\n");
	printf("Calculated: Q = %g\n",Q);
	printf("Exact: Q = %g\n",exact);
	printf("Number of calls: %d\n",calls);
	printf("Estimated error = %g\n",acc+fabs(Q)*eps);
	printf("Actual error = %g\n\n",fabs(Q-exact));

	//test for exp(-sqrt(*x))
	a=0,b=1,acc=0.001,eps=0.001;
	calls = 0;
	Q = adapt(fm,a,b,acc,eps);
	exact = 2.00000000000;
	printf("Integration from 0->1 of exp(sqrt(x)) gives:\n");
	printf("Calculated: Q = %g\n",Q);
	printf("Exact: Q = %g\n",exact);
	printf("Number of calls: %d\n",calls);
	printf("Estimated error = %g\n",acc+fabs(Q)*eps);
	printf("Actual error = %g\n\n",fabs(Q-exact));

	printf("--------------------------------------------------\n");
	printf("Task B:\n\n");

	//Test for 1/sqrt(x) with Clenshaw-function
	a=0,b=1,acc=0.001,eps=0.001;
	calls = 0;
	Q = adapt(fh,a,b,acc,eps);
	exact = 2;
	printf("Integration from 0->1 of 1/sqrt(x) WITHOUT CC variable transformation gives:\n");
	printf("Calculated: Q = %g\n",Q);
	printf("Exact: Q = %g\n",exact);
	printf("Number of calls: %d\n",calls);
	printf("Estimated error = %g\n",acc+fabs(Q)*eps);
	printf("Actual error = %g\n",fabs(Q-exact));

	calls = 0;
	Q = clenshaw(fh,a,b,acc,eps);
	printf("Integration from 0->1 of 1/sqrt(x) WITH CC variable transformation gives:\n");
	printf("Calculated: Q = %g\n",Q);
	printf("Exact: Q = %g\n",exact);
	printf("Number of calls: %d\n",calls);
	printf("Estimated error = %g\n",acc+fabs(Q)*eps);
	printf("Actual error = %g\n\n",fabs(Q-exact));

	//Test for ln(x)(/sqrt(x) with Clenshaw-function
	a=0,b=1,acc=0.001,eps=0.001;
	calls = 0;
	Q = adapt(fg,a,b,acc,eps);
	exact = -4;
	printf("Integration from 0->1 of ln(x)/sqrt(x) WITHOUT CC variable transformation gives:\n");
	printf("Calculated: Q = %g\n",Q);
	printf("Exact: Q = %g\n",exact);
	printf("Number of calls: %d\n",calls);
	printf("Estimated error = %g\n",acc+fabs(Q)*eps);
	printf("Actual error = %g\n",fabs(Q-exact));

	calls = 0;
	Q = clenshaw(fg,a,b,acc,eps);
	printf("Integration from 0->1 of ln(x)/sqrt(x) WITH CC variable transformation gives:\n");
	printf("Calculated: Q = %g\n",Q);
	printf("Exact: Q = %g\n",exact);
	printf("Number of calls: %d\n",calls);
	printf("Estimated error = %g\n",acc+fabs(Q)*eps);
	printf("Actual error = %g\n\n",fabs(Q-exact));
	printf("It is seen that the amount of calls is reducded significantly when calculted with the Clenshaw-Curtis variable transformation for both of the two integrals (1/sqrt(x) and ln(x)/sqrt(x)) which both have integrable divergiencies at the end-points of the intervals.\n\n");

	printf("To see how the variable transformation affect the accuracy for an integral without such a divergencey we look at an integral we also calculated in Task a.\n");
	calls = 0;
	Q = clenshaw(fn,a,b,acc,eps);
	exact = 3.14159265359;
	printf("Integration from 0->1 of 4*sqrt(1-x*x) WITH CC variable transformation gives:\n");
	printf("Calculated: Q = %g\n",Q);
	printf("Exact: Q = %g\n",exact);
	printf("Number of calls: %d\n",calls);
	printf("Estimated error = %g\n",acc+fabs(Q)*eps);
	printf("Actual error = %g\n\n",fabs(Q-exact));
	printf("The variable transformation roughly doubled the actual error without reducing the amount of calls a lot.\n\n");

	printf("We now calculate an integral with the GSL implementation:\n");
	double r;
	double err;
	size_t size = 10000;
	gsl_function func;
	calls = 0;
	func.function = &fn_gsl;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(size);
	gsl_integration_qags(&func,a,b,acc,eps,size,w,&r,&err);
	printf("Integration from 0->1 of 1/sqrt(x) with the GSL implementation gives::\n");
	printf("Calculated: Q = %g\n",r);
	printf("Exact: Q = %g\n",exact);
	printf("Number of calls: %d\n",calls);
	printf("Estimated error = %g\n",acc+fabs(r)*eps);
	printf("Actual error = %g\n\n",fabs(r-exact));
	printf("The amount of calls is much closer to the amount of calls of our implementation WITH variable transformation than without. The actual error is the smallest with the GSL implementation\n\n");

	printf("--------------------------------------------------\n");
	printf("Task C:\n\n");

	//test for exp(-x)
	a=0,acc=0.001,eps=0.001;
	calls = 0;
	Q = inf_integral(fj,a,acc,eps);
	exact = 1.00000000000;
	printf("Integration from 0->inf of exp(-x) gives:\n");
	printf("Calculated: Q = %g\n",Q);
	printf("Exact: Q = %g\n",exact);
	printf("Number of calls: %d\n",calls);
	printf("Estimated error = %g\n",acc+fabs(Q)*eps);
	printf("Actual error = %g\n\n",fabs(Q-exact));

	calls = 0;
	func.function = &fj_gsl;
	gsl_integration_qagiu(&func,a,acc,eps,size,w,&r,&err);
	printf("Integration from 0->inf of exp(-x) with the GSL implementation gives:\n");
	printf("Calculated: Q = %g\n",r);
	printf("Exact: Q = %g\n",exact);
	printf("Number of calls: %d\n",calls);
	printf("Estimated error = %g\n",acc+fabs(r)*eps);
	printf("Actual error = %g\n\n",fabs(r-exact));
}
