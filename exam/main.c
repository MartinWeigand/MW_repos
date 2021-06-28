#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "integrator_func.h"
#include "cc_func.h"
#include "infint_func.h"
#include <complex.h>

int calls; //# of function calls

// Functions to be used for testing of implementations

double complex f1(double complex x){
	calls++;
	return creal(x)+I*cimag(x);
}

double complex f2(double complex x){
	calls++;
	return pow(creal(x),2)+I*creal(x)*cimag(x);
}

double complex f3(double complex x){
	calls++;
	return exp(x);
}

double complex f4(double complex x){
	calls++;
	return 1/sqrt(x);
}

double complex f5(double complex x){
	calls++;
	return exp(-x);
}

double complex f6(double complex x){
	calls++;
	return exp(-pow(x,2));
}

int main(){
	double acc = 1e-6; //Absolute accuracy
	double eps = 1e-6; //Relative accuracy
	printf("The following integrations is performed with absolute and relative accuracy goal on 1e-6.\n\n");

	printf("--------------------------------------------------\n");

//First we test the implementation in integrator_func.c
	printf("Test of the implementation of the adaptive integrator for a line integral of a complex-calued function:\n\n");

	//Test for f1
	double complex a = CMPLX(0,0); //Lower limit
	double complex b = CMPLX(1,1); //Upper limit
	calls = 0;
	double complex Q = adapt(f1,a,b,acc,eps); //Calculated result of complex line integral
	double complex exact = CMPLX(0,1); //Exact result of complex line integral
	printf("Integration from 0+0i->1+i of x+yi gives:\n");
	printf("Calculated: Q = %f+%f*i\n",creal(Q),cimag(Q));
	printf("Exact: Q = %f+%f*i\n",creal(exact),cimag(exact));
	printf("Estimated error = %f+%f*i\n",acc+eps*fabs(creal(Q)),acc+eps*fabs(cimag(Q)));
	printf("Actual error = %f+%f*i\n",fabs(creal(Q)-creal(exact)),fabs(cimag(Q)-cimag(exact)));
	printf("Number of function calls: %d\n\n",calls);

	//test for f2
	a = CMPLX(1,1);
	b = CMPLX(2,4);
	calls = 0;
	Q = adapt(f2,a,b,acc,eps);
	exact = CMPLX(-9.666667,11);
	printf("Integration from 1+i->2+4i of x*x+i*x*y gives:\n");
	printf("Calculated: Q = %f+%f*i\n",creal(Q),cimag(Q));
	printf("Exact: Q = %f+%f*i\n",creal(exact),cimag(exact));
	printf("Estimated error = %f+%f*i\n",acc+eps*fabs(creal(Q)),acc+eps*fabs(cimag(Q)));
	printf("Actual error = %f+%f*i\n",fabs(creal(Q)-creal(exact)),fabs(cimag(Q)-cimag(exact)));
	printf("Number of function calls: %d\n\n",calls);

	//test for f3
	a = CMPLX(0,0);
	b = CMPLX(1,1);
	calls = 0;
	Q = adapt(f3,a,b,acc,eps);
	printf("Integration from 0+0i->1+i of exp(x+iy) gives:\n");
	printf("Calculated: Q = %f+%f*i\n",creal(Q),cimag(Q));
	printf("Exact: UNKNOWN but both real and imaginary result looks like e-1=1.718281828\n");
	printf("Estimated error = %f+%f*i\n",acc+eps*fabs(creal(Q)),acc+eps*fabs(cimag(Q)));
	printf("Number of function calls: %d\n\n",calls);

	printf("--------------------------------------------------\n");

//Second we test the implementations in cc_func.c
	printf("Test of the implementation of the Clenshaw-Curtis variable transformation for complex line integrals:\n\n");

	//Test for f4 without and with Clenshaw-Curtis
	a = CMPLX(0,0);
	b = CMPLX(1,1);
	calls = 0;
	Q = adapt(f4,a,b,acc,eps);
	printf("Integration from 0+0i->1+i of 1/sqrt(x+iy) WITHOUT CC variable transformation gives:\n");
	printf("Calculated: Q = %f+%f*i\n",creal(Q),cimag(Q));
	printf("Estimated error = %f+%f*i\n",acc+eps*fabs(creal(Q)),acc+eps*fabs(cimag(Q)));
	printf("Number of calls: %d\n\n",calls);

	calls = 0;
	Q = clenshaw(f4,a,b,acc,eps);
	printf("Integration from 0+0i->1+i of 1/sqrt(x+iy) WITH CC variable transformation gives:\n");
	printf("Calculated: Q = %f+%f*i\n",creal(Q),cimag(Q));
	printf("Estimated error = %f+%f*i\n",acc+eps*fabs(creal(Q)),acc+eps*fabs(cimag(Q)));
	printf("Number of calls: %d\n\n",calls);

	//Test for f5 with Clenshaw-Curtis
	a = CMPLX(0,0);
	b = CMPLX(1,1);
	calls = 0;
	Q = adapt(f5,a,b,acc,eps);
	printf("Integration from 0+0i->1+i of x+iy WITHOUT CC variable transformation gives:\n");
	printf("Calculated: Q = %f+%f*i\n",creal(Q),cimag(Q));
	printf("Estimated error = %f+%f*i\n",acc+eps*fabs(creal(Q)),acc+eps*fabs(cimag(Q)));
	printf("Number of calls: %d\n\n",calls);

	calls = 0;
	Q = clenshaw(f5,a,b,acc,eps);
	printf("Integration from 0+0i->1+i of x+iy WITHOUT CC variable transformation gives:\n");
	printf("Calculated: Q = %f+%f*i\n",creal(Q),cimag(Q));
	printf("Estimated error = %f+%f*i\n",acc+eps*fabs(creal(Q)),acc+eps*fabs(cimag(Q)));
	printf("Number of calls: %d\n\n",calls);

	printf("--------------------------------------------------\n");

//Test of the implementations in infint_func.c
	printf("Test of the implementation of the adaptive integrator that accepts infinite limits:\n\n");

	//test for f6
	a = CMPLX(0,0);
	calls = 0;
	Q = inf_integral(f5,a,acc,eps);
	exact = CMPLX(1,0);
	printf("Integration from 0+0i->inf of exp(-x) gives:\n");
	printf("Calculated: Q = %f+%f*i\n",creal(Q),cimag(Q));
	printf("Exact: Q = %f+%f*i\n",creal(exact),cimag(exact));
	printf("Estimated error = %f+%f*i\n",acc+eps*fabs(creal(Q)),acc+eps*fabs(cimag(Q)));
	printf("Actual error = %f+%f*i\n",fabs(creal(Q)-creal(exact)),fabs(cimag(Q)-cimag(exact)));
	printf("Number of calls: %d\n\n",calls);

	//test for f7
	a = CMPLX(0,0);
	calls = 0;
	Q = inf_integral(f6,a,acc,eps);
	exact = CMPLX(0.886227,0);
	printf("Integration from 0+0i->inf of exp(-x*x) gives:\n");
	printf("Calculated: Q = %f+%f*i\n",creal(Q),cimag(Q));
	printf("Exact: Q = sqrt(pi)/2+%f*i = %f+%f*i\n",cimag(exact),creal(exact),cimag(exact));
	printf("Estimated error = %f+%f*i\n",acc+eps*fabs(creal(Q)),acc+eps*fabs(cimag(Q)));
	printf("Actual error = %f+%f*i\n",fabs(creal(Q)-creal(exact)),fabs(cimag(Q)-cimag(exact)));
	printf("Number of calls: %d",calls);
}
