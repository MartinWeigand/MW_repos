#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <complex.h>
static double complex A;

static double complex F(double complex f(double complex), double t){
	return f(A+(1-t)/t)/t/t; //variable transformation reducing an integral with an infinite limit to an integral with onlt finite limits
}

double complex wrap24(double complex f(double complex), double complex a, double complex b, double acc, double eps, double complex f2, double complex f3, int recursions){
	assert(recursions<1000000);

	double complex f1 = F(f,a+(b-a)/6);
	double complex f4 = F(f,a+5*(b-a)/6);
	double complex Q = (2*f1+1.0*f2+1.0*f3+2*f4)/6*(b-a);
	double complex q = (f1+f4+f2+f3)/4*(b-a);
	double tol = acc+eps*(fabs(creal(Q))+fabs(cimag(Q)));
	double err = fabs(creal(Q)-creal(q))+fabs(cimag(Q)-cimag(q));
	if(err < tol){
		return Q;
	}
	else{
		double complex Q1 = wrap24(f,a,(a+b)/2,acc/sqrt(2.0),eps,f1,f2,recursions+1);
		double complex Q2 = wrap24(f,(a+b)/2,b,acc/sqrt(2.0),eps,f3,f4,recursions+1);
		return Q1+Q2;
	}
}

double complex inf_integral(double complex f(double complex), double complex a, double acc, double eps){
	A = a;
	a = 0;
	double b = 1;
	double complex f2 = F(f,a+2*(b-a)/6);
	double complex f3 = F(f,a+4*(b-a)/6);
	int recursions = 0;
	return wrap24(f,a,b,2*acc,2*eps,f2,f3,recursions);
}
