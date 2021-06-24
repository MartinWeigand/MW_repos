#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>

double complex adapt24(double complex f(double complex), double complex a, double complex b, double acc, double eps, double complex f2, double complex f3, int recursions){
	assert(recursions < 1000000); //break if err<tol haven't happen after 1e6 recursions

	double complex f1 = f(a+(b-a)/6.0);
	double complex f4 = f(a+(b-a)*5.0/6.0);
	double complex Q = (b-a)*(2.0*f1+1.0*f2+1.0*f3+2.0*f4)/6.0; //Trapezium rule
	double complex q = (b-a)*(f1+f2+f3+f4)/4.0; //Rectangle rule
	double tol = acc+eps*(fabs(creal(Q))+fabs(cimag(Q))); //Accuracy goal
	double err = fabs(creal(Q)-creal(q))+fabs(cimag(Q)-cimag(q)); //Difference between the higher and lowe order (Trapezium and Rectangle rule)

	if(err<tol){ //If difference between higher and lower order is smaller than tolerance; break
		return Q;
	}
	else{ //subdivide into half-intervals and perform procedure recursively with acc=acc/sqrt(2). This continues until err<tol
		double complex Q1 = adapt24(f,a,(a+b)/2.0,acc/sqrt(2.0),eps,f1,f2,recursions+1);
		double complex Q2 = adapt24(f,(a+b)/2.0,b,acc/sqrt(2.0),eps,f3,f4,recursions+1);
		return Q1+Q2;
	}
}

double complex adapt(double complex f(double complex), double complex a, double complex b, double acc, double eps){
	double complex f2 = f(a+2*(b-a)/6.0);
	double complex f3 = f(a+4*(b-a)/6.0);
	int recursions = 0;
	return adapt24(f,a,b,acc,eps,f2,f3,recursions);
}
