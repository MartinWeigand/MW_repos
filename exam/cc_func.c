#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <complex.h>
static double complex A,B;

static double complex F(double complex f(double complex), double t){
	return f((A+B)/2+(A-B)/2*cos(t))*sin(t)*(B-A)/2; //variable transformation
	}

double complex CC24(double complex f(double complex),double complex a, double complex b, double acc, double eps, double complex f2, double complex f3, int recursion){
	assert(recursion<1000000);

	double complex f1 = F(f,a+(b-a)/6);
	double complex f4 = F(f,a+5*(b-a)/6);
	double complex Q = (2*f1+1.0*f2+1.0*f3+2*f4)/6*(b-a);
	double complex q = (f1+f4+f2+f3)/4*(b-a);
	double tol = acc+eps*(fabs(creal(Q))+fabs(cimag(Q)));
	double err = fabs(creal(Q)-creal(q))/3.0+fabs(cimag(Q)-cimag(q))/3.0;
	if(err < tol){
		return Q;
	}
	else{
		double complex Q1 = CC24(f,a,(a+b)/2,acc/sqrt(2.0),eps,f1,f2,recursion+1);
		double complex Q2 = CC24(f,(a+b)/2,b,acc/sqrt(2.0),eps,f3,f4,recursion+1);
		return Q1+Q2;
	}
}

double complex clenshaw(double complex f(double complex), double complex a, double complex b, double acc,double eps){
	A = a;
	B = b;
	a = 0;
	b = M_PI;
	double complex f2=F(f,a+2*(b-a)/6);
	double complex f3=F(f,a+4*(b-a)/6);
	int recursion = 0;
	return CC24(f,a,b,2*acc,2*eps,f2,f3,recursion);
}
