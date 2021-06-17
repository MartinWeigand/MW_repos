#include<math.h>
#include<stdio.h>
#include<assert.h>

double adapt24(double f(double), double a, double b, double acc, double eps, double f2, double f3, int nrec){
	assert(nrec < 1000000);
	double f1 = f(a+(b-a)/6.0);
	double f4 = f(a+(b-a)*5.0/6.0);
	double Q = (b-a)*(2.0*f1+1.0*f2+1.0*f3+2.0*f4)/6.0;
	double q = (b-a)*(f1+f2+f3+f4)/4.0;
	double tol = acc+eps*fabs(Q);
	double err = fabs(Q-q);
	if(err < tol){
		return Q;
	}
	else{
		double Q1 = adapt24(f,a,(a+b)/2.0,acc/sqrt(2.0),eps,f1,f2,nrec+1);
		double Q2 = adapt24(f,(a+b)/2.0,b,acc/sqrt(2.0),eps,f3,f4,nrec+1);
		return Q1+Q2;
	}
}

double adapt(double f(double), double a, double b, double acc, double eps){
	double f2 = f(a+2*(b-a)/6.0);
	double f3 = f(a+4*(b-a)/6.0);
	int nrec=0;
	return adapt24(f,a,b,acc,eps,f2,f3,nrec);
}
