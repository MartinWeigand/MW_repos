#include<math.h>
#include<assert.h>
#include<stdio.h>

static double A;

static double F(double f(double),double t){
	return f(A+(1-t)/t)/t/t;
}

double wrap24( double f(double),double a, double b, double acc, double eps, double f2, double f3, int nrec){//wrapper calls F(f,t)
	assert(nrec<99);
	double f1=F(f,a+(b-a)/6);
	double f4=F(f,a+5*(b-a)/6);
	double Q=(2*f1+f2+f3+2*f4)/6*(b-a);
	double q=(f1+f4+f2+f3)/4*(b-a);
	double tolerance=acc+eps*fabs(Q);
	double error=fabs(Q-q);
	if(error < tolerance){
		return Q;
	}
	else{
		double Q1=wrap24(f,a,(a+b)/2,acc/sqrt(2.0),eps,f1,f2,nrec+1);
		double Q2=wrap24(f,(a+b)/2,b,acc/sqrt(2.0),eps,f3,f4,nrec+1);
		return Q1+Q2;
	}
}

double inf_integral(double f(double),double a,double acc,double eps ){
	A=a;
	a=0;
	double b=1;
	double f2=F(f,a+2*(b-a)/6);
	double f3=F(f,a+4*(b-a)/6);
	int nrec=0;
	return wrap24(f,a,b,2*acc,2*eps,f2,f3,nrec);
}
