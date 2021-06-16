#include<stdio.h>
#include<math.h>
#include"ode.h"
#include<assert.h>
#define M1 1.0
#define M2 1.0
#define M3 1.0
#define GRAVITY 1

void f(int n, double x, double* y, double*dydx){
	dydx[0]=+y[1];
	dydx[1]=-y[0];
	}

void sir1(int n, double t, double* y, double*dydt){
	double Tc = 3;
	double Tr = 7;
	double N = 6000000-y[3];
	dydt[0] = -(y[1]*y[0])/(N*Tc);
	dydt[2] = y[1]/Tr;
	dydt[1] = -dydt[0]-dydt[2];
	}

void sir2(int n, double t, double* y, double*dydt){
	double Tc = 1;
	double Tr = 7;
	double N = 6000000-y[3];
	dydt[0] = -(y[1]*y[0])/(N*Tc);
	dydt[2] = y[1]/Tr;
	dydt[1] = -dydt[0]-dydt[2];
	}

void sir3(int n, double t, double* y, double*dydt){
	double Tc = 0.5;
	double Tr = 7;
	double N = 6000000-y[3];
	dydt[0] = -(y[1]*y[0])/(N*Tc);
	dydt[2] = y[1]/Tr;
	dydt[1] = -dydt[0]-dydt[2];
	}

void three_body(int n,double t,double*u,double*dudt){
	assert(n==12);
	double x1 = u[0], y1 = u[1];
	double x2 = u[2], y2 = u[3];
	double x3 = u[4], y3 = u[5];
	double vx1 = u[6], vy1 = u[7];
	double vx2 = u[8], vy2 = u[9];
	double vx3 = u[10], vy3 = u[11];
	double dx1dt = vx1, dy1dt = vy1;
	double dx2dt = vx2, dy2dt = vy2;
	double dx3dt = vx3, dy3dt = vy3;
	double dx12 = x2-x1, dy12 = y2-y1;
	double dx13 = x3-x1, dy13 = y3-y1;
	double dx23 = x3-x2, dy23 = y3-y2;
	double r12 = sqrt(dx12*dx12+dy12*dy12);
	double r13 = sqrt(dx13*dx13+dy13*dy13);
	double r23 = sqrt(dx23*dx23+dy23*dy23);
	double f12 = M1*M2*GRAVITY/r12/r12;
	double f13 = M1*M3*GRAVITY/r13/r13;
	double f23 = M2*M3*GRAVITY/r23/r23;
	double dvx1dt = 1./M1*( f12*dx12/r12+f13*dx13/r13);
	double dvy1dt = 1./M1*( f12*dy12/r12+f13*dy13/r13);
	double dvx2dt = 1./M2*(-f12*dx12/r12+f23*dx23/r23);
	double dvy2dt = 1./M2*(-f12*dy12/r12+f23*dy23/r23);
	double dvx3dt = 1./M3*(-f13*dx13/r13-f23*dx23/r23);
	double dvy3dt = 1./M3*(-f13*dy13/r13-f23*dy23/r23);
	dudt[0] = dx1dt; dudt[1] = dy1dt;
	dudt[2] = dx2dt; dudt[3] = dy2dt;
	dudt[4] = dx3dt; dudt[5] = dy3dt;
	dudt[6] = dvx1dt; dudt[7] = dvy1dt;
	dudt[8] = dvx2dt; dudt[9] = dvy2dt;
	dudt[10] = dvx3dt; dudt[11] = dvy3dt;
}

int main(int argc, char** argv){
	int kind = argc>1 ? atoi(argv[1]):1;
	if(kind==1){
		double a=0,b=2*M_PI,h=0.1,acc=1e-2,eps=1e-2;
		int n=2;
		double y[n];
		y[0]=0;
		y[1]=1;
		driver(n,f,a,y,b,h,acc,eps);
	}
	else if(kind==2){
		double a=0,b=150,h=0.01,acc=1e-4,eps=1e-4;
		int n=3;
		double y[n];
		y[2]=0;
		y[1]=10;
		y[0]=6000000;
		driver(n,sir1,a,y,b,h,acc,eps);
	}
	else if(kind==3){
		double a=0,b=150,h=0.01,acc=1e-4,eps=1e-4;
		int n=3;
		double y[n];
		y[2]=0;
		y[1]=10;
		y[0]=6000000;
		driver(n,sir2,a,y,b,h,acc,eps);
	}
	else if(kind==4){
		double a=0,b=150,h=0.01,acc=1e-4,eps=1e-4;
		int n=3;
		double y[n];
		y[2]=0;
		y[1]=10;
		y[0]=6000000;
		driver(n,sir3,a,y,b,h,acc,eps);
	}
	else{
		double u[] = {
			-0.97000436, 0.24308753,
 			0, 0,
			0.97000436, -0.24308753,
			0.4662036850, 0.4323657300,
			-0.93240737, -0.86473146,
			0.4662036850, 0.4323657300
		};
		int n = sizeof(u)/sizeof(u[0]);
		double a=0,b=6,h=0.1,acc=1e-3,eps=1e-3;
		driver(n,three_body,a,u,b,h,acc,eps);
	}
return 0;
}
