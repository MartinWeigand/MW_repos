#include<math.h>
#include<stdio.h>

void rkstep12(
	int n,
	void f(int n, double x, double y[], double dydx[]),
	double x,
	double y[],
	double h,
	double yh[],
	double yerr[]
){
	double k0[n];
	f(n,x,y,k0);
	double y1[n];
	for(int i=0; i<n; i++){
		y1[i] = y[i]+(0.5*h)*k0[i];
	}
	double k1[n];
	f(n,x+0.5*h,y1,k1);
	for(int i=0; i<n; i++){
		yh[i] = y[i]+h*k1[i];
	}
	for(int i=0; i<n; i++){
		yerr[i] = (k1[i]-k0[i])*h;
	}
}
// #define PRINT(x) fprintf(stream,"%9.3g ",x)
//#define TRACE(x,y) PRINT(x); for(int i=0; i<n; i++){PRINT(y[i]);}; fprintf(stream,"\n")

void driver(
	int n,
	void f(int n, double x, double* y, double* dydx),
	double a,
	double* y,
	double b,
	double h,
	double acc,
	double eps,
	char file[]
){
	double x=a;
	FILE* stream = fopen(file,"w");
//	fprintf(stream,"%9.3g ",x);
//	for(int i=0; i<n; i++){
//		fprintf(stream,"%9.3g ",y[i]);
//	}
//	fprintf(stream,"\n");
	while(x<b){
		if(x+h>b){
			h=b-x;
		}
		double yh[n];
		double yerr[n];
		rkstep12(n,f,x,y,h,yh,yerr);
		double sum = 0;
		for(int i=0; i<n; i++){
			sum += y[i]*y[i];
		}
		double norm_y = sqrt(sum);
		sum = 0;
		for(int i=0; i<n; i++){
			sum += yerr[i]*yerr[i];
		}
		double err = sqrt(sum);
		double tol = (acc+eps*norm_y)*sqrt(h/(b-a));
		if(err<tol){
			x = x+h;
			for(int i=0; i<n; i++){
				y[i] = yh[i];
			}
			fprintf(stream,"%9.3g ",x);
			for(int i=0; i<n; i++){
				fprintf(stream,"%9.3g ",y[i]);
			}
			fprintf(stream,"\n");
		}
		if(err>0) h *= 0.95*pow(tol/err,0.25);
		else h *= 2;
	}
	fclose(stream);
}
