#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<assert.h>
#include<gsl/gsl_sf_erf.h>
#include<gsl/gsl_interp.h>

int binary_search(int n, double* x, double z){
	assert(z>=x[0] && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int m=(i+j)/2;
		if(z>x[m]) i=m;
		else j=m;
	}
	return i;
}

typedef struct {int n; double *x,*y,*b,*c,*d;} cubic_spline;

cubic_spline* cubic_spline_alloc(int n, double *x, double *y){
	cubic_spline* s = (cubic_spline*)malloc(sizeof(cubic_spline));
	s->x = (double*)malloc(n*sizeof(double));
	s->y = (double*)malloc(n*sizeof(double));
	s->b = (double*)malloc(n*sizeof(double));
	s->c = (double*)malloc((n-1)*sizeof(double));
	s->d = (double*)malloc((n-1)*sizeof(double));
	s->n = n;
	for(int i=0; i<n; i++){
		s->x[i]=x[i];
		s->y[i]=y[i];
	}
	double h[n-1], p[n-1];
	for(int i=0; i<n-1; i++){
		h[i]=x[i+1]-x[i];
		assert(h[i]>0);
	}
	for(int i=0; i<n-1; i++){
		p[i]=(y[i+1]-y[i])/h[i];
	}
	double D[n], Q[n-1], B[n];
	D[0]=2;
	for(int i=0; i<n-2; i++){
		D[i+1]=2*h[i]/h[i+1]+2;
	}
	D[n-1]=2;
	Q[0]=1;
	for(int i=0; i<n-2; i++){
		Q[i+1]=h[i]/h[i+1];
	}
	for(int i=0; i<n-2; i++){
		B[i+1]=3*(p[i]+p[i+1]*h[i]/h[i+1]);
	}
	B[0]=3*p[0]; B[n-1]=3*p[n-2];
	for(int i=1; i<n; i++){
		D[i]-=Q[i-1]/D[i-1];
		B[i]-=B[i-1]/D[i-1];
	}
	s->b[n-1]=B[n-1]/D[n-1];
	for(int i=n-2; i>=0; i--){
		s->b[i]=(B[i]-Q[i]*s->b[i+1])/D[i];
	}
	for(int i=0; i<n-1; i++){
		s->c[i]=(-2*s->b[i]-s->b[i+1]+3*p[i])/h[i];
		s->d[i]=(s->b[i]+s->b[i+1]-2*p[i])/h[i]/h[i];
	}
	return s;
}

double cubic_spline_eval(cubic_spline *s, double z){
	int i = binary_search(s->n,s->x,z);
	double h=z-s->x[i];
	return s->y[i]+h*(s->b[i]+h*(s->c[i]+h*s->d[i]));
}

void cubic_spline_free(cubic_spline *s){
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s->d);
	free(s);
}

double cinterp_integ(cubic_spline *s, double z){ //Integral of cubic spline
	int i = binary_search(s->n,s->x,z);
	double integral=0;
	for(int k=1;k<=i;k++){
		integral+=(s->x[k]-s->x[k-1])*(s->y[k-1]+(0.5*s->b[k-1]*(s->x[k]-s->x[k-1])+(s->x[k]-s->x[k-1])*(s->c[k-1]*(s->x[k]-s->x[k-1])/3+(s->x[k]-s->x[k-1])*(0.25*s->d[k-1]*(s->x[k]-s->x[k-1])))));
	}
	integral+=(z-s->x[i])*(s->y[i]+(0.5*s->b[i]*(z-s->x[i])+(z-s->x[i])*(s->c[i]*(z-s->x[i])/3+(z-s->x[i])*(0.25*s->d[i]*(z-s->x[i])))));
	return integral;
}

double cspline_der(cubic_spline *s, double z){ //Derivative of cubic spline
	int i = binary_search(s->n,s->x,z);
	double derivative = s->b[i]+(z-s->x[i])*(2*s->c[i]+3*s->d[i]*(z-s->x[i]));
	return derivative;
}

double cinterp_gsl(cubic_spline *s, double z){ //Cubic interpolation with GSL
	gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_cspline,s->n);
	gsl_interp_init(interpolation,s->x,s->y,s->n);
	gsl_interp_accel*accelerator = gsl_interp_accel_alloc();
	return gsl_interp_eval(interpolation,s->x,s->y,z,accelerator);
}

void cinterp_integ_gsl(cubic_spline *s, double z, double *integral, double *derivative){ // Integral and derivative of cubic interpolation with GSL
	gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_cspline,s->n);
	gsl_interp_init(interpolation,s->x,s->y,s->n);
	gsl_interp_accel*accelerator = gsl_interp_accel_alloc();
	*integral = gsl_interp_eval_integ(interpolation,s->x,s->y,0,z,accelerator);
	*derivative = gsl_interp_eval_deriv(interpolation,s->x,s->y,z,accelerator);
	gsl_interp_free(interpolation);
}

int main(){
FILE* my_os=fopen("out.cubic.data.txt","w");

int num_data = 6;
double x[num_data], y[num_data];

for(int i=0;i<num_data;i++){
	//y[i]=10.0*rand()/RAND_MAX;
	//y[i]=cos(i);
	y[i]=gsl_sf_erf(i);
	x[i]=i;
	fprintf(my_os,"%10g %10g \n",x[i],y[i]);
}

int z=0;
double ex_points = 10; //numbers of interpolated data points between two tabulated data points
FILE* my_osn=fopen("out.cubic.intdata.txt","w");

cubic_spline* q=cubic_spline_alloc(num_data,x,y);
double integral=0; double derivative=0;
for(z=0;z<=ex_points*(num_data-1);z++){
	cinterp_integ_gsl(q,z/ex_points,&integral,&derivative);
	fprintf(my_osn,"%10g %10g %10g %10g %10g %10g %10g \n",z/ex_points,cubic_spline_eval(q,z/ex_points),cinterp_integ(q,z/ex_points),cspline_der(q,z/ex_points),cinterp_gsl(q,z/ex_points),integral,derivative);//,cinterp_integ_gsl(q,z/ex_points));
}
cubic_spline_free(q);

return 0;
}
