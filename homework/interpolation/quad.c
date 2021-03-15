#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>

typedef struct {int n; double* x, *y, *b, *c;} qspline;

qspline* qspline_alloc(int n, double* x, double* y){
	qspline *s=(qspline*)malloc(sizeof(qspline));
	s->b = (double*)malloc((n-1)*sizeof(double));
	s->c = (double*)malloc((n-1)*sizeof(double));
	s->x = (double*)malloc(n*sizeof(double));
	s->y = (double*)malloc(n*sizeof(double));
	s->n = n;
	for(int i=0; i<n; i++){
		s->x[i]=x[i];
		s->y[i]=y[i];
	}
	int i;
	double p[n-1], h[n-1];
	for(i=0; i<n-1; i++){
		h[i]=x[i+1]-x[i];
		p[i]=(y[i+1]-y[i])/h[i];
	}
	s->c[0]=0;
	for(i=0; i<n-2; i++){
		s->c[i+1]=(p[i+1]-p[i]-s->c[i]*h[i])/h[i+1];
	}
	s->c[n-2]/=2;
	for(i=n-3; i>=0; i--){
		s->c[i]=(p[i+1]-p[i]-s->c[i+1]*h[i+1])/h[i];
	}
	for(i=0; i<n-1; i++){
		s->b[i]=p[i]-s->c[i]*h[i];
	}
	return s;
}

double qspline_eval(qspline *s, double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int i=0, j=s->n-1;
	while(j-i>1){
		int m=(i+j)/2;
		if(z>s->x[m]) i=m;
		else j=m;
	}
	double h=z-s->x[i];
	return s->y[i]+h*(s->b[i]+h*s->c[i]);
}

void qspline_free(qspline *s){
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s);
}

double qinterp_integ(qspline *s, double z){
int i=0, j=s->n-1;
while(j-i>1){
	int mid=(i+j)/2;
	if(z>s->x[mid]) i=mid;
		else j=mid;
}
double integral=0;
for(int k=1;k<=i;k++){ //takes care of the integral from 0 and up to the datapoint just before z
	integral+=(s->x[k]-s->x[k-1])*(s->y[k-1]+(s->x[k]-s->x[k-1])*(0.5*s->b[k-1]+s->c[k-1]*(s->x[k]-s->x[k-1])/3));
	}

integral+=(z-s->x[i])*(s->y[i]+(z-s->x[i])*(0.5*s->b[i]+s->c[i]*(z-s->x[i])/3));
return integral;
}

double qspline_der(qspline *s, double z){
int i=0, j=s->n-1;
while(j-i>1){
	int mid=(i+j)/2;
	if(z>s->x[mid]) i=mid;
		else j=mid;
}
double derivative = s->b[i]+2*s->c[i]*(z-s->x[i]);
return derivative;
}

int main(){
FILE* my_os=fopen("out.quad.data.txt","w");

int num_data = 6;
double x[num_data], y[num_data];

for(int i=0;i<num_data;i++){
//y[i]=10.0*rand()/RAND_MAX;
y[i]=cos(i);
x[i]=i;
fprintf(my_os,"%10g %10g \n",x[i],y[i]);
}

int z=0;
double ex_points = 30; //numbers of interpolated data points between two tabulated data points
FILE* my_osn=fopen("out.quad.intdata.txt","w");

qspline* q=qspline_alloc(num_data,x,y);

for(z=0;z<=ex_points*(num_data-1);z++){
fprintf(my_osn,"%10g %10g %10g %10g \n",z/ex_points,qspline_eval(q,z/ex_points),qinterp_integ(q,z/ex_points),qspline_der(q,z/ex_points));
	}
qspline_free(q);

return 0;
}
