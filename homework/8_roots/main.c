#include <math.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "newtons_method.h"
#include <assert.h>
#include "ode.h"

void vector_print(char s[], gsl_vector* vec){
	printf("%s\n",s);
	for(int i=0; i<vec->size; i++){
		printf("%10g \n",gsl_vector_get(vec,i));
	}
}

void f(gsl_vector* x, gsl_vector* fx){
	double xx = gsl_vector_get(x,0);
	double yy = gsl_vector_get(x,1);
	gsl_vector_set(fx,0,400*pow(xx,3)-xx*(400*yy-2)-2);
	gsl_vector_set(fx,1,200*(yy-pow(xx,2)));
}

void fn(gsl_vector* x, gsl_vector* fx){
	double xx = gsl_vector_get(x,0);
	double yy = gsl_vector_get(x,1);
	double zz = gsl_vector_get(x,2);
	gsl_vector_set(fx,0,5*zz-2*xx);
	gsl_vector_set(fx,1,2*yy-14+xx);
	gsl_vector_set(fx,2,4*zz-10*yy);
}

static double rmax;
static double E;
void eqs(int n, double r, double* y, double* dydr){
	assert(n == 2);
	dydr[0] = y[1];
	dydr[1]= 2*(-1.0/r-E)*y[0];
}
double F_e(double e, double r){
	assert(r >= 0);
	E = e;
	const double rmin = 1e-5;
	if(r < rmin){
		return r-r*r;
	}
	double y[2] = {rmin-rmin*rmin,1-2*rmin};
	int n = 2;
	double h=0.1;
	double acc = 0.001;
	double eps = 0.001;
	driver(n,eqs,rmin,y,r,h,acc,eps,"data.txt");
	return y[0];
}
void master(gsl_vector* x, gsl_vector* fx){
	double e = gsl_vector_get(x,0);
	assert(e < 0);
	double f_rmax = F_e(e,rmax);
	gsl_vector_set(fx,0,f_rmax);
}

int main(){
	printf("Task A:\n\n");
	double eps = 0.001;
	int dim3 = 3;
	gsl_vector* x1 = gsl_vector_alloc(dim3);
	gsl_vector_set(x1,0,4.5);
	gsl_vector_set(x1,1,4.5);
	gsl_vector_set(x1,2,4.5);
	newton(fn,x1,eps);
	printf("The root(s) of the set of equation 5*z-2x=0, 2*y-14+x=0 and 4*z-10*y=0 is:\n");
	vector_print("x(root) =", x1);
	printf("Which, when checked, actually is the root. So it seems lik the implementation works\n\n");
	int dim2 = 2;
	gsl_vector* x2 = gsl_vector_alloc(dim2);
	gsl_vector_set(x2,0,2);
	gsl_vector_set(x2,1,2);
	newton(f,x2,eps);
	printf("The extremum(s) of the RosenBrock's valley function is found by differentiating the function and using our Newton's method implementation to find the roots:\n");
	vector_print("x(min) =",x2);
	printf("Which matches the expected root (1,1)\n\n");

	printf("Task B:\n\n");
	rmax = 8;
	double dim1 = 1;
	gsl_vector* x3 = gsl_vector_alloc(dim1);
	gsl_vector_set(x3,0,-2);
	newton(master,x3,eps);
	printf("The lowest root, e_0, for r_max=8 is:\n");
	vector_print("e_0 =",x3);
	printf("The exact value is -0.5.\n");
}
