#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include"qr_func.h"
#include <stdio.h>
#include <gsl/gsl_linalg.h>

int lim = 10000;
double delta = sqrt(2.22045e-16);


int newton( void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double eps){
	int n = x->size;
	int steps = 0;
	gsl_matrix* J = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* z = gsl_vector_alloc(n);
	gsl_vector* fz = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_vector* dx = gsl_vector_alloc(n);
	while(steps < lim){
		f(x,fx);
		for(int j=0; j<n; j++){
			double xj = gsl_vector_get(x,j);
			gsl_vector_set(x,j,xj+delta);
			f(x,df);
			gsl_vector_sub(df,fx);
			for(int i=0; i<n; i++){
				gsl_matrix_set(J,i,j,gsl_vector_get(df,i)/delta);
			}
			gsl_vector_set(x,j,xj);
		}
		QR_dec(J,R);
		QR_solve(J,R,fx,dx);
		gsl_vector_scale(dx,-1);
		double s = 1;
		while(s > 1.0/64){
			gsl_vector_memcpy(z,x);
			gsl_vector_add(z,dx);
			f(z,fz);
			if(gsl_blas_dnrm2(fz) < (1-s/2)*gsl_blas_dnrm2(fx)){
				break;
			}
			s *= 0.5;
			gsl_vector_scale(dx,0.5);
		}
		gsl_vector_memcpy(x,z);
		gsl_vector_memcpy(fx,fz);
		if(gsl_blas_dnrm2(dx) < delta || gsl_blas_dnrm2(fx) < eps){
			break;
		}
		steps++;
	}
	gsl_matrix_free(J);
	gsl_vector_free(fx);
	gsl_vector_free(z);
	gsl_vector_free(fz);
	gsl_vector_free(df);
	gsl_vector_free(dx);
	return steps;
}
