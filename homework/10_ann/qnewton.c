#include <float.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

static const double delta = sqrt(2.22045e-16);

void num_gradient(double F(gsl_vector*), gsl_vector* x, gsl_vector* gradient){
	double fx=F(x);
	for(int i=0; i<x->size; i++){
		double dx;
		double xi = gsl_vector_get(x,i);
		if(fabs(xi) < sqrt(delta)){
			dx = delta;
		}
		else{
			dx = fabs(xi)*delta;
		}
		gsl_vector_set(x,i,xi+dx);
		gsl_vector_set(gradient,i,(F(x)-fx)/dx);
		gsl_vector_set(x,i,xi);
	}
}

int qnewton(double F(gsl_vector*), gsl_vector* x, double acc){
	int n = x->size;
	int steps = 0;
	double fx = F(x);
	double fz;
	double lambda;
	double sTg;
	double sTy;
	double uTy;
	double gamma;

	gsl_matrix * B = gsl_matrix_alloc(n,n);
	gsl_vector* gradient_F = gsl_vector_alloc(n);
	gsl_vector* dx = gsl_vector_alloc(n);
	gsl_vector* z = gsl_vector_alloc(n);
	gsl_vector* gz = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_vector* u = gsl_vector_alloc(n);

	gsl_matrix_set_identity(B);

	num_gradient(F,x,gradient_F);

	while(steps<1000){
		steps++;
		gsl_blas_dgemv(CblasNoTrans,-1,B,gradient_F,0,dx);
		if(gsl_blas_dnrm2(dx) < delta*gsl_blas_dnrm2(x)){
			break;
		}
		if(gsl_blas_dnrm2(gradient_F) < acc){
			break;
		}
		lambda = 1;
		while(1){
			gsl_vector_memcpy(z,x);
			gsl_vector_add(z,dx);
			fz = F(z);
			gsl_blas_ddot(dx,gradient_F,&sTg);
			if(fz < fx+0.01*sTg){
				break;
			}
			if(lambda < delta){
				gsl_matrix_set_identity(B);
				break;
			}
			lambda *= 0.5;
			gsl_vector_scale(dx,0.5);
		}
		num_gradient(F,z,gz);
		gsl_vector_memcpy(y,gz);
		gsl_blas_daxpy(-1,gradient_F,y);
		gsl_vector_memcpy(u,dx);
		gsl_blas_dgemv(CblasNoTrans,-1,B,y,1,u);
		gsl_blas_ddot(dx,y,&sTy);
		if(fabs(sTy)>1e-12){
			gsl_blas_ddot(u,y,&uTy);
			gamma = uTy/2/sTy;
			gsl_blas_daxpy(-gamma,dx,u);
			gsl_blas_dger(1.0/sTy,u,dx,B);
			gsl_blas_dger(1.0/sTy,dx,u,B);
		}
		gsl_vector_memcpy(x,z);
		gsl_vector_memcpy(gradient_F,gz);
		fx = fz;
	}
	gsl_matrix_free(B);
	gsl_vector_free(gradient_F);
	gsl_vector_free(dx);
	gsl_vector_free(z);
	gsl_vector_free(gz);
	gsl_vector_free(y);
	gsl_vector_free(u);
	return steps;
}

//Downhill simplex method

void simplex_update(int dim, double simplex[][dim], double* f_values, int* hi, int* lo, double* centroid){
	*hi = 0;
	*lo = 0;
	double lowest = f_values[0];
	double highest = f_values[0];
	for(int k=1; k<dim+1; k++){
		double next = f_values[k];
		if(highest < next){
			highest = next;*hi=k;
		}
		if(next < lowest){
			lowest = next;
			*lo = k;
		}
	}
	for(int i=0; i<dim; i++){
		double sum = 0;
		for(int k=0; k<dim+1; k++){
			if(k != *hi){
				sum += simplex[k][i];
			}
		}
		centroid[i] = sum/dim;
	}
}

void simplex_init(int dim, double (*fun)(double*), double simplex[][dim], double* f_values, int* hi, int* lo, double* centroid){
	for(int k=0; k<dim+1; k++){
		f_values[k] = fun(simplex[k]);
	}
	simplex_update(dim,simplex,f_values,hi,lo,centroid);
}

void reflection(int dim, double* highest, double* centroid, double* reflected){
	for(int i=0; i<dim; i++){
		reflected[i] = 2*centroid[i]-highest[i];
	}
}

void expansion(int dim, double* highest, double* centroid, double* expanded){
	for(int i=0; i<dim; i++){
		expanded[i] = 3*centroid[i]-2*highest[i];
	}
}

void contraction(int dim, double* highest, double* centroid, double* contracted){
	for(int i=0; i<dim; i++){
		contracted[i] = 0.5*centroid[i]+0.5*highest[i];
	}
}

void reduction(int dim, double simplex[][dim], int lo){
	for(int k=0; k<dim+1; k++){
		if(k != lo);
		for(int i=0; i<dim; i++){
			simplex[k][i] = 0.5*(simplex[k][i]+simplex[lo][i]);
		}
	}
}

double distance(int dim, double* a, double* b){
	double s=0;
	for(int i=0; i<dim; i++){
		s += pow(a[i]-b[i],2);
	}
	return sqrt(s);
}

double size(int dim, double simplex[][dim], int lo){
	double s=0;
	for(int k=1; k<dim+1; k++){
		double dist = distance(dim,simplex[lo],simplex[k]);
		if(dist > s){
			s = dist;
		}
	}
	return s;
}


int downhill(int dim, double F(double*), double* start, double* step,double simplex_size_goal){
	int hi;
	int lo;
	int steps = 0;
	double simplex[dim+1][dim];
	for(int i=0; i<dim+1; i++){
		for(int k=0; k<dim; k++){
			simplex[i][k] = start[k];
		}
	}
	for(int i=0; i<dim; i++){
		simplex[i][i] += step[i];
	}
	double centroid[dim];
	double F_value[dim+1];
	double p1[dim];
	double p2[dim];
	simplex_init(dim,F,simplex,F_value,&hi,&lo,centroid);
	while(size(dim,simplex,lo) > simplex_size_goal){
		simplex_update(dim,simplex,F_value,&hi,&lo,centroid);
		reflection(dim,simplex[hi],centroid,p1);
		double f_re = F(p1);
		if(f_re < F_value[lo]){
			expansion(dim,simplex[hi],centroid,p2);
			double f_ex = F(p2);
			if(f_ex < f_re){
				for(int k=0; k<dim; k++){
					simplex[hi][k] = p2[k];
				}
				F_value[hi] = f_ex;
			}
			else{
				for(int k=0; k<dim; k++){
					simplex[hi][k] = p1[k];
				}
				F_value[hi] = f_re;
			}
		}
		else{
			if(f_re < F_value[hi]){
				for(int k=0; k<dim; k++){
					simplex[hi][k] = p1[k];
				}
				F_value[hi] = f_re;
			}
			else{
				contraction(dim,simplex[hi],centroid,p1);
				double f_co = F(p1);
				if(f_co < F_value[hi]){
					for(int k=0; k<dim; k++){
						simplex[hi][k] = p1[k];
					}
					F_value[hi] = f_co;
				}
				else {
					reduction(dim,simplex,lo);
					simplex_init(dim,F,simplex,F_value,&hi,&lo,centroid);
				}
			}
		}
		steps++;
	}
	for(int k=0; k<dim; k++){
		start[k] = simplex[lo][k];
	}
	return steps;
}
