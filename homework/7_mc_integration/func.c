#include <math.h>
#include <stdlib.h>
#define RND (double)rand()/RAND_MAX
#include <complex.h>
#include <assert.h>

complex plainmontecarlo(
	int dim,
	double f(int dim, double* x),
	double* a,
	double* b,
	int N
){
	double V = 1;
	for(int i=0; i<dim; i++){
		V *= b[i]-a[i];
	}
	double sum = 0;
	double sum2 = 0;
	double x[dim];
	for(int i=0; i<N; i++){
		for(int i=0; i<dim; i++){
			x[i] = a[i]+RND*(b[i]-a[i]);
		}
		double fx = f(dim,x);
		sum += fx;
		sum2 += fx*fx;
	}
	double mean = sum/N;
	double sigma = sqrt(sum2/N-mean*mean);
	complex result = mean*V+I*sigma*V/sqrt(N);
	return result;
}

double corput(int n, int base){
	double q=0;
	double bk = (double) 1/base;
	while(n>0){
		q += (n%base)*bk;
		n /= base;
		bk /= base;
	}
	return q;
}

void halton1(int n, int d, double *a, double *b, double *x){
	int base[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79};
	int dmax = sizeof(base)/sizeof(int);
	assert(d <= dmax);
	for(int i=0; i<d; i++){
		x[i] = a[i]+corput(n+1,base[i])*(b[i]-a[i]);
	}
}

void halton2(int n, int d, double *a, double *b, double *x){
	int base[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83};
	int dmax = sizeof(base)/sizeof(int);
	assert(d <= dmax);
	for(int i=0; i<d; i++){
		x[i] = a[i]+corput(n+1,base[i])*(b[i]-a[i]);
	}
}

complex quasimontecarlo(
	int d,
	double f(int ,double*),
	double a[],
	double b[],
	int N
){
	double x[d];
	double vol = 1;
	for(int i=0; i<d; i++){
		vol *= b[i]-a[i];
	}
	double sum1 = 0;
	double sum2 = 0;
	for(int i=0; i<N/2; i++){
		halton1(i,d,a,b,x);
		sum1 += f(d,x);
	}
	for(int i=0;i<N/2;i++){
		halton2(i,d,a,b,x);
		sum2 += f(d,x);
	}
	double integ = (sum1+sum2)/N*vol;
	double error = fabs(sum1-sum2)/N*vol;
	return integ+I*error;
}

double stratified_sam(
	int dim,
	double f(int dim, double* x),
	double* a,
	double* b,
	double acc,
	double eps,
	int n_reuse,
	double mean_reuse
){
	int N = 16*dim;
	double V = 1;
	for(int k=0; k<dim; k++){
		V *= b[k]-a[k];
	}
	int n_l[dim];
	double n_r[dim];
	double x[dim];
	double mean_l[dim];
	double mean_r[dim];
	double mean=0;
	for(int k=0; k<dim; k++){
		mean_l[k] = 0;
		mean_r[k] = 0;
		n_l[k] = 0;
		n_r[k] = 0;
	}
	for(int i=0; i<N; i++){
		for(int k=0; k<dim; k++){
			x[k] = a[k]+RND*(b[k]-a[k]);
		}
		double fx = f(dim,x);
		for(int k=0; k<dim; k++){
			if(x[k] > (a[k]+b[k])/2){
				n_r[k]++;
				mean_r[k] += fx;
			}
			else{
				n_l[k]++;
				mean_l[k] += fx;
			}
		}
		mean += fx;
	}
	mean/=N;
	for(int k=0; k<dim; k++){
		mean_l[k] /= n_l[k];
		mean_r[k] /= n_r[k];
	}
	int kdiv = 0;
	double maxvar = 0;
	for(int k=0; k<dim; k++){
		double var = fabs(mean_r[k]-mean_l[k]);
		if(var > maxvar){
			maxvar = var;
			kdiv=k;
		}
	}
	double integ = (mean*N+mean_reuse*n_reuse)/(N+n_reuse)*V;
	double error = fabs(mean_reuse-mean)*V;
	double toler = acc+fabs(integ)*eps;
	if(error<toler){
		return integ;
	}
	double a2[dim];
	double b2[dim];
	for(int k=0; k<dim; k++){
		a2[k] = a[k];
		b2[k] = b[k];
	}
	a2[kdiv] = (a[kdiv]+b[kdiv])/2;
	b2[kdiv] = (a[kdiv]+b[kdiv])/2;
	double integ_l = stratified_sam(dim,f,a,b2,acc/sqrt(2.0),eps,n_l[kdiv],mean_l[kdiv]);
	double integ_r = stratified_sam(dim,f,a2,b,acc/sqrt(2.0),eps,n_r[kdiv],mean_r[kdiv]);
	return integ_l+integ_r;
}
