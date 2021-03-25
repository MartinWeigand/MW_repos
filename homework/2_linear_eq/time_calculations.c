#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_linalg.h>
#include"functions.h"
#include<time.h>

int main(){
clock_t my_time=clock();
clock_t gsl_time=clock();
double my_time_c;
double gsl_time_c;
for(int n=300; n<=600; n+=10){
	gsl_matrix* A=gsl_matrix_alloc(n,n);
	gsl_matrix* R=gsl_matrix_alloc(n,n);
	gsl_vector* R_gsl=gsl_vector_alloc(n);
	for(int i=0; i<n; i++){
        	for(int j=0; j<n; j++){
                	gsl_matrix_set(A,i,j,((double) rand()/RAND_MAX));
        	}
	}
	my_time=clock();
	QR_dec(A,R);
	my_time=clock()-my_time;
	my_time_c=((double) my_time)/CLOCKS_PER_SEC;

	gsl_time=clock();
	gsl_linalg_QR_decomp(A,R_gsl);
	gsl_time=clock()-gsl_time;
	gsl_time_c=((double) gsl_time)/CLOCKS_PER_SEC;
	printf("%10i %10g %10g\n",n,my_time_c,gsl_time_c);
	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_vector_free(R_gsl);
}
return 0;
}
