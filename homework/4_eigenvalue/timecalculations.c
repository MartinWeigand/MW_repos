#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<time.h>
#define RND ((double)rand()/RAND_MAX-0.5)*2

int jacobi_diag(gsl_matrix* A, gsl_matrix* V);
int jacobi2(gsl_matrix* A, gsl_vector* e, gsl_matrix* V);

int main(){
clock_t my_time=clock();
clock_t gsl_time=clock();
clock_t my_time_opt=clock();
double my_time_c;
double gsl_time_c;
double my_time_opt_c;

for(int n=5; n<=100; n+=5){
	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* V = gsl_matrix_alloc(n,n);
	gsl_matrix* AA = gsl_matrix_alloc(n,n);
	gsl_matrix* AAA = gsl_matrix_alloc(n,n);
	gsl_vector* vec = gsl_vector_alloc(n);
	gsl_vector* e = gsl_vector_alloc(n);
	gsl_matrix_set_identity(V);

	for(int i=0; i<n; i++){
		for(int j=i; j<n; j++){
			double Aij=RND;
			gsl_matrix_set(A,i,j,Aij);
			gsl_matrix_set(A,j,i,Aij);
		}
	}

	gsl_matrix_memcpy(AA,A);
	gsl_matrix_memcpy(AAA,A);

	my_time=clock();
	jacobi_diag(A,V);
	my_time=clock()-my_time;
	my_time_c=((double) my_time)/CLOCKS_PER_SEC;

	gsl_matrix_set_identity(V);

	gsl_time=clock();
	gsl_linalg_SV_decomp_jacobi(AA,V,vec);
	gsl_time=clock()-gsl_time;
	gsl_time_c=((double) gsl_time)/CLOCKS_PER_SEC;

	gsl_matrix_set_identity(V);

	my_time_opt=clock();
	jacobi2(AAA,e,V);
	my_time_opt=clock()-my_time_opt;
	my_time_opt_c=((double) my_time_opt)/CLOCKS_PER_SEC;

	printf("%10i %10g %10g %10g\n",n,my_time_c,gsl_time_c,my_time_opt_c);

	gsl_matrix_free(A);
	gsl_matrix_free(AA);
	gsl_matrix_free(AAA);
	gsl_matrix_free(V);
	gsl_vector_free(vec);
}
return 0;
}
