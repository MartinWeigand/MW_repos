#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#define RND ((double)rand()/RAND_MAX-0.5)*2
#define FMT "%7.3f"

void print_matrix(gsl_matrix *A){
	for(int r=0; r<A->size1; r++){
		for(int c=0; c<A->size2; c++){
			printf(FMT,gsl_matrix_get(A,r,c));
		}
	printf("\n");
	}
}

int jacobi_diag(gsl_matrix* A, gsl_matrix* V);

int main(){
size_t n=4;

gsl_matrix* A = gsl_matrix_alloc(n,n);
gsl_matrix* V = gsl_matrix_alloc(n,n);
gsl_matrix* AA = gsl_matrix_alloc(n,n);
gsl_matrix* X = gsl_matrix_alloc(n,n);
gsl_matrix* Y = gsl_matrix_alloc(n,n);

gsl_matrix_set_identity(V);

for(int i=0; i<n; i++){
	for(int j=i; j<n; j++){
		double Aij=RND;
		gsl_matrix_set(A,i,j,Aij);
		gsl_matrix_set(A,j,i,Aij);
	}
}

gsl_matrix_memcpy(AA,A);

printf("Task A:\n");
printf("The randomly generated symmetric matrix A is:\n");
print_matrix(A);
jacobi_diag(A,V);
printf("After matrix diagonalization with the Jacobi algorithm, the matrix A is:\n");
print_matrix(A);
printf("It is diagonal as expected\n");

printf("To prove that our implementation works we check if V^T*A*V=D, V*D*V^T=A and V^T*V=1, where V is the orthogonal matrix of eigenvectors and D is the diagonal matrix with the eigenvalues. First, calculating V^T*A*V gives:\n");
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,AA,0,X);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,X,V,0,Y);
print_matrix(Y);
printf("Which equals the matrix A after matrix diagonalization as expected\n");

printf("Second, calculating V*D*V^T gives:\n");
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,A,0,X);
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,X,V,0,Y);
print_matrix(Y);
printf("Which equals the matrix A as expected\n");

printf("Third, calculating V^T*V gives:\n");
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,V,0,X);
print_matrix(X);
printf("Which equals the identity matrix as expected\n");

printf("\n");

printf("Task B:\n");
int N=20;
double s=1.0/(N+1);
gsl_matrix* H = gsl_matrix_alloc(N,N);
for(int i=0; i<N-1; i++){
	gsl_matrix_set(H,i,i,-2);
	gsl_matrix_set(H,i,i+1,1);
	gsl_matrix_set(H,i+1,i,1);
}
gsl_matrix_set(H,N-1,N-1,-2);
gsl_matrix_scale(H,-1/s/s);
gsl_matrix* VV = gsl_matrix_alloc(N,N);
jacobi_diag(H,VV);

printf("Energy values in units of hbar²/2mL²:\n");
printf("Numerical  Exact\n");

for(int k=0; k<N/3; k++){
	double real_e = M_PI*M_PI*(k+1)*(k+1);
	double calculated_e = gsl_matrix_get(H,k,k);
	printf("%i %g %g\n",k,calculated_e,real_e);
}
printf("The energies seems correct\n");

FILE* my_os=fopen("states.data.txt","w");
fprintf(my_os,"0 0 0 0 0 0\n");
for(int i=0; i<N; i++){
	fprintf(my_os,"%g %g %g %g\n",(i+1.0)/(N+1),gsl_matrix_get(VV,i,0),gsl_matrix_get(VV,i,1),gsl_matrix_get(VV,i,2));
}
fprintf(my_os,"1 0 1 0 1 0\n");

gsl_matrix_free(A);
gsl_matrix_free(V);
gsl_matrix_free(AA);
gsl_matrix_free(X);
gsl_matrix_free(Y);
gsl_matrix_free(H);
gsl_matrix_free(VV);
return 0;
}
