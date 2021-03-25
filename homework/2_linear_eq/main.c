#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_linalg.h>
#include"functions.h"
#define FMT "%7.3f"

int main(){
size_t n=5, m=4;
gsl_matrix* A=gsl_matrix_calloc(n,m);
gsl_matrix* R=gsl_matrix_alloc(m,m);
for(int i=0; i<n; i++){
	for(int j=0; j<m;j++){
		gsl_matrix_set(A,i,j,((double) rand()/RAND_MAX));
	}
}
printf("Task A.1:\n");
printf("Our random-generated matrix, A, is:\n");
printm(A);
QR_dec(A,R);
printf("The orthogonal matrix from the QR decomposition, Q, is:\n");
printm(A);
printf("The upper triangular matrix from the QR decomposition, R, is:\n");
printm(R);
printf("So R is upper triangular.\n");
gsl_matrix *Q_trans_Q=gsl_matrix_calloc(m,m);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0.0,Q_trans_Q);
printf("Q^T*Q, which should be the identity matrix, is:\n");
printm(Q_trans_Q);
printf("Which exactly is the identity matrix.\n");
gsl_matrix *QR =gsl_matrix_calloc(n,m);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,R,0.0,QR);
printf("Q*R, which should be equal to A, is:\n");
printm(QR);
printf("This equals A.\n");
printf("\n");
printf("Task A.2:\n");
gsl_matrix *A_new = gsl_matrix_calloc(m,m);
gsl_matrix *Q_new = gsl_matrix_alloc(m,m);
for(int i=0; i<m; i++){
	for(int j=0; j<m;j++){
		gsl_matrix_set(A_new,i,j,((double) rand()/RAND_MAX));
	}
}
gsl_matrix_memcpy(Q_new,A_new);
printf("Our random-generated matrix, A, is:\n");
printm(A_new);
gsl_vector *b =gsl_vector_alloc(m);
for(int i=0; i<m; i++){
	gsl_vector_set(b,i,((double) rand()/RAND_MAX));
}
printf("Our random-generated vector, b, is:\n");
gsl_vector_fprintf(stdout,b,FMT);
printf("We now wanna solve Ax=b for x.\n");
gsl_matrix *R_new = gsl_matrix_calloc(m,m);
QR_dec(Q_new,R_new);
gsl_vector *x = gsl_vector_alloc(m);
QR_solve(Q_new,R_new,b,x);
printf("Using our implementation we get x to be:\n");
gsl_vector_fprintf(stdout,x,FMT);
gsl_blas_dgemv(CblasNoTrans,1.0,A_new,x,0.0,b);
printf("Using GSL's implementation we get that Ax is equal to:\n");
gsl_vector_fprintf(stdout,b,FMT);
printf("which equals b, so our implementation works.\n");

printf("\n");
printf("Task B:\n");
printf("We reuse the matrix and decomposition from task A.2\n");
gsl_matrix *B = gsl_matrix_calloc(m,m);
QR_inv(Q_new,R_new,B);
printf("Using our implementation, we get that the inverse of A, i.e. B, is:\n");
printm(B);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A_new,B,0.0,R_new);
printf("Using GSL's implementations we get that AB is:\n");
printm(R_new);
printf("which is the identity matrix as expected.\n");
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,B,A_new,0.0,R_new);
printf("Calculating BA gives:\n");
printm(R_new);
printf("which also is the identity matrix. So it seems like our inverse-function works.\n");
printf("\n");
printf("Task C:\n");
printf("In the out.times.png the processing time to QR-factorize a NxN matrix is plotted as a function of N. This is done for both my own implementation and GSLs implementation. We would expect my own implementation to go like O(N^3), since in my QR_decom function is calculates a dot product and a linear combination of two columns. Each of these operations contributes with O(n). So I have fitted a y=a*x^3 curve to the data of my own implementation and it looks like it actually evolves like N^3, expect for some weird 'bumps' at specfic N-values.");
gsl_matrix_free(A);
gsl_matrix_free(Q_trans_Q);
gsl_matrix_free(R);
gsl_vector_free(b);
gsl_matrix_free(A_new);
gsl_matrix_free(Q_new);
gsl_matrix_free(R_new);
gsl_matrix_free(QR);
gsl_vector_free(x);
gsl_matrix_free(B);

return 0;
}
