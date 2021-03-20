#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#define FMT "%7.3f"

void printm(gsl_matrix *A){
	for(int r=0;r<A->size1;r++){
		for(int c=0; c<A->size2;c++){
			printf(FMT,gsl_matrix_get(A,r,c));
		}
		printf("\n");
	}
}

void vector_print(char s[],gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++){
		printf("%10g ",gsl_vector_get(v,i));
	}
	printf("\n");
}

void QR_dec(gsl_matrix* A, gsl_matrix* R){ //This function replaces A with Q and determine R and put into the argument R
	int m = A->size2; //size2=columns --> m=#columns in A
	for(int i=0;i<m;i++){
		gsl_vector_view e = gsl_matrix_column(A,i); //Seperates the i'th column in A
		double r = gsl_blas_dnrm2(&e.vector); //Norm of the i'th column in A
		gsl_matrix_set(R,i,i,r);
		gsl_vector_scale(&e.vector,1/r); //Normalize e_i
		for(int j=i+1;j<m;j++){
			gsl_vector_view q =gsl_matrix_column(A,j); //Seperates columns in A that is to the right of the i'th column
			double s=0;
			gsl_blas_ddot(&e.vector,&q.vector,&s); //Dotproduct of e^T and q, returning it to s
			gsl_blas_daxpy(-s,&e.vector,&q.vector); //-s*e+q-->q
			gsl_matrix_set(R,i,j,s);
			gsl_matrix_set(R,j,i,0);
		}
	}
}

void QR_backsub(gsl_matrix *R, gsl_vector *x){ //leaves the result in x
	int m= R->size1;
	for(int i=m-1;i>=0;i--){
		double s=0;
		for(int k=i+1;k<m;k++){
			s+=gsl_matrix_get(R,i,k)*gsl_vector_get(x,k);
		}
		gsl_vector_set(x,i,(gsl_vector_get(x,i)-s)/gsl_matrix_get(R,i,i));
	}
}

void QR_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x){
	gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,x);  //Calculates Q^T*b and leaves it in x (because Rx=Q^T*b)
	QR_backsub(R,x); // Turns x into the solution
}

void QR_inv(gsl_matrix *Q, gsl_matrix *R, gsl_matrix *B){
	int n = R->size1;
	gsl_vector *b = gsl_vector_calloc(n);
	gsl_vector *x = gsl_vector_calloc(n);
	for(int i=0;i<n;i++){
		gsl_vector_set(b,i,1.0);
		QR_solve(Q,R,b,x);
		gsl_vector_set(b,i,0.0);
		gsl_matrix_set_col(B,i,x);
	}
	gsl_vector_free(b);
	gsl_vector_free(x);
}

int main(int argc, char** argv){
size_t n=5, m=4;
if(argc>1)n=atoi(argv[1]);
if(argc>2)m=atoi(argv[2]);
if(n<7){
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
}
else{
printf("n=%li\n",n);
gsl_matrix* A=gsl_matrix_alloc(n,n);
gsl_matrix* R=gsl_matrix_alloc(n,n);
for(int i=0; i<m; i++){
	for(int j=0; j<m; j++){
		gsl_matrix_set(A,i,j,((double) rand()/RAND_MAX));
	}
}
QR_dec(A,R);
gsl_matrix_free(A);
}
}
