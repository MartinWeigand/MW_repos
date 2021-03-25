#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#define FMT "%7.3f"

void printm(gsl_matrix *A){
	for(int r=0;r<A->size1;r++){
		for(int c=0; c<A->size2;c++){
			printf(FMT,gsl_matrix_get(A,r,c));
		}
		printf("\n");
	}
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
