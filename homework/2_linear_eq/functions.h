void printm(gsl_matrix *A);
void vector_print(char s[],gsl_vector* v);
void QR_dec(gsl_matrix* A, gsl_matrix* R);
void QR_backsub(gsl_matrix *R, gsl_vector *x);
void QR_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x);
void QR_inv(gsl_matrix *Q, gsl_matrix *R, gsl_matrix *B);

