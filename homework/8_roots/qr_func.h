#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

void QR_dec(gsl_matrix* A, gsl_matrix* R);
void QR_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x);
