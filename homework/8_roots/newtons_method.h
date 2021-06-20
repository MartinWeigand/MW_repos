#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

int newton( void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double eps);
