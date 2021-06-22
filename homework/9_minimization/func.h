#include <gsl/gsl_vector.h>

int qnewton(double cost(gsl_vector* x), gsl_vector* x, double acc);
int downhill(int dim, double F(double*), double* start, double* step,double simplex_size_goal);
