#include <complex.h>

double complex adapt24(double complex f(double complex), double complex a, double complex b, double acc, double eps, double complex f2, double complex f3, int recursions);
double complex adapt(double complex f(double complex), double complex a, double complex b, double acc, double eps);
