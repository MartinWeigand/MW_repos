#include <complex.h>
complex plainmontecarlo(int dim, double f(int dim, double* x), double* a, double* b, int N);
complex quasimontecarlo(int d, double f(int, double*), double a[], double b[], int N);
double stratified_sam(int dim, double f(int dim, double* x), double* a, double* b, double acc, double eps, int nreuse, double mreuse);
