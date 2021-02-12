#include "stdio.h"
#include "math.h"
#include "complex.h"
//# define M_E 2.7182818284590452354
//# define M_PI 3.1415926535897932384

int main(void){
	double gamma_5 = tgamma(5);
	double bessel_05 = j1(0.5);
	double complex z = csqrt(-2);
	double complex z1 =cexp(I*M_PI);
	double complex z2 = cpow(I,cexp(1));
	double complex z3 = cexp(I);
	double complex z4 = cpow(I,I);
	printf("\u0393=%g\n",gamma_5);
	printf("J_1(0.5)=%g\n",bessel_05);
	printf("sqrt(-2)=%g + %gI\n",creal(z),cimag(z));
	printf("exp(i\u03c0)=%f +%fI\n",creal(z1),cimag(z1)); 
	printf("i^e = %f + %fI\n",creal(z2),cimag(z2));
	printf("e^i = %f + %fI\n",creal(z3),cimag(z3));
	printf("i^i = %f + %fI\n",creal(z4),cimag(z4));
	float x_float = 1.f/9;
	double x_double = 1./9;
	long double x_long_double = 1.L/9;
	printf("The output of printing 1/9 depends on the type. A float gives %.25g, a double will give %.25lg and a long double will give %.25Lg\n",x_float,x_double,x_long_double);
return 0;
}
