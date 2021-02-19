#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_gamma.h>
#include"mygamma.h"

int main(){
	double xmin=0.01,xmax=5;
	for(double x=xmin;x<=xmax;x+=1.0/1001){
		printf("%10g %10g\n",x,exp(gsl_sf_lngamma(x)));
		}
return 0;
}
