#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_gamma.h>
#include"mygamma.h"

int main(){
	double xmin=-5,xmax=5;
	for(double x=xmin;x<=xmax;x+=1.0/1001){
		printf("%10g %10g %10g\n",x,tgamma(x),mygamma(x));
		}
return 0;
}
