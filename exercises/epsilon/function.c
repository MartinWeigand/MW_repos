#include <stdio.h>
#include <math.h>

int equal(double a, double b, double tau, double epsilon){
	if(fabs(a-b)<tau){return 1;}
	if(fabs(a-b)/(fabs(a)+fabs(b))<(epsilon/2)){return 1;}
	else{return 0;}
}

//int main(){
//int n=equal(2,1,0.1,0.0001);
//printf("Equal is %g\n",n);
//return 0;
//}
