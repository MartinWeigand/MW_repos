#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "function.c"

int main(){
int j = INT_MAX;

int i=1;
while (i+1>i){
	i++;
}
int k=1;
for (k;k+1>k;){
	k++;
}

int l=1;
do{l++;} while(l+1>l);
printf("Exercise 1.i:\n");
printf("While: My max int = %i\n",i);
printf("For: My max int = %i\n",k);
printf("Do-while: My max int = %i\n",l);
printf("INT_MAX = %i\n\n",j);

int jj=INT_MIN;
int ii=1;
while (ii-1<ii){
	ii--;
}

int kk=1;
for (kk;kk-1<kk;){
	kk--;
}

int ll=1;
do{ll--;} while(ll-1<ll);

printf("Exercise 1.ii:\n");
printf("While: My min int = %i\n",ii);
printf("For: My min int = %i\n",kk);
printf("Do-while: My min int = %i\n",ll);
printf("INT_MIN = %i\n\n",jj);

float y=1;
while(1+y!=1){y/=2;} y*=2;

double x=1;
while(1+x!=1){x/=2;} x*=2;

long double z=1;
while(1+z!=1){z/=2;} z*=2;

printf("Exercise 1.iii:\n");
printf("Float, while: Epsilon = %g\n", y);
printf("Double, while: Epsilon = %g\n",x);
printf("Long double, while: Epsilon = %Lg\n",z);

float yy; for(yy=1; 1+yy!=1; yy/=2){} yy*=2;
double xx; for(xx=1; 1+xx!=1; xx/=2){} xx*=2;
long double zz; for(zz=1; 1+zz!=1; zz/=2){} zz*=2;

printf("Float, for: Epsilon = %g\n",yy);
printf("Double, for: Epsilon = %g\n",xx);
printf("Long double, for: Epsilon = %Lg\n",zz);

float yyy=1; do{yyy/=2;} while(1+yyy!=1);
double xxx=1; do{xxx/=2;} while(1+xxx!=1);
long double zzz=1; do{zzz/=2;} while(1+zzz!=1);

printf("Float, do-while: Epsilon = %g\n",yyy);
printf("Double, do-while: Epsilon = %g\n",xxx);
printf("Long double, do-while: Epsilon = %Lg\n",zzz);

long double p = FLT_EPSILON;
long double pp = DBL_EPSILON;
long double ppp = LDBL_EPSILON;
printf("FLT_EPSILON = %Lg\n",p);
printf("DBL_EPSLION = %Lg\n",pp);
printf("LDBL_EPSILON = %Lg\n\n",ppp);

int max=INT_MAX/2;
float sum=0;
for(int q=1;q<=max;q++){sum=sum+pow(q,-1);}
printf("Exercise 2.i:\n");
printf("Sum_up_float=%g\n",sum);

float sum1=0;
for(int q=max;q>=1;q--){sum1=sum1+pow(q,-1);}
printf("Sum_down_float=%g\n",sum1);

double sumd=0;
for(int q=1;q<=max;q++){sumd=sumd+pow(q,-1);}
printf("Sum_up_double=%g\n",sumd);

double sumd1=0;
for(int q=max;q>=1;q--){sumd1=sumd1+pow(q,-1);}
printf("Sum_down_double=%g\n\n",sumd1);

printf("Exercise 3:\n");
int n = equal(2,1,0.1,0.0001);
int nn = equal(1,1,0.1,0.0001);
printf("For a=2, b=1, tau=0.1 and epsilon=0.0001, one get %i\n",n);
printf("For a=1, b=1, tau=0.1 and epsilon=0.0001, one get %i\n",nn);
return 0;
}
