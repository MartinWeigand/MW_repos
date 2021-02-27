#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>

void* pi_MC(void* arg){
	double* count = (double*) arg;
	double x=0; double y=0;
	unsigned int seedx=1; unsigned int seedy=2;
	for(int i=0;i<1e8;i++){
		x =rand_r(&seedx);
		x/=RAND_MAX;
		y =rand_r(&seedy);
		y/=RAND_MAX;
		if((pow(x,2)+pow(y,2))<=1){
		*count = *count+1;
		}
	}
	return NULL;
}

int main(){
double c1=0; double c2=0; double c3=0;

pthread_t t1,t2;
pthread_attr_t* attributes = NULL;

//pi_MC((void*)&c1);
//pi_MC((void*)&c2);
pthread_create(&t1,attributes,pi_MC,(void*)&c1);
pthread_create(&t2,attributes,pi_MC,(void*)&c2);
pi_MC((void*)&c3);

void* returnvalue=NULL;
pthread_join(t1,returnvalue);
pthread_join(t2,returnvalue);

double pi = 4 * (c1+c2+c3)/(3e8);
printf("Estimate of pi with 3e8 iterations is %g\n",pi);
return 0;
}
