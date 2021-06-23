#include <stdio.h>
#include <math.h>
#include "ann.h"
#include <gsl/gsl_vector.h>

double act_func(double x){
	return x*exp(-x*x);
}

double func_to_fit(double x){
	return sin(x)*exp(-x);
}

double func_to_fit_deri(double x){
	return (cos(x)-sin(x))*exp(-x);
}

double func_to_fit_antideri(double x){
	return -0.5*(cos(x)+sin(x))*exp(-x);
}

int main(){
	int n = 8; //# of neurons
	ann* network = ann_alloc(n,act_func);
	double a = -3;
	double b = 0;
	int nx = 100;
	gsl_vector* vx = gsl_vector_alloc(nx);
	gsl_vector* vy = gsl_vector_alloc(nx);
	for(int i=0; i<nx; i++){
		double x = a+(b-a)*i/(nx-1);
		double f = func_to_fit(x);
		gsl_vector_set(vx,i,x);
		gsl_vector_set(vy,i,f);
	}
	for(int i=0; i<network->n; i++){
		gsl_vector_set(network->params,3*i,a+(b-a)*i/(network->n-1));
		gsl_vector_set(network->params,3*i+1,1);
		gsl_vector_set(network->params,3*i+2,1);
	}
	ann_train(network,vx,vy);
	FILE* stream = fopen("data.txt","w");
	for(int i=0; i<vx->size; i++){
		double x = gsl_vector_get(vx,i);
		double f = gsl_vector_get(vy,i);
		fprintf(stream,"%g %g\n",x,f);
	}
	fprintf(stream,"\n\n");
	for(double z=a; z<=b; z+=1.0/10){
		double y = ann_response(network,z);
		fprintf(stream,"%g %g\n",z,y);
	}
	fclose(stream);
	printf("Task A and B:\n\n");
	printf("The artifical neural network is implemented and tested to interpolate a function which is illustrated in plot.png\n");
	printf("%i neurons was used and the i'th neurons parameters are:\n",n);
	for(int i=0; i<network->n; i++){
		double ai = gsl_vector_get(network->params,3*i);
		double bi = gsl_vector_get(network->params,3*i+1);
		double wi = gsl_vector_get(network->params,3*i+2);
		printf("i = %i , ai = %g , bi = %g , wi = %g\n",i,ai,bi,wi);
	}
	printf("The derivative and antiderivative is also approimated with the nerual network as shown in plot.png. The derivative match the exact quite well, while the behaviour of the antiderivative is also correct but is numerical off due to some integration constant.");
	FILE* stream2 = fopen("data2.txt","w");
	for(double i=a; i<=b; i+=1.0/10){
		double deri = ann_derivative(network,i);
		double antideri = ann_antideri(network,i);
		double exact_deri = func_to_fit_deri(i);
		double exact_antideri = func_to_fit_antideri(i);
		fprintf(stream2,"%g %g %g %g %g \n",i,deri,antideri,exact_deri,exact_antideri);
	}
	fclose(stream2);
	gsl_vector_free(vx);
	gsl_vector_free(vy);
	ann_free(network);
	return 0;
}
