#include<gsl/gsl_vector.h>
#include<stdio.h>
#include<math.h>
#include"ann.h"

int qnewton(double f(gsl_vector* x), gsl_vector* x, double eps);

ann* ann_alloc(int n, double(*f)(double)){
	ann* network = malloc(sizeof(ann));
	network->n = n;
	network->f = f;
	network->params = gsl_vector_alloc(3*n);
	return network;
}

double ann_derivative(ann* network, double x){
	double s = 0;
	for(int i=0; i<network->n; i++){
		double a = gsl_vector_get(network->params,3*i);
		double b = gsl_vector_get(network->params,3*i+1);
		double w = gsl_vector_get(network->params,3*i+2);
		s += w*(exp(-pow(x-a,2)/pow(b,2))/b-2.0*pow(x-a,2)*exp(-pow(x-a,2)/pow(b,2))/pow(b,3));
	}
	return s;
}

double ann_antideri(ann* network, double x){
	double s = 0;
	for(int i=0; i<network->n; i++){
		double a = gsl_vector_get(network->params,3*i);
		double b = gsl_vector_get(network->params,3*i+1);
		double w = gsl_vector_get(network->params,3*i+2);
		s += -0.5*b*exp(-pow(a-x,2)/pow(b,2))*w;
	}
	return s;
}

void ann_free(ann* network){
	gsl_vector_free(network->params);
	free(network);
}

double ann_response(ann* network, double x){
	double s = 0;
	for(int i=0; i<network->n; i++){
		double a = gsl_vector_get(network->params,3*i);
		double b = gsl_vector_get(network->params,3*i+1);
		double w = gsl_vector_get(network->params,3*i+2);
		s += network->f((x-a)/b)*w;
	}
	return s;
}

void ann_train(ann* network, gsl_vector* xs, gsl_vector* ys){
	double cost_function(gsl_vector* p){
		gsl_vector_memcpy(network->params,p);
		double sum = 0;
		for(int i=0; i<xs->size; i++){
			double xi = gsl_vector_get(xs,i);
			double yi = gsl_vector_get(ys,i);
			double fi = ann_response(network,xi);
			sum += fabs(fi-yi);
		}
		return sum/xs->size;
	}

	gsl_vector* p = gsl_vector_alloc(network->params->size);
	gsl_vector_memcpy(p,network->params);
	qnewton(cost_function,p,1e-3);
	gsl_vector_memcpy(network->params,p);
	gsl_vector_free(p);
}
