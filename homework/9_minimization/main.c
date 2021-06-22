#include <math.h>
#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "func.h"

void vector_print(char s[], gsl_vector* vec){
	printf("%s\n",s);
	for(int i=0; i<vec->size; i++){
		printf("%10g \n",gsl_vector_get(vec,i));
	}
}

double Rosenbrock(gsl_vector* v){
	double x = gsl_vector_get(v,0);
	double y = gsl_vector_get(v,1);
	return pow(1-x,2)+100*pow(y-x*x,2);
}

double Himmelblau(gsl_vector* v){
	double x = gsl_vector_get(v,0);
	double y = gsl_vector_get(v,1);
	return pow(x*x+y-11,2)+pow(x+y*y-7,2);
}

double Himmelblau1(double* v){
	double x = v[0];
	double y = v[1];
	return pow(x*x+y-11,2)+pow(x+y*y-7,2);
}

double f(gsl_vector* v){
	double x = gsl_vector_get(v,0);
	return pow(x,7)-8*pow(x,6);
}

double Breit_Wigner(gsl_vector* v){
	double mass = gsl_vector_get(v,0);
	double Gamma = gsl_vector_get(v,1);
	double A = gsl_vector_get(v,2);

	int data_points = 30;
	gsl_vector* energy = gsl_vector_alloc(data_points);
	gsl_vector* cross_section = gsl_vector_alloc(data_points);
	gsl_vector* err = gsl_vector_alloc(data_points);

	FILE* stream = fopen("Higgs_data.txt","r");
	int row = 0;
	double x,y,z;
	int i;
	while((i = fscanf(stream,"%lg %lg %lg",&x,&y,&z)) !=EOF){
		gsl_vector_set(energy,row,x);
		gsl_vector_set(cross_section,row,y);
		gsl_vector_set(err,row,z);
		row++;
	}
	fclose(stream);

	double deviation_function = 0;
	for(int i=0; i<data_points; i++){
		deviation_function += pow(A/(pow((gsl_vector_get(energy,i)-mass),2)+pow(Gamma,2)/4.0)-gsl_vector_get(cross_section,i),2)/pow(gsl_vector_get(err,i),2);
	}
	return deviation_function;
}
int main(){
	int n1 = 1;
	int n2 = 2;
	int n3 = 3;
	int steps;
	gsl_vector* x1 = gsl_vector_alloc(n1);
	gsl_vector* x2 = gsl_vector_alloc(n2);
	gsl_vector* x3 = gsl_vector_alloc(n3);
	double acc = 1e-3;

	printf("Task A:\n");

	printf("     Minimization of (x*x*x+y*y-z-9)^2\n");
	gsl_vector_set(x1,0,5);
	vector_print("Starting point:",x1);
	steps = qnewton(f,x1,acc);
	vector_print("Minimum is at:",x1);
	printf("Iteration steps = %i\n",steps);
	printf("Exact result is at 6.857142857\n\n");

	printf("     Minimization of Rosenbrock's function\n");
	gsl_vector_set(x2,0,20);
	gsl_vector_set(x2,1,10);
	vector_print("Starting point:",x2);
	steps = qnewton(Rosenbrock,x2,acc);
	vector_print("Minimum is at:",x2);
	printf("Iteration steps = %i\n",steps);
	printf("Exact result is (1,1)\n\n");

	printf("     Minimization of Himmelblau's function;\n");
	gsl_vector_set(x2,0,20);
	gsl_vector_set(x2,1,10);
	vector_print("Starting point:",x2);
	steps = qnewton(Himmelblau,x2,acc);
	vector_print("Minimum is at:",x2);
	printf("Iteration steps = %i\n\n",steps);

	gsl_vector_set(x2,0,5);
	gsl_vector_set(x2,1,-5);
	vector_print("Starting point:",x2);
	steps = qnewton(Himmelblau,x2,acc);
	vector_print("Minimum is at:",x2);
	printf("Iteration steps = %i\n\n",steps);


	gsl_vector_set(x2,0,-5);
	gsl_vector_set(x2,1,-5);
	vector_print("Starting point:",x2);
	steps = qnewton(Himmelblau,x2,acc);
	vector_print("Minimum is at:",x2);
	printf("Iteration steps = %i\n\n",steps);

	printf("Exact result is (3,2), (-2.805118,3.131312), (-3.779310,-3.283186) and (3.584428,-1.848126)\n\n");

	printf("Task B:\n\n");
	gsl_vector_set(x3,0,125);
	gsl_vector_set(x3,1,0.2);
	gsl_vector_set(x3,2,2);
	qnewton(Breit_Wigner,x3,acc);
	vector_print("The Breit-Wigner function is fitted to the data by minimizing the deviation function which leads to the parameters:",x3);
	printf("So according to our implementation the Higgs-boson has a mass of %g +/- %g Gev/c^2.\n\n",gsl_vector_get(x3,0),gsl_vector_get(x3,1));

	printf("Task C:\n\n");
	double x[] = {-2.0,2.0};
	double h[] = {1, 1};
	printf("We try to use the downhill simplex method on Himmelblau's function.\n\n");
	printf("Starting point: (%g,%g)\n",x[0],x[1]);
	int iterations = downhill(n2,Himmelblau1,x,h,1e-5);
	printf("Minimum is at: (%g,%g)\n",x[0],x[1]);
	printf("Iterations = %i\n",iterations);
}
