#include<math.h>
#include<gsl/gsl_matrix.h>

void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c = cos(theta), s = sin(theta);
	for(int i=0; i<A->size1; i++){
		double new_Aip = c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		double new_Aiq = s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,new_Aip);
		gsl_matrix_set(A,i,q,new_Aiq);
	}
}

void Jtimes(gsl_matrix* A, int p, int q, double theta){
	double c = cos(theta), s = sin(theta);
	for(int j=0; j<A->size2; j++){
		double new_Apj = c*gsl_matrix_get(A,p,j)+s*gsl_matrix_get(A,q,j);
		double new_Aqj = -s*gsl_matrix_get(A,p,j)+c*gsl_matrix_get(A,q,j);
		gsl_matrix_set(A,p,j,new_Apj);
		gsl_matrix_set(A,q,j,new_Aqj);
	}
}

int jacobi_diag(gsl_matrix* A, gsl_matrix* V){
	gsl_matrix_set_identity(V);
	int n=A->size1, changed, sweeps=0;
	do{
		changed=0;
		sweeps++;
		for(int p=0; p<n-1; p++){
			for(int q=p+1; q<n; q++){
				double Apq=gsl_matrix_get(A,p,q);
				double App=gsl_matrix_get(A,p,p);
				double Aqq=gsl_matrix_get(A,q,q);
				double theta=0.5*atan2(2*Apq,Aqq-App);
				double c=cos(theta),s=sin(theta);
				double new_App=c*c*App-2*s*c*Apq+s*s*Aqq;
				double new_Aqq=s*s*App+2*s*c*Apq+c*c*Aqq;
				if(new_App!=App || new_Aqq!=Aqq){
					changed=1;
					timesJ(A,p,q,theta);
					Jtimes(A,p,q,-theta);
					timesJ(V,p,q,theta);
				}
			}
		}
	}while(changed!=0);
	return sweeps;
}
