#include<math.h>
#include<stdio.h>
#include<assert.h>
#include<stdlib.h>
#include<gsl/gsl_interp.h>

double linterp(int n, double* x, double* y, double z){ //interpolates at z
assert(n>1 && z>=x[0] && z<=x[n-1]); //z should be inside the tabulated data
int i=0, j=n-1;
while(j-i>1){ //binary search
	int mid=(i+j)/2; //splits the data in half
	if(z>x[mid]) i=mid; //if z is larger than the middle element, foucs only on the data above
		else j=mid; //opposite than above
}
double dy=y[i+1]-y[i];
double dx=x[i+1]-x[i];
assert(dx>0);
return y[i]+dy/dx*(z-x[i]); //the interpolated value at z is returned
}

double linterp_integ(int n, double* x, double* y, double z){
int i=0, j=n-1;
while(j-i>1){
	int mid=(i+j)/2;
	if(z>x[mid]) i=mid;
		else j=mid;
}
double integral=0;
for(int k=1;k<=i;k++){ //takes care of the integral from 0 and up to the datapoint just before z
	integral+=y[k-1]*(x[k]-x[k-1])+0.5*(y[k]-y[k-1])*(x[k]-x[k-1]);
	}

double y_at_z = linterp(n,x,y,z);
integral+=y[i]*(z-x[i])+0.5*(y_at_z-y[i])*(z-x[i]); //integral from the last datapoint before z, to z.
return integral;
}

double linterp_gsl(int n, double* x, double* y, double z){
gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_linear,n);
gsl_interp_init(interpolation,x,y,n);
gsl_interp_accel*accelerator = gsl_interp_accel_alloc();
//gsl_interp_free(interpolation);
return gsl_interp_eval(interpolation,x,y,z,accelerator);
}

double linterp_integ_gsl(int n, double* x, double* y, double z){
int i=0, j=n-1;
while(j-i>1){
	int mid=(i+j)/2;
	if(z>x[mid]) i=mid;
		else j=mid;
}
gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_linear,n);
gsl_interp_init(interpolation,x,y,n);
gsl_interp_accel*accelerator = gsl_interp_accel_alloc();
double integral=0;
integral+= gsl_interp_eval_integ(interpolation,x,y,0,z,accelerator);

//double y_at_z = linterp_gsl(n,x,y,z);
//integral+= gsl_interp_eval_integ(interpolation,x,y,i,y_at_z,accelerator);
gsl_interp_free(interpolation);
return  integral;
}

int main(){
FILE* my_os=fopen("out.data.txt","w");

int num_data = 6;
double x[num_data], y[num_data];

for(int i=0;i<num_data;i++){
y[i]=-2*pow(i,2)+10*i;
x[i]=i;
fprintf(my_os,"%10g %10g \n",x[i],y[i]);
}

int z=0;
double ex_points = 5; //numbers of interpolated data points between two tabulated data points
FILE* my_osn=fopen("out.intdata.txt","w");
for(z=0;z<=ex_points*(num_data-1);z++){
	fprintf(my_osn,"%10g %10g %10g %10g %10g \n",z/ex_points,linterp(num_data,x,y,z/ex_points),linterp_integ(num_data,x,y,z/ex_points),linterp_gsl(num_data,x,y,z/ex_points),linterp_integ_gsl(num_data,x,y,z/ex_points));
	}


return 0;
}

