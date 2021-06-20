void driver(
	int n,
	void f(int n, double x, double* y, double* dydx),
	double a,
	double* y,
	double b,
	double h,
	double acc,
	double eps,
	char file[]
);
