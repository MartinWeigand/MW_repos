------------------------------The project------------------------------
My studentnumber is 201705072. Since 72 % 22 = 6, I'm gonna make exam
project number 6:
Adaptive integration of complex-valued functions: Implement an adaptive
integrator which calculates the integral of a complex-valued function 
f(z) of a complex variable z along a straight line between two points 
in the complex plane.

-------------------------------Overview-------------------------------
integrator_func.c : Hold the functions that implement the adaptive
		     integrator of a complex line integral
cc_func.c         : Hold the functions that improve the adaptive 
		     integrator with the Chenshaw-Curtis variable
		     transformation
infint_func.c     : Hold the functions that improve the adaptive 
		     integrator to also accept an infinite limit
main.c            : The main-file that tests the different 
		     implementations on different functions and write 
		     the results into results.txt
results.txt       : The output of the tests

---------------Discussion of implementation and results---------------
In general an adaptive integrator is an algorithm that recursively
subdivide the integration interval into adaptively refined 
subintervals until an accuracy goal is met. This is done by 
approximating the integral with a higher and lower order rule (or
quadrature) and letting the difference between the two rules determine
if the accuracy goal is met.

For this project, I chose to follow the tasks given for homework 6 and 
try to generalize the implementation that I made for that homework to 
also work for complex line integrals. As in homework 6, I base my 
adaptive integrator on the rectangle and trapezium rule as the lower
and higher order rule.

The main thing that has been done is adding the headerfile complex.h 
and changing the different functions and variables from 'double' to 
'double complex'. Looking at the description of the headerfile I 
found the three functions that I would gonna need: CMPLX(x,y), 
creal(z) and cimag(z), where the first function sets a complex 
variable to x+y*i and the two last functions gives the real and 
imaginary part of z, respectively.

For calculating the tolerence and error to determine when the 
accuracy goal was met, I chose to linearly add the tolerence of the
real part together with the tolerence of the imaginary part and the
same for the error. Another option was to calculate it as the 
average tolerance of the real and imaginary part and the same for
the error. However, I thought the former made more sense.

In results.txt the output of the runned tests are shown. Firstly, I
conclude that my adaptive integrator works for complex line integrals 
which is seen by a zero actual error for both x+yi and x*x+i*x*y, 
which I knew the analytical result of. Secondly, I see that the 
amount of function calls (and therefore also the amount of 
recursions) is reduced significantly when complex line integrals 
with integrable divergiencies inside the integration interval is 
calculated with the Clenshaw-Curtis variable transformation. If 
there is no such divergence the transformation doesn't reduce the 
amount of function calls significantly as seen for x+yi. Thirdly, 
I conclude that I have succesfully implemented an adaptive 
integrator that can compute a complex line integral with an infinite 
limit. The exact results for these integrals were found with
Wolframalpha.com.

A last note: The results of the calculations are surprisingly 
accurate compared to the exact results. When only dealing with 
integrals of real-valued function, as in homework 6, the actual 
errors were much higher. This I found a bit weird as I would expect 
the calculations of integrals of real-valued functions to be more 
accurate than those of complex-valued functions.
