CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: results.txt time.data.txt eigenfunctions.png timeplot.png

timeplot.png: time.data.txt makefile
	echo '\
	set terminal png;\
	set title "Matrix diagonalization times as function of matrix size";\
	set xlabel "Matrix size N";\
	set ylabel "Time, sec";\
	set output "$@";\
	set key left top;\
	b1=b2=1 ;\
	g1(x)=b1*x**3;\
	g2(x)=b2*x**3;\
	fit g1(x) "$<" using 1:2 via b1;\
	fit g2(x) "$<" using 1:4 via b2;\
	plot \
	"$<" using 1:2 with linespoints title "My implementation"\
	,"$<" using 1:3 with linespoints title "GSLs implementation"\
	,"$<" using 1:4 with linespoints title "My optimized implementation"\
	, g1(x) with lines title "x³ fit to my implementation"\
	, g2(x) with lines title "x³ fit to my optimized implementation";\
	'| gnuplot

eigenfunctions.png: states.data.txt makefile
	echo '\
	set terminal png;\
	set title "Eigenfunctions of infinite well";\
	set output "$@";\
	set key right bottom;\
	set xlabel "xi";\
	set ylabel "psi";\
	a1=a2=a3=1;\
	f1(x)=a1*sin(3.1459*x);\
	f2(x)=a2*sin(2*3.1459*x);\
	f3(x)=a3*sin(3*3.1459*x);\
	fit f1(x) "$<" using 1:2 via a1;\
	fit f2(x) "$<" using 1:3 via a2;\
	fit f3(x) "$<" using 1:4 via a3;\
	plot \
	"$<" using 1:2 with linespoints title "Numerical - State 1" \
	,"$<" using 1:3 with linespoints title "Numerical - State 2" \
	,"$<" using 1:4 with linespoints title "Numerical - State 3" \
	, f1(x) with lines title "Analytical - State 1" \
	, f2(x) with lines title "Analytical - State 2" \
	, f3(x) with lines title "Analytical - State 3" \
	'| gnuplot

results.txt: main
	./$< > $@

main: main.o functions.o

time.data.txt: timecalculations
	./$< > $@

timecalculations: timecalculations.o functions.o jacobi_optimized.o

.PHONEY: clean
clean:
	$(RM) *.o *.txt

