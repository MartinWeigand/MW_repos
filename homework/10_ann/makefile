CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: plot.png results.txt

plot.png: data.txt data2.txt makefile
	echo '\
	set terminal png ;\
	set out "$@" ;\
	set tics out ;\
	set title "Curve fitting of f(x)=sin(x)*exp(-x) with neural network" ;\
	set xlabel "x" ;\
	set ylabel "f(x)" ;\
	plot \
	 "$<" index 0 with points title "Exact" \
	,"$<" index 1 with lines title "Fit" \
	,"data2.txt" using 1:2 with lines title "Derivative" \
	,"data2.txt" using 1:4 with points title "Exact derivative" \
	,"data2.txt" using 1:3 with lines title "Antiderivative" \
	,"data2.txt" using 1:5 with points title "Exact antiderivative" \
	'| gnuplot

results.txt: main
	./$< > $@

main: main.o ann.o qnewton.o

.PHONEY: clean
clean:
	$(RM) *.o results.txt *.png log
