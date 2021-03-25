CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: results.txt times.txt out.times.png

out.times.png: times.txt makefile
	echo '\
	set terminal png ;\
	set title "QR-decomposition times as function of matrix size" ;\
	set xlabel "matrix size N" ;\
	set ylabel "time, sec" ;\
	set output "$@" ;\
	a=1 ;\
	f(x)=a*x**3 ;\
	fit f(x) "$<" using 1:2 via a ;\
	plot \
	"$<" using 1:2 with linespoints title "My QR-decom function" \
	,"$<" using 1:3 with linespoints title "GSLs QR-decom function"\
	, f(x) with lines title "const*NÂ³";\
	'| gnuplot


results.txt: main
	./$< > $@

main: main.o functions.o

times.txt: time_calculations
	./$< >$@

time_calculations: time_calculations.o functions.o

.PHONEY: clean
clean:
	$(RM) *.o *.txt *.log *.png
