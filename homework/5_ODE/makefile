CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: out.sin.data out.SIR1.data out.SIR2.data out.SIR3.data out.3body.data sin.plot.png SIR.plot.png SIR.comparison.plot.png 3body.plot.png

sin.plot.png: out.sin.data makefile
	echo '\
	set terminal png;\
	set title "u´´=-u, u(0)=0, u´(0)=1" ;\
	set out "$@" ;\
	set key out ;\
	set grid ;\
	set xlabel "x/pi" ;\
	set ylabel "y" ;\
	plot \
	"$<" using ($$1/pi):2 with linespoints pt 4 ps 0.5 title "u"\
	,"$<" using ($$1/pi):3 with linespoints pt 5 ps 0.5 title "u´";\
	'| gnuplot

out.sin.data: main
	./main 1 2>$@

SIR.plot.png: out.SIR1.data makefile
	echo '\
	set terminal png;\
	set title "SIR-model (T_c = 3 days)" ;\
	set out "$@" ;\
	set key out ;\
	set grid ;\
	set xlabel "Time [days]" ;\
	set ylabel "y" ;\
	plot \
	"$<" using 1:2 with linespoints pt 4 ps 0.5 title "Susceptible"\
	,"$<" using 1:3 with linespoints pt 5 ps 0.5 title "Infectious"\
	,"$<" using 1:4 with linespoints pt 6 ps 0.5 title "Removed";\
	'| gnuplot

out.SIR1.data: main
	./main 2 2>$@

SIR.comparison.plot.png: out.SIR1.data out.SIR2.data out.SIR3.data makefile
	echo '\
	set terminal png;\
	set title "SIR-model for different values of T_c" ;\
	set out "$@" ;\
	set key out ;\
	set grid ;\
	set xlabel "Time [days]" ;\
	set ylabel "y" ;\
	plot \
	"$<" using 1:3 with linespoints pt 4 ps 0.5 title "T_c = 3 days"\
	,"out.SIR2.data" using 1:3 with linespoints pt 5 ps 0.5 title "T_c = 1 day"\
	,"out.SIR3.data" using 1:3 with linespoints pt 6 ps 0.5 title "T_c = 0.5 day";\
	'| gnuplot

out.SIR2.data: main
	./main 3 2>$@

out.SIR3.data: main
	./main 4 2>$@

3body.plot.png: out.3body.data makefile
	echo '\
	set terminal png;\
	set title "Paths for the stable solution of the three-body problem" ;\
	set out "$@" ;\
	set key out ;\
	set grid ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	plot \
	"$<" using 2:3 with lines title "Body 1"\
	,"$<" using 4:5 with lines title "Body 2"\
	,"$<" using 6:7 with lines title "Body 3";\
	'| gnuplot


out.3body.data: main
	./main 5 2>$@

main: main.o ode.o

.PHONEY: clean
clean:
	$(RM) *.o *.data *.png
