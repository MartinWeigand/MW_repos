CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: results.txt plot.svg

plot.svg: results.txt makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "Time (days)";\
	set ylabel "ln(activity)";\
	set tics out;\
	set title "least sq. fit (linear) to the logarithm of the activity data of Rutherford and Soddy";\
	plot "$<" index 0 with errorbars title "data"\
	,"$<" index 1 with lines title "Least square fit"\
	,"$<" index 2 using 1:2 with lines title "fit with c0+deltac0"\
	,"$<" index 2 using 1:3 with lines title "fit with c0-deltac0"\
	,"$<" index 3 using 1:2 with lines title "fit with c1+deltac1"\
	,"$<" index 3 using 1:3 with lines title "fit with c1-deltac1"\
	'|gnuplot

results.txt: main
	./$< > $@

main: main.o functions.o

.PHONEY: clean
clean:
	$(RM) *.o *.txt *.svg

