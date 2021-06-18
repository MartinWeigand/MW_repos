CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: results_taskAB.txt errordata.txt error.plot.png stra_sam_data.txt stra_sam.plot.png

error.plot.png: errordata.txt makefile
	echo '\
	set terminal png;\
	set title "Error comparison of Montecarlo integration of the given integral in Task A" ;\
	set out "$@" ;\
	set xlabel "Number of points" ;\
	set ylabel "Actual error" ;\
	set format x "%g" ;\
	set format y "%.0e" ;\
	set xtics 25000 ;\
	set mxtics 5 ;\
	a=1; b=-0.5; f(x)=a*x**b ;\
	c=1; d=-1; h(x)=c*x**d ;\
	fit f(x) "$<" using 1:2 via a,b ;\
	fit h(x) "$<" using 1:3 via c,d ;\
	plot \
	"$<" using 1:2 with points pt 5 title "Pseduo-random" \
	, f(x) with lines title sprintf("%3.1f*n^{%3.1f}",a,b) \
	, "$<" using 1:3 with points pt 5 title "Quasi-random" \
	, h(x) with lines title sprintf("%3.1f*n^{%3.1f}",c,d) \
	'| gnuplot

stra_sam.plot.png: stra_sam_data.txt makefile
	echo '\
	set term png ;\
	set out "$@" ;\
	set title "Stratified samling of x*x+y*y<0.8*0.8 (1 if true, otherwise 0)" ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set size square ;\
	plot \
	"$<" with  points pointtype 7 pointsize 0.2 notitle \
	'| gnuplot

main: main.o func.o

results_taskAB.txt: main
	./$< > $@

ns=$(shell seq 1000 1250 200000)
errordata.txt: main
	>$@
	@for n in $(ns); do echo n=$$n; ./main $$n >>$@; done

main1: main1.o func.o

stra_sam_data.txt: main1
	./$< > $@

.PHONEY: clean
clean:
	$(RM) *.o *.txt *.log *.png
