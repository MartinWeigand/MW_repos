CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: results.txt out.times.txt out.times.png

out.times.png: out.times.txt makefile
	echo '\
	set term png ;\
	set title "QR-decomposition times as function of matrix size" ;\
	set xlabel "matrix size N" ;\
	set ylabel "time, sec" ;\
	set out "$@" ;\
	a=1 ;\
	f(x)=a*x**3 ;\
	fit f(x) "$<" via a ;\
	plot \
	 "$<" with linespoints title "measurement" \
	, f(x) with lines title "const*NÂ³";\
	'| gnuplot

out.times.txt: main makefile
	>$@
	for n in $$(seq 300 25 500);\
	do \time --format "$$n %U" --append --output $@ ./$< $$n ;\
	done

results.txt: main
	./$< > $@

.PHONEY: clean
clean:
	$(RM) *.o *.txt
