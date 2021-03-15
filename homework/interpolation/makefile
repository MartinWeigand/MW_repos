CC = gcc
CFLAGS = -Wall -Werror -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: out.data.txt out.intdata.txt out.quad.data.txt out.quad.intdata.txt lin_int.plot.png quad_int.plot.png

quad_int.plot.png: out.quad.intdata.txt out.quad.data.txt makefile
	echo '\
	set term pngcairo noenhanced ;\
	set title "Quadratic interpolation" ;\
	set key top left ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set out "$@" ;\
	plot \
	"out.quad.data.txt" with points pointtype 4 linecolor rgb "black" title "data" \
	,"$<" using 1:2 with lines linecolor rgb "red" title "Quadratic-spline" \
	,"$<" using 1:3 with lines linecolor rgb "blue" title "Integral(Quadratic-spline)" \
	,"$<" using 1:4 with lines linecolor rgb "green" title "Derivative(Quadratic-spline)" \
' | gnuplot

lin_int.plot.png: out.data.txt out.intdata.txt makefile
	echo '\
	set term pngcairo noenhanced ;\
	set title "Linear interpolation" ;\
	set key top left ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set out "$@" ;\
	plot \
	"$<" with points pointtype 4 linecolor rgb "black" title "data" \
	,"out.intdata.txt" using 1:2 with lines linecolor rgb "red" title "Linear-spline" \
	,"out.intdata.txt" using 1:3 with lines linecolor rgb "red" title "Integral(linear-spline)" \
	,"out.intdata.txt" using 1:4 with lines linecolor rgb "blue" title "GSL: Linear interpolation" \
	,"out.intdata.txt" using 1:5 with lines linecolor rgb "blue" title "GSL: Integral (linear spline)" \
	' | gnuplot

out.data.txt: linear
	./$< > $@

out.dataint.txt: linear
	./$< > $@

out.quad.data.txt: quad
	./$< > $@

out.quad.intdata.txt: quad
	./$< >$@

.PHONEY: clean
clean:
	$(RM) *.txt *.o *.png
