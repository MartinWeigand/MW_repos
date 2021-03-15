CC = gcc
CFLAGS = -Wall -Werror -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: out.linear.data.txt out.linear.intdata.txt out.quad.data.txt out.quad.intdata.txt out.cubic.data.txt out.cubic.intdata.txt lin_int.plot.png quad_int.plot.png cubic_int.plot.png

cubic_int.plot.png: out.cubic.intdata.txt out.cubic.data.txt makefile
	echo '\
	set term pngcairo noenhanced ;\
	set title "Cubic interpolation" ;\
	set key at graph 0.6, 0.98 ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set out "$@" ;\
	plot \
	"out.cubic.data.txt" with points pointtype 4 linecolor rgb "black" title "Data" \
	,"$<" using 1:2 with lines linecolor rgb "red" title "Cubic-spline" \
	,"$<" using 1:3 with lines linecolor rgb "blue" title "Integral(Cubic-spline)" \
	,"$<" using 1:4 with lines linecolor rgb "green" title "Derivative(Cubic-spline)" \
	,"$<" using 1:5 with lines linecolor rgb "red" title "GSL: Cubic-spline" \
	,"$<" using 1:6 with lines linecolor rgb "blue" title "GSL: Integral(Cubic-spline)" \
	,"$<" using 1:7 with lines linecolor rgb "green" title "GSL: Derivative(Cubic-spline)" \
' | gnuplot

quad_int.plot.png: out.quad.intdata.txt out.quad.data.txt makefile
	echo '\
	set term pngcairo noenhanced ;\
	set title "Quadratic interpolation" ;\
	set key at graph 0.6, 0.98 ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set out "$@" ;\
	plot \
	"out.quad.data.txt" with points pointtype 4 linecolor rgb "black" title "Data" \
	,"$<" using 1:2 with lines linecolor rgb "red" title "Quadratic-spline" \
	,"$<" using 1:3 with lines linecolor rgb "blue" title "Integral(Quadratic-spline)" \
	,"$<" using 1:4 with lines linecolor rgb "green" title "Derivative(Quadratic-spline)" \
	' | gnuplot

lin_int.plot.png: out.linear.intdata.txt out.linear.data.txt makefile
	echo '\
	set term pngcairo noenhanced ;\
	set title "Linear interpolation" ;\
	set key at graph 0.6, 0.98 ; \
	set xlabel "x" ;\
	set ylabel "y" ;\
	set out "$@" ;\
	plot \
	"out.linear.data" with points pointtype 4 linecolor rgb "black" title "Data" \
	,"$<" using 1:2 with lines linecolor rgb "red" title "Linear-spline" \
	,"$<" using 1:3 with lines linecolor rgb "blue" title "Integral(Linear-spline)" \
	,"$<" using 1:4 with lines linecolor rgb "red" title "GSL: Linear spline" \
	,"$<" using 1:5 with lines linecolor rgb "blue" title "GSL: Integral (Linear spline)" \
	' | gnuplot

out.linear.data.txt: linear
	./$< > $@

out.linear.dataint.txt: linear
	./$< > $@

out.quad.data.txt: quad
	./$< > $@

out.quad.intdata.txt: quad
	./$< >$@

out.cubic.data.txt: cubic
	./$< >$@

out.cubic.intdata.txt: cubic
	./$< >$@

.PHONEY: clean
clean:
	$(RM) *.txt *.o *.png
