CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: results.txt data.txt s_hydrogen_function.png convergence_con1_data.txt convergence_con2_data.txt convergence_comparison.png

s_hydrogen_function.png: data.txt makefile
	echo '\
	set terminal png ;\
	set title "Radial s-wavefunction of hydrogen for r(max)=8" ;\
	set out "$@" ;\
	set xlabel "$$r$$ [$$a_\mathrm{Bohr}$$]" ;\
	set ylabel "$$F_{epsilon}(r)$$" ;\
	set tics out ;\
	plot \
	"$<" using 1:2 with points title "Calculated data" \
	, x*exp(-x) with lines title "Exact" \
	'| tee log.gpi | pyxplot

convergence_comparison.png: convergence_con1_data.txt convergence_con2_data.txt makefile
	echo '\
	set terminal png ;\
	set title "Comparison of convergence solution based on two boundary conditions" ;\
	set out "$@" ;\
	set xlabel "r(max)" ;\
	set ylabel "Epsilon" ;\
	set tics out ;\
	plot \
	"$<" using 1:2 with points title "Condition 1" \
	,"convergence_con2_data.txt" using 1:2 with points title "Condition 2" \
	, -0.5 with lines title "Exact" \
	'| tee tmp.gpi | pyxplot

main: main.o qr_func.o newtons_method.o ode.o

results.txt: main
	./$< > $@

data.txt: main

convergence_con1_data.txt: main

convergence_con2_data.txt: main

.PHONEY: clean
clean:
	$(RM) *.o *.txt *.png *.gpi
