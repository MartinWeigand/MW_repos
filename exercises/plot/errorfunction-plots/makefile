CFLAGS = -O -Wall -std=gnu11
#CFLAGS += `gsl-config --cflags`
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY:default
default: erf.pyxplot.png

erf.pyxplot.png: data.txt makefile tabulated.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "Approximation"\
	,"tabulated.data.txt" using 1:2 with points pointtype 2 title "Tabulated data"\
	' | pyxplot

data.txt: main
	./$< > $@

main: main.o myerf.o
main.o myerf.o: myerf.h

.PHONEY:clean
clean:
	$(RM) *.o main *.png data.txt

.PHONEY:test
test:
	echo $(CFLAGS)
	echo $(LDLIBS)
