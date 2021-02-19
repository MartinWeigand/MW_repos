CFLAGS = -O -Wall -std=gnu11
#CFLAGS += `gsl-config --cflags`
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY:default
default: gamma.pyxplot.png

gamma.pyxplot.png: data.txt data1.txt makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key bottom right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set yrange [-5:5];\
	set xrange [-5:5];\
	set title "Gamma function";\
	plot \
	 "$<" using 1:2 with line title "tgamma from math.h"\
	,"data1.txt" using 1:2 with line title "gsl-sf-gamma"\
	,"$<" using 1:3 with line title "Approximation"\
	' | pyxplot

data.txt: main
	./$< > $@

main: main.o mygamma.o
main.o mygamma.o: mygamma.h

data1.txt: main1
	./$< > $@

main1: main1.o mygamma.o
main1.o mygamma.o: mygamma.h

.PHONEY:clean
clean:
	$(RM) *.o main main1 *.png data.txt data1.txt

.PHONEY:test
test:
	echo $(CFLAGS)
	echo $(LDLIBS)
