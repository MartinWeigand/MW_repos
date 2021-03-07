CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: out.data.txt ex.png ms.pdf

ms.pdf: ms.tex ex.png
	pdflatex $<
	pdflatex $<

out.data.txt: main
	./$< > $@

main: main.o

ex.png: out.data.txt makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "pow(e,x) ";\
	plot \
	"$<" using 1:3 with line title "Exponential function from math.h"\
	,"$<" using 1:2 with line title "Our implemented exponential function"\
' | pyxplot

clean:
	$(RM) ms.pdf *.log *.aux out* main *.o *.png fig*.tex
