CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: results.txt

main: main.o integrator_func.o cc_func.o infint_func.o

results.txt: main
	./$< > $@

.PHONEY: clean
clean:
	$(RM) *.o *.txt
