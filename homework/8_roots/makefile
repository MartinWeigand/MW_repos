CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: results.txt data.txt

main: main.o qr_func.o newtons_method.o ode.o

results.txt: main
	./$< > $@

data.txt: main

.PHONEY: clean
clean:
	$(RM) *.o *.txt
