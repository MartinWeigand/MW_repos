CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: results.txt

results.txt: main
	./$< > $@

main: main.o func.o

.PHONEY: clean
clean:
	$(RM) *.o results.txt
