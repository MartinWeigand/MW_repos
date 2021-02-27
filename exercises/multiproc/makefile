CFLAGS = -Wall -o
CFLAGS += -pthread
CFLAGS += -fopenmp
LDLIBS = -lm -lgomp

default: out.txt
	cat ./out.txt

out.txt: main
	./$< > $@

.PHONEY: clean
clean:
	$(RM) main out.txt
