CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

.PHONEY: default
default: out.txt
	cat ./out.txt
out.txt: main
	./$< > $@

.PHONEY: clean
clean: 
	$(RM) main out.txt
