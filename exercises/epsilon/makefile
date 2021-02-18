CC = gcc
CFLAGS = -Wall -fwrapv
LDLIBS = -lm

-PHONEY: default
default: out.txt
	cat $<

out.txt: main
	./$< > $@
main: main.o function.o
	$(CC) -o main main.o $(LDLIBS)
main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@
function.o: function.c
	$(CC) $(CFLAGS) -c $< -o $@
.PHONEY: clean
clean: 
	$(RM) main out.txt
