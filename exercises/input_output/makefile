CFLAGS = -O -Wall -std=gnu11
LDLIBS = -lm

.PHONEY:default
default:out.cmdline.txt

out.cmdline.txt:cmdline input.txt
	./cmdline $$(cat input.txt) > $@

input.txt: makefile
	cat /dev/null > $@
	for x in $$(seq 0 0.5236 6.2831); do echo "$$x" >> $@; done

.PHONEY:clean
clean:
	$(RM) input.txt cmdline out.cmdline.txt
