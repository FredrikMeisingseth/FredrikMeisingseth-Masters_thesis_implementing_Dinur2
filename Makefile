
CC=gcc
CFLAGS= -Wall -pedantic -g -O2
INC=    -I.
LIBS=   -lm
DEPS=   $(wildcard *.h)
OBJ=    $(patsubst %.c,%.o,$(wildcard *.c))


%.o: %.c $(DEPS)
	$(CC) $(INC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	$(CC) $(INC) -o $@ $^ $(CFLAGS) $(LIBS)


run: main
	./main

.PHONY: clean

clean:
	rm -f main *.o
