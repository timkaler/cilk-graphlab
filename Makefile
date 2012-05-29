CC=icpc
CFLAGS= -O3
CILK=icpc
CILKFLAGS= -Wall -Werror -O3
LDFLAGS= -L$(CURDIR)
AR=ar

all: main

main : main.cpp Graph.cpp Graph.h bag.cpp bag.h scheduler.cpp scheduler.h Makefile
	$(CILK) $(CILKFLAGS) $@.cpp $(LDFLAGS) -o $@

clean :
	rm -f main *~
