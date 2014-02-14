CC=icc
CFLAGS= -O3 -g
CILK=icc
CILKFLAGS= -Wall -Werror -O3 -g
LDFLAGS= -L$(CURDIR)
AR=ar

all: pagerank

pagerank : apps/pagerank.cpp Graph.cpp Graph.h multibag.h scheduler.cpp scheduler.h engine.cpp engine.h Makefile
	$(CILK) $(CILKFLAGS) apps/$@.cpp $(LDFLAGS) -o $@

clean :
	rm -f main pagerank *~

lint :
	python cpplint.py *.cpp *.h

run_pagerank : pagerank
	touch V$(V)_D$(D).graph
	rm V$(V)_D$(D).graph
	python graph_gen.py $(V) $(D) >> V$(V)_D$(D).graph
	./pagerank V$(V)_D$(D).graph
	rm V$(V)_D$(D).graph
