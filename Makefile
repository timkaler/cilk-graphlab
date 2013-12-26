CC=icpc
CFLAGS= -O3
CILK=icpc
CILKFLAGS= -Wall -Werror -O3
LDFLAGS= -L$(CURDIR) -lopencv_highgui -lopencv_imgproc -lopencv_core
AR=ar

OPENCV=/usr/local/include/opencv
all: main

main : main.cpp Graph.cpp Graph.h bag.cpp bag.h scheduler.cpp scheduler.h engine.cpp engine.h Makefile
	$(CILK) $(CILKFLAGS) -I$(OPENCV) $@.cpp $(LDFLAGS) -o $@

clean :
	rm -f main *~

run :
	touch V$(V)_D$(D).graph
	rm V$(V)_D$(D).graph
	python graph_gen.py $(V) $(D) >> V$(V)_D$(D).graph
	./main V$(V)_D$(D).graph
	rm V$(V)_D$(D).graph
