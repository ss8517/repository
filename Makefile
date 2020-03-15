CXX=g++
CXXFLAGS=-std=c++11 -Wall -O2
HDRS=LidDrivenCavity.h
OBJS=Courseworktest.o LidDrivenCavity.o
LDLIBS=-llapack -lcblas -lblas -lgfortran

%.o : %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

Courseworktest: $(OBJS)
	$(CXX) -o $@ $^ $(LDLIBS)

all: Courseworktest

PHONY: clean # Specify that ’clean’ is not a real file target


default: Courseworktest

clean:
	-rm -f *.o Courseworktest # Clean up (and ignore any errors)

run:
	./Courseworktest

