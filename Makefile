CXX=g++
CXXFLAGS=-std=c++11 -Wall -O2
HDRS=LidDrivenCavity.h
OBJS=LidDrivenCavitySolver.o LidDrivenCavity.o
LDLIBS=-llapack -lcblas -lblas -lgfortran

%.o : %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

LidDrivenCavitySolver: $(OBJS)
	$(CXX) -o $@ $^ $(LDLIBS)

all: LidDrivenCavitySolver

PHONY: clean # Specify that ’clean’ is not a real file target


default: LidDrivenCavitySolver

clean:
	-rm -f *.o LidDrivenCavitySolver # Clean up (and ignore any errors)

run:
	./LidDrivenCavitySolver

