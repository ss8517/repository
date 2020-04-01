CXX = mpicxx
CXXFLAGS = -std=c++11 -Wall -O2
HDRS = LidDrivenCavity.h PoissonSolver.h
OBJS = main.o LidDrivenCavity.o PoissonSolver.o
LDLIBS =/usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so -lcblas -lblas -lgfortran

%.o : %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

main: $(OBJS)
	$(CXX) -o $@ $^ $(LDLIBS)

all: main

PHONY: clean 


default: main 

clean:
	-rm -f *.o main # Clean up (and ignore any errors)

run:
	mpiexec -np $(np) ./main 
