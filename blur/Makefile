# Author: David Holmqvist <daae19@student.bth.se>

CXX=g++
CXXFLAGS=-std=c++17 -O3 -Wunused -Wall -Wunused

all: blur blur_par

blur: matrix ppm filters blur.cpp
	$(CXX) $(CXXFLAGS) blur.cpp matrix.o ppm.o filters.o -o blur

blur_par: matrix ppm filters_par blur_par.cpp
	$(CXX) $(CXXFLAGS) blur_par.cpp matrix.o ppm.o filters_par.o -o blur_par

filters: matrix filters.hpp filters.cpp
	$(CXX) $(CXXFLAGS) -c filters.cpp -o filters.o

filters_par: matrix filters_par.hpp filters_par.cpp
	$(CXX) $(CXXFLAGS) -c filters_par.cpp -o filters_par.o -lpthreads

matrix: matrix.hpp matrix.cpp
	$(CXX) $(CXXFLAGS) -c matrix.cpp -o matrix.o

ppm: ppm.hpp ppm.cpp
	$(CXX) $(CXXFLAGS) -c ppm.cpp -o ppm.o

clean:
	rm -rf blur *.ppm *.o *.dSYM 2> /dev/null
