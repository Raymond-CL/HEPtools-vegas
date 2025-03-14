# compilers
FC = gfortran
CC = gcc
CXX = g++

# rules for compiling the programs
all: cuba gsl legacy nr2c nr2cpp nr3 nrf77 nrf90

cuba: 
	$(FC) -o $@.exe $@/cuba.f90 -lcuba

gsl:
	$(CXX) -Wno-write-strings -o $@.exe $@/gsl.cpp -lgsl

legacy:
	$(FC) -std=legacy -o $@.exe $@/vegas.f $@/xvegas.f

nr2c:
	$(CC) -o $@.exe $@/nr.h $@/nrutil.h $@/nrutil.c $@/ran2.c $@/rebin.c $@/vegas.c $@/xvegas.c -lm

nr2cpp:
	$(CXX) -o $@.exe $@/nr.h $@/nrtypes_nr.h $@/nrtypes.h $@/nrutil_nr.h $@/nrutil.h $@/ran2.cpp $@/rebin.cpp $@/vegas.cpp $@/xvegas.cpp

nr3:
	$(CXX) -o $@.exe $@/xvegas.cpp

nrf77:
	$(FC) -std=legacy -fdefault-real-8 -o $@.exe $@/ran2.for $@/rebin.for $@/vegas.for $@/xvegas.for

nrf90:
	$(FC) -fno-strict-overflow -J $@/ -o $@.exe $@/nrtype.f90 $@/nr.f90 $@/nrutil.f90 $@/ran_state.f90 $@/ran1.f90 $@/vegas.f90 $@/xvegas.f90
	rm -f $@/*.mod

# rules for cleaning the directory
clean:
	rm -f *.exe *.mod *.o

# rules for phony targets
.PHONY: all clean cuba gsl legacy nr2c nr2cpp nr3 nrf77 nrf90