#include Makefile from upper-level directory
#we get everything defined so far

# Nachiket Gokhale gokhalen@bu.edu 

include ../Makefile

TESTEXECS=testgauss1d testgauss2d
TESTOBJS=testgauss1d.o testgauss2d.o

#objects are not made, executables are made
alltst:${TESTOBJS}

testgauss1d.o:testgauss1d.cpp
	${CLINKER} -o testgauss1d  ../../datastruct/shapestruct.cpp ../gauss1d.cpp testgauss1d.cpp

testgauss2d.o:testgauss2d.cpp
	${CLINKER} -o testgauss2d  ../../datastruct/shapestruct.cpp ../gauss2d.cpp testgauss2d.cpp

#running tests
rungauss1d:
	./rungauss1d
