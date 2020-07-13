#include Makefile from upper-level directory
#we get everything defined so far

# this makefile can be made better!
# Nachiket Gokhale gokhalen@bu.edu 

#include ../Makefile

TESTEXECS=testshape1d
TESTOBJS=testshape1d.o testshape2d.o
CLINKER=g++

#objects are not made, executables are made
alltst:${TESTOBJS}



testshape1d.o:testshape1d.cpp
	${CLINKER} -o -g testshape1d  ../../datastruct/shapestruct.cpp ../shape1d.cpp ../../integlib/gauss1d.cpp ./testshape1d.cpp ../checkjacobiandet.cpp

testshape2d.o:testshape2d.cpp
	${CLINKER} -o -g testshape2d  ../../datastruct/shapestruct.cpp ../shape2d.cpp ../../integlib/gauss2d.cpp ./testshape2d.cpp ../checkjacobiandet.cpp

testshape3d.o:testshape3d.cpp
#	${CLINKER} -o testshape3d testshape3d.cpp
	${CLINKER} -o testshape3d  -g ../../datastruct/shapestruct.cpp ../shape3d.cpp ../../integlib/gauss3d.cpp ./testshape3d.cpp ../checkjacobiandet.cpp

#running tests
runshape1d:
	./runshape1d
runshape2d:
	./runshape2d
runshape3d:
	./runshape3d
