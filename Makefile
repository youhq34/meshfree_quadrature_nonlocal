CC=gcc
CPP=g++
CFLAGS=-g -Wall

PDDIR=./bond-based-pd
NLDIFFDIR=./nonlocal_diffusion
KWDIR=./Kalthoff-Winkler-test

LAPACKFLAG=-llapack 

pd:
	$(CPP) $(PDDIR)/PMB_2Dweight.cpp $(LAPACKFLAG) -o PMB_2d.ex 
Nldiff:
	$(CPP) $(NLDIFFDIR)/nonlocaldiff.cpp $(LAPACKFLAG) -o nldiff.ex
KW:
	$(CPP) $(KWDIR)/KW_2Dweight_dynamic.cpp $(LAPACKFLAG) -o KW.ex
clean:
	rm -f *.ex

