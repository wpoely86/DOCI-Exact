# makefile

CPPSRC= helpers.cpp\
	Permutation.cpp\
	Molecule.cpp\
	Hamiltonian.cpp\
	SparseMatrix_CCS.cpp\

CUDASRC=

OBJ=$(CPPSRC:.cpp=.o) $(CUDASRC:.cu=.o)

EXE=main

CC=gcc
CXX=clang++ #-fsanitize=address -fno-omit-frame-pointer

CFLAGS=-Iinclude -g -Wall -O2 -march=native -std=c++11 # -Wno-sign-compare -Wunknown-pragmas #-fopenmp # -DNDEBUG
CPPFLAGS=$(CFLAGS)
LDFLAGS=-g -O2 -Wall -march=native #-fopenmp
NVFLAGS=-g -O2 --ptxas-options=-v -arch=sm_13

INCLUDE=
#LIBS=-L. -L/usr/lib64/atlas -lblas -llapack -larpack -lhdf5
LIBS=-lblas -llapack -larpack -lhdf5
#LIBS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread -lhdf5 -larpack


%.o:    %.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(@:.o=.c) -o $@

%.o:    %.cpp
	$(CXX) -c $(CPPFLAGS) $(INCLUDE) $(@:.o=.cpp) -o $@

%.o:    %.cu
	nvcc -c $(NVFLAGS) $(INCLUDE) $(@:.o=.cu) -o $@
#	nvcc -cuda $(NVFLAGS) $(INCLUDE) $(@:.o=.cu) -o $(@:.o=.cu.ii)
#	$(CXX) -c $(CPPFLAGS) $(INCLUDE) $(@:.o=.cu.ii) -o $@


all: main

main: $(OBJ) main.o
	$(CXX) $(LDFLAGS) -o main $(OBJ) main.o $(LIBS)

doc: $(CPPSRC) doc-config
	doxygen doc-config

.PHONY: clean
clean:
	rm -f $(OBJ) main.o
