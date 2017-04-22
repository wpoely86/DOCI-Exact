# Makefile: adjuct the CXX and LIBS variables as needed

CPPSRC= helpers.cpp\
	Permutation.cpp\
	Molecule.cpp\
	DOCIHamiltonian.cpp\
	SparseMatrix_CRS.cpp\
	DM2.cpp\
	SymMolecule.cpp\
	SimulatedAnnealing.cpp\
	LocalMinimizer.cpp\

OBJ=$(CPPSRC:.cpp=.o)

# name of exe
EXE=doci

# The compilers
ifeq ($(origin CC), default)
    CC = clang
endif

ifeq ($(origin CXX), default)
    CXX = g++
endif

# compile and link flags
CFLAGS=-Iinclude -Iextern/include -g -Wall -O2 -march=native -std=c++11 -fopenmp -Wno-sign-compare # -DNDEBUG
CPPFLAGS=$(CFLAGS)
LDFLAGS=-g -O2 -Wall -march=native -fopenmp

# location of headers and libraries
INCLUDE=
LIBS=-lblas -llapack -larpack -lhdf5 -Lextern -lsimanneal
#LIBS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread -lhdf5 -larpack


# You shouldn't have to change anything past this point

%.o:    %.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(@:.o=.c) -o $@

%.o:    %.cpp
	$(CXX) -c $(CPPFLAGS) $(INCLUDE) $(@:.o=.cpp) -o $@

all: extern $(EXE)

.PHONY: extern
extern:
	$(MAKE) -C extern

$(EXE): $(OBJ) $(EXE).o
	$(CXX) $(LDFLAGS) -o $(EXE) $(OBJ) $(EXE).o $(LIBS)

doc: $(CPPSRC) doc-config
	doxygen doc-config

.PHONY: clean
clean:
	rm -f $(OBJ) $(EXE).o
