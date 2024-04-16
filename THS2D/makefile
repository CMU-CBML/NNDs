UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	Compiler=g++
endif
ifeq ($(UNAME_S),Darwin)
	Compiler=/opt/homebrew/Cellar/gcc/13.2.0/bin/aarch64-apple-darwin23-g++-13
endif

EIGEN_DIR = ../eigen-3.4.0

CFLAGS= ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS}  ${PSOURCECXX} --fast-math -O3 -std=c++11 -fopenmp -Wall -I ${EIGEN_DIR}

all: TTSP2D

TTSP2D: BasicDataStructure.o BSplineBasis.o kernel.o KnotInsertion.o Laplace.o LeastSquare.o main.o SingularEval.o TTSP_3D.o TTSP_2D.o
	$(Compiler) BasicDataStructure.o BSplineBasis.o kernel.o KnotInsertion.o Laplace.o LeastSquare.o main.o SingularEval.o TTSP_3D.o TTSP_2D.o -o TTSP2D ${PETSC_LIB} $(CFLAGS)
	
BasicDataStructure.o: BasicDataStructure.cpp
	$(Compiler)	-c BasicDataStructure.cpp $(CFLAGS)
	
BSplineBasis.o: BSplineBasis.cpp
	$(Compiler)	-c BSplineBasis.cpp $(CFLAGS)
	
Elasiticity_3D.o: Elasiticity_3D.cpp
	$(Compiler)	-c Elasiticity_3D.cpp $(CFLAGS)
	
kernel.o: kernel.cpp
	$(Compiler)	-c kernel.cpp $(CFLAGS)
	
KnotInsertion.o: KnotInsertion.cpp
	$(Compiler)	-c KnotInsertion.cpp $(CFLAGS)

Laplace.o: Laplace.cpp
	$(Compiler)	-c Laplace.cpp $(CFLAGS)

LeastSquare.o: LeastSquare.cpp
	$(Compiler)	-c LeastSquare.cpp $(CFLAGS)
	
main.o: main.cpp
	$(Compiler)	-c main.cpp $(CFLAGS)
	
SingularEval.o: SingularEval.cpp
	$(Compiler)	-c SingularEval.cpp $(CFLAGS)

TTSP_2D.o: TTSP_2D.cpp
	$(Compiler)	-c TTSP_2D.cpp $(CFLAGS)

TTSP_3D.o: TTSP_3D.cpp
	$(Compiler)	-c TTSP_3D.cpp $(CFLAGS)

clean:
	rm *o TTSP2D

