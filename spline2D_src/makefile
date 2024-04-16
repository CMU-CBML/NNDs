UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	Compiler=g++
endif
ifeq ($(UNAME_S),Darwin)
	Compiler=/opt/homebrew/Cellar/gcc/13.2.0/bin/aarch64-apple-darwin23-g++-13
endif

EIGEN_DIR = ../eigen-3.4.0

CFLAGS= --fast-math -O3 -c -std=c++11 -I ${EIGEN_DIR} -Wall 
#-DMKL_ILP64
FLAGS_LINK = -fopenmp -lpthread -lm -ldl

all: spline

spline: BasicDataStructure.o BSplineBasis.o CreateProblem.o KnotInsertion.o Laplace.o LeastSquare.o mesh.o SingularEval.o Truncated_Tspline.o main.o 
	$(Compiler) BasicDataStructure.o BSplineBasis.o CreateProblem.o KnotInsertion.o Laplace.o LeastSquare.o mesh.o SingularEval.o Truncated_Tspline.o main.o  -o spline $(FLAGS_LINK)

BasicDataStructure.o: BasicDataStructure.cpp
	$(Compiler)	$(CFLAGS) BasicDataStructure.cpp $(FLAGS_LINK)
	
BSplineBasis.o: BSplineBasis.cpp
	$(Compiler)	$(CFLAGS) BSplineBasis.cpp $(FLAGS_LINK)

CreateProblem.o: CreateProblem.cpp
	$(Compiler)	$(CFLAGS) CreateProblem.cpp $(FLAGS_LINK)

KnotInsertion.o: KnotInsertion.cpp
	$(Compiler)	$(CFLAGS) KnotInsertion.cpp $(FLAGS_LINK)

Laplace.o: Laplace.cpp
	$(Compiler)	$(CFLAGS) Laplace.cpp $(FLAGS_LINK)

LeastSquare.o: LeastSquare.cpp
	$(Compiler)	$(CFLAGS) LeastSquare.cpp $(FLAGS_LINK)

mesh.o: mesh.cpp
	$(Compiler)	$(CFLAGS) mesh.cpp $(FLAGS_LINK)

SingularEval.o: SingularEval.cpp
	$(Compiler)	$(CFLAGS) SingularEval.cpp $(FLAGS_LINK)

Truncated_Tspline.o: Truncated_Tspline.cpp
	$(Compiler)	$(CFLAGS) Truncated_Tspline.cpp $(FLAGS_LINK)

main.o: main.cpp
	$(Compiler)	$(CFLAGS) main.cpp $(FLAGS_LINK)

clean:
	rm *o ./spline

