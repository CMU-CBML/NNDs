###################################################################################################	
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	PETSC_DIR = /home/kuanrenqian/Research/petsc
	Compiler=/home/kuanrenqian/Research/petsc/arch-linux-c-opt/bin/mpicxx
endif
ifeq ($(UNAME_S),Darwin)
	PETSC_DIR = /Users/kuanrenqian/petsc
	Compiler=/opt/homebrew/Cellar/gcc/13.2.0/bin/aarch64-apple-darwin23-g++-13
endif

EIGEN_DIR = ../eigen-3.4.0

# PETSC_DIR = /ocean/projects/eng170006p/ussqww/petsc
# Compiler = /ocean/projects/eng170006p/ussqww/petsc/arch-linux-c-opt/bin/mpicxx

include ${PETSC_DIR}/lib/petsc/conf/variables
CFLAGS= ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS} ${PSOURCECXX} -ffast-math -O3 -march=native -std=c++17 -fopenmp -Wall -I ${EIGEN_DIR}

all: 2DNG

2DNG: BasicDataStructure.o NeuronGrowth.o utils.o main.o 
	$(Compiler) BasicDataStructure.o NeuronGrowth.o utils.o main.o -o 2DNG ${PETSC_LIB} $(CFLAGS)
	
BasicDataStructure.o: BasicDataStructure.cpp
	$(Compiler)	-c BasicDataStructure.cpp $(CFLAGS)

NeuronGrowth.o: NeuronGrowth.cpp
	$(Compiler)	-c NeuronGrowth.cpp $(CFLAGS)

utils.o: utils.cpp
	$(Compiler)	-c utils.cpp $(CFLAGS)

main.o: main.cpp
	$(Compiler)	-c main.cpp $(CFLAGS)

clean:
	rm *o 2DNG