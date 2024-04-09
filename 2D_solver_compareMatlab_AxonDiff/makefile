###################################################################################################	
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	PETSC_DIR = /home/kuanrenqian/Research/petsc
	EIGEN_DIR = /home/kuanrenqian/Research/eigen-3.4.0
	Compiler=/home/kuanrenqian/Research/petsc/arch-linux-c-opt/bin/mpicxx
endif
ifeq ($(UNAME_S),Darwin)
	PETSC_DIR = /Users/kuanrenqian/petsc
	EIGEN_DIR = /Users/kuanrenqian/eigen-3.4.0
	Compiler=/opt/homebrew/Cellar/gcc/13.2.0/bin/aarch64-apple-darwin23-g++-13
endif

# PETSC_DIR = /ocean/projects/eng170006p/ussqww/petsc
# EIGEN_DIR = /ocean/projects/eng170006p/ussqww/eigen-3.4.0
# Compiler = /ocean/projects/eng170006p/ussqww/petsc/arch-linux-c-opt/bin/mpicxx

include ${PETSC_DIR}/lib/petsc/conf/variables
CFLAGS= ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS} ${PSOURCECXX} -ffast-math -O3 -march=native -std=c++17 -fopenmp -Wall -lpng -I ${EIGEN_DIR}

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
	
# UNAME_S := $(shell uname -s)
# ifeq ($(UNAME_S),Linux)
# 	PETSC_DIR = /home/kuanrenqian/Research/petsc
# 	EIGEN_DIR = /home/kuanrenqian/Research/eigen-3.4.0
# 	Compiler=/home/kuanrenqian/Research/petsc/arch-linux-c-opt/bin/mpicxx
# endif
# ifeq ($(UNAME_S),Darwin)
# 	PETSC_DIR = /Users/kuanrenqian/petsc
# 	EIGEN_DIR = /Users/kuanrenqian/eigen-3.4.0
# 	Compiler=/opt/homebrew/Cellar/gcc/13.2.0/bin/aarch64-apple-darwin23-g++-13
# endif

# include ${PETSC_DIR}/lib/petsc/conf/variables

# CFLAGS= ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS} ${PSOURCECXX} -ffast-math -O3 -march=native -std=c++11 -fopenmp -Wall -I ${EIGEN_DIR}

# all: 2DNG

# 2DNG: SimpleBasicDataStructure.o UserSetting.o NeuronGrowth.o utils.o BasicDataStructure.o BSplineBasis.o kernel.o KnotInsertion.o Laplace.o LeastSquare.o SingularEval.o TTSP_3D.o TTSP_2D.o main.o 
# 	$(Compiler) SimpleBasicDataStructure.o UserSetting.o NeuronGrowth.o utils.o BasicDataStructure.o BSplineBasis.o kernel.o KnotInsertion.o Laplace.o LeastSquare.o SingularEval.o TTSP_3D.o TTSP_2D.o main.o -o 2DNG ${PETSC_LIB} $(CFLAGS)
	
# SimpleBasicDataStructure.o: SimpleBasicDataStructure.cpp
# 	$(Compiler)	-c SimpleBasicDataStructure.cpp $(CFLAGS)
	
# UserSetting.o: UserSetting.cpp
# 	$(Compiler)	-c UserSetting.cpp $(CFLAGS)
	
# NeuronGrowth.o: NeuronGrowth.cpp
# 	$(Compiler)	-c NeuronGrowth.cpp $(CFLAGS)

# utils.o: utils.cpp
# 	$(Compiler)	-c utils.cpp $(CFLAGS)

# main.o: main.cpp
# 	$(Compiler)	-c main.cpp $(CFLAGS)

# clean:
# 	rm *o 2DNG ../io/* ../io/outputs/*vtk

# # all: THS3D

# # THS3D: BasicDataStructure.o BSplineBasis.o kernel.o KnotInsertion.o Laplace.o LeastSquare.o main.o SingularEval.o TTSP_3D.o TTSP_2D.o
# # 	$(Compiler) BasicDataStructure.o BSplineBasis.o kernel.o KnotInsertion.o Laplace.o LeastSquare.o main.o SingularEval.o TTSP_3D.o TTSP_2D.o -o THS3D ${PETSC_LIB} $(CFLAGS)
	
# BasicDataStructure.o: ../THS3D/BasicDataStructure.cpp
# 	$(Compiler)	-c ../THS3D/BasicDataStructure.cpp $(CFLAGS)
	
# BSplineBasis.o: ../THS3D/BSplineBasis.cpp
# 	$(Compiler)	-c ../THS3D/BSplineBasis.cpp $(CFLAGS)
	
# Elasiticity_3D.o: ../THS3D/Elasiticity_3D.cpp
# 	$(Compiler)	-c ../THS3D/Elasiticity_3D.cpp $(CFLAGS)
	
# kernel.o: ../THS3D/kernel.cpp
# 	$(Compiler)	-c ../THS3D/kernel.cpp $(CFLAGS)
	
# KnotInsertion.o: ../THS3D/KnotInsertion.cpp
# 	$(Compiler)	-c ../THS3D/KnotInsertion.cpp $(CFLAGS)

# Laplace.o: ../THS3D/Laplace.cpp
# 	$(Compiler)	-c ../THS3D/Laplace.cpp $(CFLAGS)

# LeastSquare.o: ../THS3D/LeastSquare.cpp
# 	$(Compiler)	-c ../THS3D/LeastSquare.cpp $(CFLAGS)
	
# # main.o: ../THS3D/main.cpp
# # 	$(Compiler)	-c ../THS3D/main.cpp $(CFLAGS)
	
# SingularEval.o: ../THS3D/SingularEval.cpp
# 	$(Compiler)	-c ../THS3D/SingularEval.cpp $(CFLAGS)

# TTSP_2D.o: ../THS3D/TTSP_2D.cpp
# 	$(Compiler)	-c ../THS3D/TTSP_2D.cpp $(CFLAGS)

# TTSP_3D.o: ../THS3D/TTSP_3D.cpp
# 	$(Compiler)	-c ../THS3D/TTSP_3D.cpp $(CFLAGS)