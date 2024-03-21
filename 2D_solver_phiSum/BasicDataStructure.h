#ifndef BASIC_DATA_STRUCTURE_H
#define BASIC_DATA_STRUCTURE_H

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include <petsc.h>
#include <petscksp.h>

#include "petscsys.h"   
#include "petscmat.h"
using namespace std;


const int bzpt_num = 16;
const int degree = 3;
const int dim = 2;
/////////////////////////////////
class Vertex2D
{
public:
	float coor[2];
	int label; //ID for inlet and outlet
	Vertex2D();
};


class Element2D
{
public:
	int degree;
	int order;
	int nbf;
	int type;//0 for interior and 1 for boundary, for visualization purpose
	int bzflag;//0 for spline element, 1 for Bezier element

	vector<int> IEN;
	vector<int> IENb;
	vector<array<float, 16>> cmat;
	vector<array<float, 2>> pts;//tmp
	vector<int> BC_order; // save order for boundary condition
	vector<float> ele_mat_bcvalue; // save the matrix column value for bc pts
	float velocity[2];

	Element2D(int p = 3);
	void BezierPolyn(float u, vector<float>& Nu, vector<float>& dNdu) const;
	void Basis(float u, float v, vector<float>& Nt, vector<array<float, 2>>& dNdt) const;
	void Para2Phys(float u, float v, float pt[2]) const;
	
};

//mesh format conversion

void Raw2Vtk_hex(string fn);

void Rawn2Vtk_hex(string fn);

#endif