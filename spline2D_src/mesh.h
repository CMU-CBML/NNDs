#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <utility>

using namespace std;

void CreateHemiSphere(string fn);

void CreateCylinder(double radius, double len, int nel_u, int nel_v, string fn);

void CreateCylinder_Irr(string fn_in, string fn_out, double radius, double len);

void ReadAbaqusInp_Quad(string fn, vector<array<double,3>>& pts, vector<array<int,4>>& cnct, 
	vector<int>& indxp, vector<int>& indxe);

void WriteVtk_Quad(string fn, const vector<array<double, 3>>& pts, const vector<array<int, 4>>& cnct);

void Inp2Vtk_Quad(string fn_in, string fn_out);

void ReadVtk_Quad(string fn, vector<array<double, 3>>& pts, vector<array<int, 4>>& cnct);

void RemoveDuplicatePoints(string fn_in, string fn_out);

#endif
