#ifndef UTILS_H
#define UTILS_H

#include <array>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <vector>
#include "BasicDataStructure.h"

#include <png.h>
#include <stdexcept>
#include <iostream>

using namespace std;

// Creating 2D mesh (incrementing from lo to hi)
void gen2Dmesh(int Nx, int Ny, vector<vector<float>>& vertices, vector<vector<int>>& elements);
void gen2Dmesh(int originX, int originY, int Nx, int Ny, vector<vector<float>>& vertices, vector<vector<int>>& elements);

// Export quad mesh to vtk for visualization
void write_quad_toVTK(const char* qs, vector<vector<float>>& vertices, vector<vector<int>>& elements);

// Creating 3D mesh (incrementing from lo to hi)
void gen3Dmesh(int Nx, int Ny, int Nz, vector<vector<float>>& vertices, vector<vector<int>>& elements);

// Export hex mesh to vtk for visualization
void write_hex_toVTK(const char* qs, vector<vector<float>>& vertices, vector<vector<int>>& elements);

void PrintVec2TXT(vector<float>& v, string fn, bool visualization); // print out to commandline for debugging

// generating 2D bezier mesh using spline2D_src
void bzmesh2D(string path_in);
// generating 3D bezier mesh using spline_src
void bzmesh3D();

// partitioning mesh using mpmetis
void mpmetis(int n_process, string path_in);

void THS2D(string path_in, vector<int> rfid, vector<int> rftype);

void InitializeSoma(int numNeuron, vector<array<float, 2>> &seed, int &NX, int &NY);

// void InitializeNGvars(vector<Vertex2D> cpts, vector<array<int, 2>> seed, vector<vector<float>> &NGvars, vector<float> &phi_0, vector<float> &tub_0, vector<float> &Mphi, int n);

// Local refinement based on phi interface
vector<float> ConvertTo1DFloatVector(const vector<vector<int>> input);
vector<float> ConvertTo1DFloatVector(const vector<vector<float>> input);
bool SearchPair(const vector<Vertex2D> prev_cpts, float targetX, float targetY, int &ind);
vector<float> InterpolateVars(vector<vector<int>> input, vector<Vertex2D> cpts_initial, vector<Vertex2D> cpts, int type);
vector<float> InterpolateVars(vector<float> input, vector<Vertex2D> cpts_initial, vector<Vertex2D> cpts, int type);
void ObtainRefineID(vector<float> phi, vector<Vertex2D> cpts, int NX, int NY, int originX, int originY, vector<int> &rfid, vector<int> &rftype);

void VisualizeVTK_ControlMesh(const vector<Vertex2D> &spt, const vector<Element2D> &mesh, int step, string fn, vector<float> var, string varName);

void ReadMesh(string fn, vector<Vertex2D>& pts, vector<Element2D>& mesh); //need vtk file with point label

void AssignProcessor(string fn, int &n_bzmesh, vector<vector<int>> &ele_process);

void WriteMatrixToPNG(const std::vector<float>& phi, int width, int height, const char* filename);

#endif