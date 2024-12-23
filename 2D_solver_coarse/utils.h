#ifndef UTILS_H
#define UTILS_H

#include <array>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <vector>
#include "BasicDataStructure.h"

using namespace std;

// Creating 2D mesh (incrementing from lo to hi)
void gen2Dmesh(int Nx, int Ny, vector<vector<float>>& vertices, vector<vector<int>>& elements);
void gen2Dmesh(int originX, int originY, int Nx, int Ny, vector<vector<float>>& vertices, vector<vector<int>>& elements, float ratio);
void gen2Dmesh_new(int originX, int originY, int Nx, int Ny, vector<vector<float>>& vertices, vector<vector<int>>& elements, float ratio);
// Export quad mesh to vtk for visualization
void write_quad_toVTK(const char* qs, vector<vector<float>>& vertices, vector<vector<int>>& elements);

// Creating 3D mesh (incrementing from lo to hi)
void gen3Dmesh(int Nx, int Ny, int Nz, vector<vector<float>>& vertices, vector<vector<int>>& elements);
// Export hex mesh to vtk for visualization
void write_hex_toVTK(const char* qs, vector<vector<float>>& vertices, vector<vector<int>>& elements);

void PrintVec2TXT(const vector<float>& v, string fn, bool visualization); // print out to commandline for debugging

// generating bezier mesh
void bzmesh2D(const string& path_in);
void bzmesh3D();

// partitioning mesh using mpmetis
void mpmetis(int n_process, const string& path_in);

void THS2D(const string& path_in, const vector<int>& rfid, const vector<int>& rftype);

void InitializeSoma(const int& numNeuron, vector<array<float, 2>> &seed, int &NX, int &NY);
void InitializeSoma_customizedCases(int& numNeuron, vector<array<float, 2>> &seed, int &NX, int &NY);
void InitializeRandomSoma(const int& numNeuron, vector<array<float, 2>>& seed, int& NX, int& NY);

vector<float> ConvertTo1DFloatVector(const vector<vector<int>>& input);
vector<float> ConvertTo1DFloatVector(const vector<vector<float>>& input);
vector<float> ConvertTo1DFloatVector(const vector<int>& input);
vector<int> ConvertTo1DIntVector(const vector<float>& input);

bool SearchPair(const vector<Vertex2D> prev_cpts, float targetX, float targetY, int &ind);
vector<float> InterpolateVars(vector<vector<int>> input, vector<Vertex2D> cpts_initial, vector<Vertex2D> cpts, int type);
vector<float> InterpolateVars(vector<float> input, vector<Vertex2D> cpts_initial, vector<Vertex2D> cpts, int type);
float round5(float input);
vector<float> InterpolateVars_coarse(vector<float> input, vector<Vertex2D> cpts_initial, const vector<Vertex2D>& cpts, int type);
vector<float> interpolateValues_closest(const vector<float>& phi, const vector<Vertex2D>& cpt, const vector<Vertex2D>& cpt_out);

float distanceTo(const Vertex2D& a, const Vertex2D& b);

void ObtainRefineID_coarse(vector<float> phi, vector<Vertex2D> cpts, int NX, int NY, int originX, int originY, vector<int> &rfid, vector<int> &rftype);

void ReadMesh(string fn, vector<Vertex2D>& pts, vector<Element2D>& mesh); //need vtk file with point label

void AssignProcessor(string fn, int &n_bzmesh, vector<vector<int>> &ele_process);

#endif