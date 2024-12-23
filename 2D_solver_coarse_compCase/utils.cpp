#include "utils.h"
#include <iostream>
#include <algorithm>
#include "BasicDataStructure.h"
#include <cmath>

#include <random>	// for random generator

#include <queue>

// Creating 2D mesh (incrementing from lo to hi)
void gen2Dmesh(int Nx, int Ny, vector<vector<float>>& vertices, vector<vector<int>>& elements)
{
	cout << "******************************************************************************" << endl;
	cout << "Generating 2D mesh (Nx by Ny): " << Nx << " x " << Ny << endl;
	vertices.clear(); elements.clear();
	vector<float> tmp_vtx;
	vector<int> tmp_ele;
	for(int i = 0; i <= Nx; i++)
	{
		for(int j = 0; j <= Ny; j++)
		{
			tmp_vtx.clear();
			tmp_vtx.push_back(i);
			tmp_vtx.push_back(j);
			// tmp_vtx.push_back((float)(i)/(float)Nx);
			// tmp_vtx.push_back((float)(j)/(float)Ny);
			tmp_vtx.push_back(0);
			vertices.push_back(tmp_vtx);
		}
	}

	int tl_pt;
	for(int i = 0; i < Nx; i++)
	{
		for(int j = 0; j < Ny; j++)
		{
			tl_pt = i*(Ny+1)+j;
			tmp_ele.clear();
			tmp_ele.push_back(tl_pt);
			tmp_ele.push_back(tl_pt+1);
			tmp_ele.push_back(tl_pt+Ny+2);
			tmp_ele.push_back(tl_pt+Ny+1);
			elements.push_back(tmp_ele);
		}
	}
}

void gen2Dmesh(int originX, int originY, int Nx, int Ny, vector<vector<float>>& vertices, vector<vector<int>>& elements, float ratio)
{
	cout << "******************************************************************************" << endl;
	cout << "Generating 2D structured initial mesh" << endl;
	cout << "-----------------------------------------------------------------------------" << endl;
	cout << "Nx by Ny: " << Nx << " x " << Ny << " | origin: " << originX << "," << originY << endl;
	vertices.clear(); elements.clear();
	vector<float> tmp_vtx;
	vector<int> tmp_ele;
	for(int i = 0; i <= Nx; i++)
	{
		for(int j = 0; j <= Ny; j++)
		{
			tmp_vtx.clear();
			tmp_vtx.push_back(i*2/ratio+originX);
			tmp_vtx.push_back(j*2/ratio+originY);
			// tmp_vtx.push_back((float)(i)/(float)Nx);
			// tmp_vtx.push_back((float)(j)/(float)Ny);
			tmp_vtx.push_back(0);
			vertices.push_back(tmp_vtx);
		}
	}

	int tl_pt;
	for(int i = 0; i < Nx; i++)
	{
		for(int j = 0; j < Ny; j++)
		{
			tl_pt = i*(Ny+1)+j;
			tmp_ele.clear();
			tmp_ele.push_back(tl_pt);
			tmp_ele.push_back(tl_pt+1);
			tmp_ele.push_back(tl_pt+Ny+2);
			tmp_ele.push_back(tl_pt+Ny+1);
			elements.push_back(tmp_ele);
		}
	}
}

void gen2Dmesh_new(int originX, int originY, int Nx, int Ny, vector<vector<float>>& vertices, vector<vector<int>>& elements, float ratio)
{
	cout << "******************************************************************************" << endl;
	cout << "Generating 2D structured initial mesh" << endl;
	cout << "-----------------------------------------------------------------------------" << endl;
	cout << "Nx by Ny: " << Nx << " x " << Ny << " | origin: " << originX << "," << originY << endl;
	vertices.clear(); elements.clear();
	vector<float> tmp_vtx;
	vector<int> tmp_ele;
	for(int i = originX+1; i < (originX + Nx); i++)
	{
		for(int j = originY+1; j < (originY + Ny); j++)
		{
			tmp_vtx.clear();
			tmp_vtx.push_back(i*2/ratio);
			tmp_vtx.push_back(j*2/ratio);
			// tmp_vtx.push_back((float)(i)/(float)Nx);
			// tmp_vtx.push_back((float)(j)/(float)Ny);
			tmp_vtx.push_back(0);
			vertices.push_back(tmp_vtx);
		}
	}

	int tl_pt;
	for(int i = 0; i < Nx-1; i++)
	{
		for(int j = 0; j < Ny-1; j++)
		{
			tl_pt = i*(Ny+1)+j;
			tmp_ele.clear();
			tmp_ele.push_back(tl_pt);
			tmp_ele.push_back(tl_pt+1);
			tmp_ele.push_back(tl_pt+Ny+2);
			tmp_ele.push_back(tl_pt+Ny+1);
			elements.push_back(tmp_ele);
		}
	}
}

// Export quad mesh to vtk for visualization
void write_quad_toVTK(const char* qs, vector<vector<float>>& vertices, vector<vector<int>>& elements)
{
	FILE* fp;
	fp = fopen(qs, "w");
	int nv, nquad, i, j;

	nv = vertices.size();
	nquad = elements.size();

	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "2DmeshGen\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

	fprintf(fp, "POINTS %d float\n", nv);
	for (i = 0; i < nv; i++) {
		fprintf(fp, "%.4f %.4f %.4f\n", vertices[i][0], vertices[i][1], vertices[i][2]);
	}

	fprintf(fp, "\nCELLS %d %d\n", nquad, nquad * 5);

	for (i = 0; i < nquad; i++) {
		fprintf(fp, "4 %d %d %d %d\n", elements[i][0], elements[i][1], elements[i][2], elements[i][3]);
	}

	fprintf(fp, "\nCELL_TYPES %d\n", nquad);
	for (i = 0; i < nquad; i++) {
		fprintf(fp, "9\n");

	}
	fclose(fp);
}

// Creating 2D mesh (incrementing from lo to hi)
void gen3Dmesh(int Nx, int Ny, int Nz, vector<vector<float>>& vertices, vector<vector<int>>& elements)
{
	// float x_coord,y_coord;
	// vector<float,2> pt_coord;
	vector<float> tmp_vtx;
	vector<int> tmp_ele;
	for(int i = 0; i < Nx; i++)
	{
		for(int j = 0; j < Ny; j++)
		{
			for(int k = 0; k < Nz; k++)
			{
				// pt_coord = {i,j};
				tmp_vtx.clear();
				tmp_vtx.push_back(i);
				tmp_vtx.push_back(j);
				// tmp_vtx.push_back((float)(i+1)/(float)Nx);
				// tmp_vtx.push_back((float)(j+1)/(float)Ny);
				tmp_vtx.push_back(1);
				vertices.push_back(tmp_vtx);
			}
		}
	}

	int tl_pt;
	// cout << "Testing tl_pt" << endl;
	for(int i = 0; i < Nx-1; i++)
	{
		for(int j = 0; j < Ny-1; j++)
		{
			for(int k = 0; k < Nz-1; k++)
			{
				tl_pt = i*(Nx)+j*(Ny)+k;
				// cout<< tl_pt << endl;
				tmp_ele.clear();
				tmp_ele.push_back(tl_pt);
				tmp_ele.push_back(tl_pt+1);
				tmp_ele.push_back(tl_pt+Nx+1);
				tmp_ele.push_back(tl_pt+Nx);
				tmp_ele.push_back(tl_pt+Ny);
				tmp_ele.push_back(tl_pt+Ny+1);
				tmp_ele.push_back(tl_pt+Ny+Nx+1);
				tmp_ele.push_back(tl_pt+Ny+Nx);
				elements.push_back(tmp_ele);
			}
		}
	}

}

void PrintVec2TXT(const vector<float>& v, string fn, bool visualization)
{
	ofstream fout;

	fout.open(fn);

	fout << setprecision(2) << fixed;

	int v_size = v.size();
	if (visualization == 0) {
		for (size_t i = 0; i < v.size(); i++) {
			fout << v[i] << endl;
		}
	} else {
		int sq_sz = (int)sqrt(v_size);
		int ind = 0;
		for (size_t i = 0; i < sq_sz; i++) {
			for (size_t j = 0; j < sq_sz; j++) {
				if (round(v[ind]) == 0) {
					fout << "     ";
				} else {
					fout << round(v[ind]) << " ";
				}
				ind += 1;
			}
		fout << endl;
		}
	}
	fout.close();
}

// Export hex mesh to vtk for visualization
void write_hex_toVTK(const char* qs, vector<vector<float>>& vertices, vector<vector<int>>& elements)
{
	FILE* fp;
	fp = fopen(qs, "w");
	int nv, nquad, i, j;

	nv = vertices.size();
	nquad = elements.size();

	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "2DmeshGen\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

	fprintf(fp, "POINTS %d float\n", nv);
	for (i = 0; i < nv; i++) {
		fprintf(fp, "%.2f %.2f %.2f\n", vertices[i][0], vertices[i][1], vertices[i][2]);
	}

	fprintf(fp, "\nCELLS %d %d\n", nquad, nquad * 5);

	for (i = 0; i < nquad; i++) {
		fprintf(fp, "4 %d %d %d %d\n", elements[i][0], elements[i][1], elements[i][2], elements[i][3]);
	}

	fprintf(fp, "\nCELL_TYPES %d\n", nquad);
	for (i = 0; i < nquad; i++) {
		fprintf(fp, "9\n");

	}
	fclose(fp);
}

// generating 2D bezier mesh using spline2D_src
void bzmesh2D(const string& path_in){
	cout << "******************************************************************************" << endl;
	string spline_cmd = "../spline2D_src/spline -i " + path_in;
	if (system(spline_cmd.c_str()) != 0) {
		cerr << "Command failed to execute." << endl;
	}
}

// generating 3D bezier mesh using spline_src
void bzmesh3D(){
	cout << "******************************************************************************" << endl;
	string spline_cmd = "../NeuronTransportIGA/spline_src/spline ../io/3DNG/";
	if (system(spline_cmd.c_str()) != 0) {
		cerr << "Command failed to execute." << endl;
	}
}

// partitioning mesh using mpmetis
void mpmetis(int n_process, const string& path_in){
	string mpmetis_cmd = "mpmetis " + path_in + "bzmeshinfo.txt " + to_string(n_process);
	if (system(mpmetis_cmd.c_str()) != 0) {
		cerr << "Command failed to execute." << endl;
	}
}

// // partitioning mesh using mpmetis
void THS2D(const string& path_in, const vector<int>& rfid, const vector<int>& rftype) {
	cout << "******************************************************************************" << endl;
	cout << "Local refinement based on Xiaodong's THS3D code ... " << endl;
	cout << "  - see: Truncated T-splines: Fundamentals and methods (2017)" << endl << endl;
	// cout << "-----------------------------------------------------------------------------" << endl;
	// cout << "Calling command | input mesh directory | refine ID | refine element type" << endl << endl;

	string ths2d_cmd_tmp("../THS2D/TTSP2D " + path_in + " ");
	for (size_t i = 0; i < rfid.size(); ++i) {
		ths2d_cmd_tmp += to_string(rfid[i]) + " ";
	}
	for (size_t i = 0; i < rftype.size(); ++i) {
		ths2d_cmd_tmp += to_string(rftype[i]) + " ";
	}
	// (Optional) output local refine element ID and type
	// cout << ths2d_cmd_tmp << endl;

	if (system(ths2d_cmd_tmp.c_str()) != 0) {
		cerr << "Command failed to execute." << endl;
	}
}

void InitializeSoma(const int& numNeuron, vector<array<float, 2>> &seed, int &NX, int &NY){
	seed.resize(numNeuron);
	// 2D neuron soma initialization
	switch (numNeuron) {
	case 1:
		NX = 30;
		NY = 30;
		seed = {{30, 30}};
		break;
        case 2:
		// NX = 100;
		// NY = 30;
		// seed = {{30, 30}, {170, 30}};
		NX = 50;
		NY = 30;
		seed = {{30, 30}, {70, 30}};
		break;
        // case 3:
	// 	NX = 100;
	// 	NY = 100;
	// 	seed = {{42, 56}, {151, 147}, {34, 158}};
	// 	break;
	// case 4:
	// 	NX = 100;
	// 	NY = 100;
	// 	seed = {{30, 30}, {170, 30}, {30, 170}, {170, 170}};
	// 	break;
	// case 5:
	// 	NX = 100;
	// 	NY = 100;
	// 	seed = {{30, 30}, {170, 30}, {30, 170}, {170, 170}, {100, 100}};
	// 	break;
	// case 6:
	// 	NX = 120;
	// 	NY = 120;
	// 	seed = {{30, 60}, {120, 30}, {210, 60}, {210, 180}, {120, 210}, {30, 180}};
	// 	break;
	// case 7:
	// 	NX = 120;
	// 	NY = 120;
	// 	seed = {{30, 60}, {120, 30}, {210, 60}, {210, 180}, {120, 210}, {30, 180}, {120, 120}};
	// 	break;
	}
}

// void InitializeSoma_customizedCases(int& numNeuron, vector<array<float, 2>> &seed, int &NX, int &NY){
// 	seed.resize(numNeuron);
// 	// 2D neuron soma initialization
// 	switch (numNeuron) {
// 	// case 1:
// 	// 	NX = 30;
// 	// 	NY = 30;
// 	// 	seed = {{30, 30}};
// 	// 	break;
//         // case 2: // 2 neursons
// 	// 	NX = 50;
// 	// 	NY = 50;
// 	// 	seed = {{32, 55}, {72, 30}};
// 	// 	break;
//         case 3: // 3 neurons - 1
// 		NX = 60;
// 		NY = 95;
// 		seed = {{98, 50}, {29, 166}, {80, 160}};
// 		break;
// 	case 4: // 3 neurons - 2 (just using numNeuron case easy case selection here)
// 		NX = 70;
// 		NY = 129;
// 		seed = {{55, 44}, {107, 130}, {34, 228}};
// 		numNeuron = 3;
// 		break;
// 	case 5: // 5 neurons
// 		NX = 110;
// 		NY = 200;
// 		seed = {{192, 32}, {35, 137}, {57, 239}, {190, 208}, {183, 376}};
// 		break;
// 	}
// }

void InitializeSoma_customizedCases(const string& caseType, vector<array<float, 2>> &seed, int &NX, int &NY){
	// 2D neuron soma initialization
	if (caseType == "K") {
		NX = 50;
		NY = 50;
		seed = {{40, 65}, {72, 30}};
	} else if (caseType == "L") {
		NX = 100;
		NY = 50;
		seed = {{30, 72}, {155, 32}};
	} else if (caseType == "M") {
		NX = 50;
		NY = 120;
		seed = {{30, 200}, {75, 32}};
	} else if (caseType == "N") {
		NX = 130;
		NY = 65;
		seed = {{30, 50}, {125, 100}, {220, 30}};
	} else if (caseType == "O") {
		NX = 120;
		NY = 80;
		seed = {{30, 130}, {130, 30}, {140, 85}};
	}
	std::cout << "Running comparison case: " << caseType << std::endl;
}

void InitializeRandomSoma(const int& numNeuron, vector<array<float, 2>>& seed, int& NX, int& NY) {
	seed.clear(); // Clear existing seeds if any
	random_device rd;
	default_random_engine generator(rd());

	// Set domain sizes based on the number of neurons
	switch (numNeuron) {
	// case 1:
	// 	seed = {{30, 30}}; // no need for random placement
	// 	NX = 30; NY = 30; break;
	// case 2:
	// 	seed = {{30, 30}, {170, 30}}; // no need for random placement
	// 	NX = 100; NY = 30; break;
	case 3:
		NX = 100; NY = 100; break;
	case 4:
		NX = 100; NY = 100; break;
	case 5:
		NX = 120; NY = 100; break;
	case 6:
		NX = 120; NY = 120; break;
	case 7:
		NX = 120; NY = 120; break;
	default:
		throw runtime_error("Unsupported number of neurons");
	}

	uniform_real_distribution<float> dist_x(30, 2 * NX - 30);
	uniform_real_distribution<float> dist_y(30, 2 * NY - 30);
	int max_attempts = 1000;

	for (int i = 0; i < numNeuron; ++i) {
		int attempts = 0;
		bool valid;
		array<float, 2> newSeed;
		do {
			if (attempts++ > max_attempts) {
				throw runtime_error("Unable to place all neurons without overlapping after maximum attempts.");
			}
			newSeed = {dist_x(generator), dist_y(generator)};
			valid = true;
			for (auto& s : seed) {
				float dx = s[0] - newSeed[0];
				float dy = s[1] - newSeed[1];
				if (dx * dx + dy * dy < 50*50) {
					valid = false;
					break;
				}
			}
		} while (!valid);
		seed.push_back(newSeed);
	}

	// Adjust seed positions to minimize empty space
	float min_x = numeric_limits<float>::max(), min_y = numeric_limits<float>::max();
	float max_x = 0, max_y = 0;
	for (auto& s : seed) {
		min_x = min(min_x, s[0]);
		min_y = min(min_y, s[1]);
		max_x = max(max_x, s[0]);
		max_y = max(max_y, s[1]);
	}

	// Reposition seeds to start at least 30 units from the left and bottom
	float shift_x = min_x - 30;
	float shift_y = min_y - 30;
	for (auto& s : seed) {
		s[0] -= shift_x;
		s[1] -= shift_y;
	}

	// Recalculate NX and NY based on the maximum coordinates
	NX = static_cast<int>((max_x - shift_x + 30)/2); // Add margin and convert back
	NY = static_cast<int>((max_y - shift_y + 30)/2); // Add margin and convert back
}

vector<float> ConvertTo1DFloatVector(const vector<vector<int>>& input) {
	vector<float> output;

	// Calculate total size needed
	size_t totalSize = 0;
	for (const auto& row : input) {
		totalSize += row.size();
	}

	// Reserve memory to prevent multiple reallocations
	output.reserve(totalSize);

	// Convert all integers to floats and add to the output vector
	for (const auto& row : input) {
		for (int value : row) {
			output.push_back(static_cast<float>(value));
		}
	}

	return output;
}

vector<float> ConvertTo1DFloatVector(const vector<vector<float>>& input) {
	vector<float> output;

	// Calculate total size needed
	size_t totalSize = 0;
	for (const auto& row : input) {
		totalSize += row.size();
	}

	// Reserve memory to prevent multiple reallocations
	output.reserve(totalSize);

	// Convert all integers to floats and add to the output vector
	for (const auto& row : input) {
		for (int value : row) {
			output.push_back(static_cast<float>(value));
		}
	}

	return output;
}

// Function to convert a vector of floats to a vector of integers
vector<float> ConvertTo1DFloatVector(const vector<int>& input) {
	vector<float> output;
	output.reserve(input.size());  // Pre-allocate memory for efficiency

	// Convert each float to an int by rounding
	for (float value : input) {
		output.push_back(static_cast<float>(round(value)));
	}

	return output;
}

// Function to convert a vector of floats to a vector of integers
vector<int> ConvertTo1DIntVector(const vector<float>& input) {
	vector<int> output;
	output.reserve(input.size());  // Pre-allocate memory for efficiency

	// Convert each float to an int by rounding
	for (float value : input) {
		output.push_back(static_cast<int>(round(value)));
	}

	return output;
}

// Function to search for a particular x and y in the vector of Vertex2D
bool SearchPair(const vector<Vertex2D> prev_cpts, float targetX, float targetY, int &ind) {
	float x, y;
	for (size_t i = 0; i < prev_cpts.size(); i++) {

		x = prev_cpts[i].coor[0];
		if (abs(remainder(x,1)) > 0.5) {
			x = floorf(x)+1;
		}
		if (abs(remainder(x,1)) <= 0.5) {
			x = floorf(x);
		}
		y = prev_cpts[i].coor[1];
		if (abs(remainder(y,1)) > 0.5) {
			y = floorf(y)+1;
		}
		if (abs(remainder(y,1)) <= 0.5) {
			y = floorf(y);
		}

		if (x == targetX && y == targetY) {
			ind = i;
			return true; // Found the pair (targetX, targetY) in the vector
		}
	}

	return false; // Pair not found in the vector
}

vector<float> InterpolateVars(vector<vector<int>> input, vector<Vertex2D> cpts_initial, vector<Vertex2D> cpts, int type) 
{	
	vector<float> output;
	output.resize(cpts.size());
	vector<float> tmp = ConvertTo1DFloatVector(input);

	for (size_t i = 0; i < cpts.size(); i++) {
		// float x = cpts[i].coor[0];
		// float y = cpts[i].coor[1];

		float x = cpts[i].coor[0];
		if (abs(remainder(x,1)) != 0.5) {
			x = round(x);
		}
		float y = cpts[i].coor[1];
		if (abs(remainder(y,1)) != 0.5) {
			y = round(y);
		}

		int ind;
		if (SearchPair(cpts_initial, round(x), round(y), ind)) {
		// if (SearchPair(cpts_initial, x, y, ind)) {
			output[i] = tmp[ind];
		} 
		else {
			int indDown, indUp, indLeft, indRight;
			if ((abs(remainder(x,1)) == 0.5) && (abs(remainder(y,1)) != 0.5)) {
				if (SearchPair(cpts_initial, floorf(x), round(y), indDown) &&
					SearchPair(cpts_initial, floorf(x)+1, round(y), indUp)) {
					if (type == 0) {
						output[i] = max(tmp[indDown], tmp[indUp]);
					} else if (type == 1) {
						output[i] = (tmp[indDown] + tmp[indUp])/2;
					} else if (type == 2) {
						output[i] = 0;
					}
				}
				else {
					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck0)!\n");
				}
			} else if ((abs(remainder(x,1)) != 0.5) && (abs(remainder(y,1)) == 0.5)) {
				if (SearchPair(cpts_initial, round(x), floorf(y), indLeft) &&
					SearchPair(cpts_initial, round(x), floorf(y)+1, indRight)) {
					if (type == 0) {
						output[i] = max(tmp[indLeft], tmp[indRight]);
					} else if (type == 1) {
						output[i] = (tmp[indLeft] + tmp[indRight])/2;
					} else if (type == 2) {
						output[i] = 0;
					}
				} 
				else {
					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck1)!\n");
				}
			} else if ((abs(remainder(x,1)) == 0.5) && (abs(remainder(y,1)) == 0.5)) {
				if (SearchPair(cpts_initial, floorf(x), floorf(y), indDown) &&
					SearchPair(cpts_initial, floorf(x)+1, floorf(y), indUp) &&
					SearchPair(cpts_initial, floorf(x), floorf(y), indLeft) &&
					SearchPair(cpts_initial, floorf(x), floorf(y)+1, indRight)) {
					if (type == 0) {
						output[i] = max(max(tmp[indDown], tmp[indUp]), max(tmp[indLeft], tmp[indRight]));
					} else if (type == 1) {
						output[i] = (tmp[indDown] + tmp[indUp] + tmp[indLeft] + tmp[indRight])/4;
					} else if (type == 2) {
						output[i] = 0;
					}
				} 
				else {
					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck2)!\n");
				}
			}
		}			
	}
	return output;
}

vector<float> InterpolateVars(vector<float> input, vector<Vertex2D> cpts_initial, vector<Vertex2D> cpts, int type) 
{	
	vector<float> output;
	output.resize(cpts.size());

	for (size_t i = 0; i < cpts_initial.size(); i++) {
		if (abs(remainder(cpts_initial[i].coor[0],1)) != 0.5)
			cpts_initial[i].coor[0] = round(cpts_initial[i].coor[0]);
		if (abs(remainder(cpts_initial[i].coor[1],1)) != 0.5)
			cpts_initial[i].coor[1] = round(cpts_initial[i].coor[1]);
	}
	
	for (size_t i = 0; i < cpts.size(); i++) {
		// float x = cpts[i].coor[0];
		// float y = cpts[i].coor[1];

		float x = cpts[i].coor[0];
		if (abs(remainder(x,1)) != 0.5) {
			x = round(x);
		}
		float y = cpts[i].coor[1];
		if (abs(remainder(y,1)) != 0.5) {
			y = round(y);
		}

		int ind;
		if (SearchPair(cpts_initial, round(x), round(y), ind)) {
		// if (SearchPair(cpts_initial, x, y, ind)) {
			output[i] = input[ind];
		} 
		else {
			int indDown, indUp, indLeft, indRight;
			if ((abs(remainder(x,1)) == 0.5) && (abs(remainder(y,1)) != 0.5)) {
				if (SearchPair(cpts_initial, floorf(x), round(y), indDown) &&
					SearchPair(cpts_initial, floorf(x)+1, round(y), indUp)) {
					if (type == 0) {
						output[i] = max(input[indDown], input[indUp]);
					} else if (type == 1) {
						output[i] = (input[indDown] + input[indUp])/2;
					} else if (type == 2) {
						output[i] = 0;
					}
				} else {
					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck0)!\n");
				}
			} else if ((abs(remainder(x,1)) != 0.5) && (abs(remainder(y,1)) == 0.5)) {
				if (SearchPair(cpts_initial, round(x), floorf(y), indLeft) &&
					SearchPair(cpts_initial, round(x), floorf(y)+1, indRight)) {
					if (type == 0) {
						output[i] = max(input[indLeft], input[indRight]);
					} else if (type == 1) {
						output[i] = (input[indLeft] + input[indRight])/2;
					} else if (type == 2) {
						output[i] = 0;
					}
				} else {
					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck1)!\n");
				}
			} else if ((abs(remainder(x,1)) == 0.5) && (abs(remainder(y,1)) == 0.5)) {
				if (SearchPair(cpts_initial, floorf(x), floorf(y), indDown) &&
					SearchPair(cpts_initial, floorf(x)+1, floorf(y), indUp) &&
					SearchPair(cpts_initial, floorf(x), floorf(y), indLeft) &&
					SearchPair(cpts_initial, floorf(x), floorf(y)+1, indRight)) {
					if (type == 0) {
						output[i] = max(max(input[indDown], input[indUp]), max(input[indLeft], input[indRight]));
					} else if (type == 1) {
						output[i] = (input[indDown] + input[indUp] + input[indLeft] + input[indRight])/4;
					} else if (type == 2) {
						output[i] = 0;
					}
				} else {
					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck2)!\n");
				}
			}
		}			
	}
	return output;
}

float round5(float value) {
	// Extract the fractional part and the whole part of the value
	float wholePart = floor(value);
	float fractionalPart = value - wholePart;

	// Check if the fractional part is exactly 0.5
	if (fractionalPart == 0.5 || fractionalPart == -0.5) {
		return wholePart; // Round down to 0 for 0.5
	} else {
		return round(value); // Use standard rounding for other cases
	}
}

vector<float> InterpolateVars_coarse(vector<float> input, vector<Vertex2D> cpts_initial, const vector<Vertex2D>& cpts, int type) 
{	
	vector<float> output;
	output.resize(cpts.size(), 0);

	for (size_t i = 0; i < cpts_initial.size(); i++) {
		if (abs(remainder(cpts_initial[i].coor[0],2)) != 1)
			cpts_initial[i].coor[0] = round5(cpts_initial[i].coor[0]);
		if (abs(remainder(cpts_initial[i].coor[1],2)) != 1)
			cpts_initial[i].coor[1] = round5(cpts_initial[i].coor[1]);
	}
	
	for (size_t i = 0; i < cpts.size(); i++) {
		float x = round5(cpts[i].coor[0]);
		float y = round5(cpts[i].coor[1]);
		int ind;
		if (SearchPair(cpts_initial, x, y, ind)) {
			output[i] = input[ind];
		} 
		else {
			int indDown, indUp, indLeft, indRight;
			if ((abs(remainder(x,2)) == 1) && (abs(remainder(y,2)) != 1)) {
				if (SearchPair(cpts_initial, x-1, y, indDown) &&
					SearchPair(cpts_initial, x+1, y, indUp)) {
					if (type == 0) {
						output[i] = max(input[indDown], input[indUp]);
					} else if (type == 1) {
						output[i] = (input[indDown] + input[indUp])/2;
					} else if (type == 2) {
						output[i] = 0;
					}
				} else {
					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck0)! x: %f y: %f\n", x, y);
			}
			} else if ((abs(remainder(x,2)) != 1) && (abs(remainder(y,2)) == 1)) {
				if (SearchPair(cpts_initial, x, y-1, indLeft) &&
					SearchPair(cpts_initial, x, y+1, indRight)) {
					if (type == 0) {
						output[i] = max(input[indLeft], input[indRight]);
					} else if (type == 1) {
						output[i] = (input[indLeft] + input[indRight])/2;
					} else if (type == 2) {
						output[i] = 0;
					}
				} else {
					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck1)! x: %f y: %f\n", x, y);
				}
			} else if ((abs(remainder(x,2)) == 1) && (abs(remainder(y,2)) == 1)) {
				if (SearchPair(cpts_initial, x-1, y-1, indDown) &&
					SearchPair(cpts_initial, x+1, y+1, indUp) &&
					SearchPair(cpts_initial, x-1, y+1, indLeft) &&
					SearchPair(cpts_initial, x+1, y-1, indRight)) {
					if (type == 0) {
						output[i] = max(max(input[indDown], input[indUp]), max(input[indLeft], input[indRight]));
					} else if (type == 1) {
						output[i] = (input[indDown] + input[indUp] + input[indLeft] + input[indRight])/4;
					} else if (type == 2) {
						output[i] = 0;
					}
				} else {
					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck2)! x: %f y: %f | %f %f\n", x, y);
				}
			}
		}			
	}
	return output;
}

float distanceTo(const Vertex2D& a, const Vertex2D& b) {
	// return sqrt((a.coor[0] - b.coor[0]) * (a.coor[0] - b.coor[0]) + (a.coor[1] - b.coor[1]) * (a.coor[1] - b.coor[1]));
	return abs(a.coor[0] - b.coor[0]) + abs(a.coor[1] - b.coor[1]);
}

void ObtainRefineID_coarse(vector<float> phi, vector<Vertex2D> cpts, int NX, int NY, int originX, int originY, vector<int> &rfid, vector<int> &rftype){

	rfid.clear();
	rftype.clear();
	int ind = 0;
	// int offset = 1 + (NY+1); // offset for using phi on cpts to determine local refinements on elements
	int offset = 2 + 2*(NY+1); // offset for using phi on cpts to determine local refinements on elements

	for (size_t i = offset; i < cpts.size() - offset; i++) {
		float x = cpts[i].coor[0];
		float y = cpts[i].coor[1];
		// cout << x << " ";
		rfid.push_back(ind - round(x/2) + originX + offset); // push back element id (calculated based on vertex id, ind)
		// float averagePhi = (phi[i-(NY+1)+1] + phi[i-(NY+1)] + phi[i-(NY+1)-1] + phi[i-1] + phi[i] + phi[i+1] + phi[i+(NY+1)+1] + phi[i+(NY+1)] + phi[i+(NY+1)-1])/9;
		float averagePhi = (phi[i-2*(NY+1)+2] + phi[i-2*(NY+1)+1] + phi[i-2*(NY+1)] + phi[i-2*(NY+1)-1] + phi[i-2*(NY+1)-2] 
			+ phi[i-(NY+1)+2] + phi[i-(NY+1)+1] + phi[i-(NY+1)] + phi[i-(NY+1)-1] + phi[i-(NY+1)-2] 
			+ phi[i-2] + phi[i-1] + phi[i] + phi[i+1] + phi[i+2]
			+ phi[i+(NY+1)+2] + phi[i+(NY+1)+1] + phi[i+(NY+1)] + phi[i+(NY+1)-1] + phi[i+(NY+1)-2]
			+ phi[i+2*(NY+1)+2] + phi[i+2*(NY+1)+1] + phi[i+2*(NY+1)] + phi[i+2*(NY+1)-1] + phi[i+2*(NY+1)-2])/16;
		
		if ((averagePhi < 0.99) && (averagePhi > 0.01)) {
			rftype.push_back(0); // default to type 0 local refinement
		} else { // maintain domain element strcuture, needed for later local refinement limitation (face-face intersection)
			rftype.push_back(5); // 5 for elemetmpnts without local refinement
		}
		ind += 1;
	}

	// PrintVec2TXT(tmp, "../io2D/tmp_0.txt", 1);

	// handle local refinement limitations (no face-face intersection, check Xiaodong paper - truncated t-spline)
	int ring(1), rf1(1), rf2(1);
	switch (ring) {
		case 1: // apply supporting local refinement to 1-ring of elements
			for (size_t i = (NY+1)+1; i < rftype.size()-(NY+1)-1; i++) {
				if (rftype[i] == 0) {
					rftype[i-(NY+1)-1] = min(rftype[i-(NY+1)-1], rf1);
					rftype[i-(NY+1)] = min(rftype[i-(NY+1)], rf1);
					rftype[i-(NY+1)+1] = min(rftype[i-(NY+1)+1], rf1);
					rftype[i-1] = min(rftype[i-1], rf1);
					rftype[i+1] = min(rftype[i+1], rf2);
					rftype[i+(NY+1)-1] = min(rftype[i+(NY+1)-1], rf2);
					rftype[i+(NY+1)] = min(rftype[i+(NY+1)], rf2);
					rftype[i+(NY+1)+1] = min(rftype[i+(NY+1)+1], rf2);
				} 
			}
			break;
		case 2: // apply supporting local refinement to 2-ring of elements
			for (size_t i = 2*(NY+1)+2; i < rftype.size()-2*(NY+1)-2; i++) {
				if (rftype[i] == 0) {
					rftype[i-2*(NY+1)-2] = min(rftype[i-2*(NY+1)-2], rf1);
					rftype[i-2*(NY+1)-1] = min(rftype[i-2*(NY+1)-1], rf1);
					rftype[i-2*(NY+1)] = min(rftype[i-2*(NY+1)], rf1);
					rftype[i-2*(NY+1)+1] = min(rftype[i-2*(NY+1)+1], rf1);
					rftype[i-2*(NY+1)+2] = min(rftype[i-2*(NY+1)+2], rf1);
					
					rftype[i-(NY+1)-2] = min(rftype[i-(NY+1)-2], rf1);
					rftype[i-(NY+1)-1] = min(rftype[i-(NY+1)-1], rf1);
					rftype[i-(NY+1)] = min(rftype[i-(NY+1)], rf1);
					rftype[i-(NY+1)+1] = min(rftype[i-(NY+1)+1], rf1);
					rftype[i-(NY+1)+2] = min(rftype[i-(NY+1)+2], rf1);

					rftype[i-2] = min(rftype[i-2], rf1);
					rftype[i-1] = min(rftype[i-1], rf1);
					rftype[i+1] = min(rftype[i+1], rf2);
					rftype[i+2] = min(rftype[i+2], rf2);
					
					rftype[i+(NY+1)-2] = min(rftype[i+(NY+1)-2], rf2);
					rftype[i+(NY+1)-1] = min(rftype[i+(NY+1)-1], rf2);
					rftype[i+(NY+1)] = min(rftype[i+(NY+1)], rf2);
					rftype[i+(NY+1)+1] = min(rftype[i+(NY+1)+1], rf2);
					rftype[i+(NY+1)+2] = min(rftype[i+(NY+1)+2], rf2);

					rftype[i+2*(NY+1)-2] = min(rftype[i+2*(NY+1)-2], rf2);
					rftype[i+2*(NY+1)-1] = min(rftype[i+2*(NY+1)-1], rf2);
					rftype[i+2*(NY+1)] = min(rftype[i+2*(NY+1)], rf2);
					rftype[i+2*(NY+1)+1] = min(rftype[i+2*(NY+1)+1], rf2);
					rftype[i+2*(NY+1)+2] = min(rftype[i+2*(NY+1)+2], rf2);
				}
			}
	}

	for (size_t i = 0; i < rftype.size(); i++) {
		if ((rftype[i] == 5) && (rftype[i-1] != 5) && (rftype[i+1] != 5)) {
			rftype[i] = 1;
		}
	}

	for (size_t i = 0; i < rftype.size();) {
		if (rftype[i] == 5) {
			swap(rfid[i], rfid.back());
			rfid.pop_back();
			swap(rftype[i], rftype.back());
			rftype.pop_back();
		} else {
			i++;
		}
	}

	cout << endl;
}

void ReadMesh(string fn, vector<Vertex2D>& pts, vector<Element2D>& mesh)//need vtk file with point label
{
	string fname(fn), stmp;
	int npts, neles, itmp;
	double tmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (size_t i = 0; i < 4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		pts.resize(npts);
		for (size_t i = 0; i < npts; i++)
		{
			fin >> pts[i].coor[0] >> pts[i].coor[1] >> tmp; // pts[i].coor[2];
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		mesh.resize(neles);
		for (size_t i = 0; i < neles; i++)
		{
			fin >> itmp >> mesh[i].IEN[0] >> mesh[i].IEN[1] >> mesh[i].IEN[2] >> mesh[i].IEN[3]; /* >>
				mesh[i].IEN[4] >> mesh[i].IEN[5] >> mesh[i].IEN[6] >> mesh[i].IEN[7]; */
			for (size_t j = 0; j < 4; j++)
			{
				mesh[i].pts[j][0] = pts[mesh[i].IEN[j]].coor[0];
				mesh[i].pts[j][1] = pts[mesh[i].IEN[j]].coor[1];
				// mesh[i].pts[j][2] = pts[mesh[i].IEN[j]].coor[2];
			}

		}
		for (size_t i = 0; i < neles + 5; i++) getline(fin, stmp);//skip lines
		fin.close();
		PetscPrintf(PETSC_COMM_WORLD, "Mesh Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}

void AssignProcessor(string fn, int &n_bzmesh, vector<vector<int>> &ele_process)
{
	int tmp;
	int i = 0;
	string fname(fn);
	ifstream fin;
	fin.open(fname, ios::in);
	if (fin.is_open())
	{
		while (!fin.eof() && fin.peek() != EOF)
		{
			fin >> tmp;
			ele_process[tmp].push_back(i);
			i++;
			fin.get();
		}
		n_bzmesh = i;
		//PetscPrintf(PETSC_COMM_WORLD, "Mesh partition finished!\n");
		fin.close();
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}