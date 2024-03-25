#include "utils.h"
#include <iostream>
#include <algorithm>
#include "BasicDataStructure.h"
#include <cmath>

#include <png.h>
#include <stdexcept>
#include <iostream>

#include <queue>

// Creating 2D mesh (incrementing from lo to hi)
void gen2Dmesh(int Nx, int Ny, vector<vector<float>>& vertices, vector<vector<int>>& elements)
{
	std::cout << "******************************************************************************" << std::endl;
	std::cout << "Generating 2D mesh (Nx by Ny): " << Nx << " x " << Ny << std::endl;
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
	std::cout << "******************************************************************************" << std::endl;
	std::cout << "Generating 2D structured initial mesh" << std::endl;
	std::cout << "-----------------------------------------------------------------------------" << std::endl;
	std::cout << "Nx by Ny: " << Nx << " x " << Ny << " | origin: " << originX << "," << originY << std::endl;
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
	std::cout << "******************************************************************************" << std::endl;
	std::cout << "Generating 2D structured initial mesh" << std::endl;
	std::cout << "-----------------------------------------------------------------------------" << std::endl;
	std::cout << "Nx by Ny: " << Nx << " x " << Ny << " | origin: " << originX << "," << originY << std::endl;
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
	// std::cout << "Testing tl_pt" << std::endl;
	for(int i = 0; i < Nx-1; i++)
	{
		for(int j = 0; j < Ny-1; j++)
		{
			for(int k = 0; k < Nz-1; k++)
			{
				tl_pt = i*(Nx)+j*(Ny)+k;
				// std::cout<< tl_pt << std::endl;
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

void PrintVec2TXT(const std::vector<float>& v, std::string fn, bool visualization)
{
	std::ofstream fout;

	fout.open(fn);

	fout << std::setprecision(2) << std::fixed;

	int v_size = v.size();
	if (visualization == 0) {
		for (int i = 0; i < v.size(); i++) {
			fout << v[i] << std::endl;
		}
	} else {
		int sq_sz = (int)sqrt(v_size);
		int ind = 0;
		for (int i = 0; i < sq_sz; i++) {
			for (int j = 0; j < sq_sz; j++) {
				if (round(v[ind]) == 0) {
					fout << "     ";
				} else {
					fout << round(v[ind]) << " ";
				}
				ind += 1;
			}
		fout << std::endl;
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
void bzmesh2D(string path_in){
	std::cout << "******************************************************************************" << std::endl;

	string spline_cmd_tmd("../spline2D_src/spline -i " + path_in);
	const char* spline_cmd = spline_cmd_tmd.c_str();
	system(spline_cmd);
} 

// generating 3D bezier mesh using spline_src
void bzmesh3D(){
	std::cout << "******************************************************************************" << std::endl;
	string spline_cmd_tmd = "../NeuronTransportIGA/spline_src/spline ../io/3DNG/";
	const char* spline_cmd = spline_cmd_tmd.c_str();
	system(spline_cmd);
} 

// partitioning mesh using mpmetis
void mpmetis(int n_process, string path_in){
	string mpmetis_cmd_tmp("mpmetis " + path_in + "bzmeshinfo.txt " + to_string(n_process));
	const char* mpmetis_cmd = mpmetis_cmd_tmp.c_str();
	system(mpmetis_cmd);
	}

// partitioning mesh using mpmetis
void THS2D(string path_in, vector<int> rfid, vector<int> rftype){
	std::cout << "******************************************************************************" << std::endl;
	std::cout << "Local refinement based on Xiaodong's THS3D code ... " << std::endl;
	std::cout << "  - see: Truncated T-splines: Fundamentals and methods (2017)" << std::endl << std::endl;
	std::cout << "-----------------------------------------------------------------------------" << std::endl;
	std::cout << "Calling command | input mesh directory | refine ID | refine element type" << std::endl << std::endl;
	string ths2d_cmd_tmp("../THS2D/TTSP2D " + path_in + " ");
	for (int i = 0; i < rfid.size(); i++) {
		ths2d_cmd_tmp = ths2d_cmd_tmp + std::to_string(rfid[i]) + " ";
	}
	for (int i = 0; i < rftype.size(); i++) {
		ths2d_cmd_tmp = ths2d_cmd_tmp + std::to_string(rftype[i]) + " ";
	}
	std::cout << ths2d_cmd_tmp << std::endl;
	const char* ths2d_cmd = ths2d_cmd_tmp.c_str();
	system(ths2d_cmd);
}

//
void InitializeSoma(int numNeuron, vector<array<float, 2>> &seed, int &NX, int &NY){
    seed.resize(numNeuron);
    // 2D neuron soma initialization
    switch (numNeuron) {
        case 1:
            NX = 20;
            NY = 20;
	//     NX = 40;
        //     NY = 40;
            seed[0][0] = 20;    seed[0][1] = 20;
            break;
        case 2:
            NX = 100;
            NY = 45;
            seed[0][0] = 19;    seed[0][1] = 21;
            seed[1][0] = 79;   seed[1][1] = 22;
            break;
        case 3:
            NX = 140;
            NY = 130;
            seed[0][0] = 35;    seed[0][1] = 35;
            seed[1][0] = 105;   seed[1][1] = 35;
            seed[2][0] = 70;    seed[2][1] = 95;
            break;
        case 4:
            NX = 140;
            NY = 140;
            seed[0][0] = 35;    seed[0][1] = 35;
            seed[1][0] = 105;   seed[1][1] = 35;
            seed[2][0] = 35;    seed[2][1] = 105;
            seed[3][0] = 105;   seed[3][1] = 105;
            break;
        case 5:
            NX = 140;
            NY = 140;
            seed[0][0] = 35;    seed[0][1] = 35;
            seed[1][0] = 105;   seed[1][1] = 35;
            seed[2][0] = 35;    seed[2][1] = 105;
            seed[3][0] = 105;   seed[3][1] = 105;
            seed[4][0] = 70;    seed[4][1] = 70;
            break;
    }

}

vector<float> ConvertTo1DFloatVector(const vector<vector<int>> input) 
{
	vector<float> output;

	for (const auto& row : input) {
		for (int value : row) {
			output.push_back(static_cast<float>(value)); // Cast int value to float and add to the 1D vector
		}
	}

	return output;
}

vector<float> ConvertTo1DFloatVector(const vector<vector<float>> input) 
{
	vector<float> output;

	for (const auto& row : input) {
		for (float value : row) {
			output.push_back(static_cast<float>(value)); // Cast int value to float and add to the 1D vector
		}
	}

	return output;
	}


// Function to search for a particular x and y in the vector of Vertex2D
bool SearchPair(const vector<Vertex2D> prev_cpts, float targetX, float targetY, int &ind) {
	float x, y;
	for (int i = 0; i < prev_cpts.size(); i++) {

		// for (int i = 0; i < cpts_initial.size(); i++) {
		// 	if (abs(remainder(cpts_initial[i].coor[0],2)) != 1)
		// 		cpts_initial[i].coor[0] = round(cpts_initial[i].coor[0]);
		// 	if (abs(remainder(cpts_initial[i].coor[1],2)) != 1)
		// 		cpts_initial[i].coor[1] = round(cpts_initial[i].coor[1]);
		// }

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

		// if (abs(remainder(x,1)) == 0.5) {
		// 	x = round(x) - 1;
		// }
		// if (abs(remainder(y,1)) == 0.5) {
		// 	y = round(y) - 1;
		// }
		
		// if (prev_cpts[i].coor[0] == targetX && prev_cpts[i].coor[1] == targetY) {
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

	for (int i = 0; i < cpts.size(); i++) {
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
				} else {
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
				} else {
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
				} else {
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

	for (int i = 0; i < cpts_initial.size(); i++) {
		if (abs(remainder(cpts_initial[i].coor[0],1)) != 0.5)
			cpts_initial[i].coor[0] = round(cpts_initial[i].coor[0]);
		if (abs(remainder(cpts_initial[i].coor[1],1)) != 0.5)
			cpts_initial[i].coor[1] = round(cpts_initial[i].coor[1]);
	}
	
	for (int i = 0; i < cpts.size(); i++) {
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
	float wholePart = std::floor(value);
	float fractionalPart = value - wholePart;

	// Check if the fractional part is exactly 0.5
	if (fractionalPart == 0.5 || fractionalPart == -0.5) {
		return wholePart; // Round down to 0 for 0.5
	} else {
		return std::round(value); // Use standard rounding for other cases
	}
}

vector<float> InterpolateVars_coarse(vector<float> input, vector<Vertex2D> cpts_initial, vector<Vertex2D> cpts, int type) 
{	
	vector<float> output;
	output.resize(cpts.size());

	// float max_x(-1e5), max_y(-1e5), min_x(1e5), min_y(1e5);
	for (int i = 0; i < cpts_initial.size(); i++) {
		if (abs(remainder(cpts_initial[i].coor[0],2)) != 1)
			cpts_initial[i].coor[0] = round5(cpts_initial[i].coor[0]);
		if (abs(remainder(cpts_initial[i].coor[1],2)) != 1)
			cpts_initial[i].coor[1] = round5(cpts_initial[i].coor[1]);
	}
	
	for (int i = 0; i < cpts.size(); i++) {
		float x = round5(cpts[i].coor[0]);
		float y = round5(cpts[i].coor[1]);
		int ind;
		if (SearchPair(cpts_initial, round5(x), round5(y), ind)) {
			output[i] = input[ind];
		} 
		else {
			int indDown, indUp, indLeft, indRight;
			if ((abs(remainder(x,2)) == 1) && (abs(remainder(y,2)) != 1)) {
				if (SearchPair(cpts_initial, floorf(x)-1, round5(y), indDown) &&
					SearchPair(cpts_initial, floorf(x)+1, round5(y), indUp)) {
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
				if (SearchPair(cpts_initial, round5(x), floorf(y)-1, indLeft) &&
					SearchPair(cpts_initial, round5(x), floorf(y)+1, indRight)) {
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
				if (SearchPair(cpts_initial, floorf(x)-1, floorf(y)-1, indDown) &&
					SearchPair(cpts_initial, floorf(x)+1, floorf(y)+1, indUp) &&
					SearchPair(cpts_initial, floorf(x)-1, floorf(y)+1, indLeft) &&
					SearchPair(cpts_initial, floorf(x)+1, floorf(y)-1, indRight)) {
					if (type == 0) {
						output[i] = max(max(input[indDown], input[indUp]), max(input[indLeft], input[indRight]));
					} else if (type == 1) {
						output[i] = (input[indDown] + input[indUp] + input[indLeft] + input[indRight])/4;
					} else if (type == 2) {
						output[i] = 0;
					}
				} else {
					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck2)! x: %f y: %f | %f %f\n", x, y, floorf(x)-1, floorf(y)-1);
				}
			}
		}			
	}
	return output;
}

float distance_d(const Vertex2D& a, const Vertex2D& b) {
	return std::sqrt((a.coor[0] - b.coor[0]) * (a.coor[0] - b.coor[0]) + (a.coor[1] - b.coor[1]) * (a.coor[1] - b.coor[1]));
}

// Function to interpolate values for new control points
std::vector<float> interpolateValues_closest(
	const std::vector<float>& phi,
	const std::vector<Vertex2D>& cpt,
	const std::vector<Vertex2D>& cpt_out) {

	std::vector<float> interpolatedValues(cpt_out.size());

	for (size_t i = 0; i < cpt_out.size(); ++i) {
		float minDistance = std::numeric_limits<float>::max();
		size_t closestIndex = 0;

		// Find the closest point in cpt for each point in cpt_out
		for (size_t j = 0; j < cpt.size(); ++j) {
			float distance = distance_d(cpt_out[i], cpt[j]);
			if (distance < minDistance) {
				minDistance = distance;
				closestIndex = j;
			}
		}

		// Assume the value of the closest point
		// if (abs(round(phi[closestIndex])) > 0.5) {
		// 	interpolatedValues[i] = 1;
		// } else {
		// 	interpolatedValues[i] = 0;
		// }
		interpolatedValues[i] = phi[closestIndex];
	}

	return interpolatedValues;
}

// Function to interpolate values for new control points using the average of 6 closest points
std::vector<float> interpolateValues_averageN(
	const std::vector<float>& phi,
	const std::vector<Vertex2D>& cpt,
	const std::vector<Vertex2D>& cpt_out,
	float numAverage) {

	std::vector<float> interpolatedValues(cpt_out.size());

	for (size_t i = 0; i < cpt_out.size(); ++i) {
		// Use a max heap to keep track of the closest points
		// The heap will store pairs of (distance, index), where index is the index in `cpt`
		std::priority_queue<std::pair<float, size_t>> heap;

		for (size_t j = 0; j < cpt.size(); ++j) {
			float dist = distance_d(cpt_out[i], cpt[j]);
			// Keep the heap size to 6
			if (heap.size() < numAverage) {
				heap.push({dist, j});
			} else if (dist < heap.top().first) {
				heap.pop();
				heap.push({dist, j});
			}
		}

		// Calculate the average of the values of the 6 closest points
		float sum = 0.0;
		while (!heap.empty()) {
			sum += phi[heap.top().second];
			heap.pop();
		}
		interpolatedValues[i] = sum / numAverage;
	}

	return interpolatedValues;
}

// Function to interpolate values for new control points within a specified radius
std::vector<float> interpolateValuesWithinRadius(
	const std::vector<float>& phi,
	const std::vector<Vertex2D>& cpt,
	const std::vector<Vertex2D>& cpt_out,
	float radius) {

	std::vector<float> interpolatedValues(cpt_out.size());

	for (size_t i = 0; i < cpt_out.size(); ++i) {
		float sum = 0.0;
		size_t count = 0;

		for (size_t j = 0; j < cpt.size(); ++j) {
			float dist = distance_d(cpt_out[i], cpt[j]);
			if (dist <= radius) {
				sum += phi[j];
				++count;
			}
		}

		// Avoid division by zero
		interpolatedValues[i] = count > 0 ? sum / count : 0;
	}

	return interpolatedValues;
}

void ObtainRefineID_coarse(vector<float> phi, vector<Vertex2D> cpts, int NX, int NY, int originX, int originY, vector<int> &rfid, vector<int> &rftype){

	rfid.clear();
	rftype.clear();
	int ind = 0;
	int offset = 1 + (NY+1); // offset for using phi on cpts to determine local refinements on elements

	// vector<float> tmp;

	for (int i = offset; i < cpts.size() - offset; i++) {
		float x = cpts[i].coor[0];
		float y = cpts[i].coor[1];
		// std::cout << x << " ";
		rfid.push_back(ind - round(x/2) + originX + offset); // push back element id (calculated based on vertex id, ind)
		float averagePhi = (phi[i-(NY+1)+1] + phi[i-(NY+1)] + phi[i-(NY+1)-1] + phi[i-1] + phi[i] + phi[i+1] + phi[i+(NY+1)+1] + phi[i+(NY+1)] + phi[i+(NY+1)-1])/9;
		// tmp.push_back(averagePhi);
		
		if ((averagePhi < 0.95) && (averagePhi > 0.01)) {
		// if ((averagePhi < 0.95) && (averagePhi > 0.25)) {
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
			for (int i = (NY+1)+1; i < rftype.size()-(NY+1)-1; i++) {
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
			for (int i = 2*(NY+1)+2; i < rftype.size()-2*(NY+1)-2; i++) {
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

	for (int i = 0; i < rftype.size(); i++) {
		if ((rftype[i] == 5) && (rftype[i-1] != 5) && (rftype[i+1] != 5)) {
			rftype[i] = 1;
		}
	}

	for (int i = 0; i < rftype.size();) {
		if (rftype[i] == 5) {
			std::swap(rfid[i], rfid.back());
			rfid.pop_back();
			std::swap(rftype[i], rftype.back());
			rftype.pop_back();
		} else {
			i++;
		}
	}

	std::cout << std::endl;
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
		for (int i = 0; i < 4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		pts.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> pts[i].coor[0] >> pts[i].coor[1] >> tmp; // pts[i].coor[2];
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		mesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> mesh[i].IEN[0] >> mesh[i].IEN[1] >> mesh[i].IEN[2] >> mesh[i].IEN[3]; /* >>
				mesh[i].IEN[4] >> mesh[i].IEN[5] >> mesh[i].IEN[6] >> mesh[i].IEN[7]; */
			for (int j = 0; j < 4; j++)
			{
				mesh[i].pts[j][0] = pts[mesh[i].IEN[j]].coor[0];
				mesh[i].pts[j][1] = pts[mesh[i].IEN[j]].coor[1];
				// mesh[i].pts[j][2] = pts[mesh[i].IEN[j]].coor[2];
			}

		}
		for (int i = 0; i < neles + 5; i++) getline(fin, stmp);//skip lines
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

void WriteMatrixToPNG(const std::vector<float>& phi, int width, int height, const char* filename) {
	FILE *fp = fopen(filename, "wb");
	if (!fp) throw std::runtime_error("Failed to open file for writing");

	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
	if (!png_ptr) throw std::runtime_error("Failed to create PNG write structure");

	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
		png_destroy_write_struct(&png_ptr, nullptr);
		throw std::runtime_error("Failed to create PNG info structure");
	}

	if (setjmp(png_jmpbuf(png_ptr))) {
		png_destroy_write_struct(&png_ptr, &info_ptr);
		fclose(fp);
		throw std::runtime_error("PNG writing error");
	}

	png_init_io(png_ptr, fp);
	png_set_IHDR(png_ptr, info_ptr, width, height,
			8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	png_write_info(png_ptr, info_ptr);

	std::vector<png_bytep> row_pointers(height);
	for (int y = 0; y < height; ++y) {
		row_pointers[y] = (png_bytep) malloc(png_get_rowbytes(png_ptr, info_ptr));
		for (int x = 0; x < width; ++x) {
			float value = round(phi[y * width + x]);
			// Assuming the float values are normalized between 0 and 1
			row_pointers[y][x] = static_cast<png_byte>(value * 255);
		}
	}

	png_write_image(png_ptr, row_pointers.data());
	png_write_end(png_ptr, nullptr);

	for (int y = 0; y < height; y++) {
		free(row_pointers[y]);
	}

	png_destroy_write_struct(&png_ptr, &info_ptr);
	fclose(fp);
}