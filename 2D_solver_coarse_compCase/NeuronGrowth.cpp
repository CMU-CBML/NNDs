#include "NeuronGrowth.h"
#include "utils.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <stack>	// for tic toc
#include <ctime>	// for tic toc
#include <time.h>	// for srand
#include <queue>	// for geodesic distance
#include <cfloat>	// for FLT_MAX
#include <limits>	// for numeric_limits
#include "../nanoflann/1.5.5/include/nanoflann.hpp" // for KD tree
#include <map>		// for map in tip detection
#include <random>

using namespace std;
typedef unsigned int uint;
const float PI = 3.1415926;
const int INF = 1e9;

stack<clock_t> tictoc_stack;
float t_phi, t_synTub, t_syn, t_tub, t_tip, t_write(0), t_total(0);

void tic() 
{
	tictoc_stack.push(clock());
}
void toc(float& t) 
{
	t = ((float)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC;
	tictoc_stack.pop();
}

float MatrixDet(float dxdt[2][2])
{
	float det = dxdt[0][0] * dxdt[1][1] - dxdt[0][1] * dxdt[1][0];
	return det;
}

void Matrix2DInverse(float dxdt[2][2], float dtdx[2][2])
{
	float det = MatrixDet(dxdt);
	dtdx[0][0] = 1.0 / det * (dxdt[1][1]);
	dtdx[0][1] = 1.0 / det * (-dxdt[0][1]);
	dtdx[1][0] = 1.0 / det * (-dxdt[1][0]);
	dtdx[1][1] = 1.0 / det * (dxdt[0][0]);
}

NeuronGrowth::NeuronGrowth(){
	comm = PETSC_COMM_WORLD; //MPI_COMM_WORLD;
	mpiErr = MPI_Comm_rank(comm, &comRank);
	mpiErr = MPI_Comm_size(comm, &comSize);
	nProcess = comSize;
	judge_phi = 0;
	judge_syn = 0;
	judge_tub = 0;
	sum_grad_phi0_local = 0;
	sum_grad_phi0_global = 0;

	var_save_invl 	= 100;
	expandCK_invl	 = 100;
	aniso 		= 6;
	numNeuron 	= 1;
	gc_sz 		= 2;
	Diff 		= 4;
	M_axon 		= 100;
	M_neurite 	= 50;
	M_phi 		= 60;
	alpha 		= 0.9;
	alphaT 		= 0.001;
	betaT 		= 0.001;
	c_opt 		= 1;
	delta 		= 0.20;
	dt 		= 0.0005;
	epsilonb 	= 0.04;
	g 		= 0.1;
	gamma 		= 10;
	k2 		= 0;
	kappa 		= 2;
	Dc 		= 6;
	kp75 		= 0;
	r 		= 5;
	s_coeff 	= 0.007;
	seed_radius 	= 15;
	source_coeff 	= 15;
}

void NeuronGrowth::AssignProcessor(vector<vector<int>> &ele_proc)
{
	ele_process.clear();
	for (size_t i = 0; i < ele_proc[comRank].size(); i++)
		ele_process.push_back(ele_proc[comRank][i]);
}

void NeuronGrowth::SetVariables(string fn_par) {
	ifstream inputFile(fn_par);
	if (!inputFile.is_open()) {
		PetscPrintf(PETSC_COMM_WORLD, "Error: Unable to open file %s\n", fn_par.c_str());
		return;
	}

	string line;

	while (getline(inputFile, line)) {
		istringstream iss(line);
		string variableName;
		char equalsSign;
		if (!(iss >> variableName >> equalsSign)) {
			PetscPrintf(PETSC_COMM_WORLD, "Error: Reading variable failed!\n");
			break;
		}
		if (variableName == "expandCK_invl") iss >> expandCK_invl;
		else if (variableName == "var_save_invl") iss >> var_save_invl;
		else if (variableName == "numNeuron") iss >> numNeuron;
		else if (variableName == "gc_sz") iss >> gc_sz;
		else if (variableName == "aniso") iss >> aniso;
		else if (variableName == "gamma") iss >> gamma;
		else if (variableName == "seed_radius") iss >> seed_radius;
		else if (variableName == "kappa") iss >> kappa;
		else if (variableName == "dt") iss >> dt;
		else if (variableName == "Dc") iss >> Dc;
		else if (variableName == "kp75") iss >> kp75;
		else if (variableName == "k2") iss >> k2;
		else if (variableName == "c_opt") iss >> c_opt;
		else if (variableName == "alpha") iss >> alpha;
		else if (variableName == "M_phi") iss >> M_phi;
		else if (variableName == "M_axon") iss >> M_axon;
		else if (variableName == "M_neurite") iss >> M_neurite;
		else if (variableName == "s_coeff") iss >> s_coeff;
		else if (variableName == "delta") iss >> delta;
		else if (variableName == "epsilonb") iss >> epsilonb;
		else if (variableName == "r") iss >> r;
		else if (variableName == "g") iss >> g;
		else if (variableName == "alphaT") iss >> alphaT;
		else if (variableName == "betaT") iss >> betaT;
		else if (variableName == "Diff") iss >> Diff;
		else if (variableName == "source_coeff") iss >> source_coeff;
	}
	inputFile.close();
	PetscPrintf(PETSC_COMM_WORLD, "Parameter Loaded!\n");

	// Calculate alphaOverPi after all variables are loaded
	alphaOverPi = alpha / PI;
}

void NeuronGrowth::InitializeProblemNG(const int n_bz, vector<Vertex2D>& cpts, vector<Vertex2D> prev_cpts, const KDTree& kdTree_prev,
	vector<Vertex2D>& cpts_fine, vector<Vertex2D>& prev_cpts_fine, const KDTree& kdTree_prev_fine, vector<vector<float>> &NGvars, vector<array<float, 2>>& seed, const int& NX, const int& NY)
{
	MPI_Barrier(PETSC_COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD, "Initializing-----------------------------------------------------------------\n");

	/*Initialize parameters*/
	int cpt_sz = cpts.size();
	GaussInfo(3);
	n_bzmesh = n_bz;
	N_0.resize(cpt_sz);

	// float max_x, min_x, max_y, min_y;
	max_x = -INF; min_x = INF;
	max_y = -INF; min_y = INF;
	prev_max_x = -INF; prev_min_x = INF;
	prev_max_y = -INF; prev_min_y = INF;
	for (size_t i = 0; i<cpts.size(); i++) {
		if (n==0) {
			max_x = cpts[i].coor[0];
			min_x = cpts[i].coor[0];
			max_y = cpts[i].coor[1];
			min_y = cpts[i].coor[1];

			if (i < prev_cpts.size()) {
				prev_max_x = prev_cpts[i].coor[0];
				prev_max_x = prev_cpts[i].coor[0];
				prev_max_y = prev_cpts[i].coor[1];
				prev_min_y = prev_cpts[i].coor[1];
			}

		} else {
			max_x = max(cpts[i].coor[0],max_x);
			min_x = min(cpts[i].coor[0],min_x);
			max_y = max(cpts[i].coor[1],max_y);
			min_y = min(cpts[i].coor[1],min_y);

			if (i < prev_cpts.size()) {
				prev_max_x = max(prev_cpts[i].coor[0],prev_max_x);
				prev_min_x = min(prev_cpts[i].coor[0],prev_min_x);
				prev_max_y = max(prev_cpts[i].coor[1],prev_max_y);
				prev_min_y = min(prev_cpts[i].coor[1],prev_min_y);
			}
		}
	}

	PetscPrintf(PETSC_COMM_WORLD, "Checking maximum x and y value-----------------------------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD, "max x: %f min x: %f | max y: %f min y %f\n", max_x, min_x, max_y, min_y);
	PetscPrintf(PETSC_COMM_WORLD, "max x: %f min x: %f | max y: %f min y %f\n", prev_max_x, prev_min_x, prev_max_y, prev_min_y);

	float r, x, y; // radius, x, y
	float dx = 1;
	srand(time(NULL)); // random seed

	if (n == 0) { // Initialize
		for (size_t i = 0; i<cpts.size(); i++) {
			phi.push_back(0.f);
			tub.push_back(0.f);
			syn.push_back(0.f); // synaptogenesis
		}
		for (size_t i = 0; i<cpts_fine.size(); i++) {
			theta_fine.push_back((float)(rand()%100)/(float)100); // orientation theta
		}
		theta = InterpolateValues_closest(theta_fine, cpts_fine, cpts);

		// cout << prev_id.size() << " " << tmp.size() << endl;
		for (size_t i = 0; i<cpts.size(); i++) {
			x = cpts[i].coor[0];
			y = cpts[i].coor[1];

			// calculate new boundary label
			if (x==min_x || x==max_x || y==min_y || y==max_y) {
				cpts[i].label = 1;
			} else {
				cpts[i].label = 0;
			}   

			for (size_t j = 0; j<seed.size(); j++) { // multiple neurons
				r = sqrt((x-(float)seed[j][0])*(x-(float)seed[j][0])+(y-(float)seed[j][1])*(y-(float)seed[j][1])); // radius of initial soma
				if (r < (seed_radius*dx)) { // initial soma
					phi[i] = 1.0f;
					tub[i] = (0.5+0.5*tanh((sqrt(seed_radius*dx)-r)/2)); // based eqn in literature
				} 
			} 
		}

		phi_0 = phi; // initial phi
		tub_0 = tub; // initial tub

		vector<int> tmp((2*NX+1) * (2*NY+1), 0);
		prev_id = tmp;

	} else { // Collect from NGvars - continue simulation
		phi.clear(); 		phi.resize(cpts.size());
		syn.clear(); 		syn.resize(cpts.size());
		tub.clear(); 		tub.resize(cpts.size());
		theta.clear();		theta.resize(cpts.size());
		theta_fine.clear();	theta_fine.resize(cpts.size());
		phi_0.clear(); 		phi_0.resize(cpts.size());
		tub_0.clear(); 		tub_0.resize(cpts.size());
		prev_id.clear();

		phi = InterpolateVars_coarseKDtree(NGvars[0], prev_cpts, kdTree_prev, cpts, 1, 0);
		syn = InterpolateVars_coarseKDtree(NGvars[1], prev_cpts, kdTree_prev, cpts, 1, 0);
		tub = InterpolateVars_coarseKDtree(NGvars[2], prev_cpts, kdTree_prev, cpts, 1, 0);

		theta_fine = InterpolateVars_coarseKDtree(NGvars[3], prev_cpts_fine, kdTree_prev_fine, cpts_fine, 1, 1); // isTheta = 1 for special init
		theta = InterpolateValues_closest(theta_fine, cpts_fine, cpts);

		phi_0 = InterpolateVars_coarseKDtree(NGvars[4], prev_cpts, kdTree_prev, cpts, 1, 0);
		tub_0 = InterpolateVars_coarseKDtree(NGvars[5], prev_cpts, kdTree_prev, cpts, 1, 0);
		prev_id = ConvertTo1DIntVector(InterpolateValues_closest(NGvars[6], prev_cpts_fine, cpts_fine));
	}

	this->cpts = cpts; // for passing into formFunction

	// // writing initialized variables for debugging purposes
	// bool visualization = true;
	// PrintVec2TXT(phi, "/home/kuanrenqian/Research/3D_NeuronGrowth/io/2DNG/var_check/phi_inital.txt", visualization);
	// PrintVec2TXT(tub, "/home/kuanrenqian/Research/3D_NeuronGrowth/io/2DNG/var_check/tub_inital.txt", visualization);
	// PrintVec2TXT(theta, "/home/kuanrenqian/Research/3D_NeuronGrowth/io/2DNG/var_check/theta_inital.txt", visualization);
	// PrintVec2TXT(syn, "/home/kuanrenqian/Research/3D_NeuronGrowth/io/2DNG/var_check/syn_inital.txt", visualization);

	/*Initialize petsc vector, matrix*/
	PetscInt mat_dim = cpt_sz;

	// 250 = 2*5^3. (2 var, 5 basis, 3D) || 25 = 1*5^2. (1 var, 5 basis, 2D)	

	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_phi);
	
	ierr = MatCreate(PETSC_COMM_WORLD, &J);
	ierr = MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
	ierr = MatSetType(J, MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(J, 25, NULL, 25, NULL);
	ierr = MatSetOption(J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetUp(J);
	// MatGetOwnershipRange(J, NULL, NULL);

	ierr = MatCreate(PETSC_COMM_WORLD, &GK_phi);
	ierr = MatSetSizes(GK_phi, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
	ierr = MatSetType(GK_phi, MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(GK_phi, 25, NULL, 25, NULL);
	ierr = MatSetOption(GK_phi, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetUp(GK_phi);
	// MatGetOwnershipRange(GK_phi, NULL, NULL);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &GR_phi);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_phi);

	ierr = MatCreate(PETSC_COMM_WORLD, &GK_syn);
	ierr = MatSetSizes(GK_syn, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
	ierr = MatSetType(GK_syn, MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(GK_syn, 25, NULL, 25, NULL);
	ierr = MatSetOption(GK_syn, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetUp(GK_syn);
	// MatGetOwnershipRange(GK_syn, NULL, NULL);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &GR_syn);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_syn);

	ierr = MatCreate(PETSC_COMM_WORLD, &GK_tub);
	ierr = MatSetSizes(GK_tub, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
	ierr = MatSetType(GK_tub, MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(GK_tub, 25, NULL, 25, NULL);
	ierr = MatSetOption(GK_tub, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetUp(GK_tub);
	// MatGetOwnershipRange(GK_tub, NULL, NULL);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &GR_tub);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_tub);

	vars.reserve(3*sizeof(float));
	vars.resize(3);

	pre_EMatrixSolve.clear();
	pre_EVectorSolve.clear();
}

void NeuronGrowth::ToPETScVec(vector<float> input, Vec& petscVec)
{
	PetscInt localStart, localEnd;
	VecGetOwnershipRange(petscVec, &localStart, &localEnd);
	for (int i = localStart; i<localEnd; i++) {
		PetscScalar value = input[i]; // Adjust the index to match the localValues
		VecSetValues(petscVec, 1, &i, &value, INSERT_VALUES);
	}
	VecAssemblyBegin(petscVec);
	VecAssemblyEnd(petscVec);
}

void NeuronGrowth::ReadBezierElementProcess(string fn)
{
	string stmp;
	uint npts, neles, nfunctions, itmp;
	uint add(0);

	string fname_cmat = fn + "cmat.txt";
	ifstream fin_cmat;
	fin_cmat.clear();
	fin_cmat.open(fname_cmat);

	if (fin_cmat.is_open())
	{
		fin_cmat >> neles;
		bzmesh_process.clear();
		bzmesh_process.resize(ele_process.size());
		for (size_t i = 0; i < neles; i++)
		{
			if ((i == ele_process[add]) && (add < ele_process.size()))
			{
				fin_cmat >> itmp >> nfunctions >> bzmesh_process[add].type;
				bzmesh_process[add].cmat.resize(nfunctions);
				bzmesh_process[add].IEN.resize(nfunctions);
				for (size_t j = 0; j < nfunctions; j++)
					fin_cmat >> bzmesh_process[add].IEN[j];
				for (size_t j = 0; j < nfunctions; j++)
				{
					for (int k = 0; k < 16; k++)
					{
						fin_cmat >> bzmesh_process[add].cmat[j][k];
					}
				}
				add++;
			}
			else
			{
				fin_cmat >> stmp >> nfunctions >> itmp;
				for (size_t j = 0; j < nfunctions; j++)
					fin_cmat >> stmp;
				for (size_t j = 0; j < nfunctions; j++)
					for (size_t k = 0; k < 16; k++)
						fin_cmat >> stmp;
			}
		}
		fin_cmat.close();
		// PetscPrintf(PETSC_COMM_WORLD, "Bezier Matrices Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_cmat.c_str());
	}

	string fname_bzpt = fn + "bzpt.txt";
	ifstream fin_bzpt;
	fin_bzpt.open(fname_bzpt);
	add = 0;
	float tmp_xuan;
	if (fin_bzpt.is_open())
	{
		fin_bzpt >> npts;
		getline(fin_bzpt, stmp);
		for (size_t e = 0; e < neles; e++)
		{
			if ((e == ele_process[add]) && (add < ele_process.size()))
			{
				bzmesh_process[add].pts.resize(bzpt_num);
				for (size_t i = 0; i < bzpt_num; i++)
				{
					fin_bzpt >> bzmesh_process[add].pts[i][0] >> bzmesh_process[add].pts[i][1] >> tmp_xuan;
				}
				add++;
			}
			else
			{
				for (size_t i = 0; i < bzpt_num; i++)
					fin_bzpt >> stmp >> stmp >> stmp;
			}

		}
		fin_bzpt.close();
		// PetscPrintf(PETSC_COMM_WORLD, "Bezier Points Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_bzpt.c_str());
	}
}

void NeuronGrowth::GaussInfo(int ng)
{
	Gpt.clear();
	wght.clear();
	switch (ng)
	{
		case 2:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0] = 0.2113248654051871;
			Gpt[1] = 0.7886751345948129;
			wght[0] = 1.;
			wght[1] = 1.;
			break;
		}
		case 3:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0] = 0.1127016653792583;
			Gpt[1] = 0.5;
			Gpt[2] = 0.8872983346207417;
			wght[0] = 0.5555555555555556;
			wght[1] = 0.8888888888888889;
			wght[2] = 0.5555555555555556;
			break;
		}
		case 4:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0] = 0.06943184420297371;
			Gpt[1] = 0.33000947820757187;
			Gpt[2] = 0.6699905217924281;
			Gpt[3] = 0.9305681557970262;
			wght[0] = 0.3478548451374539;
			wght[1] = 0.6521451548625461;
			wght[2] = 0.6521451548625461;
			wght[3] = 0.3478548451374539;
			break;
		}
		case 5:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0] = 0.046910077030668;
			Gpt[1] = 0.2307653449471585;
			Gpt[2] = 0.5;
			Gpt[3] = 0.7692346550528415;
			Gpt[4] = 0.953089922969332;
			wght[0] = 0.2369268850561891;
			wght[1] = 0.4786286704993665;
			wght[2] = 0.5688888888888889;
			wght[3] = 0.4786286704993665;
			wght[4] = 0.2369268850561891;
			break;
		}
		default:
		{
			Gpt.resize(2);
			wght.resize(2);
			Gpt[0] = 0.2113248654051871;
			Gpt[1] = 0.7886751345948129;
			wght[0] = 1.;
			wght[1] = 1.;
			break;
		}
	}
}

void NeuronGrowth::BasisFunction(float u, float v, const vector<array<float, 2>>& pt, const vector<array<float, 16>> &cmat,
	vector<float> &Nx, vector<array<float, 2>>& dNdx, float dudx[2][2], float& detJ)
{
	float Nu[4] = {(1 - u) * (1 - u) * (1 - u), 3 * (1 - u) * (1 - u) * u, 3 * (1 - u) * u * u, u * u * u};
	float Nv[4] = {(1 - v) * (1 - v) * (1 - v), 3 * (1 - v) * (1 - v) * v, 3 * (1 - v) * v * v, v * v * v};
	// float Nw[4] = { (1 - w)*(1 - w)*(1 - w), 3.*(1 - w)*(1 - w)*w, 3.*(1 - w)*w*w, w*w*w };
	float dNdu[4] = {-3 * (1 - u) * (1 - u), 3 - 12 * u + 9 * u * u, 3 * (2 - 3 * u) * u, 3 * u * u};
	float dNdv[4] = {-3 * (1 - v) * (1 - v), 3 - 12 * v + 9 * v * v, 3 * (2 - 3 * v) * v, 3 * v * v};
	// float dNdw[4] = { -3.*(1 - w)*(1 - w), 3 - 12.*w + 9.*w*w, 3.*(2 - 3.*w)*w, 3.*w*w };
	float dNdt[bzpt_num][2];
	float Nx_bz[bzpt_num];
	float dNdx_bz[bzpt_num][2];

	Nx.clear();
	dNdx.clear();
	Nx.resize(cmat.size(), 0);
	dNdx.resize(cmat.size(), {0, 0});

	uint i, j, a, b, loc;
	loc = 0;
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			// for (k = 0; k < 4; k++)	{
			Nx_bz[loc] = Nu[j] * Nv[i];	// * Nw[i];
			dNdt[loc][0] = dNdu[j] * Nv[i]; // * Nw[i];
			dNdt[loc][1] = Nu[j] * dNdv[i]; // * Nw[i];
			// dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
			loc++;
			// }
		}
	}

	float dxdt[2][2] = {{0}};
	for (loc = 0; loc < bzpt_num; loc++)
		for (a = 0; a < 2; a++)
			for (b = 0; b < 2; b++)
				dxdt[a][b] += pt[loc][a] * dNdt[loc][b];

	float dtdx[2][2] = {{0}};
	Matrix2DInverse(dxdt, dtdx);

	// 1st derivatives
	for (i = 0; i < bzpt_num; i++)
	{
		dNdx_bz[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0]; // + dNdt[i][2] * dtdx[2][0];
		dNdx_bz[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1]; // + dNdt[i][2] * dtdx[2][1];
										   // dNdx_bz[i][2] = dNdt[i][0] * dtdx[0][2] + dNdt[i][1] * dtdx[1][2] + dNdt[i][2] * dtdx[2][2];
	}
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			dudx[i][j] = dtdx[i][j];

	detJ = MatrixDet(dxdt);
	detJ = 0.25 * detJ;

	for (i = 0; i < cmat.size(); i++)
	{
		for (j = 0; j < bzpt_num; j++)
		{
			Nx[i] += cmat[i][j] * Nx_bz[j];
			for (int m = 0; m < 2; m++)
			{
				dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
			}
		}
	}
}

void NeuronGrowth::BasisFunction(float u, float v, int nen, const vector<array<float, 2>>& pt, const vector<array<float, 16>> &cmat,
	vector<float> &Nx, vector<array<float, 2>>& dNdx, vector<array<array<float, 2>, 2>> &dN2dx2, float dudx[2][2], float& detJ)
{
	float Nu[4] = {(1 - u) * (1 - u) * (1 - u), 3 * (1 - u) * (1 - u) * u, 3 * (1 - u) * u * u, u * u * u};
	float Nv[4] = {(1 - v) * (1 - v) * (1 - v), 3 * (1 - v) * (1 - v) * v, 3 * (1 - v) * v * v, v * v * v};
	// float Nw[4] = { (1 - w)*(1 - w)*(1 - w), 3.*(1 - w)*(1 - w)*w, 3.*(1 - w)*w*w, w*w*w };

	float dNdu[4] = {-3 * (1 - u) * (1 - u), 3 - 12 * u + 9 * u * u, 3 * (2 - 3 * u) * u, 3 * u * u};
	float dNdv[4] = {-3 * (1 - v) * (1 - v), 3 - 12 * v + 9 * v * v, 3 * (2 - 3 * v) * v, 3 * v * v};
	// float dNdw[4] = { -3.*(1 - w)*(1 - w), 3 - 12.*w + 9.*w*w, 3.*(2 - 3.*w)*w, 3.*w*w };

	float dN2du2[4] = {6 * (1 - u), -12 + 18 * u, 6 - 18 * u, 6 * u};
	float dN2dv2[4] = {6 * (1 - v), -12 + 18 * v, 6 - 18 * v, 6 * v};
	// float dN2dw2[4] = { 6.*(1 - w),-12 + 18.*w,6 - 18.*w,6.*w };

	float dNdt[bzpt_num][2];
	float dN2dt2[bzpt_num][2][2];
	float Nx_bz[bzpt_num];
	float dNdx_bz[bzpt_num][2];
	float dN2dx2_bz[bzpt_num][2][2];

	Nx.clear();
	dNdx.clear();
	dN2dx2.clear();
	Nx.resize(nen, 0);
	dNdx.resize(nen, {0});
	dN2dx2.resize(nen, {{0}});

	int i, j, k, a, b, c, loc;
	loc = 0;
	for (j = 0; j < 4; j++)
	{
		for (k = 0; k < 4; k++)
		{
			Nx_bz[loc] = Nu[k] * Nv[j];
			dNdt[loc][0] = dNdu[k] * Nv[j];
			dNdt[loc][1] = Nu[k] * dNdv[j];
			dN2dt2[loc][0][0] = dN2du2[k] * Nv[j];
			dN2dt2[loc][0][1] = dNdu[k] * dNdv[j];
			dN2dt2[loc][1][0] = dNdu[k] * dNdv[j];
			dN2dt2[loc][1][1] = Nu[k] * dN2dv2[j];
			loc++;
		}
	}

	float dxdt[2][2] = {{0}};
	for (loc = 0; loc < bzpt_num; loc++)
	{
		for (a = 0; a < 2; a++)
		{
			for (b = 0; b < 2; b++)
			{
				dxdt[a][b] += pt[loc][a] * dNdt[loc][b];
			}
		}
	}

	float dtdx[2][2] = {{0}};
	Matrix2DInverse(dxdt, dtdx);

	// 1st derivatives
	for (i = 0; i < bzpt_num; i++)
	{
		dNdx_bz[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0];
		dNdx_bz[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1];
	}
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			dudx[i][j] = dtdx[i][j];

	detJ = MatrixDet(dxdt);
	detJ = 0.25 * detJ;

	// 2nd derivatives
	float dx2dt2[2][4] = {{0}};
	float dt2dx2[2][4] = {{0}};
	for (int l = 0; l < 2; l++)
		for (loc = 0; loc < bzpt_num; loc++)
			for (a = 0; a < 2; a++)
				for (b = 0; b < 2; b++)
					dx2dt2[l][2 * a + b] += pt[loc][l] * dN2dt2[loc][a][b];

	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			for (k = 0; k < 2; k++)
				for (a = 0; a < 2; a++)
					for (b = 0; b < 2; b++)
						for (c = 0; c < 2; c++)
							dt2dx2[c][2 * i + j] -= dx2dt2[k][2 * a + b] * dtdx[a][i] * dtdx[b][j] * dtdx[c][k];

	for (loc = 0; loc < bzpt_num; loc++)
		for (i = 0; i < 2; i++)
			for (j = 0; j < 2; j++)
				dN2dx2_bz[loc][i][j] = 0.;

	for (loc = 0; loc < bzpt_num; loc++)
	{
		for (i = 0; i < 2; i++)
		{
			for (j = 0; j < 2; j++)
			{
				for (a = 0; a < 2; a++)
				{
					for (b = 0; b < 2; b++)
					{
						dN2dx2_bz[loc][i][j] += dN2dt2[loc][a][b] * dtdx[a][i] * dtdx[b][j];
					}
					dN2dx2_bz[loc][i][j] += dNdt[loc][a] * dt2dx2[a][2 * i + j];
				}
			}
		}
	}

	for (i = 0; i < nen; i++)
	{
		for (j = 0; j < bzpt_num; j++)
		{
			Nx[i] += cmat[i][j] * Nx_bz[j];
			for (int m = 0; m < 2; m++)
			{
				dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
				for (int n = 0; n < 2; n++)
				{
					dN2dx2[i][m][n] += cmat[i][j] * dN2dx2_bz[j][m][n];
				}
			}
		}
	}
}

void NeuronGrowth::ApplyBoundaryCondition(const float bc_value, int pt_num, int variable_num, vector<vector<float>> &EMatrixSolve, vector<float> &EVectorSolve)
{
	int j;
	int nen = EVectorSolve.size();
	for (j = 0; j < nen; j++)
	{
		EVectorSolve[j] -= bc_value * EMatrixSolve[j][pt_num + variable_num * nen];
	}
	for (j = 0; j < nen; j++)
	{
		EMatrixSolve[j][pt_num + variable_num * nen] = 0.0;
		EMatrixSolve[pt_num + variable_num * nen][j] = 0.0;
	}
	EMatrixSolve[pt_num + variable_num * nen][pt_num + variable_num * nen] = 1.0;
	EVectorSolve[pt_num + variable_num * nen] = bc_value;
}

void NeuronGrowth::MatrixAssembly(vector<vector<float>> &EMatrixSolve, const vector<int> &IEN, Mat &GK)
{
	// int A;
	// int row_start, row_end, row_now;
	uint add = 0;
	uint nen = IEN.size();

	PetscInt *nodeList = new PetscInt[nen * 1];
	PetscReal *tmpGK = new PetscReal[nen * 1 * nen * 1];

	for (size_t m = 0; m < IEN.size(); m++)
	{
		int A = IEN[m];
		nodeList[1 * m] = 1 * A;
		for (size_t n = 0; n < IEN.size(); n++)
		{
				tmpGK[add] = EMatrixSolve[m][n];
				add++;
		}
	}
	MatSetValues(GK, nen * 1, nodeList, nen * 1, nodeList, tmpGK, ADD_VALUES);
	delete nodeList;
	delete tmpGK;
}

void NeuronGrowth::ResidualAssembly(vector<float> &EVectorSolve, const vector<int> &IEN, Vec &GR)
{
	int A;
	int add = 0;
	uint nen = IEN.size();

	PetscInt *nodeList = new PetscInt[nen];
	PetscReal *tmpGR = new PetscReal[nen];

	for (size_t i = 0; i < IEN.size(); i++)
	{
		A = IEN[i];
		nodeList[i * 1 + 0] = A * 1 + 0;
		tmpGR[add] = EVectorSolve[i];
		add += 1;
	}
	VecSetValues(GR, nen * 1, nodeList, tmpGR, ADD_VALUES);
	delete nodeList;
	delete tmpGR;
}

void NeuronGrowth::VisualizeVTK_ControlMesh(const vector<Vertex2D>& spt, const vector<Element2D>& mesh, int step, string fn, vector<float> var, string varName)
{
	string fname;
	stringstream ss;
	ss << setw(6) << setfill('0') << step;
	ofstream fout;
	unsigned int i;
	fname = fn + "/controlmesh_" + varName + ss.str() + ".vtk";
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i < spt.size(); i++)
		{
			fout << spt[i].coor[0] << " " << spt[i].coor[1] << " " << 0.00000 /* spt[i].coor[2] */ << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
		for (i = 0; i < mesh.size(); i++)
		{
			fout << "4 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3]
			     << " " /* << mesh[i].IEN[4] << " " << mesh[i].IEN[5] << " " << mesh[i].IEN[6] << " " << mesh[i].IEN[7] */ << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << '\n';
		for (i = 0; i < mesh.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << var.size() << "\nSCALARS AllParticles float 1\nLOOKUP_TABLE default\n";
		for (size_t i = 0; i < var.size(); i++)
		{
			fout << var[i] /* +N_plus[i]+N_minus[i] */ << "\n";
			// fout << spt[i].label /* +N_plus[i]+N_minus[i] */ << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void NeuronGrowth::ConcentrationCal_Coupling_Bezier(float u, float v, const Element2D &bzel, float pt[2], float& disp, float dudx[2], float& detJ)
{
	float dUdx[2][2];
	vector<float> Nx(bzel.IEN.size());
	vector<array<float, 2>> dNdx(bzel.IEN.size());
	bzel.Para2Phys(u, v, pt);
	BasisFunction(u, v, bzel.pts, bzel.cmat, Nx, dNdx, dUdx, detJ);
	disp = 0.;
	dudx[0] = 0.;
	dudx[1] = 0.; // dudx[2] = 0.;
	for (size_t i = 0; i < bzel.IEN.size(); i++)
	{
		disp += Nx[i] * (N_0[bzel.IEN[i]]);	    // + N_plus[bzel.IEN[i]]);
		dudx[0] += dNdx[i][0] * (N_0[bzel.IEN[i]]); // + N_plus[bzel.IEN[i]]);
		dudx[1] += dNdx[i][1] * (N_0[bzel.IEN[i]]); // + N_plus[bzel.IEN[i]]);
							    // dudx[2] += dNdx[i][2] * (N_0[bzel.IEN[i]] + N_plus[bzel.IEN[i]]);
	}
}

void NeuronGrowth::VisualizeVTK_PhysicalDomain(int step, string var, string fn)
{
	vector<array<float, 3>> spt_all; // sample points
	vector<float> sresult_all;
	vector<array<int, 4>> sele_all;
	float detJ;
	int num_bzmesh_ele = bzmesh_process.size();
	float spt_proc[num_bzmesh_ele * 4 * 3];
	float sresult_proc[num_bzmesh_ele * 4 * 1];
	int sele_proc[num_bzmesh_ele * 4];

	for (int e = 0; e < num_bzmesh_ele; e++)
	{
		int ns(2);
		vector<float> su(ns);
		for (int i = 0; i < ns; i++)
		{
			su[i] = float(i) / (float(ns) - 1.);
		}

		int loc(0);
		for (int a = 0; a < ns; a++)
		{
			for (int b = 0; b < ns; b++)
			{
				// for (int c = 0; c < ns; c++)
				// {
				float pt1[2], dudx[2];
				float result;
				ConcentrationCal_Coupling_Bezier(su[b], su[a], bzmesh_process[e], pt1, result, dudx, detJ);
				spt_proc[12 * e + loc * 3 + 0] = pt1[0];
				spt_proc[12 * e + loc * 3 + 1] = pt1[1];
				spt_proc[12 * e + loc * 3 + 2] = 0.000000; // pt1[2];
				sresult_proc[4 * e + loc] = result;
				loc++;
				// }
			}
		}
		// int nns[2] = {ns * ns, ns};
		for (int a = 0; a < ns - 1; a++)
		{
			for (int b = 0; b < ns - 1; b++)
			{
				// for (int c = 0; c < ns - 1; c++)
				// {
				sele_proc[4 * e + 0] = 4 * e + a * ns + b;
				sele_proc[4 * e + 1] = 4 * e + a * ns + b + 1;
				sele_proc[4 * e + 2] = 4 * e + (a + 1) * ns + (b + 1);
				sele_proc[4 * e + 3] = 4 * e + (a + 1) * ns + b;
				// sele_proc[8 * e + 4] = 8 * e + (a + 1)*nns[1] + b*ns + c;
				// sele_proc[8 * e + 5] = 8 * e + (a + 1)*nns[1] + b*ns + c + 1;
				// sele_proc[8 * e + 6] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
				// sele_proc[8 * e + 7] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c;
				// }
			}
		}
	}

	float *spts = NULL;
	float *sresults = NULL;
	int *seles = NULL;
	int *displs_spts = NULL;
	int *displs_sresults = NULL;
	int *displs_seles = NULL;
	int *num_bzmesh_eles = NULL;
	int *recvcounts_spts = NULL;
	int *recvcounts_sresults = NULL;
	int *recvcounts_seles = NULL;

	if (comRank == 0)
	{
		num_bzmesh_eles = (int *)malloc(sizeof(int) * nProcess);
		recvcounts_spts = (int *)malloc(sizeof(int) * nProcess);
		recvcounts_sresults = (int *)malloc(sizeof(int) * nProcess);
		recvcounts_seles = (int *)malloc(sizeof(int) * nProcess);
	}
	MPI_Gather(&num_bzmesh_ele, 1, MPI_INT, num_bzmesh_eles, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Barrier(comm);

	if (comRank == 0)
	{
		spts = (float *)malloc(sizeof(float) * 12 * n_bzmesh);
		sresults = (float *)malloc(sizeof(float) * 4 * n_bzmesh);
		seles = (int *)malloc(sizeof(int) * 4 * n_bzmesh);

		displs_spts = (int *)malloc(nProcess * sizeof(int));
		displs_sresults = (int *)malloc(nProcess * sizeof(int));
		displs_seles = (int *)malloc(nProcess * sizeof(int));
		displs_spts[0] = 0;
		displs_sresults[0] = 0;
		displs_seles[0] = 0;

		for (int i = 1; i < nProcess; i++)
		{
			displs_spts[i] = displs_spts[i - 1] + num_bzmesh_eles[i - 1] * 12;
			displs_sresults[i] = displs_sresults[i - 1] + num_bzmesh_eles[i - 1] * 4;
			displs_seles[i] = displs_seles[i - 1] + num_bzmesh_eles[i - 1] * 4;
		}

		for (int i = 0; i < nProcess; i++)
		{
			recvcounts_spts[i] = num_bzmesh_eles[i] * 12;
			recvcounts_sresults[i] = num_bzmesh_eles[i] * 4;
			recvcounts_seles[i] = num_bzmesh_eles[i] * 4;
		}
	}

	MPI_Gatherv(spt_proc, num_bzmesh_ele * 4 * 3, MPI_FLOAT, spts, recvcounts_spts, displs_spts, MPI_FLOAT, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sresult_proc, num_bzmesh_ele * 4, MPI_FLOAT, sresults, recvcounts_sresults, displs_sresults, MPI_FLOAT, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sele_proc, num_bzmesh_ele * 4, MPI_INT, seles, recvcounts_seles, displs_seles, MPI_INT, 0, PETSC_COMM_WORLD);

	if (comRank == 0)
	{
		for (int i = 0; i < n_bzmesh; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				array<float, 3> pt = {spts[i * 12 + j * 3 + 0], spts[i * 12 + j * 3 + 1], spts[i * 12 + j * 3 + 2]};
				spt_all.push_back(pt);
				sresult_all.push_back(sresults[i * 4 + j]);
			}
		}
		int sum_ele = 0;
		int pstart = 0;
		for (int i = 0; i < nProcess; i++)
		{
			for (int e = 0; e < num_bzmesh_eles[i]; e++)
			{
				array<int, 4> el;
				el[0] = pstart + seles[4 * sum_ele + 0];
				el[1] = pstart + seles[4 * sum_ele + 1];
				el[2] = pstart + seles[4 * sum_ele + 2];
				el[3] = pstart + seles[4 * sum_ele + 3];
				// el[4] = pstart + seles[8 * sum_ele + 4];
				// el[5] = pstart + seles[8 * sum_ele + 5];
				// el[6] = pstart + seles[8 * sum_ele + 6];
				// el[7] = pstart + seles[8 * sum_ele + 7];
				sele_all.push_back(el);
				sum_ele++;
			}
			pstart = pstart + num_bzmesh_eles[i] * 4;
		}
		// cout << "Visualizing in Physical Domain...\n";
		WriteVTK(spt_all, sresult_all, sele_all, step, var, fn);
	}
}

void NeuronGrowth::WriteVTK(const vector<array<float, 3>> spt, const vector<float> sdisp, const vector<array<int, 4>> sele, int step, string var, string fn)
{
	stringstream ss;
	// ss << step;
	ss << setw(6) << setfill('0') << step;
	string fname = fn + "/" + var + "physics_singleparticle_" + ss.str() + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i < spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
			     << " " /* << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] */ << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << sdisp.size() << "\nSCALARS AllParticle float 1\nLOOKUP_TABLE default\n";
		for (size_t i = 0; i < sdisp.size(); i++)
		{
			fout << sdisp[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void NeuronGrowth::CalculateVarsForOutput(vector<array<float, 3>> &spt_all, vector<float> &sresult_all, vector<array<int, 4>> &sele_all) 
{
	float detJ;
	int num_bzmesh_ele = bzmesh_process.size();
	float spt_proc[num_bzmesh_ele * 4 * 3];
	float sresult_proc[num_bzmesh_ele * 4 * 1];
	int sele_proc[num_bzmesh_ele * 4];

	for (int e = 0; e < num_bzmesh_ele; e++)
	{
		int ns(2);
		vector<float> su(ns);
		for (int i = 0; i < ns; i++)
		{
			su[i] = float(i) / (float(ns) - 1.);
		}

		int loc(0);
		for (int a = 0; a < ns; a++)
		{
			for (int b = 0; b < ns; b++)
			{
				// for (int c = 0; c < ns; c++)
				// {
				float pt1[2], dudx[2];
				float result;
				ConcentrationCal_Coupling_Bezier(su[b], su[a], bzmesh_process[e], pt1, result, dudx, detJ);
				spt_proc[12 * e + loc * 3 + 0] = pt1[0];
				spt_proc[12 * e + loc * 3 + 1] = pt1[1];
				spt_proc[12 * e + loc * 3 + 2] = 0.000000; // pt1[2];
				sresult_proc[4 * e + loc] = result;
				loc++;
				// }
			}
		}
		// int nns[2] = {ns * ns, ns};
		for (int a = 0; a < ns - 1; a++)
		{
			for (int b = 0; b < ns - 1; b++)
			{
				// for (int c = 0; c < ns - 1; c++)
				// {
				sele_proc[4 * e + 0] = 4 * e + a * ns + b;
				sele_proc[4 * e + 1] = 4 * e + a * ns + b + 1;
				sele_proc[4 * e + 2] = 4 * e + (a + 1) * ns + (b + 1);
				sele_proc[4 * e + 3] = 4 * e + (a + 1) * ns + b;
				// sele_proc[8 * e + 4] = 8 * e + (a + 1)*nns[1] + b*ns + c;
				// sele_proc[8 * e + 5] = 8 * e + (a + 1)*nns[1] + b*ns + c + 1;
				// sele_proc[8 * e + 6] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
				// sele_proc[8 * e + 7] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c;
				// }
			}
		}
	}

	float *spts = NULL;
	float *sresults = NULL;
	int *seles = NULL;
	int *displs_spts = NULL;
	int *displs_sresults = NULL;
	int *displs_seles = NULL;
	int *num_bzmesh_eles = NULL;
	int *recvcounts_spts = NULL;
	int *recvcounts_sresults = NULL;
	int *recvcounts_seles = NULL;

	if (comRank == 0)
	{
		num_bzmesh_eles = (int *)malloc(sizeof(int) * nProcess);
		recvcounts_spts = (int *)malloc(sizeof(int) * nProcess);
		recvcounts_sresults = (int *)malloc(sizeof(int) * nProcess);
		recvcounts_seles = (int *)malloc(sizeof(int) * nProcess);
	}
	MPI_Gather(&num_bzmesh_ele, 1, MPI_INT, num_bzmesh_eles, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Barrier(comm);

	if (comRank == 0)
	{
		spts = (float *)malloc(sizeof(float) * 12 * n_bzmesh);
		sresults = (float *)malloc(sizeof(float) * 4 * n_bzmesh);
		seles = (int *)malloc(sizeof(int) * 4 * n_bzmesh);

		displs_spts = (int *)malloc(nProcess * sizeof(int));
		displs_sresults = (int *)malloc(nProcess * sizeof(int));
		displs_seles = (int *)malloc(nProcess * sizeof(int));
		displs_spts[0] = 0;
		displs_sresults[0] = 0;
		displs_seles[0] = 0;

		for (int i = 1; i < nProcess; i++)
		{
			displs_spts[i] = displs_spts[i - 1] + num_bzmesh_eles[i - 1] * 12;
			displs_sresults[i] = displs_sresults[i - 1] + num_bzmesh_eles[i - 1] * 4;
			displs_seles[i] = displs_seles[i - 1] + num_bzmesh_eles[i - 1] * 4;
		}

		for (int i = 0; i < nProcess; i++)
		{
			recvcounts_spts[i] = num_bzmesh_eles[i] * 12;
			recvcounts_sresults[i] = num_bzmesh_eles[i] * 4;
			recvcounts_seles[i] = num_bzmesh_eles[i] * 4;
		}
	}

	MPI_Gatherv(spt_proc, num_bzmesh_ele * 4 * 3, MPI_FLOAT, spts, recvcounts_spts, displs_spts, MPI_FLOAT, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sresult_proc, num_bzmesh_ele * 4, MPI_FLOAT, sresults, recvcounts_sresults, displs_sresults, MPI_FLOAT, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sele_proc, num_bzmesh_ele * 4, MPI_INT, seles, recvcounts_seles, displs_seles, MPI_INT, 0, PETSC_COMM_WORLD);

	if (comRank == 0)
	{
		for (int i = 0; i < n_bzmesh; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				array<float, 3> pt = {spts[i * 12 + j * 3 + 0], spts[i * 12 + j * 3 + 1], spts[i * 12 + j * 3 + 2]};
				spt_all.push_back(pt);
				sresult_all.push_back(sresults[i * 4 + j]);
			}
		}
		int sum_ele = 0;
		int pstart = 0;
		for (int i = 0; i < nProcess; i++)
		{
			for (int e = 0; e < num_bzmesh_eles[i]; e++)
			{
				array<int, 4> el;
				el[0] = pstart + seles[4 * sum_ele + 0];
				el[1] = pstart + seles[4 * sum_ele + 1];
				el[2] = pstart + seles[4 * sum_ele + 2];
				el[3] = pstart + seles[4 * sum_ele + 3];
				// el[4] = pstart + seles[8 * sum_ele + 4];
				// el[5] = pstart + seles[8 * sum_ele + 5];
				// el[6] = pstart + seles[8 * sum_ele + 6];
				// el[7] = pstart + seles[8 * sum_ele + 7];
				sele_all.push_back(el);
				sum_ele++;
			}
			pstart = pstart + num_bzmesh_eles[i] * 4;
		}
	}
}
void NeuronGrowth::VisualizeVTK_PhysicalDomain_All(int step, string fn)
{
	vector<vector<array<float, 3>>> spt_all_5var; // sample points
	spt_all_5var.resize(5);
	vector<vector<float>> sresult_all_5var;
	sresult_all_5var.resize(5);
	vector<vector<array<int, 4>>> sele_all_5var;
	sele_all_5var.resize(5);
	N_0 = phi;
	CalculateVarsForOutput(spt_all_5var[0], sresult_all_5var[0], sele_all_5var[0]);
	N_0 = syn;
	CalculateVarsForOutput(spt_all_5var[1], sresult_all_5var[1], sele_all_5var[1]);
	N_0 = tub;
	CalculateVarsForOutput(spt_all_5var[2], sresult_all_5var[2], sele_all_5var[2]);
	N_0 = tips;
	CalculateVarsForOutput(spt_all_5var[3], sresult_all_5var[3], sele_all_5var[3]);
	N_0 = theta;
	CalculateVarsForOutput(spt_all_5var[4], sresult_all_5var[4], sele_all_5var[4]);
	if (comRank == 0)
	{
		WriteVTK_All(spt_all_5var[0], sresult_all_5var, sele_all_5var[0], step, fn);
	}
}

void NeuronGrowth::WriteVTK_All(const vector<array<float, 3>> spt, const vector<vector<float>> sdisp, const vector<array<int, 4>> sele, int step, string fn)
{
	stringstream ss;
	// ss << step;
	ss << setw(6) << setfill('0') << step;
	string fname = fn + "/physics_allparticle_" + ss.str() + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i < spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
			     << " " /* << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] */ << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << sdisp[0].size() << "\nSCALARS " << "phi " << "float 1\nLOOKUP_TABLE default\n";
		for (size_t i = 0; i < sdisp[0].size(); i++)
		{
			fout << sdisp[0][i] << "\n";
		}
		fout << "\nSCALARS " << "synaptogenesis " << "float 1\nLOOKUP_TABLE default\n";
		for (size_t i = 0; i < sdisp[1].size(); i++)
		{
			fout << sdisp[1][i] << "\n";
		}
		fout << "\nSCALARS " << "tubulin " << "float 1\nLOOKUP_TABLE default\n";
		for (size_t i = 0; i < sdisp[2].size(); i++)
		{
			fout << sdisp[2][i] << "\n";
		}
		fout << "\nSCALARS " << "tips " << "float 1\nLOOKUP_TABLE default\n";
		for (size_t i = 0; i < sdisp[3].size(); i++)
		{
			fout << sdisp[3][i] << "\n";
		}
		fout << "\nSCALARS " << "theta " << "float 1\nLOOKUP_TABLE default\n";
		for (size_t i = 0; i < sdisp[4].size(); i++)
		{
			fout << sdisp[4][i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}


void NeuronGrowth::PointFormValue(vector<float> &Nx, const vector<float> &U, float Value)
{
	Value = 0;

	for (size_t j = 0; j < Nx.size(); j++)
		Value += U[j] * Nx[j];
}

void NeuronGrowth::PointFormGrad(vector<array<float, 2>>& dNdx, const vector<float> &U, float Value[2])
{
	for (int j = 0; j < 2; j++)
		Value[j] = 0.;

	for (int j = 0; j < 2; j++)
		for (size_t k = 0; k < dNdx.size(); k++)
			Value[j] += U[k] * dNdx[k][j];
}

void NeuronGrowth::PointFormHess(vector<array<array<float, 2>, 2>> &d2Ndx2, const vector<float> &U, float Value[2][2])
{
	for (int j = 0; j < 2; j++)
	{
		for (int k = 0; k < 2; k++)
		{
			Value[j][k] = 0.;
		}
	}
	for (int j = 0; j < 2; j++)
	{
		for (int k = 0; k < 2; k++)
		{
			for (size_t l = 0; l < d2Ndx2.size(); l++)
			{
				Value[j][k] += U[l] * d2Ndx2[l][j][k];
			}
		}
	}
}

void NeuronGrowth::ElementValue(const vector<float>& Nx, const vector<float>& value_node, float& value)
{
	value = 0.;
	for (size_t i = 0; i < Nx.size(); i++)
		value += value_node[i] * Nx[i];
}

void NeuronGrowth::ElementValueAll(const vector<float>& Nx, const vector<float>& elePhiGuess, float& elePG,
	const vector<float>& elePhi, float& eleP, const vector<float>& eleSyn, float& eleS,
	const vector<float>& eleTips, float& eleTp, const vector<float>& eleTubulin, float& eleTb,
	const vector<float>& eleEpsilon, float& eleEP, const vector<float>& eleEpsilonP, float& eleEEP)
{
	elePG = 0.;
	eleP = 0.;
	eleS = 0.;
	eleTp = 0.;
	eleTb = 0.;
	eleEP = 0.;
	eleEEP = 0.;
	for (size_t i = 0; i < Nx.size(); i++) {
		elePG += elePhiGuess[i] * Nx[i];
		eleP += elePhi[i] * Nx[i];
		eleS += eleSyn[i] * Nx[i];
		eleTp += eleTips[i] * Nx[i];
		eleTb += eleTubulin[i] * Nx[i];
		eleEP += eleEpsilon[i] * Nx[i];
		eleEEP += eleEpsilonP[i] * Nx[i];
	}
}

void NeuronGrowth::ElementDeriv(const uint& nen, vector<array<float, 2>>& dNdx, const vector<float>& value_node, float& dVdx, float& dVdy)
{
	dVdx = 0;
	dVdy = 0.;
	for (size_t i = 0; i < nen; i++) {
		dVdx += value_node[i] * dNdx[i][0];
		dVdy += value_node[i] * dNdx[i][1];
	}
}


void NeuronGrowth::ElementDerivAll(const uint& nen, vector<array<float, 2>>& dNdx,
	const vector<float>& elePhiGuess, float& dPGdx, float& dPGdy,
	const vector<float>& eleTheta, float& dThedx, float& dThedy, const vector<float>& eleEpsilon,
	float& dAdx, float& dAdy, const vector<float>& eleEpsilonP, float& dAPdx, float& dAPdy)
{	
	dPGdx = 0; dPGdy = 0.;
	dThedx = 0; dThedy = 0.;
	dAdx = 0; dAdy = 0.;
	dAPdx = 0; dAPdy = 0.;
	for (size_t i = 0; i < nen; i++) {
		dPGdx += elePhiGuess[i] * dNdx[i][0];	dPGdy += elePhiGuess[i] * dNdx[i][1];
		dThedx += eleTheta[i] * dNdx[i][0];	dThedy += eleTheta[i] * dNdx[i][1];
		dAdx += eleEpsilon[i] * dNdx[i][0];	dAdy += eleEpsilon[i] * dNdx[i][1];
		dAPdx += eleEpsilonP[i] * dNdx[i][0];	dAPdy += eleEpsilonP[i] * dNdx[i][1];
	}
}

void NeuronGrowth::ElementEvaluationAll_phi(const uint& nen, const vector<float>& Nx, vector<array<float, 2>>& dNdx,
	const vector<float>& elePhiGuess, float& elePG, const vector<float>& elePhi, float& eleP,
	const vector<float>& eleEpsilon, float& eleEP, const vector<float>& eleEpsilonP, float& eleEEP,
	float& dPGdx, float& dPGdy, float& dAdx, float& dAdy, float& dAPdx, float& dAPdy)
{
	elePG = 0.;
	eleP = 0.;
	eleEP = 0.;
	eleEEP = 0.;

	dPGdx = 0; dPGdy = 0.;
	dAdx = 0; dAdy = 0.;
	dAPdx = 0; dAPdy = 0.;

	for (size_t i = 0; i < nen; i++) {
		elePG += elePhiGuess[i] * Nx[i];
		eleP += elePhi[i] * Nx[i];
		eleEP += eleEpsilon[i] * Nx[i];
		eleEEP += eleEpsilonP[i] * Nx[i];

		dPGdx += elePhiGuess[i] * dNdx[i][0];	dPGdy += elePhiGuess[i] * dNdx[i][1];
		dAdx += eleEpsilon[i] * dNdx[i][0];	dAdy += eleEpsilon[i] * dNdx[i][1];
		dAPdx += eleEpsilonP[i] * dNdx[i][0];	dAPdy += eleEpsilonP[i] * dNdx[i][1];
	}
}


void NeuronGrowth::ElementEvaluationAll_phi(const uint& nen, const vector<float>& Nx, vector<array<float, 2>>& dNdx, vector<vector<float>>& eleVal, vector<float>& vars)
{
	// eleVal Array Definitions:
	// Index	Description
	// -----	-----------
	// 0		elePhiGuess	
	// 1		elePhi		
	// 2		eleTheta	
	// 3		eleEpsilon	
	// 4		eleEpsilonP	

	// Variables Definition:
	// Index   Description
	// -----   -----------
	// 0       elePG    - element value for PhiGuess
	// 1       dPGdx    - derivative of PhiGuess w.r.t. x
	// 2       dPGdy    - derivative of PhiGuess w.r.t. y
	// 3       eleP     - element value for Phi
	// 4       eleEP    - element value for epsilonP
	// 5       eleEEP   - element value for epsilonEP
	// 6       dAdx     - derivative of A w.r.t. x
	// 7       dAdy     - derivative of A w.r.t. y
	
	vars[0] = 0.;
	vars[3] = 0.;
	vars[4] = 0.;
	vars[5] = 0.;

	vars[1] = 0; vars[2] = 0.;
	vars[6] = 0; vars[7] = 0.;

	for (size_t i = 0; i < nen; i++) {
		vars[0] += eleVal[0][i] * Nx[i];
		vars[3] += eleVal[1][i] * Nx[i];
		vars[4] += eleVal[3][i] * Nx[i];
		vars[5] += eleVal[4][i] * Nx[i];

		vars[1] += eleVal[0][i] * dNdx[i][0];	vars[2] += eleVal[0][i] * dNdx[i][1];
		vars[6] += eleVal[3][i] * dNdx[i][0];	vars[7] += eleVal[3][i] * dNdx[i][1];
	}
}

void NeuronGrowth::ElementEvaluationAll_phi(const uint& nen, const vector<float>& Nx, vector<array<float, 2>>& dNdx, vector<float>& elePhiGuess, vector<float>& vars)
{
	// eleVal Array Definitions:
	// Index	Description
	// -----	-----------
	// 0		elePhiGuess	
	// 1		elePhi		
	// 2		eleTheta	
	// 3		eleEpsilon	
	// 4		eleEpsilonP	

	// Variables Definition:
	// Index   Description
	// -----   -----------
	// 0       elePG    - element value for PhiGuess
	// 1       dPGdx    - derivative of PhiGuess w.r.t. x
	// 2       dPGdy    - derivative of PhiGuess w.r.t. y
	// 3       eleP     - element value for Phi
	// 4       eleEP    - element value for epsilonP
	// 5       eleEEP   - element value for epsilonEP
	// 6       dAdx     - derivative of A w.r.t. x
	// 7       dAdy     - derivative of A w.r.t. y
	
	vars[0] = 0.;
	// vars[3] = 0.;
	// vars[4] = 0.;
	// vars[5] = 0.;

	vars[1] = 0; vars[2] = 0.;
	// vars[6] = 0; vars[7] = 0.;

	for (size_t i = 0; i < nen; i++) {
		vars[0] += elePhiGuess[i] * Nx[i];
		// vars[3] += eleVal[1][i] * Nx[i];
		// vars[4] += eleEpsilon[i] * Nx[i];
		// vars[5] += eleVal[4][i] * Nx[i];

		vars[1] += elePhiGuess[i] * dNdx[i][0];	vars[2] += elePhiGuess[i] * dNdx[i][1];
		// vars[6] += eleEpsilon[i] * dNdx[i][0];	vars[7] += eleEpsilon[i] * dNdx[i][1];
	}
}

void NeuronGrowth::ElementEvaluationAll_syn_tub(const uint& nen, const vector<float>& Nx, vector<array<float, 2>>& dNdx,
	const vector<float>& elePhiDiff, float& elePf, const vector<float>& eleSyn, float& eleS,
	const vector<float>& elePhi, float& eleP, const vector<float>& elePhiPrev, float& elePprev,
	const vector<float>& eleConct, float& eleC, float& dPdx, float& dPdy)
{
	elePf = 0.;
	eleS = 0.;
	eleP = 0.;
	elePprev = 0.;
	eleC = 0.;

	dPdx = 0; dPdy = 0.;

	for (size_t i = 0; i < nen; i++) {
		elePf += elePhiDiff[i] * Nx[i];
		eleS += eleSyn[i] * Nx[i];

		eleP += elePhi[i] * Nx[i];
		elePprev += elePhiPrev[i] * Nx[i];
		eleC += eleConct[i] * Nx[i];

		dPdx += elePhi[i] * dNdx[i][0];	dPdy += elePhi[i] * dNdx[i][1];
	}
}

void NeuronGrowth::ElementEvaluationAll_syn_tub(const uint& nen, const vector<float>& Nx, vector<array<float, 2>>& dNdx, vector<vector<float>>& eleVal, vector<float>& vars)
{
	// eleVal Array Definitions:
	// Index    Description
	// -----    -----------
	// 0        elePhiDiff  // Difference in Phi values
	// 1        eleSyn      // Synapse/electrical signal strength or similar context
	// 2        elePhi      // Current Phi value
	// 3        elePhiPrev  // Previous Phi value
	// 4        eleConct    // Concentration or a similar quantity

	// var Array Definitions:
	// Index    Variable        Description
	// -----    --------        -----------
	// 0        elePf           // Placeholder for explanation
	// 1        eleS            // Placeholder for explanation
	// 2        eleP            // Placeholder for explanation
	// 3        dPdx            // Derivative of P with respect to x
	// 4        dPdy            // Derivative of P with respect to y
	// 5        elePprev        // Previous P value
	// 6        eleC            // Placeholder for explanation
	// 7        mag_grad_phi0   // Magnitude of the gradient of phi0
	// 8        term_diff       // Diffusion term or similar
	// 9        term_alph       // Alpha term, specific to the model's context
	// 10       term_beta       // Beta term, specific to the model's context
	// 11       term_source     // Source term, specific to the model's context

	vars[0] = 0.;
	vars[1] = 0.;
	vars[2] = 0.;
	vars[5] = 0.;
	vars[6] = 0.;

	vars[3] = 0; vars[4] = 0.;

	for (size_t i = 0; i < nen; i++) {
		vars[0] += eleVal[0][i] * Nx[i];
		vars[1] += eleVal[1][i] * Nx[i];

		vars[2] += eleVal[2][i] * Nx[i];
		vars[5] += eleVal[3][i] * Nx[i];
		vars[6] += eleVal[4][i] * Nx[i];

		vars[3] += eleVal[2][i] * dNdx[i][0];	vars[4] += eleVal[2][i] * dNdx[i][1];
	}
}

void NeuronGrowth::prepareBasis() 
{
	/*Build linear system in each process*/
	float detJ;
	float dudx[2][2];
	
	float dThedx, dThedy;
	sum_grad_phi0_local = 0;

	for (size_t e = 0; e < bzmesh_process.size(); e++) {
		size_t nen = bzmesh_process[e].IEN.size();

		elePhi0.resize(nen);
		eleTheta.resize(nen);
		
		for (size_t i = 0; i < nen; i++) {
			elePhi0[i] = phi_0[bzmesh_process[e].IEN[i]];
			eleTheta[i] = theta[bzmesh_process[e].IEN[i]];
		}

		for (size_t i = 0; i < Gpt.size(); i++) {
			for (size_t j = 0; j < Gpt.size(); j++) {
				BasisFunction(Gpt[i], Gpt[j], bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dudx, detJ);
				detJ = wght[i] * wght[j] * detJ;

				pre_Nx.push_back(Nx);
				pre_dNdx.push_back(dNdx);
				pre_detJ.push_back(detJ);

				ElementDeriv(nen, dNdx, eleTheta, dThedx, dThedy);
				pre_C0.push_back((0.5 + 6 * s_coeff * sqrt(pow(dThedx, 2) + pow(dThedy, 2))));

				ElementDeriv(nen, dNdx, elePhi0, dP0dx, dP0dy);

				pre_mag_grad_phi0.push_back(pow(dP0dx, 2) + pow(dP0dy, 2));
				sum_grad_phi0_local += pow(dP0dx, 2) + pow(dP0dy, 2);

				// ElementValue(Nx, elePhi0, eleP0);
				// pre_term_source.push_back(source_coeff * eleP0);

			}
		}
	}
	MPI_Barrier(PETSC_COMM_WORLD);
}

void NeuronGrowth::preparePhaseField() 
{
	pre_eleEP.clear();
	pre_eleEEP.clear();
	pre_dAdx.clear();
	pre_dAdy.clear();
	pre_eleP.clear();
	pre_eleTh.clear();
	pre_eleMp.clear();
	pre_C1.clear();

	float eleEP(0), eleEEP(0), dAdx(0), dAdy(0), 
		eleP(0), eleTh(0), eleS(0), eleTb(0), eleTp(0), eleE(0);
	uint ind(0);

	for (size_t e = 0; e < bzmesh_process.size(); e++) {
		size_t nen = bzmesh_process[e].IEN.size();

		vector<float> elePhi(nen, 0), eleTheta(nen, 0), eleEpsilon(nen, 0), eleEpsilonP(nen, 0);
		vector<float> eleSyn(nen, 0), eleTubulin(nen, 0), eleTips(nen, 0);

		for (size_t i = 0; i < nen; i++) {

			elePhi[i] = phi[bzmesh_process[e].IEN[i]];			// elePhi
			eleTheta[i] = theta[bzmesh_process[e].IEN[i]];			// eleTheta
			eleSyn[i] = syn[bzmesh_process[e].IEN[i]];			// eleSyn
			eleTubulin[i] = tub[bzmesh_process[e].IEN[i]];			// eleTubulin
			eleTips[i] = tips[bzmesh_process[e].IEN[i]];			// eleTips
		}

		for (size_t i = 0; i < Gpt.size(); i++) {
			for (size_t j = 0; j < Gpt.size(); j++) {

				EvaluateOrientation(nen, pre_Nx[ind], pre_dNdx[ind], elePhi, eleTheta, eleEpsilon, eleEpsilonP);
				ElementValue(pre_Nx[ind], eleEpsilon, eleEP);
				pre_eleEP.push_back(eleEP);
				ElementValue(pre_Nx[ind], eleEpsilonP, eleEEP);
				pre_eleEEP.push_back(eleEEP);
				ElementDeriv(nen, pre_dNdx[ind], eleEpsilonP, dAdx, dAdy);
				pre_dAdx.push_back(dAdx);
				pre_dAdy.push_back(dAdy);

				ElementValue(pre_Nx[ind], elePhi, eleP);
				pre_eleP.push_back(eleP);

				ElementValue(pre_Nx[ind], eleTheta, eleTh);
				pre_eleTh.push_back(eleTh);

				ElementValue(pre_Nx[ind], eleSyn, eleS);
				ElementValue(pre_Nx[ind], eleTubulin, eleTb);
				ElementValue(pre_Nx[ind], eleTips, eleTp);

				// adjust rg (assembly rate) and sg (disassembly rate) based on detected tips
				if (n < 50) {
					eleE = alphaOverPi*atan(gamma * (c_opt - eleS));
					pre_eleMp.push_back(M_neurite);
					// pre_eleMp.push_back(50);
				} else {
					if (eleTp != 0) {
						eleE = alphaOverPi*atan(gamma * Regular_Heiviside_fun(50 * eleTb - 0) * (c_opt - eleS));
						if (eleTp < 0) {
							pre_eleMp.push_back(M_axon);
							// pre_eleMp.push_back(100);
						} else {
							pre_eleMp.push_back(M_neurite);
							// pre_eleMp.push_back(50);
						}
					} else {
						eleE = alphaOverPi*atan(gamma * Regular_Heiviside_fun(r * eleTb - g) * (c_opt - eleS));
						pre_eleMp.push_back(5);
					}
				}

				// calculate C1 variable for phase field energy term
				float C1 = eleE - pre_C0[ind]; // C1 = E - (0.5 + 6 * user->s_coeff * sqrt(pow(dThedx, 2) + pow(dThedy, 2));
				pre_C1.push_back(C1);
										
				ind += 1;
			}
		}
	}
}

void NeuronGrowth::prepareSourceSum() 
{
	/*Build linear system in each process*/
	sum_grad_phi0_local = 0;

	float dPdx(0), dPdy(0);
	// eleP(0)
	uint e, ind(0);
	// cout << comRank << " " << bzmesh_process.size() << endl;

	vector<float> elePhi;
	for (e = 0; e < bzmesh_process.size(); e++) {
		uint nen = bzmesh_process[e].IEN.size();

		elePhi.resize(nen);
		
		for (size_t i = 0; i < nen; i++) {
			elePhi[i] = phi[bzmesh_process[e].IEN[i]];
		}

		for (size_t i = 0; i < Gpt.size(); i++) {
			for (size_t j = 0; j < Gpt.size(); j++) {
				ElementDeriv(nen, pre_dNdx[ind], elePhi, dPdx, dPdy);

				pre_mag_grad_phi0.push_back(pow(dPdx, 2) + pow(dPdy, 2));
				sum_grad_phi0_local += pow(dPdx, 2) + pow(dPdy, 2);

				ind += 1;
			}
		}
	}
}

void NeuronGrowth::prepareTerm_source() 
{
	uint ind(0);
	pre_term_source.clear();
	for (size_t e = 0; e < bzmesh_process.size(); e++) {
		for (size_t i = 0; i < Gpt.size(); i++) {
			for (size_t j = 0; j < Gpt.size(); j++) {
				pre_term_source.push_back(source_coeff * pre_mag_grad_phi0[ind] / sum_grad_phi0_global);
			}
		}
	}
}

float NeuronGrowth::Regular_Heiviside_fun(const float& x) 
{
	// float epsilon = 0.0001; // the number is not fixed.
	float epsilon = 1e-15; // the number is not fixed.
	return 0.5*(1+(2/PI)*atan(x/epsilon));
}

void NeuronGrowth::EvaluateOrientation(const uint& nen, const vector<float>& Nx, const vector<array<float, 2>>& dNdx,
	const vector<float>& elePhi, const vector<float>& eleTheta,  vector<float>& eleEpsilon, vector<float>& eleEpsilonP)
{
	for (size_t i = 0; i < nen; i++) {
		eleEpsilon[i] += epsilonb * (1.0 + delta * cos(aniso * (atan2(elePhi[i] * dNdx[i][0], elePhi[i] * dNdx[i][1]) - eleTheta[i] * Nx[i])));
		eleEpsilonP[i] += -epsilonb * (aniso * delta * sin(aniso * (atan2(elePhi[i] * dNdx[i][0], elePhi[i] * dNdx[i][1]) - eleTheta[i] * Nx[i])));
	}
}

void NeuronGrowth::BuildLinearSystemProcessNG_phi()
{
	/*Build linear system in each process*/
	int ind(0); // pre-calculated variable index
	for (size_t e = 0; e < bzmesh_process.size(); e++) {
		uint nen = bzmesh_process[e].IEN.size(); // 16 supporting cp for 2D case

		vector<vector<float>> EMatrixSolve;	EMatrixSolve.clear();
		vector<float> EVectorSolve;		EVectorSolve.clear();
		EMatrixSolve.resize(nen);
		EVectorSolve.resize(nen);

		vector<float> elePhiGuess(nen, 0);

		for (size_t i = 0; i < nen * 1; i++) {
			// EMatrixSolve[i].reserve(nen*sizeof(float));
			EMatrixSolve[i].resize(nen);
		}
		
		// preparing element stiffness matrix and vector, extract values from control mesh
		for (size_t i = 0; i < nen; i++) {

			EVectorSolve[i] = 0.0;
			for (size_t j = 0; j < nen; j++) {
				EMatrixSolve[i][j] = 0.0;
			}
			elePhiGuess[i] = phi[bzmesh_process[e].IEN[i]];		// elePhiGuess
		}

		// loop through gaussian quadrature points
		for (size_t i = 0; i < Gpt.size(); i++) {
			for (size_t j = 0; j < Gpt.size(); j++) {				
				ElementEvaluationAll_phi(nen, pre_Nx[ind], pre_dNdx[ind], elePhiGuess, vars);
				
				// loop through control points
				for (size_t m = 0; m < nen; m++) {
					// terma2 = - eleEP * eleEP * (dPGdx * dNdx[m][0] + dPGdy * dNdx[m][1]);
					// termadx = - dAdx * eleEEP * dPGdy * dNdx[m][0];
					// termady = - dAdy * eleEEP * dPGdx * dNdx[m][1];
					// termdbl = (- elePG * elePG * elePG + (1 - C1) * elePG * elePG + C1 * elePG) * Nx[m];
					// EVectorSolve[m] += (elePG * Nx[m] - dt * M_phi * (terma2 - termadx + termady + termdbl) - eleP * Nx[m]) * detJ;

					EVectorSolve[m] += (vars[0] * pre_Nx[ind][m] - dt * pre_eleMp[ind] *\
						((- pre_eleEP[ind] * pre_eleEP[ind] * (vars[1] * pre_dNdx[ind][m][0] + vars[2] * pre_dNdx[ind][m][1])) -\
						(- pre_dAdx[ind] * pre_eleEEP[ind] * vars[2] * pre_dNdx[ind][m][0]) +\
						(- pre_dAdy[ind] * pre_eleEEP[ind] * vars[1] * pre_dNdx[ind][m][1]) +\
						(- vars[0] * vars[0] * vars[0] + (1 - pre_C1[ind]) * vars[0] * vars[0] + pre_C1[ind] * vars[0]) * pre_Nx[ind][m])
						- pre_eleP[ind] * pre_Nx[ind][m]) * pre_detJ[ind];

					// loop through 16 control points
					for (size_t n = 0; n < nen; n++) {
						// terma2 = - eleEP * eleEP * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1]);
						// termadx = - dAdx * eleEEP * dNdx[m][1] * dNdx[n][0];
						// termady = - dAdy * eleEEP * dNdx[m][0] * dNdx[n][1];
						// termdbl = (- 3 * elePG * elePG + 2 * (1 - C1) * elePG + C1 * Nx[m]) * Nx[n];
						// EMatrixSolve[m][n] += (Nx[m] * Nx[n] - dt * M_phi * (terma2 - termadx + termady + termdbl)) * detJ;
						
						EMatrixSolve[m][n] += (pre_Nx[ind][m] * pre_Nx[ind][n] - dt * pre_eleMp[ind] *
							((- pre_eleEP[ind] * pre_eleEP[ind] * (pre_dNdx[ind][m][0] * pre_dNdx[ind][n][0] + pre_dNdx[ind][m][1] * pre_dNdx[ind][n][1])) // terma2
							- (- pre_dAdx[ind] * pre_eleEEP[ind] * pre_dNdx[ind][m][1] * pre_dNdx[ind][n][0]) // termadx
							+ (- pre_dAdy[ind] * pre_eleEEP[ind] * pre_dNdx[ind][m][0] * pre_dNdx[ind][n][1]) // termady
							+ (- 3 * vars[0] * vars[0] + 2 * (1 - pre_C1[ind]) * vars[0] + pre_C1[ind] * pre_Nx[ind][m]) * pre_Nx[ind][n]) // termdbl
							) * pre_detJ[ind];
						// EMatrixSolve[m][n] += (pre_Nx[ind][m] * pre_Nx[ind][n]) * pre_detJ[ind];
					}
				}
				ind += 1; // incrementing index for extracting pre-calculated variables
			}
		}

		/*Apply Boundary Condition*/
		for (size_t i = 0; i < nen; i++) {
			int A = bzmesh_process[e].IEN[i];
			if (cpts[A].label == 1) { // domain boundary
				ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
			}
		}

		ResidualAssembly(EVectorSolve, bzmesh_process[e].IEN, GR_phi);
		MatrixAssembly(EMatrixSolve, bzmesh_process[e].IEN, GK_phi);
	}

	if (comRank == 0) {
		check_itr += 1;
	}

	VecAssemblyBegin(GR_phi);
	MatAssemblyBegin(GK_phi, MAT_FINAL_ASSEMBLY);
}

void NeuronGrowth::BuildLinearSystemProcessNG_syn(const vector<Element2D>& tmesh, const vector<Vertex2D>& cpts)
{
	/*Build linear system in each process*/
	for (size_t e = 0; e < bzmesh_process.size(); e++) {
		float detJ;
		float dudx[2][2];
		uint nen = bzmesh_process[e].IEN.size();

		vector<vector<float>> EMatrixSolve;
		vector<float> EVectorSolve;
		vector<float> Nx;
		vector<array<float, 2>> dNdx;

		EMatrixSolve.clear();
		EVectorSolve.clear();
		Nx.clear();
		dNdx.clear();

		EMatrixSolve.resize(nen * 1);
		EVectorSolve.resize(nen * 1);
		Nx.resize(nen);
		dNdx.resize(nen);

		for (size_t i = 0; i < nen; i++) {
			EMatrixSolve[i].resize(nen * 1, 0.0);
		}

		for (size_t i = 0; i < nen * 1; i++) {
			for (size_t j = 0; j < nen * 1; j++)
			{
				EMatrixSolve[i][j] = 0.0;
			}
			EVectorSolve[i] = 0.0;
		}

		vector<float> elePhiDiff, eleSyn;
		elePhiDiff.resize(nen);
		eleSyn.resize(nen);
		float elePf, eleS;

		for (size_t i = 0; i < nen; i++) {
			elePhiDiff[i] = abs(phi[bzmesh_process[e].IEN[i]] - phi_prev[bzmesh_process[e].IEN[i]]);
			eleSyn[i] = syn[bzmesh_process[e].IEN[i]];
		}

		for (size_t i = 0; i < Gpt.size(); i++) {
			for (size_t j = 0; j < Gpt.size(); j++) {
				BasisFunction(Gpt[i], Gpt[j], bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dudx, detJ);
				detJ = wght[i] * wght[j] * detJ;
				ElementValue(Nx, elePhiDiff, elePf);
				ElementValue(Nx, eleSyn, eleS);
				if (judge_syn == 0)
					Tangent_syn(nen, Nx, dNdx, detJ, EMatrixSolve);
				Residual_syn(nen, Nx, dNdx, detJ, elePf, eleS, EVectorSolve);
			}
		}

		/*Apply Boundary Condition*/
		for (size_t i = 0; i < nen; i++) {
			int A = bzmesh_process[e].IEN[i];
			if (cpts[A].label == 1)
				ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
		}

		ResidualAssembly(EVectorSolve, bzmesh_process[e].IEN, GR_syn);
		if (judge_syn == 0) // The matrix is the same, so only need to assembly once
			MatrixAssembly(EMatrixSolve, bzmesh_process[e].IEN, GK_syn);
	}
	VecAssemblyBegin(GR_syn);
	if (judge_syn == 0)
		MatAssemblyBegin(GK_syn, MAT_FINAL_ASSEMBLY);
}

void NeuronGrowth::Tangent_syn(const uint nen, vector<float> &Nx, vector<array<float, 2>>& dNdx, float detJ, vector<vector<float>> &EMatrixSolve)
{
	/*calculate tangent matrix*/
	for (size_t i = 0; i < nen; i++) {
		for (size_t j = 0; j < nen; j++) {
			EMatrixSolve[i][j] += (Nx[i] * Nx[j] + dt * Dc * (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1])) * detJ;
			// EMatrixSolve[i][j] += (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1]) * detJ;	// steady state heat transfer
		}
	}
}

void NeuronGrowth::Residual_syn(const uint nen, vector<float> &Nx, vector<array<float, 2>>& dNdx, float detJ, float elePf, float eleS, vector<float> &EVectorSolve)
{
	/*calculate residual of the equation*/
	for (size_t i = 0; i < nen; i++) {
		EVectorSolve[i] += Nx[i] * (eleS + kappa * elePf) * detJ;
		// EVectorSolve[i] += Nx[i] * (eleS) * detJ; // transient heat transfer
		// EVectorSolve[i] += 0.00; // steady state heat transfer
	}
}

void NeuronGrowth::CalculateSumGradPhi0(const vector<Element2D>& tmesh, const vector<Vertex2D>& cpts)
{
	/*Build linear system in each process*/
	for (size_t e = 0; e < bzmesh_process.size(); e++) {
		float detJ;
		float dudx[2][2];
		uint nen = bzmesh_process[e].IEN.size();

		vector<float> Nx;
		vector<array<float, 2>> dNdx;

		elePhi0.resize(nen);
		for (size_t i = 0; i < nen; i++) {
			elePhi0[i] = phi_0[bzmesh_process[e].IEN[i]];
		}

		for (size_t i = 0; i < Gpt.size(); i++) {
			for (size_t j = 0; j < Gpt.size(); j++) {
				BasisFunction(Gpt[i], Gpt[j], bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dudx, detJ);
				detJ = wght[i] * wght[j] * detJ;
				ElementDeriv(nen, dNdx, elePhi0, dP0dx, dP0dy);
				sum_grad_phi0_local += pow(dP0dx, 2) + pow(dP0dy, 2);
			}
		}
	}
}

void NeuronGrowth::BuildLinearSystemProcessNG_tub(const vector<Element2D>& tmesh, const vector<Vertex2D>& cpts)
{
	/*Build linear system in each process*/
	for (size_t e = 0; e < bzmesh_process.size(); e++) {
		float detJ;
		float dudx[2][2];
		uint nen = bzmesh_process[e].IEN.size();

		vector<vector<float>> EMatrixSolve;
		vector<float> EVectorSolve;
		vector<float> Nx;
		vector<array<float, 2>> dNdx;

		EMatrixSolve.clear();
		EVectorSolve.clear();
		Nx.clear();
		dNdx.clear();

		EMatrixSolve.resize(nen);
		EVectorSolve.resize(nen);
		Nx.resize(nen);
		dNdx.resize(nen);

		for (size_t i = 0; i < nen; i++) {
			EMatrixSolve[i].resize(nen, 0.0);
		}

		for (size_t i = 0; i < nen * 1; i++) {
			for (size_t j = 0; j < nen * 1; j++) {
				EMatrixSolve[i][j] = 0.0;
			}
			EVectorSolve[i] = 0.0;
		}

		vector<float> elePhi, elePhiPrev, eleConct;
		elePhi.resize(nen);
		elePhiPrev.resize(nen);
		eleConct.resize(nen);
		float eleP, dPdx, dPdy, elePprev, eleC, mag_grad_phi0;

		for (size_t i = 0; i < nen; i++) {
			elePhi[i] = phi[bzmesh_process[e].IEN[i]];
			elePhiPrev[i] = phi_prev[bzmesh_process[e].IEN[i]];
			eleConct[i] = tub[bzmesh_process[e].IEN[i]];
		}

		for (size_t i = 0; i < Gpt.size(); i++) {
			for (size_t j = 0; j < Gpt.size(); j++) {
				BasisFunction(Gpt[i], Gpt[j], bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dudx, detJ);
				detJ = wght[i] * wght[j] * detJ;
				ElementValue(Nx, elePhi, eleP);
				ElementValue(Nx, elePhiPrev, elePprev);
				ElementValue(Nx, eleConct, eleC);
				ElementDeriv(nen, dNdx, elePhi, dPdx, dPdy);

				ElementDeriv(nen, dNdx, elePhi0, dP0dx, dP0dy);
				mag_grad_phi0 = pow(dP0dx, 2) + pow(dP0dy, 2);
				
				Tangent_tub(nen, Nx, dNdx, detJ, eleC, eleP, dPdx, dPdy, elePprev, mag_grad_phi0, EMatrixSolve);
				Residual_tub(nen, Nx, dNdx, detJ, eleC, eleP, elePprev, mag_grad_phi0, EVectorSolve);
			}
		}

		/*Apply Boundary Condition*/
		for (size_t i = 0; i < nen; i++) {
			int A = bzmesh_process[e].IEN[i];
			if (cpts[A].label == 1)
				ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
		}

		ResidualAssembly(EVectorSolve, bzmesh_process[e].IEN, GR_tub);
		// if (judge_tub == 0) // The matrix is the same, so only need to assembly once
		MatrixAssembly(EMatrixSolve, bzmesh_process[e].IEN, GK_tub);
	}

	VecAssemblyBegin(GR_tub);
	// if (judge_tub == 0)
	MatAssemblyBegin(GK_tub, MAT_FINAL_ASSEMBLY);
}

void NeuronGrowth::Tangent_tub(const uint nen, vector<float> &Nx, vector<array<float, 2>>& dNdx, float detJ,
	float eleC, float eleP, float dPdx, float dPdy, float elePprev, float mag_grad_phi0, vector<vector<float>> &EMatrixSolve)
{
	float term_diff, term_alph, term_beta;
	/*calculate tangent matrix*/
	for (size_t i = 0; i < nen; i++) {
		for (size_t j = 0; j < nen; j++) {
			term_diff = - Diff * (eleP * (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1]));
			term_alph = alphaT * (eleP * (dNdx[i][0] + dNdx[i][1]) + Nx[i] * (dPdx + dPdy)) * Nx[j];
			term_beta = betaT * (eleP * Nx[i]) * Nx[j];
			
			EMatrixSolve[i][j] += (Nx[i] * Nx[j] - dt / eleP * (term_diff - term_alph - term_beta)) * detJ;
			// EMatrixSolve[i][j] += ((2 * eleP - elePprev) * Nx[i] * Nx[j] - dt/2 * (term_diff - term_alph - term_beta + term_source)) * detJ;
		}
	}
}

void NeuronGrowth::Residual_tub(const uint nen, vector<float> &Nx, vector<array<float, 2>>& dNdx, float detJ,
	float eleC, float eleP, float elePprev, float mag_grad_phi0, vector<float> &EVectorSolve)
{
	float term_source;
	/*calculate residual of the equation*/
	for (size_t i = 0; i < nen; i++) {
		term_source = source_coeff * mag_grad_phi0 / sum_grad_phi0_global;
		EVectorSolve[i] += (dt / eleP * term_source + eleC) * Nx[i] * detJ;
		// EVectorSolve[i] += eleP * eleC * Nx[i] * detJ;
	}
}

void NeuronGrowth::BuildLinearSystemProcessNG_syn_tub(const vector<Element2D>& tmesh, const vector<Vertex2D>& cpts)
{
	/*Build linear system in each process*/
	int ind(0);
	for (size_t e = 0; e < bzmesh_process.size(); e++) {
		uint nen = bzmesh_process[e].IEN.size();

		vector<vector<float>> EMatrixSolve_syn, EMatrixSolve_tub;
		vector<float> EVectorSolve_syn, EVectorSolve_tub;

		EMatrixSolve_syn.clear();
		EVectorSolve_syn.clear();
		EMatrixSolve_tub.clear();
		EVectorSolve_tub.clear();

		EMatrixSolve_syn.resize(nen * 1);
		EVectorSolve_syn.resize(nen * 1);
		EMatrixSolve_tub.resize(nen);
		EVectorSolve_tub.resize(nen);

		vector<vector<float>> eleVal_st;
		eleVal_st.clear();
		eleVal_st.resize(5);

		for (size_t i = 0; i < nen * 1; i++) {
			EMatrixSolve_syn[i].resize(nen * 1, 0.0);
			EMatrixSolve_tub[i].resize(nen, 0.0);
			for (size_t j = 0; j < nen * 1; j++)
			{
				EMatrixSolve_syn[i][j] = 0.0;
				EMatrixSolve_tub[i][j] = 0.0;
			}
			EVectorSolve_syn[i] = 0.0;
			EVectorSolve_tub[i] = 0.0;

			if (i < 5) 
				eleVal_st[i].resize(nen);
		}
		
		for (size_t i = 0; i < nen * 1; i++) {
			eleVal_st[0][i] = phi[bzmesh_process[e].IEN[i]] - phi_prev[bzmesh_process[e].IEN[i]];		// elePhiDiff
			eleVal_st[1][i] = syn[bzmesh_process[e].IEN[i]];						// eleSyn
			eleVal_st[2][i] = phi[bzmesh_process[e].IEN[i]];						// elePhi
			eleVal_st[3][i] = phi_prev[bzmesh_process[e].IEN[i]];						// elePhiPrev
			eleVal_st[4][i] = tub[bzmesh_process[e].IEN[i]];						// eleConct
		}


		for (size_t i = 0; i < Gpt.size(); i++) {
			for (size_t j = 0; j < Gpt.size(); j++) {
				
				vector<float> vars_st;
				vars_st.clear();
				vars_st.resize(12);
				//        0      1     2      3	  4      5        6          7          8          9         10          11
				// float elePf, eleS, eleP, dPdx, dPdy, elePprev, eleC, mag_grad_phi0, term_diff, term_alph, term_beta, term_source;

				float detJ;
				vector<float> Nx;
				vector<array<float, 2>> dNdx;
				Nx.clear();
				dNdx.clear();
				Nx.resize(nen);
				dNdx.resize(nen);

				Nx = pre_Nx[ind];
				dNdx = pre_dNdx[ind];
				detJ = pre_detJ[ind];
				vars_st[7] = pre_mag_grad_phi0[ind];
				vars_st[11] = pre_term_source[ind];
				// vars_st[11] = 0.5;
				ind += 1;

				ElementEvaluationAll_syn_tub(nen, Nx, dNdx, eleVal_st, vars_st);

				for (size_t m = 0; m < nen; m++) {
					// EVectorSolve_syn[m] += Nx[m] * (vars_st[1] + kappa * vars_st[0]) * detJ;
					EVectorSolve_syn[m] += Nx[m] * (vars_st[1] + (kappa - kp75 * vars_st[1]) * vars_st[0] - k2 * vars_st[1]) * detJ;

					// EVectorSolve_tub[m] += (dt / vars_st[2] * vars_st[11] + vars_st[6]) * Nx[m] * detJ;
					// EVectorSolve_tub[m] += (dt / vars_st[2] + vars_st[6]) * Nx[m] * detJ;
					EVectorSolve_tub[m] += (vars_st[2] * vars_st[6] + (dt/4) * vars_st[11]) * Nx[m] * detJ;

					for (size_t n = 0; n < nen; n++) {
						if (judge_syn == 0)

						EMatrixSolve_syn[m][n] += (Nx[m] * Nx[n] + (dt/4) * Dc * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1])) * detJ;

						// EMatrixSolve_tub[m][n] += (Nx[m] * Nx[n] - dt / vars_st[2] * (
						// 	(- Diff * (vars_st[2] * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1])))
						// 	- (alphaT * (vars_st[2] * (dNdx[m][0] + dNdx[n][1]) + Nx[m] * (vars_st[3] + vars_st[4])) * Nx[n])
						// 	- (betaT * (vars_st[2] * Nx[m]) * Nx[n]))
						// 	) * detJ;

						EMatrixSolve_tub[m][n] += ((vars_st[0] * Nx[m] + vars_st[2] * Nx[m]) * Nx[n] - (dt/4) * (
							(- Diff * (vars_st[2] * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1])))
							- (alphaT * (vars_st[2] * (dNdx[m][0] + dNdx[n][1]) + Nx[m] * (vars_st[3] + vars_st[4])) * Nx[n])
							- (betaT * (vars_st[2] * Nx[m]) * Nx[n]))
							) * detJ;

						// EMatrixSolve_tub[m][n] += (Nx[m] * Nx[n] - dt * (
						// 	(- Diff * (vars_st[2] * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1])))
						// 	- (alphaT * (vars_st[2] * (dNdx[m][0] + dNdx[n][1]) + Nx[m] * (vars_st[3] + vars_st[4])) * Nx[n])
						// 	- (betaT * (vars_st[2] * Nx[m]) * Nx[n]))
						// 	) * detJ;
						
					}
				}
			}
		}

		/*Apply Boundary Condition*/
		for (size_t i = 0; i < nen; i++) {
			int A = bzmesh_process[e].IEN[i];
			if (cpts[A].label == 1) // domain boundary
				ApplyBoundaryCondition(0, i, 0, EMatrixSolve_syn, EVectorSolve_syn);
			// if (round(phi_0[A]) > 0)
			// 	ApplyBoundaryCondition(tub[A] + 2e-5 * tub[A], i, 0, EMatrixSolve_tub, EVectorSolve_tub);
				// ApplyBoundaryCondition(tub_0[A], i, 0, EMatrixSolve_tub, EVectorSolve_tub);
			if (round(phi[A]) == 0) // outside soma bc for tubulin
				ApplyBoundaryCondition(0, i, 0, EMatrixSolve_tub, EVectorSolve_tub);
		}

		ResidualAssembly(EVectorSolve_syn, bzmesh_process[e].IEN, GR_syn);
		if (judge_syn == 0) // The matrix is the same, so only need to assembly once
			MatrixAssembly(EMatrixSolve_syn, bzmesh_process[e].IEN, GK_syn);

		ResidualAssembly(EVectorSolve_tub, bzmesh_process[e].IEN, GR_tub);
		MatrixAssembly(EMatrixSolve_tub, bzmesh_process[e].IEN, GK_tub);
	}

	VecAssemblyBegin(GR_syn);
	if (judge_syn == 0)
		MatAssemblyBegin(GK_syn, MAT_FINAL_ASSEMBLY);

	VecAssemblyBegin(GR_tub);
	MatAssemblyBegin(GK_tub, MAT_FINAL_ASSEMBLY);
}

int NeuronGrowth::CheckExpansion(vector<float> input, int NX, int NY) 
{
	// 0 - left | 1 - top | 2 - right | 3 - bottom | 4 - no action
	for (int i = 0; i < 10; i++) {
		for (int j = 1; j < NY; j++) {
			if (round(input[i * (NY+1) + j]) > 0) {
				return 0;
			} else if (round(input[(NX+1 - i - 1) * (NY+1) + j]) > 0) {
				return 2;
			}
		}
		for (int j = 1; j < NX; j++) {
			if (round(input[j * (NY+1) + i]) > 0) {
				return 3;
			} else if (round(input[j * (NY+1) - i]) > 0) {
				return 1;
			}
		}
	}
	return 4;
}

void NeuronGrowth::PopulateRandom(vector<float> &input) 
{

	uint length = input.size();
	for (size_t i = 0; i < length; i++) {
		if (input[i] == 0)
			input[i] = (float)(rand()%100)/(float)100;
	}
}

// Function to interpolate values for new control points
vector<float> NeuronGrowth::InterpolateValues_closest(
	const vector<float>& phi_in,
	const vector<Vertex2D>& cpt,
	const vector<Vertex2D>& cpt_out) 
{

	vector<float> interpolatedValues(cpt_out.size());

	for (size_t i = 0; i < cpt_out.size(); ++i) {
		float minDistance = numeric_limits<float>::max();
		size_t closestIndex = 0;

		// Find the closest point in cpt for each point in cpt_out
		for (size_t j = 0; j < cpt.size(); ++j) {
			float distance = distanceTo(cpt_out[i], cpt[j]);
			if (distance < minDistance) {
				minDistance = distance;
				closestIndex = j;
			}
		}

		// Assume the value of the closest point
		// if (abs(round(phi_in[closestIndex])) >= 0.5) {
		// 	interpolatedValues[i] = 1;
		// } else {
		// 	interpolatedValues[i] = 0;
		// }
		interpolatedValues[i] = phi_in[closestIndex];
	}

	return interpolatedValues;
}

vector<float> NeuronGrowth::InterpolateValues_closest(const vector<float>& input, const KDTree& kdTree, const vector<Vertex2D>& cpt_out) 
{

	vector<float> interpolatedValues(cpt_out.size());
	for (size_t i = 0; i < cpt_out.size(); ++i) {
		const float query_pt[2] = {cpt_out[i].coor[0], cpt_out[i].coor[1]};
		size_t closestIndex;
		float out_dist_sqr;
		nanoflann::KNNResultSet<float> resultSet(1);
		resultSet.init(&closestIndex, &out_dist_sqr);

		// Performing the search with correct query points
		kdTree.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(0));

		// Validate closestIndex is within bounds
		if(closestIndex >= 0 && closestIndex < input.size()) {
			interpolatedValues[i] = input[closestIndex];
			// CellBoundary(input[closestIndex], 0.25);
		} else {
			// Handle error or unexpected case
			interpolatedValues[i] = 0; // Define defaultValue appropriately
			PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck0)! x: %f y: %f\n", cpt_out[i].coor[0], cpt_out[i].coor[1]);

		}
	}

	return interpolatedValues;
}

vector<float> NeuronGrowth::InterpolateValues_closest_pushBoundary(const vector<float>& input, const KDTree& kdTree, const vector<Vertex2D>& cpt_out) 
{

	vector<float> interpolatedValues(cpt_out.size());
	for (size_t i = 0; i < cpt_out.size(); ++i) {
		const float query_pt[2] = {cpt_out[i].coor[0], cpt_out[i].coor[1]};
		size_t closestIndex;
		float out_dist_sqr;
		nanoflann::KNNResultSet<float> resultSet(1);
		resultSet.init(&closestIndex, &out_dist_sqr);

		// Performing the search with correct query points
		kdTree.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(0));

		// Validate closestIndex is within bounds
		if(closestIndex >= 0 && closestIndex < input.size()) {
			interpolatedValues[i] = CellBoundary(input[closestIndex], 0.25);
			// CellBoundary(input[closestIndex], 0.25);
		} else {
			// Handle error or unexpected case
			interpolatedValues[i] = 0; // Define defaultValue appropriately
			PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck0)! x: %f y: %f\n", cpt_out[i].coor[0], cpt_out[i].coor[1]);

		}
	}

	return interpolatedValues;
}

vector<float> NeuronGrowth::InterpolateVars_coarseKDtree(vector<float> input, vector<Vertex2D> cpts_initial, const KDTree& kdTree_initial, const vector<Vertex2D>& cpts, int type, int isTheta) 
{	
	vector<float> output;
	output.resize(cpts.size(), 0);
	
	for (size_t i = 0; i < cpts.size(); i++) {
		float x = round5(cpts[i].coor[0]);
		float y = round5(cpts[i].coor[1]);
		int ind;
		if (KD_SearchPair(cpts_initial, kdTree_initial, x, y, ind)) {
			output[i] = input[ind];
		} else if ((x > prev_max_x) || (x < prev_min_x) || (y > prev_max_y) || (y < prev_max_y)) {
			if (isTheta != 1) {
				output[i] = 0;
			} else {
				output[i] = (float)(rand()%100)/(float)100; // for theta
			}
		} else {
			int indDown, indUp, indLeft, indRight;
			if ((abs(remainder(x,2)) == 1) && (abs(remainder(y,2)) != 1)) {
				if (KD_SearchPair(cpts_initial, kdTree_initial, floorf(x)-1, round5(y), indDown) &&
				KD_SearchPair(cpts_initial, kdTree_initial, floorf(x)+1, round5(y), indUp)) {
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
				if (KD_SearchPair(cpts_initial, kdTree_initial, round5(x), floorf(y)-1, indLeft) &&
				KD_SearchPair(cpts_initial, kdTree_initial, round5(x), floorf(y)+1, indRight)) {
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
				if (KD_SearchPair(cpts_initial, kdTree_initial, floorf(x)-1, floorf(y)-1, indDown) &&
				KD_SearchPair(cpts_initial, kdTree_initial, floorf(x)+1, floorf(y)+1, indUp) &&
				KD_SearchPair(cpts_initial, kdTree_initial, floorf(x)-1, floorf(y)+1, indLeft) &&
				KD_SearchPair(cpts_initial, kdTree_initial, floorf(x)+1, floorf(y)-1, indRight)) {
					if (type == 0) {
						output[i] = max(max(input[indDown], input[indUp]), max(input[indLeft], input[indRight]));
					} else if (type == 1) {
						output[i] = (input[indDown] + input[indUp] + input[indLeft] + input[indRight])/4;
					} else if (type == 2) {
						output[i] = 0;
					}
				} else {
					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck2)! x: %f y: %f\n", x, y);
				}
			}
		}			
	}
	return output;
}

bool NeuronGrowth::KD_SearchPair(const vector<Vertex2D> prev_cpts, const KDTree& kdTree, float targetX, float targetY, int &ind) 
{

	const float query_pt[2] = {targetX, targetY};
	size_t closestIndex;
	float out_dist_sqr;
	nanoflann::KNNResultSet<float> resultSet(1);
	resultSet.init(&closestIndex, &out_dist_sqr);

	// Performing the search with correct query points
	kdTree.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(0));

	float x = prev_cpts[closestIndex].coor[0];
	float y = prev_cpts[closestIndex].coor[1];
	
	// if ((max(abs(x-targetX), abs(y-targetY)) <= 1) || (targetX < prev_min_x) || (targetY < prev_min_y) || (targetX > prev_max_x) || (targetY > prev_max_y)) {
	if (max(abs(x-targetX), abs(y-targetY)) <= 1) {
		ind = static_cast<int>(closestIndex);
		return true; // Found the pair (targetX, targetY) in the vector
	} else {
		// PetscPrintf(PETSC_COMM_WORLD, "Failed search! x: %f y: %f | %f %f\n", x, y, targetX, targetY);
		return false; // Pair not found in the vector
	}
}

float NeuronGrowth::RmOutlier(const vector<float>& data) 
{
	float sum = 0.0, mean, standardDeviation = 0.0;

	for (float val : data) {
		sum += val;
	}

	mean = sum / data.size();

	for (float val : data) {
		standardDeviation += pow(val - mean, 2);
	}

	standardDeviation = sqrt(standardDeviation / data.size());

	// // cout << mean + 2*standardDeviation << " " << maxVal << endl;
	// for (size_t i = 0; i < data.size(); i++) {
	// 	if (data[i] > (mean + 3.29*standardDeviation)) {
	// 		// data[i] = mean + 3.29*standardDeviation;
	// 	} else {
	// 		data[i] = 0;
	// 	}
	// }

	return (mean + 3.29 * standardDeviation);
}

float NeuronGrowth::CellBoundary(const float& phi, const float& threshold) 
{
	return (phi > threshold) ? 1.0f : 0.0f;
}

// Check if point coordinates are within the bounds defined by center and dx dy dz
bool NeuronGrowth::isInBox(const Vertex2D& point, const Vertex2D& center, float dx, float dy) 
{
	return point.coor[0] >= (center.coor[0] - dx/2) && point.coor[0] <= (center.coor[0] + dx/2) &&
		point.coor[1] >= (center.coor[1] - dy/2) && point.coor[1] <= (center.coor[1] + dy/2);
}

// Calculates the sum of phi values around each vertex in a given set, applying
// a threshold to identify significant points, potentially indicating neuron growth tips.
vector<float> NeuronGrowth::calculatePhiSum(const vector<Vertex2D>& cpts, float dx, float dy, vector<float> id) 
{
	vector<float> tp(cpts.size(), 0); // Initialize the result vector with zeros
	tp.clear(); tp.resize(cpts.size());
	float threshold = 0.9; // Threshold for filtering tp
	float maxVal = 0; // Track the maximum value of tp for normalization

	// Iterate over each center point
	for (size_t i = 0; i < cpts.size(); ++i) {
		const auto& center = cpts[i];

		tp[i] = 0;

		// Sum phi values for points within the box centered at 'center'
		for (size_t j = 0; j < phi.size(); ++j) {
			if (isInBox(cpts[j], center, dx, dy)) {
				tp[i] += CellBoundary(phi[j], 0.5);
				// , id[j];
			}
		}

		tp[i] = CellBoundary(phi[i], 0.5) / tp[i] * CellBoundary(phi[i], 0.5);
		// Handle NaN cases, ensuring a valid tp value
		if (isnan(tp[i])) {
			tp[i] = 0;
		}
		
		if ((center.coor[0] <= min_x+1) || (center.coor[0] >= max_x-1) 
			|| (center.coor[1] <= min_y+1) || (center.coor[1] >= max_y-1)) {
			tp[i] = 0;
		}

		// Update maxVal to the highest tp value found
		if (CellBoundary(phi[i], 0.5) > 0)
			maxVal = max(tp[i], maxVal);
	}

	for (size_t i = 0; i < tp.size(); i++) {
		tp[i] = tp[i]/maxVal;
	}

	// maxVal = RmOutlier(tp);

	// Thresholding and setting tp values based on the maximum value found
	for (size_t i = 0; i < tp.size(); ++i) {
		// if (tp[i] > (threshold * maxVal)) {
		if (tp[i] > (threshold)) {
		// if (tp[i] > (0.018)) {
			tp[i] = 1; // Set tp below threshold to 0
		} else {
			tp[i] = 0; // Set tp above threshold to 1
		}
	}
	return tp;
}

void NeuronGrowth::DetectTipsMulti(const vector<float>& phi_fine, const vector<int>& id, const int& numNeuron, vector<float>& phiSum, const int& rows, const int& cols) 
{
	// const int length = (NX + 1) * (NY + 1);
	// phiSum.assign(length, 0); // Clear and resize with zero initialization

	// const float threshold = 0.25; // threshould for where to detect neurites
	// const float intensityThreshold = 0.1; // threshould for whether to use the value for detection
	// // const float threshold = 0.5;
	// // const float intensityThreshold = 0.5;
	// const int offsetY = NY + 1; // offset to avoid index issue
	// const float phiThreshold = 0; // threshould for scaling

	// // Lambda to replace CellBoundary and reduce condition checks - mostly for optimization
	// auto cellBoundary = [](float phi, float threshold) -> float {
	// 	return phi > threshold ? 1.0f : 0.0f;
	// };

	// // Iterate with precomputed values and reduced condition checks
	// for (int i = (5 * NY + 5); i < (length - 4 * NY - 4); ++i) {
	// 	if (cellBoundary(phi_fine[i], threshold) > 0 && round(id[i]) != 9) {
	// 		for (int j = -4; j <= 4; ++j) {
	// 			for (int k = -4; k <= 4; ++k) {
	// 				int index = i + j * offsetY + k;
	// 				float roundedId = round(id[index]);
	// 				// if (roundedId == round(id[i]) || roundedId == 0 || roundedId == 9) {
	// 				if (roundedId == round(id[i])) {
	// 					phiSum[i] += cellBoundary(phi_fine[index], intensityThreshold);
	// 				}
	// 			}
	// 		}
	// 		phiSum[i] = phiSum[i] > 0 ? cellBoundary(phi_fine[i], phiThreshold) / phiSum[i] : 0;
	// 		if (isnan(phiSum[i])) phiSum[i] = 0; // Handle NaN explicitly, though it should not occur now
	// 	}
	// }

	// // // Iterate with precomputed values and reduced condition checks
	// // for (int i = (5 * NY + 5); i < (length - 4 * NY - 4); ++i) {
	// // 	if (cellBoundary(phi_fine[i], threshold) > 0 && round(id[i]) != 9) {
	// // 		const int radius = 4; // Define the radius of the circle
	// // 		for (int j = -radius; j <= radius; ++j) {
	// // 			for (int k = -radius; k <= radius; ++k) {
	// // 				if (j * j + k * k <= radius * radius) { // Check if the point (j, k) is inside the circle
	// // 					int index = i + j * offsetY + k;
	// // 					if (index >= 0 && index < length) { // Ensure index is within bounds
	// // 						int roundedId = round(id[index]);
	// // 						if (roundedId == round(id[i])) {
	// // 							phiSum[i] += cellBoundary(phi_fine[index], intensityThreshold);
	// // 						}
	// // 					}
	// // 				}
	// // 			}
	// // 		}
	// // 		phiSum[i] = phiSum[i] > 0 ? cellBoundary(phi_fine[i], phiThreshold) / phiSum[i] : 0;
	// // 		if (isnan(phiSum[i])) phiSum[i] = 0; // Handle NaN explicitly, though it should not occur now
	// // 	}
	// // }

	// // Normalize phiSum by the standard deviation value
	// float stdVal = RmOutlier(phiSum);
	// transform(phiSum.begin(), phiSum.end(), phiSum.begin(), [stdVal](float val) { return val / stdVal; });

	// // Find the maximum value in phiSum for potential use
	// float maxVal = *max_element(phiSum.begin(), phiSum.end());

	// // Filter based on a normalized threshold
	// // const float normalizedThreshold = 1.4;
	// // const float normalizedThreshold = 1.3; 
	// // const float normalizedThreshold = 1.2;//
	// // const float normalizedThreshold = 1.1;
	// // const float normalizedThreshold = 1.0;
	// const float normalizedThreshold = 0.8*maxVal;
	// // const float normalizedThreshold = 0.75*maxVal;
	// // cout << normalizedThreshold << endl;

	// // float normalizedThreshold;
	// // if (n<10000) {
	// // 	normalizedThreshold = 0.8*maxVal;
	// // } else {
	// // 	normalizedThreshold = 0.6*maxVal;
	// // }
	// transform(phiSum.begin(), phiSum.end(), phiSum.begin(), [normalizedThreshold](float val) {
	// 	return val < normalizedThreshold ? 0.0f : val;
	// });

	const int length = rows * cols;
	phiSum.assign(length, 0); // Clear and resize with zero initialization

	const float threshold = 0.25; // Threshold for where to detect neurites
	const float intensityThreshold = 0.1; // Threshold for whether to use the value for detection
	const int offsetY = cols; // Offset to avoid index issue
	const float phiThreshold = 0; // Threshold for scaling

	// Lambda to replace CellBoundary and reduce condition checks - mostly for optimization
	auto cellBoundary = [](float phi, float threshold) -> float {
		return phi > threshold ? 1.0f : 0.0f;
	};

	// Using a map to store max values for each ID
	map<int, float> maxValues;

	// Iterate with precomputed values and reduced condition checks
	for (int i = (5 * cols); i < (length - 4 * cols); ++i) {
		if (cellBoundary(phi_fine[i], threshold) > 0 && round(id[i]) != 9) {
			for (int j = -4; j <= 4; ++j) {
				for (int k = -4; k <= 4; ++k) {
					int index = i + j * offsetY + k;
					int roundedId = round(id[index]);
					if (roundedId == round(id[i]) || roundedId < 0) {
						phiSum[i] += cellBoundary(phi_fine[index], intensityThreshold);
					}
				}
			}
			phiSum[i] = phiSum[i] > 0 ? cellBoundary(phi_fine[i], phiThreshold) / phiSum[i] : 0;
			if (isnan(phiSum[i])) phiSum[i] = 0; // Handle NaN explicitly, though it should not occur now

			// Update max value for each ID
			int idKey = static_cast<int>(round(id[i]));
			if (maxValues.find(idKey) == maxValues.end() || maxValues[idKey] < phiSum[i]) {
				maxValues[idKey] = phiSum[i];
			}
		}
	}

	// Normalize phiSum by the max value for each ID and apply thresholding
	for (int i = 0; i < length; ++i) {
		int idKey = static_cast<int>(round(id[i]));
		if (maxValues.find(idKey) != maxValues.end() && maxValues[idKey] != 0) {
			// float normalizedThreshold = 0.75 * maxValues[idKey]; // for external cue guiding
			float normalizedThreshold = 0.6 * maxValues[idKey]; // for external cue guiding
			phiSum[i] = (phiSum[i] < normalizedThreshold) ? 0.0f : phiSum[i] / maxValues[idKey];
		}
	}
}

void NeuronGrowth::DetectTipsMulti(const vector<float>& phi_fine, const vector<int>& id, const int& numNeuron, vector<float>& phiSum, const int& rows, const int& cols, vector<int>& nonZeroIndices) {
	const int length = rows * cols;
	phiSum.assign(length, 0);  // Clear and resize with zero initialization
	nonZeroIndices.clear();  // Clear the indices vector to ensure it's fresh for this call

	const float threshold = 0.25; // Threshold for where to detect neurites
	const float intensityThreshold = 0.1; // Threshold for whether to use the value for detection
	const int offsetY = cols; // Offset to avoid index issues
	const float phiThreshold = 0; // Threshold for scaling

	// Lambda to replace CellBoundary and reduce condition checks - mostly for optimization
	auto cellBoundary = [](float phi, float threshold) -> float {
		return phi > threshold ? 1.0f : 0.0f;
	};

	// Using a map to store max values for each ID
	map<int, float> maxValues;

	// Iterate with precomputed values and reduced condition checks
	for (int i = (5 * cols); i < (length - 4 * cols); ++i) {
	if (cellBoundary(phi_fine[i], threshold) > 0 && round(id[i]) != 9) {
		for (int j = -4; j <= 4; ++j) {
			for (int k = -4; k <= 4; ++k) {
				int index = i + j * offsetY + k;
				if (index >= 0 && index < length && (round(id[index]) == round(id[i]) || round(id[index]) < 0)) {
					phiSum[i] += cellBoundary(phi_fine[index], intensityThreshold);
				}
			}
		}
		phiSum[i] = phiSum[i] > 0 ? cellBoundary(phi_fine[i], phiThreshold) / phiSum[i] : 0;
		if (isnan(phiSum[i])) phiSum[i] = 0;  // Handle NaN explicitly, though it should not occur now

		// Update max value for each ID
		int idKey = static_cast<int>(round(id[i]));
		if (maxValues.find(idKey) == maxValues.end() || maxValues[idKey] < phiSum[i]) {
			maxValues[idKey] = phiSum[i];
		}

		// Only store index if phiSum[i] is non-zero after processing
		if (phiSum[i] != 0) {
			nonZeroIndices.push_back(i);
		}
	}
	}

	// Normalize phiSum by the max value for each ID and apply thresholding
	for (int i = 0; i < length; ++i) {
		int idKey = static_cast<int>(round(id[i]));
		if (maxValues.find(idKey) != maxValues.end() && maxValues[idKey] != 0) {
			float normalizedThreshold = 0.6 * maxValues[idKey];
			//  : 0.825 * maxValues[idKey];
			phiSum[i] = (phiSum[i] < normalizedThreshold) ? 0.0f : phiSum[i] / maxValues[idKey];

			// Update non-zero indices if phiSum[i] changes to zero
			if (phiSum[i] == 0) {
				nonZeroIndices.erase(remove(nonZeroIndices.begin(), nonZeroIndices.end(), i), nonZeroIndices.end());
			}
		}
	}
}

// Function to perform Breadth-First Search (BFS) for clustering
float NeuronGrowth::bfs(const vector<float>& matrix, const int& rows, const int& cols, const int& row, const int& col,
	vector<bool>& visited, vector<pair<int, int>>& cluster) 
{
	
	static const int dx[] = {-1, 0, 1, 0, -1, -1, 1, 1}; // Including diagonal directions
	static const int dy[] = {0, -1, 0, 1, -1, 1, -1, 1}; // Including diagonal directions
	// static const int dx[] = {-1, 0, 1, 0}; // Including diagonal directions
	// static const int dy[] = {0, -1, 0, 1}; // Including diagonal directions
	
	float clusterValue = 0.0f;
	queue<pair<int, int>> q;
	q.push(make_pair(row, col));
	visited[row * cols + col] = true;
	cluster.push_back(make_pair(row, col));

	while (!q.empty()) {
		int r = q.front().first;
		int c = q.front().second;
		q.pop();

		for (int i = 0; i < 8; ++i) { // Iterating through all 8 directions
		// for (int i = 0; i < 4; ++i) { // Iterating through all 8 directions
			int newRow = r + dx[i];
			int newCol = c + dy[i];

			if (newRow >= 0 && newRow < rows && newCol >= 0 && newCol < cols &&
			matrix[newRow * cols + newCol] != 0 && !visited[newRow * cols + newCol]) {
				visited[newRow * cols + newCol] = true;
				q.push(make_pair(newRow, newCol));
				cluster.push_back(make_pair(newRow, newCol));
				clusterValue += matrix[newRow * cols + newCol]; // Increment cluster value with the current cell value
			}
		}
	}
	return clusterValue;
}

// Function to find connected clusters in the matrix
vector<vector<pair<int, int>>> NeuronGrowth::FindClusters(const vector<float>& matrix, const int& rows, const int& cols) 
{
    vector<bool> visited(rows * cols, false);
    vector<vector<pair<int, int>>> clusters;

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			if (matrix[i * cols + j] != 0 && !visited[i * cols + j]) {
				vector<pair<int, int>> cluster;
				float clusterValue = bfs(matrix, rows, cols, i, j, visited, cluster);

				if (!cluster.empty()) {
					clusters.push_back(cluster);
				}
			}
		}
	}
	return clusters;
}

// Function to find local maxima within connected clusters in the matrix
vector<float> NeuronGrowth::FindLocalMaximaInClusters(const vector<float>& matrix, const int& rows, const int& cols) 
{
	vector<vector<pair<int, int>>> clusters = FindClusters(matrix, rows, cols);
	vector<float> localMaxima(rows * cols, 0.0f);

	for (const auto& cluster : clusters) {
		float maxVal = -numeric_limits<float>::max();
		pair<int, int> maxPos = make_pair(-1, -1);

		for (const auto& pos : cluster) {
			int row = pos.first;
			int col = pos.second;

			if (matrix[row * cols + col] > maxVal) {
				maxVal = matrix[row * cols + col];
				maxPos = pos;
			}
		}

		if (maxPos.first != -1 && maxPos.second != -1) {
			localMaxima[maxPos.first * cols + maxPos.second] = 1.0f;
			
			// localMaxima[maxPos.first * cols + maxPos.second - cols + 1] = 1.0f;
			// localMaxima[maxPos.first * cols + maxPos.second - cols] = 1.0f;
			// localMaxima[maxPos.first * cols + maxPos.second - cols - 1] = 1.0f;
			// localMaxima[maxPos.first * cols + maxPos.second - 1] = 1.0f;
			// localMaxima[maxPos.first * cols + maxPos.second + 1] = 1.0f;
			// localMaxima[maxPos.first * cols + maxPos.second + cols + 1] = 1.0f;
			// localMaxima[maxPos.first * cols + maxPos.second + cols] = 1.0f;
			// localMaxima[maxPos.first * cols + maxPos.second + cols - 1] = 1.0f;
		}
	}
	return localMaxima;
}

vector<float> NeuronGrowth::KeepOneClusterWithMaxValue(const vector<float>& matrix, const int& rows, const int& cols) 
{
	vector<bool> visited(rows * cols, false);
	vector<float> localMaxima(rows * cols, 0.0f);
	float maxClusterValue = -numeric_limits<float>::max();
	vector<pair<int, int>> maxCluster;

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			if (matrix[i * cols + j] != 0 && !visited[i * cols + j]) {
				vector<pair<int, int>> cluster;
				float clusterValue = bfs(matrix, rows, cols, i, j, visited, cluster);

				if (clusterValue > maxClusterValue) {
					maxClusterValue = clusterValue;
					maxCluster = cluster;
				}
			}
		}
	}

	for (const auto& pos : maxCluster) {
		int row = pos.first;
		int col = pos.second;
		localMaxima[row * cols + col] = 1.0f;
	}

	// // Calculate centroid of local maxima
	// float sumRow = 0, sumCol = 0;
	// for (const auto& maxPos : maxCluster) {
	// 	sumRow += maxPos.first;
	// 	sumCol += maxPos.second;
	// }

	// int centroidRow = static_cast<int>(round(sumRow / maxCluster.size()));
	// int centroidCol = static_cast<int>(round(sumCol / maxCluster.size()));

	// localMaxima[centroidRow * cols + centroidCol] = 1.0f;


	return localMaxima;
}

// Function to find centroids of connected clusters in the matrix
vector<float> NeuronGrowth::FindCentroidsInClusters(const vector<float>& matrix, const int& rows, const int& cols) 
{
	vector<vector<pair<int, int>>> clusters = FindClusters(matrix, rows, cols);
	vector<float> centroids(rows * cols, 0.0f);

	for (const auto& cluster : clusters) {
		if (cluster.empty()) continue;

		float sumRow = 0.0f, sumCol = 0.0f;
		for (const auto& pos : cluster) {
			sumRow += pos.first;
			sumCol += pos.second;
		}
		int centroidRow = static_cast<int>(round(sumRow / cluster.size()));
		int centroidCol = static_cast<int>(round(sumCol / cluster.size()));

		// Ensure the centroid position is within bounds
		centroidRow = max(0, min(centroidRow, rows - 1));
		centroidCol = max(0, min(centroidCol, cols - 1));

		centroids[centroidRow * cols + centroidCol] = 1.0f;

		// centroids[centroidRow.first * cols + centroidRow.second - cols + 1] = 1.0f;
		// centroids[centroidRow.first * cols + centroidRow.second - cols] = 1.0f;
		// centroids[centroidRow.first * cols + centroidRow.second - cols - 1] = 1.0f;
		// centroids[centroidRow.first * cols + centroidRow.second - 1] = 1.0f;
		// centroids[centroidRow.first * cols + centroidRow.second + 1] = 1.0f;
		// centroids[centroidRow.first * cols + centroidRow.second + cols + 1] = 1.0f;
		// centroids[centroidRow.first * cols + centroidRow.second + cols] = 1.0f;
		// centroids[centroidRow.first * cols + centroidRow.second + cols - 1] = 1.0f;
	}
	return centroids;
}

bool NeuronGrowth::IsLocalMaximum(const vector<float>& matrix, const int& rows, const int& cols, const int& x, const int& y) 
{
	float value = matrix[x * cols + y];
	for (int dx = -1; dx <= 1; ++dx) {
		for (int dy = -1; dy <= 1; ++dy) {
			int nx = x + dx, ny = y + dy;
			if (nx >= 0 && nx < rows && ny >= 0 && ny < cols && (dx != 0 || dy != 0)) {
				if (matrix[nx * cols + ny] > value) {
					return false;
				}
			}
		}
	}
	return true;
}

vector<float> NeuronGrowth::FindCentroidsOfLocalMaximaClusters(const vector<float>& matrix, const int& rows, const int& cols) 
{
	vector<vector<pair<int, int>>> clusters = FindClusters(matrix, rows, cols);
	vector<float> centroids(rows * cols, 0.0f);

	for (const auto& cluster : clusters) {
		vector<pair<int, int>> localMaxima;

		// Find all local maxima in the cluster
		for (const auto& pos : cluster) {
			if (IsLocalMaximum(matrix, rows, cols, pos.first, pos.second)) {
				localMaxima.push_back(pos);
			}
		}

		// Calculate centroid of local maxima
		float sumRow = 0, sumCol = 0;
		for (const auto& maxPos : localMaxima) {
			sumRow += maxPos.first;
			sumCol += maxPos.second;
		}

		int centroidRow = static_cast<int>(round(sumRow / localMaxima.size()));
		int centroidCol = static_cast<int>(round(sumCol / localMaxima.size()));

		centroids[centroidRow * cols + centroidCol] = 5.0f;

		// centroids[centroidRow * cols + centroidCol - cols + 1] = 1.0f;
		// centroids[centroidRow * cols + centroidCol- cols] = 1.0f;
		// centroids[centroidRow * cols + centroidCol - cols - 1] = 1.0f;
		// centroids[centroidRow * cols + centroidCol - 1] = 1.0f;
		// centroids[centroidRow * cols + centroidCol + 1] = 1.0f;
		// centroids[centroidRow * cols + centroidCol + cols + 1] = 1.0f;
		// centroids[centroidRow * cols + centroidCol + cols] = 1.0f;
		// centroids[centroidRow * cols + centroidCol + cols - 1] = 1.0f;
	}

	return centroids;
}

vector<float> NeuronGrowth::FindCentroidsOfLocalMaximaClusters(const vector<float>& matrix, const int& rows, const int& cols, vector<int>& centroidIndices)
{
	vector<vector<pair<int, int>>> clusters = FindClusters(matrix, rows, cols);
	vector<float> centroids(rows * cols, 0.0f);
	centroidIndices.clear();  // Clearing existing indices to ensure fresh data on each call

	for (const auto& cluster : clusters) {
		vector<pair<int, int>> localMaxima;

		// Find all local maxima in the cluster
		for (const auto& pos : cluster) {
			if (IsLocalMaximum(matrix, rows, cols, pos.first, pos.second)) {
				localMaxima.push_back(pos);
			}
		}

		if (!localMaxima.empty()) {  // Ensure there are local maxima to calculate centroid
			// Calculate centroid of local maxima
			float sumRow = 0, sumCol = 0;
			for (const auto& maxPos : localMaxima) {
				sumRow += maxPos.first;
				sumCol += maxPos.second;
			}

			int centroidRow = static_cast<int>(round(sumRow / localMaxima.size()));
			int centroidCol = static_cast<int>(round(sumCol / localMaxima.size()));
			int centroidIndex = centroidRow * cols + centroidCol;

			centroids[centroidIndex] = 5.0f;

			// centroids[centroidRow * cols + centroidCol - cols + 1] = 1.0f;
			// centroids[centroidRow * cols + centroidCol- cols] = 1.0f;
			// centroids[centroidRow * cols + centroidCol - cols - 1] = 1.0f;
			// centroids[centroidRow * cols + centroidCol - 1] = 1.0f;
			// centroids[centroidRow * cols + centroidCol + 1] = 1.0f;
			// centroids[centroidRow * cols + centroidCol + cols + 1] = 1.0f;
			// centroids[centroidRow * cols + centroidCol + cols] = 1.0f;
			// centroids[centroidRow * cols + centroidCol + cols - 1] = 1.0f;

			centroidIndices.push_back(centroidIndex);  // Store the index of the centroid
		}
	}

	return centroids;
}

vector<float> NeuronGrowth::ComputeMaxFilter(const vector<float>& geodist, int numRows, int numCols, int windowRadius) {
	vector<float> maxGeodist(geodist.size(), numeric_limits<float>::lowest());

	for (int row = 0; row < numRows; ++row) {
		for (int col = 0; col < numCols; ++col) {
			int startIndex = row * numCols + col;
			int startRow = max(0, row - windowRadius);
			int endRow = min(numRows - 1, row + windowRadius);
			int startCol = max(0, col - windowRadius);
			int endCol = min(numCols - 1, col + windowRadius);

			for (int r = startRow; r <= endRow; ++r) {
				for (int c = startCol; c <= endCol; ++c) {
					int index = r * numCols + c;
					maxGeodist[startIndex] = max(maxGeodist[startIndex], geodist[index]);
				}
			}
		}
	}
	return maxGeodist;
}

// vector<vector<int>> NeuronGrowth::ConvertTo2DIntVector_PushBoundary(const vector<float>& input, const int& NX, const int& NY) 
// {
// 	vector<vector<int>> output;
// 	output.reserve(NX + 1); // Reserve space for rows

// 	int k = 0;
// 	for (int i = 0; i < NX + 1; i++) {
// 		vector<int> row;
// 		row.reserve(NY + 1); // Reserve space for columns

// 		for (int j = 0; j < NY + 1; j++, k++) { // Increment k in the loop header
// 			// Assuming CellBoundary returns 1 if condition is met, 0 otherwise
// 			row.push_back(CellBoundary(abs(input[k]), 0.25) * 9);
// 		}
// 		output.push_back(move(row)); // Use move to avoid copying the row
// 	}

// 	return output;
// }
vector<int> NeuronGrowth::ConvertTo1DIntVector_PushBoundary(const vector<float>& input, const int& NX, const int& NY) {
	vector<int> output((NX + 1) * (NY + 1));  // Create a 1D vector with enough space for all elements

	int k = 0;
	for (int i = 0; i < NX + 1; i++) {
		for (int j = 0; j < NY + 1; j++, k++) {
			output[i * (NY + 1) + j] = CellBoundary(abs(input[k]), 0.25) * 9;
		}
	}

	return output;
}


vector<vector<int>> NeuronGrowth::ConvertTo2DIntVector(const vector<float>& input, int NX, int NY) 
{
	vector<vector<int>> output;
	output.reserve(NX + 1); // Preallocate memory for rows

	size_t k = 0; // Use size_t for consistency with standard library sizes
	for (int i = 0; i < NX + 1; i++) {
		vector<int> row;
		row.reserve(NY + 1); // Preallocate memory for columns

		for (int j = 0; j < NY + 1; j++, k++) {
			row.push_back(static_cast<int>(input[k])); // Explicit type conversion shown
		}
		output.push_back(row);
	}
	return output;
}

vector<vector<float>> NeuronGrowth::ConvertTo2DFloatVector(const vector<float>& input, int NX, int NY) 
{
	vector<vector<float>> output;
	output.reserve(NX); // Preallocate memory for rows

	size_t k = 0; // Use size_t for consistency with standard library sizes
	for (int i = 0; i < NX; i++) {
		vector<float> row;
		row.reserve(NY); // Preallocate memory for columns

		for (int j = 0; j < NY; j++, k++) {
			row.push_back(input[k]);
		}
		output.push_back(row);
	}
	return output;
}

// // Function to perform flood fill
// void NeuronGrowth::FloodFill(vector<vector<int>>& image, int x, int y, int newColor, int originalColor, const vector<vector<int>>& prev_id) 
// {
// 	// Define the directions: up, down, left, right
// 	int dx[] = {0, 0, -1, 1};
// 	int dy[] = {-1, 1, 0, 0};

// 	if (x < 0 || x >= image.size() || y < 0 || y >= image[0].size() || image[x][y] != originalColor || image[x][y] == newColor) {
// 		return;
// 	}

// 	image[x][y] = newColor;

// 	// Apply flood fill in all four directions
// 	for (int i = 0; i < 4; ++i) {
// 		if (prev_id[x + dx[i]][y + dy[i]] == (newColor) || prev_id[x + dx[i]][y + dy[i]] == 0 || prev_id[x + dx[i]][y + dy[i]] == 9) {
// 			FloodFill(image, x + dx[i], y + dy[i], newColor, originalColor, prev_id);
// 		}
// 	}
// }
void NeuronGrowth::FloodFill(vector<int>& image, int x, int y, int newColor, int originalColor, const vector<int>& prev_id, int width) {
	int dx[] = {0, 0, -1, 1};
	int dy[] = {-1, 1, 0, 0};
	int index = x * width + y;

	if (x < 0 || x >= image.size() / width || y < 0 || y >= width || image[index] != originalColor || image[index] == newColor) {
		return;
	}

	image[index] = newColor;

	for (int i = 0; i < 4; ++i) {
		int nx = x + dx[i];
		int ny = y + dy[i];
		int nIndex = nx * width + ny;
		if (nx >= 0 && ny >= 0 && nx < image.size() / width && ny < width &&
			(prev_id[nIndex] == newColor || prev_id[nIndex] == 0 || prev_id[nIndex] == 9)) {
			FloodFill(image, nx, ny, newColor, originalColor, prev_id, width);
		}
	}
}


// // Function to perform flood fill
// void NeuronGrowth::IdentifyNeurons(vector<float>& phi_in, vector<vector<int>>& neurons, const vector<vector<int>>& prev_id, vector<array<float, 2>> seed, const int& NX, const int& NY, const int& originX, const int& originY) 
// {
// 	neurons = ConvertTo2DIntVector_PushBoundary(phi_in, NX, NY);
// 	for (size_t i = 0; i < seed.size(); i++) {
// 		int startX = seed[i][0] - originX;
// 		int startY = seed[i][1] - originY;
// 		int newColor = i+1;
// 		int originalColor = neurons[startX][startY];
// 		FloodFill(neurons, startX, startY, newColor, originalColor, prev_id);
// 	}
// }

void NeuronGrowth::IdentifyNeurons(vector<float>& phi_in, vector<int>& neurons, const vector<int>& prev_id, vector<array<float, 2>> seed, const int& NX, const int& NY, const int& originX, const int& originY) {
	neurons = ConvertTo1DIntVector_PushBoundary(phi_in, NX, NY);
	int width = NY + 1; // Define the width for correct indexing in the 1D vector

	for (size_t i = 0; i < seed.size(); i++) {
		int startX = seed[i][0] - originX;
		int startY = seed[i][1] - originY;
		int newColor = i + 1;
		int originalColor = neurons[startX * width + startY];
		FloodFill(neurons, startX, startY, newColor, originalColor, prev_id, width);
	}
}


// Helper function to check if a cell is within the grid and not already processed as a connection
bool isValid_test(int x, int y, int rows, int cols, const vector<vector<int>>& grid) {
    return x >= 0 && x < rows && y >= 0 && y < cols && grid[x][y] != -1;
}

// // Function to modify the grid where clusters are connected
// void NeuronGrowth::DetectConnections(vector<vector<int>>& grid) {

// 	const vector<pair<int, int>> directions = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}};

// 	int rows = grid.size(), cols = grid[0].size();
// 	vector<vector<bool>> visited(rows, vector<bool>(cols, false));

// 	for (int i = 0; i < rows; ++i) {
// 		for (int j = 0; j < cols; ++j) {
// 			if (grid[i][j] > 0 && !visited[i][j]) {
// 				int clusterId = grid[i][j];
// 				queue<pair<int, int>> q;
// 				q.push({i, j});
// 				visited[i][j] = true;

// 				while (!q.empty()) {
// 					auto [x, y] = q.front();
// 					q.pop();

// 					// Explore the four possible directions
// 					for (auto& dir : directions) {
// 						int nx = x + dir.first;
// 						int ny = y + dir.second;

// 						if (isValid_test(nx, ny, rows, cols, grid)) {
// 							if (grid[nx][ny] == clusterId) {
// 								if (!visited[nx][ny]) {
// 									visited[nx][ny] = true;
// 									q.push({nx, ny});
// 								}
// 							} else if (grid[nx][ny] > 0 && grid[nx][ny] != clusterId) {
// 								// Mark the current and adjacent cells as connection points
// 								grid[nx][ny] = -1;
// 								grid[x][y] = -1;
// 								visited[nx][ny] = true;  // Prevent reprocessing the marked connection
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
// }

void NeuronGrowth::DetectConnections(vector<int>& grid, int rows, int cols) {
	const vector<pair<int, int>> directions = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}};
	vector<bool> visited(rows * cols, false);

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			int index = i * cols + j; // Calculate index for flat vector
			if (grid[index] > 0 && !visited[index]) {
				int clusterId = grid[index];
				queue<pair<int, int>> q;
				q.push({i, j});
				visited[index] = true;

				while (!q.empty()) {
					auto [x, y] = q.front();
					q.pop();
					int currentIndex = x * cols + y;

					// Explore the four possible directions
					for (auto& dir : directions) {
						int nx = x + dir.first;
						int ny = y + dir.second;
						int nextIndex = nx * cols + ny;

						if (nx >= 0 && ny >= 0 && nx < rows && ny < cols) { // isValid simplified
							if (grid[nextIndex] == clusterId) {
								if (!visited[nextIndex]) {
									visited[nextIndex] = true;
									q.push({nx, ny});
								}
							} else if (grid[nextIndex] > 0 && grid[nextIndex] != clusterId) {
								// Mark the current and adjacent cells as connection points
								grid[nextIndex] = -1;
								grid[currentIndex] = -1;
								visited[nextIndex] = true;  // Prevent reprocessing the marked connection
							}
						}
					}
				}
			}
		}
	}
}

// Check if a point is valid in the grid
bool NeuronGrowth::isValid(const int& x, const int& y, const int& rows, const int& cols)
{
    return (x >= 0 && x < rows && y >= 0 && y < cols);
}

vector<vector<vector<float>>> NeuronGrowth::CalculateQuasiEuclideanDistanceFromPoint(const vector<vector<int>>& neurons, vector<array<float, 2>>& seed, const int& originX, const int& originY) 
{
	int rows = neurons.size();
	if (rows == 0) return {};

	int cols = neurons[0].size();

	// Initialize distances matrix with zeros
	vector<vector<vector<float>>> distances(seed.size(), vector<vector<float>>(rows, vector<float>(cols, -1.0f)));

	// Define the direction arrays
	static const int dx[] = {-1, 1, 0, 0};
	static const int dy[] = {0, 0, -1, 1};

	for (size_t i = 0; i < seed.size(); i++) {
		int startX = seed[i][0] - originX;
		int startY = seed[i][1] - originY;

		// Mark the starting point as visited
		vector<vector<bool>> visited(rows, vector<bool>(cols, false));
		visited[startX][startY] = true;

		// Initialize the BFS queue with the starting point
		queue<pair<int, int>> q;
		q.push({startX, startY});
		distances[i][startX][startY] = 0.0f;

		// BFS loop
		while (!q.empty()) {
			pair<int, int> current = q.front();
			q.pop();

			int x = current.first;
			int y = current.second;

			// Check neighboring cells
			for (int j = 0; j < 4; ++j) {
				int newX = x + dx[j];
				int newY = y + dy[j];

				if (isValid(newX, newY, rows, cols) && neurons[newX][newY] == neurons[startX][startY] && !visited[newX][newY]) {
					visited[newX][newY] = true;

					// Calculate Quasi-Euclidean distance
					float dx = abs(newX - startX);
					float dy = abs(newY - startY);
					float distance = distances[i][x][y] + sqrt(dx * dx + dy * dy);

					// Update distance if it's smaller than the previous or uninitialized
					if (distances[i][newX][newY] == -1.0f || distance < distances[i][newX][newY])
						distances[i][newX][newY] = distance;

					// Add the neighboring cell to the BFS queue
					q.push({newX, newY});
				}
			}
		}
	}

	return distances;
}

vector<vector<vector<float>>> NeuronGrowth::CalculateGeodesicDistanceFromPoint(const vector<vector<int>>& neurons, vector<array<float, 2>>& seed, const int& originX, const int& originY)
{
	int rows = neurons.size();
	if (rows == 0) return {};

	int cols = neurons[0].size();

	// Initialize distances matrix with zeros
	vector<vector<vector<float>>> distances(seed.size(), vector<vector<float>>(rows, vector<float>(cols, -1.0f)));

	// Define the direction arrays
	static const int dx[] = {-1, 1, 0, 0};
	static const int dy[] = {0, 0, -1, 1};

	for (size_t i = 0; i < seed.size(); i++) {
		int startX = seed[i][0] - originX;
		int startY = seed[i][1] - originY;

		// Mark the starting point as visited
		vector<vector<bool>> visited(rows, vector<bool>(cols, false));
		visited[startX][startY] = true;

		// Initialize the BFS queue with the starting point
		queue<pair<int, int>> q;
		q.push({startX, startY});
		distances[i][startX][startY] = 0.0f;

		// BFS loop
		while (!q.empty()) {
			pair<int, int> current = q.front();
			q.pop();

			int x = current.first;
			int y = current.second;

			// Check neighboring cells
			for (int j = 0; j < 4; ++j) {
				int newX = x + dx[j];
				int newY = y + dy[j];

				if (isValid(newX, newY, rows, cols) && neurons[newX][newY] == neurons[startX][startY] && !visited[newX][newY]) {
					visited[newX][newY] = true;

					// Update geodesic distance
					distances[i][newX][newY] = distances[i][x][y] + 1.0f;

					// Add the neighboring cell to the BFS queue
					q.push({newX, newY});
				}
			}
		}
	}

	return distances;
}

vector<vector<float>> NeuronGrowth::ExploreGridAndCalculateDistances(vector<int>& grid, const vector<array<float, 2>>& seed, const int& originX, const int& originY, const int& rows, const int& cols) {
	// Direction vectors for moving in the grid (up, down, left, right)
	const vector<pair<int, int>> directions = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}};
	vector<vector<float>> distances(seed.size(), vector<float>(rows * cols, -1.0f)); // Each 2D matrix is now a single vector<float>
	vector<vector<bool>> visited(rows, vector<bool>(cols, false));

	for (size_t seedIndex = 0; seedIndex < seed.size(); seedIndex++) {
		int startX = seed[seedIndex][0] - originX;
		int startY = seed[seedIndex][1] - originY;

		if (!isValid(startX, startY, rows, cols)) continue;

		queue<pair<int, int>> q;
		q.push({startX, startY});
		visited[startX][startY] = true;
		distances[seedIndex][startX * cols + startY] = 0.0f; // Access using flat index

		while (!q.empty()) {
			auto [x, y] = q.front();
			q.pop();

			for (auto& dir : directions) {
				int nx = x + dir.first;
				int ny = y + dir.second;

				if (isValid(nx, ny, rows, cols)) {
					int currIndex = x * cols + y;
					int nextIndex = nx * cols + ny;

					if (!visited[nx][ny] && grid[nextIndex] == grid[currIndex]) {
						visited[nx][ny] = true;

						// Calculate Quasi-Euclidean distance
						float dx = abs(nx - startX);
						float dy = abs(ny - startY);
						float distance = distances[seedIndex][currIndex] + sqrt(dx * dx + dy * dy);

						if (distances[seedIndex][nextIndex] == -1.0f || distance < distances[seedIndex][nextIndex])
							distances[seedIndex][nextIndex] = distance;

						q.push({nx, ny});
					} else if (grid[nextIndex] > 0 && grid[nextIndex] != grid[currIndex]) {
					grid[nextIndex] = -1; // Mark the connection point
					grid[currIndex] = -1; // Also mark the current cell as a connection point
					}
				}
			}
		}
	}

	return distances;
}


// // A combined function to explore the grid for connections between clusters (neurons) and calculate
// // Quasi-Euclidean distances from seed points within the same cluster.
// vector<vector<vector<float>>> NeuronGrowth::ExploreGridAndCalculateDistances(vector<vector<int>>& grid, const vector<array<float, 2>>& seed, const int& originX, const int& originY) 
// {
// 	// Direction vectors for moving in the grid (up, down, left, right)
// 	const vector<pair<int, int>> directions = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}};
// 	int rows = grid.size(), cols = grid[0].size();
// 	vector<vector<vector<float>>> distances(seed.size(), vector<vector<float>>(rows, vector<float>(cols, -1.0f)));
// 	vector<vector<bool>> visited(rows, vector<bool>(cols, false));

// 	for (size_t seedIndex = 0; seedIndex < seed.size(); seedIndex++) {
// 		int startX = seed[seedIndex][0] - originX;
// 		int startY = seed[seedIndex][1] - originY;

// 		if (!isValid(startX, startY, rows, cols)) continue;

// 		queue<pair<int, int>> q;
// 		q.push({startX, startY});
// 		visited[startX][startY] = true;
// 		distances[seedIndex][startX][startY] = 0.0f;

// 		while (!q.empty()) {
// 			auto [x, y] = q.front();
// 			q.pop();

// 			for (auto& dir : directions) {
// 				int nx = x + dir.first;
// 				int ny = y + dir.second;

// 				if (isValid(nx, ny, rows, cols)) {
// 					if (!visited[nx][ny] && grid[nx][ny] == grid[startX][startY]) {
// 						visited[nx][ny] = true;

// 						// Calculate Quasi-Euclidean distance
// 						float dx = abs(nx - startX);
// 						float dy = abs(ny - startY);
// 						float distance = distances[seedIndex][x][y] + sqrt(dx * dx + dy * dy);

// 						if (distances[seedIndex][nx][ny] == -1.0f || distance < distances[seedIndex][nx][ny])
// 							distances[seedIndex][nx][ny] = distance;

// 						q.push({nx, ny});
// 					} else if (grid[nx][ny] > 0 && grid[nx][ny] != grid[startX][startY]) {
// 						grid[nx][ny] = -1; // Mark the connection point
// 						grid[x][y] = -1;   // Also mark the current cell as a connection point
// 					}
// 				}
// 			}
// 		}
// 	}

// 	return distances;
// }

void NeuronGrowth::MaskSoma(vector<float>& phi, int width, int height, vector<int> somaCenter, float somaRadius) {
	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			if (sqrt(pow(x - somaCenter[0], 2) + pow(y - somaCenter[1], 2)) <= somaRadius) {
				phi[y * width + x] = 0; // Mask out soma
			}
		}
	}
}

vector<vector<vector<int>>> NeuronGrowth::TraceNeurites(vector<float>& phi, const vector<float>& geodesicDistance, int width, int height, int generation) {
	vector<vector<vector<int>>> traces;
	vector<bool> visited(width * height, false);

	while (true) {
		float maxDist = -1;
		int maxIdx = -1;
		// Find maximum unvisited geodist and untraced point
		for (int i = 0; i < geodesicDistance.size(); ++i) {
			if (phi[i] == 1 && !visited[i] && geodesicDistance[i] > maxDist) {
				maxDist = geodesicDistance[i];
				maxIdx = i;
			}
		}

		if (maxIdx == -1) break; // No more neurites to trace

		vector<vector<int>> trace;
		int idx = maxIdx;
		while (idx != -1 && !visited[idx]) {
			int x = idx % width;
			int y = idx / width;
			trace.push_back({x, y});
			visited[idx] = true;

			// Find next point with lower geodist
			int nextIdx = -1;
			float minDist = geodesicDistance[idx];
			for (int dy = -1; dy <= 1; ++dy) {
				for (int dx = -1; dx <= 1; ++dx) {
					int nx = x + dx;
					int ny = y + dy;
					int nIdx = ny * width + nx;
					if (isValid(nx, ny, width, height) && phi[nIdx] == 1 && geodesicDistance[nIdx] < minDist && !visited[nIdx]) {
						minDist = geodesicDistance[nIdx];
						nextIdx = nIdx;
					}
				}
			}
			idx = nextIdx;
		}
		traces.push_back(trace);
	}

	// Remove traced neurites from phi to not retrace them
	for (const auto& trace : traces) {
		for (const auto& p : trace) {
			phi[p[1] * width + p[0]] = 0;
		}
	}

	cout << "Generation " << generation << ": Traced " << traces.size() << " neurites." << endl;
	return traces;
}

// void NeuronGrowth::AdjustNearestTip(vector<float>& localMaximaMatrix, int width, int height, const vector<vector<int>>& cues, double adjustmentFactor) {
// 	vector<float> tempMatrix(localMaximaMatrix.size(), 0); // Start with current tip states

// 	// std::cout << cues.size() << std::endl;
// 	for (const auto& cue : cues) {
// 		double minDist = numeric_limits<double>::max();
// 		int nearestTipIndex = -1;

// 		// Traverse through the matrix to find the nearest tip to the current cue
// 		for (int y = 0; y < height; ++y) {
// 			for (int x = 0; x < width; ++x) {
// 				int idx = y * width + x;
// 				if (localMaximaMatrix[idx] > 0) { // Check if this is a tip
// 					double dist = pow(x - cue[0], 2) + pow(y - cue[1], 2);
// 					if (dist < minDist) {
// 						minDist = dist;
// 						nearestTipIndex = idx;
// 					}
// 				}
// 			}
// 		}

// 		// Adjust the found tip towards the cue (if any)
// 		if (nearestTipIndex != -1) {
// 			int tipX = nearestTipIndex % width;
// 			int tipY = nearestTipIndex / width;
// 			// Calculate new positions, adjusted towards the cue
// 			int newX = max(0, min(width - 1, tipX + static_cast<int>(adjustmentFactor * (cue[0] - tipX))));
// 			int newY = max(0, min(height - 1, tipY + static_cast<int>(adjustmentFactor * (cue[1] - tipY))));
// 			// int newX = tipX + static_cast<int>(adjustmentFactor * (cue[0] - tipX));
// 			// int newY = tipY + static_cast<int>(adjustmentFactor * (cue[1] - tipY));
// 			int newIdx = newY * width + newX;
// 			// tempMatrix[nearestTipIndex] = 0; // Remove the old position
// 			tempMatrix[newIdx] = 5; // Mark the new position as a tip
// 		}
// 	}

// 	localMaximaMatrix = tempMatrix; // Update the localMaximaMatrix with the adjusted tip positions
// }
void NeuronGrowth::AdjustNearestTip(vector<float>& localMaximaMatrix, int width, int height, const vector<vector<int>>& cues) {

	vector<float> tempMatrix(localMaximaMatrix.size(), 0); // Start with current tip states

	for (const auto& cue : cues) {
		double minDist = numeric_limits<double>::max();
		int nearestTipIndex = -1;

		// Traverse through the matrix to find the nearest tip to the current cue
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				int idx = y * width + x;
				if (localMaximaMatrix[idx] > 0) { // Check if this is a tip
					double dist = pow(x - cue[0], 2) + pow(y - cue[1], 2);
					if (dist < minDist) {
						minDist = dist;
						nearestTipIndex = idx;
					}
				}
			}
		}

		// Adjust the found tip towards the cue (if any)
		if (nearestTipIndex != -1) {
			int tipX = nearestTipIndex % width;
			int tipY = nearestTipIndex / width;
			
			// Calculate direction vector
			double dist = sqrt(pow(tipX - cue[0], 2) + pow(tipY - cue[1], 2));
			int moveX = dist > 0 ? static_cast<int>((cue[0] - tipX) / dist * 1) : 0; // Move approximately n pixels
			int moveY = dist > 0 ? static_cast<int>((cue[1] - tipY) / dist * 1) : 0;

			// Calculate new positions, adjusted towards the cue
			int newX = max(0, min(width - 1, tipX + moveX));
			int newY = max(0, min(height - 1, tipY + moveY));

			int newIdx = newY * width + newX;
			tempMatrix[newIdx] = 5; // Mark the new position as a tip

			// Optional debugging output
			// cout << "Tip moved from (" << tipX << ", " << tipY << ") to (" << newX << ", " << newY << ")" << endl;
		}
	}

	localMaximaMatrix = tempMatrix; // Update the localMaximaMatrix with the adjusted tip positions
}

vector<float> NeuronGrowth::PickNearestTip(vector<float>& tip, int width, int height, const vector<vector<int>>& cues, const vector<int>& centroidsIndex) {
	vector<float> tempMatrix(tip.size(), 0); // Initialize output matrix with zeros

	for (const auto& cue : cues) {
		if (cue.empty()) {
			// cout << "Empty cue encountered, skipping." << endl;
			continue;
		}

		int cue_x = cue[0] + rand() % 7 - 3;
		int cue_y = cue[1] + rand() % 7 - 3;
		
		// std::cout << cue[1] << " " << width << " " << cue[0] << std::endl;;
		float minDist = numeric_limits<float>::max();
		int nearestTipIndex = -1;

		// // Iterate only over indices that are non-zero tips
		// for (const int& idx : centroidsIndex) {
		// 	int x = idx % height;
		// 	int y = idx / height;

		// 	if (tip[idx] != 0) { // Check if this is a tip pixel
		// 		// Calculate Manhattan distance to the cue
		// 		float dist = abs(x - cue[0]) + abs(y - cue[1]);
		// 		if (dist < minDist) {
		// 			minDist = dist;
		// 			nearestTipIndex = idx;
		// 		}
		// 	}
		// }

		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				int idx = x * height + y;
				if (tip[idx] != 0) { // Check if this is a tip pixel
					float dist = abs(x - cue_x) + abs(y - cue_y); // Manhattan distance
					if (dist < minDist) {
						minDist = dist;
						nearestTipIndex = idx;
					}
				}
			}
		}

		// Set the nearest tip pixel to the cue in the output matrix
		if (nearestTipIndex != -1) {
			tempMatrix[nearestTipIndex] = 5;
		}

		// Find the index of the cue directly since its position is known
		int cueId = cue[0] * height + cue[1];
		// Also mark the cue position - for debugging
		if (cueId >= 0 && cueId < tempMatrix.size()) {
			tempMatrix[cueId] = -5;
		} else {
			std::cout << "cueId: " << cueId << " " << cue_x << " " << cue_y << std::endl;
 		}
	}

	return tempMatrix;
}

vector<float> NeuronGrowth::PickNearestTip(vector<float>& tip, const vector<int>& id, int width, int height, const vector<vector<vector<int>>>& Allcues, const vector<int>& centroidsIndex) {
	vector<float> tempMatrix(tip.size(), 0); // Initialize output matrix with zeros

	for (int i = 0; i < Allcues.size(); i++) {
	// std::cout << i << " " << Allcues.size() << std::endl;
		auto cues = Allcues[i];
		for (const auto& cue : cues) {
			if (cue.empty()) {
				// cout << "Empty cue encountered, skipping." << endl;
				continue;
			}
			int cue_x = cue[0] + rand() % 7 - 3;
			int cue_y = cue[1] + rand() % 7 - 3;
			// std::cout << cue[1] << " " << width << " " << cue[0] << std::endl;;
			float minDist = numeric_limits<float>::max();
			int nearestTipIndex = -1;

			// // Iterate only over indices that are non-zero tips
			// for (const int& idx : centroidsIndex) {
			// 	int x = idx % height;
			// 	int y = idx / height;

			// 	// Calculate Manhattan distance to the cue
			// 	float dist = abs(x - cue[0]) + abs(y - cue[1]);
			// 	if (dist < minDist) {
			// 		minDist = dist;
			// 		nearestTipIndex = idx;
			// 	}
			// }

			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					int idx = x * height + y;
					if (tip[idx] != 0 && id[idx] == i+1) { // Check if this is a tip pixel
						float dist = abs(x - cue_x) + abs(y - cue_y); // Manhattan distance
						if (dist < minDist) {
							minDist = dist;
							nearestTipIndex = idx;
						}
					}
				}
			}

			// Set the nearest tip pixel to the cue in the output matrix
			if (nearestTipIndex != -1) {
				tempMatrix[nearestTipIndex] = 5;
			}

			// Find the index of the cue directly since its position is known
			int cueId = cue_x * height + cue_y;
			// Also mark the cue position - for debugging
			if (cueId >= 0 && cueId < tempMatrix.size()) {
				tempMatrix[cueId] = -5;
			} else {
				std::cout << "cueId: " << cueId << " " << cue_x << " " << cue[1] << std::endl;
			}
		}
	}

	return tempMatrix;
}

// Function to calculate geodesic distance from a point
void NeuronGrowth::SaveNGvars(vector<vector<float>>& NGvars, int NX, int NY, string fn) 
{
	// writing initialized variables for debugging purposes
	bool visualization = true;
	PrintVec2TXT(NGvars[0], fn + "/phi_" + to_string(n) + ".txt", visualization);
	PrintVec2TXT(NGvars[1], fn + "/syn_" + to_string(n) + ".txt", visualization);
	PrintVec2TXT(NGvars[2], fn + "/tub_" + to_string(n) + ".txt", visualization);
	PrintVec2TXT(NGvars[3], fn + "/theta_" + to_string(n) + ".txt", visualization);
	PrintVec2TXT(NGvars[4], fn + "/phi_0_" + to_string(n) + ".txt", visualization);
	PrintVec2TXT(NGvars[5], fn + "/tub_0_" + to_string(n) + ".txt", visualization);
}

void NeuronGrowth::PrintOutNeurons(const vector<int>& neurons, int NX_fine, int NY_fine) 
{   
	int dwnRatio = 2; // Assuming line length is 78 characters, this could be adjusted for different sizes

	for (size_t i = 0; i < 78; i++) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "-");
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");

	// Adjust for the actual size of neurons in a flat structure
	for (size_t j = 0; j < NY_fine; j += 2*dwnRatio) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "| ");
		for (size_t i = 0; i < NX_fine; i += dwnRatio) {
			// Calculate the flat index from 2D indices
			int index = i * NY_fine + j;
			if (index < neurons.size()) {
				if (neurons[index] == 0) {
					ierr = PetscPrintf(PETSC_COMM_WORLD, " ");
				} else {
					ierr = PetscPrintf(PETSC_COMM_WORLD, "%d", neurons[index]);
				}
			}
		}
		ierr = PetscPrintf(PETSC_COMM_WORLD, "|\n");
	}

	for (size_t i = 0; i < 78; i++) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "-");
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");
}

void NeuronGrowth::PrintOutNeurons(const vector<int>& neurons, const vector<float>& localMax, int NX_fine, int NY_fine) 
{   
	PetscErrorCode ierr;
	int poolSize = 2; // Adjust this as needed

	// First, ensure that dimensions are divisible by poolSize or adjust them appropriately
	int adjustedNX = (NX_fine % poolSize == 0) ? NX_fine : (NX_fine / poolSize + 1) * poolSize;
	int adjustedNY = (NY_fine % poolSize == 0) ? NY_fine : (NY_fine / poolSize + 1) * poolSize;

	std::vector<float> combined(adjustedNX * adjustedNY, 0); // Create a larger buffer if needed
	for (int i = 0; i < NX_fine; i++) {
		for (int j = 0; j < NY_fine; j++) {
			combined[i * adjustedNY + j] = std::abs(neurons[i * NY_fine + j] + localMax[i * NY_fine + j]);
		}
	}

	// Perform max pooling
	std::vector<float> pooled;
	for (int y = 0; y < adjustedNY; y += poolSize) {
		for (int x = 0; x < adjustedNX; x += poolSize) {
			float maxVal = 0;
			for (int dy = 0; dy < poolSize && (y + dy) < adjustedNY; dy++) {
				for (int dx = 0; dx < poolSize && (x + dx) < adjustedNX; dx++) {
				maxVal = std::max(maxVal, combined[(y + dy) * adjustedNY + (x + dx)]);
				}
			}
			pooled.push_back(maxVal);
		}
	}

	// Calculate the dimensions of the pooled output
	int pooledHeight = std::ceil(static_cast<float>(adjustedNY) / poolSize);
	int pooledWidth = std::ceil(static_cast<float>(adjustedNX) / poolSize);

	// Print out the pooled matrix
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");
	for (size_t i = 0; i < 78; i++) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "-");
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");

	for (int y = 0; y < pooledHeight; y++) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "| ");
		for (int x = 0; x < pooledWidth; x++) {
			int index = y * pooledWidth + x;
			if (pooled[index] == 0) {
				ierr = PetscPrintf(PETSC_COMM_WORLD, " ");
			} else {
				ierr = PetscPrintf(PETSC_COMM_WORLD, "%d", static_cast<int>(pooled[index]));
			}
		}
		ierr = PetscPrintf(PETSC_COMM_WORLD, " |\n");
	}

	for (size_t i = 0; i < 78; i++) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "-");
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");
}


bool NeuronGrowth::ReadPointData(const string& filename, vector<float>& dataVector1, vector<float>& dataVector2, vector<float>& dataVector3, vector<float>& dataVector4, vector<float>& dataVector5) {
	ifstream file(filename);
	if (!file.is_open()) {
		cerr << "Failed to open file: " << filename << endl;
		return false;
	}

	string line;
	while (getline(file, line)) {
		// Look for the POINT_DATA section
		if (line.find("POINT_DATA") != string::npos) {
			return ParsePointData(file, dataVector1, dataVector2, dataVector3, dataVector4, dataVector5);
		}
	}

	cerr << "POINT_DATA section not found." << endl;
	return false;
}

bool NeuronGrowth::ParsePointData(ifstream& file, vector<float>& dataVector1, vector<float>& dataVector2, vector<float>& dataVector3, vector<float>& dataVector4, vector<float>& dataVector5) {
	string line;
	int vectorCounter = 0;

	while (getline(file, line) && vectorCounter < 5) {
		stringstream ss(line);
		string word;
		ss >> word;

		if (word == "SCALARS" || word == "VECTORS") {
			string fieldName, fieldType;
			ss >> fieldName >> fieldType;

			// Skip the lookup table line if it's a scalar
			if (word == "SCALARS") {
				getline(file, line);  // This should read the LOOKUP_TABLE line
			}

			// Read and store the data
			vector<float> tempData;
			while (getline(file, line) && !line.empty()) {
				istringstream iss(line);
				float value;
				while (iss >> value) {
					tempData.push_back(value);
				}
			}

			// Assign to the correct vector based on the counter
			switch (vectorCounter) {
				case 0:
					dataVector1 = tempData;
					break;
				case 1:
					dataVector2 = tempData;
					break;
				case 2:
					dataVector3 = tempData;
					break;
				case 3:
					dataVector4 = tempData;
					break;
				case 4:
					dataVector5 = tempData;
					break;
			}
			vectorCounter++;
		}
	}

	return true;
}

void NeuronGrowth::SetupExternalCue(const string& caseType, vector<vector<vector<int>>>& externalCues, const int& NX_fine, const int& NY_fine) {

	if (caseType == "A") {
		if (n < 35000) {
			externalCues = {{{NX_fine/2, 0+5}, {NX_fine/2, NY_fine-5}}};
		} else if (n < 60000) {
			externalCues = {{{}, {NX_fine-5, NY_fine-5}}};
		} else {
			externalCues = {{{}, {NX_fine*2/3, NY_fine-5}}};
		}
	} else if (caseType == "B") {
		if (n < 50000) {
			externalCues = {{{NX_fine*2/3, NY_fine-5}}};
		} else if (n < 90000) {
			externalCues = {{{NX_fine-5, NY_fine-5}}};
		} else if (n < 120000) {
			externalCues = {{{NX_fine*2/3, NY_fine-5}}};
		} else {
			externalCues = {{{NX_fine-5, NY_fine-5}}};
		}
	} else if (caseType == "C") {
		if (n < 40000) {
			externalCues = {{{NX_fine*2/3, 0+5}, {NX_fine-5, NY_fine*4/5}, {0+5, NY_fine*2/3}}};
		} else {
			externalCues = {{{NX_fine*3/4, 0+5}, {NX_fine-5, NY_fine-5}, {0+5, NY_fine}}};
		}
	} else if (caseType == "D") {
		if (n < 50000) {
			externalCues = {{{NX_fine/3, 0+5}, {NX_fine-5, 0+5}, {NX_fine/3, NY_fine-5}}};
		} else if (n < 90000) {
			externalCues = {{{NX_fine/2, 0+5}, {0+5, NY_fine-5}}};
		} else {
			externalCues = {{{NX_fine/2, 0+5}, {NX_fine*2/5, NY_fine-5}}};
		}
	} else if (caseType == "E") {
		if (n < 80000) {
			externalCues = {{{0+5, NY_fine/3}, {NX_fine-5, NY_fine*2/3}}};
		} else if (n < 110000) {
			externalCues = {{{0+5, NY_fine/3}, {NX_fine-5, NY_fine/3}}};
		} else {
			externalCues = {{{0+5, NY_fine-5}, {NX_fine-5, NY_fine*2/3}}};
		}
	} else if (caseType == "K") {
		if (n < 20000) {
			externalCues = {{{NX_fine/3, NY_fine-5}},
					{{NX_fine-5, NY_fine/4}}};
		} else if (n < 80000) {
			externalCues = {{{0+5, NY_fine*3/4}},
					{{NX_fine-5, NY_fine/4}}};
		} else if (n < 90000) {
			externalCues = {{{0+5, NY_fine*3/4}},
					{{NX_fine-5, NY_fine/5}, {NX_fine-5, NY_fine/2}}};
		} else if (n < 110000) {
			externalCues = {{{0+5, NY_fine*3/4}},
					{{NX_fine*4/5, NY_fine/5}, {NX_fine-5, NY_fine/3}}};
		} else {
			externalCues = {{{0+5, NY_fine*3/4}},
					{{NX_fine*4/5, 0+5}, {NX_fine-5, NY_fine/3}}};
		}
		seed_radius = 12;
	} else if (caseType == "L") {
		if (n < 50000) {
			externalCues = {{{NX_fine-5, NY_fine/2}},
					{{NX_fine*2/3, 0+5}}};
		} else if (n < 10000) {
			externalCues = {{{NX_fine-5, NY_fine/2}, {NX_fine*2/3, NY_fine-5}},
					{{NX_fine*2/3, 0+5}}};
		} else {
			externalCues = {{{NX_fine-5, NY_fine*2/5}},
					{{NX_fine*2/3, 0+5}}};
		}
		seed_radius = 12;
	} else if (caseType == "M") {
		if (n < 50000) {
			externalCues = {{{NX_fine/2, 0+5}},
					{{0+5, NY_fine/6}}};
		} else if (n < 80000) {
			externalCues = {{{NX_fine-5, 0+5}},
					{{0+5, NY_fine/3}, {NX_fine/4, 0+5}}};
		} else if (n < 120000) {
			externalCues = {{{NX_fine*3/5, 0+5}},
					{{0+5, NY_fine/4}, {NX_fine/4, 0+5}}};
		} else {
			externalCues = {{{NX_fine*3/5, 0+5}},
					{{0+5, NY_fine/5}, {NX_fine/4, 0+5}}};
		}
		seed_radius = 12;
	} else if (caseType == "N") {
		if (n < 60000) {
			externalCues = {{{0+5, 0+5}, {NX_fine-5, NY_fine-5}},
					{{0+5, NY_fine-5}, {NX_fine*2/3, 0+5}},
					{{NX_fine*4/5, 0+5}, {NX_fine-5, NY_fine-5}}};
		} else if (n < 80000){
			externalCues = {{{0+5, 0+5}, {NX_fine-5, NY_fine-5}},
					{{0+5, NY_fine-5}, {NX_fine*2/3, 0+5}},
					{{NX_fine-5, 0+5}, {0+5, 0+5}, {NX_fine-5, NY_fine-5}}};
		} else if (n < 130000){
			externalCues = {{{NX_fine/5, 0+5}},
					{{0+5, NY_fine-5}},
					{{NX_fine-5, 0+5}, {0+5, 0+5}, {NX_fine-5, NY_fine-5}}};
		} else if (n < 150000){
			externalCues = {{{NX_fine/3, 0+5}},
					{{}},
					{{NX_fine-5, 0+5}, {0+5, 0+5}, {NX_fine-5, NY_fine-5}}};
		} else {
			externalCues = {{{NX_fine/2, 0+5}},
					{{}},
					{{NX_fine-5, 0+5}, {0+5, 0+5}, {NX_fine-5, NY_fine-5}}};
		}
		seed_radius = 12;
	} else if (caseType == "O") {
		if (n < 50000) {
			externalCues = {{{0+5, NY_fine*4/5}, {NX_fine/2, 0+5}, {NX_fine-5, NY_fine*3/4}},
					{{0+5, 0+5}},
					{{NX_fine*2/3, NY_fine-5}}};
		} else if (n < 100000){
			externalCues = {{{0+5, NY_fine*3/5}, {NX_fine*2/3, 0+5}, {NX_fine-5, NY_fine*2/3}},
					{{0+5, NY_fine/4}},
					{{NX_fine*2/3, NY_fine-5}}};
		} else {
			externalCues = {{{0+5, NY_fine-5}},
					{{0+5, NY_fine/5}},
					{{NX_fine/2, NY_fine-5}}};
		}
		seed_radius = 12;
	} else if (caseType == "P") {
		if (n < 60000) {
			externalCues = {{{NX_fine*2/3, NY_fine-5}}};
		} else if (n < 90000) {
			externalCues = {{{NX_fine/3, NY_fine-5}, {NX_fine-5, NY_fine}}};
		} else {
			externalCues = {{{0+5, NY_fine-5}, {NX_fine-5, NY_fine*2/3}}};
		}
	} else if (caseType == "Q") {
		if (n < 20000) {
			externalCues = {{{0+5, NY_fine-5}}};
		} else if (n < 80000) {
			externalCues = {{{NX_fine/2, NY_fine-5}}};
		} else if (n < 100000) {
			externalCues = {{{NX_fine-5, NY_fine-5}}};
		} else {
			externalCues = {{{0+5, NY_fine*3/4}, {NX_fine-5, NY_fine-5}}};
		}
	} else if (caseType == "R") {
		if (n < 20000) {
			externalCues = {{{NX_fine-5, NY_fine-5}}};
		} else if (n < 80000) {
			externalCues = {{{NX_fine-5, NY_fine/2}}};
		} else if (n < 110000) {
			externalCues = {{{NX_fine-5, 0}, {NX_fine-5, NY_fine-5}}};
		} else {
			externalCues = {{{NX_fine*3/4, 0}, {NX_fine-5, NY_fine/2}, {NX_fine-5, NY_fine-5}}};
		}
	} else if (caseType == "S") {
		if (n < 60000) {
			externalCues = {{{NX_fine-5, 0+5}, {NX_fine/2, NY_fine-5}}};
		} else {
			externalCues = {{{NX_fine-5, 0+5}, {NX_fine/3, NY_fine-5}, {NX_fine-5, NY_fine-5}}};
		}
	} else if (caseType == "T") {
		if (n < 60000) {
			externalCues = {{{NX_fine*2/3, NY_fine-5}}};
		} else {
			externalCues = {{{NX_fine*4/5, NY_fine-5}}};
		}
	} else if (caseType == "U") {
		if (n < 80000) {
			externalCues = {{{NX_fine-5, NY_fine-5}}};
		} else if (n < 100000) {
			externalCues = {{{NX_fine*3/4, NY_fine-5}}};
		} else {
			externalCues = {{{NX_fine/3, NY_fine-5}, {NX_fine-5, NY_fine-5}}};
		}
	} else if (caseType == "V") {
		if (n < 50000) {
			externalCues = {{{NX_fine/2, NY_fine-5}}};
		} else {
			externalCues = {{{NX_fine*3/5, NY_fine-5}}};
		}
	} else if (caseType == "W") {
		if (n < 100000) {
			externalCues = {{{NX_fine/2, 0+5}, {NX_fine/2, NY_fine-5}}};
		} else {
			externalCues = {{{NX_fine/2+5, 0+5}, {NX_fine/2, NY_fine-5}}};
		}
	} else if (caseType == "X") {
		if (n < 80000) {
			externalCues = {{{NX_fine-5, NY_fine/2}}};
		} else if (n < 110000) {
			externalCues = {{{NX_fine-5, NY_fine*3/5}}};
		} else {
			externalCues = {{{NX_fine-5, NY_fine-5}}};
		}
	} else if (caseType == "Y") {
		if (n < 40000) {
			externalCues = {{{NX_fine/2, NY_fine-5}}};
		} else if (n < 80000) {
			externalCues = {{{0+5, NY_fine-5}}};
		} else if (n < 150000) {
			externalCues = {{{0+5, NY_fine/3}}};
		} else if (n < 200000) {
			externalCues = {{{0+5, NY_fine/3}}};
		} else {
			externalCues = {{{0+5, NY_fine-5}}};
		}
	} else if (caseType == "F") {
		if (n < 20000) {
			externalCues = {{{NX_fine/2, NY_fine-5}}};
		} else if (n < 100000) {
			externalCues = {{{0+5, NY_fine-5}}};
		} else {
			externalCues = {{{0+5, NY_fine-5}, {NX_fine/2, NY_fine-5}}};
		}
	} else if (caseType == "G") {
		if (n < 20000) {
			externalCues = {{{0+5, NY_fine/2}, {0+5, 0+5}, {NX_fine-5, 0+5}}};
		} else if (n < 40000) {
			externalCues = {{{0+5, NY_fine/2}, {0+5, 0+5}, {NX_fine/2, 0+5}, {NX_fine-5, 0+5}}};
		} else {
			externalCues = {{{0+5, NY_fine/2}, {0+5, 0+5}, {NX_fine/2, 0+5}, {NX_fine-5, 0+5}}};
		}
	} else {}
	
	if (n > 250000) {
		c_opt = 0;
	} else {
		c_opt = 1;
	}
}

PetscErrorCode FormFunction_phi(SNES snes, Vec x, Vec F, void *ctx)
{
	PetscErrorCode ierr;
	Vec P_seq;
	PetscScalar *Parray;
	VecScatter scatter_ctx1;
	ierr = VecScatterCreateToAll(x, &scatter_ctx1, &P_seq); CHKERRQ(ierr);
	ierr = VecScatterBegin(scatter_ctx1, x, P_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecScatterEnd(scatter_ctx1, x, P_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGetArray(P_seq, &Parray); CHKERRQ(ierr);
	ierr = VecSet(F, 0.0); CHKERRQ(ierr);

	NeuronGrowth *user = (NeuronGrowth *)ctx;

	/*Build linear system in each process*/
	uint ind(0); // pre-calculated variable index
	for (size_t e = 0; e < user->bzmesh_process.size(); e++) {
		uint nen = user->bzmesh_process[e].IEN.size(); // 16 supporting cp for 2D case

		vector<vector<float>> EMatrixSolve;	EMatrixSolve.clear();
		vector<float> EVectorSolve;		EVectorSolve.clear();
		EMatrixSolve.resize(nen);
		EVectorSolve.resize(nen);

		user->eleVal.clear();
		user->eleVal.resize(8);
		for (size_t i = 0; i < nen * 1; i++) {
			EMatrixSolve[i].resize(nen);
			if (i < 8)
				user->eleVal[i].resize(nen);
		}

		// preparing element stiffness matrix and vector, extract values from control mesh
		for (size_t i = 0; i < nen; i++) {

			EVectorSolve[i] = 0.0;
			for (size_t j = 0; j < nen; j++) {
				EMatrixSolve[i][j] = 0.0;
			}

			user->eleVal[0][i] = Parray[user->bzmesh_process[e].IEN[i]];		// elePhiGuess
			user->eleVal[1][i] = user->phi[user->bzmesh_process[e].IEN[i]];		// elePhi
			user->eleVal[4][i] = user->theta[user->bzmesh_process[e].IEN[i]];	// eleTheta
			user->eleVal[6][i] = 0;							// eleEpsilon
			user->eleVal[7][i] = 0;							// eleEpsilonP
		}

		// loop through gaussian quadrature points
		for (size_t i = 0; i < user->Gpt.size(); i++) {
			for (size_t j = 0; j < user->Gpt.size(); j++) {
				user->EvaluateOrientation(nen, user->pre_Nx[ind], user->pre_dNdx[ind], user->eleVal[1], user->eleVal[4], user->eleVal[6], user->eleVal[7]);
				user->ElementEvaluationAll_phi(nen, user->pre_Nx[ind], user->pre_dNdx[ind], user->eleVal, user->vars);
				// vars  0   1    2     3      4      5     6     7      8     9      10      11      12     13    14    15      16     17     
				// float C1, C0, elePG, dPGdx, dPGdy, eleP, eleS, eleTb, eleE, eleTp, dThedx, dThedy, eleEP, dAdx, dAdy, eleEEP, dAPdx, dAPdy;
				
				// loop through control points
				for (size_t m = 0; m < nen; m++) {
					// terma2 = - eleEP * eleEP * (dPGdx * dNdx[m][0] + dPGdy * dNdx[m][1]);
					// termadx = - dAdx * eleEEP * dPGdy * dNdx[m][0];
					// termady = - dAdy * eleEEP * dPGdx * dNdx[m][1];
					// termdbl = (- elePG * elePG * elePG + (1 - C1) * elePG * elePG + C1 * elePG) * Nx[m];
					// EVectorSolve[m] += (elePG * Nx[m] - user->dt * user->M_phi * (terma2 - termadx + termady + termdbl) - eleP * Nx[m]) * detJ;
					EVectorSolve[m] += (user->vars[2] * user->pre_Nx[ind][m] - user->dt * user->pre_eleMp[ind] *\
						((- user->vars[12] * user->vars[12] * (user->vars[3] * user->pre_dNdx[ind][m][0] + user->vars[4] * user->pre_dNdx[ind][m][1])) -\
						(- user->vars[13] * user->vars[15] * user->vars[4] * user->pre_dNdx[ind][m][0]) +\
						(- user->vars[14] * user->vars[15] * user->vars[3] * user->pre_dNdx[ind][m][1]) +\
						(- user->vars[2] * user->vars[2] * user->vars[2] + (1 - user->pre_C1[ind]) * user->vars[2] * user->vars[2] + user->pre_C1[ind] * user->vars[2]) * user->pre_Nx[ind][m])
						- user->vars[5] * user->pre_Nx[ind][m]) * user->pre_detJ[ind];
					// loop through 16 control points
					for (size_t n = 0; n < nen; n++) {
						// terma2 = - eleEP * eleEP * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1]);
						// termadx = - dAdx * eleEEP * dNdx[m][1] * dNdx[n][0];
						// termady = - dAdy * eleEEP * dNdx[m][0] * dNdx[n][1];
						// termdbl = (- 3 * elePG * elePG + 2 * (1 - C1) * elePG + C1 * Nx[m]) * Nx[n];
						// EMatrixSolve[m][n] += (Nx[m] * Nx[n] - user->dt * user->M_phi * (terma2 - termadx + termady + termdbl)) * detJ;
						EMatrixSolve[m][n] += (user->pre_Nx[ind][m] * user->pre_Nx[ind][n] - user->dt * user->pre_eleMp[ind] *
							((- user->vars[12] * user->vars[12] * (user->pre_dNdx[ind][m][0] * user->pre_dNdx[ind][n][0] + user->pre_dNdx[ind][m][1] * user->pre_dNdx[ind][n][1])) // terma2
							- (- user->vars[13] * user->vars[15] * user->pre_dNdx[ind][m][1] * user->pre_dNdx[ind][n][0]) // termadx
							+ (- user->vars[14] * user->vars[15] * user->pre_dNdx[ind][m][0] * user->pre_dNdx[ind][n][1]) // termady
							+ (- 3 * user->vars[2] * user->vars[2] + 2 * (1 - user->pre_C1[ind]) * user->vars[2] + user->pre_C1[ind] * user->pre_Nx[ind][m]) * user->pre_Nx[ind][n]) // termdbl
							) * user->pre_detJ[ind];
					}
				}
				ind += 1; // incrementing index for extracting pre-calculated variables
			}
		}

		/*Apply Boundary Condition*/
		for (size_t i = 0; i < nen; i++) {
			int A = user->bzmesh_process[e].IEN[i];
			if (user->cpts[A].label == 1) {// domain boundary
				user->ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
			}
		}

		user->ResidualAssembly(EVectorSolve, user->bzmesh_process[e].IEN, F);
		user->MatrixAssembly(EMatrixSolve, user->bzmesh_process[e].IEN, user->J);
	}

	ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

	ierr = MatAssemblyBegin(user->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(user->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	ierr = VecRestoreArray(P_seq, &Parray); CHKERRQ(ierr);
	ierr = VecScatterDestroy(&scatter_ctx1); CHKERRQ(ierr);
	ierr = VecDestroy(&P_seq); CHKERRQ(ierr);

	return 0;
}


PetscErrorCode FormFunction_phi_test(SNES snes, Vec x, Vec F, void *ctx)
{
	PetscErrorCode ierr;
	
	/*To get the global x needed for IGA, each element needs control points around that element*/
	/*very challenging if just access local vector*/
	Vec P_seq;
	PetscScalar *Parray;
	VecScatter scatter_ctx1;
	ierr = VecScatterCreateToAll(x, &scatter_ctx1, &P_seq); CHKERRQ(ierr);
	ierr = VecScatterBegin(scatter_ctx1, x, P_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecScatterEnd(scatter_ctx1, x, P_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGetArray(P_seq, &Parray); CHKERRQ(ierr);

	ierr = VecSet(F, 0.0); CHKERRQ(ierr);

	NeuronGrowth *user = (NeuronGrowth *)ctx;

	/*Build linear system in each process*/
	int ind(0); // pre-calculated variable index
	for (size_t e = 0; e < user->bzmesh_process.size(); e++) {
		uint nen = user->bzmesh_process[e].IEN.size(); // 16 supporting cp for 2D case

		vector<vector<float>> EMatrixSolve;	EMatrixSolve.clear();
		vector<float> EVectorSolve;		EVectorSolve.clear();
		EMatrixSolve.resize(nen);
		EVectorSolve.resize(nen);

		for (size_t i = 0; i < nen * 1; i++) {
			EMatrixSolve[i].resize(nen);
		}
		
		// preparing element stiffness matrix and vector, extract values from control mesh
		vector<float> elePhiGuess(nen, 0);
		for (size_t i = 0; i < nen; i++) {
			EVectorSolve[i] = 0.0;
			for (size_t j = 0; j < nen; j++) {
				EMatrixSolve[i][j] = 0.0;
			}
			elePhiGuess[i] = Parray[user->bzmesh_process[e].IEN[i]];		// elePhiGuess
		}

		// loop through gaussian quadrature points
		for (size_t i = 0; i < user->Gpt.size(); i++) {
			for (size_t j = 0; j < user->Gpt.size(); j++) {				
				user->ElementEvaluationAll_phi(nen, user->pre_Nx[ind], user->pre_dNdx[ind], elePhiGuess, user->vars);
				
				// loop through control points
				for (size_t m = 0; m < nen; m++) {
					// terma2 = - eleEP * eleEP * (dPGdx * dNdx[m][0] + dPGdy * dNdx[m][1]);
					// termadx = - dAdx * eleEEP * dPGdy * dNdx[m][0];
					// termady = - dAdy * eleEEP * dPGdx * dNdx[m][1];
					// termdbl = (- elePG * elePG * elePG + (1 - C1) * elePG * elePG + C1 * elePG) * Nx[m];
					// EVectorSolve[m] += (elePG * Nx[m] - user->dt * user->M_phi * (terma2 - termadx + termady + termdbl) - eleP * Nx[m]) * detJ;
			
					EVectorSolve[m] += (user->vars[0] * user->pre_Nx[ind][m] - user->dt * user->pre_eleMp[ind] *\
						((- user->pre_eleEP[ind] * user->pre_eleEP[ind] * (user->vars[1] * user->pre_dNdx[ind][m][0] + user->vars[2] * user->pre_dNdx[ind][m][1])) -\
						(- user->pre_dAdx[ind] * user->pre_eleEEP[ind] * user->vars[2] * user->pre_dNdx[ind][m][0]) +\
						(- user->pre_dAdy[ind] * user->pre_eleEEP[ind] * user->vars[1] * user->pre_dNdx[ind][m][1]) +\
						(- user->vars[0] * user->vars[0] * user->vars[0] + (1 - user->pre_C1[ind]) * user->vars[0] * user->vars[0] + user->pre_C1[ind] * user->vars[0]) * user->pre_Nx[ind][m])
						- user->pre_eleP[ind] * user->pre_Nx[ind][m]) * user->pre_detJ[ind];

					// loop through 16 control points
					for (size_t n = 0; n < nen; n++) {
						// terma2 = - eleEP * eleEP * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1]);
						// termadx = - dAdx * eleEEP * dNdx[m][1] * dNdx[n][0];
						// termady = - dAdy * eleEEP * dNdx[m][0] * dNdx[n][1];
						// termdbl = (- 3 * elePG * elePG + 2 * (1 - C1) * elePG + C1 * Nx[m]) * Nx[n];
						// EMatrixSolve[m][n] += (Nx[m] * Nx[n] - user->dt * user->M_phi * (terma2 - termadx + termady + termdbl)) * detJ;
						
						EMatrixSolve[m][n] += (user->pre_Nx[ind][m] * user->pre_Nx[ind][n] - user->dt * user->pre_eleMp[ind] *
							((- user->pre_eleEP[ind] * user->pre_eleEP[ind] * (user->pre_dNdx[ind][m][0] * user->pre_dNdx[ind][n][0] + user->pre_dNdx[ind][m][1] * user->pre_dNdx[ind][n][1])) // terma2
							- (- user->pre_dAdx[ind] * user->pre_eleEEP[ind] * user->pre_dNdx[ind][m][1] * user->pre_dNdx[ind][n][0]) // termadx
							+ (- user->pre_dAdy[ind] * user->pre_eleEEP[ind] * user->pre_dNdx[ind][m][0] * user->pre_dNdx[ind][n][1]) // termady
							+ (- 3 * user->vars[0] * user->vars[0] + 2 * (1 - user->pre_C1[ind]) * user->vars[0] + user->pre_C1[ind] * user->pre_Nx[ind][m]) * user->pre_Nx[ind][n]) // termdbl
							) * user->pre_detJ[ind];
					}
				}
				ind += 1; // incrementing index for extracting pre-calculated variables
			}
		}

		/*Apply Boundary Condition*/
		for (size_t i = 0; i < nen; i++) {
			int A = user->bzmesh_process[e].IEN[i];
			if (user->cpts[A].label == 1) { // domain boundary
				user->ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
			}
		}

		user->ResidualAssembly(EVectorSolve, user->bzmesh_process[e].IEN, F);
		user->MatrixAssembly(EMatrixSolve, user->bzmesh_process[e].IEN, user->J);
	}

	ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

	ierr = MatAssemblyBegin(user->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(user->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	if (user->comRank == 0) {
		user->check_itr += 1;
	}
	ierr = VecRestoreArray(P_seq, &Parray); CHKERRQ(ierr);
	ierr = VecScatterDestroy(&scatter_ctx1); CHKERRQ(ierr);
	ierr = VecDestroy(&P_seq); CHKERRQ(ierr);

	return 0;
}

PetscErrorCode FormJacobian_phi(SNES snes, Vec x, Mat J, Mat P, void *ctx) 
{
	NeuronGrowth *user = (NeuronGrowth *)ctx;
	PetscErrorCode ierr;
	ierr = MatCopy(user->J, J, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
	return 0;
}

PetscErrorCode MySNESMonitor(SNES snes, PetscInt its, PetscReal fnorm, PetscViewerAndFormat *vf)
{
	PetscFunctionBeginUser;
	SNESMonitorDefaultShort(snes, its, fnorm, vf);
	
	// Set up the KSP solver to monitor iterations
	KSP ksp;
	SNESGetKSP(snes, &ksp);

	// PetscOptionsSetValue(NULL, "-ksp_monitor", "");  // Enables KSP iteration monitoring
	PetscOptionsSetValue(NULL, "-ksp_monitor_singular_value", "");  // Enables KSP iteration monitoring
	// PetscOptionsSetValue(NULL, "-ksp_monitor_solution", "");  // Enables monitoring of the solution

	// Get the number of iterations taken by the KSP solver
	PetscInt numIterations;
	// KSPGetIterationNumber(ksp, &numIterations);
	KSPGetTotalIterations(ksp, &numIterations);
	PetscPrintf(PETSC_COMM_WORLD, "     - KSP Iterations: %d\n", numIterations);

	PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode CleanUpSolvers(NeuronGrowth& NG)
{
	// NG.ierr = SNESDestroy(&NG.snes_phi); CHKERRQ(NG.ierr);
	// NG.ierr = MatDestroy(&NG.J); CHKERRQ(NG.ierr);
	// NG.ierr = VecDestroy(&NG.temp_phi); CHKERRQ(NG.ierr);

	NG.ierr = KSPDestroy(&NG.ksp_phi); CHKERRQ(NG.ierr);
	NG.ierr = MatDestroy(&NG.GK_phi); CHKERRQ(NG.ierr);
	NG.ierr = VecDestroy(&NG.GR_phi); CHKERRQ(NG.ierr);
	NG.ierr = VecDestroy(&NG.temp_phi); CHKERRQ(NG.ierr);
	
	NG.ierr = KSPDestroy(&NG.ksp_syn); CHKERRQ(NG.ierr);
	NG.ierr = MatDestroy(&NG.GK_syn); CHKERRQ(NG.ierr);
	NG.ierr = VecDestroy(&NG.GR_syn); CHKERRQ(NG.ierr);
	NG.ierr = VecDestroy(&NG.temp_syn); CHKERRQ(NG.ierr);

	NG.ierr = KSPDestroy(&NG.ksp_tub); CHKERRQ(NG.ierr);	
	NG.ierr = MatDestroy(&NG.GK_tub); CHKERRQ(NG.ierr);
	NG.ierr = VecDestroy(&NG.GR_tub); CHKERRQ(NG.ierr);
	NG.ierr = VecDestroy(&NG.temp_tub); CHKERRQ(NG.ierr);

	return NG.ierr;
}

int RunNG(int& n_bzmesh, vector<vector<int>> ele_process_in, vector<Vertex2D>& cpts_initial, const vector<Element2D>& tmesh_initial,
	vector<Vertex2D>& cpts_fine, const vector<Element2D>& tmesh_fine, vector<Vertex2D>& cpts,
	vector<Vertex2D>& prev_cpts, vector<Vertex2D>& prev_cpts_fine, const vector<Element2D>& tmesh,
	string path_in, string path_out, int& iter, int end_iter_in, vector<vector<float>> &NGvars, 
	int &NX, int &NY, vector<array<float, 2>>& seed, int &originX, int &originY,
	vector<int> &rfid, vector<int> &rftype, bool &localRefine, const string& caseType)
{	
	/*==============================================================================*/
	// Declare and initialize Neuron Growth simulation
	NeuronGrowth NG;
	NG.SetVariables(path_in + "simulation_parameters.txt");

	// Set the number of iterations, number of neurons, and the end iteration from provided variables
	NG.n = iter; 			// Assign iteration count
	NG.numNeuron = seed.size();  	// Set the number of neurons based on the size of 'seed'
	NG.end_iter = end_iter_in;  	// Set the ending iteration

	// Initialize vertex clouds for the current, fine, and previous configurations
	Vertex2DCloud cloud(cpts);        		// Cloud for current points
	Vertex2DCloud cloud_fine(cpts_fine); 		// Cloud for finer resolution points
	Vertex2DCloud cloud_prev(prev_cpts);  		// Cloud for previous points
	Vertex2DCloud cloud_prev_fine(prev_cpts_fine);  // Cloud for previous finer resolution points
	// Initialize KD-Trees for the current, fine, and previous vertex clouds
	KDTree kdTree(2 /* dim */, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	KDTree kdTree_fine(2 /* dim */, cloud_fine, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	KDTree kdTree_prev(2 /* dim */, cloud_prev, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	KDTree kdTree_prev_fine(2 /* dim */, cloud_prev_fine, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	// Build indexes for KD-Trees to optimize search operations
	kdTree.buildIndex();
	kdTree_fine.buildIndex();
	kdTree_prev.buildIndex();
	kdTree_prev_fine.buildIndex();

	/*==============================================================================*/
	// Pre-allocate variables for tip detections (will be used every iter)
	int NX_fine(2 * NX + 1), NY_fine(2 * NY + 1);
	int totalSize = NX_fine * NY_fine;  // Total number of elements in the flattened 2D grid
	// Allocate neurons matrix as a single flat vector
	vector<int> neurons(totalSize, 0);
	// Allocate distances 3D matrix: vector of 2D matrices, each represented by a flat vector<float>
	vector<vector<float>> distances(NG.numNeuron, vector<float>(totalSize, 0.0f));
	// Allocate vectors for phi_fine, id, tip, and localMaximaMatrix
	vector<float> phi_fine(cpts_fine.size(), 0), id(cpts_fine.size(), 0), tip(cpts_fine.size(), 0), localMaximaMatrix(cpts_fine.size(), 0);
	// Allocate geodist and axonTip matrices
	vector<vector<float>> geodist(NG.numNeuron, vector<float>(totalSize, 0));
	PetscPrintf(PETSC_COMM_WORLD, "Pre-allocated vectors!-------------------------------------------------------\n");
	// For experimental comparison cases. "No" is not running compare cases.
	vector<vector<vector<int>>> externalCues;
	NG.SetupExternalCue(caseType, externalCues, NX_fine, NY_fine);

	// Initialize the Neuron Growth problem with the provided parameters and KDTree for previous points
	NG.InitializeProblemNG(n_bzmesh, cpts, prev_cpts, kdTree_prev, cpts_fine, prev_cpts_fine, kdTree_prev_fine, NGvars, seed, NX, NY);
	// Assign processing elements for the simulation
	NG.AssignProcessor(ele_process_in);
	// Synchronize and check MPI element assignments, print out information in order
	for (int i = 0; i < NG.nProcess; i++) {
		if (i == NG.comRank) {
			cout << "comRank: " << NG.comRank << "/" << NG.nProcess << " - with element process size: " << NG.ele_process.size() << endl;
			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
		} else {
			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
		}
	}
	// Read Bezier mesh for the simulation setup
	NG.ReadBezierElementProcess(path_in);
	PetscPrintf(PETSC_COMM_WORLD, "Read bzmesh!-----------------------------------------------------------------\n");
	// (Optional) Set initial guess for nonlinear solver based on the current configuration
	// NG.ToPETScVec(NG.phi, NG.temp_phi);
	// PetscPrintf(PETSC_COMM_WORLD, "Set initial guess!-----------------------------------------------------------\n");
	
	/*==============================================================================*/
	// Write initial variables
	if (NG.n == 0) {
		NG.VisualizeVTK_PhysicalDomain_All(0, path_out);
		PetscPrintf(PETSC_COMM_WORLD, "Saving all variables!--------------------------------------------------------\n");	
	}

	/*==============================================================================*/
	// Pre-calculation phase to optimize computational cost
	NG.prepareBasis();
	PetscPrintf(PETSC_COMM_WORLD, "Prepared basis!--------------------------------------------------------------\n");
	// Aggregate gradient of phi0 across processes
	NG.ierr = MPI_Allreduce(&NG.sum_grad_phi0_local, &NG.sum_grad_phi0_global, 1, MPI_FLOAT, MPI_SUM, PETSC_COMM_WORLD);
	CHKERRQ(NG.ierr); // Check for errors during MPI call
	PetscPrintf(PETSC_COMM_WORLD, "Calculated sum grad phi0!----------------------------------------------------\n");
	// Prepare source terms for computation
	NG.prepareTerm_source();
	PetscPrintf(PETSC_COMM_WORLD, "Prepared source terms!-------------------------------------------------------\n");

	/*==============================================================================*/
	// Main time iterations
	while (iter <= NG.end_iter) {
		NG.n = iter;
		tic();

		/*==============================================================================*/
		// Neuron identification and tip detection
		phi_fine = NG.InterpolateValues_closest(NG.phi, kdTree, cpts_fine);

		NG.IdentifyNeurons(phi_fine, neurons, NG.prev_id, seed, NX*2, NY*2, originX, originY);
		// Apply the transformation directly to NG.prev_id
		for (size_t i = 0; i < NG.prev_id.size(); i++) {
			int prevVal = NG.prev_id[i];
			int neuronVal = neurons[i];
			NG.prev_id[i] = (prevVal == 0 ? neuronVal : min(prevVal, neuronVal));
		}
		NG.DetectConnections(neurons, NX_fine, NY_fine);

		vector<int> centroidIndices;

		distances = NG.ExploreGridAndCalculateDistances(neurons, seed, originX, originY, NX_fine, NY_fine);
		// NG.DetectTipsMulti(phi_fine, neurons, NG.numNeuron, tip, NX_fine, NY_fine);

		if (caseType != "No") {
			NG.DetectTipsMulti(phi_fine, neurons, NG.numNeuron, tip, NX_fine, NY_fine, centroidIndices);
			// NG.AdjustNearestTip(tip, NX_fine, NY_fine, externalCues[0]);
			// localMaximaMatrix = tip;
			// localMaximaMatrix = NG.PickNearestTip(tip, NX_fine, NY_fine, externalCues[0], centroidIndices);
			localMaximaMatrix = NG.PickNearestTip(tip, neurons, NX_fine, NY_fine, externalCues, centroidIndices);
		} else {
			NG.DetectTipsMulti(phi_fine, neurons, NG.numNeuron, tip, NX_fine, NY_fine); // potential further opt can be done here
			localMaximaMatrix = NG.FindCentroidsOfLocalMaximaClusters(tip, NX_fine, NY_fine, centroidIndices);
			for (size_t i = 0; i < distances.size(); i++) {
				geodist[i] = distances[i];
				vector<float> maxGeodist = NG.ComputeMaxFilter(geodist[i], NX_fine, NY_fine, 5);
				float maxVal = -numeric_limits<float>::max();
				int maxInd = 0;
				for (size_t j = 0; j < centroidIndices.size(); ++j) {
					if (maxGeodist[centroidIndices[j]] > maxVal) {
						maxVal = maxGeodist[centroidIndices[j]];
						maxInd = centroidIndices[j];
					}
				}

				if (maxInd != 0) {
					localMaximaMatrix[maxInd] = -5;
				
				}
			}
		}

		NG.tips = NG.InterpolateValues_closest(localMaximaMatrix, kdTree_fine, cpts);

		toc(t_tip);
		tic();
		t_total += t_tip;
		
		/*==============================================================================*/
		/*Implcit Non-liear Newton Raphson solver for Phase field equation*/
		NG.phi_prev = NG.phi;	
		NG.preparePhaseField();

		float tol = 1e-4;
		float residual = tol + 1;
		int NR_itr = 0;
		while (residual > tol) {
			NG.ierr = MatZeroEntries(NG.GK_phi); CHKERRQ(NG.ierr);
			NG.ierr = VecSet(NG.GR_phi, 0); CHKERRQ(NG.ierr);
			NG.BuildLinearSystemProcessNG_phi();
			NG.ierr = VecAssemblyEnd(NG.GR_phi); CHKERRQ(NG.ierr);
			NG.ierr = MatAssemblyEnd(NG.GK_phi, MAT_FINAL_ASSEMBLY); CHKERRQ(NG.ierr);

			/*Petsc solver setting for phi*/
			if (NG.judge_phi == 0) {
				NG.ierr = KSPCreate(PETSC_COMM_WORLD, &NG.ksp_phi); CHKERRQ(NG.ierr);
				NG.ierr = KSPSetOperators(NG.ksp_phi, NG.GK_phi, NG.GK_phi); CHKERRQ(NG.ierr);
				NG.ierr = KSPGetPC(NG.ksp_phi, &NG.pc_phi); CHKERRQ(NG.ierr);
				NG.ierr = PCSetType(NG.pc_phi, PCBJACOBI); CHKERRQ(NG.ierr);
				NG.ierr = KSPSetType(NG.ksp_phi, KSPGMRES); CHKERRQ(NG.ierr);
				NG.ierr = KSPGMRESSetRestart(NG.ksp_phi, 100); CHKERRQ(NG.ierr);
				NG.ierr = KSPSetInitialGuessNonzero(NG.ksp_phi, PETSC_TRUE); CHKERRQ(NG.ierr);
				NG.ierr = KSPSetTolerances(NG.ksp_phi, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 100000); CHKERRQ(NG.ierr);
				NG.ierr = KSPSetFromOptions(NG.ksp_phi); CHKERRQ(NG.ierr);
				NG.ierr = KSPSetUp(NG.ksp_phi); CHKERRQ(NG.ierr);
				if (NG.n == 0) {
					NG.ierr = KSPView(NG.ksp_phi, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(NG.ierr);
				}
				NG.judge_phi = 1;
			}

			NG.ierr = KSPSolve(NG.ksp_phi, NG.GR_phi, NG.temp_phi); CHKERRQ(NG.ierr);
			PetscInt its_phi;
			NG.ierr = KSPGetIterationNumber(NG.ksp_phi, &its_phi); CHKERRQ(NG.ierr);
			KSPConvergedReason reason_phi;
			NG.ierr = KSPGetConvergedReason(NG.ksp_phi, &reason_phi); CHKERRQ(NG.ierr);
			if (reason_phi < 0) {
				PetscPrintf(PETSC_COMM_WORLD, "KSP phi not converging  | KSPreason:  %d | KSPiter: %d \n", reason_phi, its_phi); CHKERRQ(NG.ierr);
				return 3;	
			}

			NR_itr += 1;

			/*==============================================================================*/
			/*Collecting scattered Phi variable from all processors*/
			Vec temp_phi_seq;
			VecScatter ctx_phi;
			PetscScalar *_p;

			NG.ierr = VecScatterCreateToAll(NG.temp_phi, &ctx_phi, &temp_phi_seq); CHKERRQ(NG.ierr);
			NG.ierr = VecScatterBegin(ctx_phi, NG.temp_phi, temp_phi_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
			NG.ierr = VecScatterEnd(ctx_phi, NG.temp_phi, temp_phi_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
			NG.ierr = VecGetArray(temp_phi_seq, &_p);

			residual = 0;
			for (size_t i = 0; i < NG.phi.size(); i++) {
				NG.phi[i] = NG.phi[i] - PetscRealPart(_p[i]);
				residual = PetscMax(residual, abs(PetscRealPart(_p[i])));

				if (abs(NG.phi[i]) > 3) {
					PetscPrintf(PETSC_COMM_WORLD, "Incorrect diverging phi!-----------------------------------------------------\n");
					return 3;
				}
			}

			NG.ierr = VecRestoreArray(temp_phi_seq, &_p); CHKERRQ(NG.ierr);
			NG.ierr = VecScatterDestroy(&ctx_phi); CHKERRQ(NG.ierr);
			NG.ierr = VecDestroy(&temp_phi_seq); CHKERRQ(NG.ierr);

		}

		toc(t_phi);
		tic();
		t_total += t_phi;

		// /*==============================================================================*/
		// /*Implcit Non-liear SNES solver for Phase field equation*/
		// NG.phi_prev = NG.phi;	
		// NG.preparePhaseField();

		// NG.check_itr = 0;
		// if (NG.judge_phi == 0) { 
		// 	NG.ierr = SNESCreate(PETSC_COMM_WORLD, &NG.snes_phi); CHKERRQ(NG.ierr);
		// 	NG.ierr = SNESSetType(NG.snes_phi, SNESNEWTONLS); CHKERRQ(NG.ierr);
		// 	NG.ierr = SNESSetTolerances(NG.snes_phi, 1e-4, 1e-6, 1e-8, 100, 1000); CHKERRQ(NG.ierr);
		// 	PetscPrintf(PETSC_COMM_WORLD, "Set Tolerance!---------------------------------------------------------------\n");

		// 	SNESLineSearch linesearch;
		// 	NG.ierr = SNESGetLineSearch(NG.snes_phi, &linesearch); CHKERRQ(NG.ierr);
		// 	NG.ierr = SNESLineSearchSetType(linesearch, SNESLINESEARCHCP); CHKERRQ(NG.ierr);
		// 	PetscPrintf(PETSC_COMM_WORLD, "Set LineSearch!--------------------------------------------------------------\n");	
		// 	// NG.ierr = SNESSetFunction(NG.snes_phi, NULL, FormFunction_phi, &NG); CHKERRQ(NG.ierr);
		// 	NG.ierr = SNESSetFunction(NG.snes_phi, NULL, FormFunction_phi_test, &NG); CHKERRQ(NG.ierr);
		// 	NG.ierr = SNESSetJacobian(NG.snes_phi, NULL, NULL, FormJacobian_phi, &NG); CHKERRQ(NG.ierr);
		// 	PetscPrintf(PETSC_COMM_WORLD, "Set FormFunction and FormJacobian!-------------------------------------------\n");

		// 	// PetscViewerAndFormat *vf;
		// 	// PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf);
		// 	// SNESMonitorSet(NG.snes_phi, (PetscErrorCode(*)(SNES, PetscInt, PetscReal, void *))MySNESMonitor, vf, (PetscErrorCode(*)(void **))PetscViewerAndFormatDestroy);
		// 	// PetscPrintf(PETSC_COMM_WORLD, "Set Monitor!-----------------------------------------------------------------\n");
			
		// 	if (NG.n == 0)
		// 		NG.ierr = SNESView(NG.snes_phi, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(NG.ierr);
		// 	NG.judge_phi = 1;
		// }


		// NG.ierr = SNESSolve(NG.snes_phi, NULL, NG.temp_phi); CHKERRQ(NG.ierr);
		// PetscInt its_phi;
		// NG.ierr = SNESGetIterationNumber(NG.snes_phi, &its_phi); CHKERRQ(NG.ierr);
		// SNESConvergedReason reason_phi;
		// NG.ierr = SNESGetConvergedReason(NG.snes_phi, &reason_phi); CHKERRQ(NG.ierr);
		// if (reason_phi < 0) {
		// 	PetscPrintf(PETSC_COMM_WORLD, "SNES phi not converging ........ | SNESreason: %d | SNESiter: %d \n", reason_phi, its_phi); CHKERRQ(NG.ierr);	
		// }

		// // if (NG.comRank == 0) {
		// // 	cout << "Checking snes calls: " << NG.check_itr << endl;
		// // }

		// /*==============================================================================*/
		// /*Collecting scattered Phi variable from all processors*/
		// Vec temp_phi_seq;
		// VecScatter ctx_phi;
		// PetscScalar *_p;

		// NG.ierr = VecScatterCreateToAll(NG.temp_phi, &ctx_phi, &temp_phi_seq); CHKERRQ(NG.ierr);
		// NG.ierr = VecScatterBegin(ctx_phi, NG.temp_phi, temp_phi_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
		// NG.ierr = VecScatterEnd(ctx_phi, NG.temp_phi, temp_phi_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
		// NG.ierr = VecGetArray(temp_phi_seq, &_p);

		// for (int i = 0; i < NG.phi.size(); i++) {
		// 	NG.phi[i] = PetscRealPart(_p[i]);
		// }

		// NG.ierr = VecRestoreArray(temp_phi_seq, &_p); CHKERRQ(NG.ierr);
		// NG.ierr = VecScatterDestroy(&ctx_phi); CHKERRQ(NG.ierr);
		// NG.ierr = VecDestroy(&temp_phi_seq); CHKERRQ(NG.ierr);

		// toc(t_phi);
		// tic();
		// t_total += t_phi;
		
		/*==============================================================================*/
		/*Synaptogenesis and Tubulin equation*/
		/*Build Linear System*/
		NG.prepareSourceSum();
		NG.sum_grad_phi0_global = 0;
		NG.ierr = MPI_Allreduce(&NG.sum_grad_phi0_local, &NG.sum_grad_phi0_global, 1, MPI_FLOAT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
		NG.prepareTerm_source();
		if (NG.judge_syn == 0) {
			NG.ierr = MatZeroEntries(NG.GK_syn); CHKERRQ(NG.ierr);
		}
		NG.ierr = VecSet(NG.GR_syn, 0); CHKERRQ(NG.ierr);
		NG.ierr = MatZeroEntries(NG.GK_tub); CHKERRQ(NG.ierr);
		NG.ierr = VecSet(NG.GR_tub, 0); CHKERRQ(NG.ierr);
		NG.BuildLinearSystemProcessNG_syn_tub(tmesh, cpts); // build both Synaptogenesis and tubulin in the same loop
		NG.ierr = VecAssemblyEnd(NG.GR_syn); CHKERRQ(NG.ierr);
		if (NG.judge_syn == 0) {
			NG.ierr = MatAssemblyEnd(NG.GK_syn, MAT_FINAL_ASSEMBLY); CHKERRQ(NG.ierr);
		}
		NG.ierr = VecAssemblyEnd(NG.GR_tub); CHKERRQ(NG.ierr);
		NG.ierr = MatAssemblyEnd(NG.GK_tub, MAT_FINAL_ASSEMBLY); CHKERRQ(NG.ierr);

		/*Petsc solver setting for Synaptogenesis*/
		if (NG.judge_syn == 0) {
			NG.ierr = KSPCreate(PETSC_COMM_WORLD, &NG.ksp_syn); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetOperators(NG.ksp_syn, NG.GK_syn, NG.GK_syn); CHKERRQ(NG.ierr);
			NG.ierr = KSPGetPC(NG.ksp_syn, &NG.pc_syn); CHKERRQ(NG.ierr);
			NG.ierr = PCSetType(NG.pc_syn, PCBJACOBI); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetType(NG.ksp_syn, KSPGMRES); CHKERRQ(NG.ierr);
			NG.ierr = KSPGMRESSetRestart(NG.ksp_syn, 100); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetInitialGuessNonzero(NG.ksp_syn, PETSC_TRUE); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetTolerances(NG.ksp_syn, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 100000); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetFromOptions(NG.ksp_syn); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetUp(NG.ksp_syn); CHKERRQ(NG.ierr);
			if (NG.n == 0) {
				NG.ierr = KSPView(NG.ksp_syn, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(NG.ierr);
			}
			NG.judge_syn = 1;
		}

		NG.ierr = KSPSolve(NG.ksp_syn, NG.GR_syn, NG.temp_syn); CHKERRQ(NG.ierr);
		PetscInt its_syn;
		NG.ierr = KSPGetIterationNumber(NG.ksp_syn, &its_syn); CHKERRQ(NG.ierr);
		KSPConvergedReason reason_syn;
		NG.ierr = KSPGetConvergedReason(NG.ksp_syn, &reason_syn); CHKERRQ(NG.ierr);
		if (reason_syn < 0) {
			PetscPrintf(PETSC_COMM_WORLD, "KSP Synaptogenesis not converging  | KSPreason:  %d | KSPiter: %d \n", reason_syn, its_syn); CHKERRQ(NG.ierr);	
			return 3;
		}

		/*Petsc solver setting for Tubulin*/
		if (NG.judge_tub == 0) {
			NG.ierr = KSPCreate(PETSC_COMM_WORLD, &NG.ksp_tub); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetOperators(NG.ksp_tub, NG.GK_tub, NG.GK_tub); CHKERRQ(NG.ierr);
			NG.ierr = KSPGetPC(NG.ksp_tub, &NG.pc_tub); CHKERRQ(NG.ierr);
			NG.ierr = PCSetType(NG.pc_tub, PCBJACOBI); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetType(NG.ksp_tub, KSPGMRES); CHKERRQ(NG.ierr);
			NG.ierr = KSPGMRESSetRestart(NG.ksp_tub, 100); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetInitialGuessNonzero(NG.ksp_tub, PETSC_TRUE); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetTolerances(NG.ksp_tub, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 100000); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetFromOptions(NG.ksp_tub); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetUp(NG.ksp_tub); CHKERRQ(NG.ierr);
			if (NG.n == 0) {
				NG.ierr = KSPView(NG.ksp_tub, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(NG.ierr);
			}
			NG.judge_tub = 1;
		}

		NG.ierr = KSPSolve(NG.ksp_tub, NG.GR_tub, NG.temp_tub); CHKERRQ(NG.ierr);
		PetscInt its_tub;
		NG.ierr = KSPGetIterationNumber(NG.ksp_tub, &its_tub); CHKERRQ(NG.ierr);
		KSPConvergedReason reason_tub;
		NG.ierr = KSPGetConvergedReason(NG.ksp_tub, &reason_tub); CHKERRQ(NG.ierr);
		if (reason_tub < 0) {
			PetscPrintf(PETSC_COMM_WORLD, "KSP tubulin not converging ..... | KSPreason:  %d | KSPiter: %d \n", reason_tub, its_tub); CHKERRQ(NG.ierr);	
			return 3;	
		}

		/*==============================================================================*/
		/*Collecting scattered Synaptogenesis and Tubulin variables from all processors*/
		Vec temp_syn_seq, temp_tub_seq;
		VecScatter ctx_syn, ctx_tub;
		PetscScalar *_s, *_t;

		NG.ierr = VecScatterCreateToAll(NG.temp_syn, &ctx_syn, &temp_syn_seq); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterBegin(ctx_syn, NG.temp_syn, temp_syn_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterEnd(ctx_syn, NG.temp_syn, temp_syn_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
		NG.ierr = VecGetArray(temp_syn_seq, &_s); CHKERRQ(NG.ierr);

		NG.ierr = VecScatterCreateToAll(NG.temp_tub, &ctx_tub, &temp_tub_seq); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterBegin(ctx_tub, NG.temp_tub, temp_tub_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterEnd(ctx_tub, NG.temp_tub, temp_tub_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
		NG.ierr = VecGetArray(temp_tub_seq, &_t); CHKERRQ(NG.ierr);

		for (size_t i = 0; i < NG.syn.size(); i++) {
			NG.syn[i] = PetscRealPart(_s[i]);
			NG.tub[i] = PetscRealPart(_t[i]) * round(NG.phi[i]);
		}

		NG.ierr = VecRestoreArray(temp_syn_seq, &_s); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterDestroy(&ctx_syn); CHKERRQ(NG.ierr);
		NG.ierr = VecDestroy(&temp_syn_seq); CHKERRQ(NG.ierr);

		NG.ierr = VecRestoreArray(temp_tub_seq, &_t); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterDestroy(&ctx_tub); CHKERRQ(NG.ierr);
		NG.ierr = VecDestroy(&temp_tub_seq); CHKERRQ(NG.ierr);

		toc(t_synTub);
		tic();
		t_total += t_synTub;

		/*==============================================================================*/
		// Obtain initial local refinement information, the very first 5 iterations are purely used 
		// for getting diffused interface for applying local refinements (phi initialization is binary) 
		if ((NG.n == 5) && (localRefine == false)) {
			NG.PrintOutNeurons(neurons, NX_fine, NY_fine);
			NGvars.clear(); NGvars.resize(7);
			NGvars[0] = NG.phi;
			NGvars[1] = NG.syn;
			NGvars[2] = NG.tub;
			NGvars[3] = NG.theta_fine;
			NGvars[4] = NG.phi_0;
			NGvars[5] = NG.tub_0;
			NGvars[6] = NG.InterpolateVars_coarseKDtree(ConvertTo1DFloatVector(NG.prev_id), cpts_fine, kdTree_fine, cpts, 0, 0);

			CleanUpSolvers(NG); // Destroy solvers
			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
			iter = 0; // reset to 0 (beginning of the simulation)
			localRefine = true;
			return 2;
		}
		
		/*==============================================================================*/
		// Command line update
		PetscPrintf(PETSC_COMM_WORLD, "Step:%d/%d | Phi:[%d] %.3fs | Syn:%d[%d] Tub:%d[%d] %.3fs | Tip:%.3fs | DOFs:%d |\n",\
		NG.n, NG.end_iter, NR_itr, t_phi, reason_syn, its_syn, reason_tub, its_tub, t_synTub, t_tip, NG.phi.size()); CHKERRQ(NG.ierr);
		// PetscPrintf(PETSC_COMM_WORLD, "Step:%d/%d | Phi:%d[%d] %.3fs | Syn:%d[%d] Tub:%d[%d] %.3fs | Tip:%.3fs | DOFs:%d |\n",\
		// NG.n, NG.end_iter, reason_phi, its_phi, t_phi, reason_syn, its_syn, reason_tub, its_tub, t_synTub, t_tip, NG.phi.size()); CHKERRQ(NG.ierr);

		/*==============================================================================*/
		// Writing results to files
		if ((NG.n % NG.var_save_invl == 0) && (NG.n != 0)) {	
			NG.VisualizeVTK_PhysicalDomain_All(NG.n, path_out);
			
			// Make sure to have more than 6 threads
			string varName;	
			if (NG.comRank == 0) {			
				NG.PrintOutNeurons(neurons, 2*NX+1, 2*NY+1);
				// NG.PrintOutNeurons(neurons, localMaximaMatrix, 2*NX+1, 2*NY+1);
			// } else if (NG.comRank == 1) {			
			// 	vector<float> geoDist_all(cpts_fine.size(), 0);
			// 	for (size_t i = 0; i < distances.size(); i++) {
			// 		for (size_t j = 0; j < geodist[i].size(); j++) {
			// 			geoDist_all[j] = max(geoDist_all[j], geodist[i][j]);
			// 		}
			// 	}
			// 	varName = "geoDist_running_";
			// 	NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, NG.n, path_out, geoDist_all, varName); // solution on control points
			} else if (NG.comRank == 2) {
				varName = "localMax_running_";
				NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, NG.n, path_out, localMaximaMatrix, varName); // solution on control points
			} else if (NG.comRank == 3) {
				varName = "theta_running_"; // 0 since only theta is constant
				NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, 0, path_out, NG.theta_fine, varName); // solution on control points
			} else if (NG.comRank == 4) {
				varName = "tip_running_";
				NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, NG.n, path_out, tip, varName); // solution on control points
			// } else if (NG.comRank == 5) {
			// 	varName = "id_running_";
			// 	NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, NG.n, path_out, id, varName); // solution on control points
			} else {}

			toc(t_write);
			tic();
			t_total += t_write;

			// NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
			PetscPrintf(PETSC_COMM_WORLD, "| Wrote Physical Domain! %fs | Average time %fs | Total time: %f |\n", 
				t_write, t_total/NG.n, t_total); CHKERRQ(NG.ierr);
			PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------\n");
		}

		/*==============================================================================*/
		// Domain expansion and variable passing - back to main.cpp
		if ((NG.n % NG.expandCK_invl == 0) && (NG.n != 0)) {
			int expd_dir_local = NG.CheckExpansion(NG.phi, NX, NY);
			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
			int expd_dir_global = 4; // 4 - No expansion
			NG.ierr = MPI_Allreduce(&expd_dir_local, &expd_dir_global, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
			int expd_dir = expd_dir_global/NG.comSize; // approximate choice, different thread could have different choice 

			if (expd_dir <= 3) {
				switch (expd_dir) {
					case 0: // left
						NX += 3;	originX -= 6;		break;
					case 1: // top
						NY += 3;				break;
					case 2: // right
						NX += 3;				break;
					case 3: // bottom
						NY += 3;	originY -= 6;		break;
				}
			}

			NGvars.clear(); NGvars.resize(7);
			NGvars[0] = NG.phi;	
			NGvars[1] = NG.syn;
			NGvars[2] = NG.tub;
			NGvars[3] = NG.theta_fine;
			NGvars[4] = NG.phi_0;
			NGvars[5] = NG.tub_0;
			NGvars[6] = ConvertTo1DFloatVector(NG.prev_id);

			CleanUpSolvers(NG); // Destroy solvers
			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
			iter += 1;

			return 2;
		}
		iter += 1;
	}

	/*==============================================================================*/
	CleanUpSolvers(NG); // Destroy solvers
	NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);

	return 0; // ending simulation
}