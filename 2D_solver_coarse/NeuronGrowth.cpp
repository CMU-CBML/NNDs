// #include "NeuronGrowth.h"
// #include "utils.h"
// #include <cmath>
// #include <fstream>
// #include <iostream>
// #include <iomanip>
// #include <algorithm>
// #include <numeric>
// #include <stack>	// for tic toc
// #include <ctime>	// for tic toc
// #include <time.h>	// for srand
// #include <queue>	// for geodesic distance

// #include <cfloat>	// for FLT_MAX
// #include <limits>	// for numeric_limits

// #include <nanoflann.hpp> // for nearest points search with KD-tree
// // https://github.com/jlblancoc/nanoflann
// // sudo apt install libnanoflann-dev

// using namespace std;
// typedef unsigned int uint;
// const float PI = 3.1415926;

// std::stack<clock_t> tictoc_stack;
// float t_phi, t_synTub, t_syn, t_tub, t_tip, t_write(0), t_total(0);

// const int INF = 1e9;

// void tic() 
// {
// 	tictoc_stack.push(clock());
// }
// void toc(float &t) 
// {
// 	t = ((float)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC;
// 	tictoc_stack.pop();
// }

// float MatrixDet(float dxdt[2][2])
// {
// 	float det = dxdt[0][0] * dxdt[1][1] - dxdt[0][1] * dxdt[1][0];
// 	return det;
// }

// void Matrix2DInverse(float dxdt[2][2], float dtdx[2][2])
// {
// 	float det = MatrixDet(dxdt);
// 	dtdx[0][0] = 1.0 / det * (dxdt[1][1]);
// 	dtdx[0][1] = 1.0 / det * (-dxdt[0][1]);
// 	dtdx[1][0] = 1.0 / det * (-dxdt[1][0]);
// 	dtdx[1][1] = 1.0 / det * (dxdt[0][0]);
// }


// // // Dataset adaptor for nanoflann
// // struct Vertex2DCloud {
// //     const std::vector<Vertex2D>& pts;

// //     Vertex2DCloud(const std::vector<Vertex2D>& pts) : pts(pts) {}

// //     // Must return the number of data points
// //     inline size_t kdtree_get_point_count() const { return pts.size(); }

// //     // Returns the dim'th component of the idx'th point in the class
// //     inline float kdtree_get_pt(const size_t idx, const size_t dim) const {
// //         return dim == 0 ? pts[idx].coor[0] : pts[idx].coor[1];
// //     }

// //     // Optional: bounding-box computation
// //     template <class BBOX>
// //     bool kdtree_get_bbox(BBOX&) const { return false; }
// // };

// // using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
// //     nanoflann::L2_Simple_Adaptor<float, Vertex2DCloud>,
// //     Vertex2DCloud, 2 /* dim */>;
    

// NeuronGrowth::NeuronGrowth(){
// 	comm = PETSC_COMM_WORLD; //MPI_COMM_WORLD;
// 	mpiErr = MPI_Comm_rank(comm, &comRank);
// 	mpiErr = MPI_Comm_size(comm, &comSize);
// 	nProcess = comSize;
// 	judge_phi = 0;
// 	judge_syn = 0;
// 	judge_tub = 0;
// 	sum_grad_phi0_local = 0;
// 	sum_grad_phi0_global = 0;

// 	// integer variable setup
// 	expandCK_invl		= 100; 		// var_save_invl
// 	var_save_invl		= 100; 		// var_save_invl
// 	numNeuron 		= 1;	     	// numNeuron
// 	gc_sz			= 2;	     	// gc_sz
// 	aniso 			= 6;   		// aniso
// 	gamma 			= 10;  		// gamma
// 	// seed_radius 		= 15;		// seed radius
// 	seed_radius 		= 10;		// seed radius

// 	// variable setup
// 	kappa			= 2.0;		// kappa;
// 	dt			= 1e-3;		// time step
// 	Dc			= 3;		// syn D
// 	alpha			= 0.9;		// alpha
// 	alphaOverPi		= alpha / PI; 	// alphOverPix
// 	M_phi			= 60;		// M_phi
// 	s_coeff			= 0.007;	// s_coeff
// 	delta			= 0.20;		// delta
// 	epsilonb		= 0.04;		// epsilonb
// 	r			= 5;		// r
// 	g			= 0.1;		// g
// 	alphaT 			= 0.001;			// alpha_t
// 	betaT			= 0.001;	// beta_t
// 	Diff			= 4;		// Diff
// 	source_coeff		= 15;		// source_coeff
// }

// void NeuronGrowth::AssignProcessor(vector<vector<int>> &ele_proc)
// {
// 	ele_process.clear();
// 	for (int i = 0; i < ele_proc[comRank].size(); i++)
// 		ele_process.push_back(ele_proc[comRank][i]);
// }

// void NeuronGrowth::SetVariables(string fn_par)
// {
// 	string fname(fn_par), stmp;
// 	stringstream ss;
// 	ifstream fin;
// 	fin.open(fname);
// 	if (fin.is_open()) {
// 		fin >> stmp >> var_save_invl;
// 		fin >> stmp >> expandCK_invl;
// 		fin >> stmp >> gc_sz;
// 		fin >> stmp >> aniso;
// 		fin >> stmp >> gamma;
// 		fin >> stmp >> seed_radius;
// 		fin >> stmp >> kappa;
// 		fin >> stmp >> dt;
// 		fin >> stmp >> Dc;
// 		fin >> stmp >> alpha;
// 		alphaOverPi = alpha / PI; 	// alphOverPi
// 		fin >> stmp >> M_phi;
// 		fin >> stmp >> s_coeff;
// 		fin >> stmp >> delta;
// 		fin >> stmp >> epsilonb;
// 		fin >> stmp >> r;
// 		fin >> stmp >> g;
// 		fin >> stmp >> alphaT;
// 		fin >> stmp >> betaT;
// 		fin >> stmp >> Diff;
// 		fin >> stmp >> source_coeff;

// 		fin.close();
// 		PetscPrintf(PETSC_COMM_WORLD, "Parameter Loaded!\n");
// 	} else {
// 		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
// 	}
// }

// void NeuronGrowth::InitializeProblemNG(const int n_bz, vector<Vertex2D>& cpts, vector<Vertex2D> prev_cpts, const KDTree& kdTree_prev, vector<vector<float>> &NGvars, vector<array<float, 2>> &seed)
// {
// 	MPI_Barrier(PETSC_COMM_WORLD);
// 	PetscPrintf(PETSC_COMM_WORLD, "Initializing-----------------------------------------------------------------\n");

// 	/*Initialize parameters*/
// 	int cpt_sz = cpts.size();
// 	GaussInfo(3);
// 	n_bzmesh = n_bz;
// 	N_0.resize(cpt_sz);

// 	// float max_x, min_x, max_y, min_y;
// 	for (int i = 0; i<cpts.size(); i++) {
// 		if (n==0) {
// 			max_x = cpts[i].coor[0];
// 			min_x = cpts[i].coor[0];
// 			max_y = cpts[i].coor[1];
// 			min_y = cpts[i].coor[1];

// 			if (i < prev_cpts.size()) {
// 				prev_max_x = prev_cpts[i].coor[0];
// 				prev_min_x = prev_cpts[i].coor[0];
// 				prev_max_y = prev_cpts[i].coor[1];
// 				prev_min_y = prev_cpts[i].coor[1];
// 			}

// 		} else {
// 			max_x = max(cpts[i].coor[0],max_x);
// 			min_x = min(cpts[i].coor[0],min_x);
// 			max_y = max(cpts[i].coor[1],max_y);
// 			min_y = min(cpts[i].coor[1],min_y);

// 			if (i < prev_cpts.size()) {
// 				prev_max_x = max(prev_cpts[i].coor[0],prev_max_x);
// 				prev_min_x = min(prev_cpts[i].coor[0],prev_min_x);
// 				prev_max_y = max(prev_cpts[i].coor[1],prev_max_y);
// 				prev_min_y = min(prev_cpts[i].coor[1],prev_min_y);
// 			}
// 		}
// 	}

// 	PetscPrintf(PETSC_COMM_WORLD, "Checking maximum x and y value-----------------------------------------------\n");
// 	PetscPrintf(PETSC_COMM_WORLD, "max x: %f min x: %f | max y: %f min y %f\n", max_x, min_x, max_y, min_y);

// 	float r, x, y; // radius, x, y
// 	int nx = sqrt(cpts.size()); // only used at the beginning (square mesh)
// 	float dx = max_x/(float)nx; // only used at the beginning (square mesh)
// 	dx = 1;
// 	srand(time(NULL)); // random seed

// 	if (n == 0) { // Initialize
// 		for (int i = 0; i<cpts.size(); i++) {
// 			phi.push_back(0.f);
// 			tub.push_back(0.f);
// 			theta.push_back((float)(rand()%100)/(float)100); // orientation theta
// 			syn.push_back(0.f); // synaptogenesis
// 			Mphi.push_back(M_phi);
// 		}
// 		for (int i = 0; i<cpts.size(); i++) {
// 			x = cpts[i].coor[0];
// 			y = cpts[i].coor[1];

// 			// calculate new boundary label
// 			if (x==min_x || x==max_x || y==min_y || y==max_y) {
// 				cpts[i].label = 1;
// 			} else {
// 				cpts[i].label = 0;
// 			}   

// 			for (int j = 0; j<seed.size(); j++) { // multiple neurons
// 				r = sqrt((x-(float)seed[j][0])*(x-(float)seed[j][0])+(y-(float)seed[j][1])*(y-(float)seed[j][1])); // radius of initial soma
// 				if (r < (seed_radius*dx)) { // initial soma
// 					// std::cout << seed_radius << " " << dx << " " << seed_radius*dx << std::endl;
// 					phi[i] = 1.0f;
// 					tub[i] = (0.5+0.5*tanh((sqrt(seed_radius*dx)-r)/2)); // based eqn in literature
// 				} 
// 			} 
// 		}
// 		phi_0 = phi; // initial phi
// 		tub_0 = tub; // initial tub
// 	} else { // Collect from NGvars - continue simulation
// 		// phi.clear(); 	phi.resize(cpts.size());
// 		// syn.clear(); 	syn.resize(cpts.size());
// 		// tub.clear(); 	tub.resize(cpts.size());
// 		// theta.clear();	theta.resize(cpts.size());
// 		// phi_0.clear(); 	phi_0.resize(cpts.size());
// 		// tub_0.clear(); 	tub_0.resize(cpts.size());
// 		// for (int i = 0; i < prev_cpts.size(); i++) {
// 		// 	if (abs(remainder(prev_cpts[i].coor[0],2)) != 1)
// 		// 		prev_cpts[i].coor[0] = round(prev_cpts[i].coor[0]);
// 		// 	if (abs(remainder(prev_cpts[i].coor[1],2)) != 1)
// 		// 		prev_cpts[i].coor[1] = round(prev_cpts[i].coor[1]);
// 		// }

// 		// for (int i = 0; i < cpts.size(); i++) {
// 		// 	x = cpts[i].coor[0];
// 		// 	if (abs(remainder(x,2)) != 1) {
// 		// 		x = round(x);
// 		// 	}
// 		// 	y = cpts[i].coor[1];
// 		// 	if (abs(remainder(y,2)) != 1) {
// 		// 		y = round(y);
// 		// 	}

// 		// 	// calculate new boundary label
// 		// 	if (x==min_x || x==max_x || y==min_y || y==max_y) {
// 		// 		cpts[i].label = 1;
// 		// 	} else {
// 		// 		cpts[i].label = 0;
// 		// 	}   

// 		// 	int ind, indDown, indUp, indLeft, indRight;
// 		// 	if (SearchPair(prev_cpts, round5(x), round5(y), ind)) {
// 		// 	// if (SearchPair(prev_cpts, cpts[i].coor[0], cpts[i].coor[1], ind)) {
// 		// 		phi[i] = NGvars[0][ind];
// 		// 		syn[i] = NGvars[1][ind];
// 		// 		tub[i] = NGvars[2][ind];
// 		// 		theta[i] = NGvars[3][ind];
// 		// 		phi_0[i] = NGvars[4][ind];
// 		// 		tub_0[i] = NGvars[5][ind];
// 		// 	} else if ((abs(remainder(x,2)) == 1) && (abs(remainder(y,2)) != 1) 
// 		// 		&& SearchPair(prev_cpts, floorf(x)-1, round5(y), indDown)
// 		// 		&& SearchPair(prev_cpts, floorf(x)+1, round5(y), indUp)) {
// 		// 			phi[i] = (NGvars[0][indDown] + NGvars[0][indUp])/2;
// 		// 			syn[i] = (NGvars[1][indDown] + NGvars[1][indUp])/2;
// 		// 			tub[i] = (NGvars[2][indDown] + NGvars[2][indUp])/2;
// 		// 			// theta[i] = max(NGvars[3][indDown], NGvars[3][indUp]);
// 		// 			theta[i] = (float)(rand()%100)/(float)100;
// 		// 			phi_0[i] = (NGvars[4][indDown] + NGvars[4][indUp])/2;
// 		// 			tub_0[i] = (NGvars[5][indDown] + NGvars[5][indUp])/2;
// 		// 	} else if ((abs(remainder(x,2)) != 1) && (abs(remainder(y,2)) == 1)
// 		// 		&& SearchPair(prev_cpts, round5(x), floorf(y)-1, indDown)
// 		// 		&& SearchPair(prev_cpts, round5(x), floorf(y)+1, indUp)){
// 		// 			phi[i] = (NGvars[0][indDown] + NGvars[0][indUp])/2;
// 		// 			syn[i] = (NGvars[1][indDown] + NGvars[1][indUp])/2;
// 		// 			tub[i] = (NGvars[2][indDown] + NGvars[2][indUp])/2;
// 		// 			// theta[i] = max(NGvars[3][indDown], NGvars[3][indUp]);
// 		// 			theta[i] = (float)(rand()%100)/(float)100;
// 		// 			phi_0[i] = (NGvars[4][indDown] + NGvars[4][indUp])/2;
// 		// 			tub_0[i] = (NGvars[5][indDown] + NGvars[5][indUp])/2;
// 		// 	} else if ((abs(remainder(x,2)) == 1) && (abs(remainder(y,2)) == 1)
// 		// 		&& SearchPair(prev_cpts, floorf(x)-1, floorf(y)-1, indDown)
// 		// 		&& SearchPair(prev_cpts, floorf(x)+1, floorf(y)+1, indUp)
// 		// 		&& SearchPair(prev_cpts, floorf(x)-1, floorf(y)+1, indLeft)
// 		// 		&& SearchPair(prev_cpts, floorf(x)+1, floorf(y)-1, indRight)) {
// 		// 			phi[i] = (NGvars[0][indDown] + NGvars[0][indUp] + NGvars[0][indLeft] + NGvars[0][indRight])/4;
// 		// 			syn[i] = (NGvars[1][indDown] + NGvars[1][indUp] + NGvars[1][indLeft] + NGvars[1][indRight])/4;
// 		// 			tub[i] = (NGvars[2][indDown] + NGvars[2][indUp] + NGvars[2][indLeft] + NGvars[2][indRight])/4;
// 		// 			// theta[i] = max(max(NGvars[3][indDown], NGvars[3][indUp]), max(NGvars[3][indLeft], NGvars[3][indRight]));
// 		// 			theta[i] = (float)(rand()%100)/(float)100;
// 		// 			phi_0[i] = (NGvars[4][indDown] + NGvars[4][indUp] + NGvars[4][indLeft] + NGvars[4][indRight])/4;
// 		// 			tub_0[i] = (NGvars[5][indDown] + NGvars[5][indUp] + NGvars[5][indLeft] + NGvars[5][indRight])/4;
// 		// 	} else {
// 		// 		phi[i] = 0;
// 		// 		syn[i] = 0;
// 		// 		tub[i] = 0;
// 		// 		theta[i] = (float)(rand()%100)/(float)100;
// 		// 		phi_0[i] = 0;
// 		// 		tub_0[i] = 0;
// 		// 	}		
// 		// 	Mphi.push_back(M_phi);
// 		// }

// 		phi.clear(); 	phi.resize(cpts.size());
// 		syn.clear(); 	syn.resize(cpts.size());
// 		tub.clear(); 	tub.resize(cpts.size());
// 		theta.clear();	theta.resize(cpts.size());
// 		phi_0.clear(); 	phi_0.resize(cpts.size());
// 		tub_0.clear(); 	tub_0.resize(cpts.size());

// 		// phi = InterpolateVars_coarse(NGvars[0], prev_cpts, cpts, 1);
// 		// syn = InterpolateVars_coarse(NGvars[1], prev_cpts, cpts, 1);
// 		// tub = InterpolateVars_coarse(NGvars[2], prev_cpts, cpts, 1);
// 		// theta = InterpolateVars_coarse(NGvars[3], prev_cpts, cpts, 1);
// 		// phi_0 = InterpolateVars_coarse(NGvars[4], prev_cpts, cpts, 1);
// 		// tub_0 = InterpolateVars_coarse(NGvars[5], prev_cpts, cpts, 1);

// 		phi = InterpolateVars_coarse1(NGvars[0], prev_cpts, kdTree_prev, cpts, 1, 0);
// 		syn = InterpolateVars_coarse1(NGvars[1], prev_cpts, kdTree_prev, cpts, 1, 0);
// 		tub = InterpolateVars_coarse1(NGvars[2], prev_cpts, kdTree_prev, cpts, 1, 0);
// 		theta = InterpolateVars_coarse1(NGvars[3], prev_cpts, kdTree_prev, cpts, 1, 1); // isTheta = 1 for special init
// 		phi_0 = InterpolateVars_coarse1(NGvars[4], prev_cpts, kdTree_prev, cpts, 1, 0);
// 		tub_0 = InterpolateVars_coarse1(NGvars[5], prev_cpts, kdTree_prev, cpts, 1, 0);

// 		for (int i = 0; i < cpts.size(); i++) {	
// 			Mphi.push_back(M_phi);
// 		}

// 	}

// 	this->cpts = cpts; // for passing into formFunction

// 	// // writing initialized variables for debugging purposes
// 	// bool visualization = true;
// 	// PrintVec2TXT(phi, "/home/kuanrenqian/Research/3D_NeuronGrowth/io/2DNG/var_check/phi_inital.txt", visualization);
// 	// PrintVec2TXT(tub, "/home/kuanrenqian/Research/3D_NeuronGrowth/io/2DNG/var_check/tub_inital.txt", visualization);
// 	// PrintVec2TXT(theta, "/home/kuanrenqian/Research/3D_NeuronGrowth/io/2DNG/var_check/theta_inital.txt", visualization);
// 	// PrintVec2TXT(syn, "/home/kuanrenqian/Research/3D_NeuronGrowth/io/2DNG/var_check/syn_inital.txt", visualization);

// 	/*Initialize petsc vector, matrix*/
// 	PetscInt mat_dim = cpt_sz;

// 	// 250 = 2*5^3. (2 var, 5 basis, 3D) || 25 = 1*5^2. (1 var, 5 basis, 2D)	

// 	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_phi);
	
// 	ierr = MatCreate(PETSC_COMM_WORLD, &J);
// 	ierr = MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
// 	ierr = MatSetType(J, MATMPIAIJ);
// 	ierr = MatMPIAIJSetPreallocation(J, 25, NULL, 25, NULL);
// 	ierr = MatSetOption(J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
// 	ierr = MatSetUp(J);
// 	// MatGetOwnershipRange(J, NULL, NULL);

// 	ierr = MatCreate(PETSC_COMM_WORLD, &GK_phi);
// 	ierr = MatSetSizes(GK_phi, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
// 	ierr = MatSetType(GK_phi, MATMPIAIJ);
// 	ierr = MatMPIAIJSetPreallocation(GK_phi, 25, NULL, 25, NULL);
// 	ierr = MatSetOption(GK_phi, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
// 	ierr = MatSetUp(GK_phi);
// 	// MatGetOwnershipRange(GK_phi, NULL, NULL);
// 	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &GR_phi);
// 	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_phi);

// 	ierr = MatCreate(PETSC_COMM_WORLD, &GK_syn);
// 	ierr = MatSetSizes(GK_syn, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
// 	ierr = MatSetType(GK_syn, MATMPIAIJ);
// 	ierr = MatMPIAIJSetPreallocation(GK_syn, 25, NULL, 25, NULL);
// 	ierr = MatSetOption(GK_syn, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
// 	ierr = MatSetUp(GK_syn);
// 	// MatGetOwnershipRange(GK_syn, NULL, NULL);
// 	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &GR_syn);
// 	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_syn);

// 	ierr = MatCreate(PETSC_COMM_WORLD, &GK_tub);
// 	ierr = MatSetSizes(GK_tub, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
// 	ierr = MatSetType(GK_tub, MATMPIAIJ);
// 	ierr = MatMPIAIJSetPreallocation(GK_tub, 25, NULL, 25, NULL);
// 	ierr = MatSetOption(GK_tub, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
// 	ierr = MatSetUp(GK_tub);
// 	// MatGetOwnershipRange(GK_tub, NULL, NULL);
// 	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &GR_tub);
// 	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_tub);

// 	vars.reserve(3*sizeof(float));
// 	vars.resize(3);

// 	pre_EMatrixSolve.clear();
// 	pre_EVectorSolve.clear();
// }

// void NeuronGrowth::ToPETScVec(vector<float> input, Vec& petscVec)
// {
// 	PetscInt localStart, localEnd;
// 	VecGetOwnershipRange(petscVec, &localStart, &localEnd);
// 	for (int i = localStart; i<localEnd; i++) {
// 		PetscScalar value = input[i]; // Adjust the index to match the localValues
// 		VecSetValues(petscVec, 1, &i, &value, INSERT_VALUES);
// 	}
// 	VecAssemblyBegin(petscVec);
// 	VecAssemblyEnd(petscVec);
// }

// void NeuronGrowth::ReadBezierElementProcess(string fn)
// {
// 	string stmp;
// 	int npts, neles, nfunctions, itmp, itmp1;
// 	int add(0);

// 	string fname_cmat = fn + "cmat.txt";
// 	ifstream fin_cmat;
// 	fin_cmat.clear();
// 	fin_cmat.open(fname_cmat);

// 	if (fin_cmat.is_open())
// 	{
// 		fin_cmat >> neles;
// 		bzmesh_process.clear();
// 		bzmesh_process.resize(ele_process.size());
// 		for (int i = 0; i < neles; i++)
// 		{
// 			if ((i == ele_process[add]) && (add < ele_process.size()))
// 			{
// 				fin_cmat >> itmp >> nfunctions >> bzmesh_process[add].type;
// 				bzmesh_process[add].cmat.resize(nfunctions);
// 				bzmesh_process[add].IEN.resize(nfunctions);
// 				for (int j = 0; j < nfunctions; j++)
// 					fin_cmat >> bzmesh_process[add].IEN[j];
// 				for (int j = 0; j < nfunctions; j++)
// 				{
// 					for (int k = 0; k < 16; k++)
// 					{
// 						fin_cmat >> bzmesh_process[add].cmat[j][k];
// 					}
// 				}
// 				add++;
// 			}
// 			else
// 			{
// 				fin_cmat >> stmp >> nfunctions >> itmp;
// 				for (int j = 0; j < nfunctions; j++)
// 					fin_cmat >> stmp;
// 				for (int j = 0; j < nfunctions; j++)
// 					for (int k = 0; k < 16; k++)
// 						fin_cmat >> stmp;
// 			}
// 		}
// 		fin_cmat.close();
// 		// PetscPrintf(PETSC_COMM_WORLD, "Bezier Matrices Loaded!\n");
// 	}
// 	else
// 	{
// 		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_cmat.c_str());
// 	}

// 	string fname_bzpt = fn + "bzpt.txt";
// 	ifstream fin_bzpt;
// 	fin_bzpt.open(fname_bzpt);
// 	add = 0;
// 	float tmp_xuan;
// 	if (fin_bzpt.is_open())
// 	{
// 		fin_bzpt >> npts;
// 		getline(fin_bzpt, stmp);
// 		for (int e = 0; e < neles; e++)
// 		{
// 			if ((e == ele_process[add]) && (add < ele_process.size()))
// 			{
// 				bzmesh_process[add].pts.resize(bzpt_num);
// 				for (int i = 0; i < bzpt_num; i++)
// 				{
// 					fin_bzpt >> bzmesh_process[add].pts[i][0] >> bzmesh_process[add].pts[i][1] >> tmp_xuan;
// 				}
// 				add++;
// 			}
// 			else
// 			{
// 				for (int i = 0; i < bzpt_num; i++)
// 					fin_bzpt >> stmp >> stmp >> stmp;
// 			}

// 		}
// 		fin_bzpt.close();
// 		// PetscPrintf(PETSC_COMM_WORLD, "Bezier Points Loaded!\n");
// 	}
// 	else
// 	{
// 		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_bzpt.c_str());
// 	}
// }

// void NeuronGrowth::GaussInfo(int ng)
// {
// 	Gpt.clear();
// 	wght.clear();
// 	switch (ng)
// 	{
// 		case 2:
// 		{
// 			Gpt.resize(ng);
// 			wght.resize(ng);
// 			Gpt[0] = 0.2113248654051871;
// 			Gpt[1] = 0.7886751345948129;
// 			wght[0] = 1.;
// 			wght[1] = 1.;
// 			break;
// 		}
// 		case 3:
// 		{
// 			Gpt.resize(ng);
// 			wght.resize(ng);
// 			Gpt[0] = 0.1127016653792583;
// 			Gpt[1] = 0.5;
// 			Gpt[2] = 0.8872983346207417;
// 			wght[0] = 0.5555555555555556;
// 			wght[1] = 0.8888888888888889;
// 			wght[2] = 0.5555555555555556;
// 			break;
// 		}
// 		case 4:
// 		{
// 			Gpt.resize(ng);
// 			wght.resize(ng);
// 			Gpt[0] = 0.06943184420297371;
// 			Gpt[1] = 0.33000947820757187;
// 			Gpt[2] = 0.6699905217924281;
// 			Gpt[3] = 0.9305681557970262;
// 			wght[0] = 0.3478548451374539;
// 			wght[1] = 0.6521451548625461;
// 			wght[2] = 0.6521451548625461;
// 			wght[3] = 0.3478548451374539;
// 			break;
// 		}
// 		case 5:
// 		{
// 			Gpt.resize(ng);
// 			wght.resize(ng);
// 			Gpt[0] = 0.046910077030668;
// 			Gpt[1] = 0.2307653449471585;
// 			Gpt[2] = 0.5;
// 			Gpt[3] = 0.7692346550528415;
// 			Gpt[4] = 0.953089922969332;
// 			wght[0] = 0.2369268850561891;
// 			wght[1] = 0.4786286704993665;
// 			wght[2] = 0.5688888888888889;
// 			wght[3] = 0.4786286704993665;
// 			wght[4] = 0.2369268850561891;
// 			break;
// 		}
// 		default:
// 		{
// 			Gpt.resize(2);
// 			wght.resize(2);
// 			Gpt[0] = 0.2113248654051871;
// 			Gpt[1] = 0.7886751345948129;
// 			wght[0] = 1.;
// 			wght[1] = 1.;
// 			break;
// 		}
// 	}
// }

// void NeuronGrowth::BasisFunction(float u, float v, const vector<array<float, 2>> &pt, const vector<array<float, 16>> &cmat, vector<float> &Nx, vector<array<float, 2>> &dNdx, float dudx[2][2], float &detJ)
// {
// 	float Nu[4] = {(1. - u) * (1. - u) * (1. - u), 3. * (1. - u) * (1. - u) * u, 3. * (1. - u) * u * u, u * u * u};
// 	float Nv[4] = {(1. - v) * (1. - v) * (1. - v), 3. * (1. - v) * (1. - v) * v, 3. * (1. - v) * v * v, v * v * v};
// 	// float Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
// 	float dNdu[4] = {-3. * (1. - u) * (1. - u), 3. - 12. * u + 9. * u * u, 3. * (2. - 3. * u) * u, 3. * u * u};
// 	float dNdv[4] = {-3. * (1. - v) * (1. - v), 3. - 12. * v + 9. * v * v, 3. * (2. - 3. * v) * v, 3. * v * v};
// 	// float dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
// 	float dNdt[bzpt_num][2];
// 	float Nx_bz[bzpt_num];
// 	float dNdx_bz[bzpt_num][2];

// 	Nx.clear();
// 	dNdx.clear();
// 	Nx.resize(cmat.size(), 0);
// 	dNdx.resize(cmat.size(), {0, 0});

// 	int i, j, a, b, c, loc;
// 	loc = 0;
// 	for (i = 0; i < 4; i++)
// 	{
// 		for (j = 0; j < 4; j++)
// 		{
// 			// for (k = 0; k < 4; k++)	{
// 			Nx_bz[loc] = Nu[j] * Nv[i];	// * Nw[i];
// 			dNdt[loc][0] = dNdu[j] * Nv[i]; // * Nw[i];
// 			dNdt[loc][1] = Nu[j] * dNdv[i]; // * Nw[i];
// 			// dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
// 			loc++;
// 			// }
// 		}
// 	}

// 	float dxdt[2][2] = {{0}};
// 	for (loc = 0; loc < bzpt_num; loc++)
// 		for (a = 0; a < 2; a++)
// 			for (b = 0; b < 2; b++)
// 				dxdt[a][b] += pt[loc][a] * dNdt[loc][b];

// 	float dtdx[2][2] = {{0}};
// 	Matrix2DInverse(dxdt, dtdx);

// 	// 1st derivatives
// 	for (i = 0; i < bzpt_num; i++)
// 	{
// 		dNdx_bz[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0]; // + dNdt[i][2] * dtdx[2][0];
// 		dNdx_bz[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1]; // + dNdt[i][2] * dtdx[2][1];
// 										   // dNdx_bz[i][2] = dNdt[i][0] * dtdx[0][2] + dNdt[i][1] * dtdx[1][2] + dNdt[i][2] * dtdx[2][2];
// 	}
// 	for (int i = 0; i < 2; i++)
// 		for (int j = 0; j < 2; j++)
// 			dudx[i][j] = dtdx[i][j];

// 	detJ = MatrixDet(dxdt);
// 	detJ = 0.25 * detJ;

// 	for (i = 0; i < cmat.size(); i++)
// 	{
// 		for (j = 0; j < bzpt_num; j++)
// 		{
// 			Nx[i] += cmat[i][j] * Nx_bz[j];
// 			for (int m = 0; m < 2; m++)
// 			{
// 				dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
// 			}
// 		}
// 	}
// }

// void NeuronGrowth::BasisFunction(float u, float v, int nen, const vector<array<float, 2>> &pt, const vector<array<float, 16>> &cmat, vector<float> &Nx, vector<array<float, 2>> &dNdx, vector<array<array<float, 2>, 2>> &dN2dx2, float dudx[2][2], float &detJ)
// {
// 	float Nu[4] = {(1. - u) * (1. - u) * (1. - u), 3. * (1. - u) * (1. - u) * u, 3. * (1. - u) * u * u, u * u * u};
// 	float Nv[4] = {(1. - v) * (1. - v) * (1. - v), 3. * (1. - v) * (1. - v) * v, 3. * (1. - v) * v * v, v * v * v};
// 	// float Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };

// 	float dNdu[4] = {-3. * (1. - u) * (1. - u), 3. - 12. * u + 9. * u * u, 3. * (2. - 3. * u) * u, 3. * u * u};
// 	float dNdv[4] = {-3. * (1. - v) * (1. - v), 3. - 12. * v + 9. * v * v, 3. * (2. - 3. * v) * v, 3. * v * v};
// 	// float dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };

// 	float dN2du2[4] = {6. * (1. - u), -12. + 18. * u, 6. - 18. * u, 6. * u};
// 	float dN2dv2[4] = {6. * (1. - v), -12. + 18. * v, 6. - 18. * v, 6. * v};
// 	// float dN2dw2[4] = { 6.*(1. - w),-12. + 18.*w,6. - 18.*w,6.*w };

// 	float dNdt[bzpt_num][2];
// 	float dN2dt2[bzpt_num][2][2];
// 	float Nx_bz[bzpt_num];
// 	float dNdx_bz[bzpt_num][2];
// 	float dN2dx2_bz[bzpt_num][2][2];

// 	Nx.clear();
// 	dNdx.clear();
// 	dN2dx2.clear();
// 	Nx.resize(nen, 0);
// 	dNdx.resize(nen, {0});
// 	dN2dx2.resize(nen, {{0}});

// 	int i, j, k, a, b, c, loc;
// 	loc = 0;
// 	for (j = 0; j < 4; j++)
// 	{
// 		for (k = 0; k < 4; k++)
// 		{
// 			Nx_bz[loc] = Nu[k] * Nv[j];
// 			dNdt[loc][0] = dNdu[k] * Nv[j];
// 			dNdt[loc][1] = Nu[k] * dNdv[j];
// 			dN2dt2[loc][0][0] = dN2du2[k] * Nv[j];
// 			dN2dt2[loc][0][1] = dNdu[k] * dNdv[j];
// 			dN2dt2[loc][1][0] = dNdu[k] * dNdv[j];
// 			dN2dt2[loc][1][1] = Nu[k] * dN2dv2[j];
// 			loc++;
// 		}
// 	}

// 	float dxdt[2][2] = {{0}};
// 	for (loc = 0; loc < bzpt_num; loc++)
// 	{
// 		for (a = 0; a < 2; a++)
// 		{
// 			for (b = 0; b < 2; b++)
// 			{
// 				dxdt[a][b] += pt[loc][a] * dNdt[loc][b];
// 			}
// 		}
// 	}

// 	float dtdx[2][2] = {{0}};
// 	Matrix2DInverse(dxdt, dtdx);

// 	// 1st derivatives
// 	for (i = 0; i < bzpt_num; i++)
// 	{
// 		dNdx_bz[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0];
// 		dNdx_bz[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1];
// 	}
// 	for (int i = 0; i < 2; i++)
// 		for (int j = 0; j < 2; j++)
// 			dudx[i][j] = dtdx[i][j];

// 	detJ = MatrixDet(dxdt);
// 	detJ = 0.25 * detJ;

// 	// 2nd derivatives
// 	float dx2dt2[2][4] = {{0}};
// 	float dt2dx2[2][4] = {{0}};
// 	for (int l = 0; l < 2; l++)
// 		for (loc = 0; loc < bzpt_num; loc++)
// 			for (a = 0; a < 2; a++)
// 				for (b = 0; b < 2; b++)
// 					dx2dt2[l][2 * a + b] += pt[loc][l] * dN2dt2[loc][a][b];

// 	for (i = 0; i < 2; i++)
// 		for (j = 0; j < 2; j++)
// 			for (k = 0; k < 2; k++)
// 				for (a = 0; a < 2; a++)
// 					for (b = 0; b < 2; b++)
// 						for (c = 0; c < 2; c++)
// 							dt2dx2[c][2 * i + j] -= dx2dt2[k][2 * a + b] * dtdx[a][i] * dtdx[b][j] * dtdx[c][k];

// 	for (loc = 0; loc < bzpt_num; loc++)
// 		for (i = 0; i < 2; i++)
// 			for (j = 0; j < 2; j++)
// 				dN2dx2_bz[loc][i][j] = 0.;

// 	for (loc = 0; loc < bzpt_num; loc++)
// 	{
// 		for (i = 0; i < 2; i++)
// 		{
// 			for (j = 0; j < 2; j++)
// 			{
// 				for (a = 0; a < 2; a++)
// 				{
// 					for (b = 0; b < 2; b++)
// 					{
// 						dN2dx2_bz[loc][i][j] += dN2dt2[loc][a][b] * dtdx[a][i] * dtdx[b][j];
// 					}
// 					dN2dx2_bz[loc][i][j] += dNdt[loc][a] * dt2dx2[a][2 * i + j];
// 				}
// 			}
// 		}
// 	}

// 	for (i = 0; i < nen; i++)
// 	{
// 		for (j = 0; j < bzpt_num; j++)
// 		{
// 			Nx[i] += cmat[i][j] * Nx_bz[j];
// 			for (int m = 0; m < 2; m++)
// 			{
// 				dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
// 				for (int n = 0; n < 2; n++)
// 				{
// 					dN2dx2[i][m][n] += cmat[i][j] * dN2dx2_bz[j][m][n];
// 				}
// 			}
// 		}
// 	}
// }

// void NeuronGrowth::WeightingFunction(const float velocity[2], const float &s, const float &tau, const vector<float> &Nx, const vector<array<float, 2>> &dNdx, vector<float> &Wx)
// {
// 	/*Calculate weighting function using SUPG parameter Tau*/
// 	float U = sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] /* + velocity[2] * velocity[2] */);
// 	for (int loc = 0; loc < Nx.size(); loc++)
// 		Wx[loc] = Nx[loc] + tau * (velocity[0] * dNdx[loc][0] + velocity[1] * dNdx[loc][1] /* + velocity[2] * dNdx[loc][2] */);
// }

// void NeuronGrowth::ApplyBoundaryCondition(const float bc_value, int pt_num, int variable_num, vector<vector<float>> &EMatrixSolve, vector<float> &EVectorSolve)
// {
// 	int j, k;
// 	int nen = EVectorSolve.size();
// 	for (j = 0; j < nen; j++)
// 	{
// 		EVectorSolve[j] -= bc_value * EMatrixSolve[j][pt_num + variable_num * nen];
// 	}
// 	for (j = 0; j < nen; j++)
// 	{
// 		EMatrixSolve[j][pt_num + variable_num * nen] = 0.0;
// 		EMatrixSolve[pt_num + variable_num * nen][j] = 0.0;
// 	}
// 	EMatrixSolve[pt_num + variable_num * nen][pt_num + variable_num * nen] = 1.0;
// 	EVectorSolve[pt_num + variable_num * nen] = bc_value;
// }

// void NeuronGrowth::MatrixAssembly(vector<vector<float>> &EMatrixSolve, const vector<int> &IEN, Mat &GK)
// {
// 	int i, j, A, B, m, n;
// 	int row_start, row_end, row_now;
// 	int add = 0;
// 	int nen = IEN.size();

// 	PetscInt *nodeList = new PetscInt[nen * 1];
// 	PetscReal *tmpGK = new PetscReal[nen * 1 * nen * 1];

// 	for (m = 0; m < IEN.size(); m++)
// 	{
// 		A = IEN[m];
// 		nodeList[1 * m] = 1 * A;
// 		for (n = 0; n < IEN.size(); n++)
// 		{
// 				tmpGK[add] = EMatrixSolve[m][n];
// 				add++;
// 		}
// 	}
// 	MatSetValues(GK, nen * 1, nodeList, nen * 1, nodeList, tmpGK, ADD_VALUES);
// 	delete nodeList;
// 	delete tmpGK;
// }

// void NeuronGrowth::ResidualAssembly(vector<float> &EVectorSolve, const vector<int> &IEN, Vec &GR)
// {
// 	int i, m, A;
// 	int add = 0;
// 	int nen = IEN.size();

// 	PetscInt *nodeList = new PetscInt[nen];
// 	PetscReal *tmpGR = new PetscReal[nen];

// 	for (i = 0; i < IEN.size(); i++)
// 	{
// 		A = IEN[i];
// 		nodeList[i * 1 + 0] = A * 1 + 0;
// 		tmpGR[add] = EVectorSolve[i];
// 		add += 1;
// 	}
// 	VecSetValues(GR, nen * 1, nodeList, tmpGR, ADD_VALUES);
// 	delete nodeList;
// 	delete tmpGR;
// }

// void NeuronGrowth::MatrixAssembly_insert(vector<vector<float>> &EMatrixSolve, const vector<int> &IEN, Mat &GK)
// {
// 	int i, j, A, B, m, n;
// 	int row_start, row_end, row_now;
// 	int add = 0;
// 	int nen = IEN.size();

// 	PetscInt *nodeList = new PetscInt[nen * 1];
// 	PetscReal *tmpGK = new PetscReal[nen * 1 * nen * 1];

// 	for (m = 0; m < IEN.size(); m++)
// 	{
// 		A = IEN[m];
// 		nodeList[1 * m] = 1 * A;
// 		for (n = 0; n < IEN.size(); n++)
// 		{
// 				tmpGK[add] = EMatrixSolve[m][n];
// 				add++;
// 		}
// 	}
// 	MatSetValues(GK, nen * 1, nodeList, nen * 1, nodeList, tmpGK, INSERT_VALUES);
// 	delete nodeList;
// 	delete tmpGK;
// }

// void NeuronGrowth::ResidualAssembly_insert(vector<float> &EVectorSolve, const vector<int> &IEN, Vec &GR)
// {
// 	int i, m, A;
// 	int add = 0;
// 	int nen = IEN.size();

// 	PetscInt *nodeList = new PetscInt[nen];
// 	PetscReal *tmpGR = new PetscReal[nen];

// 	for (i = 0; i < IEN.size(); i++)
// 	{
// 		A = IEN[i];
// 		nodeList[i * 1 + 0] = A * 1 + 0;
// 		tmpGR[add] = EVectorSolve[i];
// 		add += 1;
// 	}
// 	VecSetValues(GR, nen * 1, nodeList, tmpGR, INSERT_VALUES);
// 	delete nodeList;
// 	delete tmpGR;
// }

// void NeuronGrowth::MatrixAssembly_2var(vector<vector<float>>& EMatrixSolve, const vector<int>& IEN, Mat& GK)
// {
// 	int i, j, A, B, m, n;
// 	int row_start, row_end, row_now;
// 	int add = 0;
// 	int nen = IEN.size();


// 	PetscInt *nodeList = new PetscInt[nen * 2];	
// 	PetscReal *tmpGK = new PetscReal[nen * 2 * nen * 2];
	
// 	for (m = 0; m<IEN.size(); m++){
// 		A = IEN[m];
// 		for (i = 0; i < 2; i++)	{
// 			nodeList[2 * m + i] = 2 * A + i;
// 			for (n = 0; n<IEN.size(); n++){
// 				for (j = 0; j < 2; j++)	{
// 					tmpGK[add] = EMatrixSolve[m + i*nen][n + j*nen];
// 					add++;
// 				}   
// 			}
// 		}
// 	}	
// 	MatSetValues(GK, nen * 2, nodeList, nen * 2, nodeList, tmpGK, ADD_VALUES);
// 	delete nodeList;
// 	delete tmpGK;
// }

// void NeuronGrowth::ResidualAssembly_2var(vector<float>& EVectorSolve, const vector<int>& IEN, Vec& GR)
// {
// 	int i, m, A;
// 	int add = 0;
// 	int nen = IEN.size();

// 	PetscInt *nodeList = new PetscInt[nen * 2];
// 	PetscReal *tmpGR = new PetscReal[nen * 2];
	
// 	for (i = 0; i<IEN.size(); i++){
// 		A = IEN[i];
// 		nodeList[i * 2 + 0] = A * 2 + 0;
// 		nodeList[i * 2 + 1] = A * 2 + 1;
// 		tmpGR[add] = EVectorSolve[i];
// 		tmpGR[add + 1] = EVectorSolve[i + nen];
// 		add+=2;
// 	}
// 	VecSetValues(GR, nen * 2, nodeList, tmpGR, ADD_VALUES);
// 	delete nodeList;
// 	delete tmpGR;
// }

// void NeuronGrowth::VisualizeVTK_ControlMesh(const vector<Vertex2D> &spt, const vector<Element2D> &mesh, int step, string fn, vector<float> var, string varName)
// {
// 	string fname;
// 	stringstream ss;
// 	ss << std::setw(6) << std::setfill('0') << step;
// 	ofstream fout;
// 	unsigned int i;
// 	fname = fn + "/controlmesh_" + varName + ss.str() + ".vtk";
// 	fout.open(fname.c_str());
// 	if (fout.is_open())
// 	{
// 		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
// 		fout << "POINTS " << spt.size() << " float\n";
// 		for (i = 0; i < spt.size(); i++)
// 		{
// 			fout << spt[i].coor[0] << " " << spt[i].coor[1] << " " << 0.00000 /* spt[i].coor[2] */ << "\n";
// 		}
// 		fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
// 		for (i = 0; i < mesh.size(); i++)
// 		{
// 			fout << "4 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3]
// 			     << " " /* << mesh[i].IEN[4] << " " << mesh[i].IEN[5] << " " << mesh[i].IEN[6] << " " << mesh[i].IEN[7] */ << '\n';
// 		}
// 		fout << "\nCELL_TYPES " << mesh.size() << '\n';
// 		for (i = 0; i < mesh.size(); i++)
// 		{
// 			fout << "9\n";
// 		}
// 		fout << "POINT_DATA " << var.size() << "\nSCALARS AllParticles float 1\nLOOKUP_TABLE default\n";
// 		for (uint i = 0; i < var.size(); i++)
// 		{
// 			fout << var[i] /* +N_plus[i]+N_minus[i] */ << "\n";
// 			// fout << spt[i].label /* +N_plus[i]+N_minus[i] */ << "\n";
// 		}
// 		fout.close();
// 	}
// 	else
// 	{
// 		cout << "Cannot open " << fname << "!\n";
// 	}
// }

// void NeuronGrowth::ConcentrationCal_Coupling_Bezier(float u, float v, const Element2D &bzel, float pt[2], float &disp, float dudx[2], float &detJ)
// {
// 	float dUdx[2][2];
// 	vector<float> Nx(bzel.IEN.size());
// 	vector<array<float, 2>> dNdx(bzel.IEN.size());
// 	bzel.Para2Phys(u, v, pt);
// 	BasisFunction(u, v, bzel.pts, bzel.cmat, Nx, dNdx, dUdx, detJ);
// 	disp = 0.;
// 	dudx[0] = 0.;
// 	dudx[1] = 0.; // dudx[2] = 0.;
// 	for (uint i = 0; i < bzel.IEN.size(); i++)
// 	{
// 		disp += Nx[i] * (N_0[bzel.IEN[i]]);	    // + N_plus[bzel.IEN[i]]);
// 		dudx[0] += dNdx[i][0] * (N_0[bzel.IEN[i]]); // + N_plus[bzel.IEN[i]]);
// 		dudx[1] += dNdx[i][1] * (N_0[bzel.IEN[i]]); // + N_plus[bzel.IEN[i]]);
// 							    // dudx[2] += dNdx[i][2] * (N_0[bzel.IEN[i]] + N_plus[bzel.IEN[i]]);
// 	}
// }

// void NeuronGrowth::VisualizeVTK_PhysicalDomain(int step, string var, string fn)
// {
// 	vector<array<float, 3>> spt_all; // sample points
// 	vector<float> sresult_all;
// 	vector<array<int, 4>> sele_all;
// 	float detJ;
// 	int num_bzmesh_ele = bzmesh_process.size();
// 	float spt_proc[num_bzmesh_ele * 4 * 3];
// 	float sresult_proc[num_bzmesh_ele * 4 * 1];
// 	int sele_proc[num_bzmesh_ele * 4];

// 	for (unsigned int e = 0; e < num_bzmesh_ele; e++)
// 	{
// 		int ns(2);
// 		vector<float> su(ns);
// 		for (int i = 0; i < ns; i++)
// 		{
// 			su[i] = float(i) / (float(ns) - 1.);
// 		}

// 		int loc(0);
// 		for (int a = 0; a < ns; a++)
// 		{
// 			for (int b = 0; b < ns; b++)
// 			{
// 				// for (int c = 0; c < ns; c++)
// 				// {
// 				float pt1[2], dudx[2];
// 				float result;
// 				ConcentrationCal_Coupling_Bezier(su[b], su[a], bzmesh_process[e], pt1, result, dudx, detJ);
// 				spt_proc[12 * e + loc * 3 + 0] = pt1[0];
// 				spt_proc[12 * e + loc * 3 + 1] = pt1[1];
// 				spt_proc[12 * e + loc * 3 + 2] = 0.000000; // pt1[2];
// 				sresult_proc[4 * e + loc] = result;
// 				loc++;
// 				// }
// 			}
// 		}
// 		int nns[2] = {ns * ns, ns};
// 		for (int a = 0; a < ns - 1; a++)
// 		{
// 			for (int b = 0; b < ns - 1; b++)
// 			{
// 				// for (int c = 0; c < ns - 1; c++)
// 				// {
// 				sele_proc[4 * e + 0] = 4 * e + a * ns + b;
// 				sele_proc[4 * e + 1] = 4 * e + a * ns + b + 1;
// 				sele_proc[4 * e + 2] = 4 * e + (a + 1) * ns + (b + 1);
// 				sele_proc[4 * e + 3] = 4 * e + (a + 1) * ns + b;
// 				// sele_proc[8 * e + 4] = 8 * e + (a + 1)*nns[1] + b*ns + c;
// 				// sele_proc[8 * e + 5] = 8 * e + (a + 1)*nns[1] + b*ns + c + 1;
// 				// sele_proc[8 * e + 6] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
// 				// sele_proc[8 * e + 7] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c;
// 				// }
// 			}
// 		}
// 	}

// 	float *spts = NULL;
// 	float *sresults = NULL;
// 	int *seles = NULL;
// 	int *displs_spts = NULL;
// 	int *displs_sresults = NULL;
// 	int *displs_seles = NULL;
// 	int *num_bzmesh_eles = NULL;
// 	int *recvcounts_spts = NULL;
// 	int *recvcounts_sresults = NULL;
// 	int *recvcounts_seles = NULL;

// 	if (comRank == 0)
// 	{
// 		num_bzmesh_eles = (int *)malloc(sizeof(int) * nProcess);
// 		recvcounts_spts = (int *)malloc(sizeof(int) * nProcess);
// 		recvcounts_sresults = (int *)malloc(sizeof(int) * nProcess);
// 		recvcounts_seles = (int *)malloc(sizeof(int) * nProcess);
// 	}
// 	MPI_Gather(&num_bzmesh_ele, 1, MPI_INT, num_bzmesh_eles, 1, MPI_INT, 0, PETSC_COMM_WORLD);
// 	MPI_Barrier(comm);

// 	if (comRank == 0)
// 	{
// 		spts = (float *)malloc(sizeof(float) * 12 * n_bzmesh);
// 		sresults = (float *)malloc(sizeof(float) * 4 * n_bzmesh);
// 		seles = (int *)malloc(sizeof(int) * 4 * n_bzmesh);

// 		displs_spts = (int *)malloc(nProcess * sizeof(int));
// 		displs_sresults = (int *)malloc(nProcess * sizeof(int));
// 		displs_seles = (int *)malloc(nProcess * sizeof(int));
// 		displs_spts[0] = 0;
// 		displs_sresults[0] = 0;
// 		displs_seles[0] = 0;

// 		for (int i = 1; i < nProcess; i++)
// 		{
// 			displs_spts[i] = displs_spts[i - 1] + num_bzmesh_eles[i - 1] * 12;
// 			displs_sresults[i] = displs_sresults[i - 1] + num_bzmesh_eles[i - 1] * 4;
// 			displs_seles[i] = displs_seles[i - 1] + num_bzmesh_eles[i - 1] * 4;
// 		}

// 		for (int i = 0; i < nProcess; i++)
// 		{
// 			recvcounts_spts[i] = num_bzmesh_eles[i] * 12;
// 			recvcounts_sresults[i] = num_bzmesh_eles[i] * 4;
// 			recvcounts_seles[i] = num_bzmesh_eles[i] * 4;
// 		}
// 	}

// 	MPI_Gatherv(spt_proc, num_bzmesh_ele * 4 * 3, MPI_FLOAT, spts, recvcounts_spts, displs_spts, MPI_FLOAT, 0, PETSC_COMM_WORLD);
// 	MPI_Gatherv(sresult_proc, num_bzmesh_ele * 4, MPI_FLOAT, sresults, recvcounts_sresults, displs_sresults, MPI_FLOAT, 0, PETSC_COMM_WORLD);
// 	MPI_Gatherv(sele_proc, num_bzmesh_ele * 4, MPI_INT, seles, recvcounts_seles, displs_seles, MPI_INT, 0, PETSC_COMM_WORLD);

// 	if (comRank == 0)
// 	{
// 		for (int i = 0; i < n_bzmesh; i++)
// 		{
// 			for (int j = 0; j < 4; j++)
// 			{
// 				array<float, 3> pt = {spts[i * 12 + j * 3 + 0], spts[i * 12 + j * 3 + 1], spts[i * 12 + j * 3 + 2]};
// 				spt_all.push_back(pt);
// 				sresult_all.push_back(sresults[i * 4 + j]);
// 			}
// 		}
// 		int sum_ele = 0;
// 		int pstart = 0;
// 		for (int i = 0; i < nProcess; i++)
// 		{
// 			for (int e = 0; e < num_bzmesh_eles[i]; e++)
// 			{
// 				array<int, 4> el;
// 				el[0] = pstart + seles[4 * sum_ele + 0];
// 				el[1] = pstart + seles[4 * sum_ele + 1];
// 				el[2] = pstart + seles[4 * sum_ele + 2];
// 				el[3] = pstart + seles[4 * sum_ele + 3];
// 				// el[4] = pstart + seles[8 * sum_ele + 4];
// 				// el[5] = pstart + seles[8 * sum_ele + 5];
// 				// el[6] = pstart + seles[8 * sum_ele + 6];
// 				// el[7] = pstart + seles[8 * sum_ele + 7];
// 				sele_all.push_back(el);
// 				sum_ele++;
// 			}
// 			pstart = pstart + num_bzmesh_eles[i] * 4;
// 		}
// 		// cout << "Visualizing in Physical Domain...\n";
// 		WriteVTK(spt_all, sresult_all, sele_all, step, var, fn);
// 	}
// }

// void NeuronGrowth::WriteVTK(const vector<array<float, 3>> spt, const vector<float> sdisp, const vector<array<int, 4>> sele, int step, string var, string fn)
// {
// 	stringstream ss;
// 	// ss << step;
// 	ss << std::setw(6) << std::setfill('0') << step;
// 	string fname = fn + "/" + var + "physics_singleparticle_" + ss.str() + ".vtk";
// 	ofstream fout;
// 	fout.open(fname.c_str());
// 	unsigned int i;
// 	if (fout.is_open())
// 	{
// 		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
// 		fout << "POINTS " << spt.size() << " float\n";
// 		for (i = 0; i < spt.size(); i++)
// 		{
// 			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
// 		}
// 		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
// 		for (i = 0; i < sele.size(); i++)
// 		{
// 			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
// 			     << " " /* << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] */ << '\n';
// 		}
// 		fout << "\nCELL_TYPES " << sele.size() << '\n';
// 		for (i = 0; i < sele.size(); i++)
// 		{
// 			fout << "9\n";
// 		}
// 		fout << "POINT_DATA " << sdisp.size() << "\nSCALARS AllParticle float 1\nLOOKUP_TABLE default\n";
// 		for (uint i = 0; i < sdisp.size(); i++)
// 		{
// 			fout << sdisp[i] << "\n";
// 		}
// 		fout.close();
// 	}
// 	else
// 	{
// 		cout << "Cannot open " << fname << "!\n";
// 	}
// }

// void NeuronGrowth::CalculateVarsForOutput(vector<array<float, 3>> &spt_all, vector<float> &sresult_all, vector<array<int, 4>> &sele_all) {
// 	float detJ;
// 	int num_bzmesh_ele = bzmesh_process.size();
// 	float spt_proc[num_bzmesh_ele * 4 * 3];
// 	float sresult_proc[num_bzmesh_ele * 4 * 1];
// 	int sele_proc[num_bzmesh_ele * 4];

// 	for (unsigned int e = 0; e < num_bzmesh_ele; e++)
// 	{
// 		int ns(2);
// 		vector<float> su(ns);
// 		for (int i = 0; i < ns; i++)
// 		{
// 			su[i] = float(i) / (float(ns) - 1.);
// 		}

// 		int loc(0);
// 		for (int a = 0; a < ns; a++)
// 		{
// 			for (int b = 0; b < ns; b++)
// 			{
// 				// for (int c = 0; c < ns; c++)
// 				// {
// 				float pt1[2], dudx[2];
// 				float result;
// 				ConcentrationCal_Coupling_Bezier(su[b], su[a], bzmesh_process[e], pt1, result, dudx, detJ);
// 				spt_proc[12 * e + loc * 3 + 0] = pt1[0];
// 				spt_proc[12 * e + loc * 3 + 1] = pt1[1];
// 				spt_proc[12 * e + loc * 3 + 2] = 0.000000; // pt1[2];
// 				sresult_proc[4 * e + loc] = result;
// 				loc++;
// 				// }
// 			}
// 		}
// 		int nns[2] = {ns * ns, ns};
// 		for (int a = 0; a < ns - 1; a++)
// 		{
// 			for (int b = 0; b < ns - 1; b++)
// 			{
// 				// for (int c = 0; c < ns - 1; c++)
// 				// {
// 				sele_proc[4 * e + 0] = 4 * e + a * ns + b;
// 				sele_proc[4 * e + 1] = 4 * e + a * ns + b + 1;
// 				sele_proc[4 * e + 2] = 4 * e + (a + 1) * ns + (b + 1);
// 				sele_proc[4 * e + 3] = 4 * e + (a + 1) * ns + b;
// 				// sele_proc[8 * e + 4] = 8 * e + (a + 1)*nns[1] + b*ns + c;
// 				// sele_proc[8 * e + 5] = 8 * e + (a + 1)*nns[1] + b*ns + c + 1;
// 				// sele_proc[8 * e + 6] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
// 				// sele_proc[8 * e + 7] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c;
// 				// }
// 			}
// 		}
// 	}

// 	float *spts = NULL;
// 	float *sresults = NULL;
// 	int *seles = NULL;
// 	int *displs_spts = NULL;
// 	int *displs_sresults = NULL;
// 	int *displs_seles = NULL;
// 	int *num_bzmesh_eles = NULL;
// 	int *recvcounts_spts = NULL;
// 	int *recvcounts_sresults = NULL;
// 	int *recvcounts_seles = NULL;

// 	if (comRank == 0)
// 	{
// 		num_bzmesh_eles = (int *)malloc(sizeof(int) * nProcess);
// 		recvcounts_spts = (int *)malloc(sizeof(int) * nProcess);
// 		recvcounts_sresults = (int *)malloc(sizeof(int) * nProcess);
// 		recvcounts_seles = (int *)malloc(sizeof(int) * nProcess);
// 	}
// 	MPI_Gather(&num_bzmesh_ele, 1, MPI_INT, num_bzmesh_eles, 1, MPI_INT, 0, PETSC_COMM_WORLD);
// 	MPI_Barrier(comm);

// 	if (comRank == 0)
// 	{
// 		spts = (float *)malloc(sizeof(float) * 12 * n_bzmesh);
// 		sresults = (float *)malloc(sizeof(float) * 4 * n_bzmesh);
// 		seles = (int *)malloc(sizeof(int) * 4 * n_bzmesh);

// 		displs_spts = (int *)malloc(nProcess * sizeof(int));
// 		displs_sresults = (int *)malloc(nProcess * sizeof(int));
// 		displs_seles = (int *)malloc(nProcess * sizeof(int));
// 		displs_spts[0] = 0;
// 		displs_sresults[0] = 0;
// 		displs_seles[0] = 0;

// 		for (int i = 1; i < nProcess; i++)
// 		{
// 			displs_spts[i] = displs_spts[i - 1] + num_bzmesh_eles[i - 1] * 12;
// 			displs_sresults[i] = displs_sresults[i - 1] + num_bzmesh_eles[i - 1] * 4;
// 			displs_seles[i] = displs_seles[i - 1] + num_bzmesh_eles[i - 1] * 4;
// 		}

// 		for (int i = 0; i < nProcess; i++)
// 		{
// 			recvcounts_spts[i] = num_bzmesh_eles[i] * 12;
// 			recvcounts_sresults[i] = num_bzmesh_eles[i] * 4;
// 			recvcounts_seles[i] = num_bzmesh_eles[i] * 4;
// 		}
// 	}

// 	MPI_Gatherv(spt_proc, num_bzmesh_ele * 4 * 3, MPI_FLOAT, spts, recvcounts_spts, displs_spts, MPI_FLOAT, 0, PETSC_COMM_WORLD);
// 	MPI_Gatherv(sresult_proc, num_bzmesh_ele * 4, MPI_FLOAT, sresults, recvcounts_sresults, displs_sresults, MPI_FLOAT, 0, PETSC_COMM_WORLD);
// 	MPI_Gatherv(sele_proc, num_bzmesh_ele * 4, MPI_INT, seles, recvcounts_seles, displs_seles, MPI_INT, 0, PETSC_COMM_WORLD);

// 	if (comRank == 0)
// 	{
// 		for (int i = 0; i < n_bzmesh; i++)
// 		{
// 			for (int j = 0; j < 4; j++)
// 			{
// 				array<float, 3> pt = {spts[i * 12 + j * 3 + 0], spts[i * 12 + j * 3 + 1], spts[i * 12 + j * 3 + 2]};
// 				spt_all.push_back(pt);
// 				sresult_all.push_back(sresults[i * 4 + j]);
// 			}
// 		}
// 		int sum_ele = 0;
// 		int pstart = 0;
// 		for (int i = 0; i < nProcess; i++)
// 		{
// 			for (int e = 0; e < num_bzmesh_eles[i]; e++)
// 			{
// 				array<int, 4> el;
// 				el[0] = pstart + seles[4 * sum_ele + 0];
// 				el[1] = pstart + seles[4 * sum_ele + 1];
// 				el[2] = pstart + seles[4 * sum_ele + 2];
// 				el[3] = pstart + seles[4 * sum_ele + 3];
// 				// el[4] = pstart + seles[8 * sum_ele + 4];
// 				// el[5] = pstart + seles[8 * sum_ele + 5];
// 				// el[6] = pstart + seles[8 * sum_ele + 6];
// 				// el[7] = pstart + seles[8 * sum_ele + 7];
// 				sele_all.push_back(el);
// 				sum_ele++;
// 			}
// 			pstart = pstart + num_bzmesh_eles[i] * 4;
// 		}
// 	}
// }
// void NeuronGrowth::VisualizeVTK_PhysicalDomain_All(int step, string fn)
// {
// 	vector<vector<array<float, 3>>> spt_all_4var; // sample points
// 	spt_all_4var.resize(4);
// 	vector<vector<float>> sresult_all_4var;
// 	sresult_all_4var.resize(4);
// 	vector<vector<array<int, 4>>> sele_all_4var;
// 	sele_all_4var.resize(4);
// 	N_0 = phi;
// 	CalculateVarsForOutput(spt_all_4var[0], sresult_all_4var[0], sele_all_4var[0]);
// 	N_0 = syn;
// 	// N_0 = theta;
// 	CalculateVarsForOutput(spt_all_4var[1], sresult_all_4var[1], sele_all_4var[1]);
// 	N_0 = tub;
// 	CalculateVarsForOutput(spt_all_4var[2], sresult_all_4var[2], sele_all_4var[2]);
// 	N_0 = tips;
// 	CalculateVarsForOutput(spt_all_4var[3], sresult_all_4var[3], sele_all_4var[3]);
// 	if (comRank == 0)
// 	{
// 		WriteVTK_All(spt_all_4var[0], sresult_all_4var, sele_all_4var[0], step, fn);
// 	}
// }

// void NeuronGrowth::WriteVTK_All(const vector<array<float, 3>> spt, const vector<vector<float>> sdisp, const vector<array<int, 4>> sele, int step, string fn)
// {
// 	stringstream ss;
// 	// ss << step;
// 	ss << std::setw(6) << std::setfill('0') << step;
// 	string fname = fn + "/physics_allparticle_" + ss.str() + ".vtk";
// 	ofstream fout;
// 	fout.open(fname.c_str());
// 	unsigned int i;
// 	if (fout.is_open())
// 	{
// 		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
// 		fout << "POINTS " << spt.size() << " float\n";
// 		for (i = 0; i < spt.size(); i++)
// 		{
// 			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
// 		}
// 		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
// 		for (i = 0; i < sele.size(); i++)
// 		{
// 			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
// 			     << " " /* << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] */ << '\n';
// 		}
// 		fout << "\nCELL_TYPES " << sele.size() << '\n';
// 		for (i = 0; i < sele.size(); i++)
// 		{
// 			fout << "9\n";
// 		}
// 		fout << "POINT_DATA " << sdisp[0].size() << "\nSCALARS " << "phi " << "float 1\nLOOKUP_TABLE default\n";
// 		for (uint i = 0; i < sdisp[0].size(); i++)
// 		{
// 			fout << sdisp[0][i] << "\n";
// 		}
// 		fout << "\nSCALARS " << "synaptogenesis " << "float 1\nLOOKUP_TABLE default\n";
// 		for (uint i = 0; i < sdisp[1].size(); i++)
// 		{
// 			fout << sdisp[1][i] << "\n";
// 		}
// 		fout << "\nSCALARS " << "tubulin " << "float 1\nLOOKUP_TABLE default\n";
// 		for (uint i = 0; i < sdisp[2].size(); i++)
// 		{
// 			fout << sdisp[2][i] << "\n";
// 		}
// 		fout << "\nSCALARS " << "tips " << "float 1\nLOOKUP_TABLE default\n";
// 		for (uint i = 0; i < sdisp[3].size(); i++)
// 		{
// 			fout << sdisp[3][i] << "\n";
// 		}
// 		fout.close();
// 	}
// 	else
// 	{
// 		cout << "Cannot open " << fname << "!\n";
// 	}
// }


// void NeuronGrowth::PointFormValue(vector<float> &Nx, const vector<float> &U, float Value)
// {
// 	Value = 0;

// 	for (int j = 0; j < Nx.size(); j++)
// 		Value += U[j] * Nx[j];
// }

// void NeuronGrowth::PointFormGrad(vector<array<float, 2>> &dNdx, const vector<float> &U, float Value[2])
// {
// 	for (int j = 0; j < 2; j++)
// 		Value[j] = 0.;

// 	for (int j = 0; j < 2; j++)
// 		for (int k = 0; k < dNdx.size(); k++)
// 			Value[j] += U[k] * dNdx[k][j];
// }

// void NeuronGrowth::PointFormHess(vector<array<array<float, 2>, 2>> &d2Ndx2, const vector<float> &U, float Value[2][2])
// {
// 	for (int j = 0; j < 2; j++)
// 	{
// 		for (int k = 0; k < 2; k++)
// 		{
// 			Value[j][k] = 0.;
// 		}
// 	}
// 	for (int j = 0; j < 2; j++)
// 	{
// 		for (int k = 0; k < 2; k++)
// 		{
// 			for (int l = 0; l < d2Ndx2.size(); l++)
// 			{
// 				Value[j][k] += U[l] * d2Ndx2[l][j][k];
// 			}
// 		}
// 	}
// }

// void NeuronGrowth::ElementValue(const vector<float> &Nx, const vector<float> value_node, float &value)
// {
// 	value = 0.;
// 	for (int i = 0; i < Nx.size(); i++)
// 		value += value_node[i] * Nx[i];
// }

// void NeuronGrowth::ElementValueAll(const vector<float> &Nx, const vector<float> elePhiGuess, float &elePG, const vector<float> elePhi, float &eleP, const vector<float> eleSyn, float &eleS, const vector<float> eleTips, float &eleTp, const vector<float> eleTubulin, float &eleTb, const vector<float> eleEpsilon, float &eleEP, const vector<float> eleEpsilonP, float &eleEEP)
// {
// 	elePG = 0.;
// 	eleP = 0.;
// 	eleS = 0.;
// 	eleTp = 0.;
// 	eleTb = 0.;
// 	eleEP = 0.;
// 	eleEEP = 0.;
// 	for (int i = 0; i < Nx.size(); i++) {
// 		elePG += elePhiGuess[i] * Nx[i];
// 		eleP += elePhi[i] * Nx[i];
// 		eleS += eleSyn[i] * Nx[i];
// 		eleTp += eleTips[i] * Nx[i];
// 		eleTb += eleTubulin[i] * Nx[i];
// 		eleEP += eleEpsilon[i] * Nx[i];
// 		eleEEP += eleEpsilonP[i] * Nx[i];
// 	}
// }

// void NeuronGrowth::ElementDeriv(const int nen, vector<array<float, 2>> &dNdx, const vector<float> value_node, float &dVdx, float &dVdy)
// {
// 	dVdx = 0;
// 	dVdy = 0.;
// 	for (int i = 0; i < nen; i++) {
// 		dVdx += value_node[i] * dNdx[i][0];
// 		dVdy += value_node[i] * dNdx[i][1];
// 	}
// }


// void NeuronGrowth::ElementDerivAll(const int nen, vector<array<float, 2>> &dNdx, const vector<float> elePhiGuess, float &dPGdx, float &dPGdy, const vector<float> eleTheta, float &dThedx, float &dThedy, const vector<float> eleEpsilon, float &dAdx, float &dAdy, const vector<float> eleEpsilonP, float &dAPdx, float &dAPdy)
// {	
// 	dPGdx = 0; dPGdy = 0.;
// 	dThedx = 0; dThedy = 0.;
// 	dAdx = 0; dAdy = 0.;
// 	dAPdx = 0; dAPdy = 0.;
// 	for (int i = 0; i < nen; i++) {
// 		dPGdx += elePhiGuess[i] * dNdx[i][0];	dPGdy += elePhiGuess[i] * dNdx[i][1];
// 		dThedx += eleTheta[i] * dNdx[i][0];	dThedy += eleTheta[i] * dNdx[i][1];
// 		dAdx += eleEpsilon[i] * dNdx[i][0];	dAdy += eleEpsilon[i] * dNdx[i][1];
// 		dAPdx += eleEpsilonP[i] * dNdx[i][0];	dAPdy += eleEpsilonP[i] * dNdx[i][1];
// 	}
// }

// void NeuronGrowth::ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 2>> &dNdx, const vector<float> elePhiGuess, float &elePG, const vector<float> elePhi, float &eleP, const vector<float> eleEpsilon, float &eleEP, const vector<float> eleEpsilonP, float &eleEEP, float &dPGdx, float &dPGdy, float &dAdx, float &dAdy, float &dAPdx, float &dAPdy)
// {
// 	elePG = 0.;
// 	eleP = 0.;
// 	eleEP = 0.;
// 	eleEEP = 0.;

// 	dPGdx = 0; dPGdy = 0.;
// 	dAdx = 0; dAdy = 0.;
// 	dAPdx = 0; dAPdy = 0.;

// 	for (int i = 0; i < nen; i++) {
// 		elePG += elePhiGuess[i] * Nx[i];
// 		eleP += elePhi[i] * Nx[i];
// 		eleEP += eleEpsilon[i] * Nx[i];
// 		eleEEP += eleEpsilonP[i] * Nx[i];

// 		dPGdx += elePhiGuess[i] * dNdx[i][0];	dPGdy += elePhiGuess[i] * dNdx[i][1];
// 		dAdx += eleEpsilon[i] * dNdx[i][0];	dAdy += eleEpsilon[i] * dNdx[i][1];
// 		dAPdx += eleEpsilonP[i] * dNdx[i][0];	dAPdy += eleEpsilonP[i] * dNdx[i][1];
// 	}
// }


// void NeuronGrowth::ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 2>> &dNdx, vector<vector<float>> &eleVal, vector<float> &vars)
// {
// 	// eleVal Array Definitions:
// 	// Index	Description
// 	// -----	-----------
// 	// 0		elePhiGuess	
// 	// 1		elePhi		
// 	// 2		eleTheta	
// 	// 3		eleEpsilon	
// 	// 4		eleEpsilonP	

// 	// Variables Definition:
// 	// Index   Description
// 	// -----   -----------
// 	// 0       elePG    - element value for PhiGuess
// 	// 1       dPGdx    - derivative of PhiGuess w.r.t. x
// 	// 2       dPGdy    - derivative of PhiGuess w.r.t. y
// 	// 3       eleP     - element value for Phi
// 	// 4       eleEP    - element value for epsilonP
// 	// 5       eleEEP   - element value for epsilonEP
// 	// 6       dAdx     - derivative of A w.r.t. x
// 	// 7       dAdy     - derivative of A w.r.t. y
	
// 	vars[0] = 0.;
// 	vars[3] = 0.;
// 	vars[4] = 0.;
// 	vars[5] = 0.;

// 	vars[1] = 0; vars[2] = 0.;
// 	vars[6] = 0; vars[7] = 0.;

// 	for (int i = 0; i < nen; i++) {
// 		vars[0] += eleVal[0][i] * Nx[i];
// 		vars[3] += eleVal[1][i] * Nx[i];
// 		vars[4] += eleVal[3][i] * Nx[i];
// 		vars[5] += eleVal[4][i] * Nx[i];

// 		vars[1] += eleVal[0][i] * dNdx[i][0];	vars[2] += eleVal[0][i] * dNdx[i][1];
// 		vars[6] += eleVal[3][i] * dNdx[i][0];	vars[7] += eleVal[3][i] * dNdx[i][1];
// 	}
// }

// // void NeuronGrowth::ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 2>> &dNdx, vector<vector<float>> &eleVal, vector<float> &vars)
// void NeuronGrowth::ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 2>> &dNdx, vector<float> &elePhiGuess, vector<float> &vars)
// {
// 	// eleVal Array Definitions:
// 	// Index	Description
// 	// -----	-----------
// 	// 0		elePhiGuess	
// 	// 1		elePhi		
// 	// 2		eleTheta	
// 	// 3		eleEpsilon	
// 	// 4		eleEpsilonP	

// 	// Variables Definition:
// 	// Index   Description
// 	// -----   -----------
// 	// 0       elePG    - element value for PhiGuess
// 	// 1       dPGdx    - derivative of PhiGuess w.r.t. x
// 	// 2       dPGdy    - derivative of PhiGuess w.r.t. y
// 	// 3       eleP     - element value for Phi
// 	// 4       eleEP    - element value for epsilonP
// 	// 5       eleEEP   - element value for epsilonEP
// 	// 6       dAdx     - derivative of A w.r.t. x
// 	// 7       dAdy     - derivative of A w.r.t. y
	
// 	vars[0] = 0.;
// 	// vars[3] = 0.;
// 	// vars[4] = 0.;
// 	// vars[5] = 0.;

// 	vars[1] = 0; vars[2] = 0.;
// 	// vars[6] = 0; vars[7] = 0.;

// 	for (int i = 0; i < nen; i++) {
// 		vars[0] += elePhiGuess[i] * Nx[i];
// 		// vars[3] += eleVal[1][i] * Nx[i];
// 		// vars[4] += eleEpsilon[i] * Nx[i];
// 		// vars[5] += eleVal[4][i] * Nx[i];

// 		vars[1] += elePhiGuess[i] * dNdx[i][0];	vars[2] += elePhiGuess[i] * dNdx[i][1];
// 		// vars[6] += eleEpsilon[i] * dNdx[i][0];	vars[7] += eleEpsilon[i] * dNdx[i][1];
// 	}
// }

// void NeuronGrowth::ElementEvaluationAll_syn_tub(const int nen, const vector<float> &Nx, vector<array<float, 2>> &dNdx, const vector<float> elePhiDiff, float &elePf, const vector<float> eleSyn, float &eleS, const vector<float> elePhi, float &eleP, const vector<float> elePhiPrev, float &elePprev, const vector<float> eleConct, float &eleC, float &dPdx, float &dPdy)
// {
// 	elePf = 0.;
// 	eleS = 0.;
// 	eleP = 0.;
// 	elePprev = 0.;
// 	eleC = 0.;

// 	dPdx = 0; dPdy = 0.;

// 	for (int i = 0; i < nen; i++) {
// 		elePf += elePhiDiff[i] * Nx[i];
// 		eleS += eleSyn[i] * Nx[i];

// 		eleP += elePhi[i] * Nx[i];
// 		elePprev += elePhiPrev[i] * Nx[i];
// 		eleC += eleConct[i] * Nx[i];

// 		dPdx += elePhi[i] * dNdx[i][0];	dPdy += elePhi[i] * dNdx[i][1];
// 	}
// }

// void NeuronGrowth::ElementEvaluationAll_syn_tub(const int nen, const vector<float> &Nx, vector<array<float, 2>> &dNdx, vector<vector<float>> &eleVal, vector<float> &vars)
// {
// 	// eleVal Array Definitions:
// 	// Index    Description
// 	// -----    -----------
// 	// 0        elePhiDiff  // Difference in Phi values
// 	// 1        eleSyn      // Synapse/electrical signal strength or similar context
// 	// 2        elePhi      // Current Phi value
// 	// 3        elePhiPrev  // Previous Phi value
// 	// 4        eleConct    // Concentration or a similar quantity

// 	// var Array Definitions:
// 	// Index    Variable        Description
// 	// -----    --------        -----------
// 	// 0        elePf           // Placeholder for explanation
// 	// 1        eleS            // Placeholder for explanation
// 	// 2        eleP            // Placeholder for explanation
// 	// 3        dPdx            // Derivative of P with respect to x
// 	// 4        dPdy            // Derivative of P with respect to y
// 	// 5        elePprev        // Previous P value
// 	// 6        eleC            // Placeholder for explanation
// 	// 7        mag_grad_phi0   // Magnitude of the gradient of phi0
// 	// 8        term_diff       // Diffusion term or similar
// 	// 9        term_alph       // Alpha term, specific to the model's context
// 	// 10       term_beta       // Beta term, specific to the model's context
// 	// 11       term_source     // Source term, specific to the model's context

// 	vars[0] = 0.;
// 	vars[1] = 0.;
// 	vars[2] = 0.;
// 	vars[5] = 0.;
// 	vars[6] = 0.;

// 	vars[3] = 0; vars[4] = 0.;

// 	for (int i = 0; i < nen; i++) {
// 		vars[0] += eleVal[0][i] * Nx[i];
// 		vars[1] += eleVal[1][i] * Nx[i];

// 		vars[2] += eleVal[2][i] * Nx[i];
// 		vars[5] += eleVal[3][i] * Nx[i];
// 		vars[6] += eleVal[4][i] * Nx[i];

// 		vars[3] += eleVal[2][i] * dNdx[i][0];	vars[4] += eleVal[2][i] * dNdx[i][1];
// 	}
// }

// void NeuronGrowth::prepareBasis() {

// 	/*Build linear system in each process*/
// 	float detJ;
// 	float dudx[2][2];
	
// 	float dThedx, dThedy;
// 	sum_grad_phi0_local = 0;

// 	float eleP0(0);
// 	int e;
// 	// std::cout << comRank << " " << bzmesh_process.size() << std::endl;

// 	for (e = 0; e < bzmesh_process.size(); e++) {
// 		int nen = bzmesh_process[e].IEN.size();

// 		// if (comRank == 1)
// 			// std::cout << comRank << " " << e << "/" << bzmesh_process.size() << std::endl;

// 		elePhi0.resize(nen);
// 		eleTheta.resize(nen);
		
// 		for (int i = 0; i < nen; i++) {
// 			elePhi0[i] = phi_0[bzmesh_process[e].IEN[i]];
// 			eleTheta[i] = theta[bzmesh_process[e].IEN[i]];
// 		}

// 		for (int i = 0; i < Gpt.size(); i++) {
// 			for (int j = 0; j < Gpt.size(); j++) {
// 				BasisFunction(Gpt[i], Gpt[j], bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dudx, detJ);
// 				detJ = wght[i] * wght[j] * detJ;

// 				pre_Nx.push_back(Nx);
// 				pre_dNdx.push_back(dNdx);
// 				pre_detJ.push_back(detJ);

// 				ElementDeriv(nen, dNdx, eleTheta, dThedx, dThedy);
// 				pre_C0.push_back((0.5 + 6 * s_coeff * sqrt(pow(dThedx, 2) + pow(dThedy, 2))));

// 				ElementDeriv(nen, dNdx, elePhi0, dP0dx, dP0dy);

// 				pre_mag_grad_phi0.push_back(pow(dP0dx, 2) + pow(dP0dy, 2));
// 				sum_grad_phi0_local += pow(dP0dx, 2) + pow(dP0dy, 2);

// 				// ElementValue(Nx, elePhi0, eleP0);
// 				// pre_term_source.push_back(source_coeff * eleP0);

// 			}
// 		}
// 	}
// 	MPI_Barrier(PETSC_COMM_WORLD);
// }

// void NeuronGrowth::preparePhaseField() {

// 	pre_eleEP.clear();
// 	pre_eleEEP.clear();
// 	pre_dAdx.clear();
// 	pre_dAdy.clear();
// 	pre_eleP.clear();
// 	// pre_dPdx.clear();
// 	// pre_dPdy.clear();
// 	pre_eleTh.clear();
// 	pre_eleMp.clear();
// 	pre_C1.clear();

// 	float eleEP(0), eleEEP(0), dAdx(0), dAdy(0), 
// 		eleP(0), eleTh(0), eleS(0), eleTb(0), eleTp(0), eleE(0), eleMp(0); 
// 	int e, ind(0);

// 	for (e = 0; e < bzmesh_process.size(); e++) {
// 		int nen = bzmesh_process[e].IEN.size();

// 		vector<float> elePhi(nen, 0), eleTheta(nen, 0), eleEpsilon(nen, 0), eleEpsilonP(nen, 0);
// 		vector<float> eleSyn(nen, 0), eleTubulin(nen, 0), eleTips(nen, 0), eleMphi(nen, 0);
// 		// elePhi.resize(nen);
// 		// eleTheta.resize(nen);
// 		// eleEpsilon.resize(nen);
// 		// eleEpsilonP.resize(nen);
// 		// eleSyn.resize(nen);
// 		// eleTubulin.resize(nen);
// 		// eleTips.resize(nen);
// 		// eleMphi.resize(nen);
		
// 		// std::cout << Mphi.size() << std::endl;

// 		for (int i = 0; i < nen; i++) {

// 			elePhi[i] = phi[bzmesh_process[e].IEN[i]];			// elePhi
// 			eleTheta[i] = theta[bzmesh_process[e].IEN[i]];			// eleTheta
// 			// eleEpsilon[i] = 0;						// eleEpsilon
// 			// eleEpsilonP[i] = 0;						// eleEpsilonP

// 			eleSyn[i] = syn[bzmesh_process[e].IEN[i]];			// eleSyn
// 			eleTubulin[i] = tub[bzmesh_process[e].IEN[i]];			// eleTubulin
// 			eleTips[i] = tips[bzmesh_process[e].IEN[i]];			// eleTips
// 			// if (n < 1) {
// 			// 	eleMphi[i] = 60;
// 			// } else {
// 			// 	eleMphi[i] = Mphi[bzmesh_process[e].IEN[i]];
// 			// }
// 		}

// 		for (int i = 0; i < Gpt.size(); i++) {
// 			for (int j = 0; j < Gpt.size(); j++) {

// 				EvaluateOrientation(nen, pre_Nx[ind], pre_dNdx[ind], elePhi, eleTheta, eleEpsilon, eleEpsilonP);
// 				ElementValue(pre_Nx[ind], eleEpsilon, eleEP);
// 				pre_eleEP.push_back(eleEP);
// 				ElementValue(pre_Nx[ind], eleEpsilonP, eleEEP);
// 				pre_eleEEP.push_back(eleEEP);
// 				ElementDeriv(nen, pre_dNdx[ind], eleEpsilonP, dAdx, dAdy);
// 				pre_dAdx.push_back(dAdx);
// 				pre_dAdy.push_back(dAdy);

// 				ElementValue(pre_Nx[ind], elePhi, eleP);
// 				pre_eleP.push_back(eleP);
// 				// ElementDeriv(nen, pre_Nx[ind], elePhi, dPdx, dPdy);
// 				// pre_dPdx.push_back(dPdx);
// 				// pre_dPdy.push_back(dPdy);

// 				ElementValue(pre_Nx[ind], eleTheta, eleTh);
// 				pre_eleTh.push_back(eleTh);

// 				ElementValue(pre_Nx[ind], eleSyn, eleS);
// 				ElementValue(pre_Nx[ind], eleTubulin, eleTb);
// 				ElementValue(pre_Nx[ind], eleTips, eleTp);

// 				// float eleMp;
// 				// ElementValue(pre_Nx[ind], eleMphi, eleMp);
// 				// pre_eleMp.push_back(eleMp);

// 				// adjust rg (assembly rate) and sg (disassembly rate) based on detected tips
// 				if (n < 1) {
// 					eleE = alphaOverPi*atan(gamma * (1 - eleS));
// 					pre_eleMp.push_back(50);
// 				} else {
// 					if (eleTp > 0) {
// 						eleE = alphaOverPi*atan(gamma * Regular_Heiviside_fun(50 * eleTb - 0) * (1 - eleS));
// 						if (eleTp > 2) {
// 							// pre_eleMp.push_back(120);
// 							pre_eleMp.push_back(50);
// 						} else {
// 							// pre_eleMp.push_back(20);
// 							pre_eleMp.push_back(50);
// 						}
// 					} else {
// 						eleE = alphaOverPi*atan(gamma * Regular_Heiviside_fun(r * eleTb - g) * (1 - eleS));
// 						// pre_eleMp.push_back(1);
// 						pre_eleMp.push_back(5);
// 					}
// 				}

// 				// calculate C1 variable for phase field energy term
// 				float C1 = eleE - pre_C0[ind]; // C1 = E - (0.5 + 6 * user->s_coeff * sqrt(pow(dThedx, 2) + pow(dThedy, 2));
// 				pre_C1.push_back(C1);
										
// 				ind += 1;
// 			}
// 		}
// 	}
// 	// std::cout << pre_eleMp.size() << std::endl;
// 	// MPI_Barrier(PETSC_COMM_WORLD);
// }

// void NeuronGrowth::prepareSourceSum() {

// 	/*Build linear system in each process*/
// 	sum_grad_phi0_local = 0;

// 	float eleP(0), dPdx(0), dPdy(0);
// 	int e, ind(0);
// 	// std::cout << comRank << " " << bzmesh_process.size() << std::endl;

// 	vector<float> elePhi;
// 	for (e = 0; e < bzmesh_process.size(); e++) {
// 		int nen = bzmesh_process[e].IEN.size();

// 		elePhi.resize(nen);
		
// 		for (int i = 0; i < nen; i++) {
// 			elePhi[i] = phi[bzmesh_process[e].IEN[i]];
// 		}

// 		for (int i = 0; i < Gpt.size(); i++) {
// 			for (int j = 0; j < Gpt.size(); j++) {
// 				ElementDeriv(nen, pre_dNdx[ind], elePhi, dPdx, dPdy);

// 				pre_mag_grad_phi0.push_back(pow(dPdx, 2) + pow(dPdy, 2));
// 				sum_grad_phi0_local += pow(dPdx, 2) + pow(dPdy, 2);

// 				ind += 1;
// 			}
// 		}
// 	}
// 	// MPI_Barrier(PETSC_COMM_WORLD);
// }

// void NeuronGrowth::prepareTerm_source() {
// 	int ind(0);
// 	pre_term_source.clear();
// 	for (int e = 0; e < bzmesh_process.size(); e++) {
// 		for (int i = 0; i < Gpt.size(); i++) {
// 			for (int j = 0; j < Gpt.size(); j++) {
// 				pre_term_source.push_back(source_coeff * pre_mag_grad_phi0[ind] / sum_grad_phi0_global);
// 			}
// 		}
// 	}
// }

// void NeuronGrowth::prepareEE()
// {
// 	/*Build linear system in each process*/
// 	PetscInt ind(0), e;
// 	for (e = 0; e < bzmesh_process.size(); e++) {
// 		PetscInt nen = bzmesh_process[e].IEN.size();

// 		for (PetscInt i = 0; i < nen * 1; i++) {
// 			EVectorSolve[i] = 0.0;
// 			for (PetscInt j = 0; j < nen * 1; j++) {
// 				EMatrixSolve[i][j] = 0.0;
// 			}
// 		}
	
// 		for (PetscInt i = 0; i < nen * 1; i++) {
// 			eleVal[0][i] = phi[bzmesh_process[e].IEN[i]];		// elePhiGuess
// 			// eleVal[1][i] = phi[bzmesh_process[e].IEN[i]];		// elePhi
// 			eleVal[0][i] = 0;		// elePhiGuess
// 			eleVal[1][i] = 0;		// elePhi

// 			eleVal[2][i] = syn[bzmesh_process[e].IEN[i]];		// eleSyn
// 			eleVal[3][i] = tub[bzmesh_process[e].IEN[i]];		// eleTubulin
// 			eleVal[4][i] = theta[bzmesh_process[e].IEN[i]];	// eleTheta
// 			eleVal[5][i] = tips[bzmesh_process[e].IEN[i]];	// eleTips
// 			eleVal[6][i] = 0;							// eleEpsilon
// 			eleVal[7][i] = 0;							// eleEpsilonP
// 		}

// 		for (PetscInt i = 0; i < Gpt.size(); i++) {
// 			for (PetscInt j = 0; j < Gpt.size(); j++) {

// 				PetscReal detJ;

// 				vector<float> Nx;
// 				vector<array<float, 2>> dNdx;
// 				Nx.clear();
// 				dNdx.clear();
// 				Nx.resize(nen);
// 				dNdx.resize(nen);

// 				// use pre-calculated values (save computational cost)
// 				Nx = pre_Nx[ind];
// 				dNdx = pre_dNdx[ind];
// 				detJ = pre_detJ[ind];
// 				vars[1] = pre_C0[ind]; // C0
// 				ind += 1;

// 				EvaluateOrientation(nen, Nx, dNdx, eleVal[1], eleVal[4], eleVal[6], eleVal[7]);
// 				ElementEvaluationAll_phi(nen, Nx, dNdx, eleVal, vars);
// 				//	 0   1    2     3      4      5     6     7      8     9      10      11      12     13    14    15      16     17     
// 				// float C1, C0, elePG, dPGdx, dPGdy, eleP, eleS, eleTb, eleE, eleTp, dThedx, dThedy, eleEP, dAdx, dAdy, eleEEP, dAPdx, dAPdy;
				
// 				vars[8] = alphaOverPi*atan(gamma * (1 - vars[6]));

// 				vars[0] = vars[8] - vars[1]; // E - (0.5 + 6 * s_coeff * sqrt(pow(dThedx, 2) + pow(dThedy, 2));
									
// 				for (PetscInt m = 0; m < nen; m++) {
// 						EVectorSolve[m] += (vars[2] * Nx[m] - dt * M_phi *\
// 						((- vars[12] * vars[12] * (vars[3] * dNdx[m][0] + vars[4] * dNdx[m][1])) -\
// 						(- vars[13] * vars[15] * vars[4] * dNdx[m][0]) +\
// 						(- vars[14] * vars[15] * vars[3] * dNdx[m][1]) +\
// 						(- vars[2] * vars[2] * vars[2] + (1 - vars[0]) * vars[2] * vars[2] + vars[0] * vars[2]) * Nx[m])
// 						- vars[5] * Nx[m]) * detJ;
// 					for (PetscInt n = 0; n < nen; n++) {
// 						EMatrixSolve[m][n] += (Nx[m] * Nx[n] - dt * M_phi *
// 							((- vars[12] * vars[12] * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1])) // terma2
// 							- (- vars[13] * vars[15] * dNdx[m][1] * dNdx[n][0]) // termadx
// 							+ (- vars[14] * vars[15] * dNdx[m][0] * dNdx[n][1]) // termady
// 							+ (- 3 * vars[2] * vars[2] + 2 * (1 - vars[0]) * vars[2] + vars[0] * Nx[m]) * Nx[n]) // termdbl
// 							) * detJ;
// 					}
// 				}
// 			}
// 		}

// 		/*Apply Boundary Condition*/
// 		for (PetscInt i = 0; i < nen; i++) {
// 			PetscInt A = bzmesh_process[e].IEN[i];
// 			if (cpts[A].label == 1) // domain boundary
// 				ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
// 		}
		
// 		pre_EMatrixSolve.push_back(EMatrixSolve);
// 		pre_EVectorSolve.push_back(EVectorSolve);
// 	}
// }

// void NeuronGrowth::EvaluateEnergy(const int nen, const vector<float> &Nx, const vector<float> eleSyn, vector<float>& E)
// {
// 	for (int i = 0; i < nen; i++) {
// 		E[i] = alphaOverPi*atan(gamma*(1-eleSyn[i]));
// 	}

// }

// float NeuronGrowth::Regular_Heiviside_fun(float x) 
// {
// 	// float epsilon = 0.0001; // the number is not fixed.
// 	float epsilon = 1e-15; // the number is not fixed.
// 	return 0.5*(1+(2/PI)*atan(x/epsilon));
// }

// void NeuronGrowth::EvaluateOrientation(const int nen, const vector<float> &Nx, const vector<array<float, 2>> &dNdx, const vector<float> elePhi, const vector<float> eleTheta,  vector<float>& eleEpsilon, vector<float>& eleEpsilonP)
// {
// 	for (int i = 0; i < nen; i++) {
// 		eleEpsilon[i] += epsilonb * (1.0 + delta * cos(aniso * (atan2(elePhi[i] * dNdx[i][0], elePhi[i] * dNdx[i][1]) - eleTheta[i] * Nx[i])));
// 		eleEpsilonP[i] += -epsilonb * (aniso * delta * sin(aniso * (atan2(elePhi[i] * dNdx[i][0], elePhi[i] * dNdx[i][1]) - eleTheta[i] * Nx[i])));
// 	}
// }

// // void NeuronGrowth::BuildLinearSystemProcessNG_phi()
// // {
// // 	/*Build linear system in each process*/
// // 	int ind(0); // pre-calculated variable index
// // 	for (int e = 0; e < bzmesh_process.size(); e++) {
// // 		int nen = bzmesh_process[e].IEN.size(); // 16 supporting cp for 2D case

// // 		vector<vector<float>> EMatrixSolve;	EMatrixSolve.clear();
// // 		vector<float> EVectorSolve;		EVectorSolve.clear();
// // 		EMatrixSolve.resize(nen);
// // 		EVectorSolve.resize(nen);

// // 		vector<float> elePhiGuess(nen, 0);

// // 		for (int i = 0; i < nen * 1; i++) {
// // 			// EMatrixSolve[i].reserve(nen*sizeof(float));
// // 			EMatrixSolve[i].resize(nen);
// // 		}
		
// // 		// preparing element stiffness matrix and vector, extract values from control mesh
// // 		for (int i = 0; i < nen; i++) {

// // 			EVectorSolve[i] = 0.0;
// // 			for (int j = 0; j < nen; j++) {
// // 				EMatrixSolve[i][j] = 0.0;
// // 			}
// // 			elePhiGuess[i] = phi[bzmesh_process[e].IEN[i]];		// elePhiGuess
// // 		}

// // 		// loop through gaussian quadrature points
// // 		for (int i = 0; i < Gpt.size(); i++) {
// // 			for (int j = 0; j < Gpt.size(); j++) {				
// // 				ElementEvaluationAll_phi(nen, pre_Nx[ind], pre_dNdx[ind], elePhiGuess, vars);
				
// // 				// loop through control points
// // 				for (int m = 0; m < nen; m++) {
// // 					// terma2 = - eleEP * eleEP * (dPGdx * dNdx[m][0] + dPGdy * dNdx[m][1]);
// // 					// termadx = - dAdx * eleEEP * dPGdy * dNdx[m][0];
// // 					// termady = - dAdy * eleEEP * dPGdx * dNdx[m][1];
// // 					// termdbl = (- elePG * elePG * elePG + (1 - C1) * elePG * elePG + C1 * elePG) * Nx[m];
// // 					// EVectorSolve[m] += (elePG * Nx[m] - dt * M_phi * (terma2 - termadx + termady + termdbl) - eleP * Nx[m]) * detJ;
			
// // 					EVectorSolve[m] += ( - vars[0] * pre_Nx[ind][m] - dt * pre_eleMp[ind] *\
// // 						(
// // 						// (- pre_eleEP[ind] * pre_eleEP[ind] * (vars[1] * pre_dNdx[ind][m][0] + vars[2] * pre_dNdx[ind][m][1])) -\
// // 						// (- pre_dAdx[ind] * pre_eleEEP[ind] * vars[2] * pre_dNdx[ind][m][0]) +\
// // 						// (- pre_dAdy[ind] * pre_eleEEP[ind] * vars[1] * pre_dNdx[ind][m][1]) +\

// // 						( - vars[0] * vars[0] * vars[0] + (1 - pre_C1[ind]) * vars[0] * vars[0] + pre_C1[ind] * vars[0]) * pre_Nx[ind][m])
// // 						+ pre_eleP[ind] * pre_Nx[ind][m]) * pre_detJ[ind];

// // 					// loop through 16 control points
// // 					for (int n = 0; n < nen; n++) {
// // 						// terma2 = - eleEP * eleEP * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1]);
// // 						// termadx = - dAdx * eleEEP * dNdx[m][1] * dNdx[n][0];
// // 						// termady = - dAdy * eleEEP * dNdx[m][0] * dNdx[n][1];
// // 						// termdbl = (- 3 * elePG * elePG + 2 * (1 - C1) * elePG + C1 * Nx[m]) * Nx[n];
// // 						// EMatrixSolve[m][n] += (Nx[m] * Nx[n] - dt * M_phi * (terma2 - termadx + termady + termdbl)) * detJ;
						
// // 						EMatrixSolve[m][n] += ( - pre_Nx[ind][m] * pre_Nx[ind][n] - dt * pre_eleMp[ind] *
// // 							(
// // 							// (- pre_eleEP[ind] * pre_eleEP[ind] * (pre_dNdx[ind][m][0] * pre_dNdx[ind][n][0] + pre_dNdx[ind][m][1] * pre_dNdx[ind][n][1])) // terma2
// // 							// - (- pre_dAdx[ind] * pre_eleEEP[ind] * pre_dNdx[ind][m][1] * pre_dNdx[ind][n][0]) // termadx
// // 							// + (- pre_dAdy[ind] * pre_eleEEP[ind] * pre_dNdx[ind][m][0] * pre_dNdx[ind][n][1]) // termady
							
// // 							( - 3 * vars[0] * vars[0] + 2 * (1 - pre_C1[ind]) * vars[0] + pre_C1[ind] * pre_Nx[ind][m]) * pre_Nx[ind][n]) // termdbl
// // 							) * pre_detJ[ind];
// // 					}
// // 				}
// // 				ind += 1; // incrementing index for extracting pre-calculated variables
// // 			}
// // 		}

// // 		/*Apply Boundary Condition*/
// // 		for (int i = 0; i < nen; i++) {
// // 			int A = bzmesh_process[e].IEN[i];
// // 			if (cpts[A].label == 1) { // domain boundary
// // 				ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
// // 			}
// // 			// if (tips[A] == 0) { // domain boundary
// // 			// 	ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
// // 			// }
// // 		}

// // 		ResidualAssembly(EVectorSolve, bzmesh_process[e].IEN, GR_phi);
// // 		MatrixAssembly(EMatrixSolve, bzmesh_process[e].IEN, GK_phi);
// // 	}

// // 	if (comRank == 0) {
// // 		check_itr += 1;
// // 	}

// // 	VecAssemblyBegin(GR_phi);
// // 	MatAssemblyBegin(GK_phi, MAT_FINAL_ASSEMBLY);
// // }

// void NeuronGrowth::BuildLinearSystemProcessNG_phi()
// {
// 	/*Build linear system in each process*/
// 	int ind(0); // pre-calculated variable index
// 	for (int e = 0; e < bzmesh_process.size(); e++) {
// 		int nen = bzmesh_process[e].IEN.size(); // 16 supporting cp for 2D case

// 		vector<vector<float>> EMatrixSolve;	EMatrixSolve.clear();
// 		vector<float> EVectorSolve;		EVectorSolve.clear();
// 		EMatrixSolve.resize(nen);
// 		EVectorSolve.resize(nen);

// 		vector<float> elePhiGuess(nen, 0);

// 		for (int i = 0; i < nen * 1; i++) {
// 			// EMatrixSolve[i].reserve(nen*sizeof(float));
// 			EMatrixSolve[i].resize(nen);
// 		}
		
// 		// preparing element stiffness matrix and vector, extract values from control mesh
// 		for (int i = 0; i < nen; i++) {

// 			EVectorSolve[i] = 0.0;
// 			for (int j = 0; j < nen; j++) {
// 				EMatrixSolve[i][j] = 0.0;
// 			}
// 			elePhiGuess[i] = phi[bzmesh_process[e].IEN[i]];		// elePhiGuess
// 		}

// 		// loop through gaussian quadrature points
// 		for (int i = 0; i < Gpt.size(); i++) {
// 			for (int j = 0; j < Gpt.size(); j++) {				
// 				ElementEvaluationAll_phi(nen, pre_Nx[ind], pre_dNdx[ind], elePhiGuess, vars);
				
// 				// loop through control points
// 				for (int m = 0; m < nen; m++) {
// 					// terma2 = - eleEP * eleEP * (dPGdx * dNdx[m][0] + dPGdy * dNdx[m][1]);
// 					// termadx = - dAdx * eleEEP * dPGdy * dNdx[m][0];
// 					// termady = - dAdy * eleEEP * dPGdx * dNdx[m][1];
// 					// termdbl = (- elePG * elePG * elePG + (1 - C1) * elePG * elePG + C1 * elePG) * Nx[m];
// 					// EVectorSolve[m] += (elePG * Nx[m] - dt * M_phi * (terma2 - termadx + termady + termdbl) - eleP * Nx[m]) * detJ;
			
// 					EVectorSolve[m] += (vars[0] * pre_Nx[ind][m] - dt * pre_eleMp[ind] *\
// 						((- pre_eleEP[ind] * pre_eleEP[ind] * (vars[1] * pre_dNdx[ind][m][0] + vars[2] * pre_dNdx[ind][m][1])) -\
// 						(- pre_dAdx[ind] * pre_eleEEP[ind] * vars[2] * pre_dNdx[ind][m][0]) +\
// 						(- pre_dAdy[ind] * pre_eleEEP[ind] * vars[1] * pre_dNdx[ind][m][1]) +\
// 						(- vars[0] * vars[0] * vars[0] + (1 - pre_C1[ind]) * vars[0] * vars[0] + pre_C1[ind] * vars[0]) * pre_Nx[ind][m])
// 						- pre_eleP[ind] * pre_Nx[ind][m]) * pre_detJ[ind];

// 					// loop through 16 control points
// 					for (int n = 0; n < nen; n++) {
// 						// terma2 = - eleEP * eleEP * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1]);
// 						// termadx = - dAdx * eleEEP * dNdx[m][1] * dNdx[n][0];
// 						// termady = - dAdy * eleEEP * dNdx[m][0] * dNdx[n][1];
// 						// termdbl = (- 3 * elePG * elePG + 2 * (1 - C1) * elePG + C1 * Nx[m]) * Nx[n];
// 						// EMatrixSolve[m][n] += (Nx[m] * Nx[n] - dt * M_phi * (terma2 - termadx + termady + termdbl)) * detJ;
						
// 						EMatrixSolve[m][n] += (pre_Nx[ind][m] * pre_Nx[ind][n] - dt * pre_eleMp[ind] *
// 							((- pre_eleEP[ind] * pre_eleEP[ind] * (pre_dNdx[ind][m][0] * pre_dNdx[ind][n][0] + pre_dNdx[ind][m][1] * pre_dNdx[ind][n][1])) // terma2
// 							- (- pre_dAdx[ind] * pre_eleEEP[ind] * pre_dNdx[ind][m][1] * pre_dNdx[ind][n][0]) // termadx
// 							+ (- pre_dAdy[ind] * pre_eleEEP[ind] * pre_dNdx[ind][m][0] * pre_dNdx[ind][n][1]) // termady
// 							+ (- 3 * vars[0] * vars[0] + 2 * (1 - pre_C1[ind]) * vars[0] + pre_C1[ind] * pre_Nx[ind][m]) * pre_Nx[ind][n]) // termdbl
// 							) * pre_detJ[ind];
// 						// EMatrixSolve[m][n] += (pre_Nx[ind][m] * pre_Nx[ind][n]) * pre_detJ[ind];
// 					}
// 				}
// 				ind += 1; // incrementing index for extracting pre-calculated variables
// 			}
// 		}

// 		/*Apply Boundary Condition*/
// 		for (int i = 0; i < nen; i++) {
// 			int A = bzmesh_process[e].IEN[i];
// 			if (cpts[A].label == 1) { // domain boundary
// 				ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
// 			}
// 		}

// 		ResidualAssembly(EVectorSolve, bzmesh_process[e].IEN, GR_phi);
// 		MatrixAssembly(EMatrixSolve, bzmesh_process[e].IEN, GK_phi);
// 	}

// 	if (comRank == 0) {
// 		check_itr += 1;
// 	}

// 	VecAssemblyBegin(GR_phi);
// 	MatAssemblyBegin(GK_phi, MAT_FINAL_ASSEMBLY);
// }

// void NeuronGrowth::BuildLinearSystemProcessNG_syn(const vector<Element2D> &tmesh, const vector<Vertex2D> &cpts)
// {
// 	/*Build linear system in each process*/
// 	for (int e = 0; e < bzmesh_process.size(); e++) {
// 		float detJ;
// 		float dudx[2][2];
// 		int nen = bzmesh_process[e].IEN.size();

// 		vector<vector<float>> EMatrixSolve;
// 		vector<float> EVectorSolve;
// 		vector<float> Nx;
// 		vector<array<float, 2>> dNdx;

// 		EMatrixSolve.clear();
// 		EVectorSolve.clear();
// 		Nx.clear();
// 		dNdx.clear();

// 		EMatrixSolve.resize(nen * 1);
// 		EVectorSolve.resize(nen * 1);
// 		Nx.resize(nen);
// 		dNdx.resize(nen);

// 		for (int i = 0; i < nen; i++) {
// 			EMatrixSolve[i].resize(nen * 1, 0.0);
// 		}

// 		for (int i = 0; i < nen * 1; i++) {
// 			for (int j = 0; j < nen * 1; j++)
// 			{
// 				EMatrixSolve[i][j] = 0.0;
// 			}
// 			EVectorSolve[i] = 0.0;
// 		}

// 		vector<float> elePhiDiff, eleSyn;
// 		elePhiDiff.resize(nen);
// 		eleSyn.resize(nen);
// 		float elePf, eleS;

// 		for (int i = 0; i < nen; i++) {
// 			elePhiDiff[i] = abs(phi[bzmesh_process[e].IEN[i]] - phi_prev[bzmesh_process[e].IEN[i]]);
// 			eleSyn[i] = syn[bzmesh_process[e].IEN[i]];
// 		}

// 		for (int i = 0; i < Gpt.size(); i++) {
// 			for (int j = 0; j < Gpt.size(); j++) {
// 				BasisFunction(Gpt[i], Gpt[j], bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dudx, detJ);
// 				detJ = wght[i] * wght[j] * detJ;
// 				ElementValue(Nx, elePhiDiff, elePf);
// 				ElementValue(Nx, eleSyn, eleS);
// 				if (judge_syn == 0)
// 					Tangent_syn(nen, Nx, dNdx, detJ, EMatrixSolve);
// 				Residual_syn(nen, Nx, dNdx, detJ, elePf, eleS, EVectorSolve);
// 			}
// 		}

// 		/*Apply Boundary Condition*/
// 		for (int i = 0; i < nen; i++) {
// 			int A = bzmesh_process[e].IEN[i];
// 			if (cpts[A].label == 1)
// 				ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
// 		}

// 		ResidualAssembly(EVectorSolve, bzmesh_process[e].IEN, GR_syn);
// 		if (judge_syn == 0) // The matrix is the same, so only need to assembly once
// 			MatrixAssembly(EMatrixSolve, bzmesh_process[e].IEN, GK_syn);
// 	}
// 	VecAssemblyBegin(GR_syn);
// 	if (judge_syn == 0)
// 		MatAssemblyBegin(GK_syn, MAT_FINAL_ASSEMBLY);
// }

// void NeuronGrowth::Tangent_syn(const int nen, vector<float> &Nx, vector<array<float, 2>> &dNdx, float detJ, vector<vector<float>> &EMatrixSolve)
// {
// 	/*calculate tangent matrix*/
// 	for (int i = 0; i < nen; i++) {
// 		for (int j = 0; j < nen; j++) {
// 			EMatrixSolve[i][j] += (Nx[i] * Nx[j] + dt * Dc * (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1])) * detJ;
// 			// EMatrixSolve[i][j] += (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1]) * detJ;	// steady state heat transfer
// 		}
// 	}
// }

// void NeuronGrowth::Residual_syn(const int nen, vector<float> &Nx, vector<array<float, 2>> &dNdx, float detJ, float elePf, float eleS, vector<float> &EVectorSolve)
// {
// 	/*calculate residual of the equation*/
// 	for (int i = 0; i < nen; i++) {
// 		EVectorSolve[i] += Nx[i] * (eleS + kappa * elePf) * detJ;
// 		// EVectorSolve[i] += Nx[i] * (eleS) * detJ; // transient heat transfer
// 		// EVectorSolve[i] += 0.00; // steady state heat transfer
// 	}
// }

// void NeuronGrowth::CalculateSumGradPhi0(const vector<Element2D> &tmesh, const vector<Vertex2D> &cpts)
// {
// 	/*Build linear system in each process*/
// 	for (int e = 0; e < bzmesh_process.size(); e++) {
// 		float detJ;
// 		float dudx[2][2];
// 		int nen = bzmesh_process[e].IEN.size();

// 		vector<float> Nx;
// 		vector<array<float, 2>> dNdx;

// 		elePhi0.resize(nen);
// 		for (int i = 0; i < nen; i++) {
// 			elePhi0[i] = phi_0[bzmesh_process[e].IEN[i]];
// 		}

// 		for (int i = 0; i < Gpt.size(); i++) {
// 			for (int j = 0; j < Gpt.size(); j++) {
// 				BasisFunction(Gpt[i], Gpt[j], bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dudx, detJ);
// 				detJ = wght[i] * wght[j] * detJ;
// 				ElementDeriv(nen, dNdx, elePhi0, dP0dx, dP0dy);
// 				sum_grad_phi0_local += pow(dP0dx, 2) + pow(dP0dy, 2);
// 			}
// 		}
// 	}
// }

// void NeuronGrowth::BuildLinearSystemProcessNG_tub(const vector<Element2D> &tmesh, const vector<Vertex2D> &cpts)
// {
// 	/*Build linear system in each process*/
// 	for (int e = 0; e < bzmesh_process.size(); e++) {
// 		float detJ;
// 		float dudx[2][2];
// 		int nen = bzmesh_process[e].IEN.size();

// 		vector<vector<float>> EMatrixSolve;
// 		vector<float> EVectorSolve;
// 		vector<float> Nx;
// 		vector<array<float, 2>> dNdx;

// 		EMatrixSolve.clear();
// 		EVectorSolve.clear();
// 		Nx.clear();
// 		dNdx.clear();

// 		EMatrixSolve.resize(nen);
// 		EVectorSolve.resize(nen);
// 		Nx.resize(nen);
// 		dNdx.resize(nen);

// 		for (int i = 0; i < nen; i++) {
// 			EMatrixSolve[i].resize(nen, 0.0);
// 		}

// 		for (int i = 0; i < nen * 1; i++) {
// 			for (int j = 0; j < nen * 1; j++) {
// 				EMatrixSolve[i][j] = 0.0;
// 			}
// 			EVectorSolve[i] = 0.0;
// 		}

// 		vector<float> elePhi, elePhiPrev, eleConct;
// 		elePhi.resize(nen);
// 		elePhiPrev.resize(nen);
// 		eleConct.resize(nen);
// 		float eleP, dPdx, dPdy, elePprev, eleC, mag_grad_phi0;

// 		for (int i = 0; i < nen; i++) {
// 			elePhi[i] = phi[bzmesh_process[e].IEN[i]];
// 			elePhiPrev[i] = phi_prev[bzmesh_process[e].IEN[i]];
// 			eleConct[i] = tub[bzmesh_process[e].IEN[i]];
// 		}

// 		for (int i = 0; i < Gpt.size(); i++) {
// 			for (int j = 0; j < Gpt.size(); j++) {
// 				BasisFunction(Gpt[i], Gpt[j], bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dudx, detJ);
// 				detJ = wght[i] * wght[j] * detJ;
// 				ElementValue(Nx, elePhi, eleP);
// 				ElementValue(Nx, elePhiPrev, elePprev);
// 				ElementValue(Nx, eleConct, eleC);
// 				ElementDeriv(nen, dNdx, elePhi, dPdx, dPdy);

// 				ElementDeriv(nen, dNdx, elePhi0, dP0dx, dP0dy);
// 				mag_grad_phi0 = pow(dP0dx, 2) + pow(dP0dy, 2);
				
// 				Tangent_tub(nen, Nx, dNdx, detJ, eleC, eleP, dPdx, dPdy, elePprev, mag_grad_phi0, EMatrixSolve);
// 				Residual_tub(nen, Nx, dNdx, detJ, eleC, eleP, elePprev, mag_grad_phi0, EVectorSolve);
// 			}
// 		}

// 		/*Apply Boundary Condition*/
// 		for (int i = 0; i < nen; i++) {
// 			int A = bzmesh_process[e].IEN[i];
// 			if (cpts[A].label == 1)
// 				ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
// 		}

// 		ResidualAssembly(EVectorSolve, bzmesh_process[e].IEN, GR_tub);
// 		// if (judge_tub == 0) // The matrix is the same, so only need to assembly once
// 		MatrixAssembly(EMatrixSolve, bzmesh_process[e].IEN, GK_tub);
// 	}

// 	VecAssemblyBegin(GR_tub);
// 	// if (judge_tub == 0)
// 	MatAssemblyBegin(GK_tub, MAT_FINAL_ASSEMBLY);
// }

// void NeuronGrowth::Tangent_tub(const int nen, vector<float> &Nx, vector<array<float, 2>> &dNdx, float detJ, float eleC, float eleP, float dPdx, float dPdy, float elePprev, float mag_grad_phi0, vector<vector<float>> &EMatrixSolve)
// {
// 	float term_diff, term_alph, term_beta;
// 	/*calculate tangent matrix*/
// 	for (int i = 0; i < nen; i++) {
// 		for (int j = 0; j < nen; j++) {
// 			term_diff = - Diff * (eleP * (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1]));
// 			term_alph = alphaT * (eleP * (dNdx[i][0] + dNdx[i][1]) + Nx[i] * (dPdx + dPdy)) * Nx[j];
// 			term_beta = betaT * (eleP * Nx[i]) * Nx[j];
			
// 			EMatrixSolve[i][j] += (Nx[i] * Nx[j] - dt / eleP * (term_diff - term_alph - term_beta)) * detJ;
// 			// EMatrixSolve[i][j] += ((2 * eleP - elePprev) * Nx[i] * Nx[j] - dt/2 * (term_diff - term_alph - term_beta + term_source)) * detJ;
// 		}
// 	}
// }

// void NeuronGrowth::Residual_tub(const int nen, vector<float> &Nx, vector<array<float, 2>> &dNdx, float detJ, float eleC, float eleP, float elePprev, float mag_grad_phi0, vector<float> &EVectorSolve)
// {
// 	float term_source;
// 	/*calculate residual of the equation*/
// 	for (int i = 0; i < nen; i++) {
// 		term_source = source_coeff * mag_grad_phi0 / sum_grad_phi0_global;
// 		EVectorSolve[i] += (dt / eleP * term_source + eleC) * Nx[i] * detJ;
// 		// EVectorSolve[i] += eleP * eleC * Nx[i] * detJ;
// 	}
// }

// void NeuronGrowth::BuildLinearSystemProcessNG_syn_tub(const vector<Element2D> &tmesh, const vector<Vertex2D> &cpts)
// {
// 	/*Build linear system in each process*/
// 	int ind(0);
// 	for (int e = 0; e < bzmesh_process.size(); e++) {
// 		int nen = bzmesh_process[e].IEN.size();

// 		vector<vector<float>> EMatrixSolve_syn, EMatrixSolve_tub;
// 		vector<float> EVectorSolve_syn, EVectorSolve_tub;

// 		EMatrixSolve_syn.clear();
// 		EVectorSolve_syn.clear();
// 		EMatrixSolve_tub.clear();
// 		EVectorSolve_tub.clear();

// 		EMatrixSolve_syn.resize(nen * 1);
// 		EVectorSolve_syn.resize(nen * 1);
// 		EMatrixSolve_tub.resize(nen);
// 		EVectorSolve_tub.resize(nen);

// 		vector<vector<float>> eleVal_st;
// 		eleVal_st.clear();
// 		eleVal_st.resize(5);

// 		for (int i = 0; i < nen * 1; i++) {
// 			EMatrixSolve_syn[i].resize(nen * 1, 0.0);
// 			EMatrixSolve_tub[i].resize(nen, 0.0);
// 			for (int j = 0; j < nen * 1; j++)
// 			{
// 				EMatrixSolve_syn[i][j] = 0.0;
// 				EMatrixSolve_tub[i][j] = 0.0;
// 			}
// 			EVectorSolve_syn[i] = 0.0;
// 			EVectorSolve_tub[i] = 0.0;

// 			if (i < 5) 
// 				eleVal_st[i].resize(nen);
// 		}
		
// 		for (int i = 0; i < nen * 1; i++) {
// 			eleVal_st[0][i] = phi[bzmesh_process[e].IEN[i]] - phi_prev[bzmesh_process[e].IEN[i]];		// elePhiDiff
// 			eleVal_st[1][i] = syn[bzmesh_process[e].IEN[i]];						// eleSyn
// 			eleVal_st[2][i] = phi[bzmesh_process[e].IEN[i]];						// elePhi
// 			eleVal_st[3][i] = phi_prev[bzmesh_process[e].IEN[i]];						// elePhiPrev
// 			eleVal_st[4][i] = tub[bzmesh_process[e].IEN[i]];						// eleConct
// 		}


// 		for (int i = 0; i < Gpt.size(); i++) {
// 			for (int j = 0; j < Gpt.size(); j++) {
				
// 				vector<float> vars_st;
// 				vars_st.clear();
// 				vars_st.resize(12);
// 				//        0      1     2      3	  4      5        6          7          8          9         10          11
// 				// float elePf, eleS, eleP, dPdx, dPdy, elePprev, eleC, mag_grad_phi0, term_diff, term_alph, term_beta, term_source;

// 				float detJ;
// 				vector<float> Nx;
// 				vector<array<float, 2>> dNdx;
// 				Nx.clear();
// 				dNdx.clear();
// 				Nx.resize(nen);
// 				dNdx.resize(nen);

// 				Nx = pre_Nx[ind];
// 				dNdx = pre_dNdx[ind];
// 				detJ = pre_detJ[ind];
// 				vars_st[7] = pre_mag_grad_phi0[ind];
// 				vars_st[11] = pre_term_source[ind];
// 				// vars_st[11] = 0.5;
// 				ind += 1;

// 				ElementEvaluationAll_syn_tub(nen, Nx, dNdx, eleVal_st, vars_st);

// 				for (int m = 0; m < nen; m++) {
// 					EVectorSolve_syn[m] += Nx[m] * (vars_st[1] + kappa * vars_st[0]) * detJ;

// 					// EVectorSolve_tub[m] += (dt / vars_st[2] * vars_st[11] + vars_st[6]) * Nx[m] * detJ;
// 					// EVectorSolve_tub[m] += (dt / vars_st[2] + vars_st[6]) * Nx[m] * detJ;
// 					EVectorSolve_tub[m] += (vars_st[2] * vars_st[6] + (dt/4) * vars_st[11]) * Nx[m] * detJ;

// 					for (int n = 0; n < nen; n++) {
// 						if (judge_syn == 0)
// 							EMatrixSolve_syn[m][n] += (Nx[m] * Nx[n] + dt/4 * Dc * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1])) * detJ;

// 						// EMatrixSolve_tub[m][n] += (Nx[m] * Nx[n] - dt / vars_st[2] * (
// 						// 	(- Diff * (vars_st[2] * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1])))
// 						// 	- (alphaT * (vars_st[2] * (dNdx[m][0] + dNdx[n][1]) + Nx[m] * (vars_st[3] + vars_st[4])) * Nx[n])
// 						// 	- (betaT * (vars_st[2] * Nx[m]) * Nx[n]))
// 						// 	) * detJ;

// 						EMatrixSolve_tub[m][n] += ((vars_st[0] * Nx[m] + vars_st[2] * Nx[m]) * Nx[n] - (dt/4) * (
// 							(- Diff * (vars_st[2] * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1])))
// 							- (alphaT * (vars_st[2] * (dNdx[m][0] + dNdx[n][1]) + Nx[m] * (vars_st[3] + vars_st[4])) * Nx[n])
// 							- (betaT * (vars_st[2] * Nx[m]) * Nx[n]))
// 							) * detJ;

// 						// EMatrixSolve_tub[m][n] += (Nx[m] * Nx[n] - dt * (
// 						// 	(- Diff * (vars_st[2] * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1])))
// 						// 	- (alphaT * (vars_st[2] * (dNdx[m][0] + dNdx[n][1]) + Nx[m] * (vars_st[3] + vars_st[4])) * Nx[n])
// 						// 	- (betaT * (vars_st[2] * Nx[m]) * Nx[n]))
// 						// 	) * detJ;
						
// 					}
// 				}
// 			}
// 		}

// 		/*Apply Boundary Condition*/
// 		for (int i = 0; i < nen; i++) {
// 			int A = bzmesh_process[e].IEN[i];
// 			if (cpts[A].label == 1) // domain boundary
// 				ApplyBoundaryCondition(0, i, 0, EMatrixSolve_syn, EVectorSolve_syn);
// 			// if (round(phi_0[A]) > 0)
// 			// 	ApplyBoundaryCondition(tub[A] + 2e-5 * tub[A], i, 0, EMatrixSolve_tub, EVectorSolve_tub);
// 				// ApplyBoundaryCondition(tub_0[A], i, 0, EMatrixSolve_tub, EVectorSolve_tub);
// 			if (round(phi[A]) == 0) // outside soma bc for tubulin
// 				ApplyBoundaryCondition(0, i, 0, EMatrixSolve_tub, EVectorSolve_tub);
// 		}

// 		ResidualAssembly(EVectorSolve_syn, bzmesh_process[e].IEN, GR_syn);
// 		if (judge_syn == 0) // The matrix is the same, so only need to assembly once
// 			MatrixAssembly(EMatrixSolve_syn, bzmesh_process[e].IEN, GK_syn);

// 		ResidualAssembly(EVectorSolve_tub, bzmesh_process[e].IEN, GR_tub);
// 		MatrixAssembly(EMatrixSolve_tub, bzmesh_process[e].IEN, GK_tub);
// 	}

// 	VecAssemblyBegin(GR_syn);
// 	if (judge_syn == 0)
// 		MatAssemblyBegin(GK_syn, MAT_FINAL_ASSEMBLY);

// 	VecAssemblyBegin(GR_tub);
// 	MatAssemblyBegin(GK_tub, MAT_FINAL_ASSEMBLY);
// }

// int NeuronGrowth::CheckExpansion(vector<float> input) {

// 	int length = input.size();
// 	int sz = sqrt(length);

// 	// 0 - left | 1 - top | 2 - right | 3 - bottom | 4 - no action
// 	for (int i = 0; i < 2; i++) {
// 		for (int j = 1; j < (sz-1); j++) {
// 			if (round(input[i * sz + j]) > 0) {
// 				return 0;
// 			} else if (round(input[j * sz + i]) > 0) {
// 				return 3;
// 			} else if (round(input[(sz - i - 1) * sz + j]) > 0) {
// 				return 2;
// 			} else if (round(input[j * sz - i - 1]) > 0) {
// 				return 1;
// 			}
// 		}
// 	}
// 	return 4;
// }

// int NeuronGrowth::CheckExpansion(vector<float> input, int NX, int NY) 
// {
// 	// 0 - left | 1 - top | 2 - right | 3 - bottom | 4 - no action
// 	for (int i = 0; i < 10; i++) {
// 		for (int j = 1; j < NY; j++) {
// 			if (round(input[i * (NY+1) + j]) > 0) {
// 				return 0;
// 			} else if (round(input[(NX+1 - i - 1) * (NY+1) + j]) > 0) {
// 				return 2;
// 			}
// 		}
// 		for (int j = 1; j < NX; j++) {
// 			if (round(input[j * (NY+1) + i]) > 0) {
// 				return 3;
// 			} else if (round(input[j * (NY+1) - i]) > 0) {
// 				return 1;
// 			}
// 		}
// 	}
// 	return 4;
// }

// void NeuronGrowth::ExpandDomain(vector<float> input, vector<float> &expd_var, int edge) 
// {
// 	int length = input.size();
// 	int sz = sqrt(length);

// 	int expd_sz = 10; // directional expanding size
// 	int new_sz = sz + expd_sz;

// 	expd_var.clear(); expd_var.resize(pow(new_sz,2));
// 	for (int i = 0; i < new_sz; i++)
// 		expd_var[i] = 0;
	
// 	int ind;
// 	switch (edge) {
// 		case 0: // left
// 			ind = 0;
// 			for (int i = expd_sz; i < new_sz-1; i++) {
// 				for (int j = (int)(expd_sz/2); j < (new_sz - (int)(expd_sz/2)); j++) {
// 					expd_var[i * new_sz + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding left!" << std::endl;
// 			break;
// 		case 1: // top
// 			ind = 0;
// 			for (int i = (int)(expd_sz/2); i < (new_sz - (int)(expd_sz/2)); i++) {
// 				for (int j = 0; j < new_sz-expd_sz; j++) {
// 					expd_var[i * new_sz + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding top!" << std::endl;
// 			break;
// 		case 2: // right
// 			ind = 0;
// 			for (int i = 0; i < new_sz-1-expd_sz; i++) {
// 				for (int j = (int)(expd_sz/2); j < (new_sz - (int)(expd_sz/2)); j++) {
// 					expd_var[i * new_sz + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding right!" << std::endl;
// 			break;
// 		case 3: // bottom
// 			ind = 0;
// 			for (int i = (int)(expd_sz/2); i < (new_sz - (int)(expd_sz/2)); i++) {
// 				for (int j = expd_sz; j < new_sz; j++) {
// 					expd_var[i * new_sz + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding bottom!" << std::endl;
// 			break;
// 		case 4: // all direction
// 			ind = 0;
// 			for (int i = (int)(expd_sz/2); i < (new_sz - (int)(expd_sz/2)); i++) {
// 				for (int j = (int)(expd_sz/2); j < (new_sz - (int)(expd_sz/2)); j++) {
// 					expd_var[i * new_sz + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding all direction!" << std::endl;
// 			break;
// 	}
// 	// input.swap(expd_var);
// }

// void NeuronGrowth::ExpandDomain(vector<float> input, vector<float> &expd_var, int edge, int NX, int NY) 
// {
// 	int expd_sz = 10; // directional expanding size
	
// 	expd_var.clear(); expd_var.resize((NX+1) * (NY+1));
// 	for (int i = 0; i < (expd_var.size()); i++)
// 		expd_var[i] = 0;
	
// 	// (0-left|1-top|2-right|3-bottom)
// 	int ind;
// 	switch (edge) {
// 		case 0: // left
// 			ind = 0;
// 			for (int i = expd_sz; i < NX; i++) {
// 				for (int j = 0; j <= NY; j++) {
// 					expd_var[i * (NY+1) + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding left!" << std::endl;
// 			break;
// 		case 1: // top
// 			ind = 0;
// 			for (int i = 0; i < NX; i++) {
// 				for (int j = 0; j <= NY-expd_sz; j++) {
// 					expd_var[i * (NY+1) + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding top!" << std::endl;
// 			break;
// 		case 2: // right
// 			ind = 0;
// 			for (int i = 0; i < NX-expd_sz; i++) {
// 				for (int j = 0; j <= NY; j++) {
// 					expd_var[i * (NY+1) + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding right!" << std::endl;
// 			break;
// 		case 3: // bottom
// 			ind = 0;
// 			for (int i = 0; i < NX; i++) {
// 				for (int j = expd_sz; j <= NY; j++) {
// 					expd_var[i * (NY+1) + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding bottom!" << std::endl;
// 			break;
// 	}
// }

// void NeuronGrowth::PopulateRandom(vector<float> &input) {

// 	int length = input.size();
// 	for (int i = 0; i < length; i++) {
// 		if (input[i] == 0)
// 			input[i] = (float)(rand()%100)/(float)100;
// 	}
// }


// // Function to interpolate values for new control points
// std::vector<float> NeuronGrowth::InterpolateValues_closest(
// 	const std::vector<float>& phi_in,
// 	const std::vector<Vertex2D>& cpt,
// 	const std::vector<Vertex2D>& cpt_out) {

// 	std::vector<float> interpolatedValues(cpt_out.size());

// 	for (size_t i = 0; i < cpt_out.size(); ++i) {
// 		float minDistance = std::numeric_limits<float>::max();
// 		size_t closestIndex = 0;

// 		// Find the closest point in cpt for each point in cpt_out
// 		for (size_t j = 0; j < cpt.size(); ++j) {
// 			float distance = distance_d(cpt_out[i], cpt[j]);
// 			if (distance < minDistance) {
// 				minDistance = distance;
// 				closestIndex = j;
// 			}
// 		}

// 		// Assume the value of the closest point
// 		// if (abs(round(phi_in[closestIndex])) >= 0.5) {
// 		// 	interpolatedValues[i] = 1;
// 		// } else {
// 		// 	interpolatedValues[i] = 0;
// 		// }
// 		interpolatedValues[i] = phi_in[closestIndex];
// 	}

// 	return interpolatedValues;
// }

// vector<float> NeuronGrowth::InterpolateValues_closest(const vector<float>& input, const KDTree& kdTree, const vector<Vertex2D>& cpt_out) {

// 	std::vector<float> interpolatedValues(cpt_out.size());
// 	for (size_t i = 0; i < cpt_out.size(); ++i) {
// 		const float query_pt[2] = {cpt_out[i].coor[0], cpt_out[i].coor[1]};
// 		size_t closestIndex;
// 		float out_dist_sqr;
// 		nanoflann::KNNResultSet<float> resultSet(1);
// 		resultSet.init(&closestIndex, &out_dist_sqr);

// 		// Performing the search with correct query points
// 		kdTree.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(0));

// 		// Validate closestIndex is within bounds
// 		if(closestIndex >= 0 && closestIndex < input.size()) {
// 			interpolatedValues[i] = input[closestIndex];
// 			// CellBoundary(input[closestIndex], 0.25);
// 		} else {
// 			// Handle error or unexpected case
// 			interpolatedValues[i] = 0; // Define defaultValue appropriately
// 			PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck0)! x: %f y: %f\n", cpt_out[i].coor[0], cpt_out[i].coor[1]);

// 		}
// 	}

// 	return interpolatedValues;
// }

// vector<float> NeuronGrowth::InterpolateVars_coarse1(vector<float> input, vector<Vertex2D> cpts_initial, const KDTree& kdTree_initial, const vector<Vertex2D>& cpts, int type, int isTheta) 
// {	
// 	vector<float> output;
// 	output.resize(cpts.size(), 0);
	
// 	for (int i = 0; i < cpts.size(); i++) {
// 		float x = round5(cpts[i].coor[0]);
// 		float y = round5(cpts[i].coor[1]);
// 		int ind;
// 		if (KD_SearchPair(cpts_initial, kdTree_initial, x, y, ind)) {
// 			output[i] = input[ind];
// 		} 
// 		else {
// 			int indDown, indUp, indLeft, indRight;
// 			if ((abs(remainder(x,2)) == 1) && (abs(remainder(y,2)) != 1)) {
// 				if (KD_SearchPair(cpts_initial, kdTree_initial, floorf(x)-1, round5(y), indDown) &&
// 					KD_SearchPair(cpts_initial, kdTree_initial, floorf(x)+1, round5(y), indUp)) {
// 					if (type == 0) {
// 						output[i] = max(input[indDown], input[indUp]);
// 					} else if (type == 1) {
// 						output[i] = (input[indDown] + input[indUp])/2;
// 					} else if (type == 2) {
// 						output[i] = 0;
// 					}
// 				} else {
// 					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck0)! x: %f y: %f\n", x, y);
// 			}
// 			} else if ((abs(remainder(x,2)) != 1) && (abs(remainder(y,2)) == 1)) {
// 				if (KD_SearchPair(cpts_initial, kdTree_initial, round5(x), floorf(y)-1, indLeft) &&
// 					KD_SearchPair(cpts_initial, kdTree_initial, round5(x), floorf(y)+1, indRight)) {
// 					if (type == 0) {
// 						output[i] = max(input[indLeft], input[indRight]);
// 					} else if (type == 1) {
// 						output[i] = (input[indLeft] + input[indRight])/2;
// 					} else if (type == 2) {
// 						output[i] = 0;
// 					}
// 				} else {
// 					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck1)! x: %f y: %f\n", x, y);
// 				}
// 			} else if ((abs(remainder(x,2)) == 1) && (abs(remainder(y,2)) == 1)) {
// 				if (KD_SearchPair(cpts_initial, kdTree_initial, floorf(x)-1, floorf(y)-1, indDown) &&
// 					KD_SearchPair(cpts_initial, kdTree_initial, floorf(x)+1, floorf(y)+1, indUp) &&
// 					KD_SearchPair(cpts_initial, kdTree_initial, floorf(x)-1, floorf(y)+1, indLeft) &&
// 					KD_SearchPair(cpts_initial, kdTree_initial, floorf(x)+1, floorf(y)-1, indRight)) {
// 					if (type == 0) {
// 						output[i] = max(max(input[indDown], input[indUp]), max(input[indLeft], input[indRight]));
// 					} else if (type == 1) {
// 						output[i] = (input[indDown] + input[indUp] + input[indLeft] + input[indRight])/4;
// 					} else if (type == 2) {
// 						output[i] = 0;
// 					}
// 				} else {
// 					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck2)! x: %f y: %f\n", x, y);
// 				}
// 			} else {
// 				if (isTheta != 1) {
// 					output[i] = 0;
// 				} else {
// 					output[i] = (float)(rand()%100)/(float)100; // for theta
// 				}
// 			}
// 		}			
// 	}
// 	return output;
// }

// bool NeuronGrowth::KD_SearchPair(const vector<Vertex2D> prev_cpts, const KDTree& kdTree, float targetX, float targetY, int &ind) {

// 	const float query_pt[2] = {targetX, targetY};
// 	size_t closestIndex;
// 	float out_dist_sqr;
// 	nanoflann::KNNResultSet<float> resultSet(1);
// 	resultSet.init(&closestIndex, &out_dist_sqr);

// 	// Performing the search with correct query points
// 	kdTree.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(0));

// 	float x = prev_cpts[closestIndex].coor[0];
// 	float y = prev_cpts[closestIndex].coor[1];
	
// 	if ((max(abs(x-targetX), abs(y-targetY)) <= 1) || (targetX < prev_min_x) || (targetY < prev_min_y) || (targetX > prev_max_x) || (targetY > prev_max_y)) {
// 		ind = static_cast<int>(closestIndex);
// 		return true; // Found the pair (targetX, targetY) in the vector
// 	} else {
// 		PetscPrintf(PETSC_COMM_WORLD, "Failed search! x: %f y: %f | %f %f\n", x, y, targetX, targetY);
// 		return false; // Pair not found in the vector
// 	}
// }

// float NeuronGrowth::RmOutlier(vector<float> &data) 
// {
// 	float sum = 0.0, mean, standardDeviation = 0.0;

// 	for(int i = 0; i < data.size(); i++) {
// 		sum += data[i];
// 	}

// 	mean = sum / data.size();

// 	for(int i = 0; i < data.size(); i++) {
// 		standardDeviation += pow(data[i] - mean, 2);
// 	}

// 	standardDeviation = sqrt(standardDeviation / data.size());

// 	// // std::cout << mean + 2*standardDeviation << " " << maxVal << std::endl;
// 	// for(int i = 0; i < data.size(); i++) {
// 	// 	if (data[i] > (mean + 3.29*standardDeviation)) {
// 	// 		// data[i] = mean + 3.29*standardDeviation;
// 	// 	} else {
// 	// 		data[i] = 0;
// 	// 	}
// 	// }
// 	return (mean + 3.29*standardDeviation);
// }

// float NeuronGrowth::CellBoundary(float phi, float threshold) {
// 	float P;
// 	if (phi > threshold) {
// 		P = 1;
// 	} else {
// 		P = 0;
// 	}
// 	return P;
// }

// // Function to apply a simple smoothing operation to a 2D binary variable
// std::vector<float> NeuronGrowth::SmoothBinary2D(const std::vector<float>& binaryData, int rows, int cols) {
// 	// Create a new vector to store the smoothed data
// 	std::vector<float> smoothedData(rows * cols, 0.0f);

// 	// Define a 5x5 smoothing kernel
// 	const std::vector<float> kernel = {
// 		1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
// 		1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
// 		1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
// 		1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
// 	1.0f, 1.0f, 1.0f, 1.0f, 1.0f
// 	};

// 	// Apply the smoothing operation using the defined kernel
// 	for (int i = 2; i < rows - 2; ++i) {
// 		for (int j = 2; j < cols - 2; ++j) {
// 			float sum = 0.0f;
// 			for (int k = -2; k <= 2; ++k) {
// 				for (int l = -2; l <= 2; ++l) {
// 					sum += binaryData[(i + k) * cols + (j + l)] * kernel[(k + 2) * 5 + (l + 2)];
// 				}
// 			}
// 			smoothedData[i * cols + j] = (sum > 12.0f) ? 1.0f : 0.0f; // Apply threshold to convert sum to binary value
// 		}
// 	}

// 	return smoothedData;
// }


// // Check if point coordinates are within the bounds defined by center and dx dy dz
// bool NeuronGrowth::isInBox(const Vertex2D& point, const Vertex2D& center, float dx, float dy) {
// 	return point.coor[0] >= (center.coor[0] - dx/2) && point.coor[0] <= (center.coor[0] + dx/2) &&
// 		point.coor[1] >= (center.coor[1] - dy/2) && point.coor[1] <= (center.coor[1] + dy/2);
// }

// /**
//  * Calculates the sum of phi values around each vertex in a given set, applying
//  * a threshold to identify significant points, potentially indicating neuron growth tips.
//  * 
//  * This function iterates over a collection of vertices (`cpts`), summing the phi values
//  * within a defined vicinity around each vertex. The vicinity is determined by the `dx`, and `dy`
//  * parameters, which define the dimensions of a box centered on each vertex. Points within
//  * this box contribute to the sum. After calculating the sums, the function applies a threshold
//  * to these values, normalizing them against the highest value found. Points with summed values
//  * above this threshold are marked as potential growth tips by setting their corresponding value
//  * in the output vector to 1; all others are set to 0.
//  * 
//  * @param cpts A vector of Vertex2D objects representing the neuron's vertices.
//  * @param dx The delta in the x-direction to define the vicinity around a point.
//  * @param dy The delta in the y-direction to define the vicinity around a point.
//  * @return std::vector<float> A vector of the same size as `cpts`, where each element is either 0
//  *         (indicating the corresponding vertex is not a tip) or 1 (indicating a potential tip),
//  *         based on the thresholding of summed phi values.
//  * 
//  * Note: The function assumes that `phi` is a pre-defined vector accessible within the class
//  *       that contains phi values corresponding to each Vertex3D in `cpts`. The function
//  *       `isInBox` checks whether a point is within the specified vicinity of another,
//  *       and `CellBoundary` computes a boundary-related value for a given phi, with
//  *       the second argument presumably allowing for further customization.
//  */
// vector<float> NeuronGrowth::calculatePhiSum(const std::vector<Vertex2D>& cpts, float dx, float dy, vector<float> id) {
// 	std::vector<float> tp(cpts.size(), 0); // Initialize the result vector with zeros
// 	tp.clear(); tp.resize(cpts.size());
// 	float threshold = 0.9; // Threshold for filtering tp
// 	float maxVal = 0; // Track the maximum value of tp for normalization

// 	// Iterate over each center point
// 	for (int i = 0; i < cpts.size(); ++i) {
// 		const auto& center = cpts[i];

// 		tp[i] = 0;

// 		// Sum phi values for points within the box centered at 'center'
// 		for (int j = 0; j < phi.size(); ++j) {
// 			if (isInBox(cpts[j], center, dx, dy)) {
// 				tp[i] += CellBoundary(phi[j], 0.5), id[j];
// 			}
// 		}

// 		tp[i] = CellBoundary(phi[i], 0.5) / tp[i] * CellBoundary(phi[i], 0.5);
// 		// Handle NaN cases, ensuring a valid tp value
// 		if (std::isnan(tp[i])) {
// 			tp[i] = 0;
// 		}
		
// 		if ((center.coor[0] <= min_x+1) || (center.coor[0] >= max_x-1) 
// 			|| (center.coor[1] <= min_y+1) || (center.coor[1] >= max_y-1)) {
// 			tp[i] = 0;
// 		}

// 		// Update maxVal to the highest tp value found
// 		if (CellBoundary(phi[i], 0.5) > 0)
// 			maxVal = max(tp[i], maxVal);
// 	}

// 	for (int i = 0; i < tp.size(); i++) {
// 		tp[i] = tp[i]/maxVal;
// 	}

// 	// maxVal = RmOutlier(tp);

// 	// Thresholding and setting tp values based on the maximum value found
// 	for (int i = 0; i < tp.size(); ++i) {
// 		// if (tp[i] > (threshold * maxVal)) {
// 		if (tp[i] > (threshold)) {
// 		// if (tp[i] > (0.018)) {
// 			tp[i] = 1; // Set tp below threshold to 0
// 		} else {
// 			tp[i] = 0; // Set tp above threshold to 1
// 		}
// 	}
// 	return tp;
// }

// void NeuronGrowth::DetectTipsMulti(const vector<float>& phi_fine, const vector<float>& id, const int& numNeuron, vector<float>& phiSum, const int& NX, const int& NY)
// {
// 	// float threshold(0.9997), maxVal(0);
// 	float threshold(0.9), stdVal(0), maxVal(0);
// 	int ind, length((NX+1)*(NY+1));
// 	phiSum.clear();
// 	phiSum.resize(length, 0);

// 	for (int i = (5*NY+5); i < (length-4*NY-4); i++) {
// 		// phiSum[i] = 0;
// 		if (CellBoundary(phi_fine[i], 0.25) > 0) {
// 		// if (CellBoundary(phi_fine[i], 0.5) > 0) {
// 			// int baseID = static_cast<int>(round(id[i]));
// 			for (int j = -4; j < 5; j++) {
// 				for (int k = -4; k < 5; k++) {
// 			// for (int j = -6; j < 7; j++) {
// 			// 	for (int k = -6; k < 7; k++) {
// 					// for (int l = 0; l < numNeuron; l++) {
// 						// int inID = static_cast<int>(round(id[i+j*(NY+1)+k]));
// 						// std::cout << baseID << " " << inID << std::endl;
// 						if ((round(id[i]) == round(id[i+j*(NY+1)+k])) || (round(id[i+j*(NY+1)+k]) == 0)) {
// 						// if (inID == l+1) { // not the same neuron
// 						// if ((l+1) == static_cast<int>(round(id[i+j*(NY+1)+k]))) {
// 							// phiSum[i] += CellBoundary(phi_fine[i+j*(NY+1)+k], 0.1);
// 							phiSum[i] += CellBoundary(phi_fine[i+j*(NY+1)+k], 0.25);
// 							// phiSum[i] += CellBoundary(phi_fine[i+j*(NY+1)+k], 0.5);
// 						}
// 					// }
// 				}
// 			}	
// 			if (phiSum[i] > 0)
// 				phiSum[i] = CellBoundary(phi_fine[i], 0) / phiSum[i];
// 			if (std::isnan(phiSum[i]))
// 				phiSum[i] = 0;
// 			// if (phiSum[i] > maxVal)
// 			// 	maxVal = phiSum[i];
// 		}
// 	}

// 	stdVal = RmOutlier(phiSum);

// 	for (int i = 0; i < phiSum.size(); i++) {
// 		phiSum[i] = phiSum[i]/stdVal;
// 		maxVal = max(phiSum[i], maxVal);
// 	}

// 	// std::cout << maxVal << std::endl;

// 	for (int i = 1+NY; i < length-NY-1; i++) {
// 		// if (phiSum[i] < (threshold * maxVal)) {
// 		// if (phiSum[i] < (threshold)) {
// 		if (phiSum[i] < (1.3)) {
// 			phiSum[i] = 0;
// 		} 
// 		// else {
// 		// 	tip[i] = 1;
// 		// }
// 	}
// }

// void NeuronGrowth::DetectTipsMulti_new(const std::vector<float>& phi_fine, const std::vector<float>& id, int numNeuron, std::vector<float>& phiSum, int NX, int NY) {
// 	float threshold = 0.9997f;
// 	int length = (NX + 1) * (NY + 1);
// 	phiSum.clear();
// 	phiSum.resize(length, 0.0f);

// 	// Kernel for convolution; adjust size/shape for your specific needs
// 	const int kernelSize = 5; // 5x5 kernel
// 	const int kernelOffset = kernelSize / 2;

// 	for (int x = kernelOffset; x < NX - kernelOffset; ++x) {
// 		for (int y = kernelOffset; y < NY - kernelOffset; ++y) {
// 			int index = x * (NY + 1) + y;
// 			if (CellBoundary(phi_fine[index], 0.25) > 0) {
// 				for (int dx = -kernelOffset; dx <= kernelOffset; ++dx) {
// 					for (int dy = -kernelOffset; dy <= kernelOffset; ++dy) {
// 						int neighborIndex = index + dx * (NY + 1) + dy;
// 						if (static_cast<int>(round(id[neighborIndex])) == numNeuron) {
// 							phiSum[index] += CellBoundary(phi_fine[neighborIndex], 0.1);
// 						}
// 					}
// 				}
// 				phiSum[index] = phiSum[index] > 0 ? CellBoundary(phi_fine[index], 0) / phiSum[index] : 0;
// 			}
// 		}
// 	}

// 	// Assuming RmOutlier modifies phiSum to remove outliers and returns new maxVal
// 	float maxVal = RmOutlier(phiSum);

// 	for (float& value : phiSum) {
// 		value /= maxVal;
// 		if (value < threshold) {
// 			value = 0;
// 		}
// 	}
// }

// // Function to perform Breadth-First Search (BFS) for clustering
// void NeuronGrowth::bfs(const vector<float>& matrix, int rows, int cols, int row, int col,
//          vector<bool>& visited, vector<pair<int, int>>& cluster) {
// 	static const int dx[] = {-1, 0, 1, 0, -1, -1, 1, 1}; // Including diagonal directions
// 	static const int dy[] = {0, -1, 0, 1, -1, 1, -1, 1}; // Including diagonal directions
// 	// static const int dx[] = {-1, 0, 1, 0}; // Including diagonal directions
// 	// static const int dy[] = {0, -1, 0, 1}; // Including diagonal directions
	
// 	queue<pair<int, int>> q;
// 	q.push(make_pair(row, col));
// 	visited[row * cols + col] = true;
// 	cluster.push_back(make_pair(row, col));

// 	while (!q.empty()) {
// 		int r = q.front().first;
// 		int c = q.front().second;
// 		q.pop();

// 		for (int i = 0; i < 8; ++i) { // Iterating through all 8 directions
// 		// for (int i = 0; i < 4; ++i) { // Iterating through all 8 directions
// 			int newRow = r + dx[i];
// 			int newCol = c + dy[i];

// 			if (newRow >= 0 && newRow < rows && newCol >= 0 && newCol < cols &&
// 			matrix[newRow * cols + newCol] != 0 && !visited[newRow * cols + newCol]) {
// 				visited[newRow * cols + newCol] = true;
// 				q.push(make_pair(newRow, newCol));
// 				cluster.push_back(make_pair(newRow, newCol));
// 			}
// 		}
// 	}
// }

// // Function to find connected clusters in the matrix
// vector<vector<pair<int, int>>> NeuronGrowth::FindClusters(const vector<float>& matrix, int rows, int cols) {
//     vector<bool> visited(rows * cols, false);
//     vector<vector<pair<int, int>>> clusters;

// 	for (int i = 0; i < rows; ++i) {
// 		for (int j = 0; j < cols; ++j) {
// 			if (matrix[i * cols + j] != 0 && !visited[i * cols + j]) {
// 				vector<pair<int, int>> cluster;
// 				bfs(matrix, rows, cols, i, j, visited, cluster);

// 				if (!cluster.empty()) {
// 					clusters.push_back(cluster);
// 				}
// 			}
// 		}
// 	}
// 	return clusters;
// }

// // Function to find local maxima within connected clusters in the matrix
// vector<float> NeuronGrowth::FindLocalMaximaInClusters(const vector<float>& matrix, int rows, int cols) {
// 	vector<vector<pair<int, int>>> clusters = FindClusters(matrix, rows, cols);
// 	vector<float> localMaxima(rows * cols, 0.0f);

// 	for (const auto& cluster : clusters) {
// 		float maxVal = -numeric_limits<float>::max();
// 		pair<int, int> maxPos = make_pair(-1, -1);

// 		for (const auto& pos : cluster) {
// 			int row = pos.first;
// 			int col = pos.second;

// 			if (matrix[row * cols + col] > maxVal) {
// 				maxVal = matrix[row * cols + col];
// 				maxPos = pos;
// 			}
// 		}

// 		if (maxPos.first != -1 && maxPos.second != -1) {
// 			localMaxima[maxPos.first * cols + maxPos.second] = 1.0f;
			
// 			// localMaxima[maxPos.first * cols + maxPos.second - cols + 1] = 1.0f;
// 			// localMaxima[maxPos.first * cols + maxPos.second - cols] = 1.0f;
// 			// localMaxima[maxPos.first * cols + maxPos.second - cols - 1] = 1.0f;
// 			// localMaxima[maxPos.first * cols + maxPos.second - 1] = 1.0f;
// 			// localMaxima[maxPos.first * cols + maxPos.second + 1] = 1.0f;
// 			// localMaxima[maxPos.first * cols + maxPos.second + cols + 1] = 1.0f;
// 			// localMaxima[maxPos.first * cols + maxPos.second + cols] = 1.0f;
// 			// localMaxima[maxPos.first * cols + maxPos.second + cols - 1] = 1.0f;
// 		}
// 	}
// 	return localMaxima;
// }

// // Function to find centroids of connected clusters in the matrix
// vector<float> NeuronGrowth::FindCentroidsInClusters(const vector<float>& matrix, int rows, int cols) {
// 	vector<vector<pair<int, int>>> clusters = FindClusters(matrix, rows, cols);
// 	vector<float> centroids(rows * cols, 0.0f);

// 	for (const auto& cluster : clusters) {
// 		if (cluster.empty()) continue;

// 		float sumRow = 0.0f, sumCol = 0.0f;
// 		for (const auto& pos : cluster) {
// 			sumRow += pos.first;
// 			sumCol += pos.second;
// 		}
// 		int centroidRow = static_cast<int>(round(sumRow / cluster.size()));
// 		int centroidCol = static_cast<int>(round(sumCol / cluster.size()));

// 		// Ensure the centroid position is within bounds
// 		centroidRow = std::max(0, std::min(centroidRow, rows - 1));
// 		centroidCol = std::max(0, std::min(centroidCol, cols - 1));

// 		centroids[centroidRow * cols + centroidCol] = 1.0f;

// 		// centroids[centroidRow.first * cols + centroidRow.second - cols + 1] = 1.0f;
// 		// centroids[centroidRow.first * cols + centroidRow.second - cols] = 1.0f;
// 		// centroids[centroidRow.first * cols + centroidRow.second - cols - 1] = 1.0f;
// 		// centroids[centroidRow.first * cols + centroidRow.second - 1] = 1.0f;
// 		// centroids[centroidRow.first * cols + centroidRow.second + 1] = 1.0f;
// 		// centroids[centroidRow.first * cols + centroidRow.second + cols + 1] = 1.0f;
// 		// centroids[centroidRow.first * cols + centroidRow.second + cols] = 1.0f;
// 		// centroids[centroidRow.first * cols + centroidRow.second + cols - 1] = 1.0f;
// 	}
// 	return centroids;
// }

// bool NeuronGrowth::IsLocalMaximum(const vector<float>& matrix, int rows, int cols, int x, int y) {
//     float value = matrix[x * cols + y];
//     for (int dx = -1; dx <= 1; ++dx) {
//         for (int dy = -1; dy <= 1; ++dy) {
//             int nx = x + dx, ny = y + dy;
//             if (nx >= 0 && nx < rows && ny >= 0 && ny < cols && (dx != 0 || dy != 0)) {
//                 if (matrix[nx * cols + ny] > value) {
//                     return false;
//                 }
//             }
//         }
//     }
//     return true;
// }

// vector<float> NeuronGrowth::FindCentroidsOfLocalMaximaClusters(const vector<float>& matrix, int rows, int cols) {
//     vector<vector<pair<int, int>>> clusters = FindClusters(matrix, rows, cols);
//     vector<float> centroids(rows * cols, 0.0f);

//     for (const auto& cluster : clusters) {
//         vector<pair<int, int>> localMaxima;

//         // Find all local maxima in the cluster
//         for (const auto& pos : cluster) {
//             if (IsLocalMaximum(matrix, rows, cols, pos.first, pos.second)) {
//                 localMaxima.push_back(pos);
//             }
//         }

//         // Calculate centroid of local maxima
//         float sumRow = 0, sumCol = 0;
//         for (const auto& maxPos : localMaxima) {
//             sumRow += maxPos.first;
//             sumCol += maxPos.second;
//         }

//         int centroidRow = static_cast<int>(round(sumRow / localMaxima.size()));
//         int centroidCol = static_cast<int>(round(sumCol / localMaxima.size()));

//         centroids[centroidRow * cols + centroidCol] = 1.0f;

// 	// centroids[centroidRow * cols + centroidCol - cols + 1] = 1.0f;
// 	// centroids[centroidRow * cols + centroidCol- cols] = 1.0f;
// 	// centroids[centroidRow * cols + centroidCol - cols - 1] = 1.0f;
// 	// centroids[centroidRow * cols + centroidCol - 1] = 1.0f;
// 	// centroids[centroidRow * cols + centroidCol + 1] = 1.0f;
// 	// centroids[centroidRow * cols + centroidCol + cols + 1] = 1.0f;
// 	// centroids[centroidRow * cols + centroidCol + cols] = 1.0f;
// 	// centroids[centroidRow * cols + centroidCol + cols - 1] = 1.0f;
//     }

//     return centroids;
// }

// vector<vector<int>> NeuronGrowth::ConvertTo2DIntVector(const vector<float> input, int NX, int NY) 
// {
// 	vector<vector<int>> output;

// 	int k = 0;
// 	for (int i = 0; i < NX+1; i++) {
// 		vector<int> row;
// 		for (int j = 0; j < NY+1; j++) {
// 			// row.push_back(CellBoundary(abs(input[k]), 0.01));
// 			row.push_back(CellBoundary(abs(input[k]), 0.25));
// 			// row.push_back(CellBoundary(abs(input[k]), 0.75));
// 			// row.push_back(CellBoundary(input[k], 0.5));
// 			k++;
// 		}
// 		output.push_back(row);
// 	}
// 	return output;
// }

// vector<vector<float>> NeuronGrowth::ConvertTo2DFloatVector(const vector<float> input, int NX, int NY) 
// {
// 	vector<vector<float>> output;
// 	int k = 0;
// 	for (int i = 0; i < NX; i++) {
// 		vector<float> row;
// 		for (int j = 0; j < NY; j++) {
// 			row.push_back(input[k]);
// 			k++;
// 		}
// 		output.push_back(row);
// 	}
// 	return output;
// }

// // Function to perform flood fill
// void NeuronGrowth::FloodFill(vector<vector<int>>& image, int x, int y, int newColor, int originalColor) 
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
// 		FloodFill(image, x + dx[i], y + dy[i], newColor, originalColor);
// 	}
// }

// // // Function to perform flood fill
// // void NeuronGrowth::ExtractPhi(vector<Vertex2D> cpts) 
// // {
// // 	neurons = ConvertTo2DIntVector(phi, NX, NY);   
// // 	for (int i = 0; i < seed.size(); i++) {
// // 		int startX = seed[i][0] - originX;
// // 		int startY = seed[i][1] - originY;
// // 		// std::cout << " startX " << startX << " startY " << startY << std::endl;
// // 		int newColor = i+1;
// // 		int originalColor = neurons[startX][startY];
// // 		FloodFill(neurons, startX, startY, newColor, originalColor);
// // 	}
// // }

// // Function to perform flood fill
// void NeuronGrowth::IdentifyNeurons(vector<float>& phi_in, vector<vector<int>>& neurons, vector<array<float, 2>> seed, int NX, int NY, int originX, int originY) 
// {
// 	neurons = ConvertTo2DIntVector(phi_in, NX, NY);   
// 	// PrintOutNeurons(neurons);
// 	// std::cout << neurons.size() << " " << NX << " " << originX << " " << NY << " " << originY << std::endl;
// 	for (int i = 0; i < seed.size(); i++) {
// 		int startX = (seed[i][0] - originX);
// 		int startY = (seed[i][1] - originY);
// 		// std::cout << " startX " << startX << " startY " << startY << std::endl;
// 		int newColor = i+1;
// 		int originalColor = neurons[startX][startY];
// 		FloodFill(neurons, startX, startY, newColor, originalColor);
// 	}
// }

// // Check if a point is valid in the grid
// bool NeuronGrowth::isValid(int x, int y, int rows, int cols) 
// {
//     return (x >= 0 && x < rows && y >= 0 && y < cols);
// }

// // Function to calculate geodesic distance from a point
// vector<vector<int>> NeuronGrowth::CalculateGeodesicDistanceFromPoint(vector<vector<int>> neurons, int startX, int startY) 
// {
// 	int rows = neurons.size();
// 	if (rows == 0) return {};

// 	int cols = neurons[0].size();

// 	vector<vector<int>> distances(rows, vector<int>(cols, INF));
// 	vector<vector<bool>> visited(rows, vector<bool>(cols, false));

// 	distances[startX][startY] = 0;
// 	visited[startX][startY] = true;

// 	queue<pair<int, int>> q;
// 	q.push({startX, startY});

// 	int dx[] = {-1, 1, 0, 0};
// 	int dy[] = {0, 0, -1, 1};

// 	while (!q.empty()) {
// 		pair<int, int> current = q.front();
// 		q.pop();

// 		int x = current.first;
// 		int y = current.second;

// 		for (int i = 0; i < 4; ++i) {
// 			int newX = x + dx[i];
// 			int newY = y + dy[i];

// 			if (isValid(newX, newY, rows, cols) && neurons[newX][newY] == 1 && !visited[newX][newY]) {
// 				visited[newX][newY] = true;
// 				distances[newX][newY] = distances[x][y] + 1;
// 				q.push({newX, newY});
// 			}
// 		}
// 	}
// 	return distances;
// }

// // Function to calculate geodesic distance from a point
// vector<vector<int>> NeuronGrowth::CalculateGeodesicDistanceFromPoint(vector<vector<int>> neurons, vector<array<float, 2>> &seed, int originX, int originY) 
// {
// 	int rows = neurons.size();
// 	if (rows == 0) return {};

// 	int cols = neurons[0].size();

// 	vector<vector<int>> distances(rows, vector<int>(cols, INF));
// 	vector<vector<bool>> visited(rows, vector<bool>(cols, false));

// 	int startX, startY;
// 	for (int i = 0; i < seed.size(); i++) {
// 		startX = seed[i][0] - originX;
// 		startY = seed[i][1] - originY;

// 		distances[startX][startY] = 0;
// 		visited[startX][startY] = true;

// 		queue<pair<int, int>> q;
// 		q.push({startX, startY});

// 		int dx[] = {-1, 1, 0, 0};
// 		int dy[] = {0, 0, -1, 1};

// 		while (!q.empty()) {
// 			pair<int, int> current = q.front();
// 			q.pop();

// 			int x = current.first;
// 			int y = current.second;

// 			for (int i = 0; i < 4; ++i) {
// 				int newX = x + dx[i];
// 				int newY = y + dy[i];

// 				if (isValid(newX, newY, rows, cols) && neurons[newX][newY] == 1 && !visited[newX][newY]) {
// 					visited[newX][newY] = true;
// 					distances[newX][newY] = distances[x][y] + 1;
// 					q.push({newX, newY});
// 				}
// 			}
// 		}
// 	}
// 	return distances;
// }

// // Function to calculate geodesic distance from a point
// vector<vector<array<int, 2>>> NeuronGrowth::NeuriteTracing(vector<vector<double>> distance) 
// {
// 	vector<vector<array<int, 2>>> traces;
// 	return traces;
// }

// // Function to calculate geodesic distance from a point
// void NeuronGrowth::SaveNGvars(vector<vector<float>> &NGvars, int NX, int NY, string fn) 
// {
// 	// writing initialized variables for debugging purposes
// 	bool visualization = true;
// 	PrintVec2TXT(NGvars[0], fn + "/phi_" + std::to_string(n) + ".txt", visualization);
// 	PrintVec2TXT(NGvars[1], fn + "/syn_" + std::to_string(n) + ".txt", visualization);
// 	PrintVec2TXT(NGvars[2], fn + "/tub_" + std::to_string(n) + ".txt", visualization);
// 	PrintVec2TXT(NGvars[3], fn + "/theta_" + std::to_string(n) + ".txt", visualization);
// 	PrintVec2TXT(NGvars[4], fn + "/phi_0_" + std::to_string(n) + ".txt", visualization);
// 	PrintVec2TXT(NGvars[5], fn + "/tub_0_" + std::to_string(n) + ".txt", visualization);
// }

// void NeuronGrowth::PrintOutNeurons(vector<vector<int>> neurons) 
// {
// 	ierr = PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------------\n");
// 	int dwnRatio = 1;
// 	for (int i = 0; i < neurons.size(); i+=2*dwnRatio) {
// 		ierr = PetscPrintf(PETSC_COMM_WORLD, "| ");
// 		for (int j = 0; j < neurons[i].size(); j+=dwnRatio) {
// 			if (round(neurons[i][j]) == 0) {
// 				ierr = PetscPrintf(PETSC_COMM_WORLD, " ");
// 			} else {
// 				ierr = PetscPrintf(PETSC_COMM_WORLD, "#");
// 			}
// 		}
// 		ierr = PetscPrintf(PETSC_COMM_WORLD, "|\n");
// 	}
// 	ierr = PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------------\n");
// }

// PetscErrorCode FormFunction_phi(SNES snes, Vec x, Vec F, void *ctx)
// {
// 	PetscErrorCode ierr;
// 	Vec P_seq;
// 	PetscScalar *Parray;
// 	VecScatter scatter_ctx1;
// 	ierr = VecScatterCreateToAll(x, &scatter_ctx1, &P_seq); CHKERRQ(ierr);
// 	ierr = VecScatterBegin(scatter_ctx1, x, P_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
// 	ierr = VecScatterEnd(scatter_ctx1, x, P_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
// 	ierr = VecGetArray(P_seq, &Parray); CHKERRQ(ierr);
// 	ierr = VecSet(F, 0.0); CHKERRQ(ierr);

// 	NeuronGrowth *user = (NeuronGrowth *)ctx;

// 	/*Build linear system in each process*/
// 	int ind(0); // pre-calculated variable index
// 	for (int e = 0; e < user->bzmesh_process.size(); e++) {
// 		int nen = user->bzmesh_process[e].IEN.size(); // 16 supporting cp for 2D case

// 		vector<vector<float>> EMatrixSolve;	EMatrixSolve.clear();
// 		vector<float> EVectorSolve;		EVectorSolve.clear();
// 		EMatrixSolve.resize(nen);
// 		EVectorSolve.resize(nen);

// 		user->eleVal.clear();
// 		user->eleVal.resize(8);
// 		for (int i = 0; i < nen * 1; i++) {
// 			// EMatrixSolve[i].reserve(nen*sizeof(float));
// 			EMatrixSolve[i].resize(nen);
// 			if (i < 8)
// 				user->eleVal[i].resize(nen);
// 		}

// 		// preparing element stiffness matrix and vector, extract values from control mesh
// 		for (int i = 0; i < nen; i++) {

// 			EVectorSolve[i] = 0.0;
// 			for (int j = 0; j < nen; j++) {
// 				EMatrixSolve[i][j] = 0.0;
// 			}

// 			user->eleVal[0][i] = Parray[user->bzmesh_process[e].IEN[i]];		// elePhiGuess
// 			user->eleVal[1][i] = user->phi[user->bzmesh_process[e].IEN[i]];		// elePhi
// 			user->eleVal[4][i] = user->theta[user->bzmesh_process[e].IEN[i]];	// eleTheta
// 			user->eleVal[6][i] = 0;							// eleEpsilon
// 			user->eleVal[7][i] = 0;							// eleEpsilonP
// 		}

// 		// loop through gaussian quadrature points
// 		for (int i = 0; i < user->Gpt.size(); i++) {
// 			for (int j = 0; j < user->Gpt.size(); j++) {
// 				user->EvaluateOrientation(nen, user->pre_Nx[ind], user->pre_dNdx[ind], user->eleVal[1], user->eleVal[4], user->eleVal[6], user->eleVal[7]);
// 				user->ElementEvaluationAll_phi(nen, user->pre_Nx[ind], user->pre_dNdx[ind], user->eleVal, user->vars);
// 				// vars  0   1    2     3      4      5     6     7      8     9      10      11      12     13    14    15      16     17     
// 				// float C1, C0, elePG, dPGdx, dPGdy, eleP, eleS, eleTb, eleE, eleTp, dThedx, dThedy, eleEP, dAdx, dAdy, eleEEP, dAPdx, dAPdy;
				
// 				// loop through control points
// 				for (int m = 0; m < nen; m++) {
// 					// terma2 = - eleEP * eleEP * (dPGdx * dNdx[m][0] + dPGdy * dNdx[m][1]);
// 					// termadx = - dAdx * eleEEP * dPGdy * dNdx[m][0];
// 					// termady = - dAdy * eleEEP * dPGdx * dNdx[m][1];
// 					// termdbl = (- elePG * elePG * elePG + (1 - C1) * elePG * elePG + C1 * elePG) * Nx[m];
// 					// EVectorSolve[m] += (elePG * Nx[m] - user->dt * user->M_phi * (terma2 - termadx + termady + termdbl) - eleP * Nx[m]) * detJ;
// 					EVectorSolve[m] += (user->vars[2] * user->pre_Nx[ind][m] - user->dt * user->pre_eleMp[ind] *\
// 						((- user->vars[12] * user->vars[12] * (user->vars[3] * user->pre_dNdx[ind][m][0] + user->vars[4] * user->pre_dNdx[ind][m][1])) -\
// 						(- user->vars[13] * user->vars[15] * user->vars[4] * user->pre_dNdx[ind][m][0]) +\
// 						(- user->vars[14] * user->vars[15] * user->vars[3] * user->pre_dNdx[ind][m][1]) +\
// 						(- user->vars[2] * user->vars[2] * user->vars[2] + (1 - user->pre_C1[ind]) * user->vars[2] * user->vars[2] + user->pre_C1[ind] * user->vars[2]) * user->pre_Nx[ind][m])
// 						- user->vars[5] * user->pre_Nx[ind][m]) * user->pre_detJ[ind];
// 					// loop through 16 control points
// 					for (int n = 0; n < nen; n++) {
// 						// terma2 = - eleEP * eleEP * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1]);
// 						// termadx = - dAdx * eleEEP * dNdx[m][1] * dNdx[n][0];
// 						// termady = - dAdy * eleEEP * dNdx[m][0] * dNdx[n][1];
// 						// termdbl = (- 3 * elePG * elePG + 2 * (1 - C1) * elePG + C1 * Nx[m]) * Nx[n];
// 						// EMatrixSolve[m][n] += (Nx[m] * Nx[n] - user->dt * user->M_phi * (terma2 - termadx + termady + termdbl)) * detJ;
// 						EMatrixSolve[m][n] += (user->pre_Nx[ind][m] * user->pre_Nx[ind][n] - user->dt * user->pre_eleMp[ind] *
// 							((- user->vars[12] * user->vars[12] * (user->pre_dNdx[ind][m][0] * user->pre_dNdx[ind][n][0] + user->pre_dNdx[ind][m][1] * user->pre_dNdx[ind][n][1])) // terma2
// 							- (- user->vars[13] * user->vars[15] * user->pre_dNdx[ind][m][1] * user->pre_dNdx[ind][n][0]) // termadx
// 							+ (- user->vars[14] * user->vars[15] * user->pre_dNdx[ind][m][0] * user->pre_dNdx[ind][n][1]) // termady
// 							+ (- 3 * user->vars[2] * user->vars[2] + 2 * (1 - user->pre_C1[ind]) * user->vars[2] + user->pre_C1[ind] * user->pre_Nx[ind][m]) * user->pre_Nx[ind][n]) // termdbl
// 							) * user->pre_detJ[ind];
// 					}
// 				}
// 				ind += 1; // incrementing index for extracting pre-calculated variables
// 			}
// 		}

// 		/*Apply Boundary Condition*/
// 		for (int i = 0; i < nen; i++) {
// 			int A = user->bzmesh_process[e].IEN[i];
// 			if (user->cpts[A].label == 1) {// domain boundary
// 				user->ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
// 			}
// 		}

// 		user->ResidualAssembly(EVectorSolve, user->bzmesh_process[e].IEN, F);
// 		user->MatrixAssembly(EMatrixSolve, user->bzmesh_process[e].IEN, user->J);
// 	}

// 	ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
// 	ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

// 	ierr = MatAssemblyBegin(user->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
// 	ierr = MatAssemblyEnd(user->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

// 	ierr = VecRestoreArray(P_seq, &Parray); CHKERRQ(ierr);
// 	ierr = VecScatterDestroy(&scatter_ctx1); CHKERRQ(ierr);
// 	ierr = VecDestroy(&P_seq); CHKERRQ(ierr);

// 	return 0;
// }


// PetscErrorCode FormFunction_phi_test(SNES snes, Vec x, Vec F, void *ctx)
// {
// 	PetscErrorCode ierr;
	
// 	/*To get the global x needed for IGA, each element needs control points around that element*/
// 	/*very challenging if just access local vector*/
// 	Vec P_seq;
// 	PetscScalar *Parray;
// 	VecScatter scatter_ctx1;
// 	ierr = VecScatterCreateToAll(x, &scatter_ctx1, &P_seq); CHKERRQ(ierr);
// 	ierr = VecScatterBegin(scatter_ctx1, x, P_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
// 	ierr = VecScatterEnd(scatter_ctx1, x, P_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
// 	ierr = VecGetArray(P_seq, &Parray); CHKERRQ(ierr);

// 	ierr = VecSet(F, 0.0); CHKERRQ(ierr);

// 	NeuronGrowth *user = (NeuronGrowth *)ctx;

// 	/*Build linear system in each process*/
// 	int ind(0); // pre-calculated variable index
// 	for (int e = 0; e < user->bzmesh_process.size(); e++) {
// 		int nen = user->bzmesh_process[e].IEN.size(); // 16 supporting cp for 2D case

// 		vector<vector<float>> EMatrixSolve;	EMatrixSolve.clear();
// 		vector<float> EVectorSolve;		EVectorSolve.clear();
// 		EMatrixSolve.resize(nen);
// 		EVectorSolve.resize(nen);

// 		for (int i = 0; i < nen * 1; i++) {
// 			// EMatrixSolve[i].reserve(nen*sizeof(float));
// 			EMatrixSolve[i].resize(nen);
// 		}
		
// 		// preparing element stiffness matrix and vector, extract values from control mesh
// 		vector<float> elePhiGuess(nen, 0);
// 		for (int i = 0; i < nen; i++) {
// 			EVectorSolve[i] = 0.0;
// 			for (int j = 0; j < nen; j++) {
// 				EMatrixSolve[i][j] = 0.0;
// 			}
// 			elePhiGuess[i] = Parray[user->bzmesh_process[e].IEN[i]];		// elePhiGuess
// 		}

// 		// loop through gaussian quadrature points
// 		for (int i = 0; i < user->Gpt.size(); i++) {
// 			for (int j = 0; j < user->Gpt.size(); j++) {				
// 				user->ElementEvaluationAll_phi(nen, user->pre_Nx[ind], user->pre_dNdx[ind], elePhiGuess, user->vars);
				
// 				// loop through control points
// 				for (int m = 0; m < nen; m++) {
// 					// terma2 = - eleEP * eleEP * (dPGdx * dNdx[m][0] + dPGdy * dNdx[m][1]);
// 					// termadx = - dAdx * eleEEP * dPGdy * dNdx[m][0];
// 					// termady = - dAdy * eleEEP * dPGdx * dNdx[m][1];
// 					// termdbl = (- elePG * elePG * elePG + (1 - C1) * elePG * elePG + C1 * elePG) * Nx[m];
// 					// EVectorSolve[m] += (elePG * Nx[m] - user->dt * user->M_phi * (terma2 - termadx + termady + termdbl) - eleP * Nx[m]) * detJ;
			
// 					EVectorSolve[m] += (user->vars[0] * user->pre_Nx[ind][m] - user->dt * user->pre_eleMp[ind] *\
// 						((- user->pre_eleEP[ind] * user->pre_eleEP[ind] * (user->vars[1] * user->pre_dNdx[ind][m][0] + user->vars[2] * user->pre_dNdx[ind][m][1])) -\
// 						(- user->pre_dAdx[ind] * user->pre_eleEEP[ind] * user->vars[2] * user->pre_dNdx[ind][m][0]) +\
// 						(- user->pre_dAdy[ind] * user->pre_eleEEP[ind] * user->vars[1] * user->pre_dNdx[ind][m][1]) +\
// 						(- user->vars[0] * user->vars[0] * user->vars[0] + (1 - user->pre_C1[ind]) * user->vars[0] * user->vars[0] + user->pre_C1[ind] * user->vars[0]) * user->pre_Nx[ind][m])
// 						- user->pre_eleP[ind] * user->pre_Nx[ind][m]) * user->pre_detJ[ind];

// 					// loop through 16 control points
// 					for (int n = 0; n < nen; n++) {
// 						// terma2 = - eleEP * eleEP * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1]);
// 						// termadx = - dAdx * eleEEP * dNdx[m][1] * dNdx[n][0];
// 						// termady = - dAdy * eleEEP * dNdx[m][0] * dNdx[n][1];
// 						// termdbl = (- 3 * elePG * elePG + 2 * (1 - C1) * elePG + C1 * Nx[m]) * Nx[n];
// 						// EMatrixSolve[m][n] += (Nx[m] * Nx[n] - user->dt * user->M_phi * (terma2 - termadx + termady + termdbl)) * detJ;
						
// 						EMatrixSolve[m][n] += (user->pre_Nx[ind][m] * user->pre_Nx[ind][n] - user->dt * user->pre_eleMp[ind] *
// 							((- user->pre_eleEP[ind] * user->pre_eleEP[ind] * (user->pre_dNdx[ind][m][0] * user->pre_dNdx[ind][n][0] + user->pre_dNdx[ind][m][1] * user->pre_dNdx[ind][n][1])) // terma2
// 							- (- user->pre_dAdx[ind] * user->pre_eleEEP[ind] * user->pre_dNdx[ind][m][1] * user->pre_dNdx[ind][n][0]) // termadx
// 							+ (- user->pre_dAdy[ind] * user->pre_eleEEP[ind] * user->pre_dNdx[ind][m][0] * user->pre_dNdx[ind][n][1]) // termady
// 							+ (- 3 * user->vars[0] * user->vars[0] + 2 * (1 - user->pre_C1[ind]) * user->vars[0] + user->pre_C1[ind] * user->pre_Nx[ind][m]) * user->pre_Nx[ind][n]) // termdbl
// 							) * user->pre_detJ[ind];
// 					}
// 				}
// 				ind += 1; // incrementing index for extracting pre-calculated variables
// 			}
// 		}

// 		/*Apply Boundary Condition*/
// 		for (int i = 0; i < nen; i++) {
// 			int A = user->bzmesh_process[e].IEN[i];
// 			if (user->cpts[A].label == 1) { // domain boundary
// 				user->ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
// 			}
// 		}

// 		user->ResidualAssembly(EVectorSolve, user->bzmesh_process[e].IEN, F);
// 		user->MatrixAssembly(EMatrixSolve, user->bzmesh_process[e].IEN, user->J);
// 	}

// 	ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
// 	ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

// 	ierr = MatAssemblyBegin(user->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
// 	ierr = MatAssemblyEnd(user->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

// 	// if (user->n % 50 == 0) {
// 	// 	PetscInt rows, cols;
// 	// 	ierr = MatGetSize(user->J, &rows, &cols);
// 	// 	PetscPrintf(PETSC_COMM_WORLD, "Global matrix size is %D x %D\n", rows, cols);
// 	// }

// 	if (user->comRank == 0) {
// 		user->check_itr += 1;
// 	}
// 	ierr = VecRestoreArray(P_seq, &Parray); CHKERRQ(ierr);
// 	ierr = VecScatterDestroy(&scatter_ctx1); CHKERRQ(ierr);
// 	ierr = VecDestroy(&P_seq); CHKERRQ(ierr);

// 	return 0;
// }

// PetscErrorCode FormJacobian_phi(SNES snes, Vec x, Mat J, Mat P, void *ctx) 
// {
// 	NeuronGrowth *user = (NeuronGrowth *)ctx;
// 	PetscErrorCode ierr;
// 	ierr = MatCopy(user->J, J, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
// 	return 0;
// }

// PetscErrorCode MySNESMonitor(SNES snes, PetscInt its, PetscReal fnorm, PetscViewerAndFormat *vf)
// {
// 	PetscFunctionBeginUser;
// 	SNESMonitorDefaultShort(snes, its, fnorm, vf);
	
// 	// Set up the KSP solver to monitor iterations
// 	KSP ksp;
// 	SNESGetKSP(snes, &ksp);

// 	// PetscOptionsSetValue(NULL, "-ksp_monitor", "");  // Enables KSP iteration monitoring
// 	PetscOptionsSetValue(NULL, "-ksp_monitor_singular_value", "");  // Enables KSP iteration monitoring
// 	// PetscOptionsSetValue(NULL, "-ksp_monitor_solution", "");  // Enables monitoring of the solution

// 	// Get the number of iterations taken by the KSP solver
// 	PetscInt numIterations;
// 	// KSPGetIterationNumber(ksp, &numIterations);
// 	KSPGetTotalIterations(ksp, &numIterations);
// 	PetscPrintf(PETSC_COMM_WORLD, "     - KSP Iterations: %d\n", numIterations);

// 	PetscFunctionReturn(PETSC_SUCCESS);
// }

// PetscErrorCode CleanUpSolvers(NeuronGrowth &NG)
// {
// 	// NG.ierr = SNESDestroy(&NG.snes_phi); CHKERRQ(NG.ierr);
// 	// NG.ierr = MatDestroy(&NG.J); CHKERRQ(NG.ierr);
// 	// NG.ierr = VecDestroy(&NG.temp_phi); CHKERRQ(NG.ierr);

// 	NG.ierr = KSPDestroy(&NG.ksp_phi); CHKERRQ(NG.ierr);
// 	NG.ierr = MatDestroy(&NG.GK_phi); CHKERRQ(NG.ierr);
// 	NG.ierr = VecDestroy(&NG.GR_phi); CHKERRQ(NG.ierr);
// 	NG.ierr = VecDestroy(&NG.temp_phi); CHKERRQ(NG.ierr);
	
// 	NG.ierr = KSPDestroy(&NG.ksp_syn); CHKERRQ(NG.ierr);
// 	NG.ierr = MatDestroy(&NG.GK_syn); CHKERRQ(NG.ierr);
// 	NG.ierr = VecDestroy(&NG.GR_syn); CHKERRQ(NG.ierr);
// 	NG.ierr = VecDestroy(&NG.temp_syn); CHKERRQ(NG.ierr);

// 	NG.ierr = KSPDestroy(&NG.ksp_tub); CHKERRQ(NG.ierr);	
// 	NG.ierr = MatDestroy(&NG.GK_tub); CHKERRQ(NG.ierr);
// 	NG.ierr = VecDestroy(&NG.GR_tub); CHKERRQ(NG.ierr);
// 	NG.ierr = VecDestroy(&NG.temp_tub); CHKERRQ(NG.ierr);

// 	return NG.ierr;
// }

// int RunNG(int& n_bzmesh, vector<vector<int>> ele_process_in, vector<Vertex2D>& cpts_initial, const vector<Element2D>& tmesh_initial, vector<Vertex2D>& cpts_fine, const vector<Element2D>& tmesh_fine, vector<Vertex2D>& cpts, vector<Vertex2D>& prev_cpts, const vector<Element2D>& tmesh, string path_in, string path_out, int& iter, int end_iter_in,
// 	vector<vector<float>> &NGvars, int &NX, int &NY, vector<array<float, 2>> &seed, int &originX, int &originY, vector<int> &rfid, vector<int> &rftype, bool &localRefine)
// {	
// 	// Declare and initialize Neuron Growth simulation
// 	NeuronGrowth NG;

// 	// Set the number of iterations, number of neurons, and the end iteration from provided variables
// 	NG.n = iter;  // Assign iteration count
// 	NG.numNeuron = seed.size();  // Set the number of neurons based on the size of 'seed'
// 	NG.end_iter = end_iter_in;  // Set the ending iteration

// 	// Initialize vertex clouds for the current, fine, and previous configurations
// 	Vertex2DCloud cloud(cpts);        // Cloud for current points
// 	Vertex2DCloud cloud_fine(cpts_fine);  // Cloud for finer resolution points
// 	Vertex2DCloud cloud_prev(prev_cpts);  // Cloud for previous points

// 	// Initialize KD-Trees for the current, fine, and previous vertex clouds
// 	KDTree kdTree(2 /* dim */, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
// 	KDTree kdTree_fine(2 /* dim */, cloud_fine, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
// 	KDTree kdTree_prev(2 /* dim */, cloud_prev, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));

// 	// Build indexes for KD-Trees to optimize search operations
// 	kdTree.buildIndex();
// 	kdTree_fine.buildIndex();
// 	kdTree_prev.buildIndex();

// 	// Initialize the Neuron Growth problem with the provided parameters and KDTree for previous points
// 	NG.InitializeProblemNG(n_bzmesh, cpts, prev_cpts, kdTree_prev, NGvars, seed);

// 	// Assign processing elements for the simulation
// 	NG.AssignProcessor(ele_process_in);

// 	// Synchronize and check MPI element assignments, print out information in order
// 	for (int i = 0; i < NG.nProcess; i++) {
// 	if (i == NG.comRank) {
// 		std::cout << "comRank: " << NG.comRank << "/" << NG.nProcess << " - with element process size: " << NG.ele_process.size() << std::endl;
// 		NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
// 	} else {
// 		NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
// 	}
// 	}

// 	// Read Bezier mesh for the simulation setup
// 	NG.ReadBezierElementProcess(path_in);
// 	PetscPrintf(PETSC_COMM_WORLD, "Read bzmesh!-----------------------------------------------------------------\n");

// 	// (Optional) Set initial guess for SNES based on the current configuration
// 	// NG.ToPETScVec(NG.phi, NG.temp_phi); // Commented out, indicating it's optional

// 	PetscPrintf(PETSC_COMM_WORLD, "Set initial guess!-----------------------------------------------------------\n");

// 	// /*========================================================*/
// 	// // Initializations
// 	// NeuronGrowth NG;
// 	// // NG.SetVariables("simulation_parameters.txt");
// 	// NG.n = iter;
// 	// NG.numNeuron = seed.size();
// 	// NG.end_iter = end_iter_in;

// 	// Vertex2DCloud cloud(cpts), cloud_fine(cpts_fine), cloud_prev(prev_cpts);
// 	// KDTree kdTree(2 /* dim */, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
// 	// KDTree kdTree_fine(2 /* dim */, cloud_fine, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
// 	// KDTree kdTree_prev(2 /* dim */, cloud_prev, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
// 	// kdTree.buildIndex();
// 	// kdTree_fine.buildIndex();
// 	// kdTree_prev.buildIndex();

// 	// NG.InitializeProblemNG(n_bzmesh, cpts, prev_cpts, kdTree_prev, NGvars, seed);
// 	// NG.AssignProcessor(ele_process_in);
// 	// // Check MPI element assignments, and print out in orders
// 	// for (int i = 0; i < NG.nProcess; i++){
// 	// 	if (i == NG.comRank) {
// 	// 		std::cout << "comRank: " << NG.comRank << "/" << NG.nProcess << " - with element process size: " << NG.ele_process.size() << std::endl;
// 	// 		NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
// 	// 	} else {
// 	// 		NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
// 	// 	}
// 	// }				
// 	// // Read bezier mesh and prepare SNES initial guess
// 	// NG.ReadBezierElementProcess(path_in);
// 	// PetscPrintf(PETSC_COMM_WORLD, "Read bzmesh!-----------------------------------------------------------------\n");	
// 	// // NG.ToPETScVec(NG.phi, NG.temp_phi); // initial guess for SNES (optional)
// 	// PetscPrintf(PETSC_COMM_WORLD, "Set initial guess!-----------------------------------------------------------\n");	
	
// 	/*========================================================*/
// 	// Initial neuron identifications and geodesic distance calculation
// 	vector<vector<int>> neurons, distances;
// 	vector<float> phi_fine, id, geodist, tip, localMaximaMatrix, Mphi;
// 	// phi_fine = InterpolateVars_coarse(NG.phi, cpts, cpts_fine, 1);
// 	// phi_fine = NG.InterpolateVars_coarse1(NG.phi, cpts, kdTree, cpts_fine, 0, 0);
// 	// phi_fine = NG.InterpolateValues_closest(NG.phi, cpts, cpts_fine);
// 	phi_fine = NG.InterpolateValues_closest(NG.phi, kdTree, cpts_fine);
// 	NG.IdentifyNeurons(phi_fine, neurons, seed, NX*2, NY*2, originX, originY);
// 	distances = NG.CalculateGeodesicDistanceFromPoint(neurons, seed, originX, originY);
// 	geodist = ConvertTo1DFloatVector(distances);
// 	PetscPrintf(PETSC_COMM_WORLD, "Calculated geodesic distances!-----------------------------------------------\n");
// 	id = ConvertTo1DFloatVector(neurons);
// 	NG.DetectTipsMulti(phi_fine, id, NG.numNeuron, tip, NX*2, NY*2);
// 	// NG.DetectTipsMulti_new(phi_fine, id, NG.numNeuron, tip, NX*2, NY*2);
// 	localMaximaMatrix = NG.FindCentroidsOfLocalMaximaClusters(tip, NX*2+1, NY*2+1);
// 	NG.tips.clear();
// 	NG.tips.resize(NG.phi.size(),0);
// 	// NG.tips = InterpolateVars_coarse(localMaximaMatrix, cpts_initial, cpts, 0);
// 	// NG.tips = NG.InterpolateVars_coarse1(localMaximaMatrix, cpts_fine, kdTree_fine, cpts, 2, 0);
// 	// NG.tips = NG.InterpolateValues_closest(localMaximaMatrix, cpts_fine, cpts);	
// 	NG.tips = NG.InterpolateValues_closest(localMaximaMatrix, kdTree_fine, cpts);
// 	PetscPrintf(PETSC_COMM_WORLD, "Detected initial tips!-------------------------------------------------------\n");

// 	/*========================================================*/
// 	// Write initial variables
// 	string varName;	
// 	if (NG.n == 0) {
// 		NG.VisualizeVTK_PhysicalDomain_All(0, path_out);
// 		PetscPrintf(PETSC_COMM_WORLD, "Saving all variables!--------------------------------------------------------\n");	
// 	}

// 	/*========================================================*/
// 	// calculate some variable in advance to save computational cost
// 	NG.prepareBasis();
// 	PetscPrintf(PETSC_COMM_WORLD, "Prepared basis!--------------------------------------------------------------\n");	
// 	NG.ierr = MPI_Allreduce(&NG.sum_grad_phi0_local, &NG.sum_grad_phi0_global, 1, MPI_FLOAT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
// 	PetscPrintf(PETSC_COMM_WORLD, "Calculated sum grad phi0!----------------------------------------------------\n");
// 	NG.prepareTerm_source();
// 	// NG.prepareEE();
// 	PetscPrintf(PETSC_COMM_WORLD, "Prepared variables!----------------------------------------------------------\n");
	
// 	while (iter <= NG.end_iter) {
// 		NG.n = iter;
// 		tic();		

// 		/*========================================================*/
// 		/*Implcit Non-liear Newton Raphson solver for Phase field equation*/
// 		NG.phi_prev = NG.phi;	
	
// 		NG.preparePhaseField();

// 		float tol = 1e-4;
// 		float residual = tol + 1;
// 		int NR_itr = 0;
// 		while (residual > tol) {
// 			// PetscPrintf(PETSC_COMM_WORLD, "residual: %f | %f\n", residual, tol);

// 			NG.ierr = MatZeroEntries(NG.GK_phi); CHKERRQ(NG.ierr);
// 			NG.ierr = VecSet(NG.GR_phi, 0); CHKERRQ(NG.ierr);
// 			NG.BuildLinearSystemProcessNG_phi();
// 			NG.ierr = VecAssemblyEnd(NG.GR_phi); CHKERRQ(NG.ierr);
// 			NG.ierr = MatAssemblyEnd(NG.GK_phi, MAT_FINAL_ASSEMBLY); CHKERRQ(NG.ierr);

// 			/*Petsc solver setting for phi*/
// 			if (NG.judge_phi == 0) {
// 				NG.ierr = KSPCreate(PETSC_COMM_WORLD, &NG.ksp_phi); CHKERRQ(NG.ierr);
// 				NG.ierr = KSPSetOperators(NG.ksp_phi, NG.GK_phi, NG.GK_phi); CHKERRQ(NG.ierr);
// 				NG.ierr = KSPGetPC(NG.ksp_phi, &NG.pc_phi); CHKERRQ(NG.ierr);
// 				NG.ierr = PCSetType(NG.pc_phi, PCBJACOBI); CHKERRQ(NG.ierr);
// 				NG.ierr = KSPSetType(NG.ksp_phi, KSPGMRES); CHKERRQ(NG.ierr);
// 				NG.ierr = KSPGMRESSetRestart(NG.ksp_phi, 100); CHKERRQ(NG.ierr);
// 				NG.ierr = KSPSetInitialGuessNonzero(NG.ksp_phi, PETSC_TRUE); CHKERRQ(NG.ierr);
// 				NG.ierr = KSPSetTolerances(NG.ksp_phi, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 100000); CHKERRQ(NG.ierr);
// 				NG.ierr = KSPSetFromOptions(NG.ksp_phi); CHKERRQ(NG.ierr);
// 				NG.ierr = KSPSetUp(NG.ksp_phi); CHKERRQ(NG.ierr);
// 				if (NG.n == 0)
// 					NG.ierr = KSPView(NG.ksp_phi, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(NG.ierr);
// 				NG.judge_phi = 1;
// 			}

// 			NG.ierr = KSPSolve(NG.ksp_phi, NG.GR_phi, NG.temp_phi); CHKERRQ(NG.ierr);
// 			PetscInt its_phi;
// 			NG.ierr = KSPGetIterationNumber(NG.ksp_phi, &its_phi); CHKERRQ(NG.ierr);
// 			KSPConvergedReason reason_phi;
// 			NG.ierr = KSPGetConvergedReason(NG.ksp_phi, &reason_phi); CHKERRQ(NG.ierr);
// 			if (reason_phi < 0) {
// 				PetscPrintf(PETSC_COMM_WORLD, "KSP phi not converging  | KSPreason:  %d | KSPiter: %d \n", reason_phi, its_phi); CHKERRQ(NG.ierr);
// 				return 3;	
// 			}

// 			NR_itr += 1;

// 			/*========================================================*/
// 			/*Collecting scattered Phi variable from all processors*/
// 			Vec temp_phi_seq;
// 			VecScatter ctx_phi;
// 			PetscScalar *_p;

// 			NG.ierr = VecScatterCreateToAll(NG.temp_phi, &ctx_phi, &temp_phi_seq); CHKERRQ(NG.ierr);
// 			NG.ierr = VecScatterBegin(ctx_phi, NG.temp_phi, temp_phi_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
// 			NG.ierr = VecScatterEnd(ctx_phi, NG.temp_phi, temp_phi_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
// 			NG.ierr = VecGetArray(temp_phi_seq, &_p);

// 			residual = 0;
// 			for (int i = 0; i < NG.phi.size(); i++) {
// 				NG.phi[i] = NG.phi[i] - PetscRealPart(_p[i]);
// 				residual = PetscMax(residual, abs(PetscRealPart(_p[i])));

// 				if (abs(NG.phi[i]) > 3) {
// 					PetscPrintf(PETSC_COMM_WORLD, "Incorrect diverging phi!-----------------------------------------------------\n");
// 					return 3;
// 				}
// 			}

// 			NG.ierr = VecRestoreArray(temp_phi_seq, &_p); CHKERRQ(NG.ierr);
// 			NG.ierr = VecScatterDestroy(&ctx_phi); CHKERRQ(NG.ierr);
// 			NG.ierr = VecDestroy(&temp_phi_seq); CHKERRQ(NG.ierr);

// 		}

// 		toc(t_phi);
// 		tic();
// 		t_write += t_phi;
// 		t_total += t_phi;

// 		// /*========================================================*/
// 		// /*Implcit Non-liear SNES solver for Phase field equation*/
// 		// NG.phi_prev = NG.phi;	
	
// 		// NG.preparePhaseField();

// 		// NG.check_itr = 0;
// 		// if (NG.judge_phi == 0) { 
// 		// 	NG.ierr = SNESCreate(PETSC_COMM_WORLD, &NG.snes_phi); CHKERRQ(NG.ierr);
// 		// 	NG.ierr = SNESSetType(NG.snes_phi, SNESNEWTONLS); CHKERRQ(NG.ierr);
// 		// 	NG.ierr = SNESSetTolerances(NG.snes_phi, 1e-4, 1e-6, 1e-8, 100, 1000); CHKERRQ(NG.ierr);
// 		// 	PetscPrintf(PETSC_COMM_WORLD, "Set Tolerance!---------------------------------------------------------------\n");

// 		// 	SNESLineSearch linesearch;
// 		// 	NG.ierr = SNESGetLineSearch(NG.snes_phi, &linesearch); CHKERRQ(NG.ierr);
// 		// 	NG.ierr = SNESLineSearchSetType(linesearch, SNESLINESEARCHCP); CHKERRQ(NG.ierr);
// 		// 	PetscPrintf(PETSC_COMM_WORLD, "Set LineSearch!--------------------------------------------------------------\n");	
// 		// 	// NG.ierr = SNESSetFunction(NG.snes_phi, NULL, FormFunction_phi, &NG); CHKERRQ(NG.ierr);
// 		// 	NG.ierr = SNESSetFunction(NG.snes_phi, NULL, FormFunction_phi_test, &NG); CHKERRQ(NG.ierr);
// 		// 	NG.ierr = SNESSetJacobian(NG.snes_phi, NULL, NULL, FormJacobian_phi, &NG); CHKERRQ(NG.ierr);
// 		// 	PetscPrintf(PETSC_COMM_WORLD, "Set FormFunction and FormJacobian!-------------------------------------------\n");

// 		// 	// PetscViewerAndFormat *vf;
// 		// 	// PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf);
// 		// 	// SNESMonitorSet(NG.snes_phi, (PetscErrorCode(*)(SNES, PetscInt, PetscReal, void *))MySNESMonitor, vf, (PetscErrorCode(*)(void **))PetscViewerAndFormatDestroy);
// 		// 	// PetscPrintf(PETSC_COMM_WORLD, "Set Monitor!-----------------------------------------------------------------\n");
			
// 		// 	if (NG.n == 0)
// 		// 		NG.ierr = SNESView(NG.snes_phi, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(NG.ierr);
// 		// 	NG.judge_phi = 1;
// 		// }


// 		// NG.ierr = SNESSolve(NG.snes_phi, NULL, NG.temp_phi); CHKERRQ(NG.ierr);
// 		// PetscInt its_phi;
// 		// NG.ierr = SNESGetIterationNumber(NG.snes_phi, &its_phi); CHKERRQ(NG.ierr);
// 		// SNESConvergedReason reason_phi;
// 		// NG.ierr = SNESGetConvergedReason(NG.snes_phi, &reason_phi); CHKERRQ(NG.ierr);
// 		// if (reason_phi < 0) {
// 		// 	PetscPrintf(PETSC_COMM_WORLD, "SNES phi not converging ........ | SNESreason: %d | SNESiter: %d \n", reason_phi, its_phi); CHKERRQ(NG.ierr);	
// 		// }

// 		// // if (NG.comRank == 0) {
// 		// // 	std::cout << "Checking snes calls: " << NG.check_itr << std::endl;
// 		// // }

// 		// /*========================================================*/
// 		// /*Collecting scattered Phi variable from all processors*/
// 		// Vec temp_phi_seq;
// 		// VecScatter ctx_phi;
// 		// PetscScalar *_p;

// 		// NG.ierr = VecScatterCreateToAll(NG.temp_phi, &ctx_phi, &temp_phi_seq); CHKERRQ(NG.ierr);
// 		// NG.ierr = VecScatterBegin(ctx_phi, NG.temp_phi, temp_phi_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
// 		// NG.ierr = VecScatterEnd(ctx_phi, NG.temp_phi, temp_phi_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
// 		// NG.ierr = VecGetArray(temp_phi_seq, &_p);

// 		// for (int i = 0; i < NG.phi.size(); i++) {
// 		// 	NG.phi[i] = PetscRealPart(_p[i]);
// 		// }

// 		// NG.ierr = VecRestoreArray(temp_phi_seq, &_p); CHKERRQ(NG.ierr);
// 		// NG.ierr = VecScatterDestroy(&ctx_phi); CHKERRQ(NG.ierr);
// 		// NG.ierr = VecDestroy(&temp_phi_seq); CHKERRQ(NG.ierr);

// 		// toc(t_phi);
// 		// tic();
// 		// t_write += t_phi;
// 		// t_total += t_phi;
		
// 		/*========================================================*/
// 		/*Synaptogenesis and Tubulin equation*/
// 		/*Build Linear System*/
// 		// NG.prepareBasis();
// 		NG.prepareSourceSum();
// 		NG.sum_grad_phi0_global = 0;
// 		NG.ierr = MPI_Allreduce(&NG.sum_grad_phi0_local, &NG.sum_grad_phi0_global, 1, MPI_FLOAT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
// 		NG.prepareTerm_source();

// 		if (NG.judge_syn == 0) {
// 			NG.ierr = MatZeroEntries(NG.GK_syn); CHKERRQ(NG.ierr);
// 		}
// 		NG.ierr = VecSet(NG.GR_syn, 0); CHKERRQ(NG.ierr);
// 		NG.ierr = MatZeroEntries(NG.GK_tub); CHKERRQ(NG.ierr);
// 		NG.ierr = VecSet(NG.GR_tub, 0); CHKERRQ(NG.ierr);
// 		NG.BuildLinearSystemProcessNG_syn_tub(tmesh, cpts); // build both Synaptogenesis and tubulin in the same loop
// 		NG.ierr = VecAssemblyEnd(NG.GR_syn); CHKERRQ(NG.ierr);
// 		if (NG.judge_syn == 0) {
// 			NG.ierr = MatAssemblyEnd(NG.GK_syn, MAT_FINAL_ASSEMBLY); CHKERRQ(NG.ierr);
// 		}
// 		NG.ierr = VecAssemblyEnd(NG.GR_tub); CHKERRQ(NG.ierr);
// 		NG.ierr = MatAssemblyEnd(NG.GK_tub, MAT_FINAL_ASSEMBLY); CHKERRQ(NG.ierr);

// 		/*Petsc solver setting for Synaptogenesis*/
// 		if (NG.judge_syn == 0) {
// 			NG.ierr = KSPCreate(PETSC_COMM_WORLD, &NG.ksp_syn); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPSetOperators(NG.ksp_syn, NG.GK_syn, NG.GK_syn); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPGetPC(NG.ksp_syn, &NG.pc_syn); CHKERRQ(NG.ierr);
// 			NG.ierr = PCSetType(NG.pc_syn, PCBJACOBI); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPSetType(NG.ksp_syn, KSPGMRES); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPGMRESSetRestart(NG.ksp_syn, 100); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPSetInitialGuessNonzero(NG.ksp_syn, PETSC_TRUE); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPSetTolerances(NG.ksp_syn, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 100000); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPSetFromOptions(NG.ksp_syn); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPSetUp(NG.ksp_syn); CHKERRQ(NG.ierr);
// 			if (NG.n == 0)
// 				NG.ierr = KSPView(NG.ksp_syn, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(NG.ierr);
// 			NG.judge_syn = 1;
// 		}

// 		NG.ierr = KSPSolve(NG.ksp_syn, NG.GR_syn, NG.temp_syn); CHKERRQ(NG.ierr);
// 		PetscInt its_syn;
// 		NG.ierr = KSPGetIterationNumber(NG.ksp_syn, &its_syn); CHKERRQ(NG.ierr);
// 		KSPConvergedReason reason_syn;
// 		NG.ierr = KSPGetConvergedReason(NG.ksp_syn, &reason_syn); CHKERRQ(NG.ierr);
// 		if (reason_syn < 0) {
// 			PetscPrintf(PETSC_COMM_WORLD, "KSP Synaptogenesis not converging  | KSPreason:  %d | KSPiter: %d \n", reason_syn, its_syn); CHKERRQ(NG.ierr);	
// 			return 3;
// 		}

// 		/*Petsc solver setting for Tubulin*/
// 		if (NG.judge_tub == 0) {
// 			NG.ierr = KSPCreate(PETSC_COMM_WORLD, &NG.ksp_tub); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPSetOperators(NG.ksp_tub, NG.GK_tub, NG.GK_tub); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPGetPC(NG.ksp_tub, &NG.pc_tub); CHKERRQ(NG.ierr);
// 			NG.ierr = PCSetType(NG.pc_tub, PCBJACOBI); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPSetType(NG.ksp_tub, KSPGMRES); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPGMRESSetRestart(NG.ksp_tub, 100); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPSetInitialGuessNonzero(NG.ksp_tub, PETSC_TRUE); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPSetTolerances(NG.ksp_tub, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 100000); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPSetFromOptions(NG.ksp_tub); CHKERRQ(NG.ierr);
// 			NG.ierr = KSPSetUp(NG.ksp_tub); CHKERRQ(NG.ierr);
// 			if (NG.n == 0)
// 				NG.ierr = KSPView(NG.ksp_tub, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(NG.ierr);
// 			NG.judge_tub = 1;
// 		}

// 		NG.ierr = KSPSolve(NG.ksp_tub, NG.GR_tub, NG.temp_tub); CHKERRQ(NG.ierr);
// 		PetscInt its_tub;
// 		NG.ierr = KSPGetIterationNumber(NG.ksp_tub, &its_tub); CHKERRQ(NG.ierr);
// 		KSPConvergedReason reason_tub;
// 		NG.ierr = KSPGetConvergedReason(NG.ksp_tub, &reason_tub); CHKERRQ(NG.ierr);
// 		if (reason_tub < 0) {
// 			PetscPrintf(PETSC_COMM_WORLD, "KSP tubulin not converging ..... | KSPreason:  %d | KSPiter: %d \n", reason_tub, its_tub); CHKERRQ(NG.ierr);	
// 			return 3;	
// 		}

// 		/*========================================================*/
// 		/*Collecting scattered Synaptogenesis and Tubulin variables from all processors*/
// 		Vec temp_syn_seq, temp_tub_seq;
// 		VecScatter ctx_syn, ctx_tub;
// 		PetscScalar *_s, *_t;

// 		NG.ierr = VecScatterCreateToAll(NG.temp_syn, &ctx_syn, &temp_syn_seq); CHKERRQ(NG.ierr);
// 		NG.ierr = VecScatterBegin(ctx_syn, NG.temp_syn, temp_syn_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
// 		NG.ierr = VecScatterEnd(ctx_syn, NG.temp_syn, temp_syn_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
// 		NG.ierr = VecGetArray(temp_syn_seq, &_s); CHKERRQ(NG.ierr);

// 		NG.ierr = VecScatterCreateToAll(NG.temp_tub, &ctx_tub, &temp_tub_seq); CHKERRQ(NG.ierr);
// 		NG.ierr = VecScatterBegin(ctx_tub, NG.temp_tub, temp_tub_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
// 		NG.ierr = VecScatterEnd(ctx_tub, NG.temp_tub, temp_tub_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
// 		NG.ierr = VecGetArray(temp_tub_seq, &_t); CHKERRQ(NG.ierr);

// 		for (int i = 0; i < NG.syn.size(); i++) {
// 			NG.syn[i] = PetscRealPart(_s[i]);
// 			NG.tub[i] = PetscRealPart(_t[i]) * round(NG.phi[i]);
// 		}

// 		NG.ierr = VecRestoreArray(temp_syn_seq, &_s); CHKERRQ(NG.ierr);
// 		NG.ierr = VecScatterDestroy(&ctx_syn); CHKERRQ(NG.ierr);
// 		NG.ierr = VecDestroy(&temp_syn_seq); CHKERRQ(NG.ierr);

// 		NG.ierr = VecRestoreArray(temp_tub_seq, &_t); CHKERRQ(NG.ierr);
// 		NG.ierr = VecScatterDestroy(&ctx_tub); CHKERRQ(NG.ierr);
// 		NG.ierr = VecDestroy(&temp_tub_seq); CHKERRQ(NG.ierr);

// 		toc(t_synTub);
// 		tic();
// 		t_write += t_synTub;
// 		t_total += t_synTub;

// 		/*========================================================*/
// 		// Obtain initial local refinement information, the very first 25 iterations are purely used 
// 		// for getting diffused interface for applying local refinements (phi initialization is binary) 
// 		if ((NG.n == 25) && (localRefine == false)) {
// 			NGvars.clear(); NGvars.resize(6);
// 			NGvars[0] = NG.phi;	
// 			NGvars[1] = NG.syn;
// 			NGvars[2] = NG.tub;
// 			NGvars[3] = NG.theta;
// 			NGvars[4] = NG.phi_0;
// 			NGvars[5] = NG.tub_0;

// 			CleanUpSolvers(NG); // Destroy solvers
// 			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
// 			iter = 0; // reset to 0 (beginning of the simulation)
// 			localRefine = true;
// 			return 2;
// 		}
		
// 		/*========================================================*/
// 		// Neuron identification and tip detection
// 		// if ((NG.n % 5 == 0) && (NG.n != 0)) {
		
// 		float t_cc(0);
// 		if ((NG.n % 1 == 0) || (NG.n == 0)) {	

// 			PetscPrintf(PETSC_COMM_WORLD, "Step:%d/%d | Phi:[%d]%.3fs | Syn:%d[%d] Tub:%d[%d] %.3fs | Tip:%.3fs | DOFs:%d |\n",\
// 			NG.n, NG.end_iter, NR_itr, t_phi, reason_syn, its_syn, reason_tub, its_tub, t_synTub, t_tip, NG.phi.size()); CHKERRQ(NG.ierr);

// 			// phi_fine = InterpolateVars_coarse(NG.phi, cpts, cpts_fine, 0);
// 			// phi_fine = NG.InterpolateVars_coarse1(NG.phi, cpts, kdTree, cpts_fine, 0, 0);
// 			// phi_fine = NG.InterpolateValues_closest(NG.phi, cpts, cpts_fine);
// 			phi_fine = NG.InterpolateValues_closest(NG.phi, kdTree, cpts_fine);
// 			// PetscPrintf(PETSC_COMM_WORLD, "Got fine |");
// 			NG.IdentifyNeurons(phi_fine, neurons, seed, NX*2, NY*2, originX, originY);
// 			// PetscPrintf(PETSC_COMM_WORLD, "Got id |");
// 			distances = NG.CalculateGeodesicDistanceFromPoint(neurons, seed, originX, originY);
// 			geodist = ConvertTo1DFloatVector(distances);
// 			// NG.PrintOutNeurons(neurons);
// 			id = ConvertTo1DFloatVector(neurons);
// 			// PetscPrintf(PETSC_COMM_WORLD, "Got geo |");
// 			NG.DetectTipsMulti(phi_fine, id, NG.numNeuron, tip, NX*2, NY*2);
// 			// PetscPrintf(PETSC_COMM_WORLD, "Got tip |");
// 			// NG.DetectTipsMulti_new(phi_fine, id, NG.numNeuron, tip, NX*2, NY*2);
// 			// localMaximaMatrix = NG.FindLocalMaximaInClusters(tip, NX*2+1, NY*2+1);
// 			// localMaximaMatrix = NG.FindCentroidsInClusters(tip, NX+1, NY+1);
//  			localMaximaMatrix = NG.FindCentroidsOfLocalMaximaClusters(tip, NX*2+1, NY*2+1);
// 			// PetscPrintf(PETSC_COMM_WORLD, "Got loMax 0|");

// 			int maxGeoInd(0);
// 			float maxGeo(0);
// 			for (int i = 0; i < geodist.size(); i++) {
// 				if ((geodist[i] != INF) && (geodist[i] > maxGeo)) {
// 					maxGeo = geodist[i];
// 					maxGeoInd = i;
// 				}
// 			}
// 			if (maxGeoInd > 1) {
// 				for (int i = -1; i < 1; i++) { // increase Mphi at longest tip 
// 					for (int j = -1; j < 1; j++) {
// 							localMaximaMatrix[maxGeoInd+j*(2*NY+1)-i] = 5;
// 					}
// 				}
// 			}

// 			// localMaximaMatrix[maxGeoInd] = 5;
// 			// PetscPrintf(PETSC_COMM_WORLD, "Got loMax 1|");

// 			// NG.tips = InterpolateVars_coarse(localMaximaMatrix, cpts_fine, cpts, 0);
// 			// NG.tips = NG.InterpolateVars_coarse1(localMaximaMatrix, cpts_fine, kdTree_fine, cpts, 2, 0);
// 			// NG.tips = NG.InterpolateValues_closest(localMaximaMatrix, cpts_fine, cpts);
// 			NG.tips = NG.InterpolateValues_closest(localMaximaMatrix, kdTree_fine, cpts);
// 			// PetscPrintf(PETSC_COMM_WORLD, "Got NGtips|\n");

// 			toc(t_tip);
// 			tic();
// 			t_write += t_tip;
// 			t_total += t_tip;
// 		}

// 		if ((NG.n % NG.var_save_invl == 0) && (NG.n != 0)) {	
// 			NG.PrintOutNeurons(neurons);
// 			/*========================================================*/
// 			// Writing results to files
// 			PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------\n");
// 			NG.VisualizeVTK_PhysicalDomain_All(NG.n, path_out);
			
// 			// if (NG.comRank == 0) {
// 			// 	string outputPath = path_out + "/phi_" + to_string(NG.n) + ".png";
// 			// 	WriteMatrixToPNG(NG.phi, NX+1, NY+1, outputPath.c_str());
// 			// }
// 			varName = "phifine_running_";
// 			NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, NG.n, path_out, phi_fine, varName); // solution on control points
			
// 			varName = "geoDist_running_";
// 			NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, NG.n, path_out, geodist, varName); // solution on control points
// 			varName = "localMax_running_";
// 			NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, NG.n, path_out, localMaximaMatrix, varName); // solution on control points
// 			varName = "tip_running_";
// 			NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, NG.n, path_out, tip, varName); // solution on control points
// 			varName = "id_running_";
// 			NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, NG.n, path_out, id, varName); // solution on control points
// 			// varName = "phi_running_";
// 			// NG.VisualizeVTK_ControlMesh(cpts, tmesh, NG.n, path_out, NG.phi, varName); // solution on control points
// 			// varName = "Mphi_running_";
// 			// NG.VisualizeVTK_ControlMesh(cpts_initial, tmesh_initial, NG.n, path_out, Mphi, varName); // solution on control points
			
// 			PetscPrintf(PETSC_COMM_WORLD, "| Wrote Physical Domain! | Average time %fs | Total time: %f |\n", 
// 				t_total/NG.n, t_total); CHKERRQ(NG.ierr);
// 				// t_write/NG.var_save_invl, t_total); CHKERRQ(NG.ierr);
// 			PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------\n");
// 			t_write = 0;
// 		}

// 		/*========================================================*/
// 		// Domain expansion and variable passing - back to main.cpp
// 		if ((NG.n % NG.expandCK_invl == 0) && (NG.n != 0)) {
// 			// NG.PrintOutNeurons(neurons); 
// 			// int expd_dir_local = NG.CheckExpansion(id, NX, NY);
// 			int expd_dir_local = NG.CheckExpansion(NG.phi, NX, NY);
// 			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
// 			int expd_dir_global = 4; // 4 - No expansion
// 			NG.ierr = MPI_Allreduce(&expd_dir_local, &expd_dir_global, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
// 			int expd_dir = expd_dir_global/NG.comSize; // approximate choice, different thread could have different choice 

// 			if (expd_dir <= 3) {
// 				switch (expd_dir) {
// 					case 0: // left
// 						NX += 3;	originX -= 6;		break;
// 					case 1: // top
// 						NY += 3;				break;
// 					case 2: // right
// 						NX += 3;				break;
// 					case 3: // bottom
// 						NY += 3;	originY -= 6;		break;
// 				}
// 			}

// 			NGvars.clear(); NGvars.resize(6);
// 			NGvars[0] = NG.phi;	
// 			NGvars[1] = NG.syn;
// 			NGvars[2] = NG.tub;
// 			NGvars[3] = NG.theta;
// 			NGvars[4] = NG.phi_0;
// 			NGvars[5] = NG.tub_0;

// 			CleanUpSolvers(NG); // Destroy solvers
// 			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
// 			iter += 1;
// 			return 2;

// 		}
// 		iter += 1;		 
// 	}

// 	/*========================================================*/
// 	CleanUpSolvers(NG); // Destroy solvers
// 	NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);

// 	return 0; // ending simulation
// }

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

#include <nanoflann.hpp> // for nearest points search with KD-tree
// https://github.com/jlblancoc/nanoflann
// sudo apt install libnanoflann-dev

using namespace std;
typedef unsigned int uint;
const float PI = 3.1415926;

std::stack<clock_t> tictoc_stack;
float t_phi, t_synTub, t_syn, t_tub, t_tip, t_write(0), t_total(0);

const int INF = 1e9;

void tic() 
{
	tictoc_stack.push(clock());
}
void toc(float &t) 
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


// // Dataset adaptor for nanoflann
// struct Vertex2DCloud {
//     const std::vector<Vertex2D>& pts;

//     Vertex2DCloud(const std::vector<Vertex2D>& pts) : pts(pts) {}

//     // Must return the number of data points
//     inline size_t kdtree_get_point_count() const { return pts.size(); }

//     // Returns the dim'th component of the idx'th point in the class
//     inline float kdtree_get_pt(const size_t idx, const size_t dim) const {
//         return dim == 0 ? pts[idx].coor[0] : pts[idx].coor[1];
//     }

//     // Optional: bounding-box computation
//     template <class BBOX>
//     bool kdtree_get_bbox(BBOX&) const { return false; }
// };

// using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
//     nanoflann::L2_Simple_Adaptor<float, Vertex2DCloud>,
//     Vertex2DCloud, 2 /* dim */>;
    

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

	// // integer variable setup
	// var_save_invl		= 25; 		// var_save_invl
	// expandCK_invl		= 100; 		// var_save_invl
	// gc_sz			= 2;	     	// gc_sz
	// aniso 			= 6;   		// aniso
	// gamma 			= 10;  		// gamma
	// seed_radius 		= 10;		// seed radius

	// // variable setup
	// kappa			= 2.0;		// kappa;
	// dt			= 1e-2;		// time step
	// Dc			= 3;		// syn D
	// alpha			= 0.9;		// alpha
	// alphaOverPi		= alpha / PI; 	// alphOverPix
	// M_phi			= 10;		// M_phi
	// s_coeff			= 0.007;	// s_coeff
	// delta			= 0.20;		// delta
	// epsilonb		= 0.03;		// epsilonb
	// r			= 5;		// r
	// g			= 0.05;		// g
	// alphaT 			= 0.005;			// alpha_t
	// betaT			= 0.001;	// beta_t
	// Diff			= 4;		// Diff
	// source_coeff		= 15;		// source_coeff

	// integer variable setup
	expandCK_invl		= 100; 		// var_save_invl
	var_save_invl		= 100; 		// var_save_invl
	numNeuron 		= 1;	     	// numNeuron
	gc_sz			= 2;	     	// gc_sz
	aniso 			= 6;   		// aniso
	gamma 			= 10;  		// gamma
	// seed_radius 		= 15;		// seed radius
	seed_radius 		= 10;		// seed radius

	// variable setup
	kappa			= 2.0;		// kappa;
	dt			= 1e-3;		// time step
	Dc			= 3;		// syn D
	alpha			= 0.9;		// alpha
	alphaOverPi		= alpha / PI; 	// alphOverPix
	M_phi			= 60;		// M_phi
	s_coeff			= 0.007;	// s_coeff
	delta			= 0.20;		// delta
	epsilonb		= 0.04;		// epsilonb
	r			= 5;		// r
	g			= 0.1;		// g
	alphaT 			= 0.001;			// alpha_t
	betaT			= 0.001;	// beta_t
	Diff			= 4;		// Diff
	source_coeff		= 15;		// source_coeff
}

void NeuronGrowth::AssignProcessor(vector<vector<int>> &ele_proc)
{
	ele_process.clear();
	for (int i = 0; i < ele_proc[comRank].size(); i++)
		ele_process.push_back(ele_proc[comRank][i]);
}

void NeuronGrowth::SetVariables(string fn_par)
{
	string fname(fn_par), stmp;
	stringstream ss;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open()) {
		fin >> stmp >> var_save_invl;
		fin >> stmp >> expandCK_invl;
		fin >> stmp >> gc_sz;
		fin >> stmp >> aniso;
		fin >> stmp >> gamma;
		fin >> stmp >> seed_radius;
		fin >> stmp >> kappa;
		fin >> stmp >> dt;
		fin >> stmp >> Dc;
		fin >> stmp >> alpha;
		alphaOverPi = alpha / PI; 	// alphOverPi
		fin >> stmp >> M_phi;
		fin >> stmp >> s_coeff;
		fin >> stmp >> delta;
		fin >> stmp >> epsilonb;
		fin >> stmp >> r;
		fin >> stmp >> g;
		fin >> stmp >> alphaT;
		fin >> stmp >> betaT;
		fin >> stmp >> Diff;
		fin >> stmp >> source_coeff;

		fin.close();
		PetscPrintf(PETSC_COMM_WORLD, "Parameter Loaded!\n");
	} else {
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}

void NeuronGrowth::InitializeProblemNG(const int n_bz, vector<Vertex2D>& cpts, vector<Vertex2D> prev_cpts, const KDTree& kdTree_prev, vector<vector<float>> &NGvars, vector<array<float, 2>> &seed)
{
	MPI_Barrier(PETSC_COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD, "Initializing-----------------------------------------------------------------\n");

	/*Initialize parameters*/
	int cpt_sz = cpts.size();
	GaussInfo(3);
	n_bzmesh = n_bz;
	N_0.resize(cpt_sz);

	// float max_x, min_x, max_y, min_y;
	for (int i = 0; i<cpts.size(); i++) {
		if (n==0) {
			max_x = cpts[i].coor[0];
			min_x = cpts[i].coor[0];
			max_y = cpts[i].coor[1];
			min_y = cpts[i].coor[1];

			if (i < prev_cpts.size()) {
				prev_max_x = prev_cpts[i].coor[0];
				prev_min_x = prev_cpts[i].coor[0];
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

	// // Vertex2DCloud cloud(cpts), cloud_fine(cpts_fine), cloud_prev(prev_cpts);
	// KDTree kdTree(2 /* dim */, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	// KDTree kdTree_fine(2 /* dim */, cloud_fine, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	// KDTree kdTree_prev(2 /* dim */, cloud_prev, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	// kdTree.buildIndex();
	// kdTree_fine.buildIndex();
	// kdTree_prev.buildIndex();

	float r, x, y; // radius, x, y
	int nx = sqrt(cpts.size()); // only used at the beginning (square mesh)
	float dx = max_x/(float)nx; // only used at the beginning (square mesh)
	dx = 1;
	srand(time(NULL)); // random seed

	if (n == 0) { // Initialize
		for (int i = 0; i<cpts.size(); i++) {
			phi.push_back(0.f);
			tub.push_back(0.f);
			theta.push_back((float)(rand()%100)/(float)100); // orientation theta
			syn.push_back(0.f); // synaptogenesis
			Mphi.push_back(M_phi);
			prev_id.push_back(0.f);
		}
		for (int i = 0; i<cpts.size(); i++) {
			x = cpts[i].coor[0];
			y = cpts[i].coor[1];

			// calculate new boundary label
			if (x==min_x || x==max_x || y==min_y || y==max_y) {
				cpts[i].label = 1;
			} else {
				cpts[i].label = 0;
			}   

			for (int j = 0; j<seed.size(); j++) { // multiple neurons
				r = sqrt((x-(float)seed[j][0])*(x-(float)seed[j][0])+(y-(float)seed[j][1])*(y-(float)seed[j][1])); // radius of initial soma
				if (r < (seed_radius*dx)) { // initial soma
					// std::cout << seed_radius << " " << dx << " " << seed_radius*dx << std::endl;
					phi[i] = 1.0f;
					tub[i] = (0.5+0.5*tanh((sqrt(seed_radius*dx)-r)/2)); // based eqn in literature
				} 
			} 
		}
		phi_0 = phi; // initial phi
		tub_0 = tub; // initial tub
	} else { // Collect from NGvars - continue simulation
		// phi.clear(); 	phi.resize(cpts.size());
		// syn.clear(); 	syn.resize(cpts.size());
		// tub.clear(); 	tub.resize(cpts.size());
		// theta.clear();	theta.resize(cpts.size());
		// phi_0.clear(); 	phi_0.resize(cpts.size());
		// tub_0.clear(); 	tub_0.resize(cpts.size());
		// for (int i = 0; i < prev_cpts.size(); i++) {
		// 	if (abs(remainder(prev_cpts[i].coor[0],2)) != 1)
		// 		prev_cpts[i].coor[0] = round(prev_cpts[i].coor[0]);
		// 	if (abs(remainder(prev_cpts[i].coor[1],2)) != 1)
		// 		prev_cpts[i].coor[1] = round(prev_cpts[i].coor[1]);
		// }

		// for (int i = 0; i < cpts.size(); i++) {
		// 	x = cpts[i].coor[0];
		// 	if (abs(remainder(x,2)) != 1) {
		// 		x = round(x);
		// 	}
		// 	y = cpts[i].coor[1];
		// 	if (abs(remainder(y,2)) != 1) {
		// 		y = round(y);
		// 	}

		// 	// calculate new boundary label
		// 	if (x==min_x || x==max_x || y==min_y || y==max_y) {
		// 		cpts[i].label = 1;
		// 	} else {
		// 		cpts[i].label = 0;
		// 	}   

		// 	int ind, indDown, indUp, indLeft, indRight;
		// 	if (SearchPair(prev_cpts, round5(x), round5(y), ind)) {
		// 	// if (SearchPair(prev_cpts, cpts[i].coor[0], cpts[i].coor[1], ind)) {
		// 		phi[i] = NGvars[0][ind];
		// 		syn[i] = NGvars[1][ind];
		// 		tub[i] = NGvars[2][ind];
		// 		theta[i] = NGvars[3][ind];
		// 		phi_0[i] = NGvars[4][ind];
		// 		tub_0[i] = NGvars[5][ind];
		// 	} else if ((abs(remainder(x,2)) == 1) && (abs(remainder(y,2)) != 1) 
		// 		&& SearchPair(prev_cpts, floorf(x)-1, round5(y), indDown)
		// 		&& SearchPair(prev_cpts, floorf(x)+1, round5(y), indUp)) {
		// 			phi[i] = (NGvars[0][indDown] + NGvars[0][indUp])/2;
		// 			syn[i] = (NGvars[1][indDown] + NGvars[1][indUp])/2;
		// 			tub[i] = (NGvars[2][indDown] + NGvars[2][indUp])/2;
		// 			// theta[i] = max(NGvars[3][indDown], NGvars[3][indUp]);
		// 			theta[i] = (float)(rand()%100)/(float)100;
		// 			phi_0[i] = (NGvars[4][indDown] + NGvars[4][indUp])/2;
		// 			tub_0[i] = (NGvars[5][indDown] + NGvars[5][indUp])/2;
		// 	} else if ((abs(remainder(x,2)) != 1) && (abs(remainder(y,2)) == 1)
		// 		&& SearchPair(prev_cpts, round5(x), floorf(y)-1, indDown)
		// 		&& SearchPair(prev_cpts, round5(x), floorf(y)+1, indUp)){
		// 			phi[i] = (NGvars[0][indDown] + NGvars[0][indUp])/2;
		// 			syn[i] = (NGvars[1][indDown] + NGvars[1][indUp])/2;
		// 			tub[i] = (NGvars[2][indDown] + NGvars[2][indUp])/2;
		// 			// theta[i] = max(NGvars[3][indDown], NGvars[3][indUp]);
		// 			theta[i] = (float)(rand()%100)/(float)100;
		// 			phi_0[i] = (NGvars[4][indDown] + NGvars[4][indUp])/2;
		// 			tub_0[i] = (NGvars[5][indDown] + NGvars[5][indUp])/2;
		// 	} else if ((abs(remainder(x,2)) == 1) && (abs(remainder(y,2)) == 1)
		// 		&& SearchPair(prev_cpts, floorf(x)-1, floorf(y)-1, indDown)
		// 		&& SearchPair(prev_cpts, floorf(x)+1, floorf(y)+1, indUp)
		// 		&& SearchPair(prev_cpts, floorf(x)-1, floorf(y)+1, indLeft)
		// 		&& SearchPair(prev_cpts, floorf(x)+1, floorf(y)-1, indRight)) {
		// 			phi[i] = (NGvars[0][indDown] + NGvars[0][indUp] + NGvars[0][indLeft] + NGvars[0][indRight])/4;
		// 			syn[i] = (NGvars[1][indDown] + NGvars[1][indUp] + NGvars[1][indLeft] + NGvars[1][indRight])/4;
		// 			tub[i] = (NGvars[2][indDown] + NGvars[2][indUp] + NGvars[2][indLeft] + NGvars[2][indRight])/4;
		// 			// theta[i] = max(max(NGvars[3][indDown], NGvars[3][indUp]), max(NGvars[3][indLeft], NGvars[3][indRight]));
		// 			theta[i] = (float)(rand()%100)/(float)100;
		// 			phi_0[i] = (NGvars[4][indDown] + NGvars[4][indUp] + NGvars[4][indLeft] + NGvars[4][indRight])/4;
		// 			tub_0[i] = (NGvars[5][indDown] + NGvars[5][indUp] + NGvars[5][indLeft] + NGvars[5][indRight])/4;
		// 	} else {
		// 		phi[i] = 0;
		// 		syn[i] = 0;
		// 		tub[i] = 0;
		// 		theta[i] = (float)(rand()%100)/(float)100;
		// 		phi_0[i] = 0;
		// 		tub_0[i] = 0;
		// 	}		
		// 	Mphi.push_back(M_phi);
		// }

		phi.clear(); 	phi.resize(cpts.size());
		syn.clear(); 	syn.resize(cpts.size());
		tub.clear(); 	tub.resize(cpts.size());
		theta.clear();	theta.resize(cpts.size());
		phi_0.clear(); 	phi_0.resize(cpts.size());
		tub_0.clear(); 	tub_0.resize(cpts.size());
		prev_id.clear(); prev_id.resize(cpts.size());

		// phi = InterpolateVars_coarse(NGvars[0], prev_cpts, cpts, 1);
		// syn = InterpolateVars_coarse(NGvars[1], prev_cpts, cpts, 1);
		// tub = InterpolateVars_coarse(NGvars[2], prev_cpts, cpts, 1);
		// theta = InterpolateVars_coarse(NGvars[3], prev_cpts, cpts, 1);
		// phi_0 = InterpolateVars_coarse(NGvars[4], prev_cpts, cpts, 1);
		// tub_0 = InterpolateVars_coarse(NGvars[5], prev_cpts, cpts, 1);

		phi = InterpolateVars_coarse1(NGvars[0], prev_cpts, kdTree_prev, cpts, 1, 0);
		syn = InterpolateVars_coarse1(NGvars[1], prev_cpts, kdTree_prev, cpts, 1, 0);
		tub = InterpolateVars_coarse1(NGvars[2], prev_cpts, kdTree_prev, cpts, 1, 0);
		theta = InterpolateVars_coarse1(NGvars[3], prev_cpts, kdTree_prev, cpts, 1, 1); // isTheta = 1 for special init
		phi_0 = InterpolateVars_coarse1(NGvars[4], prev_cpts, kdTree_prev, cpts, 1, 0);
		tub_0 = InterpolateVars_coarse1(NGvars[5], prev_cpts, kdTree_prev, cpts, 1, 0);
		prev_id = InterpolateVars_coarse1(NGvars[6], prev_cpts, kdTree_prev, cpts, 1, 0);

		for (int i = 0; i < cpts.size(); i++) {	
			Mphi.push_back(M_phi);
		}

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
	int npts, neles, nfunctions, itmp, itmp1;
	int add(0);

	string fname_cmat = fn + "cmat.txt";
	ifstream fin_cmat;
	fin_cmat.clear();
	fin_cmat.open(fname_cmat);

	if (fin_cmat.is_open())
	{
		fin_cmat >> neles;
		bzmesh_process.clear();
		bzmesh_process.resize(ele_process.size());
		for (int i = 0; i < neles; i++)
		{
			if ((i == ele_process[add]) && (add < ele_process.size()))
			{
				fin_cmat >> itmp >> nfunctions >> bzmesh_process[add].type;
				bzmesh_process[add].cmat.resize(nfunctions);
				bzmesh_process[add].IEN.resize(nfunctions);
				for (int j = 0; j < nfunctions; j++)
					fin_cmat >> bzmesh_process[add].IEN[j];
				for (int j = 0; j < nfunctions; j++)
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
				for (int j = 0; j < nfunctions; j++)
					fin_cmat >> stmp;
				for (int j = 0; j < nfunctions; j++)
					for (int k = 0; k < 16; k++)
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
		for (int e = 0; e < neles; e++)
		{
			if ((e == ele_process[add]) && (add < ele_process.size()))
			{
				bzmesh_process[add].pts.resize(bzpt_num);
				for (int i = 0; i < bzpt_num; i++)
				{
					fin_bzpt >> bzmesh_process[add].pts[i][0] >> bzmesh_process[add].pts[i][1] >> tmp_xuan;
				}
				add++;
			}
			else
			{
				for (int i = 0; i < bzpt_num; i++)
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

void NeuronGrowth::BasisFunction(float u, float v, const vector<array<float, 2>> &pt, const vector<array<float, 16>> &cmat, vector<float> &Nx, vector<array<float, 2>> &dNdx, float dudx[2][2], float &detJ)
{
	float Nu[4] = {(1. - u) * (1. - u) * (1. - u), 3. * (1. - u) * (1. - u) * u, 3. * (1. - u) * u * u, u * u * u};
	float Nv[4] = {(1. - v) * (1. - v) * (1. - v), 3. * (1. - v) * (1. - v) * v, 3. * (1. - v) * v * v, v * v * v};
	// float Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
	float dNdu[4] = {-3. * (1. - u) * (1. - u), 3. - 12. * u + 9. * u * u, 3. * (2. - 3. * u) * u, 3. * u * u};
	float dNdv[4] = {-3. * (1. - v) * (1. - v), 3. - 12. * v + 9. * v * v, 3. * (2. - 3. * v) * v, 3. * v * v};
	// float dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
	float dNdt[bzpt_num][2];
	float Nx_bz[bzpt_num];
	float dNdx_bz[bzpt_num][2];

	Nx.clear();
	dNdx.clear();
	Nx.resize(cmat.size(), 0);
	dNdx.resize(cmat.size(), {0, 0});

	int i, j, a, b, c, loc;
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

void NeuronGrowth::BasisFunction(float u, float v, int nen, const vector<array<float, 2>> &pt, const vector<array<float, 16>> &cmat, vector<float> &Nx, vector<array<float, 2>> &dNdx, vector<array<array<float, 2>, 2>> &dN2dx2, float dudx[2][2], float &detJ)
{
	float Nu[4] = {(1. - u) * (1. - u) * (1. - u), 3. * (1. - u) * (1. - u) * u, 3. * (1. - u) * u * u, u * u * u};
	float Nv[4] = {(1. - v) * (1. - v) * (1. - v), 3. * (1. - v) * (1. - v) * v, 3. * (1. - v) * v * v, v * v * v};
	// float Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };

	float dNdu[4] = {-3. * (1. - u) * (1. - u), 3. - 12. * u + 9. * u * u, 3. * (2. - 3. * u) * u, 3. * u * u};
	float dNdv[4] = {-3. * (1. - v) * (1. - v), 3. - 12. * v + 9. * v * v, 3. * (2. - 3. * v) * v, 3. * v * v};
	// float dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };

	float dN2du2[4] = {6. * (1. - u), -12. + 18. * u, 6. - 18. * u, 6. * u};
	float dN2dv2[4] = {6. * (1. - v), -12. + 18. * v, 6. - 18. * v, 6. * v};
	// float dN2dw2[4] = { 6.*(1. - w),-12. + 18.*w,6. - 18.*w,6.*w };

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

void NeuronGrowth::WeightingFunction(const float velocity[2], const float &s, const float &tau, const vector<float> &Nx, const vector<array<float, 2>> &dNdx, vector<float> &Wx)
{
	/*Calculate weighting function using SUPG parameter Tau*/
	float U = sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] /* + velocity[2] * velocity[2] */);
	for (int loc = 0; loc < Nx.size(); loc++)
		Wx[loc] = Nx[loc] + tau * (velocity[0] * dNdx[loc][0] + velocity[1] * dNdx[loc][1] /* + velocity[2] * dNdx[loc][2] */);
}

void NeuronGrowth::ApplyBoundaryCondition(const float bc_value, int pt_num, int variable_num, vector<vector<float>> &EMatrixSolve, vector<float> &EVectorSolve)
{
	int j, k;
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
	int i, j, A, B, m, n;
	int row_start, row_end, row_now;
	int add = 0;
	int nen = IEN.size();

	PetscInt *nodeList = new PetscInt[nen * 1];
	PetscReal *tmpGK = new PetscReal[nen * 1 * nen * 1];

	for (m = 0; m < IEN.size(); m++)
	{
		A = IEN[m];
		nodeList[1 * m] = 1 * A;
		for (n = 0; n < IEN.size(); n++)
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
	int i, m, A;
	int add = 0;
	int nen = IEN.size();

	PetscInt *nodeList = new PetscInt[nen];
	PetscReal *tmpGR = new PetscReal[nen];

	for (i = 0; i < IEN.size(); i++)
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

void NeuronGrowth::MatrixAssembly_insert(vector<vector<float>> &EMatrixSolve, const vector<int> &IEN, Mat &GK)
{
	int i, j, A, B, m, n;
	int row_start, row_end, row_now;
	int add = 0;
	int nen = IEN.size();

	PetscInt *nodeList = new PetscInt[nen * 1];
	PetscReal *tmpGK = new PetscReal[nen * 1 * nen * 1];

	for (m = 0; m < IEN.size(); m++)
	{
		A = IEN[m];
		nodeList[1 * m] = 1 * A;
		for (n = 0; n < IEN.size(); n++)
		{
				tmpGK[add] = EMatrixSolve[m][n];
				add++;
		}
	}
	MatSetValues(GK, nen * 1, nodeList, nen * 1, nodeList, tmpGK, INSERT_VALUES);
	delete nodeList;
	delete tmpGK;
}

void NeuronGrowth::ResidualAssembly_insert(vector<float> &EVectorSolve, const vector<int> &IEN, Vec &GR)
{
	int i, m, A;
	int add = 0;
	int nen = IEN.size();

	PetscInt *nodeList = new PetscInt[nen];
	PetscReal *tmpGR = new PetscReal[nen];

	for (i = 0; i < IEN.size(); i++)
	{
		A = IEN[i];
		nodeList[i * 1 + 0] = A * 1 + 0;
		tmpGR[add] = EVectorSolve[i];
		add += 1;
	}
	VecSetValues(GR, nen * 1, nodeList, tmpGR, INSERT_VALUES);
	delete nodeList;
	delete tmpGR;
}

void NeuronGrowth::MatrixAssembly_2var(vector<vector<float>>& EMatrixSolve, const vector<int>& IEN, Mat& GK)
{
	int i, j, A, B, m, n;
	int row_start, row_end, row_now;
	int add = 0;
	int nen = IEN.size();


	PetscInt *nodeList = new PetscInt[nen * 2];	
	PetscReal *tmpGK = new PetscReal[nen * 2 * nen * 2];
	
	for (m = 0; m<IEN.size(); m++){
		A = IEN[m];
		for (i = 0; i < 2; i++)	{
			nodeList[2 * m + i] = 2 * A + i;
			for (n = 0; n<IEN.size(); n++){
				for (j = 0; j < 2; j++)	{
					tmpGK[add] = EMatrixSolve[m + i*nen][n + j*nen];
					add++;
				}   
			}
		}
	}	
	MatSetValues(GK, nen * 2, nodeList, nen * 2, nodeList, tmpGK, ADD_VALUES);
	delete nodeList;
	delete tmpGK;
}

void NeuronGrowth::ResidualAssembly_2var(vector<float>& EVectorSolve, const vector<int>& IEN, Vec& GR)
{
	int i, m, A;
	int add = 0;
	int nen = IEN.size();

	PetscInt *nodeList = new PetscInt[nen * 2];
	PetscReal *tmpGR = new PetscReal[nen * 2];
	
	for (i = 0; i<IEN.size(); i++){
		A = IEN[i];
		nodeList[i * 2 + 0] = A * 2 + 0;
		nodeList[i * 2 + 1] = A * 2 + 1;
		tmpGR[add] = EVectorSolve[i];
		tmpGR[add + 1] = EVectorSolve[i + nen];
		add+=2;
	}
	VecSetValues(GR, nen * 2, nodeList, tmpGR, ADD_VALUES);
	delete nodeList;
	delete tmpGR;
}

void NeuronGrowth::VisualizeVTK_ControlMesh(const vector<Vertex2D> &spt, const vector<Element2D> &mesh, int step, string fn, vector<float> var, string varName)
{
	string fname;
	stringstream ss;
	ss << std::setw(6) << std::setfill('0') << step;
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
		for (uint i = 0; i < var.size(); i++)
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

void NeuronGrowth::ConcentrationCal_Coupling_Bezier(float u, float v, const Element2D &bzel, float pt[2], float &disp, float dudx[2], float &detJ)
{
	float dUdx[2][2];
	vector<float> Nx(bzel.IEN.size());
	vector<array<float, 2>> dNdx(bzel.IEN.size());
	bzel.Para2Phys(u, v, pt);
	BasisFunction(u, v, bzel.pts, bzel.cmat, Nx, dNdx, dUdx, detJ);
	disp = 0.;
	dudx[0] = 0.;
	dudx[1] = 0.; // dudx[2] = 0.;
	for (uint i = 0; i < bzel.IEN.size(); i++)
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

	for (unsigned int e = 0; e < num_bzmesh_ele; e++)
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
		int nns[2] = {ns * ns, ns};
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
	ss << std::setw(6) << std::setfill('0') << step;
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
		for (uint i = 0; i < sdisp.size(); i++)
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

void NeuronGrowth::CalculateVarsForOutput(vector<array<float, 3>> &spt_all, vector<float> &sresult_all, vector<array<int, 4>> &sele_all) {
	float detJ;
	int num_bzmesh_ele = bzmesh_process.size();
	float spt_proc[num_bzmesh_ele * 4 * 3];
	float sresult_proc[num_bzmesh_ele * 4 * 1];
	int sele_proc[num_bzmesh_ele * 4];

	for (unsigned int e = 0; e < num_bzmesh_ele; e++)
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
		int nns[2] = {ns * ns, ns};
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
	vector<vector<array<float, 3>>> spt_all_4var; // sample points
	spt_all_4var.resize(4);
	vector<vector<float>> sresult_all_4var;
	sresult_all_4var.resize(4);
	vector<vector<array<int, 4>>> sele_all_4var;
	sele_all_4var.resize(4);
	N_0 = phi;
	CalculateVarsForOutput(spt_all_4var[0], sresult_all_4var[0], sele_all_4var[0]);
	N_0 = syn;
	// N_0 = theta;
	CalculateVarsForOutput(spt_all_4var[1], sresult_all_4var[1], sele_all_4var[1]);
	N_0 = tub;
	CalculateVarsForOutput(spt_all_4var[2], sresult_all_4var[2], sele_all_4var[2]);
	N_0 = tips;
	CalculateVarsForOutput(spt_all_4var[3], sresult_all_4var[3], sele_all_4var[3]);
	if (comRank == 0)
	{
		WriteVTK_All(spt_all_4var[0], sresult_all_4var, sele_all_4var[0], step, fn);
	}
}

void NeuronGrowth::WriteVTK_All(const vector<array<float, 3>> spt, const vector<vector<float>> sdisp, const vector<array<int, 4>> sele, int step, string fn)
{
	stringstream ss;
	// ss << step;
	ss << std::setw(6) << std::setfill('0') << step;
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
		for (uint i = 0; i < sdisp[0].size(); i++)
		{
			fout << sdisp[0][i] << "\n";
		}
		fout << "\nSCALARS " << "synaptogenesis " << "float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < sdisp[1].size(); i++)
		{
			fout << sdisp[1][i] << "\n";
		}
		fout << "\nSCALARS " << "tubulin " << "float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < sdisp[2].size(); i++)
		{
			fout << sdisp[2][i] << "\n";
		}
		fout << "\nSCALARS " << "tips " << "float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < sdisp[3].size(); i++)
		{
			fout << sdisp[3][i] << "\n";
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

	for (int j = 0; j < Nx.size(); j++)
		Value += U[j] * Nx[j];
}

void NeuronGrowth::PointFormGrad(vector<array<float, 2>> &dNdx, const vector<float> &U, float Value[2])
{
	for (int j = 0; j < 2; j++)
		Value[j] = 0.;

	for (int j = 0; j < 2; j++)
		for (int k = 0; k < dNdx.size(); k++)
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
			for (int l = 0; l < d2Ndx2.size(); l++)
			{
				Value[j][k] += U[l] * d2Ndx2[l][j][k];
			}
		}
	}
}

void NeuronGrowth::ElementValue(const vector<float> &Nx, const vector<float> value_node, float &value)
{
	value = 0.;
	for (int i = 0; i < Nx.size(); i++)
		value += value_node[i] * Nx[i];
}

void NeuronGrowth::ElementValueAll(const vector<float> &Nx, const vector<float> elePhiGuess, float &elePG, const vector<float> elePhi, float &eleP, const vector<float> eleSyn, float &eleS, const vector<float> eleTips, float &eleTp, const vector<float> eleTubulin, float &eleTb, const vector<float> eleEpsilon, float &eleEP, const vector<float> eleEpsilonP, float &eleEEP)
{
	elePG = 0.;
	eleP = 0.;
	eleS = 0.;
	eleTp = 0.;
	eleTb = 0.;
	eleEP = 0.;
	eleEEP = 0.;
	for (int i = 0; i < Nx.size(); i++) {
		elePG += elePhiGuess[i] * Nx[i];
		eleP += elePhi[i] * Nx[i];
		eleS += eleSyn[i] * Nx[i];
		eleTp += eleTips[i] * Nx[i];
		eleTb += eleTubulin[i] * Nx[i];
		eleEP += eleEpsilon[i] * Nx[i];
		eleEEP += eleEpsilonP[i] * Nx[i];
	}
}

void NeuronGrowth::ElementDeriv(const int nen, vector<array<float, 2>> &dNdx, const vector<float> value_node, float &dVdx, float &dVdy)
{
	dVdx = 0;
	dVdy = 0.;
	for (int i = 0; i < nen; i++) {
		dVdx += value_node[i] * dNdx[i][0];
		dVdy += value_node[i] * dNdx[i][1];
	}
}


void NeuronGrowth::ElementDerivAll(const int nen, vector<array<float, 2>> &dNdx, const vector<float> elePhiGuess, float &dPGdx, float &dPGdy, const vector<float> eleTheta, float &dThedx, float &dThedy, const vector<float> eleEpsilon, float &dAdx, float &dAdy, const vector<float> eleEpsilonP, float &dAPdx, float &dAPdy)
{	
	dPGdx = 0; dPGdy = 0.;
	dThedx = 0; dThedy = 0.;
	dAdx = 0; dAdy = 0.;
	dAPdx = 0; dAPdy = 0.;
	for (int i = 0; i < nen; i++) {
		dPGdx += elePhiGuess[i] * dNdx[i][0];	dPGdy += elePhiGuess[i] * dNdx[i][1];
		dThedx += eleTheta[i] * dNdx[i][0];	dThedy += eleTheta[i] * dNdx[i][1];
		dAdx += eleEpsilon[i] * dNdx[i][0];	dAdy += eleEpsilon[i] * dNdx[i][1];
		dAPdx += eleEpsilonP[i] * dNdx[i][0];	dAPdy += eleEpsilonP[i] * dNdx[i][1];
	}
}

void NeuronGrowth::ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 2>> &dNdx, const vector<float> elePhiGuess, float &elePG, const vector<float> elePhi, float &eleP, const vector<float> eleEpsilon, float &eleEP, const vector<float> eleEpsilonP, float &eleEEP, float &dPGdx, float &dPGdy, float &dAdx, float &dAdy, float &dAPdx, float &dAPdy)
{
	elePG = 0.;
	eleP = 0.;
	eleEP = 0.;
	eleEEP = 0.;

	dPGdx = 0; dPGdy = 0.;
	dAdx = 0; dAdy = 0.;
	dAPdx = 0; dAPdy = 0.;

	for (int i = 0; i < nen; i++) {
		elePG += elePhiGuess[i] * Nx[i];
		eleP += elePhi[i] * Nx[i];
		eleEP += eleEpsilon[i] * Nx[i];
		eleEEP += eleEpsilonP[i] * Nx[i];

		dPGdx += elePhiGuess[i] * dNdx[i][0];	dPGdy += elePhiGuess[i] * dNdx[i][1];
		dAdx += eleEpsilon[i] * dNdx[i][0];	dAdy += eleEpsilon[i] * dNdx[i][1];
		dAPdx += eleEpsilonP[i] * dNdx[i][0];	dAPdy += eleEpsilonP[i] * dNdx[i][1];
	}
}


void NeuronGrowth::ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 2>> &dNdx, vector<vector<float>> &eleVal, vector<float> &vars)
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

	for (int i = 0; i < nen; i++) {
		vars[0] += eleVal[0][i] * Nx[i];
		vars[3] += eleVal[1][i] * Nx[i];
		vars[4] += eleVal[3][i] * Nx[i];
		vars[5] += eleVal[4][i] * Nx[i];

		vars[1] += eleVal[0][i] * dNdx[i][0];	vars[2] += eleVal[0][i] * dNdx[i][1];
		vars[6] += eleVal[3][i] * dNdx[i][0];	vars[7] += eleVal[3][i] * dNdx[i][1];
	}
}

// void NeuronGrowth::ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 2>> &dNdx, vector<vector<float>> &eleVal, vector<float> &vars)
void NeuronGrowth::ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 2>> &dNdx, vector<float> &elePhiGuess, vector<float> &vars)
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

	for (int i = 0; i < nen; i++) {
		vars[0] += elePhiGuess[i] * Nx[i];
		// vars[3] += eleVal[1][i] * Nx[i];
		// vars[4] += eleEpsilon[i] * Nx[i];
		// vars[5] += eleVal[4][i] * Nx[i];

		vars[1] += elePhiGuess[i] * dNdx[i][0];	vars[2] += elePhiGuess[i] * dNdx[i][1];
		// vars[6] += eleEpsilon[i] * dNdx[i][0];	vars[7] += eleEpsilon[i] * dNdx[i][1];
	}
}

void NeuronGrowth::ElementEvaluationAll_syn_tub(const int nen, const vector<float> &Nx, vector<array<float, 2>> &dNdx, const vector<float> elePhiDiff, float &elePf, const vector<float> eleSyn, float &eleS, const vector<float> elePhi, float &eleP, const vector<float> elePhiPrev, float &elePprev, const vector<float> eleConct, float &eleC, float &dPdx, float &dPdy)
{
	elePf = 0.;
	eleS = 0.;
	eleP = 0.;
	elePprev = 0.;
	eleC = 0.;

	dPdx = 0; dPdy = 0.;

	for (int i = 0; i < nen; i++) {
		elePf += elePhiDiff[i] * Nx[i];
		eleS += eleSyn[i] * Nx[i];

		eleP += elePhi[i] * Nx[i];
		elePprev += elePhiPrev[i] * Nx[i];
		eleC += eleConct[i] * Nx[i];

		dPdx += elePhi[i] * dNdx[i][0];	dPdy += elePhi[i] * dNdx[i][1];
	}
}

void NeuronGrowth::ElementEvaluationAll_syn_tub(const int nen, const vector<float> &Nx, vector<array<float, 2>> &dNdx, vector<vector<float>> &eleVal, vector<float> &vars)
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

	for (int i = 0; i < nen; i++) {
		vars[0] += eleVal[0][i] * Nx[i];
		vars[1] += eleVal[1][i] * Nx[i];

		vars[2] += eleVal[2][i] * Nx[i];
		vars[5] += eleVal[3][i] * Nx[i];
		vars[6] += eleVal[4][i] * Nx[i];

		vars[3] += eleVal[2][i] * dNdx[i][0];	vars[4] += eleVal[2][i] * dNdx[i][1];
	}
}

void NeuronGrowth::prepareBasis() {

	/*Build linear system in each process*/
	float detJ;
	float dudx[2][2];
	
	float dThedx, dThedy;
	sum_grad_phi0_local = 0;

	float eleP0(0);
	int e;
	// std::cout << comRank << " " << bzmesh_process.size() << std::endl;

	for (e = 0; e < bzmesh_process.size(); e++) {
		int nen = bzmesh_process[e].IEN.size();

		// if (comRank == 1)
			// std::cout << comRank << " " << e << "/" << bzmesh_process.size() << std::endl;

		elePhi0.resize(nen);
		eleTheta.resize(nen);
		
		for (int i = 0; i < nen; i++) {
			elePhi0[i] = phi_0[bzmesh_process[e].IEN[i]];
			eleTheta[i] = theta[bzmesh_process[e].IEN[i]];
		}

		for (int i = 0; i < Gpt.size(); i++) {
			for (int j = 0; j < Gpt.size(); j++) {
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

void NeuronGrowth::preparePhaseField() {

	pre_eleEP.clear();
	pre_eleEEP.clear();
	pre_dAdx.clear();
	pre_dAdy.clear();
	pre_eleP.clear();
	// pre_dPdx.clear();
	// pre_dPdy.clear();
	pre_eleTh.clear();
	pre_eleMp.clear();
	pre_C1.clear();

	float eleEP(0), eleEEP(0), dAdx(0), dAdy(0), 
		eleP(0), eleTh(0), eleS(0), eleTb(0), eleTp(0), eleE(0), eleMp(0); 
	int e, ind(0);

	for (e = 0; e < bzmesh_process.size(); e++) {
		int nen = bzmesh_process[e].IEN.size();

		vector<float> elePhi(nen, 0), eleTheta(nen, 0), eleEpsilon(nen, 0), eleEpsilonP(nen, 0);
		vector<float> eleSyn(nen, 0), eleTubulin(nen, 0), eleTips(nen, 0), eleMphi(nen, 0);
		// elePhi.resize(nen);
		// eleTheta.resize(nen);
		// eleEpsilon.resize(nen);
		// eleEpsilonP.resize(nen);
		// eleSyn.resize(nen);
		// eleTubulin.resize(nen);
		// eleTips.resize(nen);
		// eleMphi.resize(nen);
		
		// std::cout << Mphi.size() << std::endl;

		for (int i = 0; i < nen; i++) {

			elePhi[i] = phi[bzmesh_process[e].IEN[i]];			// elePhi
			eleTheta[i] = theta[bzmesh_process[e].IEN[i]];			// eleTheta
			// eleEpsilon[i] = 0;						// eleEpsilon
			// eleEpsilonP[i] = 0;						// eleEpsilonP

			eleSyn[i] = syn[bzmesh_process[e].IEN[i]];			// eleSyn
			eleTubulin[i] = tub[bzmesh_process[e].IEN[i]];			// eleTubulin
			eleTips[i] = tips[bzmesh_process[e].IEN[i]];			// eleTips
			// if (n < 1) {
			// 	eleMphi[i] = 60;
			// } else {
			// 	eleMphi[i] = Mphi[bzmesh_process[e].IEN[i]];
			// }
		}

		for (int i = 0; i < Gpt.size(); i++) {
			for (int j = 0; j < Gpt.size(); j++) {

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
				// ElementDeriv(nen, pre_Nx[ind], elePhi, dPdx, dPdy);
				// pre_dPdx.push_back(dPdx);
				// pre_dPdy.push_back(dPdy);

				ElementValue(pre_Nx[ind], eleTheta, eleTh);
				pre_eleTh.push_back(eleTh);

				ElementValue(pre_Nx[ind], eleSyn, eleS);
				ElementValue(pre_Nx[ind], eleTubulin, eleTb);
				ElementValue(pre_Nx[ind], eleTips, eleTp);

				// float eleMp;
				// ElementValue(pre_Nx[ind], eleMphi, eleMp);
				// pre_eleMp.push_back(eleMp);

				// adjust rg (assembly rate) and sg (disassembly rate) based on detected tips
				if (n < 1) {
					eleE = alphaOverPi*atan(gamma * (1 - eleS));
					pre_eleMp.push_back(50);
				} else {
					if (eleTp > 0) {
						eleE = alphaOverPi*atan(gamma * Regular_Heiviside_fun(50 * eleTb - 0) * (1 - eleS));
						if (eleTp > 2) {
							// pre_eleMp.push_back(120);
							pre_eleMp.push_back(50);
						} else {
							// pre_eleMp.push_back(50);
							pre_eleMp.push_back(10);
						}
					} else {
						eleE = alphaOverPi*atan(gamma * Regular_Heiviside_fun(r * eleTb - g) * (1 - eleS));
						// pre_eleMp.push_back(5);
						pre_eleMp.push_back(1);
					}
				}

				// calculate C1 variable for phase field energy term
				float C1 = eleE - pre_C0[ind]; // C1 = E - (0.5 + 6 * user->s_coeff * sqrt(pow(dThedx, 2) + pow(dThedy, 2));
				pre_C1.push_back(C1);
										
				ind += 1;
			}
		}
	}
	// std::cout << pre_eleMp.size() << std::endl;
	// MPI_Barrier(PETSC_COMM_WORLD);
}

void NeuronGrowth::prepareSourceSum() {

	/*Build linear system in each process*/
	sum_grad_phi0_local = 0;

	float eleP(0), dPdx(0), dPdy(0);
	int e, ind(0);
	// std::cout << comRank << " " << bzmesh_process.size() << std::endl;

	vector<float> elePhi;
	for (e = 0; e < bzmesh_process.size(); e++) {
		int nen = bzmesh_process[e].IEN.size();

		elePhi.resize(nen);
		
		for (int i = 0; i < nen; i++) {
			elePhi[i] = phi[bzmesh_process[e].IEN[i]];
		}

		for (int i = 0; i < Gpt.size(); i++) {
			for (int j = 0; j < Gpt.size(); j++) {
				ElementDeriv(nen, pre_dNdx[ind], elePhi, dPdx, dPdy);

				pre_mag_grad_phi0.push_back(pow(dPdx, 2) + pow(dPdy, 2));
				sum_grad_phi0_local += pow(dPdx, 2) + pow(dPdy, 2);

				ind += 1;
			}
		}
	}
	// MPI_Barrier(PETSC_COMM_WORLD);
}

void NeuronGrowth::prepareTerm_source() {
	int ind(0);
	pre_term_source.clear();
	for (int e = 0; e < bzmesh_process.size(); e++) {
		for (int i = 0; i < Gpt.size(); i++) {
			for (int j = 0; j < Gpt.size(); j++) {
				pre_term_source.push_back(source_coeff * pre_mag_grad_phi0[ind] / sum_grad_phi0_global);
			}
		}
	}
}

void NeuronGrowth::prepareEE()
{
	/*Build linear system in each process*/
	PetscInt ind(0), e;
	for (e = 0; e < bzmesh_process.size(); e++) {
		PetscInt nen = bzmesh_process[e].IEN.size();

		for (PetscInt i = 0; i < nen * 1; i++) {
			EVectorSolve[i] = 0.0;
			for (PetscInt j = 0; j < nen * 1; j++) {
				EMatrixSolve[i][j] = 0.0;
			}
		}
	
		for (PetscInt i = 0; i < nen * 1; i++) {
			eleVal[0][i] = phi[bzmesh_process[e].IEN[i]];		// elePhiGuess
			// eleVal[1][i] = phi[bzmesh_process[e].IEN[i]];		// elePhi
			eleVal[0][i] = 0;		// elePhiGuess
			eleVal[1][i] = 0;		// elePhi

			eleVal[2][i] = syn[bzmesh_process[e].IEN[i]];		// eleSyn
			eleVal[3][i] = tub[bzmesh_process[e].IEN[i]];		// eleTubulin
			eleVal[4][i] = theta[bzmesh_process[e].IEN[i]];	// eleTheta
			eleVal[5][i] = tips[bzmesh_process[e].IEN[i]];	// eleTips
			eleVal[6][i] = 0;							// eleEpsilon
			eleVal[7][i] = 0;							// eleEpsilonP
		}

		for (PetscInt i = 0; i < Gpt.size(); i++) {
			for (PetscInt j = 0; j < Gpt.size(); j++) {

				PetscReal detJ;

				vector<float> Nx;
				vector<array<float, 2>> dNdx;
				Nx.clear();
				dNdx.clear();
				Nx.resize(nen);
				dNdx.resize(nen);

				// use pre-calculated values (save computational cost)
				Nx = pre_Nx[ind];
				dNdx = pre_dNdx[ind];
				detJ = pre_detJ[ind];
				vars[1] = pre_C0[ind]; // C0
				ind += 1;

				EvaluateOrientation(nen, Nx, dNdx, eleVal[1], eleVal[4], eleVal[6], eleVal[7]);
				ElementEvaluationAll_phi(nen, Nx, dNdx, eleVal, vars);
				//	 0   1    2     3      4      5     6     7      8     9      10      11      12     13    14    15      16     17     
				// float C1, C0, elePG, dPGdx, dPGdy, eleP, eleS, eleTb, eleE, eleTp, dThedx, dThedy, eleEP, dAdx, dAdy, eleEEP, dAPdx, dAPdy;
				
				vars[8] = alphaOverPi*atan(gamma * (1 - vars[6]));

				vars[0] = vars[8] - vars[1]; // E - (0.5 + 6 * s_coeff * sqrt(pow(dThedx, 2) + pow(dThedy, 2));
									
				for (PetscInt m = 0; m < nen; m++) {
						EVectorSolve[m] += (vars[2] * Nx[m] - dt * M_phi *\
						((- vars[12] * vars[12] * (vars[3] * dNdx[m][0] + vars[4] * dNdx[m][1])) -\
						(- vars[13] * vars[15] * vars[4] * dNdx[m][0]) +\
						(- vars[14] * vars[15] * vars[3] * dNdx[m][1]) +\
						(- vars[2] * vars[2] * vars[2] + (1 - vars[0]) * vars[2] * vars[2] + vars[0] * vars[2]) * Nx[m])
						- vars[5] * Nx[m]) * detJ;
					for (PetscInt n = 0; n < nen; n++) {
						EMatrixSolve[m][n] += (Nx[m] * Nx[n] - dt * M_phi *
							((- vars[12] * vars[12] * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1])) // terma2
							- (- vars[13] * vars[15] * dNdx[m][1] * dNdx[n][0]) // termadx
							+ (- vars[14] * vars[15] * dNdx[m][0] * dNdx[n][1]) // termady
							+ (- 3 * vars[2] * vars[2] + 2 * (1 - vars[0]) * vars[2] + vars[0] * Nx[m]) * Nx[n]) // termdbl
							) * detJ;
					}
				}
			}
		}

		/*Apply Boundary Condition*/
		for (PetscInt i = 0; i < nen; i++) {
			PetscInt A = bzmesh_process[e].IEN[i];
			if (cpts[A].label == 1) // domain boundary
				ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
		}
		
		pre_EMatrixSolve.push_back(EMatrixSolve);
		pre_EVectorSolve.push_back(EVectorSolve);
	}
}

void NeuronGrowth::EvaluateEnergy(const int nen, const vector<float> &Nx, const vector<float> eleSyn, vector<float>& E)
{
	for (int i = 0; i < nen; i++) {
		E[i] = alphaOverPi*atan(gamma*(1-eleSyn[i]));
	}

}

float NeuronGrowth::Regular_Heiviside_fun(float x) 
{
	// float epsilon = 0.0001; // the number is not fixed.
	float epsilon = 1e-15; // the number is not fixed.
	return 0.5*(1+(2/PI)*atan(x/epsilon));
}

void NeuronGrowth::EvaluateOrientation(const int nen, const vector<float> &Nx, const vector<array<float, 2>> &dNdx, const vector<float> elePhi, const vector<float> eleTheta,  vector<float>& eleEpsilon, vector<float>& eleEpsilonP)
{
	for (int i = 0; i < nen; i++) {
		eleEpsilon[i] += epsilonb * (1.0 + delta * cos(aniso * (atan2(elePhi[i] * dNdx[i][0], elePhi[i] * dNdx[i][1]) - eleTheta[i] * Nx[i])));
		eleEpsilonP[i] += -epsilonb * (aniso * delta * sin(aniso * (atan2(elePhi[i] * dNdx[i][0], elePhi[i] * dNdx[i][1]) - eleTheta[i] * Nx[i])));
	}
}

// void NeuronGrowth::BuildLinearSystemProcessNG_phi()
// {
// 	/*Build linear system in each process*/
// 	int ind(0); // pre-calculated variable index
// 	for (int e = 0; e < bzmesh_process.size(); e++) {
// 		int nen = bzmesh_process[e].IEN.size(); // 16 supporting cp for 2D case

// 		vector<vector<float>> EMatrixSolve;	EMatrixSolve.clear();
// 		vector<float> EVectorSolve;		EVectorSolve.clear();
// 		EMatrixSolve.resize(nen);
// 		EVectorSolve.resize(nen);

// 		vector<float> elePhiGuess(nen, 0);

// 		for (int i = 0; i < nen * 1; i++) {
// 			// EMatrixSolve[i].reserve(nen*sizeof(float));
// 			EMatrixSolve[i].resize(nen);
// 		}
		
// 		// preparing element stiffness matrix and vector, extract values from control mesh
// 		for (int i = 0; i < nen; i++) {

// 			EVectorSolve[i] = 0.0;
// 			for (int j = 0; j < nen; j++) {
// 				EMatrixSolve[i][j] = 0.0;
// 			}
// 			elePhiGuess[i] = phi[bzmesh_process[e].IEN[i]];		// elePhiGuess
// 		}

// 		// loop through gaussian quadrature points
// 		for (int i = 0; i < Gpt.size(); i++) {
// 			for (int j = 0; j < Gpt.size(); j++) {				
// 				ElementEvaluationAll_phi(nen, pre_Nx[ind], pre_dNdx[ind], elePhiGuess, vars);
				
// 				// loop through control points
// 				for (int m = 0; m < nen; m++) {
// 					// terma2 = - eleEP * eleEP * (dPGdx * dNdx[m][0] + dPGdy * dNdx[m][1]);
// 					// termadx = - dAdx * eleEEP * dPGdy * dNdx[m][0];
// 					// termady = - dAdy * eleEEP * dPGdx * dNdx[m][1];
// 					// termdbl = (- elePG * elePG * elePG + (1 - C1) * elePG * elePG + C1 * elePG) * Nx[m];
// 					// EVectorSolve[m] += (elePG * Nx[m] - dt * M_phi * (terma2 - termadx + termady + termdbl) - eleP * Nx[m]) * detJ;
			
// 					EVectorSolve[m] += ( - vars[0] * pre_Nx[ind][m] - dt * pre_eleMp[ind] *\
// 						(
// 						// (- pre_eleEP[ind] * pre_eleEP[ind] * (vars[1] * pre_dNdx[ind][m][0] + vars[2] * pre_dNdx[ind][m][1])) -\
// 						// (- pre_dAdx[ind] * pre_eleEEP[ind] * vars[2] * pre_dNdx[ind][m][0]) +\
// 						// (- pre_dAdy[ind] * pre_eleEEP[ind] * vars[1] * pre_dNdx[ind][m][1]) +\

// 						( - vars[0] * vars[0] * vars[0] + (1 - pre_C1[ind]) * vars[0] * vars[0] + pre_C1[ind] * vars[0]) * pre_Nx[ind][m])
// 						+ pre_eleP[ind] * pre_Nx[ind][m]) * pre_detJ[ind];

// 					// loop through 16 control points
// 					for (int n = 0; n < nen; n++) {
// 						// terma2 = - eleEP * eleEP * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1]);
// 						// termadx = - dAdx * eleEEP * dNdx[m][1] * dNdx[n][0];
// 						// termady = - dAdy * eleEEP * dNdx[m][0] * dNdx[n][1];
// 						// termdbl = (- 3 * elePG * elePG + 2 * (1 - C1) * elePG + C1 * Nx[m]) * Nx[n];
// 						// EMatrixSolve[m][n] += (Nx[m] * Nx[n] - dt * M_phi * (terma2 - termadx + termady + termdbl)) * detJ;
						
// 						EMatrixSolve[m][n] += ( - pre_Nx[ind][m] * pre_Nx[ind][n] - dt * pre_eleMp[ind] *
// 							(
// 							// (- pre_eleEP[ind] * pre_eleEP[ind] * (pre_dNdx[ind][m][0] * pre_dNdx[ind][n][0] + pre_dNdx[ind][m][1] * pre_dNdx[ind][n][1])) // terma2
// 							// - (- pre_dAdx[ind] * pre_eleEEP[ind] * pre_dNdx[ind][m][1] * pre_dNdx[ind][n][0]) // termadx
// 							// + (- pre_dAdy[ind] * pre_eleEEP[ind] * pre_dNdx[ind][m][0] * pre_dNdx[ind][n][1]) // termady
							
// 							( - 3 * vars[0] * vars[0] + 2 * (1 - pre_C1[ind]) * vars[0] + pre_C1[ind] * pre_Nx[ind][m]) * pre_Nx[ind][n]) // termdbl
// 							) * pre_detJ[ind];
// 					}
// 				}
// 				ind += 1; // incrementing index for extracting pre-calculated variables
// 			}
// 		}

// 		/*Apply Boundary Condition*/
// 		for (int i = 0; i < nen; i++) {
// 			int A = bzmesh_process[e].IEN[i];
// 			if (cpts[A].label == 1) { // domain boundary
// 				ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
// 			}
// 			// if (tips[A] == 0) { // domain boundary
// 			// 	ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
// 			// }
// 		}

// 		ResidualAssembly(EVectorSolve, bzmesh_process[e].IEN, GR_phi);
// 		MatrixAssembly(EMatrixSolve, bzmesh_process[e].IEN, GK_phi);
// 	}

// 	if (comRank == 0) {
// 		check_itr += 1;
// 	}

// 	VecAssemblyBegin(GR_phi);
// 	MatAssemblyBegin(GK_phi, MAT_FINAL_ASSEMBLY);
// }

void NeuronGrowth::BuildLinearSystemProcessNG_phi()
{
	/*Build linear system in each process*/
	int ind(0); // pre-calculated variable index
	for (int e = 0; e < bzmesh_process.size(); e++) {
		int nen = bzmesh_process[e].IEN.size(); // 16 supporting cp for 2D case

		vector<vector<float>> EMatrixSolve;	EMatrixSolve.clear();
		vector<float> EVectorSolve;		EVectorSolve.clear();
		EMatrixSolve.resize(nen);
		EVectorSolve.resize(nen);

		vector<float> elePhiGuess(nen, 0);

		for (int i = 0; i < nen * 1; i++) {
			// EMatrixSolve[i].reserve(nen*sizeof(float));
			EMatrixSolve[i].resize(nen);
		}
		
		// preparing element stiffness matrix and vector, extract values from control mesh
		for (int i = 0; i < nen; i++) {

			EVectorSolve[i] = 0.0;
			for (int j = 0; j < nen; j++) {
				EMatrixSolve[i][j] = 0.0;
			}
			elePhiGuess[i] = phi[bzmesh_process[e].IEN[i]];		// elePhiGuess
		}

		// loop through gaussian quadrature points
		for (int i = 0; i < Gpt.size(); i++) {
			for (int j = 0; j < Gpt.size(); j++) {				
				ElementEvaluationAll_phi(nen, pre_Nx[ind], pre_dNdx[ind], elePhiGuess, vars);
				
				// loop through control points
				for (int m = 0; m < nen; m++) {
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
					for (int n = 0; n < nen; n++) {
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
		for (int i = 0; i < nen; i++) {
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

void NeuronGrowth::BuildLinearSystemProcessNG_syn(const vector<Element2D> &tmesh, const vector<Vertex2D> &cpts)
{
	/*Build linear system in each process*/
	for (int e = 0; e < bzmesh_process.size(); e++) {
		float detJ;
		float dudx[2][2];
		int nen = bzmesh_process[e].IEN.size();

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

		for (int i = 0; i < nen; i++) {
			EMatrixSolve[i].resize(nen * 1, 0.0);
		}

		for (int i = 0; i < nen * 1; i++) {
			for (int j = 0; j < nen * 1; j++)
			{
				EMatrixSolve[i][j] = 0.0;
			}
			EVectorSolve[i] = 0.0;
		}

		vector<float> elePhiDiff, eleSyn;
		elePhiDiff.resize(nen);
		eleSyn.resize(nen);
		float elePf, eleS;

		for (int i = 0; i < nen; i++) {
			elePhiDiff[i] = abs(phi[bzmesh_process[e].IEN[i]] - phi_prev[bzmesh_process[e].IEN[i]]);
			eleSyn[i] = syn[bzmesh_process[e].IEN[i]];
		}

		for (int i = 0; i < Gpt.size(); i++) {
			for (int j = 0; j < Gpt.size(); j++) {
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
		for (int i = 0; i < nen; i++) {
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

void NeuronGrowth::Tangent_syn(const int nen, vector<float> &Nx, vector<array<float, 2>> &dNdx, float detJ, vector<vector<float>> &EMatrixSolve)
{
	/*calculate tangent matrix*/
	for (int i = 0; i < nen; i++) {
		for (int j = 0; j < nen; j++) {
			EMatrixSolve[i][j] += (Nx[i] * Nx[j] + dt * Dc * (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1])) * detJ;
			// EMatrixSolve[i][j] += (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1]) * detJ;	// steady state heat transfer
		}
	}
}

void NeuronGrowth::Residual_syn(const int nen, vector<float> &Nx, vector<array<float, 2>> &dNdx, float detJ, float elePf, float eleS, vector<float> &EVectorSolve)
{
	/*calculate residual of the equation*/
	for (int i = 0; i < nen; i++) {
		EVectorSolve[i] += Nx[i] * (eleS + kappa * elePf) * detJ;
		// EVectorSolve[i] += Nx[i] * (eleS) * detJ; // transient heat transfer
		// EVectorSolve[i] += 0.00; // steady state heat transfer
	}
}

void NeuronGrowth::CalculateSumGradPhi0(const vector<Element2D> &tmesh, const vector<Vertex2D> &cpts)
{
	/*Build linear system in each process*/
	for (int e = 0; e < bzmesh_process.size(); e++) {
		float detJ;
		float dudx[2][2];
		int nen = bzmesh_process[e].IEN.size();

		vector<float> Nx;
		vector<array<float, 2>> dNdx;

		elePhi0.resize(nen);
		for (int i = 0; i < nen; i++) {
			elePhi0[i] = phi_0[bzmesh_process[e].IEN[i]];
		}

		for (int i = 0; i < Gpt.size(); i++) {
			for (int j = 0; j < Gpt.size(); j++) {
				BasisFunction(Gpt[i], Gpt[j], bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dudx, detJ);
				detJ = wght[i] * wght[j] * detJ;
				ElementDeriv(nen, dNdx, elePhi0, dP0dx, dP0dy);
				sum_grad_phi0_local += pow(dP0dx, 2) + pow(dP0dy, 2);
			}
		}
	}
}

void NeuronGrowth::BuildLinearSystemProcessNG_tub(const vector<Element2D> &tmesh, const vector<Vertex2D> &cpts)
{
	/*Build linear system in each process*/
	for (int e = 0; e < bzmesh_process.size(); e++) {
		float detJ;
		float dudx[2][2];
		int nen = bzmesh_process[e].IEN.size();

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

		for (int i = 0; i < nen; i++) {
			EMatrixSolve[i].resize(nen, 0.0);
		}

		for (int i = 0; i < nen * 1; i++) {
			for (int j = 0; j < nen * 1; j++) {
				EMatrixSolve[i][j] = 0.0;
			}
			EVectorSolve[i] = 0.0;
		}

		vector<float> elePhi, elePhiPrev, eleConct;
		elePhi.resize(nen);
		elePhiPrev.resize(nen);
		eleConct.resize(nen);
		float eleP, dPdx, dPdy, elePprev, eleC, mag_grad_phi0;

		for (int i = 0; i < nen; i++) {
			elePhi[i] = phi[bzmesh_process[e].IEN[i]];
			elePhiPrev[i] = phi_prev[bzmesh_process[e].IEN[i]];
			eleConct[i] = tub[bzmesh_process[e].IEN[i]];
		}

		for (int i = 0; i < Gpt.size(); i++) {
			for (int j = 0; j < Gpt.size(); j++) {
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
		for (int i = 0; i < nen; i++) {
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

void NeuronGrowth::Tangent_tub(const int nen, vector<float> &Nx, vector<array<float, 2>> &dNdx, float detJ, float eleC, float eleP, float dPdx, float dPdy, float elePprev, float mag_grad_phi0, vector<vector<float>> &EMatrixSolve)
{
	float term_diff, term_alph, term_beta;
	/*calculate tangent matrix*/
	for (int i = 0; i < nen; i++) {
		for (int j = 0; j < nen; j++) {
			term_diff = - Diff * (eleP * (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1]));
			term_alph = alphaT * (eleP * (dNdx[i][0] + dNdx[i][1]) + Nx[i] * (dPdx + dPdy)) * Nx[j];
			term_beta = betaT * (eleP * Nx[i]) * Nx[j];
			
			EMatrixSolve[i][j] += (Nx[i] * Nx[j] - dt / eleP * (term_diff - term_alph - term_beta)) * detJ;
			// EMatrixSolve[i][j] += ((2 * eleP - elePprev) * Nx[i] * Nx[j] - dt/2 * (term_diff - term_alph - term_beta + term_source)) * detJ;
		}
	}
}

void NeuronGrowth::Residual_tub(const int nen, vector<float> &Nx, vector<array<float, 2>> &dNdx, float detJ, float eleC, float eleP, float elePprev, float mag_grad_phi0, vector<float> &EVectorSolve)
{
	float term_source;
	/*calculate residual of the equation*/
	for (int i = 0; i < nen; i++) {
		term_source = source_coeff * mag_grad_phi0 / sum_grad_phi0_global;
		EVectorSolve[i] += (dt / eleP * term_source + eleC) * Nx[i] * detJ;
		// EVectorSolve[i] += eleP * eleC * Nx[i] * detJ;
	}
}

void NeuronGrowth::BuildLinearSystemProcessNG_syn_tub(const vector<Element2D> &tmesh, const vector<Vertex2D> &cpts)
{
	/*Build linear system in each process*/
	int ind(0);
	for (int e = 0; e < bzmesh_process.size(); e++) {
		int nen = bzmesh_process[e].IEN.size();

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

		for (int i = 0; i < nen * 1; i++) {
			EMatrixSolve_syn[i].resize(nen * 1, 0.0);
			EMatrixSolve_tub[i].resize(nen, 0.0);
			for (int j = 0; j < nen * 1; j++)
			{
				EMatrixSolve_syn[i][j] = 0.0;
				EMatrixSolve_tub[i][j] = 0.0;
			}
			EVectorSolve_syn[i] = 0.0;
			EVectorSolve_tub[i] = 0.0;

			if (i < 5) 
				eleVal_st[i].resize(nen);
		}
		
		for (int i = 0; i < nen * 1; i++) {
			eleVal_st[0][i] = phi[bzmesh_process[e].IEN[i]] - phi_prev[bzmesh_process[e].IEN[i]];		// elePhiDiff
			eleVal_st[1][i] = syn[bzmesh_process[e].IEN[i]];						// eleSyn
			eleVal_st[2][i] = phi[bzmesh_process[e].IEN[i]];						// elePhi
			eleVal_st[3][i] = phi_prev[bzmesh_process[e].IEN[i]];						// elePhiPrev
			eleVal_st[4][i] = tub[bzmesh_process[e].IEN[i]];						// eleConct
		}


		for (int i = 0; i < Gpt.size(); i++) {
			for (int j = 0; j < Gpt.size(); j++) {
				
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

				for (int m = 0; m < nen; m++) {
					EVectorSolve_syn[m] += Nx[m] * (vars_st[1] + kappa * vars_st[0]) * detJ;

					// EVectorSolve_tub[m] += (dt / vars_st[2] * vars_st[11] + vars_st[6]) * Nx[m] * detJ;
					// EVectorSolve_tub[m] += (dt / vars_st[2] + vars_st[6]) * Nx[m] * detJ;
					EVectorSolve_tub[m] += (vars_st[2] * vars_st[6] + (dt/4) * vars_st[11]) * Nx[m] * detJ;

					for (int n = 0; n < nen; n++) {
						if (judge_syn == 0)
							EMatrixSolve_syn[m][n] += (Nx[m] * Nx[n] + dt/4 * Dc * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1])) * detJ;

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
		for (int i = 0; i < nen; i++) {
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

int NeuronGrowth::CheckExpansion(vector<float> input) {

	int length = input.size();
	int sz = sqrt(length);

	// 0 - left | 1 - top | 2 - right | 3 - bottom | 4 - no action
	for (int i = 0; i < 2; i++) {
		for (int j = 1; j < (sz-1); j++) {
			if (round(input[i * sz + j]) > 0) {
				return 0;
			} else if (round(input[j * sz + i]) > 0) {
				return 3;
			} else if (round(input[(sz - i - 1) * sz + j]) > 0) {
				return 2;
			} else if (round(input[j * sz - i - 1]) > 0) {
				return 1;
			}
		}
	}
	return 4;
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

void NeuronGrowth::ExpandDomain(vector<float> input, vector<float> &expd_var, int edge) 
{
	int length = input.size();
	int sz = sqrt(length);

	int expd_sz = 10; // directional expanding size
	int new_sz = sz + expd_sz;

	expd_var.clear(); expd_var.resize(pow(new_sz,2));
	for (int i = 0; i < new_sz; i++)
		expd_var[i] = 0;
	
	int ind;
	switch (edge) {
		case 0: // left
			ind = 0;
			for (int i = expd_sz; i < new_sz-1; i++) {
				for (int j = (int)(expd_sz/2); j < (new_sz - (int)(expd_sz/2)); j++) {
					expd_var[i * new_sz + j] = input[ind];
					ind += 1;
				}
			}
			// std::cout << "Expanding left!" << std::endl;
			break;
		case 1: // top
			ind = 0;
			for (int i = (int)(expd_sz/2); i < (new_sz - (int)(expd_sz/2)); i++) {
				for (int j = 0; j < new_sz-expd_sz; j++) {
					expd_var[i * new_sz + j] = input[ind];
					ind += 1;
				}
			}
			// std::cout << "Expanding top!" << std::endl;
			break;
		case 2: // right
			ind = 0;
			for (int i = 0; i < new_sz-1-expd_sz; i++) {
				for (int j = (int)(expd_sz/2); j < (new_sz - (int)(expd_sz/2)); j++) {
					expd_var[i * new_sz + j] = input[ind];
					ind += 1;
				}
			}
			// std::cout << "Expanding right!" << std::endl;
			break;
		case 3: // bottom
			ind = 0;
			for (int i = (int)(expd_sz/2); i < (new_sz - (int)(expd_sz/2)); i++) {
				for (int j = expd_sz; j < new_sz; j++) {
					expd_var[i * new_sz + j] = input[ind];
					ind += 1;
				}
			}
			// std::cout << "Expanding bottom!" << std::endl;
			break;
		case 4: // all direction
			ind = 0;
			for (int i = (int)(expd_sz/2); i < (new_sz - (int)(expd_sz/2)); i++) {
				for (int j = (int)(expd_sz/2); j < (new_sz - (int)(expd_sz/2)); j++) {
					expd_var[i * new_sz + j] = input[ind];
					ind += 1;
				}
			}
			// std::cout << "Expanding all direction!" << std::endl;
			break;
	}
	// input.swap(expd_var);
}

void NeuronGrowth::ExpandDomain(vector<float> input, vector<float> &expd_var, int edge, int NX, int NY) 
{
	int expd_sz = 10; // directional expanding size
	
	expd_var.clear(); expd_var.resize((NX+1) * (NY+1));
	for (int i = 0; i < (expd_var.size()); i++)
		expd_var[i] = 0;
	
	// (0-left|1-top|2-right|3-bottom)
	int ind;
	switch (edge) {
		case 0: // left
			ind = 0;
			for (int i = expd_sz; i < NX; i++) {
				for (int j = 0; j <= NY; j++) {
					expd_var[i * (NY+1) + j] = input[ind];
					ind += 1;
				}
			}
			// std::cout << "Expanding left!" << std::endl;
			break;
		case 1: // top
			ind = 0;
			for (int i = 0; i < NX; i++) {
				for (int j = 0; j <= NY-expd_sz; j++) {
					expd_var[i * (NY+1) + j] = input[ind];
					ind += 1;
				}
			}
			// std::cout << "Expanding top!" << std::endl;
			break;
		case 2: // right
			ind = 0;
			for (int i = 0; i < NX-expd_sz; i++) {
				for (int j = 0; j <= NY; j++) {
					expd_var[i * (NY+1) + j] = input[ind];
					ind += 1;
				}
			}
			// std::cout << "Expanding right!" << std::endl;
			break;
		case 3: // bottom
			ind = 0;
			for (int i = 0; i < NX; i++) {
				for (int j = expd_sz; j <= NY; j++) {
					expd_var[i * (NY+1) + j] = input[ind];
					ind += 1;
				}
			}
			// std::cout << "Expanding bottom!" << std::endl;
			break;
	}
}

void NeuronGrowth::PopulateRandom(vector<float> &input) {

	int length = input.size();
	for (int i = 0; i < length; i++) {
		if (input[i] == 0)
			input[i] = (float)(rand()%100)/(float)100;
	}
}


// Function to interpolate values for new control points
std::vector<float> NeuronGrowth::InterpolateValues_closest(
	const std::vector<float>& phi_in,
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
		// if (abs(round(phi_in[closestIndex])) >= 0.5) {
		// 	interpolatedValues[i] = 1;
		// } else {
		// 	interpolatedValues[i] = 0;
		// }
		interpolatedValues[i] = phi_in[closestIndex];
	}

	return interpolatedValues;
}

vector<float> NeuronGrowth::InterpolateValues_closest(const vector<float>& input, const KDTree& kdTree, const vector<Vertex2D>& cpt_out) {

	std::vector<float> interpolatedValues(cpt_out.size());
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

vector<float> NeuronGrowth::InterpolateVars_coarse1(vector<float> input, vector<Vertex2D> cpts_initial, const KDTree& kdTree_initial, const vector<Vertex2D>& cpts, int type, int isTheta) 
{	
	vector<float> output;
	output.resize(cpts.size(), 0);
	
	for (int i = 0; i < cpts.size(); i++) {
		float x = round5(cpts[i].coor[0]);
		float y = round5(cpts[i].coor[1]);
		int ind;
		if (KD_SearchPair(cpts_initial, kdTree_initial, x, y, ind)) {
			output[i] = input[ind];
		} 
		else {
			int indDown, indUp, indLeft, indRight;
			if ((abs(remainder(x,2)) == 1) && (abs(remainder(y,2)) != 1)
			&& (x <= prev_max_x) && (x >= prev_min_x) && (y <= prev_max_y) && (y >= prev_max_y)) {
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
			} else if ((abs(remainder(x,2)) != 1) && (abs(remainder(y,2)) == 1)
			&& (x <= prev_max_x) && (x >= prev_min_x) && (y <= prev_max_y) && (y >= prev_max_y)) {
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
			} else if ((abs(remainder(x,2)) == 1) && (abs(remainder(y,2)) == 1)
			&& (x <= prev_max_x) && (x >= prev_min_x) && (y <= prev_max_y) && (y >= prev_max_y)) {
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
			} else {
				if (isTheta != 1) {
					output[i] = 0;
				} else {
					output[i] = (float)(rand()%100)/(float)100; // for theta
				}
			}
		}			
	}
	return output;
}

bool NeuronGrowth::KD_SearchPair(const vector<Vertex2D> prev_cpts, const KDTree& kdTree, float targetX, float targetY, int &ind) {

	const float query_pt[2] = {targetX, targetY};
	size_t closestIndex;
	float out_dist_sqr;
	nanoflann::KNNResultSet<float> resultSet(1);
	resultSet.init(&closestIndex, &out_dist_sqr);

	// Performing the search with correct query points
	kdTree.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(0));

	float x = prev_cpts[closestIndex].coor[0];
	float y = prev_cpts[closestIndex].coor[1];
	
	if ((max(abs(x-targetX), abs(y-targetY)) <= 1) || (targetX < prev_min_x) || (targetY < prev_min_y) || (targetX > prev_max_x) || (targetY > prev_max_y)) {
		ind = static_cast<int>(closestIndex);
		return true; // Found the pair (targetX, targetY) in the vector
	} else {
		PetscPrintf(PETSC_COMM_WORLD, "Failed search! x: %f y: %f | %f %f\n", x, y, targetX, targetY);
		return false; // Pair not found in the vector
	}
}

float NeuronGrowth::RmOutlier(vector<float> &data) 
{
	float sum = 0.0, mean, standardDeviation = 0.0;

	for(int i = 0; i < data.size(); i++) {
		sum += data[i];
	}

	mean = sum / data.size();

	for(int i = 0; i < data.size(); i++) {
		standardDeviation += pow(data[i] - mean, 2);
	}

	standardDeviation = sqrt(standardDeviation / data.size());

	// // std::cout << mean + 2*standardDeviation << " " << maxVal << std::endl;
	// for(int i = 0; i < data.size(); i++) {
	// 	if (data[i] > (mean + 3.29*standardDeviation)) {
	// 		// data[i] = mean + 3.29*standardDeviation;
	// 	} else {
	// 		data[i] = 0;
	// 	}
	// }
	return (mean + 3.29*standardDeviation);
}

float NeuronGrowth::CellBoundary(float phi, float threshold) {
	float P;
	if (phi > threshold) {
		P = 1;
	} else {
		P = 0;
	}
	return P;
}

// Function to apply a simple smoothing operation to a 2D binary variable
std::vector<float> NeuronGrowth::SmoothBinary2D(const std::vector<float>& binaryData, int rows, int cols) {
	// Create a new vector to store the smoothed data
	std::vector<float> smoothedData(rows * cols, 0.0f);

	// Define a 5x5 smoothing kernel
	const std::vector<float> kernel = {
		1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
		1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
		1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
		1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
	1.0f, 1.0f, 1.0f, 1.0f, 1.0f
	};

	// Apply the smoothing operation using the defined kernel
	for (int i = 2; i < rows - 2; ++i) {
		for (int j = 2; j < cols - 2; ++j) {
			float sum = 0.0f;
			for (int k = -2; k <= 2; ++k) {
				for (int l = -2; l <= 2; ++l) {
					sum += binaryData[(i + k) * cols + (j + l)] * kernel[(k + 2) * 5 + (l + 2)];
				}
			}
			smoothedData[i * cols + j] = (sum > 12.0f) ? 1.0f : 0.0f; // Apply threshold to convert sum to binary value
		}
	}

	return smoothedData;
}


// Check if point coordinates are within the bounds defined by center and dx dy dz
bool NeuronGrowth::isInBox(const Vertex2D& point, const Vertex2D& center, float dx, float dy) {
	return point.coor[0] >= (center.coor[0] - dx/2) && point.coor[0] <= (center.coor[0] + dx/2) &&
		point.coor[1] >= (center.coor[1] - dy/2) && point.coor[1] <= (center.coor[1] + dy/2);
}

/**
 * Calculates the sum of phi values around each vertex in a given set, applying
 * a threshold to identify significant points, potentially indicating neuron growth tips.
 * 
 * This function iterates over a collection of vertices (`cpts`), summing the phi values
 * within a defined vicinity around each vertex. The vicinity is determined by the `dx`, and `dy`
 * parameters, which define the dimensions of a box centered on each vertex. Points within
 * this box contribute to the sum. After calculating the sums, the function applies a threshold
 * to these values, normalizing them against the highest value found. Points with summed values
 * above this threshold are marked as potential growth tips by setting their corresponding value
 * in the output vector to 1; all others are set to 0.
 * 
 * @param cpts A vector of Vertex2D objects representing the neuron's vertices.
 * @param dx The delta in the x-direction to define the vicinity around a point.
 * @param dy The delta in the y-direction to define the vicinity around a point.
 * @return std::vector<float> A vector of the same size as `cpts`, where each element is either 0
 *         (indicating the corresponding vertex is not a tip) or 1 (indicating a potential tip),
 *         based on the thresholding of summed phi values.
 * 
 * Note: The function assumes that `phi` is a pre-defined vector accessible within the class
 *       that contains phi values corresponding to each Vertex3D in `cpts`. The function
 *       `isInBox` checks whether a point is within the specified vicinity of another,
 *       and `CellBoundary` computes a boundary-related value for a given phi, with
 *       the second argument presumably allowing for further customization.
 */
vector<float> NeuronGrowth::calculatePhiSum(const std::vector<Vertex2D>& cpts, float dx, float dy, vector<float> id) {
	std::vector<float> tp(cpts.size(), 0); // Initialize the result vector with zeros
	tp.clear(); tp.resize(cpts.size());
	float threshold = 0.9; // Threshold for filtering tp
	float maxVal = 0; // Track the maximum value of tp for normalization

	// Iterate over each center point
	for (int i = 0; i < cpts.size(); ++i) {
		const auto& center = cpts[i];

		tp[i] = 0;

		// Sum phi values for points within the box centered at 'center'
		for (int j = 0; j < phi.size(); ++j) {
			if (isInBox(cpts[j], center, dx, dy)) {
				tp[i] += CellBoundary(phi[j], 0.5), id[j];
			}
		}

		tp[i] = CellBoundary(phi[i], 0.5) / tp[i] * CellBoundary(phi[i], 0.5);
		// Handle NaN cases, ensuring a valid tp value
		if (std::isnan(tp[i])) {
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

	for (int i = 0; i < tp.size(); i++) {
		tp[i] = tp[i]/maxVal;
	}

	// maxVal = RmOutlier(tp);

	// Thresholding and setting tp values based on the maximum value found
	for (int i = 0; i < tp.size(); ++i) {
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

void NeuronGrowth::DetectTipsMulti(const vector<float>& phi_fine, const vector<float>& id, const int& numNeuron, vector<float>& phiSum, const int& NX, const int& NY)
{
	// float threshold(0.9997), maxVal(0);
	float threshold(0.9), stdVal(0), maxVal(0);
	int ind, length((NX+1)*(NY+1));
	phiSum.clear();
	phiSum.resize(length, 0);

	for (int i = (5*NY+5); i < (length-4*NY-4); i++) {
		// phiSum[i] = 0;
		if (CellBoundary(phi_fine[i], 0.25) > 0) {
		// if (CellBoundary(phi_fine[i], 0.5) > 0) {
			// int baseID = static_cast<int>(round(id[i]));
			for (int j = -4; j < 5; j++) {
				for (int k = -4; k < 5; k++) {
			// for (int j = -6; j < 7; j++) {
			// 	for (int k = -6; k < 7; k++) {
					// for (int l = 0; l < numNeuron; l++) {
						// int inID = static_cast<int>(round(id[i+j*(NY+1)+k]));
						// std::cout << baseID << " " << inID << std::endl;
						if (round(id[i]) == round(id[i+j*(NY+1)+k])) {
						// if (inID == l+1) { // not the same neuron
						// if ((l+1) == static_cast<int>(round(id[i+j*(NY+1)+k]))) {
							// phiSum[i] += CellBoundary(phi_fine[i+j*(NY+1)+k], 0.1);
							phiSum[i] += CellBoundary(phi_fine[i+j*(NY+1)+k], 0.25);
							// phiSum[i] += CellBoundary(phi_fine[i+j*(NY+1)+k], 0.5);
						}
					// }
				}
			}	
			if (phiSum[i] > 0)
				phiSum[i] = CellBoundary(phi_fine[i], 0) / phiSum[i];
			if (std::isnan(phiSum[i]))
				phiSum[i] = 0;
			// if (phiSum[i] > maxVal)
			// 	maxVal = phiSum[i];
		}
	}

	stdVal = RmOutlier(phiSum);

	for (int i = 0; i < phiSum.size(); i++) {
		phiSum[i] = phiSum[i]/stdVal;
		maxVal = max(phiSum[i], maxVal);
	}

	// std::cout << maxVal << std::endl;

	for (int i = 1+NY; i < length-NY-1; i++) {
		// if (phiSum[i] < (threshold * maxVal)) {
		// if (phiSum[i] < (threshold)) {
		if (phiSum[i] < (1.3)) {
			phiSum[i] = 0;
		} 
		// else {
		// 	tip[i] = 1;
		// }
	}
}

void NeuronGrowth::DetectTipsMulti_new(const std::vector<float>& phi_fine, const std::vector<float>& id, int numNeuron, std::vector<float>& phiSum, int NX, int NY) {
	float threshold = 0.9997f;
	int length = (NX + 1) * (NY + 1);
	phiSum.clear();
	phiSum.resize(length, 0.0f);

	// Kernel for convolution; adjust size/shape for your specific needs
	const int kernelSize = 5; // 5x5 kernel
	const int kernelOffset = kernelSize / 2;

	for (int x = kernelOffset; x < NX - kernelOffset; ++x) {
		for (int y = kernelOffset; y < NY - kernelOffset; ++y) {
			int index = x * (NY + 1) + y;
			if (CellBoundary(phi_fine[index], 0.25) > 0) {
				for (int dx = -kernelOffset; dx <= kernelOffset; ++dx) {
					for (int dy = -kernelOffset; dy <= kernelOffset; ++dy) {
						int neighborIndex = index + dx * (NY + 1) + dy;
						if (static_cast<int>(round(id[neighborIndex])) == numNeuron) {
							phiSum[index] += CellBoundary(phi_fine[neighborIndex], 0.1);
						}
					}
				}
				phiSum[index] = phiSum[index] > 0 ? CellBoundary(phi_fine[index], 0) / phiSum[index] : 0;
			}
		}
	}

	// Assuming RmOutlier modifies phiSum to remove outliers and returns new maxVal
	float maxVal = RmOutlier(phiSum);

	for (float& value : phiSum) {
		value /= maxVal;
		if (value < threshold) {
			value = 0;
		}
	}
}

// Function to perform Breadth-First Search (BFS) for clustering
void NeuronGrowth::bfs(const vector<float>& matrix, int rows, int cols, int row, int col,
         vector<bool>& visited, vector<pair<int, int>>& cluster) {
	static const int dx[] = {-1, 0, 1, 0, -1, -1, 1, 1}; // Including diagonal directions
	static const int dy[] = {0, -1, 0, 1, -1, 1, -1, 1}; // Including diagonal directions
	// static const int dx[] = {-1, 0, 1, 0}; // Including diagonal directions
	// static const int dy[] = {0, -1, 0, 1}; // Including diagonal directions
	
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
			}
		}
	}
}

// Function to find connected clusters in the matrix
vector<vector<pair<int, int>>> NeuronGrowth::FindClusters(const vector<float>& matrix, int rows, int cols) {
    vector<bool> visited(rows * cols, false);
    vector<vector<pair<int, int>>> clusters;

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			if (matrix[i * cols + j] != 0 && !visited[i * cols + j]) {
				vector<pair<int, int>> cluster;
				bfs(matrix, rows, cols, i, j, visited, cluster);

				if (!cluster.empty()) {
					clusters.push_back(cluster);
				}
			}
		}
	}
	return clusters;
}

// Function to find local maxima within connected clusters in the matrix
vector<float> NeuronGrowth::FindLocalMaximaInClusters(const vector<float>& matrix, int rows, int cols) {
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

// Function to find centroids of connected clusters in the matrix
vector<float> NeuronGrowth::FindCentroidsInClusters(const vector<float>& matrix, int rows, int cols) {
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
		centroidRow = std::max(0, std::min(centroidRow, rows - 1));
		centroidCol = std::max(0, std::min(centroidCol, cols - 1));

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

bool NeuronGrowth::IsLocalMaximum(const vector<float>& matrix, int rows, int cols, int x, int y) {
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

vector<float> NeuronGrowth::FindCentroidsOfLocalMaximaClusters(const vector<float>& matrix, int rows, int cols) {
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

        centroids[centroidRow * cols + centroidCol] = 1.0f;

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

vector<vector<int>> NeuronGrowth::ConvertTo2DIntVector(const vector<float> input, int NX, int NY) 
{
	vector<vector<int>> output;

	int k = 0;
	for (int i = 0; i < NX+1; i++) {
		vector<int> row;
		for (int j = 0; j < NY+1; j++) {
			// row.push_back(CellBoundary(abs(input[k]), 0.01));
			row.push_back(CellBoundary(abs(input[k]), 0.25));
			// row.push_back(CellBoundary(input[k], 0.5));
			k++;
		}
		output.push_back(row);
	}
	return output;
}

vector<vector<float>> NeuronGrowth::ConvertTo2DFloatVector(const vector<float> input, int NX, int NY) 
{
	vector<vector<float>> output;
	int k = 0;
	for (int i = 0; i < NX; i++) {
		vector<float> row;
		for (int j = 0; j < NY; j++) {
			row.push_back(input[k]);
			k++;
		}
		output.push_back(row);
	}
	return output;
}

// Function to perform flood fill
void NeuronGrowth::FloodFill(vector<vector<int>>& image, int x, int y, int newColor, int originalColor) 
{
	// Define the directions: up, down, left, right
	int dx[] = {0, 0, -1, 1};
	int dy[] = {-1, 1, 0, 0};

	if (x < 0 || x >= image.size() || y < 0 || y >= image[0].size() || image[x][y] != originalColor || image[x][y] == newColor) {
		return;
	}

	image[x][y] = newColor;

	// Apply flood fill in all four directions
	for (int i = 0; i < 4; ++i) {
		FloodFill(image, x + dx[i], y + dy[i], newColor, originalColor);
	}
}

// // Function to perform flood fill
// void NeuronGrowth::ExtractPhi(vector<Vertex2D> cpts) 
// {
// 	neurons = ConvertTo2DIntVector(phi, NX, NY);   
// 	for (int i = 0; i < seed.size(); i++) {
// 		int startX = seed[i][0] - originX;
// 		int startY = seed[i][1] - originY;
// 		// std::cout << " startX " << startX << " startY " << startY << std::endl;
// 		int newColor = i+1;
// 		int originalColor = neurons[startX][startY];
// 		FloodFill(neurons, startX, startY, newColor, originalColor);
// 	}
// }

// Function to perform flood fill
void NeuronGrowth::IdentifyNeurons(vector<float>& phi_in, vector<vector<int>>& neurons, const vector<float>& prev_id, vector<array<float, 2>> seed, int NX, int NY, int originX, int originY) 
{
	neurons = ConvertTo2DIntVector(phi_in, NX, NY);   
	// PrintOutNeurons(neurons);
	// std::cout << neurons.size() << " " << NX << " " << originX << " " << NY << " " << originY << std::endl;
	for (int i = 0; i < seed.size(); i++) {
		int startX = (seed[i][0] - originX);
		int startY = (seed[i][1] - originY);
		// std::cout << " startX " << startX << " startY " << startY << std::endl;
		int newColor = i+1;
		int originalColor = neurons[startX][startY];
		FloodFill(neurons, startX, startY, newColor, originalColor);
	}
}

// Check if a point is valid in the grid
bool NeuronGrowth::isValid(int x, int y, int rows, int cols) 
{
    return (x >= 0 && x < rows && y >= 0 && y < cols);
}

// Function to calculate geodesic distance from a point
vector<vector<int>> NeuronGrowth::CalculateGeodesicDistanceFromPoint(vector<vector<int>> neurons, int startX, int startY) 
{
	int rows = neurons.size();
	if (rows == 0) return {};

	int cols = neurons[0].size();

	vector<vector<int>> distances(rows, vector<int>(cols, INF));
	vector<vector<bool>> visited(rows, vector<bool>(cols, false));

	distances[startX][startY] = 0;
	visited[startX][startY] = true;

	queue<pair<int, int>> q;
	q.push({startX, startY});

	int dx[] = {-1, 1, 0, 0};
	int dy[] = {0, 0, -1, 1};

	while (!q.empty()) {
		pair<int, int> current = q.front();
		q.pop();

		int x = current.first;
		int y = current.second;

		for (int i = 0; i < 4; ++i) {
			int newX = x + dx[i];
			int newY = y + dy[i];

			if (isValid(newX, newY, rows, cols) && neurons[newX][newY] == 1 && !visited[newX][newY]) {
				visited[newX][newY] = true;
				distances[newX][newY] = distances[x][y] + 1;
				q.push({newX, newY});
			}
		}
	}
	return distances;
}

// Function to calculate geodesic distance from a point
vector<vector<int>> NeuronGrowth::CalculateGeodesicDistanceFromPoint(vector<vector<int>> neurons, vector<array<float, 2>> &seed, int originX, int originY) 
{
	int rows = neurons.size();
	if (rows == 0) return {};

	int cols = neurons[0].size();

	vector<vector<int>> distances(rows, vector<int>(cols, INF));
	vector<vector<bool>> visited(rows, vector<bool>(cols, false));

	int startX, startY;
	for (int i = 0; i < seed.size(); i++) {
		startX = seed[i][0] - originX;
		startY = seed[i][1] - originY;

		distances[startX][startY] = 0;
		visited[startX][startY] = true;

		queue<pair<int, int>> q;
		q.push({startX, startY});

		int dx[] = {-1, 1, 0, 0};
		int dy[] = {0, 0, -1, 1};

		while (!q.empty()) {
			pair<int, int> current = q.front();
			q.pop();

			int x = current.first;
			int y = current.second;

			for (int i = 0; i < 4; ++i) {
				int newX = x + dx[i];
				int newY = y + dy[i];

				if (isValid(newX, newY, rows, cols) && neurons[newX][newY] == 1 && !visited[newX][newY]) {
					visited[newX][newY] = true;
					distances[newX][newY] = distances[x][y] + 1;
					q.push({newX, newY});
				}
			}
		}
	}
	return distances;
}

// Function to calculate geodesic distance from a point
vector<vector<array<int, 2>>> NeuronGrowth::NeuriteTracing(vector<vector<double>> distance) 
{
	vector<vector<array<int, 2>>> traces;
	return traces;
}

// Function to calculate geodesic distance from a point
void NeuronGrowth::SaveNGvars(vector<vector<float>> &NGvars, int NX, int NY, string fn) 
{
	// writing initialized variables for debugging purposes
	bool visualization = true;
	PrintVec2TXT(NGvars[0], fn + "/phi_" + std::to_string(n) + ".txt", visualization);
	PrintVec2TXT(NGvars[1], fn + "/syn_" + std::to_string(n) + ".txt", visualization);
	PrintVec2TXT(NGvars[2], fn + "/tub_" + std::to_string(n) + ".txt", visualization);
	PrintVec2TXT(NGvars[3], fn + "/theta_" + std::to_string(n) + ".txt", visualization);
	PrintVec2TXT(NGvars[4], fn + "/phi_0_" + std::to_string(n) + ".txt", visualization);
	PrintVec2TXT(NGvars[5], fn + "/tub_0_" + std::to_string(n) + ".txt", visualization);
}

void NeuronGrowth::PrintOutNeurons(vector<vector<int>> neurons) 
{
	ierr = PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------------\n");
	int dwnRatio = 1;
	for (int i = 0; i < neurons.size(); i+=2*dwnRatio) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "| ");
		for (int j = 0; j < neurons[i].size(); j+=dwnRatio) {
			if (round(neurons[i][j]) == 0) {
				ierr = PetscPrintf(PETSC_COMM_WORLD, " ");
			} else {
				ierr = PetscPrintf(PETSC_COMM_WORLD, "#");
			}
		}
		ierr = PetscPrintf(PETSC_COMM_WORLD, "|\n");
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------------\n");
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
	int ind(0); // pre-calculated variable index
	for (int e = 0; e < user->bzmesh_process.size(); e++) {
		int nen = user->bzmesh_process[e].IEN.size(); // 16 supporting cp for 2D case

		vector<vector<float>> EMatrixSolve;	EMatrixSolve.clear();
		vector<float> EVectorSolve;		EVectorSolve.clear();
		EMatrixSolve.resize(nen);
		EVectorSolve.resize(nen);

		user->eleVal.clear();
		user->eleVal.resize(8);
		for (int i = 0; i < nen * 1; i++) {
			// EMatrixSolve[i].reserve(nen*sizeof(float));
			EMatrixSolve[i].resize(nen);
			if (i < 8)
				user->eleVal[i].resize(nen);
		}

		// preparing element stiffness matrix and vector, extract values from control mesh
		for (int i = 0; i < nen; i++) {

			EVectorSolve[i] = 0.0;
			for (int j = 0; j < nen; j++) {
				EMatrixSolve[i][j] = 0.0;
			}

			user->eleVal[0][i] = Parray[user->bzmesh_process[e].IEN[i]];		// elePhiGuess
			user->eleVal[1][i] = user->phi[user->bzmesh_process[e].IEN[i]];		// elePhi
			user->eleVal[4][i] = user->theta[user->bzmesh_process[e].IEN[i]];	// eleTheta
			user->eleVal[6][i] = 0;							// eleEpsilon
			user->eleVal[7][i] = 0;							// eleEpsilonP
		}

		// loop through gaussian quadrature points
		for (int i = 0; i < user->Gpt.size(); i++) {
			for (int j = 0; j < user->Gpt.size(); j++) {
				user->EvaluateOrientation(nen, user->pre_Nx[ind], user->pre_dNdx[ind], user->eleVal[1], user->eleVal[4], user->eleVal[6], user->eleVal[7]);
				user->ElementEvaluationAll_phi(nen, user->pre_Nx[ind], user->pre_dNdx[ind], user->eleVal, user->vars);
				// vars  0   1    2     3      4      5     6     7      8     9      10      11      12     13    14    15      16     17     
				// float C1, C0, elePG, dPGdx, dPGdy, eleP, eleS, eleTb, eleE, eleTp, dThedx, dThedy, eleEP, dAdx, dAdy, eleEEP, dAPdx, dAPdy;
				
				// loop through control points
				for (int m = 0; m < nen; m++) {
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
					for (int n = 0; n < nen; n++) {
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
		for (int i = 0; i < nen; i++) {
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
	for (int e = 0; e < user->bzmesh_process.size(); e++) {
		int nen = user->bzmesh_process[e].IEN.size(); // 16 supporting cp for 2D case

		vector<vector<float>> EMatrixSolve;	EMatrixSolve.clear();
		vector<float> EVectorSolve;		EVectorSolve.clear();
		EMatrixSolve.resize(nen);
		EVectorSolve.resize(nen);

		for (int i = 0; i < nen * 1; i++) {
			// EMatrixSolve[i].reserve(nen*sizeof(float));
			EMatrixSolve[i].resize(nen);
		}
		
		// preparing element stiffness matrix and vector, extract values from control mesh
		vector<float> elePhiGuess(nen, 0);
		for (int i = 0; i < nen; i++) {
			EVectorSolve[i] = 0.0;
			for (int j = 0; j < nen; j++) {
				EMatrixSolve[i][j] = 0.0;
			}
			elePhiGuess[i] = Parray[user->bzmesh_process[e].IEN[i]];		// elePhiGuess
		}

		// loop through gaussian quadrature points
		for (int i = 0; i < user->Gpt.size(); i++) {
			for (int j = 0; j < user->Gpt.size(); j++) {				
				user->ElementEvaluationAll_phi(nen, user->pre_Nx[ind], user->pre_dNdx[ind], elePhiGuess, user->vars);
				
				// loop through control points
				for (int m = 0; m < nen; m++) {
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
					for (int n = 0; n < nen; n++) {
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
		for (int i = 0; i < nen; i++) {
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

	// if (user->n % 50 == 0) {
	// 	PetscInt rows, cols;
	// 	ierr = MatGetSize(user->J, &rows, &cols);
	// 	PetscPrintf(PETSC_COMM_WORLD, "Global matrix size is %D x %D\n", rows, cols);
	// }

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

PetscErrorCode CleanUpSolvers(NeuronGrowth &NG)
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

int RunNG(int& n_bzmesh, vector<vector<int>> ele_process_in, vector<Vertex2D>& cpts_initial, const vector<Element2D>& tmesh_initial, vector<Vertex2D>& cpts_fine, const vector<Element2D>& tmesh_fine, vector<Vertex2D>& cpts, vector<Vertex2D>& prev_cpts, const vector<Element2D>& tmesh, string path_in, string path_out, int& iter, int end_iter_in,
	vector<vector<float>> &NGvars, int &NX, int &NY, vector<array<float, 2>> &seed, int &originX, int &originY, vector<int> &rfid, vector<int> &rftype, bool &localRefine)
{	
	// Declare and initialize Neuron Growth simulation
	NeuronGrowth NG;

	// Set the number of iterations, number of neurons, and the end iteration from provided variables
	NG.n = iter;  // Assign iteration count
	NG.numNeuron = seed.size();  // Set the number of neurons based on the size of 'seed'
	NG.end_iter = end_iter_in;  // Set the ending iteration

	// Initialize vertex clouds for the current, fine, and previous configurations
	Vertex2DCloud cloud(cpts);        // Cloud for current points
	Vertex2DCloud cloud_fine(cpts_fine);  // Cloud for finer resolution points
	Vertex2DCloud cloud_prev(prev_cpts);  // Cloud for previous points

	// Initialize KD-Trees for the current, fine, and previous vertex clouds
	KDTree kdTree(2 /* dim */, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	KDTree kdTree_fine(2 /* dim */, cloud_fine, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	KDTree kdTree_prev(2 /* dim */, cloud_prev, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));

	// Build indexes for KD-Trees to optimize search operations
	kdTree.buildIndex();
	kdTree_fine.buildIndex();
	kdTree_prev.buildIndex();

	// Initialize the Neuron Growth problem with the provided parameters and KDTree for previous points
	NG.InitializeProblemNG(n_bzmesh, cpts, prev_cpts, kdTree_prev, NGvars, seed);

	// Assign processing elements for the simulation
	NG.AssignProcessor(ele_process_in);

	// Synchronize and check MPI element assignments, print out information in order
	for (int i = 0; i < NG.nProcess; i++) {
	if (i == NG.comRank) {
		std::cout << "comRank: " << NG.comRank << "/" << NG.nProcess << " - with element process size: " << NG.ele_process.size() << std::endl;
		NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
	} else {
		NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
	}
	}

	// Read Bezier mesh for the simulation setup
	NG.ReadBezierElementProcess(path_in);
	PetscPrintf(PETSC_COMM_WORLD, "Read bzmesh!-----------------------------------------------------------------\n");

	// (Optional) Set initial guess for SNES based on the current configuration
	// NG.ToPETScVec(NG.phi, NG.temp_phi); // Commented out, indicating it's optional

	PetscPrintf(PETSC_COMM_WORLD, "Set initial guess!-----------------------------------------------------------\n");

	// /*========================================================*/
	// // Initializations
	// NeuronGrowth NG;
	// // NG.SetVariables("simulation_parameters.txt");
	// NG.n = iter;
	// NG.numNeuron = seed.size();
	// NG.end_iter = end_iter_in;

	// Vertex2DCloud cloud(cpts), cloud_fine(cpts_fine), cloud_prev(prev_cpts);
	// KDTree kdTree(2 /* dim */, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	// KDTree kdTree_fine(2 /* dim */, cloud_fine, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	// KDTree kdTree_prev(2 /* dim */, cloud_prev, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	// kdTree.buildIndex();
	// kdTree_fine.buildIndex();
	// kdTree_prev.buildIndex();

	// NG.InitializeProblemNG(n_bzmesh, cpts, prev_cpts, kdTree_prev, NGvars, seed);
	// NG.AssignProcessor(ele_process_in);
	// // Check MPI element assignments, and print out in orders
	// for (int i = 0; i < NG.nProcess; i++){
	// 	if (i == NG.comRank) {
	// 		std::cout << "comRank: " << NG.comRank << "/" << NG.nProcess << " - with element process size: " << NG.ele_process.size() << std::endl;
	// 		NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
	// 	} else {
	// 		NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
	// 	}
	// }				
	// // Read bezier mesh and prepare SNES initial guess
	// NG.ReadBezierElementProcess(path_in);
	// PetscPrintf(PETSC_COMM_WORLD, "Read bzmesh!-----------------------------------------------------------------\n");	
	// // NG.ToPETScVec(NG.phi, NG.temp_phi); // initial guess for SNES (optional)
	// PetscPrintf(PETSC_COMM_WORLD, "Set initial guess!-----------------------------------------------------------\n");	
	
	/*========================================================*/
	// Initial neuron identifications and geodesic distance calculation
	vector<vector<int>> neurons, distances;
	vector<float> phi_fine, id, geodist, tip, localMaximaMatrix, Mphi;
	// phi_fine = InterpolateVars_coarse(NG.phi, cpts, cpts_fine, 1);
	// phi_fine = NG.InterpolateVars_coarse1(NG.phi, cpts, kdTree, cpts_fine, 0, 0);
	// phi_fine = NG.InterpolateValues_closest(NG.phi, cpts, cpts_fine);
	phi_fine = NG.InterpolateValues_closest(NG.phi, kdTree, cpts_fine);
	NG.IdentifyNeurons(phi_fine, neurons, NG.prev_id, seed, NX*2, NY*2, originX, originY);
	distances = NG.CalculateGeodesicDistanceFromPoint(neurons, seed, originX, originY);
	geodist = ConvertTo1DFloatVector(distances);
	PetscPrintf(PETSC_COMM_WORLD, "Calculated geodesic distances!-----------------------------------------------\n");
	id = ConvertTo1DFloatVector(neurons);
	NG.prev_id = id;
	NG.DetectTipsMulti(phi_fine, id, NG.numNeuron, tip, NX*2, NY*2);
	// NG.DetectTipsMulti_new(phi_fine, id, NG.numNeuron, tip, NX*2, NY*2);
	localMaximaMatrix = NG.FindCentroidsOfLocalMaximaClusters(tip, NX*2+1, NY*2+1);
	NG.tips.clear();
	NG.tips.resize(NG.phi.size(),0);
	// NG.tips = InterpolateVars_coarse(localMaximaMatrix, cpts_initial, cpts, 0);
	// NG.tips = NG.InterpolateVars_coarse1(localMaximaMatrix, cpts_fine, kdTree_fine, cpts, 2, 0);
	// NG.tips = NG.InterpolateValues_closest(localMaximaMatrix, cpts_fine, cpts);	
	NG.tips = NG.InterpolateValues_closest(localMaximaMatrix, kdTree_fine, cpts);
	PetscPrintf(PETSC_COMM_WORLD, "Detected initial tips!-------------------------------------------------------\n");

	/*========================================================*/
	// Write initial variables
	string varName;	
	if (NG.n == 0) {
		NG.VisualizeVTK_PhysicalDomain_All(0, path_out);
		PetscPrintf(PETSC_COMM_WORLD, "Saving all variables!--------------------------------------------------------\n");	
	}

	/*========================================================*/
	// calculate some variable in advance to save computational cost
	NG.prepareBasis();
	PetscPrintf(PETSC_COMM_WORLD, "Prepared basis!--------------------------------------------------------------\n");	
	NG.ierr = MPI_Allreduce(&NG.sum_grad_phi0_local, &NG.sum_grad_phi0_global, 1, MPI_FLOAT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
	PetscPrintf(PETSC_COMM_WORLD, "Calculated sum grad phi0!----------------------------------------------------\n");
	NG.prepareTerm_source();
	// NG.prepareEE();
	PetscPrintf(PETSC_COMM_WORLD, "Prepared variables!----------------------------------------------------------\n");
	
	while (iter <= NG.end_iter) {
		NG.n = iter;
		tic();		

		/*========================================================*/
		/*Implcit Non-liear Newton Raphson solver for Phase field equation*/
		NG.phi_prev = NG.phi;	
	
		NG.preparePhaseField();

		float tol = 1e-4;
		float residual = tol + 1;
		int NR_itr = 0;
		while (residual > tol) {
			// PetscPrintf(PETSC_COMM_WORLD, "residual: %f | %f\n", residual, tol);

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
				if (NG.n == 0)
					NG.ierr = KSPView(NG.ksp_phi, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(NG.ierr);
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

			/*========================================================*/
			/*Collecting scattered Phi variable from all processors*/
			Vec temp_phi_seq;
			VecScatter ctx_phi;
			PetscScalar *_p;

			NG.ierr = VecScatterCreateToAll(NG.temp_phi, &ctx_phi, &temp_phi_seq); CHKERRQ(NG.ierr);
			NG.ierr = VecScatterBegin(ctx_phi, NG.temp_phi, temp_phi_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
			NG.ierr = VecScatterEnd(ctx_phi, NG.temp_phi, temp_phi_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
			NG.ierr = VecGetArray(temp_phi_seq, &_p);

			residual = 0;
			for (int i = 0; i < NG.phi.size(); i++) {
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
		t_write += t_phi;
		t_total += t_phi;

		// /*========================================================*/
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
		// // 	std::cout << "Checking snes calls: " << NG.check_itr << std::endl;
		// // }

		// /*========================================================*/
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
		// t_write += t_phi;
		// t_total += t_phi;
		
		/*========================================================*/
		/*Synaptogenesis and Tubulin equation*/
		/*Build Linear System*/
		// NG.prepareBasis();
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
			if (NG.n == 0)
				NG.ierr = KSPView(NG.ksp_syn, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(NG.ierr);
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
			if (NG.n == 0)
				NG.ierr = KSPView(NG.ksp_tub, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(NG.ierr);
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

		/*========================================================*/
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

		for (int i = 0; i < NG.syn.size(); i++) {
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
		t_write += t_synTub;
		t_total += t_synTub;

		/*========================================================*/
		// Obtain initial local refinement information, the very first 25 iterations are purely used 
		// for getting diffused interface for applying local refinements (phi initialization is binary) 
		if ((NG.n == 25) && (localRefine == false)) {
			NGvars.clear(); NGvars.resize(7);
			NGvars[0] = NG.phi;	
			NGvars[1] = NG.syn;
			NGvars[2] = NG.tub;
			NGvars[3] = NG.theta;
			NGvars[4] = NG.phi_0;
			NGvars[5] = NG.tub_0;
			NGvars[6] = NG.prev_id;

			CleanUpSolvers(NG); // Destroy solvers
			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
			iter = 0; // reset to 0 (beginning of the simulation)
			localRefine = true;
			return 2;
		}
		
		/*========================================================*/
		// Neuron identification and tip detection
		// if ((NG.n % 5 == 0) && (NG.n != 0)) {
		
		float t_cc(0);
		if ((NG.n % 1 == 0) || (NG.n == 0)) {	

			PetscPrintf(PETSC_COMM_WORLD, "Step:%d/%d | Phi:[%d]%.3fs | Syn:%d[%d] Tub:%d[%d] %.3fs | Tip:%.3fs | Mesh:%d |\n",\
			NG.n, NG.end_iter, NR_itr, t_phi, reason_syn, its_syn, reason_tub, its_tub, t_synTub, t_tip, NG.phi.size()); CHKERRQ(NG.ierr);

			// phi_fine = InterpolateVars_coarse(NG.phi, cpts, cpts_fine, 0);
			// phi_fine = NG.InterpolateVars_coarse1(NG.phi, cpts, kdTree, cpts_fine, 0, 0);
			// phi_fine = NG.InterpolateValues_closest(NG.phi, cpts, cpts_fine);
			phi_fine = NG.InterpolateValues_closest(NG.phi, kdTree, cpts_fine);
			// PetscPrintf(PETSC_COMM_WORLD, "Got fine |");
			NG.IdentifyNeurons(phi_fine, neurons, NG.prev_id, seed, NX*2, NY*2, originX, originY);
			// PetscPrintf(PETSC_COMM_WORLD, "Got id |");
			distances = NG.CalculateGeodesicDistanceFromPoint(neurons, seed, originX, originY);
			geodist = ConvertTo1DFloatVector(distances);
			// NG.PrintOutNeurons(neurons);
			id = ConvertTo1DFloatVector(neurons);
			// PetscPrintf(PETSC_COMM_WORLD, "Got geo |");
			NG.DetectTipsMulti(phi_fine, id, NG.numNeuron, tip, NX*2, NY*2);
			// PetscPrintf(PETSC_COMM_WORLD, "Got tip |");
			// NG.DetectTipsMulti_new(phi_fine, id, NG.numNeuron, tip, NX*2, NY*2);
			// localMaximaMatrix = NG.FindLocalMaximaInClusters(tip, NX*2+1, NY*2+1);
			// localMaximaMatrix = NG.FindCentroidsInClusters(tip, NX+1, NY+1);
 			localMaximaMatrix = NG.FindCentroidsOfLocalMaximaClusters(tip, NX*2+1, NY*2+1);
			// PetscPrintf(PETSC_COMM_WORLD, "Got loMax 0|");

			int maxGeoInd(0);
			float maxGeo(0);
			for (int i = 0; i < geodist.size(); i++) {
				if ((geodist[i] != INF) && (geodist[i] > maxGeo)) {
					maxGeo = geodist[i];
					maxGeoInd = i;
				}
			}
			if (maxGeoInd > 1) {
				for (int i = -1; i < 1; i++) { // increase Mphi at longest tip 
					for (int j = -1; j < 1; j++) {
							localMaximaMatrix[maxGeoInd+j*(2*NY+1)-i] = 5;
					}
				}
			}

			// localMaximaMatrix[maxGeoInd] = 5;
			// PetscPrintf(PETSC_COMM_WORLD, "Got loMax 1|");

			// NG.tips = InterpolateVars_coarse(localMaximaMatrix, cpts_fine, cpts, 0);
			// NG.tips = NG.InterpolateVars_coarse1(localMaximaMatrix, cpts_fine, kdTree_fine, cpts, 2, 0);
			// NG.tips = NG.InterpolateValues_closest(localMaximaMatrix, cpts_fine, cpts);
			NG.tips = NG.InterpolateValues_closest(localMaximaMatrix, kdTree_fine, cpts);
			// PetscPrintf(PETSC_COMM_WORLD, "Got NGtips|\n");

			toc(t_tip);
			tic();
			t_write += t_tip;
			t_total += t_tip;
		}

		if ((NG.n % NG.var_save_invl == 0) && (NG.n != 0)) {	
			NG.PrintOutNeurons(neurons);
			/*========================================================*/
			// Writing results to files
			PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------\n");
			NG.VisualizeVTK_PhysicalDomain_All(NG.n, path_out);
			
			// if (NG.comRank == 0) {
			// 	string outputPath = path_out + "/phi_" + to_string(NG.n) + ".png";
			// 	WriteMatrixToPNG(NG.phi, NX+1, NY+1, outputPath.c_str());
			// }
			varName = "phifine_running_";
			NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, NG.n, path_out, phi_fine, varName); // solution on control points
			
			varName = "geoDist_running_";
			NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, NG.n, path_out, geodist, varName); // solution on control points
			varName = "localMax_running_";
			NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, NG.n, path_out, localMaximaMatrix, varName); // solution on control points
			varName = "tip_running_";
			NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, NG.n, path_out, tip, varName); // solution on control points
			varName = "id_running_";
			NG.VisualizeVTK_ControlMesh(cpts_fine, tmesh_fine, NG.n, path_out, id, varName); // solution on control points
			// varName = "phi_running_";
			// NG.VisualizeVTK_ControlMesh(cpts, tmesh, NG.n, path_out, NG.phi, varName); // solution on control points
			// varName = "Mphi_running_";
			// NG.VisualizeVTK_ControlMesh(cpts_initial, tmesh_initial, NG.n, path_out, Mphi, varName); // solution on control points
			
			PetscPrintf(PETSC_COMM_WORLD, "| Wrote Physical Domain! | Average time %fs | Total time: %f |\n", 
				t_total/NG.n, t_total); CHKERRQ(NG.ierr);
				// t_write/NG.var_save_invl, t_total); CHKERRQ(NG.ierr);
			PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------\n");
			t_write = 0;
		}

		/*========================================================*/
		// Domain expansion and variable passing - back to main.cpp
		if ((NG.n % NG.expandCK_invl == 0) && (NG.n != 0)) {
			// NG.PrintOutNeurons(neurons); 
			// int expd_dir_local = NG.CheckExpansion(id, NX, NY);
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

			NGvars.clear(); NGvars.resize(6);
			NGvars[0] = NG.phi;	
			NGvars[1] = NG.syn;
			NGvars[2] = NG.tub;
			NGvars[3] = NG.theta;
			NGvars[4] = NG.phi_0;
			NGvars[5] = NG.tub_0;

			CleanUpSolvers(NG); // Destroy solvers
			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
			iter += 1;
			return 2;

		}
		iter += 1;		 
	}

	/*========================================================*/
	CleanUpSolvers(NG); // Destroy solvers
	NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);

	return 0; // ending simulation
}