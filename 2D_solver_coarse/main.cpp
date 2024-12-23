#include <iostream>
#include <vector>
#include <array>
#include "NeuronGrowth.h"
#include <sstream>
#include <iomanip>
#include "time.h"
#include "utils.h"

using namespace std;

static char help[] = "Solve 2D Neuron Growth\n";

int main(int argc, char **argv)
{
	int rank, nProcs;
	PetscErrorCode ierr;
	/// start up petsc
	ierr = PetscInitialize(&argc, &argv, (char*)0, help); if (ierr) return ierr;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &nProcs); CHKERRQ(ierr);

	int numNeuron = atoi(argv[1]);		// user input number of neurons
	int end_iter = atoi(argv[2]);		// user input number of iterations
	string path_in = argv[3];		// user input working directory

	int NX, NY, originX(0), originY(0);
	vector<array<float, 2>> seed; 
	if (numNeuron <= 2) {
		InitializeSoma(numNeuron, seed, NX, NY); // predefined soma placements
	} else {
		if (rank == 0) {
			InitializeRandomSoma(numNeuron, seed, NX, NY); // randomized soma placements
			// InitializeSoma_customizedCases(numNeuron, seed, NX, NY); // randomized soma placements
		}
		ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);
		// Broadcast NX, NY, and seed from the root process to all other processes
		MPI_Bcast(&NX, 1, MPI_INT, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&NY, 1, MPI_INT, 0, PETSC_COMM_WORLD);
		int size = seed.size();
		// First broadcast the size of the vector
		MPI_Bcast(&size, 1, MPI_INT, 0, PETSC_COMM_WORLD);
		// Resize the vector on non-root processes
		if (rank != 0) {
			seed.resize(size);
		}
		// Broadcast the actual data
		MPI_Bcast(seed.data(), size * 2, MPI_FLOAT, 0, PETSC_COMM_WORLD);
	}
	
	/// Set simulation parameters and mesh
	string fn_mesh_initial(path_in + "controlmesh_initial.vtk");
	string fn_mesh(path_in + "controlmesh.vtk");
	string fn_mesh_fine(path_in + "controlmesh_fine.vtk");
	string fn_bz(path_in + "bzmeshinfo.txt.epart." + to_string(nProcs));
	string path_out(path_in + "outputs/");
	// string path_out(path_in);

	int n_bzmesh;
	vector<vector<float>> vertices;
	vector<vector<int>> elements;
	vector<vector<int>> ele_process;
	vector<Vertex2D> cpts_initial, cpts, cpts_fine, prev_cpts, prev_cpts_fine;
	vector<Element2D> tmesh_initial, tmesh, tmesh_fine;
	ele_process.resize(nProcs);
	PetscPrintf(PETSC_COMM_WORLD, "Initial element process size: %ld\n", ele_process.size()); CHKERRQ(ierr);

	vector<int> rfid, rftype; // list of elements to refine
	vector<vector<float>> NGvars; // for passing neuron growth variables in-between domain expansion
	NGvars.clear(); NGvars.resize(7);

	bool localRefine = false;
	int iter(0), state(1); // 0-end, 1-running, 2-expanding, 3-diverging
	while (iter <= end_iter) {
		prev_cpts = cpts; // back up old control points for later NGvars interpolations (old cpts to new cpts)
		prev_cpts_fine = cpts_fine;
		cpts_initial.clear(); tmesh_initial.clear(); 
		cpts.clear(); tmesh.clear(); 
		ele_process.clear();
		ele_process.resize(nProcs);

		if (rank == 0) {
			// to make sure correct files are generated and then read later on
			std::remove("../io/controlmesh.vtk");
			std::remove("../io/controlmesh_initial.vtk");
			std::remove("../io/controlmesh_fine.vtk");
			std::remove("../io/bzpt.txt");
			std::remove("../io/cmat.txt");
			std::remove("../io/bzmesh.vtk");
			std::remove("../io/bzmeshinfo.txt");
			string epart("../io/bzmeshinfo.txt.epart." + std::to_string(nProcs));
			string npart("../io/bzmeshinfo.txt.npart." + std::to_string(nProcs));
			std::remove(epart.c_str());
			std::remove(npart.c_str());

			gen2Dmesh(originX, originY, NX, NY, vertices, elements, 1); // Generating 2D quad mesh
			write_quad_toVTK(fn_mesh_initial.c_str(), vertices, elements);

			gen2Dmesh(originX, originY, NX*2, NY*2, vertices, elements, 2); // Generating 2D quad mesh
			write_quad_toVTK(fn_mesh_fine.c_str(), vertices, elements);
			
			// localRefine = false;
			if (localRefine != true) {
				write_quad_toVTK(fn_mesh.c_str(), vertices, elements);
				bzmesh2D(path_in); // generating 2D bezier mesh information (bzmeshinfo.txt)
			} else {
				ReadMesh(fn_mesh_initial, cpts_initial, tmesh_initial);
				vector<float> tmp = InterpolateVars_coarse(NGvars[0], prev_cpts, cpts_initial, 1);
				ObtainRefineID_coarse(tmp, cpts_initial, NX, NY, originX/2, originY/2, rfid, rftype);
				THS2D(path_in, rfid, rftype);
			}
			mpmetis(nProcs, path_in); // partitioning bzmeshinfo using mpmetis for parallelization
		}
		ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

		ReadMesh(fn_mesh_initial, cpts_initial, tmesh_initial);
		ReadMesh(fn_mesh, cpts, tmesh);
		ReadMesh(fn_mesh_fine, cpts_fine, tmesh_fine);
		AssignProcessor(fn_bz, n_bzmesh, ele_process);
		PetscPrintf(PETSC_COMM_WORLD, "Processor Assigned!----------------------------------------------------------\n");
	
		state = RunNG(n_bzmesh, ele_process, cpts_initial, tmesh_initial, cpts_fine, tmesh_fine, cpts, prev_cpts, prev_cpts_fine, tmesh, path_in, path_out,
			iter, end_iter, NGvars, NX, NY, seed, originX, originY, rfid, rftype, localRefine);

		if (state == 3) {
			PetscPrintf(PETSC_COMM_WORLD, "Simulation diverged!---------------------------------------------------------\n");
			return 0;
		}
	}

	PetscPrintf(PETSC_COMM_WORLD, "Done - Main reached end iteration!\n");
	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;
}
