#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cstdlib>

// #include "Elasiticity_3D.h"
#include "KnotInsertion.h"
#include "TTSP_2D.h"
#include "TTSP_3D.h"
#include "Laplace.h"
#include "kernel.h"

using namespace std;

static char help[] = "THS\n";

int main(int argc, char **argv)
{

	string path_in = argv[1];

	vector<int> ids;
	ids.resize(argc-2);
	for (int i = 0; i<argc-2; i++) {
		ids[i] = atoi(argv[i+2]);
		// std::cout << ids[i] << " ";
	}
	// std::cout << std::endl;
	// std::cout << "==============="<< std::endl;

	// Generating locally refined 2D control mesh
	TruncatedTspline_2D tt2;
	// tt2.runXP_complex("../iotest/controlmesh_initial");
	tt2.run_surf_XP(path_in, ids);
	// tt2.run_surf_XP("../io/surf/controlmesh");
	// tt2.run_surf_XP("./io/surf/2D_quad");
	// std::cout << "==============="<< std::endl;

	// // 3D version
	// kernel app;
	// app.run_complex_fit();
	// // // tt3.SetProblem("./io/hex_input/cube5");

	return 0;
}