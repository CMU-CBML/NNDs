#ifndef KERNEL_H
#define KERNEL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include "KnotInsertion.h"
//#include "TTSP_2D.h"
#include "TTSP_3D.h"
// #include "Elasiticity_3D.h"
#include "Laplace.h"
#include "LeastSquare.h"
#include <sstream>
#include <iomanip>
using namespace std;

class kernel
{
public:
	// void run();

	// void run_complex();
	void OutputMesh(const vector<BezierElement3D>& bzmesh, string fn);
	void run_complex_fit();
	// void run_complex_glb();
	// void run_Bspline();
	// void run_Bspline_fit();

	// void run_benchmark();
	// void run_leastsquare();
	// void run_leastsquare_1();
	// void run_pipeline();
	void output_err(string fn, const vector<int>& dof, const vector<double>& err);
	void output_err(string fn, const vector<double>& dof, const vector<double>& err);
//private:
//	TruncatedTspline_3D tt3;
//	LinearElasticity le;
//	Laplace lap;
};

#endif