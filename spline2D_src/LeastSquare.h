#ifndef LEASTSQUARE_H
#define LEASTSQUARE_H

#include <vector>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "BasicDataStructure.h"
//#include "Truncated_Tspline.h"

using namespace std;
using namespace Eigen;

class LeastSquare
{
public:
	LeastSquare();
	void SetProblem(int npt_in);
	void SetProblem(const vector<int>& IDBC_in, const vector<double>& gh_in);
	void SetSamplingPoints(int n=5);
	void GaussInfo(int ng=3);
	void BasisFunction(double u, double v, const double pt[16][3], double Nx[16], double dNdx[16][2], double& detJ);
	void BasisFunction4(double u, double v, const double pt[25][3], double Nx[25], double dNdx[25][2], double& detJ);
	void Solver(SparseMatrix<double>& GK, vector<VectorXd>& GF, vector<array<double, 3>>& cp_out);
	void Solver(SparseMatrix<double>& GK, VectorXd& GF, vector<double>& cp);
	void VisualizeVTK(string fn,const vector<BezierElement>& bzmesh);
	void Run(const vector<BezierElement>& bzmesh, string fn, vector<array<double, 3>>& cp_out);
	void Run_ScalarFitting(const vector<BezierElement>& bzmesh, string fn, vector<double>& gh_out);

	void BasisFunction_TSP(double u, double v, const BezierElement& bzel, vector<double>& Nx, vector<array<double,2>>& dNdx, double& detJ);
	void BasisFunction_TSP(double u, double v, const BezierElement& bzel,
		vector<double>& Nx, vector<array<double, 2>>& dNdx, double& detJ, double pt[3]);
	void ElementMatrix_TSP(const vector<double>& dNdx, vector<vector<double>>& EK);
	void ElementForce_TSP(const vector<double>& Nx, double Fb, vector<double>& EF);
	void ElementForce_TSP(const vector<double>& Nx, double Fb[3], vector<vector<double>>& EF);
	void GetFb(double u, double v, const BezierElement& bzel, double Fb[3]);
	void Assembly_TSP(const vector<vector<double>>& EK, const vector<vector<double>>& EF, const vector<int>& IEN, 
		SparseMatrix<double>& GK, vector<VectorXd>& GF);
	void Assembly_TSP(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN,
		SparseMatrix<double>& GK, VectorXd& GF);
	void InitializeSparseMat(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK);
	void BuildLinearSystem_TSP(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK, vector<VectorXd>& GF);
	void BuildLinearSystem_TSP(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void Para2Phys_TSP(double u, double v, const BezierElement& bzel, array<double,3>& pt);
	void DispCal_TSP(double u,double v,const BezierElement& bzel,double& disp,double& detJ);
	void Quantityh_TSP(double u,double v,const BezierElement& bzel,double val[3],double& detJ);

	double ExactValue(double x, double y);
	double ExactValue_1(double x, double y);
	double ExactValue_6(double x, double y);
	double ExactValue_7(double x, double y);

	
private:
	vector<double> Gpt;
	vector<double> wght;
	vector<int> IDBC;
	vector<int> BCList;
	vector<double> gh;
	vector<double> uh;
	int npt;
	int neq;

	//vector<array<double, 3>> cp;
};

#endif