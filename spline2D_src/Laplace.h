#ifndef LAPLACE_H
#define LAPLACE_H

#include <vector>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "BasicDataStructure.h"
#include "Truncated_Tspline.h"

using namespace std;
using namespace Eigen;

class Laplace
{
public:
	Laplace();
	void GaussInfo(int ng=3);
	void SetProblem(const TruncatedTspline& tts);
	void SetBoundary();
	void BasisFunction(double u, double v, const double pt[16][3], double Nx[16], double dNdx[16][2], double& detJ);
	void BasisFunction4(double u, double v, const double pt[25][3], double Nx[25], double dNdx[25][2], double& detJ);
	void ElementMatrix(double dNdx[16][2], double detJ, double EK[16][16]);
	void ElementForce(double Nx[16], double Jmod, double Fb[2], double EF[32]);
	void Assembly(double EK[16][16], double EF[16], const vector<array<double,16>>& cmat, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	void BuildLinearSystem(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void Solver(SparseMatrix<double>& GK, VectorXd& GF);
	void Para2Phys(double u, double v, const double cpt[25][3], double pt[3]);
	void DispCal(double u,double v,const BezierElement& bzel,double& disp,double& detJ);
	void ElementError(const BezierElement& bzel, double& L2, double& H1);
	void VisualizeVTK(string fn,const vector<BezierElement>& bzmesh);

	void DispCal_Coupling_Bezier(double u, double v, double w, const BezierElement& bzel, double pt[3], double& disp, double dudx[3], double& detJ);
	void VisualizeVTK3D(string fn, const vector<BezierElement>& bzmesh);
	void VisualizeError(string fn);
	//void ElementErrorEstimate(int eid,double& err,double& area,vector<int>& pid);
	//void TotalErrorEstimate();
	void Run(const vector<BezierElement>& bzmesh, string fn);

	void SetProblem(const vector<int>& IDBC_in, const vector<double>& gh_in);

	void BasisFunction_TSP(double u, double v, const BezierElement& bzel, vector<double>& Nx, vector<array<double,2>>& dNdx, double& detJ);
	void ElementMatrix_TSP(const vector<array<double,2>>& dNdx, double detJ, vector<vector<double>>& EK);
	void ElementForce_TSP(const vector<double>& Nx, double detJ, double Fb, vector<double>& EF);
	void Assembly_TSP(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	void BuildLinearSystem_TSP(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void Para2Phys_TSP(double u, double v, const BezierElement& bzel, array<double,3>& pt);
	void DispCal_TSP(double u,double v,const BezierElement& bzel,double& disp,double& detJ);
	void Quantityh_TSP(double u,double v,const BezierElement& bzel,double val[3],double& detJ);

	void setProblem_ptestHO(const TruncatedTspline& tts);

	void BasisFunction_Dual(double u, double v, const BezierElement& bzel, vector<double>& Nx, vector<array<double, 2>>& dNdx, double& detJ);
	void InitializeSparseMat(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK);
	void BuildLinearSystem_Dual(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void BasisFunctionDual_Irr(double u, double v, const BezierElement& bzel, vector<double>& Nx, vector<array<double, 2>>& dNdx, double& detJ);
	void Para2Phys_Dual(double u, double v, const BezierElement& bzel, const vector<double>& Nx, array<double, 3>& pt);
	void DispCal_Dual(double u, double v, const BezierElement& bzel, double& disp, double& detJ);
	void Quantityh_Dual(const vector<int>& IEN, const vector<double>& Nx, const vector<array<double, 2>>& dNdx, double val[3]);
	void SetProblem_Dual(const vector<int>& IDBC_in, const vector<double>& gh_in);
	void ElementErrorDual(const BezierElement& bzel, double& L2, double& H1);
	void VisualizeVTK_Dual(string fn, const vector<BezierElement>& bzmesh);
	void Run_Dual(const vector<BezierElement>& bzmesh, string fn);
	void VisualizeBasisFunction(int pid, const vector<BezierElement>& bzmesh, string fn);

	void BasisFunction_CC(double u, double v, const BezierElement& bzel, vector<double>& Nx, vector<array<double, 2>>& dNdx, double& detJ);
	void BuildLinearSystem_CC(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void BasisFunctionCC_Irr(double u, double v, const BezierElement& bzel, vector<double>& Nx, vector<array<double, 2>>& dNdx, double& detJ);
	void ElementErrorCC(const BezierElement& bzel, double& L2, double& H1);
	void VisualizeVTK_CC(string fn, const vector<BezierElement>& bzmesh);
	void Run_CC(const vector<BezierElement>& bzmesh, string fn);
	
private:
	vector<double> Gpt;
	vector<double> wght;
	vector<int> IDBC;
	vector<int> BCList;
	vector<double> gh;
	vector<double> uh;
	int npt;
	int neq;
};

double LDomainSolution(double x,double y);
void LDomainSolutionDeriv(double x, double y, double& u, double& ux, double& uy);

double ptestHO_sol_1(double x, double y);
double ptestHO_sol_2(double x, double y);
double ptestHO_sol_3(double x, double y);

double exact_sol(double x, double y);
double exact_sol_1(double x, double y);
double exact_sol_4(double x, double y);
double exact_sol_5(double x, double y);
double exact_sol_6(double x, double y);
double exact_sol_7(double x, double y);

void grad_u(double x, double y, double gu[2]);
void grad_u_1(double x, double y, double gu[2]);
void grad_u_4(double x, double y, double gu[2]);
void grad_u_6(double x, double y, double gu[2]);
void grad_u_7(double x, double y, double gu[2]);

double f_source(double x, double y);
double f_source_1(double x, double y);//u=x+y, f=0
double f_source_2(double x, double y);//u=x^2+y^2, f=-4
double f_source_3(double x, double y);//u=x^3+y^3, f=-6(x+y)
double f_source_4(double x, double y);//u=x^3+y^3, f=-6(x+y)
double f_source_5(double x, double y);//u=x, f=0
double f_source_6(double x, double y);//u=sin(PI x)sin(PI y)
double f_source_7(double x, double y);//u=exp((x+y)/2), f=0

#endif