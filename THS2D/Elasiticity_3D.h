#ifndef ELASITICITY_3D_H
#define ELASITICITY_3D_H

#include <vector>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "BasicDataStructure.h"
//#include "Truncated_Tspline.h"

using namespace std;
using namespace Eigen;

class LinearElasticity
{
public:
	LinearElasticity();
	void SetMaterialProp(double Ey=1., double v=0.3);
	void GaussInfo(int ng=3);
	void SetProblem(const vector<int>& IDBC_in, const vector<double>& gh_in);
	void SetBoundary();
	void BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, double Nx[64], double dNdx[64][3], double& detJ);
	void ElementMatrix(double dNdx[64][3], double detJ, double EK[192][192]);
	void ElementForce(double Nx[16], double Jmod, double Fb[2], double EF[32]);
	void Assembly(double EK[192][192], double EF[192], const vector<vector<double>>& cmat, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	void BuildLinearSystem(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void Solver(SparseMatrix<double>& GK, VectorXd& GF);
	void Para2Phys(double u, double v, const double cpt[16][3], double pt[3]);
	void DispStrainStress(double u,double v,double w, const BezierElement3D& bzel,double disp[3],double se[6],double ss[6]);
	//void StressCalculate(double u,double v,const BezierElement& bzel,double stress[6]);
	void VisualizeVTK(const vector<BezierElement3D>& bzmesh, string fn);
	//void VisualizeError(string fn);
	//void ElementErrorEstimate(int eid,double& err,double& area,vector<int>& pid);
	//void TotalErrorEstimate();
	void Run(const vector<BezierElement3D>& bzmesh, string fn);

	void BasisFunction4(double u, double v, const double pt[25][3], double Nx[25], double dNdx[25][2], double& detJ);
	void ElementMatrix4(double dNdx[25][2], double detJ, double EK[50][50]);
	void ElementForce4(double Nx[25], double Jmod, double Fb[2], double EF[50]);
	void Assembly4(double EK[50][50], double EF[50], const vector<array<double,25>>& cmat, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	void Para2Phys4(double u, double v, const double cpt[25][3], double pt[3]);
	void DispStrainStress4(double u,double v,const BezierElement3D& bzel,double disp[2],double se[3],double ss[3]);

	void BasisFunction_TSP(double u, double v, const BezierElement3D& bzel, vector<double>& Nx, vector<array<double,2>>& dNdx, double& detJ);
	void ElementMatrix_TSP(const vector<array<double,2>>& dNdx, double detJ, vector<vector<double>>& EK);
	void ElementForce_TSP(const BezierElement3D& bzel, vector<double>& EF);
	void Assembly_TSP(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	void AssemblyForce_TSP(const vector<double>& EF, const vector<int>& IEN, VectorXd& GF);
	void BuildLinearSystem_TSP(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void DispStrainStress_TSP(double u,double v,const BezierElement3D& bzel,double disp[2],double se[3],double ss[3]);
	void Para2Phys_TSP(double u, double v, const BezierElement3D& bzel, array<double,3>& pt);

	//void SetProblem_PlateHole(const TruncatedTspline& tts);
	void Quantityh_TSP(double u,double v,const BezierElement3D& bzel,double disp[2],double ss[3], double& detJ);
	void ElementError(const BezierElement3D& bzel, double& L2);
	double OverallError(const vector<BezierElement3D>& bzmesh);

	void BasisFunction_Rational(double u, double v, const BezierElement3D& bzel, vector<double>& Rx, vector<array<double,2>>& dRdx, double& detJ);
	void BasisFunction4_Rational(double u, double v, const BezierElement3D& bzel, vector<double>& Rx, vector<array<double,2>>& dRdx, double& detJ);

	void Neumann_Element(double u, double w, int ids[4], int edid, const BezierElement3D& bzel, vector<double>& EF);

	void OutputMatrix(SparseMatrix<double>& GK, string fn);
	
private:
	vector<double> Gpt;
	vector<double> wght;
	double lambda;//Lame parameters
	double mu;
	vector<int> IDBC;
	vector<int> BCList;
	vector<int> BC_Neumann;
	vector<double> gh;
	vector<double> uh;
	int npt;
	int neq;
	int dim;
};

void ExactDisp_PlateHole(double x, double y, double u[2]);
void ExactStress_PlateHole(double x, double y, double ss[3]);

#endif
