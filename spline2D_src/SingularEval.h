#ifndef SINGULAREVAL_H
#define SINGULAREVAL_H

#include <vector>
#include <Eigen/Dense>
#include "BasicDataStructure.h"
using std::vector;
using namespace Eigen;

class EigenStruct
{
public:
	vector<double> L;//eigenvalues
	vector<double> iV;//inverse of eigenvectors
	vector<vector<double>> x;//coeffs of splines
};

class SingularPatchEval
{
public:
	SingularPatchEval();
	void ReadEigenStruct();
	void Print();
	void ProjectPoints(vector<Vertex>& Cp,vector<Vertex>& C,int N);
	void ObtainIVcoef(int N, vector<double>& vec);
	void ObtainIVMat(int N, vector<vector<double>>& mat);
	void EvalSpline(const vector<double>& coef,double u,double v,int i,int N,double& val,double Dval[2]);
	//void Evaluate(Vertex& P,double u,double v,vector<Vertex>& C,int N);
	double ExplicitBasis(int bsid,double u,double v,int N);
	//void SingularBasis(double u,double v,int N,VectorXd& Nt,MatrixXd& dNdt);

	//void EvalSpline2(const vector<double>& coef,double u,double v,int i,int N,double& val,double Dval[2],double D2val[3]);
	//void SingularBasis2(double u,double v,int N,VectorXd& Nt,MatrixXd& dNdt,MatrixXd& d2Ndt);
private:
	int Nmax;
	vector<EigenStruct> eigen;
};

#endif