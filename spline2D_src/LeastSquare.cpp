#include "LeastSquare.h"
#include "BSplineBasis.h"
#include <cmath>
#include <fstream>
#include <iostream>

#define PI 3.141592654

using namespace std;

typedef unsigned int uint;

LeastSquare::LeastSquare()
{
}

void LeastSquare::SetProblem(int npt_in)
{
	//cp.resize(cp_in.size());
	//for (uint i = 0; i < cp_in.size(); i++)
	//{
	//	cp[i][0] = cp_in[i].coor[0];
	//	cp[i][1] = cp_in[i].coor[1];
	//	cp[i][2] = cp_in[i].coor[2];
	//}
	neq = npt_in;
	npt = npt_in;
	IDBC.resize(npt_in);
	for (int i = 0; i < npt_in; i++)
	{
		IDBC[i] = i;
	}
}

void LeastSquare::SetProblem(const vector<int>& IDBC_in, const vector<double>& gh_in)
{
	neq = 0;
	npt = IDBC_in.size();
	IDBC = IDBC_in;
	gh = gh_in;
	for (uint i = 0; i < IDBC.size(); i++)
	{
		if (IDBC[i] != -1)
		{
			neq++;
		}
	}
}

void LeastSquare::SetSamplingPoints(int n)
{
	Gpt.clear();
	Gpt.resize(n);
	for (int i = 0; i < n; i++)
	{
		Gpt[i] = double(i) / double(n - 1);
	}
}

void LeastSquare::GaussInfo(int ng)
{
	Gpt.clear();
	wght.clear();
	switch(ng)
	{
	case 2:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.2113248654051871;			Gpt[1]=0.7886751345948129;
			wght[0]=1.;			wght[1]=1.;
			break;
		}
	case 3:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.1127016653792583;			Gpt[1]=0.5;			Gpt[2]=0.8872983346207417;
			wght[0]=0.5555555555555556;			wght[1]=0.8888888888888889;			wght[2]=0.5555555555555556;
			break;
		}
	case 4:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.06943184420297371;			Gpt[1]=0.33000947820757187;			Gpt[2]=0.6699905217924281;			Gpt[3]=0.9305681557970262;
			wght[0]=0.3478548451374539;			wght[1]=0.6521451548625461;			wght[2]=0.6521451548625461;			wght[3]=0.3478548451374539;
			break;
		}
	default:
		{
			Gpt.resize(2);
			wght.resize(2);
			Gpt[0]=0.2113248654051871;			Gpt[1]=0.7886751345948129;
			wght[0]=1.;			wght[1]=1.;
			break;
		}
	}
}

void LeastSquare::BasisFunction(double u, double v, const double pt[16][3], double Nx[16], double dNdx[16][2], double& detJ)
{
	double Nu[4]={(1.-u)*(1.-u)*(1.-u),3.*(1.-u)*(1.-u)*u,3.*(1.-u)*u*u,u*u*u};
	double Nv[4]={(1.-v)*(1.-v)*(1.-v),3.*(1.-v)*(1.-v)*v,3.*(1.-v)*v*v,v*v*v};
	double dNdu[4]={-3.*(1.-u)*(1.-u),3.-12.*u+9.*u*u,3.*(2.-3.*u)*u,3.*u*u};
	double dNdv[4]={-3.*(1.-v)*(1.-v),3.-12.*v+9.*v*v,3.*(2.-3.*v)*v,3.*v*v};
	double dNdt[16][2];
	int i,j,a,b,loc(0);
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			Nx[loc]=Nu[j]*Nv[i];
			dNdt[loc][0]=dNdu[j]*Nv[i];
			dNdt[loc][1]=Nu[j]*dNdv[i];
			loc++;
		}
	}
	double dxdt[2][2]={{0.,0.},{0.,0.}};
	loc=0;
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			for (a=0;a<2;a++)
			{
				for(b=0;b<2;b++)
				{
					dxdt[a][b]+=pt[loc][a]*dNdt[loc][b];
				}
			}
			loc++;
		}
	}
	detJ=dxdt[0][0]*dxdt[1][1]-dxdt[0][1]*dxdt[1][0];
	double dtdx[2][2]={{dxdt[1][1]/detJ,-dxdt[0][1]/detJ},{-dxdt[1][0]/detJ,dxdt[0][0]/detJ}};
	for(i=0;i<16;i++)
	{
		dNdx[i][0]=dNdt[i][0]*dtdx[0][0]+dNdt[i][1]*dtdx[1][0];
		dNdx[i][1]=dNdt[i][0]*dtdx[0][1]+dNdt[i][1]*dtdx[1][1];
	}
	detJ=0.25*detJ;
}

void LeastSquare::BasisFunction4(double u, double v, const double pt[25][3], double Nx[25], double dNdx[25][2], double& detJ)
{
	double Nu[5]={(1.-u)*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-u)*u,6.*(1.-u)*(1.-u)*u*u,4.*(1.-u)*u*u*u,u*u*u*u};
	double Nv[5]={(1.-v)*(1.-v)*(1.-v)*(1.-v),4.*(1.-v)*(1.-v)*(1.-v)*v,6.*(1.-v)*(1.-v)*v*v,4.*(1.-v)*v*v*v,v*v*v*v};
	double dNdu[5]={-4.*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-4.*u),12.*u*(1.-3.*u+2.*u*u),4.*(3.-4.*u)*u*u,4.*u*u*u};
	double dNdv[5]={-4.*(1.-v)*(1.-v)*(1.-v),4.*(1.-v)*(1.-v)*(1.-4.*v),12.*v*(1.-3.*v+2.*v*v),4.*(3.-4.*v)*v*v,4.*v*v*v};
	double dNdt[25][2];
	int i,j,a,b,loc(0);
	for(i=0;i<5;i++)
	{
		for(j=0;j<5;j++)
		{
			Nx[loc]=Nu[j]*Nv[i];
			dNdt[loc][0]=dNdu[j]*Nv[i];
			dNdt[loc][1]=Nu[j]*dNdv[i];
			loc++;
		}
	}
	double dxdt[2][2]={{0.,0.},{0.,0.}};
	loc=0;
	for(i=0;i<5;i++)
	{
		for(j=0;j<5;j++)
		{
			for (a=0;a<2;a++)
			{
				for(b=0;b<2;b++)
				{
					dxdt[a][b]+=pt[loc][a]*dNdt[loc][b];
				}
			}
			loc++;
		}
	}
	detJ=dxdt[0][0]*dxdt[1][1]-dxdt[0][1]*dxdt[1][0];
	double dtdx[2][2]={{dxdt[1][1]/detJ,-dxdt[0][1]/detJ},{-dxdt[1][0]/detJ,dxdt[0][0]/detJ}};
	for(i=0;i<25;i++)
	{
		dNdx[i][0]=dNdt[i][0]*dtdx[0][0]+dNdt[i][1]*dtdx[1][0];
		dNdx[i][1]=dNdt[i][0]*dtdx[0][1]+dNdt[i][1]*dtdx[1][1];
	}
	detJ=0.25*detJ;
}

void LeastSquare::BasisFunction_TSP(double u, double v, const BezierElement& bzel, vector<double>& Nx, vector<array<double,2>>& dNdx, double& detJ)
{
	Nx.clear();
	dNdx.clear();
	unsigned int i,j;
	if(bzel.order==3)
	{
		double Nx0[16];
		double dNdx0[16][2];
		BasisFunction(u,v,bzel.pts,Nx0,dNdx0,detJ);
		Nx.resize(bzel.cmat.size(),0);
		dNdx.resize(bzel.cmat.size());
		for(i=0; i<bzel.cmat.size(); i++)
		{
			dNdx[i][0]=0.; dNdx[i][1]=0.;
			for(j=0; j<bzel.cmat[i].size(); j++)
			{
				Nx[i]+=bzel.cmat[i][j]*Nx0[j];
				dNdx[i][0]+=bzel.cmat[i][j]*dNdx0[j][0];
				dNdx[i][1]+=bzel.cmat[i][j]*dNdx0[j][1];
			}
		}
	}
	else if(bzel.order==4)
	{
		double Nx0[25];
		double dNdx0[25][2];
		BasisFunction4(u,v,bzel.pts4,Nx0,dNdx0,detJ);
		Nx.resize(bzel.cmat4.size(),0);
		dNdx.resize(bzel.cmat4.size());
		for(i=0; i<bzel.cmat4.size(); i++)
		{
			dNdx[i][0]=0.; dNdx[i][1]=0.;
			for(j=0; j<bzel.cmat4[i].size(); j++)
			{
				Nx[i]+=bzel.cmat4[i][j]*Nx0[j];
				dNdx[i][0]+=bzel.cmat4[i][j]*dNdx0[j][0];
				dNdx[i][1]+=bzel.cmat4[i][j]*dNdx0[j][1];
			}
		}
	}
}

void LeastSquare::BasisFunction_TSP(double u, double v, const BezierElement& bzel, 
	vector<double>& Nx, vector<array<double, 2>>& dNdx, double& detJ, double pt[3])
{
	Nx.clear();
	dNdx.clear();
	pt[0] = 0.; pt[1] = 0.; pt[2] = 0.;
	unsigned int i, j;
	if (bzel.order == 3)
	{
		double Nx0[16];
		double dNdx0[16][2];
		BasisFunction(u, v, bzel.pts, Nx0, dNdx0, detJ);
		for (i = 0; i < 16; i++)
		{
			pt[0] += Nx0[i] * bzel.pts[i][0];
			pt[1] += Nx0[i] * bzel.pts[i][1];
			pt[2] += Nx0[i] * bzel.pts[i][2];
		}
		Nx.resize(bzel.cmat.size(), 0);
		dNdx.resize(bzel.cmat.size());
		for (i = 0; i<bzel.cmat.size(); i++)
		{
			dNdx[i][0] = 0.; dNdx[i][1] = 0.;
			for (j = 0; j<bzel.cmat[i].size(); j++)
			{
				Nx[i] += bzel.cmat[i][j] * Nx0[j];
				dNdx[i][0] += bzel.cmat[i][j] * dNdx0[j][0];
				dNdx[i][1] += bzel.cmat[i][j] * dNdx0[j][1];
			}
		}
	}
	else if (bzel.order == 4)
	{
		double Nx0[25];
		double dNdx0[25][2];
		BasisFunction4(u, v, bzel.pts4, Nx0, dNdx0, detJ);
		for (i = 0; i < 25; i++)
		{
			pt[0] += Nx0[i] * bzel.pts4[i][0];
			pt[1] += Nx0[i] * bzel.pts4[i][1];
			pt[2] += Nx0[i] * bzel.pts4[i][2];
		}
		Nx.resize(bzel.cmat4.size(), 0);
		dNdx.resize(bzel.cmat4.size());
		for (i = 0; i<bzel.cmat4.size(); i++)
		{
			dNdx[i][0] = 0.; dNdx[i][1] = 0.;
			for (j = 0; j<bzel.cmat4[i].size(); j++)
			{
				Nx[i] += bzel.cmat4[i][j] * Nx0[j];
				dNdx[i][0] += bzel.cmat4[i][j] * dNdx0[j][0];
				dNdx[i][1] += bzel.cmat4[i][j] * dNdx0[j][1];
			}
		}
	}
}


void LeastSquare::ElementMatrix_TSP(const vector<double>& Nx, vector<vector<double>>& EK)
{
	unsigned int i,j;
	for(i=0; i<Nx.size(); i++)
	{
		for(j=0; j<Nx.size(); j++)
		{
			EK[i][j] += Nx[i] * Nx[j];
		}
	}
}

void LeastSquare::ElementForce_TSP(const vector<double>& Nx, double Fb[3], vector<vector<double>>& EF)
{
	for(unsigned int i=0; i<Nx.size(); i++)
	{
		EF[0][i] += Nx[i] * Fb[0];
		EF[1][i] += Nx[i] * Fb[1];
		EF[2][i] += Nx[i] * Fb[2];
	}
}

void LeastSquare::ElementForce_TSP(const vector<double>& Nx, double Fb, vector<double>& EF)
{
	for (unsigned int i = 0; i<Nx.size(); i++)
	{
		EF[i] += Nx[i] * Fb;
	}
}

void LeastSquare::GetFb(double u, double v, const BezierElement& bzel, double Fb[3])
{
	const double r(10.);//hemisphere
	//const double r(100.);//cylinder
	array<double, 3> pt;
	Para2Phys_TSP(u, v, bzel, pt);
	double len = sqrt(pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2]);//hemisphere
	//double len = sqrt(pt[0] * pt[0] + pt[1] * pt[1]);//cylinder
	//if (fabs(len) < 1.e-6)
	//{
	//	cout << u << " " << v << "\n"; 
	//	getchar();
	//}
	len = r / len;
	Fb[0] = len*pt[0];
	Fb[1] = len*pt[1];
	//Fb[2] = len*pt[2];
	//cout << Fb[0] << " " << Fb[1] << " " << Fb[2] << "\n";
	//getchar();
}

void LeastSquare::Assembly_TSP(const vector<vector<double>>& EK, const vector<vector<double>>& EF, const vector<int>& IEN, 
	SparseMatrix<double>& GK, vector<VectorXd>& GF)
{
	unsigned int i,j;
	for(i=0; i<IEN.size(); i++)
	{
		for(j=0; j<IEN.size(); j++)
		{
			if(IDBC[IEN[i]]!=-1 && IDBC[IEN[j]]!=-1)
			{
				GK.coeffRef(IDBC[IEN[i]],IDBC[IEN[j]]) += EK[i][j];
			}
			//else if(IDBC[IEN[i]]!=-1 && IDBC[IEN[j]]==-1)
			//{
			//	GF(IDBC[IEN[i]]) -= EK[i][j]*gh[IEN[j]];
			//}
		}
		//assemble element force
		if(IDBC[IEN[i]]!=-1)
		{
			GF[0](IDBC[IEN[i]]) += EF[0][i];
			GF[1](IDBC[IEN[i]]) += EF[1][i];
			GF[2](IDBC[IEN[i]]) += EF[2][i];
		}
	}
}

void LeastSquare::Assembly_TSP(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN,
	SparseMatrix<double>& GK, VectorXd& GF)
{
	unsigned int i, j;
	for (i = 0; i<IEN.size(); i++)
	{
		for (j = 0; j<IEN.size(); j++)
		{
			if (IDBC[IEN[i]] != -1 && IDBC[IEN[j]] != -1)
			{
				GK.coeffRef(IDBC[IEN[i]], IDBC[IEN[j]]) += EK[i][j];
			}
			//else if(IDBC[IEN[i]]!=-1 && IDBC[IEN[j]]==-1)
			//{
			//	GF(IDBC[IEN[i]]) -= EK[i][j]*gh[IEN[j]];
			//}
		}
		//assemble element force
		if (IDBC[IEN[i]] != -1)
		{
			GF(IDBC[IEN[i]]) += EF[i];
		}
	}
}

void LeastSquare::InitializeSparseMat(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK)
{
	cout << "Initialize sparse matrix...\n";
	cout << "# of eles: " << bzmesh.size() << "\n";
	vector<Triplet<double>> trilist;
	for (int e = 0; e < bzmesh.size(); e++)
	{
		if (e != 0 && e % 500 == 0)
		{
			cout << e << " ";
		}
		for (int i = 0; i < bzmesh[e].IEN.size(); i++)
		{
			for (int j = 0; j < bzmesh[e].IEN.size(); j++)
			{
				int A(bzmesh[e].IEN[i]), B(bzmesh[e].IEN[j]);
				if (A != -1 && B != -1)
					//if (A != -1 && B != -1 && A >= B)
				{
					//if (A >= IDBC.size() || B >= IDBC.size())
					//{
					//	cout << A << " " << B << "\n";
					//	getchar();
					//}
					if (IDBC[A] != -1 && IDBC[B] != -1)
					{
						trilist.push_back(Triplet<double>(IDBC[A], IDBC[B], 0.));
					}
				}
			}
		}
	}
	GK.setFromTriplets(trilist.begin(), trilist.end());
	GK.makeCompressed();
	cout << "done initializing\n";
	cout << "npt: " << npt << "\n";
	cout << "neq: " << neq << "\n";
}

void LeastSquare::BuildLinearSystem_TSP(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK, vector<VectorXd>& GF)
{
	InitializeSparseMat(bzmesh, GK);

	unsigned int e,i,j;
	for(e=0; e<bzmesh.size(); e++)
	{
		//cout << "eid: " << e << "\n";
		vector<vector<double>> EK(bzmesh[e].IEN.size(),vector<double>(bzmesh[e].IEN.size(),0.));
		vector<vector<double>> EF(3, vector<double>(bzmesh[e].IEN.size(), 0.));
		vector<double> Nx;
		vector<array<double,2>> dNdx;
		double detJ, Fb[3];
		array<double,3> x;
		for(i=0; i<Gpt.size(); i++)
		{
			for(j=0; j<Gpt.size(); j++)
			{
				BasisFunction_TSP(Gpt[i],Gpt[j],bzmesh[e],Nx,dNdx,detJ);
				ElementMatrix_TSP(Nx,EK);
				GetFb(Gpt[i], Gpt[j], bzmesh[e], Fb);
				ElementForce_TSP(Nx,Fb,EF);
			}
		}
		Assembly_TSP(EK, EF, bzmesh[e].IEN, GK, GF);
	}
}

void LeastSquare::BuildLinearSystem_TSP(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF)
{
	InitializeSparseMat(bzmesh, GK);

	unsigned int e, i, j;
	for (e = 0; e<bzmesh.size(); e++)
	{
		//cout << "eid: " << e << "\n";
		vector<vector<double>> EK(bzmesh[e].IEN.size(), vector<double>(bzmesh[e].IEN.size(), 0.));
		vector<double> EF(bzmesh[e].IEN.size(), 0.);
		vector<double> Nx;
		vector<array<double, 2>> dNdx;
		double detJ, Fb;
		array<double, 3> x;
		if (bzmesh[e].bc[0] == 1)
		{
			for (i = 0; i < Gpt.size(); i++)
			{
				BasisFunction_TSP(Gpt[i], 0., bzmesh[e], Nx, dNdx, detJ, x.data());
				ElementMatrix_TSP(Nx, EK);
				Fb = ExactValue(x[0], x[1]);
				ElementForce_TSP(Nx, Fb, EF);
			}
		}
		if (bzmesh[e].bc[1] == 1)
		{
			for (i = 0; i < Gpt.size(); i++)
			{
				BasisFunction_TSP(1., Gpt[i], bzmesh[e], Nx, dNdx, detJ, x.data());
				ElementMatrix_TSP(Nx, EK);
				Fb = ExactValue(x[0], x[1]);
				ElementForce_TSP(Nx, Fb, EF);
			}
		}
		if (bzmesh[e].bc[2] == 1)
		{
			for (i = 0; i < Gpt.size(); i++)
			{
				BasisFunction_TSP(Gpt[i], 1., bzmesh[e], Nx, dNdx, detJ, x.data());
				ElementMatrix_TSP(Nx, EK);
				Fb = ExactValue(x[0], x[1]);
				ElementForce_TSP(Nx, Fb, EF);
			}
		}
		if (bzmesh[e].bc[3] == 1)
		{
			for (i = 0; i < Gpt.size(); i++)
			{
				BasisFunction_TSP(0., Gpt[i], bzmesh[e], Nx, dNdx, detJ, x.data());
				ElementMatrix_TSP(Nx, EK);
				Fb = ExactValue(x[0], x[1]);
				ElementForce_TSP(Nx, Fb, EF);
			}
		}
		Assembly_TSP(EK, EF, bzmesh[e].IEN, GK, GF);
	}
}

void LeastSquare::Solver(SparseMatrix<double>& GK, vector<VectorXd>& GF, vector<array<double, 3>>& cp)
{
	SimplicialLDLT<SparseMatrix<double>> solver;
	solver.compute(GK);
	VectorXd solx = solver.solve(GF[0]);
	VectorXd soly = solver.solve(GF[1]);
	//VectorXd solz = solver.solve(GF[2]);
	//cout<<GF[0];
	//getchar();
	//uh.resize(npt);
	cp.resize(npt);
	for(int i=0; i<npt; i++)
	{
		if (IDBC[i] != -1)
		{
			//uh[i] = sol(IDBC[i]);
			cp[i][0] = solx(IDBC[i]);
			cp[i][1] = soly(IDBC[i]);
			//cp[i][2] = solz(IDBC[i]);
			//cout << cp[i][0] << " " << cp[i][1] << " " << cp[i][2] << "\n";
			//getchar();
		}
		else
		{
			//uh[i] = gh[i];
		}
		//cout << uh[i] << '\n';
	}
	//cout << GK << "\n";
	//cout << GF << "\n";
	//getchar();
}

void LeastSquare::Solver(SparseMatrix<double>& GK, VectorXd& GF, vector<double>& gh_out)
{
	SimplicialLDLT<SparseMatrix<double>> solver;
	solver.compute(GK);
	VectorXd sol = solver.solve(GF);
	//cout<<GF[0];
	//getchar();
	gh_out.resize(npt, 0.);
	for (int i = 0; i<npt; i++)
	{
		if (IDBC[i] != -1)
		{
			gh_out[i] = sol(IDBC[i]);
		}
		//cout << gh_out[i] << '\n'; getchar();
	}
	//cout << GK << "\n";
	//cout << GF << "\n";
	//getchar();
}

void LeastSquare::Para2Phys_TSP(double u, double v, const BezierElement& bzel, array<double,3>& pt)
{
	if(bzel.order==3)
	{
		double Nu[4]={(1.-u)*(1.-u)*(1.-u),3.*(1.-u)*(1.-u)*u,3.*(1.-u)*u*u,u*u*u};
		double Nv[4]={(1.-v)*(1.-v)*(1.-v),3.*(1.-v)*(1.-v)*v,3.*(1.-v)*v*v,v*v*v};
		int i,j,loc(0);
		double tmp;
		pt[0]=0.; pt[1]=0.; pt[2]=0.;
		for(i=0;i<4;i++)
		{
			for(j=0;j<4;j++)
			{
				tmp=Nu[j]*Nv[i];
				pt[0]+=tmp*bzel.pts[loc][0];
				pt[1]+=tmp*bzel.pts[loc][1];
				pt[2]+=tmp*bzel.pts[loc][2];
				loc++;
			}
		}
	}
	else if(bzel.order==4)
	{
		double Nu[5]={(1.-u)*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-u)*u,6.*(1.-u)*(1.-u)*u*u,4.*(1.-u)*u*u*u,u*u*u*u};
		double Nv[5]={(1.-v)*(1.-v)*(1.-v)*(1.-v),4.*(1.-v)*(1.-v)*(1.-v)*v,6.*(1.-v)*(1.-v)*v*v,4.*(1.-v)*v*v*v,v*v*v*v};
		int i,j,loc(0);
		double tmp;
		pt[0]=0.; pt[1]=0.; pt[2]=0.;
		for(i=0;i<5;i++)
		{
			for(j=0;j<5;j++)
			{
				tmp=Nu[j]*Nv[i];
				pt[0]+=tmp*bzel.pts4[loc][0];
				pt[1]+=tmp*bzel.pts4[loc][1];
				pt[2]+=tmp*bzel.pts4[loc][2];
				loc++;
			}
		}
	}
}

void LeastSquare::DispCal_TSP(double u,double v,const BezierElement& bzel,double& disp,double& detJ)
{
	if(bzel.order==3)
	{
		double Nx[16];
		double dNdx[16][2];
		BasisFunction(u,v,bzel.pts,Nx,dNdx,detJ);
		double uloc[16];
		unsigned int i,j;
		for(i=0;i<16;i++)
		{
			uloc[i]=0.;
			for(j=0;j<bzel.IEN.size();j++)
			{
				uloc[i]+=bzel.cmat[j][i]*uh[bzel.IEN[j]];
			}
		}
		//displacement
		disp=0.;
		for(i=0;i<16;i++)
		{
			disp+=Nx[i]*uloc[i];
		}
	}
	else if(bzel.order==4)
	{
		double Nx[25];
		double dNdx[25][2];
		BasisFunction4(u,v,bzel.pts4,Nx,dNdx,detJ);
		double uloc[25];
		unsigned int i,j;
		for(i=0;i<25;i++)
		{
			uloc[i]=0.;
			for(j=0;j<bzel.IEN.size();j++)
			{
				uloc[i]+=bzel.cmat4[j][i]*uh[bzel.IEN[j]];
			}
		}
		//displacement
		disp=0.;
		for(i=0;i<25;i++)
		{
			disp+=Nx[i]*uloc[i];
		}
	}
}

void LeastSquare::Quantityh_TSP(double u,double v,const BezierElement& bzel,double val[3],double& detJ)
{
	if(bzel.order==3)
	{
		double Nx[16];
		double dNdx[16][2];
		BasisFunction(u,v,bzel.pts,Nx,dNdx,detJ);
		double uloc[16];
		unsigned int i,j;
		for(i=0;i<16;i++)
		{
			uloc[i]=0.;
			for(j=0;j<bzel.IEN.size();j++)
			{
				uloc[i]+=bzel.cmat[j][i]*uh[bzel.IEN[j]];
			}
		}
		//displacement
		val[0]=0.; val[1]=0.; val[2]=0.;
		for(i=0;i<16;i++)
		{
			val[0]+=Nx[i]*uloc[i];
			val[1]+=dNdx[i][0]*uloc[i];
			val[2]+=dNdx[i][1]*uloc[i];
		}
	}
	else if(bzel.order==4)
	{
		double Nx[25];
		double dNdx[25][2];
		BasisFunction4(u,v,bzel.pts4,Nx,dNdx,detJ);
		double uloc[25];
		unsigned int i,j;
		for(i=0;i<25;i++)
		{
			uloc[i]=0.;
			for(j=0;j<bzel.IEN.size();j++)
			{
				uloc[i]+=bzel.cmat4[j][i]*uh[bzel.IEN[j]];
			}
		}
		//displacement
		val[0]=0.; val[1]=0.; val[2]=0.;
		for(i=0;i<25;i++)
		{
			val[0]+=Nx[i]*uloc[i];
			val[1]+=dNdx[i][0]*uloc[i];
			val[2]+=dNdx[i][1]*uloc[i];
		}
	}
}

//void LeastSquare::VisualizeVTK(string fn,const vector<BezierElement>& bzmesh)
//{
//	vector<array<double,3>> spt;//s means sample
//	vector<array<int,4>> sele;
//	vector<array<double,3>> lpt;//visulize parameter lines
//	vector<array<int,2>> led;//line connectivity
//	vector<double> sdisp;
//	vector<double> errL2;
//	int ns(5),ecount(0);
//	vector<double> su(ns);
//	for(int i=0;i<ns;i++)
//	{
//		su[i]=double(i)/(double(ns)-1.);
//	}
//
//	double L2_all(0.), H1_all(0.);
//	for(unsigned int e=0;e<bzmesh.size();e++)
//	{
//		int loc(0);
//		double L2,H1;
//		ElementError(bzmesh[e],L2,H1);
//		L2_all+=L2; H1_all+=H1;
//		errL2.push_back(sqrt(L2));
//		//errL2.push_back(sqrt(H1));
//		for(int a=0;a<ns;a++)
//		{
//			for(int b=0;b<ns;b++)
//			{
//				double disp, detJ;
//				array<double,3> pt;
//				//Para2Phys(su[b],su[a],bzmesh[e].pts,pt1);
//				//DispCal(su[b],su[a],bzmesh[e],disp,detJ);
//				Para2Phys_TSP(su[b],su[a],bzmesh[e],pt);
//				//DispCal_TSP(su[b],su[a],bzmesh[e],disp,detJ);
//				spt.push_back(pt);
//				//sdisp.push_back(disp);
//				if(a==0||a==ns-1||b==0||b==ns-1)
//				{
//					lpt.push_back(pt);
//				}
//			}
//		}
//		for(int a=0;a<ns-1;a++)
//		{
//			for(int b=0;b<ns-1;b++)
//			{
//				array<int,4> el;
//				el[0]=ecount*ns*ns+a*ns+b;
//				el[1]=ecount*ns*ns+a*ns+b+1;
//				el[2]=ecount*ns*ns+(a+1)*ns+b+1;
//				el[3]=ecount*ns*ns+(a+1)*ns+b;
//				sele.push_back(el);
//			}
//		}
//		for(int a=0;a<ns-1;a++)
//		{
//			array<int,2> lc;
//			lc[0]=ecount*4*(ns-1)+a;
//			lc[1]=ecount*4*(ns-1)+a+1;
//			led.push_back(lc);
//			lc[0]=ecount*4*(ns-1)+3*ns-4+a;
//			lc[1]=ecount*4*(ns-1)+3*ns-4+a+1;
//			led.push_back(lc);
//		}
//		for(int a=0;a<ns-2;a++)
//		{
//			array<int,2> lc;
//			lc[0]=ecount*4*(ns-1)+ns+2*a;
//			lc[1]=ecount*4*(ns-1)+ns+2*a+2;
//			led.push_back(lc);
//			lc[0]=ecount*4*(ns-1)+ns+2*a-1;
//			lc[1]=ecount*4*(ns-1)+ns+2*a+1;
//			led.push_back(lc);
//		}
//		array<int,2> lc1;
//		lc1[0]=ecount*4*(ns-1);
//		lc1[1]=ecount*4*(ns-1)+ns;
//		led.push_back(lc1);
//		lc1[0]=ecount*4*(ns-1)+3*ns-5;
//		lc1[1]=ecount*4*(ns-1)+4*ns-5;
//		led.push_back(lc1);
//		ecount++;
//	}
//
//	cout<<"total L2: "<<sqrt(L2_all)<<"\n";
//	//cout<<"total H1: "<<sqrt(H1_all)<<"\n";
//
//	string fname=fn+".vtk";
//	ofstream fout;
//	fout.open(fname.c_str());
//	unsigned int i;
//	if(fout.is_open())
//	{
//		fout<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
//		fout<<"POINTS "<<spt.size()<<" float\n";
//		for(i=0;i<spt.size();i++)
//		{
//			fout<<spt[i][0]<<" "<<spt[i][1]<<" "<<spt[i][2]<<"\n";
//		}
//		fout<<"\nCELLS "<<sele.size()<<" "<<5*sele.size()<<'\n';
//		for(i=0;i<sele.size();i++)
//		{
//			fout<<"4 "<<sele[i][0]<<" "<<sele[i][1]<<" "<<sele[i][2]<<" "<<sele[i][3]<<'\n';
//		}
//		fout<<"\nCELL_TYPES "<<sele.size()<<'\n';
//		for(i=0;i<sele.size();i++)
//		{
//			fout<<"9\n";
//		}
//		//fout<<"POINT_DATA "<<sdisp.size()<<"\nVECTORS disp float\n";
//		//for(i=0;i<sdisp.size();i++)
//		//{
//		//	fout<<sdisp[i][0]<<" "<<sdisp[i][1]<<" 0\n";
//		//}
//		//fout<<"POINT_DATA "<<sse.size()<<"\nVECTORS strain float\n";
//		//for(i=0;i<sse.size();i++)
//		//{
//		//	fout<<sse[i][0]<<" "<<sse[i][1]<<" "<<sse[i][2]<<"\n";
//		//}
//		//fout<<"POINT_DATA "<<sss.size()<<"\nVECTORS stress float\n";
//		//for(i=0;i<sss.size();i++)
//		//{
//		//	fout<<sss[i][0]<<" "<<sss[i][1]<<" "<<sss[i][2]<<"\n";
//		//}
//
//		//fout<<"POINT_DATA "<<sdisp.size()<<"\nSCALARS u float 1\nLOOKUP_TABLE default\n";
//		//for(i=0;i<sdisp.size();i++)
//		//{
//		//	fout<<sdisp[i]<<"\n";
//		//}
//
//		fout<<"\nCELL_DATA "<<(ns-1)*(ns-1)*errL2.size()<<"\nSCALARS err float 1\nLOOKUP_TABLE default\n";
//		for(i=0;i<errL2.size();i++)
//		{
//			for(int j=0; j<(ns-1)*(ns-1); j++)
//			fout<<errL2[i]<<"\n";
//		}
//		fout.close();
//	}
//	else
//	{
//		cout<<"Cannot open "<<fname<<"!\n";
//	}
//
//	string fname1(fn+"-lines.vtk");
//	ofstream fout1;
//	fout1.open(fname1.c_str());
//	if(fout1.is_open())
//	{
//		fout1<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
//		fout1<<"POINTS "<<lpt.size()<<" float\n";
//		for(uint i=0;i<lpt.size();i++)
//		{
//			fout1<<lpt[i][0]<<" "<<lpt[i][1]<<" "<<lpt[i][2]<<"\n";
//		}
//		fout1<<"\nCELLS "<<led.size()<<" "<<3*led.size()<<'\n';
//		for(uint i=0;i<led.size();i++)
//		{
//			fout1<<"2 "<<led[i][0]<<" "<<led[i][1]<<'\n';
//		}
//		fout1<<"\nCELL_TYPES "<<led.size()<<'\n';
//		for(uint i=0;i<led.size();i++)
//		{
//			fout1<<"3\n";
//		}
//		fout1.close();
//	}
//	else
//	{
//		cout<<"Cannot open "<<fname1<<"!\n";
//	}
//}

void LeastSquare::Run(const vector<BezierElement>& bzmesh, string fn, vector<array<double,3>>& cp_out)
{
	//GaussInfo(4);
	SetSamplingPoints();
	SparseMatrix<double> GK(neq,neq);
	vector<VectorXd> GF(3);
	GK.setZero();
	for (int i = 0; i < 3; i++) GF[i] = VectorXd::Zero(neq);
	cout<<"Building linear system...\n";
	BuildLinearSystem_TSP(bzmesh,GK,GF);
	Solver(GK,GF,cp_out);
	//VisualizeVTK(fn,bzmesh);
}


void LeastSquare::Run_ScalarFitting(const vector<BezierElement>& bzmesh, string fn, vector<double>& gh_out)
{
	//GaussInfo(4);
	SetSamplingPoints();
	SparseMatrix<double> GK(neq, neq);
	VectorXd GF = VectorXd::Zero(neq);
	GK.setZero();
	cout << "Building linear system...\n";
	BuildLinearSystem_TSP(bzmesh, GK, GF);
	Solver(GK, GF, gh_out);
	//VisualizeVTK(fn,bzmesh);
}













double LeastSquare::ExactValue(double x, double y)
{
	//double tmp = ExactValue_1(x, y);
	//double tmp = ExactValue_6(x,y);
	double tmp = ExactValue_7(x, y);
	return tmp;
}

double LeastSquare::ExactValue_1(double x, double y)
{
	return x + y;
}

double LeastSquare::ExactValue_6(double x, double y)
{
	return sin(PI*x) * sin(PI*y);
}

double LeastSquare::ExactValue_7(double x, double y)
{
	return exp((x + y) / 2.);
}






