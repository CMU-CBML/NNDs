#include "Laplace.h"
#include "BSplineBasis.h"
#include <cmath>
#include <fstream>
#include <iostream>

#define PI 3.141592654

using namespace std;

typedef unsigned int uint;

Laplace::Laplace()
{
}

void Laplace::GaussInfo(int ng)
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

void Laplace::BasisFunction(double u, double v, const double pt[16][3], double Nx[16], double dNdx[16][2], double& detJ)
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

void Laplace::BasisFunction4(double u, double v, const double pt[25][3], double Nx[25], double dNdx[25][2], double& detJ)
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

void Laplace::BasisFunction_TSP(double u, double v, const BezierElement& bzel, vector<double>& Nx, vector<array<double,2>>& dNdx, double& detJ)
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

void Laplace::ElementMatrix(double dNdx[16][2], double detJ, double EK[16][16])
{
	int i,j;
	for(i=0;i<16;i++)
	{
		for(j=0;j<16;j++)
		{
			EK[i][j]+=(dNdx[i][0]*dNdx[j][0]+dNdx[i][1]*dNdx[j][1])*detJ;
		}
	}
}

void Laplace::ElementMatrix_TSP(const vector<array<double,2>>& dNdx, double detJ, vector<vector<double>>& EK)
{
	unsigned int i,j;
	for(i=0; i<dNdx.size(); i++)
	{
		for(j=0; j<dNdx.size(); j++)
		{
			EK[i][j]+=(dNdx[i][0]*dNdx[j][0]+dNdx[i][1]*dNdx[j][1])*detJ;
		}
	}
}

void Laplace::ElementForce_TSP(const vector<double>& Nx, double detJ, double Fb, vector<double>& EF)
{
	for(unsigned int i=0; i<Nx.size(); i++)
	{
		EF[i]+=Nx[i]*detJ*Fb;
	}
}

void Laplace::SetProblem(const TruncatedTspline& tts)
{
	npt=0;
	unsigned int i,j,dof;
	BCList.clear();
	gh.clear();
	for(i=0; i<tts.cp.size(); i++)
	{
		if(tts.cp[i].act==1)
		{
			if((tts.cp[i].coor[0]==-1.)||(tts.cp[i].coor[1]==-1.)||
				(tts.cp[i].coor[0]==0.&&tts.cp[i].coor[1]>=0.&&tts.cp[i].coor[1]<=1.)||
				(tts.cp[i].coor[1]==0.&&tts.cp[i].coor[0]>=0.&&tts.cp[i].coor[0]<=1.)||
				(tts.cp[i].coor[1]==1.&&tts.cp[i].coor[0]>=-1.&&tts.cp[i].coor[0]<=0.)||
				(tts.cp[i].coor[0]==1.&&tts.cp[i].coor[1]>=-1.&&tts.cp[i].coor[1]<=0.))
			{
				BCList.push_back(tts.paid[i]);
				double tmp=LDomainSolution(tts.cp[i].coor[0],tts.cp[i].coor[1]);
				gh.push_back(tmp);
			}
			else
			{
				gh.push_back(0.);
			}
			npt++;
		}
	}
}

void Laplace::SetProblem(const vector<int>& IDBC_in, const vector<double>& gh_in)
{
	npt = IDBC_in.size();
	neq = 0;
	IDBC = IDBC_in;
	gh = gh_in;
	for (unsigned int i = 0; i < IDBC.size(); i++)
	{
		if (IDBC[i] != -1) neq++;
		//cout << i << ": " << IDBC[i] << "\n";
	}
}

void Laplace::SetBoundary()
{
	//Dirichilet boundary condition
	IDBC.clear();
	IDBC.resize(npt);
	neq=0;
	for(int i=0; i<npt; i++)
	{
		vector<int>::iterator it=find(BCList.begin(),BCList.end(),i);
		if(it==BCList.end())
		{
			IDBC[i]=neq;
			neq++;
		}
		else
		{
			IDBC[i]=-1;
		}
	}
}

void Laplace::Assembly(double EK1[16][16], double EF1[16], const vector<array<double,16>>& cmat, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF)
{
	unsigned int i,j,a,b,A,B;
	vector<vector<double>> EK(IEN.size(),vector<double>(IEN.size(),0.));
	for(i=0; i<IEN.size(); i++)
	{
		for(j=0; j<IEN.size(); j++)
		{
			for(a=0; a<16; a++)
			{
				for(b=0; b<16; b++)
				{
					EK[i][j] += cmat[i][a]*cmat[j][b]*EK1[a][b];
				}
			}
		}
	}

	for(i=0; i<IEN.size(); i++)
	{
		for(j=0; j<IEN.size(); j++)
		{
			if(IDBC[IEN[i]]!=-1 && IDBC[IEN[j]]!=-1)
			{
				GK.coeffRef(IDBC[IEN[i]],IDBC[IEN[j]]) += EK[i][j];
			}
			else if(IDBC[IEN[i]]!=-1 && IDBC[IEN[j]]==-1)
			{
				GF(IDBC[IEN[i]]) -= EK[i][j]*gh[IEN[j]];
			}
		}
	}
}

void Laplace::Assembly_TSP(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF)
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
			else if(IDBC[IEN[i]]!=-1 && IDBC[IEN[j]]==-1)
			{
				GF(IDBC[IEN[i]]) -= EK[i][j]*gh[IEN[j]];
			}
		}
		//assemble element force
		if(IDBC[IEN[i]]!=-1)
		{
			GF(IDBC[IEN[i]]) += EF[i];
		}
	}
}

void Laplace::BuildLinearSystem(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF)
{
	unsigned int e,i,j,a,b;
	for(e=0; e<bzmesh.size(); e++)
	{
		double EK[16][16];
		double EF[16];
		for(i=0; i<16; i++)
		{
			EF[i]=0.;
			for(j=0; j<16; j++)
			{
				EK[i][j]=0.;
			}
		}
		double Nx[16];
		double dNdx[16][2];
		double detJ;
		for(i=0; i<Gpt.size(); i++)
		{
			for(j=0; j<Gpt.size(); j++)
			{
				BasisFunction(Gpt[i],Gpt[j],bzmesh[e].pts,Nx,dNdx,detJ);
				detJ=wght[i]*wght[j]*detJ;
				ElementMatrix(dNdx,detJ,EK);
			}
		}
		Assembly(EK,EF,bzmesh[e].cmat,bzmesh[e].IEN,GK,GF);
	}

	//int num(0);
	//for(int i=0; i<neq; i++)
	//{
	//	//cout<<GK.coeffRef(i,i)<<' ';
	//	cout<<GF(i)<<' ';
	//	num++;
	//	if(num==5)
	//	{
	//		num=0; cout<<'\n';
	//	}
	//}
	//getchar();
}

void Laplace::BuildLinearSystem_TSP(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF)
{
	InitializeSparseMat(bzmesh, GK);
	unsigned int e,i,j;
	for(e=0; e<bzmesh.size(); e++)
	{
		vector<vector<double>> EK(bzmesh[e].IEN.size(),vector<double>(bzmesh[e].IEN.size(),0.));
		vector<double> Nx, EF(bzmesh[e].IEN.size(),0.);
		vector<array<double,2>> dNdx;
		double detJ, Fb;
		array<double,3> x;
		for(i=0; i<Gpt.size(); i++)
		{
			for(j=0; j<Gpt.size(); j++)
			{
				BasisFunction_TSP(Gpt[i],Gpt[j],bzmesh[e],Nx,dNdx,detJ);
				detJ=wght[i]*wght[j]*detJ;
				ElementMatrix_TSP(dNdx,detJ,EK);
				//element force
				Para2Phys_TSP(Gpt[i],Gpt[j],bzmesh[e],x);
				//Fb=f_source_1(x[0],x[1]);
				//Fb=f_source_2(x[0],x[1]);
				//Fb=f_source_3(x[0],x[1]);
				Fb = f_source(x[0], x[1]);
				ElementForce_TSP(Nx,detJ,Fb,EF);
			}
		}
		//for(unsigned int i=0; i<EF.size(); i++)
		//{
		//	cout<<EF[i]<<" ";
		//}
		//cout<<"\n";
		//getchar();
		Assembly_TSP(EK,EF,bzmesh[e].IEN,GK,GF);
	}
}

void Laplace::Solver(SparseMatrix<double>& GK, VectorXd& GF)
{
	SimplicialLDLT<SparseMatrix<double>> solver;
	VectorXd sol = solver.compute(GK).solve(GF);
	//cout<<sol;
	//getchar();
	uh.resize(npt);
	for(int i=0; i<npt; i++)
	{
		if(IDBC[i]!=-1)
			uh[i]=sol(IDBC[i]);
		else
			uh[i]=gh[i];
		//cout << i <<": "<< uh[i] << '\n';
	}
	//cout << GK << "\n";
	//cout << GF << "\n";
	//getchar();
}

void Laplace::Para2Phys(double u, double v, const double cpt[16][3], double pt[3])
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
			pt[0]+=tmp*cpt[loc][0];
			pt[1]+=tmp*cpt[loc][1];
			pt[2]+=tmp*cpt[loc][2];
			loc++;
		}
	}
}

void Laplace::Para2Phys_TSP(double u, double v, const BezierElement& bzel, array<double,3>& pt)
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

void Laplace::DispCal(double u,double v,const BezierElement& bzel,double& disp,double& detJ)
{
	double Nx[16];
	double dNdx[16][2];
	//double detJ;
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

void Laplace::DispCal_TSP(double u,double v,const BezierElement& bzel,double& disp,double& detJ)
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

void Laplace::Quantityh_TSP(double u,double v,const BezierElement& bzel,double val[3],double& detJ)
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

void Laplace::ElementError(const BezierElement& bzel, double& L2, double& H1)
{
	L2=0.; H1=0.;
	for(unsigned int gid1=0;gid1<Gpt.size();gid1++)
	{
		for(unsigned int gid2=0;gid2<Gpt.size();gid2++)
		{
			double us, ue, detJ, val[3], ux, uy;
			array<double,3> x;
			//DispCal(Gpt[gid1],Gpt[gid2],bzel,us,detJ);
			//Para2Phys(Gpt[gid1],Gpt[gid2],bzel.pts,x);
			//DispCal_TSP(Gpt[gid1],Gpt[gid2],bzel,us,detJ);
			Quantityh_TSP(Gpt[gid1],Gpt[gid2],bzel,val,detJ);
			Para2Phys_TSP(Gpt[gid1],Gpt[gid2],bzel,x);
			//ue=LDomainSolution(x[0],x[1]);
			//LDomainSolutionDeriv(x[0],x[1],ue,ux,uy);
			//ue=ptestHO_sol_1(x[0],x[1]);
			//ue=ptestHO_sol_2(x[0],x[1]);
			//ue=ptestHO_sol_3(x[0],x[1]);
			ue = exact_sol(x[0], x[1]);
			//L2+=wght[gid1]*wght[gid2]*detJ*(us-ue)*(us-ue);
			L2+=wght[gid1]*wght[gid2]*detJ*(val[0]-ue)*(val[0]-ue);
			//H1+=wght[gid1]*wght[gid2]*detJ*((val[1]-ux)*(val[1]-ux)+(val[2]-uy)*(val[2]-uy));
		}
	}
}

void Laplace::VisualizeVTK(string fn, const vector<BezierElement>& bzmesh)
{
	vector<array<double, 3>> spt;//s means sample
	vector<array<int, 4>> sele;
	vector<array<double, 3>> lpt;//visulize parameter lines
	vector<array<int, 2>> led;//line connectivity
	vector<double> sdisp;
	vector<double> errL2;
	int ns(8), ecount(0);
	vector<double> su(ns);
	for (int i = 0; i < ns; i++)
	{
		su[i] = double(i) / (double(ns) - 1.);
	}

	double L2_all(0.), H1_all(0.);
	double hmax(0.);
	for (unsigned int e = 0; e < bzmesh.size(); e++)
	{
		double pd1[2] = { bzmesh[e].pts[15][0] - bzmesh[e].pts[0][0],bzmesh[e].pts[15][1] - bzmesh[e].pts[0][1] };
		double h1 = sqrt(pd1[0] * pd1[0] + pd1[1] * pd1[1]);
		double pd2[2] = { bzmesh[e].pts[12][0] - bzmesh[e].pts[3][0],bzmesh[e].pts[12][1] - bzmesh[e].pts[3][1] };
		double h2 = sqrt(pd2[0] * pd2[0] + pd2[1] * pd2[1]);
		if (h1 < h2) h1 = h2;
		if (hmax < h1) hmax = h1;

		int loc(0);
		double L2, H1;
		//ElementError(bzmesh[e],L2,H1);
		//L2_all+=L2; H1_all+=H1;
		//errL2.push_back(sqrt(L2));
		//errL2.push_back(sqrt(H1));
		for (int a = 0; a < ns; a++)
		{
			for (int b = 0; b < ns; b++)
			{
				double disp, detJ;
				array<double, 3> pt;
				//Para2Phys(su[b],su[a],bzmesh[e].pts,pt1);
				//DispCal(su[b],su[a],bzmesh[e],disp,detJ);
				Para2Phys_TSP(su[b], su[a], bzmesh[e], pt);
				//DispCal_TSP(su[b],su[a],bzmesh[e],disp,detJ);
				spt.push_back(pt);
				//sdisp.push_back(disp);
				if (a == 0 || a == ns - 1 || b == 0 || b == ns - 1)
				{
					lpt.push_back(pt);
				}
			}
		}
		for (int a = 0; a < ns - 1; a++)
		{
			for (int b = 0; b < ns - 1; b++)
			{
				array<int, 4> el;
				el[0] = ecount * ns*ns + a * ns + b;
				el[1] = ecount * ns*ns + a * ns + b + 1;
				el[2] = ecount * ns*ns + (a + 1)*ns + b + 1;
				el[3] = ecount * ns*ns + (a + 1)*ns + b;
				sele.push_back(el);
			}
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> lc;
			lc[0] = ecount * 4 * (ns - 1) + a;
			lc[1] = ecount * 4 * (ns - 1) + a + 1;
			led.push_back(lc);
			lc[0] = ecount * 4 * (ns - 1) + 3 * ns - 4 + a;
			lc[1] = ecount * 4 * (ns - 1) + 3 * ns - 4 + a + 1;
			led.push_back(lc);
		}
		for (int a = 0; a < ns - 2; a++)
		{
			array<int, 2> lc;
			lc[0] = ecount * 4 * (ns - 1) + ns + 2 * a;
			lc[1] = ecount * 4 * (ns - 1) + ns + 2 * a + 2;
			led.push_back(lc);
			lc[0] = ecount * 4 * (ns - 1) + ns + 2 * a - 1;
			lc[1] = ecount * 4 * (ns - 1) + ns + 2 * a + 1;
			led.push_back(lc);
		}
		array<int, 2> lc1;
		lc1[0] = ecount * 4 * (ns - 1);
		lc1[1] = ecount * 4 * (ns - 1) + ns;
		led.push_back(lc1);
		lc1[0] = ecount * 4 * (ns - 1) + 3 * ns - 5;
		lc1[1] = ecount * 4 * (ns - 1) + 4 * ns - 5;
		led.push_back(lc1);
		ecount++;
	}

	cout << "hmax: " << hmax << "\n";
	cout << "total L2: " << sqrt(L2_all) << "\n";
	//cout<<"total H1: "<<sqrt(H1_all)<<"\n";

	string fname = fn + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i < spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
		}
		//fout<<"POINT_DATA "<<sdisp.size()<<"\nVECTORS disp float\n";
		//for(i=0;i<sdisp.size();i++)
		//{
		//	fout<<sdisp[i][0]<<" "<<sdisp[i][1]<<" 0\n";
		//}
		//fout<<"POINT_DATA "<<sse.size()<<"\nVECTORS strain float\n";
		//for(i=0;i<sse.size();i++)
		//{
		//	fout<<sse[i][0]<<" "<<sse[i][1]<<" "<<sse[i][2]<<"\n";
		//}
		//fout<<"POINT_DATA "<<sss.size()<<"\nVECTORS stress float\n";
		//for(i=0;i<sss.size();i++)
		//{
		//	fout<<sss[i][0]<<" "<<sss[i][1]<<" "<<sss[i][2]<<"\n";
		//}

		//fout<<"POINT_DATA "<<sdisp.size()<<"\nSCALARS u float 1\nLOOKUP_TABLE default\n";
		//for(i=0;i<sdisp.size();i++)
		//{
		//	fout<<sdisp[i]<<"\n";
		//}

		//fout<<"\nCELL_DATA "<<(ns-1)*(ns-1)*errL2.size()<<"\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		//for(i=0;i<errL2.size();i++)
		//{
		//	for(int j=0; j<(ns-1)*(ns-1); j++)
		//	fout<<errL2[i]<<"\n";
		//}

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	string fname1(fn + "-lines.vtk");
	ofstream fout1;
	fout1.open(fname1.c_str());
	if (fout1.is_open())
	{
		fout1 << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout1 << "POINTS " << lpt.size() << " float\n";
		for (uint i = 0; i < lpt.size(); i++)
		{
			fout1 << lpt[i][0] << " " << lpt[i][1] << " " << lpt[i][2] << "\n";
		}
		fout1 << "\nCELLS " << led.size() << " " << 3 * led.size() << '\n';
		for (uint i = 0; i < led.size(); i++)
		{
			fout1 << "2 " << led[i][0] << " " << led[i][1] << '\n';
		}
		fout1 << "\nCELL_TYPES " << led.size() << '\n';
		for (uint i = 0; i < led.size(); i++)
		{
			fout1 << "3\n";
		}
		fout1.close();
	}
	else
	{
		cout << "Cannot open " << fname1 << "!\n";
	}
}

//void Laplace::DispCal_Coupling_Bezier(double u, double v, double w, const BezierElement& bzel, double pt[3], double& disp, double dudx[3], double& detJ)
//{
//	if (bzel.bzflag == 0)
//	{
//		vector<double> Nx(bzel.IEN.size());
//		vector<array<double, 3>> dNdx(bzel.IEN.size());
//		bzel.Para2Phys(u, v, w, pt);
//		BasisFunction_IGA(u, v, w, bzel, Nx, dNdx, detJ);
//		disp = 0.;
//		dudx[0] = 0.; dudx[1] = 0.; dudx[2] = 0.;
//		for (uint i = 0; i < bzel.IEN.size(); i++)
//		{
//			disp += Nx[i] * uh[bzel.IEN[i]];
//			dudx[0] += dNdx[i][0] * uh[bzel.IEN[i]];
//			dudx[1] += dNdx[i][1] * uh[bzel.IEN[i]];
//			dudx[2] += dNdx[i][2] * uh[bzel.IEN[i]];
//		}
//	}
//	else
//	{
//		vector<double> Nx(bzel.IENb.size());
//		vector<array<double, 3>> dNdx(bzel.IENb.size());
//		bzel.Para2Phys(u, v, w, pt);
//		BasisFunction_Bezier(u, v, w, bzel.pts, Nx, dNdx, detJ);
//		disp = 0.;
//		dudx[0] = 0.; dudx[1] = 0.; dudx[2] = 0.;
//		for (uint i = 0; i < bzel.IENb.size(); i++)
//		{
//			disp += Nx[i] * uh[bzel.IENb[i]];
//			dudx[0] += dNdx[i][0] * uh[bzel.IENb[i]];
//			dudx[1] += dNdx[i][1] * uh[bzel.IENb[i]];
//			dudx[2] += dNdx[i][2] * uh[bzel.IENb[i]];
//		}
//	}
//}
//
//void Laplace::VisualizeVTK3D(string fn, const vector<BezierElement>& bzmesh)
//{
//	vector<array<double, 3>> spt;//sample points
//	vector<double> sdisp;
//	vector<array<int, 8>> sele;
//	vector<double> errL2;
//	vector<array<double, 3>> lpt;//visulize parameter lines
//	vector<array<int, 2>> led;//line connectivity
//	double detJ;
//
//	for (unsigned int e = 0; e < bzmesh.size(); e++)
//	{
//		int ns(4);
//		if (bzmesh[e].type == 1) ns = 5;
//		vector<double> su(ns);
//		for (int i = 0; i < ns; i++)
//		{
//			su[i] = double(i) / (double(ns) - 1.);
//		}
//
//		int loc(0);
//		int pstart = spt.size();
//		for (int a = 0; a < ns; a++)
//		{
//			for (int b = 0; b < ns; b++)
//			{
//				for (int c = 0; c < ns; c++)
//				{
//					double pt1[3], dudx[3];
//					double disp;
//					//bzmesh[e].Para2Phys(su[c], su[b], su[a], pt1);
//					DispCal_Coupling_Bezier(su[c], su[b], su[a], bzmesh[e], pt1, disp, dudx, detJ);
//					array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
//					spt.push_back(pt);
//					sdisp.push_back(disp);
//				}
//			}
//		}
//		int nns[2] = { ns*ns*ns, ns*ns };
//		for (int a = 0; a < ns - 1; a++)
//		{
//			for (int b = 0; b < ns - 1; b++)
//			{
//				for (int c = 0; c < ns - 1; c++)
//				{
//					array<int, 8> el;
//					el[0] = pstart + a * nns[1] + b * ns + c;
//					el[1] = pstart + a * nns[1] + b * ns + c + 1;
//					el[2] = pstart + a * nns[1] + (b + 1)*ns + c + 1;
//					el[3] = pstart + a * nns[1] + (b + 1)*ns + c;
//					el[4] = pstart + (a + 1)*nns[1] + b * ns + c;
//					el[5] = pstart + (a + 1)*nns[1] + b * ns + c + 1;
//					el[6] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
//					el[7] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c;
//					sele.push_back(el);
//				}
//			}
//		}
//		//edges
//		int lstart = lpt.size();
//		for (int a = 0; a < ns; a++)
//		{
//			double pt1[3];
//			bzmesh[e].Para2Phys(su[a], 0., 0., pt1);
//			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
//			lpt.push_back(pt);
//		}
//		for (int a = 0; a < ns - 1; a++)
//		{
//			array<int, 2> ed = { lstart + a, lstart + a + 1 };
//			led.push_back(ed);
//		}
//		lstart = lpt.size();
//		for (int a = 0; a < ns; a++)
//		{
//			double pt1[3];
//			bzmesh[e].Para2Phys(su[a], 1., 0., pt1);
//			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
//			lpt.push_back(pt);
//		}
//		for (int a = 0; a < ns - 1; a++)
//		{
//			array<int, 2> ed = { lstart + a, lstart + a + 1 };
//			led.push_back(ed);
//		}
//		lstart = lpt.size();
//		for (int a = 0; a < ns; a++)
//		{
//			double pt1[3];
//			bzmesh[e].Para2Phys(su[a], 0., 1., pt1);
//			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
//			lpt.push_back(pt);
//		}
//		for (int a = 0; a < ns - 1; a++)
//		{
//			array<int, 2> ed = { lstart + a, lstart + a + 1 };
//			led.push_back(ed);
//		}
//		lstart = lpt.size();
//		for (int a = 0; a < ns; a++)
//		{
//			double pt1[3];
//			bzmesh[e].Para2Phys(su[a], 1., 1., pt1);
//			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
//			lpt.push_back(pt);
//		}
//		for (int a = 0; a < ns - 1; a++)
//		{
//			array<int, 2> ed = { lstart + a, lstart + a + 1 };
//			led.push_back(ed);
//		}
//		lstart = lpt.size();
//		for (int a = 0; a < ns; a++)
//		{
//			double pt1[3];
//			bzmesh[e].Para2Phys(0., su[a], 0., pt1);
//			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
//			lpt.push_back(pt);
//		}
//		for (int a = 0; a < ns - 1; a++)
//		{
//			array<int, 2> ed = { lstart + a, lstart + a + 1 };
//			led.push_back(ed);
//		}
//		lstart = lpt.size();
//		for (int a = 0; a < ns; a++)
//		{
//			double pt1[3];
//			bzmesh[e].Para2Phys(1., su[a], 0., pt1);
//			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
//			lpt.push_back(pt);
//		}
//		for (int a = 0; a < ns - 1; a++)
//		{
//			array<int, 2> ed = { lstart + a, lstart + a + 1 };
//			led.push_back(ed);
//		}
//		lstart = lpt.size();
//		for (int a = 0; a < ns; a++)
//		{
//			double pt1[3];
//			bzmesh[e].Para2Phys(0., su[a], 1., pt1);
//			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
//			lpt.push_back(pt);
//		}
//		for (int a = 0; a < ns - 1; a++)
//		{
//			array<int, 2> ed = { lstart + a, lstart + a + 1 };
//			led.push_back(ed);
//		}
//		lstart = lpt.size();
//		for (int a = 0; a < ns; a++)
//		{
//			double pt1[3];
//			bzmesh[e].Para2Phys(1., su[a], 1., pt1);
//			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
//			lpt.push_back(pt);
//		}
//		for (int a = 0; a < ns - 1; a++)
//		{
//			array<int, 2> ed = { lstart + a, lstart + a + 1 };
//			led.push_back(ed);
//		}
//		lstart = lpt.size();
//		for (int a = 0; a < ns; a++)
//		{
//			double pt1[3];
//			bzmesh[e].Para2Phys(0., 0., su[a], pt1);
//			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
//			lpt.push_back(pt);
//		}
//		for (int a = 0; a < ns - 1; a++)
//		{
//			array<int, 2> ed = { lstart + a, lstart + a + 1 };
//			led.push_back(ed);
//		}
//		lstart = lpt.size();
//		for (int a = 0; a < ns; a++)
//		{
//			double pt1[3];
//			bzmesh[e].Para2Phys(1., 0., su[a], pt1);
//			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
//			lpt.push_back(pt);
//		}
//		for (int a = 0; a < ns - 1; a++)
//		{
//			array<int, 2> ed = { lstart + a, lstart + a + 1 };
//			led.push_back(ed);
//		}
//		lstart = lpt.size();
//		for (int a = 0; a < ns; a++)
//		{
//			double pt1[3];
//			bzmesh[e].Para2Phys(0., 1., su[a], pt1);
//			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
//			lpt.push_back(pt);
//		}
//		for (int a = 0; a < ns - 1; a++)
//		{
//			array<int, 2> ed = { lstart + a, lstart + a + 1 };
//			led.push_back(ed);
//		}
//		lstart = lpt.size();
//		for (int a = 0; a < ns; a++)
//		{
//			double pt1[3];
//			bzmesh[e].Para2Phys(1., 1., su[a], pt1);
//			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
//			lpt.push_back(pt);
//		}
//		for (int a = 0; a < ns - 1; a++)
//		{
//			array<int, 2> ed = { lstart + a, lstart + a + 1 };
//			led.push_back(ed);
//		}
//	}
//
//	string fname = fn + ".vtk";
//	ofstream fout;
//	fout.open(fname.c_str());
//	unsigned int i;
//	if (fout.is_open())
//	{
//		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
//		fout << "POINTS " << spt.size() << " float\n";
//		for (i = 0; i < spt.size(); i++)
//		{
//			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
//		}
//		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
//		for (i = 0; i < sele.size(); i++)
//		{
//			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << '\n';
//		}
//		fout << "\nCELL_TYPES " << sele.size() << '\n';
//		for (i = 0; i < sele.size(); i++)
//		{
//			fout << "9\n";
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
//		//fout<<"\nCELL_DATA "<<(ns-1)*(ns-1)*errL2.size()<<"\nSCALARS err float 1\nLOOKUP_TABLE default\n";
//		//for(i=0;i<errL2.size();i++)
//		//{
//		//	for(int j=0; j<(ns-1)*(ns-1); j++)
//		//	fout<<errL2[i]<<"\n";
//		//}
//
//		fout.close();
//	}
//	else
//	{
//		cout << "Cannot open " << fname << "!\n";
//	}
//
//	string fname1(fn + "-lines.vtk");
//	ofstream fout1;
//	fout1.open(fname1.c_str());
//	if (fout1.is_open())
//	{
//		fout1 << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
//		fout1 << "POINTS " << lpt.size() << " float\n";
//		for (uint i = 0; i < lpt.size(); i++)
//		{
//			fout1 << lpt[i][0] << " " << lpt[i][1] << " " << lpt[i][2] << "\n";
//		}
//		fout1 << "\nCELLS " << led.size() << " " << 3 * led.size() << '\n';
//		for (uint i = 0; i < led.size(); i++)
//		{
//			fout1 << "2 " << led[i][0] << " " << led[i][1] << '\n';
//		}
//		fout1 << "\nCELL_TYPES " << led.size() << '\n';
//		for (uint i = 0; i < led.size(); i++)
//		{
//			fout1 << "3\n";
//		}
//		fout1.close();
//	}
//	else
//	{
//		cout << "Cannot open " << fname1 << "!\n";
//	}
//}

// Backup
//void Laplace::VisualizeVTK(string fn,const vector<BezierElement>& bzmesh)
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
//	double hmax(0.);
//	for(unsigned int e=0;e<bzmesh.size();e++)
//	{
//		double pd1[2] = { bzmesh[e].pts[15][0] - bzmesh[e].pts[0][0],bzmesh[e].pts[15][1] - bzmesh[e].pts[0][1] };
//		double h1 = sqrt(pd1[0] * pd1[0] + pd1[1] * pd1[1]);
//		double pd2[2] = { bzmesh[e].pts[12][0] - bzmesh[e].pts[3][0],bzmesh[e].pts[12][1] - bzmesh[e].pts[3][1] };
//		double h2 = sqrt(pd2[0] * pd2[0] + pd2[1] * pd2[1]);
//		if (h1 < h2) h1 = h2;
//		if (hmax < h1) hmax = h1;
//
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
//				DispCal_TSP(su[b],su[a],bzmesh[e],disp,detJ);
//				spt.push_back(pt);
//				sdisp.push_back(disp);
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
//	cout << "hmax: " << hmax << "\n";
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
//		fout<<"POINT_DATA "<<sdisp.size()<<"\nSCALARS u float 1\nLOOKUP_TABLE default\n";
//		for(i=0;i<sdisp.size();i++)
//		{
//			fout<<sdisp[i]<<"\n";
//		}
//
//		//fout<<"\nCELL_DATA "<<(ns-1)*(ns-1)*errL2.size()<<"\nSCALARS err float 1\nLOOKUP_TABLE default\n";
//		//for(i=0;i<errL2.size();i++)
//		//{
//		//	for(int j=0; j<(ns-1)*(ns-1); j++)
//		//	fout<<errL2[i]<<"\n";
//		//}
//
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

void Laplace::Run(const vector<BezierElement>& bzmesh, string fn)
{
	GaussInfo(4);
	//SetBoundary();
	SparseMatrix<double> GK(neq,neq);
	VectorXd GF(neq);
	GK.setZero();
	GF.setZero();
	cout<<"Building linear system...\n";
	//BuildLinearSystem(bzmesh,GK,GF);
	BuildLinearSystem_TSP(bzmesh,GK,GF);
	Solver(GK,GF);
	VisualizeVTK(fn,bzmesh);
}

void Laplace::setProblem_ptestHO(const TruncatedTspline& tts)
{
	npt=0;
	unsigned int i,j,dof;
	BCList.clear();
	gh.clear();
	for(i=0; i<tts.cp.size(); i++)
	{
		if(tts.cp[i].act==1)
		{
			if((tts.cp[i].coor[0]==0.)||
				(tts.cp[i].coor[1]==0.)||
				(tts.cp[i].coor[0]==1.)||
				(tts.cp[i].coor[1]==1.))
			{
				BCList.push_back(tts.paid[i]);
				//double tmp=ptestHO_sol_1(tts.cp[i].coor[0],tts.cp[i].coor[1]);
				//double tmp=ptestHO_sol_2(tts.cp[i].coor[0],tts.cp[i].coor[1]);
				double tmp=ptestHO_sol_3(tts.cp[i].coor[0],tts.cp[i].coor[1]);
				gh.push_back(tmp);
			}
			else
			{
				gh.push_back(0.);
			}
			npt++;
		}
	}
}


////////////////////////////////////////////////////////////////////////////

void Laplace::SetProblem_Dual(const vector<int>& IDBC_in, const vector<double>& gh_in)
{
	IDBC = IDBC_in;
	gh = gh_in;
	npt = IDBC.size();
	neq = 0;
	for (uint i = 0; i < IDBC.size(); i++)
	{
		if (IDBC[i] != -1)
		{
			neq++;
		}
	}
	//cout << npt << " " << neq << "\n"; getchar();
}

void Laplace::BasisFunction_Dual(double u, double v, const BezierElement& bzel, 
	vector<double>& Nx, vector<array<double, 2>>& dNdx, double& detJ)
{
	//Nx.clear();
	//dNdx.clear();
	unsigned int i, j;
	if (bzel.type == 0)
	{
		double Nx0[16];
		double dNdx0[16][2];
		BasisFunction(u, v, bzel.pts, Nx0, dNdx0, detJ);
		//Nx.resize(bzel.cmat.size(), 0);
		//dNdx.resize(bzel.cmat.size());
		for (i = 0; i<bzel.cmat.size(); i++)
		{
			Nx[i] = 0.;
			dNdx[i][0] = 0.; dNdx[i][1] = 0.;
			for (j = 0; j<bzel.cmat[i].size(); j++)
			{
				Nx[i] += bzel.cmat[i][j] * Nx0[j];
				dNdx[i][0] += bzel.cmat[i][j] * dNdx0[j][0];
				dNdx[i][1] += bzel.cmat[i][j] * dNdx0[j][1];
			}
			//cout << dNdx[i][0] << " " << dNdx[i][1] << "\n";
		}
		//getchar();
	}
	else if (bzel.type == 4)
	{
		BasisFunctionDual_Irr(u, v, bzel, Nx, dNdx, detJ);
	}
}

void Laplace::BasisFunctionDual_Irr(double u, double v, const BezierElement& bzel,
	vector<double>& Nx, vector<array<double, 2>>& dNdx, double& detJ)
{
	double tol(1.e-10), eps(1.e-8);
	double u1 = max(tol, u);
	double v1 = max(tol, v);
	int n = floor(min(-log2(u1), -log2(v1))) + 1, subid;
	double pow2 = pow(2., double(n - 1));
	u1 *= pow2; v1 *= pow2;
	if (v1<0.5)
	{
		subid = 0; u1 = 2.*u1 - 1.; v1 = 2.*v1;
	}
	else if (u1<0.5)
	{
		subid = 2; u1 = 2.*u1; v1 = 2.*v1 - 1.;
	}
	else
	{
		subid = 1; u1 = 2.*u1 - 1.; v1 = 2.*v1 - 1.;
	}
	//cout << "subid: " << subid << "\n";
	//cout << "uv: " << u1 << " " << v1 << "\n";
	//evaluate first nv+5 and the remaining 7 functions separately
	//the last 7 functions
	vector<double> ku(5, 0.), kv(5, 0.), uval, vval;
	BSplineBasis bu, bv;
	int nv(bzel.IEN.size() - 12);
	//for (uint i = nv; i < tmesh[eid].IEN.size(); i++)
	//{
	//	ku.assign(tmesh[eid].patch_ku[i].begin(), tmesh[eid].patch_ku[i].end());
	//	kv.assign(tmesh[eid].patch_kv[i].begin(), tmesh[eid].patch_kv[i].end());
	//	bu.Set(3, ku);
	//	bv.Set(3, kv);
	//	bu.BasisFunction(0, u, 1, uval);
	//	bv.BasisFunction(0, v, 1, vval);
	//	Nt[i] = uval[0] * vval[0];
	//	dNdt[i][0] = uval[1] * vval[0];
	//	dNdt[i][1] = uval[0] * vval[1];
	//}
	//the first nv functions
	//cout << n << " " << subid << " " << u1 << " " << v1 << "\n";
	//double sum(0.);
	double Nt1[16];
	double dNdt1[16][2];
	double subku[3][8];
	double subkv[3][8];
	double knot[9] = { -2.,-1.,0.,0.,1.,2.,3.,4.,5. };
	int ist[3][2] = { { 1,0 },{ 1,1 },{ 0,1 } };
	double shift[3][2] = { { -1.,0. },{ -1.,-1. },{ 0.,-1. } };
	int loc(0);
	pow2 *= 2.;
	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int k = 0; k < 5; k++)
			{
				ku[k] = knot[ist[subid][0] + i + k] + shift[subid][0];
				kv[k] = knot[ist[subid][1] + j + k] + shift[subid][1];
			}
			bu.Set(3, ku);
			bv.Set(3, kv);
			bu.BasisFunction(0, u1, 1, uval);
			bv.BasisFunction(0, v1, 1, vval);
			Nt1[loc] = uval[0] * vval[0];
			//dNdt1[loc][0] = uval[1] * vval[0]/pow2;
			//dNdt1[loc][1] = uval[0] * vval[1]/pow2;
			dNdt1[loc][0] = uval[1] * vval[0];
			dNdt1[loc][1] = uval[0] * vval[1];
			//cout << dNdt1[loc][0] << " " << dNdt1[loc][1] << "\n"; getchar();
			loc++;
		}
	}
	int Pk[3][16] = { { nv - 1,nv,nv + 5,nv + 12,0,nv + 1,nv + 6,nv + 13,nv + 3,nv + 2,nv + 7,nv + 14,nv + 10,nv + 9,nv + 8,nv + 15 },
	{ 0,nv + 1,nv + 6,nv + 13,nv + 3,nv + 2,nv + 7,nv + 14,nv + 10,nv + 9,nv + 8,nv + 15,nv + 19,nv + 18,nv + 17,nv + 16 },
	{ 1,0,nv + 1,nv + 6,nv + 4,nv + 3,nv + 2,nv + 7,nv + 11,nv + 10,nv + 9,nv + 8,nv + 20,nv + 19,nv + 18,nv + 17 } };
	//cout << "sum: " << sum << "\n";
	MatrixXd mat(bzel.IEN.size(), 16);
	for (uint i = 0; i < bzel.IEN.size(); i++)
	{
		for (int j = 0; j < 16; j++)
		{
			mat(i, j) = bzel.abar(i, Pk[subid][j]);
		}
	}
	//cout << "eid: " << eid << "\n";
	//cout << "IEN size: " << tmesh[eid].IEN.size() << "\n";
	//cout << "nv: " << nv << "\n";
	//cout << "mat: " << mat.rows() << " " << mat.cols() << "\n";
	//cout << "smat[0]: " << tmesh[eid].smat[0].rows() << " " << tmesh[eid].smat[0].rows() << "\n";
	//getchar();
	//if (n > 10)
	//{
	//	cout << "n: " << n << "\n";
	//}
	for (int i = 0; i < n - 1; i++)
	{
		mat = bzel.amat * mat;
	}
	vector<array<double, 2>> dNdt(bzel.IEN.size());
	double dxdt[2][2] = { {0.,0.},{0.,0.} };
	for (uint i = 0; i < bzel.IEN.size(); i++)
	{
		Nx[i] = 0.;
		dNdt[i][0] = 0.; dNdt[i][1] = 0.;
		for (int j = 0; j < 16; j++)
		{
			Nx[i] += mat(i, j)*Nt1[j];
			dNdt[i][0] += mat(i, j)*dNdt1[j][0];
			dNdt[i][1] += mat(i, j)*dNdt1[j][1];
		}
		dNdt[i][0] *= pow2;
		dNdt[i][1] *= pow2;
		//cout << dNdt[i][0] << " " << dNdt[i][1] << "\n"; getchar();
		for (int a = 0; a < 2; a++)
		{
			for (int b = 0; b < 2; b++)
			{
				dxdt[a][b] += bzel.cp[i][a]*dNdt[i][b];
			}
		}
	}
	//cout << dxdt[0][0] << " " << dxdt[0][1] << "\n";
	//cout << dxdt[1][0] << " " << dxdt[1][1] << "\n";
	detJ = dxdt[0][0] * dxdt[1][1] - dxdt[0][1] * dxdt[1][0];
	//cout << detJ << "\n"; getchar();
	double dtdx[2][2] = { { dxdt[1][1] / detJ,-dxdt[0][1] / detJ },{ -dxdt[1][0] / detJ,dxdt[0][0] / detJ } };
	for (uint i = 0; i<bzel.IEN.size(); i++)
	{
		dNdx[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0];
		dNdx[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1];
		//cout << dNdx[i][0] << " " << dNdx[i][1] << "\n";
		//getchar();
	}
	detJ = 0.25*detJ;
}

void Laplace::Para2Phys_Dual(double u, double v, const BezierElement& bzel, const vector<double>& Nx, array<double, 3>& pt)
{
	if (bzel.type == 0)
	{
		double Nu[4] = { (1. - u)*(1. - u)*(1. - u),3.*(1. - u)*(1. - u)*u,3.*(1. - u)*u*u,u*u*u };
		double Nv[4] = { (1. - v)*(1. - v)*(1. - v),3.*(1. - v)*(1. - v)*v,3.*(1. - v)*v*v,v*v*v };
		int i, j, loc(0);
		double tmp;
		pt[0] = 0.; pt[1] = 0.; pt[2] = 0.;
		for (i = 0; i<4; i++)
		{
			for (j = 0; j<4; j++)
			{
				tmp = Nu[j] * Nv[i];
				pt[0] += tmp*bzel.pts[loc][0];
				pt[1] += tmp*bzel.pts[loc][1];
				pt[2] += tmp*bzel.pts[loc][2];
				loc++;
			}
		}
	}
	else if (bzel.type == 4)
	{
		pt[0] = 0.; pt[1] = 0.; pt[2] = 0.;
		for (uint i = 0; i < bzel.IEN.size(); i++)
		{
			pt[0] += Nx[i] * bzel.cp[i][0];
			pt[1] += Nx[i] * bzel.cp[i][1];
			pt[2] += Nx[i] * bzel.cp[i][2];
		}
	}
}

void Laplace::InitializeSparseMat(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK)
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

void Laplace::BuildLinearSystem_Dual(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF)
{
	InitializeSparseMat(bzmesh, GK);
	unsigned int e, i, j;
	int lev_intg(15);//levels of subdivision for integration, 0 means no subdivision
	for (e = 0; e<bzmesh.size(); e++)
	{
		if (e != 0 && e % 100 == 0)
		{
			cout << e << " ";
		}
		vector<vector<double>> EK(bzmesh[e].IEN.size(), vector<double>(bzmesh[e].IEN.size(), 0.));
		vector<double> Nx(bzmesh[e].IEN.size()), EF(bzmesh[e].IEN.size(), 0.);
		vector<array<double, 2>> dNdx(bzmesh[e].IEN.size());
		double detJ, Fb;
		array<double, 3> x;
		//cout << "eid: " << e << "\n";
		if (bzmesh[e].type == 0)
		{
			for (i = 0; i<Gpt.size(); i++)
			{
				for (j = 0; j<Gpt.size(); j++)
				{
					BasisFunction_Dual(Gpt[i], Gpt[j], bzmesh[e], Nx, dNdx, detJ);
					detJ = wght[i] * wght[j] * detJ;
					ElementMatrix_TSP(dNdx, detJ, EK);
					//element force
					Para2Phys_Dual(Gpt[i], Gpt[j], bzmesh[e], Nx, x);
					Fb = f_source(x[0], x[1]);
					ElementForce_TSP(Nx, detJ, Fb, EF);
				}
			}
		}
		else if (bzmesh[e].type == 4)
		{
			//if(lev_intg==0)
			//{
			//	for (i = 0; i<Gpt.size(); i++)
			//	{
			//		for (j = 0; j<Gpt.size(); j++)
			//		{
			//			BasisFunction_Dual(Gpt[i], Gpt[j], bzmesh[e], Nx, dNdx, detJ);
			//			detJ = wght[i] * wght[j] * detJ;
			//			ElementMatrix_TSP(dNdx, detJ, EK);
			//			//element force
			//			Para2Phys_Dual(Gpt[i], Gpt[j], bzmesh[e], Nx, x);
			//			Fb = f_source(x[0], x[1]);
			//			ElementForce_TSP(Nx, detJ, Fb, EF);
			//		}
			//	}
			//}
			//else//greater than 1
			{
				double pow2, u[2];
				for (int it = 1; it <= lev_intg; it++)
				{
					pow2 = pow(2., double(it));
					//sub-region 1
					for (i = 0; i<Gpt.size(); i++)//u
					{
						for (j = 0; j<Gpt.size(); j++)//v
						{
							u[0] = (Gpt[i] + 1.) / pow2;
							u[1] = Gpt[j] / pow2;
							BasisFunction_Dual(u[0], u[1], bzmesh[e], Nx, dNdx, detJ);
							detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
							ElementMatrix_TSP(dNdx, detJ, EK);
							//element force
							Para2Phys_Dual(u[0], u[1], bzmesh[e], Nx, x);
							Fb = f_source(x[0], x[1]);
							ElementForce_TSP(Nx, detJ, Fb, EF);
						}
					}
					//sub-region 2
					for (i = 0; i<Gpt.size(); i++)//u
					{
						for (j = 0; j<Gpt.size(); j++)//v
						{
							u[0] = (Gpt[i] + 1.) / pow2;
							u[1] = (Gpt[j] + 1.) / pow2;
							BasisFunction_Dual(u[0], u[1], bzmesh[e], Nx, dNdx, detJ);
							detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
							ElementMatrix_TSP(dNdx, detJ, EK);
							//element force
							Para2Phys_Dual(u[0], u[1], bzmesh[e], Nx, x);
							Fb = f_source(x[0], x[1]);
							ElementForce_TSP(Nx, detJ, Fb, EF);
						}
					}
					//sub-region 3
					for (i = 0; i<Gpt.size(); i++)//u
					{
						for (j = 0; j<Gpt.size(); j++)//v
						{
							u[0] = Gpt[i] / pow2;
							u[1] = (Gpt[j] + 1.) / pow2;
							BasisFunction_Dual(u[0], u[1], bzmesh[e], Nx, dNdx, detJ);
							detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
							ElementMatrix_TSP(dNdx, detJ, EK);
							//element force
							Para2Phys_Dual(u[0], u[1], bzmesh[e], Nx, x);
							Fb = f_source(x[0], x[1]);
							ElementForce_TSP(Nx, detJ, Fb, EF);
						}
					}
				}
				//it==lev_intg, sub-region 0
				pow2 = pow(2., double(lev_intg));
				for (i = 0; i<Gpt.size(); i++)
				{
					for (j = 0; j<Gpt.size(); j++)
					{
						u[0] = Gpt[i] / pow2;
						u[1] = Gpt[j] / pow2;
						BasisFunction_Dual(u[0], u[1], bzmesh[e], Nx, dNdx, detJ);
						detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
						ElementMatrix_TSP(dNdx, detJ, EK);
						//element force
						Para2Phys_Dual(u[0], u[1], bzmesh[e], Nx, x);
						Fb = f_source(x[0], x[1]);
						ElementForce_TSP(Nx, detJ, Fb, EF);
					}
				}
			}
		}
		
		Assembly_TSP(EK, EF, bzmesh[e].IEN, GK, GF);
	}
}

void Laplace::Quantityh_Dual(const vector<int>& IEN, const vector<double>& Nx, const vector<array<double,2>>& dNdx, double val[3])
{
	val[0] = 0.; val[1] = 0.; val[2] = 0.;
	for (uint i = 0; i < IEN.size(); i++)
	{
		val[0] += Nx[i] * uh[IEN[i]];
		val[1] += dNdx[i][0] * uh[IEN[i]];
		val[2] += dNdx[i][1] * uh[IEN[i]];
	}
}

void Laplace::ElementErrorDual(const BezierElement& bzel, double& L2, double& H1)
{
	L2 = 0.; H1 = 0.;
	vector<double> Nx(bzel.IEN.size());
	vector<array<double, 2>> dNdx(bzel.IEN.size());
	double detJ, val[3], ue, uxy[2];
	array<double, 3> x;
	int lev_intg(15);
	if (bzel.type == 0)
	{
		for (unsigned int gid1 = 0; gid1<Gpt.size(); gid1++)
		{
			for (unsigned int gid2 = 0; gid2<Gpt.size(); gid2++)
			{
				BasisFunction_Dual(Gpt[gid1], Gpt[gid2], bzel, Nx, dNdx, detJ);
				detJ = wght[gid1] * wght[gid2] * detJ;
				Quantityh_Dual(bzel.IEN, Nx, dNdx, val);
				Para2Phys_Dual(Gpt[gid1], Gpt[gid2], bzel, Nx, x);
				ue = exact_sol(x[0], x[1]);
				grad_u(x[0], x[1], uxy);
				L2 += detJ*(val[0] - ue)*(val[0] - ue);
				H1 += detJ*((val[1] - uxy[0])*(val[1] - uxy[0]) + (val[2] - uxy[1])*(val[2] - uxy[1]));
			}
		}
	}
	else if (bzel.type == 4)
	{
		double pow2, u[2];
		uint i, j;
		for (int it = 1; it <= lev_intg; it++)
		{
			pow2 = pow(2., double(it));
			//sub-region 1
			for (i = 0; i<Gpt.size(); i++)//u
			{
				for (j = 0; j<Gpt.size(); j++)//v
				{
					u[0] = (Gpt[i] + 1.) / pow2;
					u[1] = Gpt[j] / pow2;
					BasisFunction_Dual(u[0], u[1], bzel, Nx, dNdx, detJ);
					detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
					Quantityh_Dual(bzel.IEN, Nx, dNdx, val);
					Para2Phys_Dual(u[0], u[1], bzel, Nx, x);
					ue = exact_sol(x[0], x[1]);
					grad_u(x[0], x[1], uxy);
					L2 += detJ*(val[0] - ue)*(val[0] - ue);
					H1 += detJ*((val[1] - uxy[0])*(val[1] - uxy[0]) + (val[2] - uxy[1])*(val[2] - uxy[1]));
				}
			}
			//sub-region 2
			for (i = 0; i<Gpt.size(); i++)//u
			{
				for (j = 0; j<Gpt.size(); j++)//v
				{
					u[0] = (Gpt[i] + 1.) / pow2;
					u[1] = (Gpt[j] + 1.) / pow2;
					BasisFunction_Dual(u[0], u[1], bzel, Nx, dNdx, detJ);
					detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
					Quantityh_Dual(bzel.IEN, Nx, dNdx, val);
					Para2Phys_Dual(u[0], u[1], bzel, Nx, x);
					ue = exact_sol(x[0], x[1]);
					grad_u(x[0], x[1], uxy);
					L2 += detJ*(val[0] - ue)*(val[0] - ue);
					H1 += detJ*((val[1] - uxy[0])*(val[1] - uxy[0]) + (val[2] - uxy[1])*(val[2] - uxy[1]));
				}
			}
			//sub-region 3
			for (i = 0; i<Gpt.size(); i++)//u
			{
				for (j = 0; j<Gpt.size(); j++)//v
				{
					u[0] = Gpt[i] / pow2;
					u[1] = (Gpt[j] + 1.) / pow2;
					BasisFunction_Dual(u[0], u[1], bzel, Nx, dNdx, detJ);
					detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
					Quantityh_Dual(bzel.IEN, Nx, dNdx, val);
					Para2Phys_Dual(u[0], u[1], bzel, Nx, x);
					ue = exact_sol(x[0], x[1]);
					grad_u(x[0], x[1], uxy);
					L2 += detJ*(val[0] - ue)*(val[0] - ue);
					H1 += detJ*((val[1] - uxy[0])*(val[1] - uxy[0]) + (val[2] - uxy[1])*(val[2] - uxy[1]));
				}
			}
		}
		//it==lev_intg, sub-region 0
		pow2 = pow(2., double(lev_intg));
		for (i = 0; i<Gpt.size(); i++)
		{
			for (j = 0; j<Gpt.size(); j++)
			{
				u[0] = Gpt[i] / pow2;
				u[1] = Gpt[j] / pow2;
				BasisFunction_Dual(u[0], u[1], bzel, Nx, dNdx, detJ);
				detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
				Quantityh_Dual(bzel.IEN, Nx, dNdx, val);
				Para2Phys_Dual(u[0], u[1], bzel, Nx, x);
				ue = exact_sol(x[0], x[1]);
				grad_u(x[0], x[1], uxy);
				L2 += detJ*(val[0] - ue)*(val[0] - ue);
				H1 += detJ*((val[1] - uxy[0])*(val[1] - uxy[0]) + (val[2] - uxy[1])*(val[2] - uxy[1]));
			}
		}
	}
}

void Laplace::VisualizeVTK_Dual(string fn, const vector<BezierElement>& bzmesh)
{
	vector<array<double, 3>> spt;//s means sample
	vector<array<int, 4>> sele;
	vector<array<double, 3>> lpt;//visulize parameter lines
	vector<array<int, 2>> led;//line connectivity
	vector<double> sdisp;
	vector<double> errL2;
	int ns(2), ecount(0);
	vector<double> su(ns);
	for (int i = 0; i<ns; i++)
	{
		su[i] = double(i) / (double(ns) - 1.);
	}

	double L2_all(0.), H1_all(0.);
	for (unsigned int e = 0; e<bzmesh.size(); e++)
	{
		int loc(0);
		double L2(0.), H1(0.);
		//if (bzmesh[e].focus == 1)
		{
			ElementErrorDual(bzmesh[e], L2, H1);
		}		
		L2_all += L2; 
		H1_all += H1;
		errL2.push_back(sqrt(L2));
		//errL2.push_back(sqrt(H1));
		vector<double> Nx(bzmesh[e].IEN.size());
		vector<array<double, 2>> dNdx(bzmesh[e].IEN.size());
		double detJ, val[3];
		for (int a = 0; a<ns; a++)
		{
			for (int b = 0; b<ns; b++)
			{
				array<double, 3> pt;
				BasisFunction_Dual(su[b], su[a], bzmesh[e], Nx, dNdx, detJ);
				Para2Phys_Dual(su[b], su[a], bzmesh[e], Nx, pt);
				Quantityh_Dual(bzmesh[e].IEN, Nx, dNdx, val);
				spt.push_back(pt);
				//val[0] = 0.;
				//for (uint i = 0; i < Nx.size(); i++) val[0] += Nx[i];
				sdisp.push_back(val[0]);
				if (a == 0 || a == ns - 1 || b == 0 || b == ns - 1)
				{
					lpt.push_back(pt);
				}
			}
		}
		for (int a = 0; a<ns - 1; a++)
		{
			for (int b = 0; b<ns - 1; b++)
			{
				array<int, 4> el;
				el[0] = ecount*ns*ns + a*ns + b;
				el[1] = ecount*ns*ns + a*ns + b + 1;
				el[2] = ecount*ns*ns + (a + 1)*ns + b + 1;
				el[3] = ecount*ns*ns + (a + 1)*ns + b;
				sele.push_back(el);
			}
		}
		for (int a = 0; a<ns - 1; a++)
		{
			array<int, 2> lc;
			lc[0] = ecount * 4 * (ns - 1) + a;
			lc[1] = ecount * 4 * (ns - 1) + a + 1;
			led.push_back(lc);
			lc[0] = ecount * 4 * (ns - 1) + 3 * ns - 4 + a;
			lc[1] = ecount * 4 * (ns - 1) + 3 * ns - 4 + a + 1;
			led.push_back(lc);
		}
		for (int a = 0; a<ns - 2; a++)
		{
			array<int, 2> lc;
			lc[0] = ecount * 4 * (ns - 1) + ns + 2 * a;
			lc[1] = ecount * 4 * (ns - 1) + ns + 2 * a + 2;
			led.push_back(lc);
			lc[0] = ecount * 4 * (ns - 1) + ns + 2 * a - 1;
			lc[1] = ecount * 4 * (ns - 1) + ns + 2 * a + 1;
			led.push_back(lc);
		}
		array<int, 2> lc1;
		lc1[0] = ecount * 4 * (ns - 1);
		lc1[1] = ecount * 4 * (ns - 1) + ns;
		led.push_back(lc1);
		lc1[0] = ecount * 4 * (ns - 1) + 3 * ns - 5;
		lc1[1] = ecount * 4 * (ns - 1) + 4 * ns - 5;
		led.push_back(lc1);
		ecount++;
	}

	//hmax
	double hmax(0.), ltmp;
	int imax(0);
	for (uint i = 0; i < sele.size(); i++)
	{
		if (bzmesh[i].focus == 1)
		{
			ltmp = sqrt((spt[sele[i][2]][0] - spt[sele[i][0]][0])*(spt[sele[i][2]][0] - spt[sele[i][0]][0]) +
				(spt[sele[i][2]][1] - spt[sele[i][0]][1])*(spt[sele[i][2]][1] - spt[sele[i][0]][1]));
			if (ltmp > hmax)
			{
				hmax = ltmp;
				imax = i;
			}
			ltmp = sqrt((spt[sele[i][3]][0] - spt[sele[i][1]][0])*(spt[sele[i][3]][0] - spt[sele[i][1]][0]) +
				(spt[sele[i][3]][1] - spt[sele[i][1]][1])*(spt[sele[i][3]][1] - spt[sele[i][1]][1]));
			if (ltmp > hmax)
			{
				hmax = ltmp;
				imax = i;
			}
		}
	}
	cout << "hmax: " << hmax << "\n";
	cout << "hmax etype: " << imax << " " << bzmesh[imax].type << "\n";

	cout << "total L2: " << sqrt(L2_all) << "\n";
	cout<<"total H1: "<<sqrt(H1_all)<<"\n";

	string fname = fn + "_sol.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "9\n";
		}
		//fout<<"POINT_DATA "<<sdisp.size()<<"\nVECTORS disp float\n";
		//for(i=0;i<sdisp.size();i++)
		//{
		//	fout<<sdisp[i][0]<<" "<<sdisp[i][1]<<" 0\n";
		//}
		//fout<<"POINT_DATA "<<sse.size()<<"\nVECTORS strain float\n";
		//for(i=0;i<sse.size();i++)
		//{
		//	fout<<sse[i][0]<<" "<<sse[i][1]<<" "<<sse[i][2]<<"\n";
		//}
		//fout<<"POINT_DATA "<<sss.size()<<"\nVECTORS stress float\n";
		//for(i=0;i<sss.size();i++)
		//{
		//	fout<<sss[i][0]<<" "<<sss[i][1]<<" "<<sss[i][2]<<"\n";
		//}

		fout<<"POINT_DATA "<<sdisp.size()<<"\nSCALARS u float 1\nLOOKUP_TABLE default\n";
		for(i=0;i<sdisp.size();i++)
		{
			fout<<sdisp[i]<<"\n";
		}

		//fout << "\nCELL_DATA " << (ns - 1)*(ns - 1)*errL2.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		//for (i = 0; i<errL2.size(); i++)
		//{
		//	for (int j = 0; j<(ns - 1)*(ns - 1); j++)
		//		fout << errL2[i] << "\n";
		//}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	/*string fname1(fn + "-lines.vtk");
	ofstream fout1;
	fout1.open(fname1.c_str());
	if (fout1.is_open())
	{
		fout1 << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout1 << "POINTS " << lpt.size() << " float\n";
		for (uint i = 0; i<lpt.size(); i++)
		{
			fout1 << lpt[i][0] << " " << lpt[i][1] << " " << lpt[i][2] << "\n";
		}
		fout1 << "\nCELLS " << led.size() << " " << 3 * led.size() << '\n';
		for (uint i = 0; i<led.size(); i++)
		{
			fout1 << "2 " << led[i][0] << " " << led[i][1] << '\n';
		}
		fout1 << "\nCELL_TYPES " << led.size() << '\n';
		for (uint i = 0; i<led.size(); i++)
		{
			fout1 << "3\n";
		}
		fout1.close();
	}
	else
	{
		cout << "Cannot open " << fname1 << "!\n";
	}*/

	string fname_err(fn + "_err.txt");
	fout.open(fname_err);
	if (fout.is_open())
	{
		fout << "npt: " << npt << "\n";
		fout << "neq: " << neq << "\n";
		fout << "total L2: " << sqrt(L2_all) << "\n";
		fout << "total H1: " << sqrt(H1_all) << "\n";
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname_err << "!\n";
	}
}

void Laplace::VisualizeBasisFunction(int pid, const vector<BezierElement>& bzmesh, string fn)
{
	vector<array<double, 3>> spt;//s means sample
	vector<array<int, 4>> sele;
	vector<array<double, 3>> lpt;
	vector<array<int, 2>> led;
	//vector<double> sval;
	vector<array<double, 3>> sval;
	int ns(5), ecount(0);
	vector<double> su(ns);
	for (int i = 0; i<ns; i++)
	{
		su[i] = double(i) / (double(ns) - 1.);
	}

	for (unsigned int e = 0; e<bzmesh.size(); e++)
	{
		vector<double> Nx(bzmesh[e].IEN.size());
		vector<array<double, 2>> dNdx(bzmesh[e].IEN.size());
		double detJ;
		int pos(-1);
		vector<int>::const_iterator it = find(bzmesh[e].IEN.begin(), bzmesh[e].IEN.end(), pid);
		if (it != bzmesh[e].IEN.end()) pos = it - bzmesh[e].IEN.begin();
		for (int a = 0; a<ns; a++)
		{
			for (int b = 0; b<ns; b++)
			{
				array<double, 3> pt;
				BasisFunction_Dual(su[b], su[a], bzmesh[e], Nx, dNdx, detJ);
				Para2Phys_Dual(su[b], su[a], bzmesh[e], Nx, pt);
				spt.push_back(pt);
				array<double, 3> bfval = { 0.,0.,0. };
				if (pos != -1)
				{
					bfval[0] = Nx[pos];
					bfval[1] = dNdx[pos][0];
					bfval[2] = dNdx[pos][1];
				}
				sval.push_back(bfval);
				if (a == 0 || a == ns - 1 || b == 0 || b == ns - 1)
				{
					lpt.push_back(pt);
				}
			}
		}
		for (int a = 0; a<ns - 1; a++)
		{
			for (int b = 0; b<ns - 1; b++)
			{
				array<int, 4> el;
				el[0] = ecount*ns*ns + a*ns + b;
				el[1] = ecount*ns*ns + a*ns + b + 1;
				el[2] = ecount*ns*ns + (a + 1)*ns + b + 1;
				el[3] = ecount*ns*ns + (a + 1)*ns + b;
				sele.push_back(el);
			}
		}
		for (int a = 0; a<ns - 1; a++)
		{
			array<int, 2> lc;
			lc[0] = ecount * 4 * (ns - 1) + a;
			lc[1] = ecount * 4 * (ns - 1) + a + 1;
			led.push_back(lc);
			lc[0] = ecount * 4 * (ns - 1) + 3 * ns - 4 + a;
			lc[1] = ecount * 4 * (ns - 1) + 3 * ns - 4 + a + 1;
			led.push_back(lc);
		}
		for (int a = 0; a<ns - 2; a++)
		{
			array<int, 2> lc;
			lc[0] = ecount * 4 * (ns - 1) + ns + 2 * a;
			lc[1] = ecount * 4 * (ns - 1) + ns + 2 * a + 2;
			led.push_back(lc);
			lc[0] = ecount * 4 * (ns - 1) + ns + 2 * a - 1;
			lc[1] = ecount * 4 * (ns - 1) + ns + 2 * a + 1;
			led.push_back(lc);
		}
		array<int, 2> lc1;
		lc1[0] = ecount * 4 * (ns - 1);
		lc1[1] = ecount * 4 * (ns - 1) + ns;
		led.push_back(lc1);
		lc1[0] = ecount * 4 * (ns - 1) + 3 * ns - 5;
		lc1[1] = ecount * 4 * (ns - 1) + 4 * ns - 5;
		led.push_back(lc1);
		ecount++;
	}

	string fname = fn + to_string(pid) + "_bfval.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "9\n";
		}
		fout<<"POINT_DATA "<<sval.size()<<"\nVECTORS disp float\n";
		for(i=0;i<sval.size();i++)
		{
			fout<<sval[i][0]<<" "<<sval[i][1] << " " << sval[i][2] <<"\n";
		}
		//fout<<"POINT_DATA "<<sse.size()<<"\nVECTORS strain float\n";
		//for(i=0;i<sse.size();i++)
		//{
		//	fout<<sse[i][0]<<" "<<sse[i][1]<<" "<<sse[i][2]<<"\n";
		//}
		//fout<<"POINT_DATA "<<sss.size()<<"\nVECTORS stress float\n";
		//for(i=0;i<sss.size();i++)
		//{
		//	fout<<sss[i][0]<<" "<<sss[i][1]<<" "<<sss[i][2]<<"\n";
		//}
		//fout << "POINT_DATA " << sdisp.size() << "\nSCALARS u float 1\nLOOKUP_TABLE default\n";
		//for (i = 0; i<sdisp.size(); i++)
		//{
		//	fout << sdisp[i] << "\n";
		//}
		//fout << "\nCELL_DATA " << (ns - 1)*(ns - 1)*errL2.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		//for (i = 0; i<errL2.size(); i++)
		//{
		//	for (int j = 0; j<(ns - 1)*(ns - 1); j++)
		//		fout << errL2[i] << "\n";
		//}
		//fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	string fname1(fname + "-lines.vtk");
	ofstream fout1;
	fout1.open(fname1.c_str());
	if (fout1.is_open())
	{
		fout1 << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout1 << "POINTS " << lpt.size() << " float\n";
		for (uint i = 0; i<lpt.size(); i++)
		{
			fout1 << lpt[i][0] << " " << lpt[i][1] << " " << lpt[i][2] << "\n";
		}
		fout1 << "\nCELLS " << led.size() << " " << 3 * led.size() << '\n';
		for (uint i = 0; i<led.size(); i++)
		{
			fout1 << "2 " << led[i][0] << " " << led[i][1] << '\n';
		}
		fout1 << "\nCELL_TYPES " << led.size() << '\n';
		for (uint i = 0; i<led.size(); i++)
		{
			fout1 << "3\n";
		}
		fout1.close();
	}
	else
	{
		cout << "Cannot open " << fname1 << "!\n";
	}
}

void Laplace::Run_Dual(const vector<BezierElement>& bzmesh, string fn)
{
	//int pid(115);
	//VisualizeBasisFunction(pid, bzmesh, fn);

	GaussInfo(4);
	SparseMatrix<double> GK(neq, neq);
	VectorXd GF(neq);
	GK.setZero();
	GF.setZero();
	cout << "Building linear system...\n";
	BuildLinearSystem_Dual(bzmesh, GK, GF);
	Solver(GK, GF);
	VisualizeVTK_Dual(fn, bzmesh);
}





//////////////////////////////////////////////////////////////////////////////////////

void Laplace::BasisFunction_CC(double u, double v, const BezierElement& bzel,
	vector<double>& Nx, vector<array<double, 2>>& dNdx, double& detJ)
{
	//Nx.clear();
	//dNdx.clear();
	unsigned int i, j;
	if (bzel.type == 0)
	{
		double Nx0[16];
		double dNdx0[16][2];
		BasisFunction(u, v, bzel.pts, Nx0, dNdx0, detJ);
		//Nx.resize(bzel.cmat.size(), 0);
		//dNdx.resize(bzel.cmat.size());
		for (i = 0; i<bzel.cmat.size(); i++)
		{
			Nx[i] = 0.;
			dNdx[i][0] = 0.; dNdx[i][1] = 0.;
			for (j = 0; j<bzel.cmat[i].size(); j++)
			{
				Nx[i] += bzel.cmat[i][j] * Nx0[j];
				dNdx[i][0] += bzel.cmat[i][j] * dNdx0[j][0];
				dNdx[i][1] += bzel.cmat[i][j] * dNdx0[j][1];
			}
			//cout << dNdx[i][0] << " " << dNdx[i][1] << "\n";
		}
		//getchar();
	}
	else if (bzel.type == 4)
	{
		BasisFunctionCC_Irr(u, v, bzel, Nx, dNdx, detJ);
	}
}

void Laplace::BasisFunctionCC_Irr(double u, double v, const BezierElement& bzel,
	vector<double>& Nx, vector<array<double, 2>>& dNdx, double& detJ)
{
	double tol(1.e-10), eps(1.e-8);
	double u1 = max(tol, u);
	double v1 = max(tol, v);
	int n = floor(min(-log2(u1), -log2(v1))) + 1, subid;
	double pow2 = pow(2., double(n - 1));
	u1 *= pow2; v1 *= pow2;
	if (v1<0.5)
	{
		subid = 0; u1 = 2.*u1 - 1.; v1 = 2.*v1;
	}
	else if (u1<0.5)
	{
		subid = 2; u1 = 2.*u1; v1 = 2.*v1 - 1.;
	}
	else
	{
		subid = 1; u1 = 2.*u1 - 1.; v1 = 2.*v1 - 1.;
	}
	int nv((bzel.IEN.size() - 8)/2);
	RegularPatchBasis Nu, Nv;
	Nu.Evaluate(u1, 1., 0.);
	Nv.Evaluate(v1, 1., 0.);
	double Nt1[16];
	double dNdt1[16][2];
	int loc(0);
	pow2 *= 2.;
	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			Nt1[loc] = Nu.val[i] * Nv.val[j];
			dNdt1[loc][0] = Nu.Dval[i] * Nv.val[j];
			dNdt1[loc][1] = Nu.val[i] * Nv.Dval[j];
			loc++;
		}
	}
	int i7(7), nv2(2 * nv);
	if (nv == 3) i7 = 1;
	int Pk[3][16] = { { i7,6,nv2 + 1,nv2 + 8,0,5,nv2 + 2,nv2 + 9,3,4,nv2 + 3,nv2 + 10,nv2 + 6,nv2 + 5,nv2 + 4,nv2 + 11 },
	{ 0,5,nv2 + 2,nv2 + 9,3,4,nv2 + 3,nv2 + 10,nv2 + 6,nv2 + 5,nv2 + 4,nv2 + 11,nv2 + 15,nv2 + 14,nv2 + 13,nv2 + 12 },
	{ 1,0,5,nv2 + 2,2,3,4,nv2 + 3,nv2 + 7,nv2 + 6,nv2 + 5,nv2 + 4,nv2 + 16,nv2 + 15,nv2 + 14,nv2 + 13 } };
	MatrixXd mat(bzel.IEN.size(), 16);
	for (uint i = 0; i < bzel.IEN.size(); i++)
	{
		for (int j = 0; j < 16; j++)
		{
			mat(i, j) = bzel.abar(i, Pk[subid][j]);
		}
	}
	for (int i = 0; i < n - 1; i++)
	{
		mat = bzel.amat * mat;
	}
	vector<array<double, 2>> dNdt(bzel.IEN.size());
	double dxdt[2][2] = { { 0.,0. },{ 0.,0. } };
	for (uint i = 0; i < bzel.IEN.size(); i++)
	{
		Nx[i] = 0.;
		dNdt[i][0] = 0.; dNdt[i][1] = 0.;
		for (int j = 0; j < 16; j++)
		{
			Nx[i] += mat(i, j)*Nt1[j];
			dNdt[i][0] += mat(i, j)*dNdt1[j][0];
			dNdt[i][1] += mat(i, j)*dNdt1[j][1];
		}
		dNdt[i][0] *= pow2;
		dNdt[i][1] *= pow2;
		for (int a = 0; a < 2; a++)
		{
			for (int b = 0; b < 2; b++)
			{
				dxdt[a][b] += bzel.cp[i][a] * dNdt[i][b];
			}
		}
	}
	detJ = dxdt[0][0] * dxdt[1][1] - dxdt[0][1] * dxdt[1][0];
	double dtdx[2][2] = { { dxdt[1][1] / detJ,-dxdt[0][1] / detJ },{ -dxdt[1][0] / detJ,dxdt[0][0] / detJ } };
	for (uint i = 0; i<bzel.IEN.size(); i++)
	{
		dNdx[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0];
		dNdx[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1];
	}
	detJ = 0.25*detJ;
}

void Laplace::BuildLinearSystem_CC(const vector<BezierElement>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF)
{
	InitializeSparseMat(bzmesh, GK);
	unsigned int e, i, j;
	int lev_intg(15);//levels of subdivision for integration, 0 means no subdivision
	for (e = 0; e<bzmesh.size(); e++)
	{
		if (e != 0 && e % 100 == 0)
		{
			cout << e << " ";
		}
		vector<vector<double>> EK(bzmesh[e].IEN.size(), vector<double>(bzmesh[e].IEN.size(), 0.));
		vector<double> Nx(bzmesh[e].IEN.size()), EF(bzmesh[e].IEN.size(), 0.);
		vector<array<double, 2>> dNdx(bzmesh[e].IEN.size());
		double detJ, Fb;
		array<double, 3> x;
		//cout << "eid: " << e << "\n";
		if (bzmesh[e].type == 0)
		{
			for (i = 0; i<Gpt.size(); i++)
			{
				for (j = 0; j<Gpt.size(); j++)
				{
					BasisFunction_CC(Gpt[i], Gpt[j], bzmesh[e], Nx, dNdx, detJ);
					detJ = wght[i] * wght[j] * detJ;
					ElementMatrix_TSP(dNdx, detJ, EK);
					//element force
					Para2Phys_Dual(Gpt[i], Gpt[j], bzmesh[e], Nx, x);
					Fb = f_source(x[0], x[1]);
					ElementForce_TSP(Nx, detJ, Fb, EF);
				}
			}
		}
		else if (bzmesh[e].type == 4)
		{
			//if(lev_intg==0)
			//{
			//	for (i = 0; i<Gpt.size(); i++)
			//	{
			//		for (j = 0; j<Gpt.size(); j++)
			//		{
			//			BasisFunction_Dual(Gpt[i], Gpt[j], bzmesh[e], Nx, dNdx, detJ);
			//			detJ = wght[i] * wght[j] * detJ;
			//			ElementMatrix_TSP(dNdx, detJ, EK);
			//			//element force
			//			Para2Phys_Dual(Gpt[i], Gpt[j], bzmesh[e], Nx, x);
			//			Fb = f_source(x[0], x[1]);
			//			ElementForce_TSP(Nx, detJ, Fb, EF);
			//		}
			//	}
			//}
			//else//greater than 1
			{
				double pow2, u[2];
				for (int it = 1; it <= lev_intg; it++)
				{
					pow2 = pow(2., double(it));
					//sub-region 1
					for (i = 0; i<Gpt.size(); i++)//u
					{
						for (j = 0; j<Gpt.size(); j++)//v
						{
							u[0] = (Gpt[i] + 1.) / pow2;
							u[1] = Gpt[j] / pow2;
							BasisFunction_CC(u[0], u[1], bzmesh[e], Nx, dNdx, detJ);
							detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
							ElementMatrix_TSP(dNdx, detJ, EK);
							//element force
							Para2Phys_Dual(u[0], u[1], bzmesh[e], Nx, x);
							Fb = f_source(x[0], x[1]);
							ElementForce_TSP(Nx, detJ, Fb, EF);
						}
					}
					//sub-region 2
					for (i = 0; i<Gpt.size(); i++)//u
					{
						for (j = 0; j<Gpt.size(); j++)//v
						{
							u[0] = (Gpt[i] + 1.) / pow2;
							u[1] = (Gpt[j] + 1.) / pow2;
							BasisFunction_CC(u[0], u[1], bzmesh[e], Nx, dNdx, detJ);
							detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
							ElementMatrix_TSP(dNdx, detJ, EK);
							//element force
							Para2Phys_Dual(u[0], u[1], bzmesh[e], Nx, x);
							Fb = f_source(x[0], x[1]);
							ElementForce_TSP(Nx, detJ, Fb, EF);
						}
					}
					//sub-region 3
					for (i = 0; i<Gpt.size(); i++)//u
					{
						for (j = 0; j<Gpt.size(); j++)//v
						{
							u[0] = Gpt[i] / pow2;
							u[1] = (Gpt[j] + 1.) / pow2;
							BasisFunction_CC(u[0], u[1], bzmesh[e], Nx, dNdx, detJ);
							detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
							ElementMatrix_TSP(dNdx, detJ, EK);
							//element force
							Para2Phys_Dual(u[0], u[1], bzmesh[e], Nx, x);
							Fb = f_source(x[0], x[1]);
							ElementForce_TSP(Nx, detJ, Fb, EF);
						}
					}
				}
				//it==lev_intg, sub-region 0
				pow2 = pow(2., double(lev_intg));
				for (i = 0; i<Gpt.size(); i++)
				{
					for (j = 0; j<Gpt.size(); j++)
					{
						u[0] = Gpt[i] / pow2;
						u[1] = Gpt[j] / pow2;
						BasisFunction_CC(u[0], u[1], bzmesh[e], Nx, dNdx, detJ);
						detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
						ElementMatrix_TSP(dNdx, detJ, EK);
						//element force
						Para2Phys_Dual(u[0], u[1], bzmesh[e], Nx, x);
						Fb = f_source(x[0], x[1]);
						ElementForce_TSP(Nx, detJ, Fb, EF);
					}
				}
			}
		}

		Assembly_TSP(EK, EF, bzmesh[e].IEN, GK, GF);
	}
}

void Laplace::ElementErrorCC(const BezierElement& bzel, double& L2, double& H1)
{
	L2 = 0.; H1 = 0.;
	vector<double> Nx(bzel.IEN.size());
	vector<array<double, 2>> dNdx(bzel.IEN.size());
	double detJ, val[3], ue, uxy[2];
	array<double, 3> x;
	int lev_intg(15);
	if (bzel.type == 0)
	{
		for (unsigned int gid1 = 0; gid1<Gpt.size(); gid1++)
		{
			for (unsigned int gid2 = 0; gid2<Gpt.size(); gid2++)
			{
				BasisFunction_CC(Gpt[gid1], Gpt[gid2], bzel, Nx, dNdx, detJ);
				detJ = wght[gid1] * wght[gid2] * detJ;
				Quantityh_Dual(bzel.IEN, Nx, dNdx, val);
				Para2Phys_Dual(Gpt[gid1], Gpt[gid2], bzel, Nx, x);
				ue = exact_sol(x[0], x[1]);
				grad_u(x[0], x[1], uxy);
				L2 += detJ*(val[0] - ue)*(val[0] - ue);
				H1 += detJ*((val[1] - uxy[0])*(val[1] - uxy[0]) + (val[2] - uxy[1])*(val[2] - uxy[1]));
			}
		}
	}
	else if (bzel.type == 4)
	{
		double pow2, u[2];
		uint i, j;
		for (int it = 1; it <= lev_intg; it++)
		{
			pow2 = pow(2., double(it));
			//sub-region 1
			for (i = 0; i<Gpt.size(); i++)//u
			{
				for (j = 0; j<Gpt.size(); j++)//v
				{
					u[0] = (Gpt[i] + 1.) / pow2;
					u[1] = Gpt[j] / pow2;
					BasisFunction_CC(u[0], u[1], bzel, Nx, dNdx, detJ);
					detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
					Quantityh_Dual(bzel.IEN, Nx, dNdx, val);
					Para2Phys_Dual(u[0], u[1], bzel, Nx, x);
					ue = exact_sol(x[0], x[1]);
					grad_u(x[0], x[1], uxy);
					L2 += detJ*(val[0] - ue)*(val[0] - ue);
					H1 += detJ*((val[1] - uxy[0])*(val[1] - uxy[0]) + (val[2] - uxy[1])*(val[2] - uxy[1]));
				}
			}
			//sub-region 2
			for (i = 0; i<Gpt.size(); i++)//u
			{
				for (j = 0; j<Gpt.size(); j++)//v
				{
					u[0] = (Gpt[i] + 1.) / pow2;
					u[1] = (Gpt[j] + 1.) / pow2;
					BasisFunction_CC(u[0], u[1], bzel, Nx, dNdx, detJ);
					detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
					Quantityh_Dual(bzel.IEN, Nx, dNdx, val);
					Para2Phys_Dual(u[0], u[1], bzel, Nx, x);
					ue = exact_sol(x[0], x[1]);
					grad_u(x[0], x[1], uxy);
					L2 += detJ*(val[0] - ue)*(val[0] - ue);
					H1 += detJ*((val[1] - uxy[0])*(val[1] - uxy[0]) + (val[2] - uxy[1])*(val[2] - uxy[1]));
				}
			}
			//sub-region 3
			for (i = 0; i<Gpt.size(); i++)//u
			{
				for (j = 0; j<Gpt.size(); j++)//v
				{
					u[0] = Gpt[i] / pow2;
					u[1] = (Gpt[j] + 1.) / pow2;
					BasisFunction_CC(u[0], u[1], bzel, Nx, dNdx, detJ);
					detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
					Quantityh_Dual(bzel.IEN, Nx, dNdx, val);
					Para2Phys_Dual(u[0], u[1], bzel, Nx, x);
					ue = exact_sol(x[0], x[1]);
					grad_u(x[0], x[1], uxy);
					L2 += detJ*(val[0] - ue)*(val[0] - ue);
					H1 += detJ*((val[1] - uxy[0])*(val[1] - uxy[0]) + (val[2] - uxy[1])*(val[2] - uxy[1]));
				}
			}
		}
		//it==lev_intg, sub-region 0
		pow2 = pow(2., double(lev_intg));
		for (i = 0; i<Gpt.size(); i++)
		{
			for (j = 0; j<Gpt.size(); j++)
			{
				u[0] = Gpt[i] / pow2;
				u[1] = Gpt[j] / pow2;
				BasisFunction_CC(u[0], u[1], bzel, Nx, dNdx, detJ);
				detJ = wght[i] * wght[j] * detJ / (pow2 * pow2);
				Quantityh_Dual(bzel.IEN, Nx, dNdx, val);
				Para2Phys_Dual(u[0], u[1], bzel, Nx, x);
				ue = exact_sol(x[0], x[1]);
				grad_u(x[0], x[1], uxy);
				L2 += detJ*(val[0] - ue)*(val[0] - ue);
				H1 += detJ*((val[1] - uxy[0])*(val[1] - uxy[0]) + (val[2] - uxy[1])*(val[2] - uxy[1]));
			}
		}
	}
}

void Laplace::VisualizeVTK_CC(string fn, const vector<BezierElement>& bzmesh)
{
	vector<array<double, 3>> spt;//s means sample
	vector<array<int, 4>> sele;
	vector<array<double, 3>> lpt;//visulize parameter lines
	vector<array<int, 2>> led;//line connectivity
	vector<double> sdisp;
	vector<double> errL2;
	int ns(2), ecount(0);
	vector<double> su(ns);
	for (int i = 0; i<ns; i++)
	{
		su[i] = double(i) / (double(ns) - 1.);
	}

	double L2_all(0.), H1_all(0.);
	for (unsigned int e = 0; e<bzmesh.size(); e++)
	{
		int loc(0);
		double L2, H1;
		ElementErrorCC(bzmesh[e], L2, H1);
		L2_all += L2;
		H1_all += H1;
		errL2.push_back(sqrt(L2));
		//errL2.push_back(sqrt(H1));
		vector<double> Nx(bzmesh[e].IEN.size());
		vector<array<double, 2>> dNdx(bzmesh[e].IEN.size());
		double detJ, val[3];
		for (int a = 0; a<ns; a++)
		{
			for (int b = 0; b<ns; b++)
			{
				array<double, 3> pt;
				BasisFunction_CC(su[b], su[a], bzmesh[e], Nx, dNdx, detJ);
				Para2Phys_Dual(su[b], su[a], bzmesh[e], Nx, pt);
				Quantityh_Dual(bzmesh[e].IEN, Nx, dNdx, val);
				spt.push_back(pt);
				//val[0] = 0.;
				//for (uint i = 0; i < Nx.size(); i++) val[0] += Nx[i];
				sdisp.push_back(val[0]);
				if (a == 0 || a == ns - 1 || b == 0 || b == ns - 1)
				{
					lpt.push_back(pt);
				}
			}
		}
		for (int a = 0; a<ns - 1; a++)
		{
			for (int b = 0; b<ns - 1; b++)
			{
				array<int, 4> el;
				el[0] = ecount*ns*ns + a*ns + b;
				el[1] = ecount*ns*ns + a*ns + b + 1;
				el[2] = ecount*ns*ns + (a + 1)*ns + b + 1;
				el[3] = ecount*ns*ns + (a + 1)*ns + b;
				sele.push_back(el);
			}
		}
		for (int a = 0; a<ns - 1; a++)
		{
			array<int, 2> lc;
			lc[0] = ecount * 4 * (ns - 1) + a;
			lc[1] = ecount * 4 * (ns - 1) + a + 1;
			led.push_back(lc);
			lc[0] = ecount * 4 * (ns - 1) + 3 * ns - 4 + a;
			lc[1] = ecount * 4 * (ns - 1) + 3 * ns - 4 + a + 1;
			led.push_back(lc);
		}
		for (int a = 0; a<ns - 2; a++)
		{
			array<int, 2> lc;
			lc[0] = ecount * 4 * (ns - 1) + ns + 2 * a;
			lc[1] = ecount * 4 * (ns - 1) + ns + 2 * a + 2;
			led.push_back(lc);
			lc[0] = ecount * 4 * (ns - 1) + ns + 2 * a - 1;
			lc[1] = ecount * 4 * (ns - 1) + ns + 2 * a + 1;
			led.push_back(lc);
		}
		array<int, 2> lc1;
		lc1[0] = ecount * 4 * (ns - 1);
		lc1[1] = ecount * 4 * (ns - 1) + ns;
		led.push_back(lc1);
		lc1[0] = ecount * 4 * (ns - 1) + 3 * ns - 5;
		lc1[1] = ecount * 4 * (ns - 1) + 4 * ns - 5;
		led.push_back(lc1);
		ecount++;
	}

	//hmax
	double hmax(0.), ltmp;
	int imax(0);
	for (uint i = 0; i < sele.size(); i++)
	{
		ltmp = sqrt((spt[sele[i][2]][0] - spt[sele[i][0]][0])*(spt[sele[i][2]][0] - spt[sele[i][0]][0]) +
			(spt[sele[i][2]][1] - spt[sele[i][0]][1])*(spt[sele[i][2]][1] - spt[sele[i][0]][1]));
		if (ltmp > hmax)
		{
			hmax = ltmp;
			imax = i;
		}
		ltmp = sqrt((spt[sele[i][3]][0] - spt[sele[i][1]][0])*(spt[sele[i][3]][0] - spt[sele[i][1]][0]) +
			(spt[sele[i][3]][1] - spt[sele[i][1]][1])*(spt[sele[i][3]][1] - spt[sele[i][1]][1]));
		if (ltmp > hmax)
		{
			hmax = ltmp;
			imax = i;
		}
	}
	cout << "hmax: " << hmax << "\n";
	cout << "hmax etype: " << imax << " " << bzmesh[imax].type << "\n";

	cout << "total L2: " << sqrt(L2_all) << "\n";
	cout << "total H1: " << sqrt(H1_all) << "\n";

	string fname = fn + "_sol.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "9\n";
		}
		//fout<<"POINT_DATA "<<sdisp.size()<<"\nVECTORS disp float\n";
		//for(i=0;i<sdisp.size();i++)
		//{
		//	fout<<sdisp[i][0]<<" "<<sdisp[i][1]<<" 0\n";
		//}
		//fout<<"POINT_DATA "<<sse.size()<<"\nVECTORS strain float\n";
		//for(i=0;i<sse.size();i++)
		//{
		//	fout<<sse[i][0]<<" "<<sse[i][1]<<" "<<sse[i][2]<<"\n";
		//}
		//fout<<"POINT_DATA "<<sss.size()<<"\nVECTORS stress float\n";
		//for(i=0;i<sss.size();i++)
		//{
		//	fout<<sss[i][0]<<" "<<sss[i][1]<<" "<<sss[i][2]<<"\n";
		//}

		//fout<<"POINT_DATA "<<sdisp.size()<<"\nSCALARS u float 1\nLOOKUP_TABLE default\n";
		//for(i=0;i<sdisp.size();i++)
		//{
		//	fout<<sdisp[i]<<"\n";
		//}

		fout << "\nCELL_DATA " << (ns - 1)*(ns - 1)*errL2.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		for (i = 0; i<errL2.size(); i++)
		{
			for (int j = 0; j<(ns - 1)*(ns - 1); j++)
				fout << errL2[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	/*string fname1(fn + "-lines.vtk");
	ofstream fout1;
	fout1.open(fname1.c_str());
	if (fout1.is_open())
	{
	fout1 << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	fout1 << "POINTS " << lpt.size() << " float\n";
	for (uint i = 0; i<lpt.size(); i++)
	{
	fout1 << lpt[i][0] << " " << lpt[i][1] << " " << lpt[i][2] << "\n";
	}
	fout1 << "\nCELLS " << led.size() << " " << 3 * led.size() << '\n';
	for (uint i = 0; i<led.size(); i++)
	{
	fout1 << "2 " << led[i][0] << " " << led[i][1] << '\n';
	}
	fout1 << "\nCELL_TYPES " << led.size() << '\n';
	for (uint i = 0; i<led.size(); i++)
	{
	fout1 << "3\n";
	}
	fout1.close();
	}
	else
	{
	cout << "Cannot open " << fname1 << "!\n";
	}*/

	string fname_err(fn + "_err.txt");
	fout.open(fname_err);
	if (fout.is_open())
	{
		fout << "npt: " << npt << "\n";
		fout << "neq: " << neq << "\n";
		fout << "total L2: " << sqrt(L2_all) << "\n";
		fout << "total H1: " << sqrt(H1_all) << "\n";
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname_err << "!\n";
	}
}

void Laplace::Run_CC(const vector<BezierElement>& bzmesh, string fn)
{
	//int pid(115);
	//VisualizeBasisFunction(pid, bzmesh, fn);

	GaussInfo(4);
	SparseMatrix<double> GK(neq, neq);
	VectorXd GF(neq);
	GK.setZero();
	GF.setZero();
	cout << "Building linear system...\n";
	BuildLinearSystem_CC(bzmesh, GK, GF);
	Solver(GK, GF);
	VisualizeVTK_CC(fn, bzmesh);
}






///////////////////////////////////////////////////////////////////////////////////////////

double LDomainSolution(double x,double y)
{
	//const double PI(3.14159265359);
	double r=sqrt(x*x+y*y);
	double phi=atan2(y,x);
	if(phi<=0.) phi+=2.*PI;
	double u=pow(r,0.6666666667)*sin(2.*phi/3.-PI/3.);
	return u;
}

void LDomainSolutionDeriv(double x, double y, double& u, double& ux, double& uy)
{
	//const double PI(3.14159265359);
	double r=sqrt(x*x+y*y);
	double phi=atan2(y,x);
	if(phi<=0.) phi+=2.*PI;
	u=pow(r,0.6666666667)*sin(2.*phi/3.-PI/3.);
	ux=0.6666666667*pow(r,-0.3333333333)*(sin(2.*phi/3.-PI/3.)*cos(phi)-cos(2.*phi/3.-PI/3.)*sin(phi));
	uy=0.6666666667*pow(r,-0.3333333333)*(sin(2.*phi/3.-PI/3.)*sin(phi)+cos(2.*phi/3.-PI/3.)*cos(phi));
}

double ptestHO_sol_1(double x, double y)
{
	return (x+y);
}

double ptestHO_sol_2(double x, double y)
{
	return (x*x+y*y);
	//return (x*x);
}

double ptestHO_sol_3(double x, double y)
{
	//return (x*x*x+y*y*y);
	double u=x*(1.-x)*y*(1.-y);
	return u;
}



double exact_sol(double x, double y)
{
	//double u = exact_sol_1(x, y);
	//double u = exact_sol_4(x, y);
	//double u = exact_sol_5(x, y);
	//double u = exact_sol_6(x, y);
	double u = exact_sol_7(x, y);
	return u;
}

double exact_sol_1(double x, double y)
{
	double tmp = x + y;
	return tmp;
}

double exact_sol_4(double x, double y)
{
	double tmp = x * (1. - x)* y * (1. - y);
	return tmp;
}

double exact_sol_5(double x, double y)
{
	return x;
}

double exact_sol_6(double x, double y)
{
	return sin(PI*x) * sin(PI*y);
}

double exact_sol_7(double x, double y)
{
	return exp((x + y) / 2.);
}


void grad_u(double x, double y, double gu[2])
{
	//grad_u_1(x, y, gu);
	//grad_u_4(x, y, gu);
	//grad_u_6(x, y, gu);
	grad_u_7(x, y, gu);
}

void grad_u_1(double x, double y, double gu[2])
{
	gu[0] = 1.;
	gu[1] = 1.;
}

void grad_u_4(double x, double y, double gu[2])
{
	gu[0] = (1. - 2.*x)*(y - y*y);
	gu[1] = (1. - 2.*y)*(x - x*x);
}

void grad_u_6(double x, double y, double gu[2])
{
	gu[0] = PI*cos(PI*x)*sin(PI*y);
	gu[1] = PI*sin(PI*x)*cos(PI*y);
}

void grad_u_7(double x, double y, double gu[2])
{
	gu[0] = .5*exp((x + y) / 2.);
	gu[1] = .5*exp((x + y) / 2.);
}



double f_source(double x, double y)
{
	//double f = f_source_1(x, y);
	//double f = f_source_4(x, y);
	//double f = f_source_5(x, y);
	//double f = f_source_6(x, y);
	double f = f_source_7(x, y);
	return f;
}

double f_source_1(double x, double y)
{
	return 0.;
}

double f_source_2(double x, double y)
{
	return -4.;
	//return -2.;
}

double f_source_3(double x, double y)
{
	//return -6.*(x+y);
	double f=x+y-x*x-y*y;
	return f;
}

double f_source_4(double x, double y)
{
	double tmp = 2.*(x - x*x + y - y*y);
	return tmp;
}

double f_source_5(double x, double y)
{
	return 0.;
}

double f_source_6(double x, double y)
{
	double tmp = 2.*PI*PI*(sin(PI*x) * sin(PI*y));
	return tmp;
}

double f_source_7(double x, double y)
{
	double tmp = -0.5*exp((x + y) / 2.);
	return tmp;
}