#include "Elasiticity_3D.h"
#include <fstream>
#include <iostream>

using namespace std;

LinearElasticity::LinearElasticity()
{
	dim = 3;
}

void LinearElasticity::GaussInfo(int ng)
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
	case 5:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.046910077030668;			Gpt[1]=0.2307653449471585;			Gpt[2]=0.5;			Gpt[3]=0.7692346550528415;  Gpt[4]=0.953089922969332;
			wght[0]=0.2369268850561891;			wght[1]=0.4786286704993665;			wght[2]=0.5688888888888889;			wght[3]=0.4786286704993665; wght[4]=0.2369268850561891;
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

void LinearElasticity::SetMaterialProp(double Ey,double v)
{
	lambda=v*Ey/((1.+v)*(1.-2.*v));//plain strain, or volume
	mu=Ey/(2.*(1.+v));
	//lambda=v*Ey/(1.-v*v);//plane stress
}

void LinearElasticity::BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, double Nx[64], double dNdx[64][3], double& detJ)
{
	double Nu[4]={(1.-u)*(1.-u)*(1.-u),3.*(1.-u)*(1.-u)*u,3.*(1.-u)*u*u,u*u*u};
	double Nv[4]={(1.-v)*(1.-v)*(1.-v),3.*(1.-v)*(1.-v)*v,3.*(1.-v)*v*v,v*v*v};
	double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
	double dNdu[4]={-3.*(1.-u)*(1.-u),3.-12.*u+9.*u*u,3.*(2.-3.*u)*u,3.*u*u};
	double dNdv[4]={-3.*(1.-v)*(1.-v),3.-12.*v+9.*v*v,3.*(2.-3.*v)*v,3.*v*v};
	double dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
	double dNdt[64][3];
	int i,j,k,a,b,loc(0);
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			for (k = 0; k < 4; k++)
			{
				Nx[loc] = Nu[k] * Nv[j]*Nw[i];
				dNdt[loc][0] = dNdu[k] * Nv[j]*Nw[i];
				dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
				dNdt[loc][2] = Nu[k] * Nv[j]*dNdw[i];
				loc++;
			}
		}
	}
	Matrix3d dxdt=Matrix3d::Zero();
	loc=0;
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			for (k = 0; k < 4; k++)
			{
				for (a = 0; a<3; a++)
				{
					for (b = 0; b<3; b++)
					{
						dxdt(a,b) += pt[loc][a] * dNdt[loc][b];
					}
				}
				loc++;
			}
		}
	}
	//cout << dxdt << "\n\n";
	Matrix3d dtdx = dxdt.inverse();
	for(i=0;i<64;i++)
	{
		dNdx[i][0] = dNdt[i][0] * dtdx(0, 0) + dNdt[i][1] * dtdx(1, 0) + dNdt[i][2] * dtdx(2, 0);
		dNdx[i][1] = dNdt[i][0] * dtdx(0, 1) + dNdt[i][1] * dtdx(1, 1) + dNdt[i][2] * dtdx(2, 1);
		dNdx[i][2] = dNdt[i][0] * dtdx(0, 2) + dNdt[i][1] * dtdx(1, 2) + dNdt[i][2] * dtdx(2, 2);
	}
	detJ = dxdt.determinant();
	detJ=0.25*detJ;
}

//void LinearElasticity::BasisFunction_Rational(double u, double v, const BezierElement3D& bzel, vector<double>& Rx, vector<array<double,2>>& dRdx, double& detJ)
//{
//	double Bt[16], dBdt[16][2], wb(0.), dwbdt[2]={0.,0.};
//	bzel.Basis(u,v,Bt,dBdt);
//	for(int i=0; i<16; i++)
//	{
//		wb+=Bt[i]*bzel.w[i];
//		dwbdt[0]+=dBdt[i][0]*bzel.w[i];
//		dwbdt[1]+=dBdt[i][1]*bzel.w[i];
//	}
//	Rx.clear();
//	dRdx.clear();
//	Rx.resize(bzel.IEN.size(),0.);
//	dRdx.resize(bzel.IEN.size());
//	vector<array<double,2>> dRdt(bzel.IEN.size());
//	for(unsigned int i=0; i<bzel.IEN.size(); i++)
//	{
//		dRdt[i][0]=0.; dRdt[i][1]=0.; 
//		for(int j=0; j<16; j++)
//		{
//			Rx[i]+=bzel.w_nurbs[i]*bzel.cmat[i][j]*Bt[j]/wb;
//			dRdt[i][0]+=bzel.w_nurbs[i]*bzel.cmat[i][j]*(dBdt[j][0]/wb-dwbdt[0]*Bt[j]/(wb*wb));
//			dRdt[i][1]+=bzel.w_nurbs[i]*bzel.cmat[i][j]*(dBdt[j][1]/wb-dwbdt[1]*Bt[j]/(wb*wb));
//		}
//	}
//	//double Nx[16], dNdt[16][2], sum(0.), sum_u(0.), sum_v(0.);
//	//for(int i=0; i<16; i++)
//	//{
//	//	Bt[i]*=bzel.w[i];
//	//	dBdt[i][0]*=bzel.w[i];
//	//	dBdt[i][1]*=bzel.w[i];
//	//	sum+=Bt[i];
//	//	sum_u+=dBdt[i][0];
//	//	sum_v+=dBdt[i][1];
//	//}
//	for(int i=0; i<16; i++)
//	{
//		dBdt[i][0]=bzel.w[i]*(dBdt[i][0]*wb-Bt[i]*dwbdt[0])/(wb*wb);
//		dBdt[i][1]=bzel.w[i]*(dBdt[i][1]*wb-Bt[i]*dwbdt[1])/(wb*wb);
//	}
//	int loc(0);
//	double dxdt[2][2]={{0.,0.},{0.,0.}};
//	for(int i=0;i<16;i++)
//	{
//		for (int a=0;a<2;a++)
//		{
//			for(int b=0;b<2;b++)
//			{
//				dxdt[a][b]+=bzel.pts[i][a]*dBdt[i][b];
//			}
//		}
//	}
//
//	detJ=dxdt[0][0]*dxdt[1][1]-dxdt[0][1]*dxdt[1][0];
//	double dtdx[2][2]={{dxdt[1][1]/detJ,-dxdt[0][1]/detJ},{-dxdt[1][0]/detJ,dxdt[0][0]/detJ}};
//	for(unsigned int i=0; i<bzel.IEN.size(); i++)
//	{
//		dRdx[i][0]=dRdt[i][0]*dtdx[0][0]+dRdt[i][1]*dtdx[1][0];
//		dRdx[i][1]=dRdt[i][0]*dtdx[0][1]+dRdt[i][1]*dtdx[1][1];
//	}
//	detJ=0.25*detJ;
//
//	//double Nu[4]={(1.-u)*(1.-u)*(1.-u),3.*(1.-u)*(1.-u)*u,3.*(1.-u)*u*u,u*u*u};
//	//double Nv[4]={(1.-v)*(1.-v)*(1.-v),3.*(1.-v)*(1.-v)*v,3.*(1.-v)*v*v,v*v*v};
//	//double dNdu[4]={-3.*(1.-u)*(1.-u),3.-12.*u+9.*u*u,3.*(2.-3.*u)*u,3.*u*u};
//	//double dNdv[4]={-3.*(1.-v)*(1.-v),3.-12.*v+9.*v*v,3.*(2.-3.*v)*v,3.*v*v};
//	//double dNdt[16][2];
//	//int i,j,a,b,loc(0);
//	//double sum(0.), sum_u(0.), sum_v(0.);
//	//for(i=0;i<4;i++)
//	//{
//	//	for(j=0;j<4;j++)
//	//	{
//	//		//Nx[loc]=Nu[j]*Nv[i]*w[loc];
//	//		//dNdt[loc][0]=dNdu[j]*Nv[i]*w[loc];
//	//		//dNdt[loc][1]=Nu[j]*dNdv[i]*w[loc];
//	//		//sum+=Nx[loc];
//	//		//sum_u+=dNdt[loc][0];
//	//		//sum_v+=dNdt[loc][1];
//	//		Nx[loc]=Nu[j]*Nv[i];
//	//		dNdt[loc][0]=dNdu[j]*Nv[i];
//	//		dNdt[loc][1]=Nu[j]*dNdv[i];
//	//		loc++;
//	//	}
//	//}
//	//for(i=0; i<16; i++)
//	//{
//	//	dNdt[i][0]=(dNdt[i][0]*sum-Nx[i]*sum_u)/(sum*sum);
//	//	dNdt[i][1]=(dNdt[i][1]*sum-Nx[i]*sum_v)/(sum*sum);
//	//	Nx[i]/=sum;
//	//}
//	//double dxdt[2][2]={{0.,0.},{0.,0.}};
//	//loc=0;
//	//for(i=0;i<4;i++)
//	//{
//	//	for(j=0;j<4;j++)
//	//	{
//	//		for (a=0;a<2;a++)
//	//		{
//	//			for(b=0;b<2;b++)
//	//			{
//	//				dxdt[a][b]+=pt[loc][a]*dNdt[loc][b];
//	//			}
//	//		}
//	//		loc++;
//	//	}
//	//}
//	//detJ=dxdt[0][0]*dxdt[1][1]-dxdt[0][1]*dxdt[1][0];
//	//double dtdx[2][2]={{dxdt[1][1]/detJ,-dxdt[0][1]/detJ},{-dxdt[1][0]/detJ,dxdt[0][0]/detJ}};
//	//for(i=0;i<16;i++)
//	//{
//	//	dNdx[i][0]=dNdt[i][0]*dtdx[0][0]+dNdt[i][1]*dtdx[1][0];
//	//	dNdx[i][1]=dNdt[i][0]*dtdx[0][1]+dNdt[i][1]*dtdx[1][1];
//	//}
//	//detJ=0.25*detJ;
//}

void LinearElasticity::ElementMatrix(double dNdx[64][3], double detJ, double EK[192][192])
{
	int i,j;
	for(i=0;i<64;i++)
	{
		for(j=0;j<64;j++)
		{
			EK[dim*i][dim*j] += ((lambda + 2.*mu)*dNdx[i][0] * dNdx[j][0] + mu*(dNdx[i][1] * dNdx[j][1] + dNdx[i][2] * dNdx[j][2]))*detJ;
			EK[dim*i+1][dim*j+1] += ((lambda + 2.*mu)*dNdx[i][1] * dNdx[j][1] + mu*(dNdx[i][0] * dNdx[j][0] + dNdx[i][2] * dNdx[j][2]))*detJ;
			EK[dim*i + 2][dim*j + 2] += ((lambda + 2.*mu)*dNdx[i][2] * dNdx[j][2] + mu*(dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1]))*detJ;
			EK[dim*i][dim*j+1]+=(lambda*dNdx[i][0]*dNdx[j][1]+mu*dNdx[i][1]*dNdx[j][0])*detJ;
			EK[dim*i][dim*j + 2] += (lambda*dNdx[i][0] * dNdx[j][2] + mu*dNdx[i][2] * dNdx[j][0])*detJ;
			EK[dim*i + 1][dim*j] += (lambda*dNdx[i][1] * dNdx[j][0] + mu*dNdx[i][0] * dNdx[j][1])*detJ;
			EK[dim*i + 1][dim*j+2] += (lambda*dNdx[i][1] * dNdx[j][2] + mu*dNdx[i][2] * dNdx[j][1])*detJ;
			EK[dim*i + 2][dim*j] += (lambda*dNdx[i][2] * dNdx[j][0] + mu*dNdx[i][0] * dNdx[j][2])*detJ;
			EK[dim*i + 2][dim*j+1] += (lambda*dNdx[i][2] * dNdx[j][1] + mu*dNdx[i][1] * dNdx[j][2])*detJ;
		}
	}
}

//void LinearElasticity::ReadMeshFromVTK(string fname)
//{
//	pts.clear();
//	eles.clear();
//	string stmp;
//	int npts,neles,itmp,i;
//	vector<int> tmp;
//	ifstream fin;
//	fin.open(fname);
//	if(fin.is_open())
//	{
//		for(i=0;i<4;i++) getline(fin,stmp);//skip lines
//		fin>>stmp>>npts>>stmp;
//		pts.resize(npts);
//		tmp.resize(npts);
//		for(i=0;i<npts;i++)
//		{
//			fin>>pts[i].coor[0]>>pts[i].coor[1]>>pts[i].coor[2];
//		}
//		getline(fin,stmp);
//		fin>>stmp>>neles>>itmp;
//		eles.resize(neles);
//		for(i=0;i<neles;i++)
//		{
//			fin>>itmp>>eles[i].cnct[0]>>eles[i].cnct[1]>>eles[i].cnct[2]>>eles[i].cnct[3];
//		}
//	}
//	else
//	{
//		cerr<<"Cannot open "<<fname<<"!\n";
//	}
//	fin.close();
//}

//void LinearElasticity::SetProblem()
//{
//	//manually assign local knot intervals
//	double kiv[8]={0.,0.,0.,1.,1.,0.,0.,0.};
//	for(int i=0;i<25;i++)
//	{
//		for(int j=0;j<4;j++)
//		{
//			pts[i].kitvU[j]=kiv[i%5+j];
//			pts[i].kitvV[j]=kiv[i/5+j];
//		}
//	}
//
//	//manually construct elements
//	int ec[4][4]={{6,7,12,11},{7,8,13,12},{11,12,17,16},{12,13,18,17}};
//	int env[4][16]={{0,1,2,3,5,6,7,8,10,11,12,13,15,16,17,18},{1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19},{5,6,7,8,10,11,12,13,15,16,17,18,20,21,22,23},{6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24}};
//	eles.resize(4);
//	for(int i=0;i<4;i++)
//	{
//		for(int j=0;j<4;j++)
//		{
//			eles[i].cnct[j]=ec[i][j];
//		}
//		for(int j=0;j<16;j++)
//		{
//			eles[i].vertx_neibors[j]=env[i][j];
//		}
//	}
//}

void LinearElasticity::SetProblem(const vector<int>& IDBC_in, const vector<double>& gh_in)
{
	npt=IDBC_in.size();
	IDBC = IDBC_in;
	gh = gh_in;
	neq = 0;
	int pid(0);
	for (unsigned int i = 0; i < IDBC.size(); i++)
	{
		if (IDBC[i] != -1) neq++;
		//cout << IDBC[i] << " ";
		//if (i % 3==2)
		//{
		//	pid++;
		//	cout << "\n";
		//	cout <<"pid: "<< pid << "\n";
		//}
	}
	//getchar();

	//unsigned int i,j,dof;
	//BCList.clear();
	//for(i=0; i<tts.cp.size(); i++)
	//{
	//	if(tts.cp[i].act==1)
	//	{
	//		for(dof=0; dof<2; dof++)
	//		{
	//			if(dof==0 && tts.cp[i].coor[dof]==0.)//constrain left boundary in x-direction
	//			{
	//				BCList.push_back(2*tts.paid[i]+dof);
	//			}
	//			else if(dof==1 && tts.cp[i].coor[dof]==0.)//constrain bottom boundary in y-direction
	//			{
	//				BCList.push_back(2*tts.paid[i]+dof);
	//			}
	//			else if(dof==0 && tts.cp[i].coor[dof]==1.)//apply 0.1 disp in x-direction at right boundary
	//			{
	//				BCList.push_back(2*tts.paid[i]+dof);
	//			}
	//			npt++;
	//		}
	//	}
	//}
	//gh.clear();
	//gh.resize(npt,0.);
	//for(i=0; i<tts.cp.size(); i++)
	//{
	//	if(tts.cp[i].act==1)
	//	{
	//		for(dof=0; dof<2; dof++)
	//		{
	//			if(dof==0 && tts.cp[i].coor[dof]==1.)//apply 0.1 disp in x-direction at right boundary
	//			{
	//				gh[2*tts.paid[i]+dof]=0.1;
	//			}
	//		}
	//	}
	//}
}

void LinearElasticity::SetBoundary()
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

void LinearElasticity::Assembly(double EK1[192][192], double EF1[192], const vector<vector<double>>& cmat, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF)
{
	unsigned int i,j,a,b,A,B;
	vector<vector<double>> EK(dim*IEN.size(),vector<double>(dim*IEN.size(),0.));
	for(i=0; i<IEN.size(); i++)
	{
		for(j=0; j<IEN.size(); j++)
		{
			for(A=0; A<64; A++)
			{
				for(B=0; B<64; B++)
				{
					for (a = 0; a < dim; a++)
					{
						for (b = 0; b < dim; b++)
						{
							EK[dim*i+a][dim*j+b] += cmat[i][A] * cmat[j][B] * EK1[dim * A+a][dim * B+b];
						}
					}
				}
			}
		}
	}
	for(i=0; i<IEN.size(); i++)
	{
		for(j=0; j<IEN.size(); j++)
		{
			for(a=0; a<dim; a++)
			{
				for(b=0; b<dim; b++)
				{
					A=dim*IEN[i]+a; B=dim*IEN[j]+b;
					if(IDBC[A]!=-1 && IDBC[B]!=-1)
					{
						GK.coeffRef(IDBC[A],IDBC[B]) += EK[dim*i+a][dim*j+b];
					}
					else if(IDBC[A]!=-1 && IDBC[B]==-1)
					{
						GF(IDBC[A]) -= EK[dim*i+a][dim*j+b]*gh[B];
					}
				}
			}
		}
	}
}

void LinearElasticity::BuildLinearSystem(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF)
{
	unsigned int e,i,j,k,a,b;
	double EK[192][192];
	double EF[192];
	double Nx[64];
	double dNdx[64][3];
	double detJ;
	for(e=0; e<bzmesh.size(); e++)
	{
		//if(bzmesh[e].order==3/* || bzmesh[e].order==4*/)
		{
			cout << "element id: " << e << "/" << bzmesh.size()<< "\n";
			for(i=0; i<192; i++)
			{
				EF[i]=0.;
				for(j=0; j<192; j++)
				{
					EK[i][j]=0.;
				}
			}
			for(i=0; i<Gpt.size(); i++)
			{
				for(j=0; j<Gpt.size(); j++)
				{
					for (k = 0; k < Gpt.size(); k++)
					{
						BasisFunction(Gpt[i], Gpt[j], Gpt[k], bzmesh[e].pts, Nx, dNdx, detJ);
						//cout << detJ << "\n";
						//if (detJ < 0.) detJ = -detJ;
						detJ = wght[i] * wght[j] * wght[k] * detJ;
						ElementMatrix(dNdx, detJ, EK);
						//for (a = 0; a < 192; a++)
						//{
						//	if (EK[a][a]<0.)
						//		cout << EK[a][a] << " ";
						//}
						//cout << "\n";
						//getchar();
					}
				}
			}
			Assembly(EK,EF,bzmesh[e].cmat,bzmesh[e].IEN,GK,GF);
			//cout << GK.coeffRef(0, 0) << "\n";
			//getchar();
		}
		//else if(bzmesh[e].order==4)
		//{
		//	double EK[50][50];
		//	double EF[50];
		//	for(i=0; i<50; i++)
		//	{
		//		EF[i]=0.;
		//		for(j=0; j<50; j++)
		//		{
		//			EK[i][j]=0.;
		//		}
		//	}
		//	double Nx[25];
		//	double dNdx[25][2];
		//	double detJ;
		//	for(i=0; i<Gpt.size(); i++)
		//	{
		//		for(j=0; j<Gpt.size(); j++)
		//		{
		//			BasisFunction4(Gpt[i],Gpt[j],bzmesh[e].pts4,Nx,dNdx,detJ);
		//			detJ=wght[i]*wght[j]*detJ;
		//			ElementMatrix4(dNdx,detJ,EK);
		//		}
		//	}
		//	Assembly4(EK,EF,bzmesh[e].cmat4,bzmesh[e].IEN,GK,GF);
		//}
	}
}

void LinearElasticity::Solver(SparseMatrix<double>& GK, VectorXd& GF)
{
	SimplicialLDLT<SparseMatrix<double>> solver;
	VectorXd sol = solver.compute(GK).solve(GF);

	//ConjugateGradient<SparseMatrix<double>, Eigen::Upper> solver;
	//VectorXd sol = solver.compute(GK).solve(GF);

	uh.resize(npt);
	for(int i=0; i<npt; i++)
	{
		if (IDBC[i] != -1)
		{
			uh[i] = sol(IDBC[i]);
			//uh[i] = 1.;
		}
		else
		{
			uh[i] = gh[i];
			//uh[i] = 1.;
		}
		//cout << uh[i] << " ";
		//if (i % 3 == 2) cout << "\n";
	}
}

//void LinearElasticity::Para2Phys(double u, double v, double w, const double cpt[16][3], double pt[3])
//{
//	double Nu[4]={(1.-u)*(1.-u)*(1.-u),3.*(1.-u)*(1.-u)*u,3.*(1.-u)*u*u,u*u*u};
//	double Nv[4]={(1.-v)*(1.-v)*(1.-v),3.*(1.-v)*(1.-v)*v,3.*(1.-v)*v*v,v*v*v};
//	int i,j,loc(0);
//	double tmp;
//	pt[0]=0.; pt[1]=0.; pt[2]=0.;
//	for(i=0;i<4;i++)
//	{
//		for(j=0;j<4;j++)
//		{
//			tmp=Nu[j]*Nv[i];
//			pt[0]+=tmp*cpt[loc][0];
//			pt[1]+=tmp*cpt[loc][1];
//			pt[2]+=tmp*cpt[loc][2];
//			loc++;
//		}
//	}
//}

void LinearElasticity::DispStrainStress(double u, double v, double w, const BezierElement3D& bzel, double disp[3], double se[6], double ss[6])
{
	double Nx[64];
	double dNdx[64][3];
	double detJ;
	BasisFunction(u,v,w,bzel.pts,Nx,dNdx,detJ);
	double uloc[64][3];
	unsigned int i,j;
	for(i=0;i<64;i++)
	{
		uloc[i][0] = 0.; uloc[i][1] = 0.; uloc[i][2] = 0.;
		for(j=0;j<bzel.IEN.size();j++)
		{
			uloc[i][0]+=bzel.cmat[j][i]*uh[dim*bzel.IEN[j]];
			uloc[i][1]+=bzel.cmat[j][i]*uh[dim*bzel.IEN[j]+1];
			uloc[i][2] += bzel.cmat[j][i] * uh[dim*bzel.IEN[j] + 2];
		}
	}
	//displacement
	disp[0] = 0.; disp[1] = 0.; disp[2] = 0.;
	for(i=0;i<64;i++)
	{
		disp[0]+=Nx[i]*uloc[i][0];
		disp[1]+=Nx[i]*uloc[i][1];
		disp[2] += Nx[i] * uloc[i][2];
	}
	//strain
	se[0]=0.; se[1]=0.; se[2]=0.;
	se[3] = 0.; se[4] = 0.; se[5] = 0.;
	for(i=0;i<64;i++)
	{
		se[0]+=dNdx[i][0]*uloc[i][0];
		se[1]+=dNdx[i][1]*uloc[i][1];
		se[2] += dNdx[i][2] * uloc[i][2];
		//se[3]+=(dNdx[i][1]*uloc[i][0]+dNdx[i][0]*uloc[i][1]);
		//se[4] += (dNdx[i][1] * uloc[i][0] + dNdx[i][0] * uloc[i][1]);
	}
	//stress
	ss[0]=(lambda+2.*mu)*se[0]+lambda*(se[1]+se[2]);
	ss[1]=(lambda+2.*mu)*se[1]+lambda*(se[0]+se[2]);
	ss[2] = (lambda + 2.*mu)*se[2] + lambda*(se[1] + se[2]);
	ss[3] = 0.;
	ss[4] = 0.;
	ss[5] = 0.;
}

void LinearElasticity::VisualizeVTK(const vector<BezierElement3D>& bzmesh, string fn)
{
	vector<array<double,3>> spt;//s means sample
	vector<array<double,3>> sdisp;
	vector<array<double,3>> sdisp_err;
	vector<array<double,3>> sse;
	vector<array<double,3>> sss;
	vector<array<int,8>> sele;
	vector<double> errL2;
	vector<array<double,3>> lpt;//visulize parameter lines
	vector<array<int,2>> led;//line connectivity
	int ns(2),ecount(0);
	vector<double> su(ns);
	for(int i=0;i<ns;i++)
	{
		su[i]=double(i)/(double(ns)-1.);
	}

	for(unsigned int e=0;e<bzmesh.size();e++)
	{
		//double errtmp;
		//ElementError(bzmesh[e],errtmp);
		//errL2.push_back(sqrt(errtmp));

		int loc(0);
		for(int a=0;a<ns;a++)
		{
			for(int b=0;b<ns;b++)
			{
				for (int c = 0; c < ns; c++)
				{
					double pt1[3];
					double disp1[3], se1[6], ss1[6];
					//if (bzmesh[e].order == 3)
					{
						bzmesh[e].Para2Phys(su[c], su[b], su[a], pt1);
						DispStrainStress(su[c],su[b], su[a], bzmesh[e], disp1, se1, ss1);
						//DispStrainStress_TSP(su[b],su[a],bzmesh[e],disp1,se1,ss1);
					}
					//else if(bzmesh[e].order==4)
					//{
					//	Para2Phys4(su[b],su[a],bzmesh[e].pts4,pt1);
					//	DispStrainStress4(su[b],su[a],bzmesh[e],disp1,se1,ss1);
					//}
					array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
					array<double, 3> disp = { disp1[0], disp1[1], disp1[2] };
					array<double, 3> se = { se1[0], se1[1], se1[2] };
					array<double, 3> ss = { ss1[0], ss1[1], ss1[2] };
					spt.push_back(pt);
					sdisp.push_back(disp);
					sse.push_back(se);
					sss.push_back(ss);
					//double ue[2], ss_ex[3];
					//ExactDisp_PlateHole(pt1[0],pt1[1],ue);
					//ExactStress_PlateHole(pt1[0],pt1[1],ss_ex);
					////array<double,2> uerr={ue[0]-disp1[0],ue[1]-disp1[1]};
					//array<double,2> uerr={ss_ex[0]-ss1[0],ss_ex[1]-ss1[1]};
					//if(uerr[0]<0.) uerr[0]=-uerr[0];
					//if(uerr[1]<0.) uerr[1]=-uerr[1];
					//sdisp_err.push_back(uerr);
					if (a == 0 || a == ns - 1 || b == 0 || b == ns - 1)
					{
						lpt.push_back(pt);
					}
				}
			}
		}
		int nns[2] = {ns*ns*ns,ns*ns};
		for(int a=0;a<ns-1;a++)
		{
			for(int b=0;b<ns-1;b++)
			{
				for (int c = 0; c < ns - 1; c++)
				{
					array<int, 8> el;
					el[0] = ecount*nns[0] + a*nns[1] + b*ns+c;
					el[1] = ecount*nns[0] + a*nns[1] + b*ns + c+1;
					el[2] = ecount*nns[0] + a*nns[1] + (b+1)*ns + c+1;
					el[3] = ecount*nns[0] + a*nns[1] + (b+1)*ns + c;
					el[4] = ecount*nns[0] + (a+1)*nns[1] + b*ns + c;
					el[5] = ecount*nns[0] + (a + 1)*nns[1] + b*ns + c + 1;
					el[6] = ecount*nns[0] + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
					el[7] = ecount*nns[0] + (a + 1)*nns[1] + (b + 1)*ns + c;
					sele.push_back(el);
				}
			}
		}
		/*for(int a=0;a<ns-1;a++)
		{
			array<int,2> lc;
			lc[0]=ecount*4*(ns-1)+a;
			lc[1]=ecount*4*(ns-1)+a+1;
			led.push_back(lc);
			lc[0]=ecount*4*(ns-1)+3*ns-4+a;
			lc[1]=ecount*4*(ns-1)+3*ns-4+a+1;
			led.push_back(lc);
		}
		for(int a=0;a<ns-2;a++)
		{
			array<int,2> lc;
			lc[0]=ecount*4*(ns-1)+ns+2*a;
			lc[1]=ecount*4*(ns-1)+ns+2*a+2;
			led.push_back(lc);
			lc[0]=ecount*4*(ns-1)+ns+2*a-1;
			lc[1]=ecount*4*(ns-1)+ns+2*a+1;
			led.push_back(lc);
		}
		array<int,2> lc1;
		lc1[0]=ecount*4*(ns-1);
		lc1[1]=ecount*4*(ns-1)+ns;
		led.push_back(lc1);
		lc1[0]=ecount*4*(ns-1)+3*ns-5;
		lc1[1]=ecount*4*(ns-1)+4*ns-5;
		led.push_back(lc1);*/
		ecount++;
	}

	string fname=fn+".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if(fout.is_open())
	{
		fout<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout<<"POINTS "<<spt.size()<<" float\n";
		for(i=0;i<spt.size();i++)
		{
			fout<<spt[i][0]<<" "<<spt[i][1]<<" "<<spt[i][2]<<"\n";
		}
		fout<<"\nCELLS "<<sele.size()<<" "<<9*sele.size()<<'\n';
		for(i=0;i<sele.size();i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
				<<" " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout<<"\nCELL_TYPES "<<sele.size()<<'\n';
		for(i=0;i<sele.size();i++)
		{
			fout<<"12\n";
		}
		//fout<<"\nPOINT_DATA "<<sdisp.size()<<"\nVECTORS disp float\n";
		//for(i=0;i<sdisp.size();i++)
		//{
		//	fout << sdisp[i][0] << " " << sdisp[i][1] << " " << sdisp[i][2] << "\n";
		//}
		fout<<"\nPOINT_DATA "<<sse.size()<<"\nVECTORS strain float\n";
		for(i=0;i<sse.size();i++)
		{
			fout<<sse[i][0]<<" "<<sse[i][1]<<" "<<sse[i][2]<<"\n";
		}
		//fout<<"POINT_DATA "<<sss.size()<<"\nVECTORS stress float\n";
		//for(i=0;i<sss.size();i++)
		//{
		//	fout<<sss[i][0]<<" "<<sss[i][1]<<" "<<sss[i][2]<<"\n";
		//}
		//fout<<"POINT_DATA "<<sdisp_err.size()<<"\nVECTORS disp_err float\n";
		//for(i=0;i<sdisp_err.size();i++)
		//{
		//	fout<<sdisp_err[i][0]<<" "<<sdisp_err[i][1]<<" 0\n";
		//}

		//fout<<"POINT_DATA "<<sval.size()<<"\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		//for(unit i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i]<<"\n";
		//}

		//fout<<"\nCELL_DATA "<<errL2.size()<<"\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		//for(i=0;i<errL2.size();i++)
		//{
		//	fout<<errL2[i]<<"\n";
		//	//fout<<eles[i].type<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout<<"Cannot open "<<fname<<"!\n";
	}

	//string fname1(fn+"-lines.vtk");
	//ofstream fout1;
	//fout1.open(fname1.c_str());
	//if(fout1.is_open())
	//{
	//	fout1<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	//	fout1<<"POINTS "<<lpt.size()<<" float\n";
	//	for(i=0;i<lpt.size();i++)
	//	{
	//		fout1<<lpt[i][0]<<" "<<lpt[i][1]<<" "<<lpt[i][2]<<"\n";
	//	}
	//	fout1<<"\nCELLS "<<led.size()<<" "<<3*led.size()<<'\n';
	//	for(i=0;i<led.size();i++)
	//	{
	//		fout1<<"2 "<<led[i][0]<<" "<<led[i][1]<<'\n';
	//	}
	//	fout1<<"\nCELL_TYPES "<<led.size()<<'\n';
	//	for(i=0;i<led.size();i++)
	//	{
	//		fout1<<"3\n";
	//	}
	//	fout1.close();
	//}
	//else
	//{
	//	cout<<"Cannot open "<<fname1<<"!\n";
	//}
}

void LinearElasticity::Run(const vector<BezierElement3D>& bzmesh, string fn)
{
	GaussInfo(5);
	//ReadMeshFromVTK("Test1.vtk");
	//SetProblem(tts);
	SetMaterialProp(1.e5,.3);
	//SetBoundary();
	SparseMatrix<double> GK(neq,neq);
	//MatrixXd GK(neq,neq);
	VectorXd GF(neq);
	GK.setZero();
	GF.setZero();
	cout<<"Building linear system...\n";
	BuildLinearSystem(bzmesh,GK,GF);
	//BuildLinearSystem_TSP(bzmesh,GK,GF);
	//OutputMatrix(GK,fn+"_mat");

	cout << "Solving...\n";
	Solver(GK,GF);
	cout << "Finish solving...\n";
	VisualizeVTK(bzmesh,fn);

	//double err=OverallError(bzmesh);
	//cout<<"L2 stress err: "<<err<<"\n";
}

void LinearElasticity::OutputMatrix(SparseMatrix<double>& GK, string fn)
{


	string fname = fn + ".txt";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << GK << "\n";
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

//void LinearElasticity::BasisFunction4(double u, double v, const double pt[25][3], double Nx[25], double dNdx[25][2], double& detJ)
//{
//	double Nu[5]={(1.-u)*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-u)*u,6.*(1.-u)*(1.-u)*u*u,4.*(1.-u)*u*u*u,u*u*u*u};
//	double Nv[5]={(1.-v)*(1.-v)*(1.-v)*(1.-v),4.*(1.-v)*(1.-v)*(1.-v)*v,6.*(1.-v)*(1.-v)*v*v,4.*(1.-v)*v*v*v,v*v*v*v};
//	double dNdu[5]={-4.*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-4.*u),12.*u*(1.-3.*u+2.*u*u),4.*(3.-4.*u)*u*u,4.*u*u*u};
//	double dNdv[5]={-4.*(1.-v)*(1.-v)*(1.-v),4.*(1.-v)*(1.-v)*(1.-4.*v),12.*v*(1.-3.*v+2.*v*v),4.*(3.-4.*v)*v*v,4.*v*v*v};
//	double dNdt[25][2];
//	int i,j,a,b,loc(0);
//	for(i=0;i<5;i++)
//	{
//		for(j=0;j<5;j++)
//		{
//			Nx[loc]=Nu[j]*Nv[i];
//			dNdt[loc][0]=dNdu[j]*Nv[i];
//			dNdt[loc][1]=Nu[j]*dNdv[i];
//			loc++;
//		}
//	}
//	double dxdt[2][2]={{0.,0.},{0.,0.}};
//	loc=0;
//	for(i=0;i<5;i++)
//	{
//		for(j=0;j<5;j++)
//		{
//			for (a=0;a<2;a++)
//			{
//				for(b=0;b<2;b++)
//				{
//					dxdt[a][b]+=pt[loc][a]*dNdt[loc][b];
//				}
//			}
//			loc++;
//		}
//	}
//	detJ=dxdt[0][0]*dxdt[1][1]-dxdt[0][1]*dxdt[1][0];
//	double dtdx[2][2]={{dxdt[1][1]/detJ,-dxdt[0][1]/detJ},{-dxdt[1][0]/detJ,dxdt[0][0]/detJ}};
//	for(i=0;i<25;i++)
//	{
//		dNdx[i][0]=dNdt[i][0]*dtdx[0][0]+dNdt[i][1]*dtdx[1][0];
//		dNdx[i][1]=dNdt[i][0]*dtdx[0][1]+dNdt[i][1]*dtdx[1][1];
//	}
//	detJ=0.25*detJ;
//}
//
//void LinearElasticity::BasisFunction4_Rational(double u, double v, const BezierElement3D& bzel, vector<double>& Rx, vector<array<double,2>>& dRdx, double& detJ)
//{
//	double Bt[25], dBdt[25][2], wb(0.), dwbdt[2]={0.,0.};
//	bzel.Basis4(u,v,Bt,dBdt);
//	for(int i=0; i<25; i++)
//	{
//		wb+=Bt[i]*bzel.w4[i];
//		dwbdt[0]+=dBdt[i][0]*bzel.w4[i];
//		dwbdt[1]+=dBdt[i][1]*bzel.w4[i];
//	}
//	Rx.clear();
//	dRdx.clear();
//	Rx.resize(bzel.IEN.size(),0.);
//	dRdx.resize(bzel.IEN.size());
//	vector<array<double,2>> dRdt(bzel.IEN.size());
//	for(unsigned int i=0; i<bzel.IEN.size(); i++)
//	{
//		dRdt[i][0]=0.; dRdt[i][1]=0.; 
//		for(int j=0; j<25; j++)
//		{
//			Rx[i]+=bzel.w_nurbs[i]*bzel.cmat4[i][j]*Bt[j]/wb;
//			dRdt[i][0]+=bzel.w_nurbs[i]*bzel.cmat4[i][j]*(dBdt[j][0]/wb-dwbdt[0]*Bt[j]/(wb*wb));
//			dRdt[i][1]+=bzel.w_nurbs[i]*bzel.cmat4[i][j]*(dBdt[j][1]/wb-dwbdt[1]*Bt[j]/(wb*wb));
//		}
//	}
//	for(int i=0; i<25; i++)
//	{
//		dBdt[i][0]=bzel.w4[i]*(dBdt[i][0]*wb-Bt[i]*dwbdt[0])/(wb*wb);
//		dBdt[i][1]=bzel.w4[i]*(dBdt[i][1]*wb-Bt[i]*dwbdt[1])/(wb*wb);
//	}
//	int loc(0);
//	double dxdt[2][2]={{0.,0.},{0.,0.}};
//	for(int i=0;i<25;i++)
//	{
//		for (int a=0;a<2;a++)
//		{
//			for(int b=0;b<2;b++)
//			{
//				dxdt[a][b]+=bzel.pts4[i][a]*dBdt[i][b];
//			}
//		}
//	}
//
//	detJ=dxdt[0][0]*dxdt[1][1]-dxdt[0][1]*dxdt[1][0];
//	double dtdx[2][2]={{dxdt[1][1]/detJ,-dxdt[0][1]/detJ},{-dxdt[1][0]/detJ,dxdt[0][0]/detJ}};
//	for(unsigned int i=0; i<bzel.IEN.size(); i++)
//	{
//		dRdx[i][0]=dRdt[i][0]*dtdx[0][0]+dRdt[i][1]*dtdx[1][0];
//		dRdx[i][1]=dRdt[i][0]*dtdx[0][1]+dRdt[i][1]*dtdx[1][1];
//	}
//	detJ=0.25*detJ;
//
//	/*double Nu[5]={(1.-u)*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-u)*u,6.*(1.-u)*(1.-u)*u*u,4.*(1.-u)*u*u*u,u*u*u*u};
//	double Nv[5]={(1.-v)*(1.-v)*(1.-v)*(1.-v),4.*(1.-v)*(1.-v)*(1.-v)*v,6.*(1.-v)*(1.-v)*v*v,4.*(1.-v)*v*v*v,v*v*v*v};
//	double dNdu[5]={-4.*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-4.*u),12.*u*(1.-3.*u+2.*u*u),4.*(3.-4.*u)*u*u,4.*u*u*u};
//	double dNdv[5]={-4.*(1.-v)*(1.-v)*(1.-v),4.*(1.-v)*(1.-v)*(1.-4.*v),12.*v*(1.-3.*v+2.*v*v),4.*(3.-4.*v)*v*v,4.*v*v*v};
//	double dNdt[25][2];
//	int i,j,a,b,loc(0);
//	double sum(0.), sum_u(0.), sum_v(0.);
//	for(i=0;i<5;i++)
//	{
//		for(j=0;j<5;j++)
//		{
//			Nx[loc]=Nu[j]*Nv[i]*w[loc];
//			dNdt[loc][0]=dNdu[j]*Nv[i]*w[loc];
//			dNdt[loc][1]=Nu[j]*dNdv[i]*w[loc];
//			sum+=Nx[loc];
//			sum_u+=dNdt[loc][0];
//			sum_v+=dNdt[loc][1];
//			loc++;
//		}
//	}
//	for(i=0; i<25; i++)
//	{
//		dNdt[i][0]=(dNdt[i][0]*sum-Nx[i]*sum_u)/(sum*sum);
//		dNdt[i][1]=(dNdt[i][1]*sum-Nx[i]*sum_v)/(sum*sum);
//		Nx[i]/=sum;
//	}
//	double dxdt[2][2]={{0.,0.},{0.,0.}};
//	loc=0;
//	for(i=0;i<5;i++)
//	{
//		for(j=0;j<5;j++)
//		{
//			for (a=0;a<2;a++)
//			{
//				for(b=0;b<2;b++)
//				{
//					dxdt[a][b]+=pt[loc][a]*dNdt[loc][b];
//				}
//			}
//			loc++;
//		}
//	}
//	detJ=dxdt[0][0]*dxdt[1][1]-dxdt[0][1]*dxdt[1][0];
//	double dtdx[2][2]={{dxdt[1][1]/detJ,-dxdt[0][1]/detJ},{-dxdt[1][0]/detJ,dxdt[0][0]/detJ}};
//	for(i=0;i<25;i++)
//	{
//		dNdx[i][0]=dNdt[i][0]*dtdx[0][0]+dNdt[i][1]*dtdx[1][0];
//		dNdx[i][1]=dNdt[i][0]*dtdx[0][1]+dNdt[i][1]*dtdx[1][1];
//	}
//	detJ=0.25*detJ;*/
//}
//
//void LinearElasticity::ElementMatrix4(double dNdx[25][2], double detJ, double EK[50][50])
//{
//	int i,j;
//	for(i=0;i<25;i++)
//	{
//		for(j=0;j<25;j++)
//		{
//			EK[2*i][2*j]+=((lambda+2.*mu)*dNdx[i][0]*dNdx[j][0]+mu*dNdx[i][1]*dNdx[j][1])*detJ;
//			EK[2*i+1][2*j+1]+=((lambda+2.*mu)*dNdx[i][1]*dNdx[j][1]+mu*dNdx[i][0]*dNdx[j][0])*detJ;
//			//EK[2*i+1][2*j]+=(lambda*dNdx[i][0]*dNdx[j][1]+mu*dNdx[i][1]*dNdx[j][0])*detJ;
//			//EK[2*i][2*j+1]+=(lambda*dNdx[i][1]*dNdx[j][0]+mu*dNdx[i][0]*dNdx[j][1])*detJ;
//			EK[2*i+1][2*j]+=(lambda*dNdx[i][1]*dNdx[j][0]+mu*dNdx[i][0]*dNdx[j][1])*detJ;
//			EK[2*i][2*j+1]+=(lambda*dNdx[i][0]*dNdx[j][1]+mu*dNdx[i][1]*dNdx[j][0])*detJ;
//		}
//	}
//}
//
//void LinearElasticity::Assembly4(double EK1[50][50], double EF1[50], const vector<array<double,25>>& cmat, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF)
//{
//	unsigned int i,j,a,b,A,B;
//	vector<vector<double>> EK(2*IEN.size(),vector<double>(2*IEN.size(),0.));
//	for(i=0; i<IEN.size(); i++)
//	{
//		for(j=0; j<IEN.size(); j++)
//		{
//			for(a=0; a<25; a++)
//			{
//				for(b=0; b<25; b++)
//				{
//					EK[2*i][2*j] += cmat[i][a]*cmat[j][b]*EK1[2*a][2*b];
//					EK[2*i+1][2*j] += cmat[i][a]*cmat[j][b]*EK1[2*a+1][2*b];
//					EK[2*i][2*j+1] += cmat[i][a]*cmat[j][b]*EK1[2*a][2*b+1];
//					EK[2*i+1][2*j+1] += cmat[i][a]*cmat[j][b]*EK1[2*a+1][2*b+1];
//				}
//			}
//		}
//	}
//	for(i=0; i<IEN.size(); i++)
//	{
//		for(j=0; j<IEN.size(); j++)
//		{
//			for(a=0; a<2; a++)
//			{
//				for(b=0; b<2; b++)
//				{
//					A=2*IEN[i]+a; B=2*IEN[j]+b;
//					if(IDBC[A]!=-1 && IDBC[B]!=-1)
//					{
//						GK.coeffRef(IDBC[A],IDBC[B]) += EK[2*i+a][2*j+b];
//					}
//					else if(IDBC[A]!=-1 && IDBC[B]==-1)
//					{
//						GF(IDBC[A]) -= EK[2*i+a][2*j+b]*gh[B];
//					}
//				}
//			}
//		}
//	}
//}
//
//void LinearElasticity::DispStrainStress4(double u,double v,const BezierElement3D& bzel,double disp[2],double se[3],double ss[3])
//{
//	double Nx[25];
//	double dNdx[25][2];
//	double detJ;
//	BasisFunction4(u,v,bzel.pts4,Nx,dNdx,detJ);
//	double uloc[25][2];
//	unsigned int i,j;
//	for(i=0;i<25;i++)
//	{
//		uloc[i][0]=0.; uloc[i][1]=0.;
//		for(j=0;j<bzel.IEN.size();j++)
//		{
//			uloc[i][0]+=bzel.cmat4[j][i]*uh[2*bzel.IEN[j]];
//			uloc[i][1]+=bzel.cmat4[j][i]*uh[2*bzel.IEN[j]+1];
//		}
//	}
//	//displacement
//	disp[0]=0.; disp[1]=0.;
//	for(i=0;i<25;i++)
//	{
//		disp[0]+=Nx[i]*uloc[i][0];
//		disp[1]+=Nx[i]*uloc[i][1];
//	}
//	//strain
//	se[0]=0.; se[1]=0.; se[2]=0.;
//	for(i=0;i<25;i++)
//	{
//		se[0]+=dNdx[i][0]*uloc[i][0];
//		se[1]+=dNdx[i][1]*uloc[i][1];
//		se[2]+=(dNdx[i][1]*uloc[i][0]+dNdx[i][0]*uloc[i][1]);
//	}
//	//cout.precision(16);
//	//cout<<se[1]<<'\n';
//	//getchar();
//	//stress
//	ss[0]=(lambda+2.*mu)*se[0]+lambda*se[1];
//	ss[1]=(lambda+2.*mu)*se[1]+lambda*se[0];
//	ss[2]=mu*se[2];
//}
//
//void LinearElasticity::Para2Phys4(double u, double v, const double cpt[25][3], double pt[3])
//{
//	double Nu[5]={(1.-u)*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-u)*u,6.*(1.-u)*(1.-u)*u*u,4.*(1.-u)*u*u*u,u*u*u*u};
//	double Nv[5]={(1.-v)*(1.-v)*(1.-v)*(1.-v),4.*(1.-v)*(1.-v)*(1.-v)*v,6.*(1.-v)*(1.-v)*v*v,4.*(1.-v)*v*v*v,v*v*v*v};
//	int i,j,loc(0);
//	double tmp;
//	pt[0]=0.; pt[1]=0.; pt[2]=0.;
//	for(i=0;i<5;i++)
//	{
//		for(j=0;j<5;j++)
//		{
//			tmp=Nu[j]*Nv[i];
//			pt[0]+=tmp*cpt[loc][0];
//			pt[1]+=tmp*cpt[loc][1];
//			pt[2]+=tmp*cpt[loc][2];
//			loc++;
//		}
//	}
//}
//
//void LinearElasticity::BasisFunction_TSP(double u, double v, const BezierElement3D& bzel, vector<double>& Nx, vector<array<double,2>>& dNdx, double& detJ)
//{
//	Nx.clear();
//	dNdx.clear();
//	unsigned int i,j;
//	if(bzel.order==3 /*|| bzel.order==4*/)
//	{
//		double Nx0[16];
//		double dNdx0[16][2];
//		BasisFunction(u,v,bzel.pts,Nx0,dNdx0,detJ);
//		Nx.resize(bzel.cmat.size(),0);
//		dNdx.resize(bzel.cmat.size());
//		for(i=0; i<bzel.cmat.size(); i++)
//		{
//			dNdx[i][0]=0.; dNdx[i][1]=0.;
//			for(j=0; j<bzel.cmat[i].size(); j++)
//			{
//				Nx[i]+=bzel.cmat[i][j]*Nx0[j];
//				dNdx[i][0]+=bzel.cmat[i][j]*dNdx0[j][0];
//				dNdx[i][1]+=bzel.cmat[i][j]*dNdx0[j][1];
//			}
//		}
//		//BasisFunction_Rational(u,v,bzel,Nx,dNdx,detJ);
//	}
//	else if(bzel.order==4)
//	{
//		double Nx0[25];
//		double dNdx0[25][2];
//		BasisFunction4(u,v,bzel.pts4,Nx0,dNdx0,detJ);
//		Nx.resize(bzel.cmat4.size(),0);
//		dNdx.resize(bzel.cmat4.size());
//		for(i=0; i<bzel.cmat4.size(); i++)
//		{
//			dNdx[i][0]=0.; dNdx[i][1]=0.;
//			for(j=0; j<bzel.cmat4[i].size(); j++)
//			{
//				Nx[i]+=bzel.cmat4[i][j]*Nx0[j];
//				dNdx[i][0]+=bzel.cmat4[i][j]*dNdx0[j][0];
//				dNdx[i][1]+=bzel.cmat4[i][j]*dNdx0[j][1];
//			}
//		}
//		//BasisFunction4_Rational(u,v,bzel,Nx,dNdx,detJ);
//	}
//}
//
//void LinearElasticity::ElementMatrix_TSP(const vector<array<double,2>>& dNdx, double detJ, vector<vector<double>>& EK)
//{
//	unsigned int i,j;
//	for(i=0; i<dNdx.size(); i++)
//	{
//		for(j=0; j<dNdx.size(); j++)
//		{
//			EK[2*i][2*j]+=((lambda+2.*mu)*dNdx[i][0]*dNdx[j][0]+mu*dNdx[i][1]*dNdx[j][1])*detJ;
//			EK[2*i+1][2*j+1]+=((lambda+2.*mu)*dNdx[i][1]*dNdx[j][1]+mu*dNdx[i][0]*dNdx[j][0])*detJ;
//			EK[2*i+1][2*j]+=(lambda*dNdx[i][1]*dNdx[j][0]+mu*dNdx[i][0]*dNdx[j][1])*detJ;
//			EK[2*i][2*j+1]+=(lambda*dNdx[i][0]*dNdx[j][1]+mu*dNdx[i][1]*dNdx[j][0])*detJ;
//		}
//	}
//}
//
//void LinearElasticity::ElementForce_TSP(const BezierElement3D& bzel, vector<double>& EF)
//{
//	unsigned int i, j, k;
//	EF.clear();
//	EF.resize(2*bzel.IEN.size(),0.);
//	for(i=0; i<bzel.neum_edge.size(); i++)
//	{
//		if(bzel.neum_edge[i]==0)
//		{
//			int ids[4]={0,1,2,3};
//			for(j=0; j<Gpt.size(); j++)
//			{
//				Neumann_Element(Gpt[j],wght[j],ids,i,bzel,EF);
//			}
//		}
//		else if(bzel.neum_edge[i]==1)
//		{
//			int ids[4]={3,7,11,15};
//			for(j=0; j<Gpt.size(); j++)
//			{
//				Neumann_Element(Gpt[j],wght[j],ids,i,bzel,EF);
//			}
//		}
//		else if(bzel.neum_edge[i]==2)
//		{
//			int ids[4]={12,13,14,15};
//			for(j=0; j<Gpt.size(); j++)
//			{
//				Neumann_Element(Gpt[j],wght[j],ids,i,bzel,EF);
//			}
//		}
//		else if(bzel.neum_edge[i]==3)
//		{
//			int ids[4]={0,4,8,12};
//			for(j=0; j<Gpt.size(); j++)
//			{
//				Neumann_Element(Gpt[j],wght[j],ids,i,bzel,EF);
//			}
//		}
//	}
//}
//
//void LinearElasticity::Neumann_Element(double u, double w, int ids[4], int edid, const BezierElement3D& bzel, vector<double>& EF)
//{
//	unsigned int j, k;
//	double Nu[4]={(1.-u)*(1.-u)*(1.-u),3.*(1.-u)*(1.-u)*u,3.*(1.-u)*u*u,u*u*u};
//	double dNdu[4]={-3.*(1.-u)*(1.-u),3.-12.*u+9.*u*u,3.*(2.-3.*u)*u,3.*u*u};
//	double x[2]={0.,0.};
//	double dxdt[2]={0.,0.};
//	for(j=0; j<4; j++)
//	{
//		x[0]+=bzel.pts[ids[j]][0]*Nu[j];
//		x[1]+=bzel.pts[ids[j]][1]*Nu[j];
//		dxdt[0]+=bzel.pts[ids[j]][0]*dNdu[j];
//		dxdt[1]+=bzel.pts[ids[j]][1]*dNdu[j];
//	}
//	double ds=sqrt(dxdt[0]*dxdt[0]+dxdt[1]*dxdt[1]);
//	double nm[2]={-dxdt[1]/ds,dxdt[0]/ds};
//	if(bzel.neum_edge[edid]==0 || bzel.neum_edge[edid]==1)
//	{
//		nm[0]=-nm[0]; nm[1]=-nm[1];
//	}
//	//cout<<nm[0]<<" "<<nm[1]<<"\n";
//	//getchar();
//	ds*=0.5;
//	double ss[3], h[2];
//	ExactStress_PlateHole(x[0],x[1],ss);
//	h[0]=ss[0]*nm[0]+ss[2]*nm[1];
//	h[1]=ss[2]*nm[0]+ss[1]*nm[1];
//	double Ntmp(0.);
//	for(j=0; j<bzel.neum_ID[edid].size(); j++)
//	{
//		Ntmp=0.;
//		for(k=0; k<4; k++)
//		{
//			Ntmp+=bzel.cmat[bzel.neum_ID[edid][j]][ids[k]]*Nu[k];
//		}
//		EF[2*bzel.neum_ID[edid][j]]+=h[0]*Ntmp*ds*w;
//		EF[2*bzel.neum_ID[edid][j]+1]+=h[1]*Ntmp*ds*w;
//	}
//}
//
//void LinearElasticity::Assembly_TSP(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF)
//{
//	unsigned int i,j,a,b,A,B;
//	for(i=0; i<IEN.size(); i++)
//	{
//		for(j=0; j<IEN.size(); j++)
//		{
//			for(a=0; a<2; a++)
//			{
//				for(b=0; b<2; b++)
//				{
//					A=2*IEN[i]+a; B=2*IEN[j]+b;
//					if(IDBC[A]!=-1 && IDBC[B]!=-1)
//					{
//						GK.coeffRef(IDBC[A],IDBC[B]) += EK[2*i+a][2*j+b];
//					}
//					else if(IDBC[A]!=-1 && IDBC[B]==-1)
//					{
//						GF(IDBC[A]) -= EK[2*i+a][2*j+b]*gh[B];
//					}
//				}
//			}
//		}
//	}
//}
//
//void LinearElasticity::AssemblyForce_TSP(const vector<double>& EF, const vector<int>& IEN, VectorXd& GF)
//{
//	unsigned int i,a,A;
//	for(i=0; i<IEN.size(); i++)
//	{
//		for(a=0; a<2; a++)
//		{
//			A=2*IEN[i]+a;
//			if(IDBC[A]!=-1)
//			{
//				GF(IDBC[A])+=EF[2*i+a];
//			}
//		}
//	}
//}
//
//void LinearElasticity::BuildLinearSystem_TSP(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF)
//{
//	unsigned int e,i,j,a,b;
//	for(e=0; e<bzmesh.size(); e++)
//	{
//		//if(bzmesh[e].order==3)
//		//{
//			vector<vector<double>> EK(2*bzmesh[e].IEN.size(),vector<double>(2*bzmesh[e].IEN.size(),0.));
//			vector<double> Nx, EF;
//			vector<array<double,2>> dNdx;
//			double detJ;
//			for(i=0; i<Gpt.size(); i++)
//			{
//				for(j=0; j<Gpt.size(); j++)
//				{
//					BasisFunction_TSP(Gpt[i],Gpt[j],bzmesh[e],Nx,dNdx,detJ);
//					detJ=wght[i]*wght[j]*detJ;
//					ElementMatrix_TSP(dNdx,detJ,EK);
//				}
//			}
//			Assembly_TSP(EK,EF,bzmesh[e].IEN,GK,GF);
//
//			if(bzmesh[e].neum_edge.size()!=0)//have neumann bc edge
//			{
//				ElementForce_TSP(bzmesh[e],EF);
//				AssemblyForce_TSP(EF,bzmesh[e].IEN,GF);
//			}
//		//}
//		//else if(bzmesh[e].order==4)
//		//{
//		//	vector<vector<double>> EK(2*bzmesh[e].IEN.size(),vector<double>(2*bzmesh[e].IEN.size(),0.));
//		//	vector<double> Nx, EF;
//		//	vector<array<double,2>> dNdx;
//		//	double detJ;
//		//	for(i=0; i<Gpt.size(); i++)
//		//	{
//		//		for(j=0; j<Gpt.size(); j++)
//		//		{
//		//			BasisFunction_TSP(Gpt[i],Gpt[j],bzmesh[e],Nx,dNdx,detJ);
//		//			detJ=wght[i]*wght[j]*detJ;
//		//			ElementMatrix_TSP(dNdx,detJ,EK);
//		//		}
//		//	}
//		//	Assembly_TSP(EK,EF,bzmesh[e].IEN,GK,GF);
//		//}
//	}
//}
//
//void LinearElasticity::DispStrainStress_TSP(double u,double v,const BezierElement3D& bzel,double disp[2],double se[3],double ss[3])
//{
//	vector<double> Nx;
//	vector<array<double,2>> dNdx;
//	double detJ;
//	BasisFunction_TSP(u,v,bzel,Nx,dNdx,detJ);
//	unsigned int i,j;
//	//displacement
//	disp[0]=0.; disp[1]=0.;
//	for(i=0; i<Nx.size(); i++)
//	{
//		disp[0]+=Nx[i]*uh[2*bzel.IEN[i]];
//		disp[1]+=Nx[i]*uh[2*bzel.IEN[i]+1];
//	}
//	//strain
//	se[0]=0.; se[1]=0.; se[2]=0.;
//	for(i=0; i<dNdx.size(); i++)
//	{
//		se[0]+=dNdx[i][0]*uh[2*bzel.IEN[i]];
//		se[1]+=dNdx[i][1]*uh[2*bzel.IEN[i]+1];
//		se[2]+=(dNdx[i][1]*uh[2*bzel.IEN[i]]+dNdx[i][0]*uh[2*bzel.IEN[i]+1]);
//	}
//	//stress
//	ss[0]=(lambda+2.*mu)*se[0]+lambda*se[1];
//	ss[1]=(lambda+2.*mu)*se[1]+lambda*se[0];
//	ss[2]=mu*se[2];
//}
//
//void LinearElasticity::Para2Phys_TSP(double u, double v, const BezierElement3D& bzel, array<double,3>& pt)
//{
//	if(bzel.order==3)
//	{
//		double Nu[4]={(1.-u)*(1.-u)*(1.-u),3.*(1.-u)*(1.-u)*u,3.*(1.-u)*u*u,u*u*u};
//		double Nv[4]={(1.-v)*(1.-v)*(1.-v),3.*(1.-v)*(1.-v)*v,3.*(1.-v)*v*v,v*v*v};
//		int i,j,loc(0);
//		double tmp;
//		pt[0]=0.; pt[1]=0.; pt[2]=0.;
//		for(i=0;i<4;i++)
//		{
//			for(j=0;j<4;j++)
//			{
//				tmp=Nu[j]*Nv[i];
//				pt[0]+=tmp*bzel.pts[loc][0];
//				pt[1]+=tmp*bzel.pts[loc][1];
//				pt[2]+=tmp*bzel.pts[loc][2];
//				loc++;
//			}
//		}
//	}
//	else if(bzel.order==4)
//	{
//		double Nu[5]={(1.-u)*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-u)*u,6.*(1.-u)*(1.-u)*u*u,4.*(1.-u)*u*u*u,u*u*u*u};
//		double Nv[5]={(1.-v)*(1.-v)*(1.-v)*(1.-v),4.*(1.-v)*(1.-v)*(1.-v)*v,6.*(1.-v)*(1.-v)*v*v,4.*(1.-v)*v*v*v,v*v*v*v};
//		int i,j,loc(0);
//		double tmp;
//		pt[0]=0.; pt[1]=0.; pt[2]=0.;
//		for(i=0;i<5;i++)
//		{
//			for(j=0;j<5;j++)
//			{
//				tmp=Nu[j]*Nv[i];
//				pt[0]+=tmp*bzel.pts4[loc][0];
//				pt[1]+=tmp*bzel.pts4[loc][1];
//				pt[2]+=tmp*bzel.pts4[loc][2];
//				loc++;
//			}
//		}
//	}
//}
//
//void LinearElasticity::SetProblem_PlateHole(const TruncatedTspline& tts)
//{
//	npt=0;
//	unsigned int i,j,dof;
//	BCList.clear();
//	for(i=0; i<tts.cp.size(); i++)
//	{
//		if(tts.cp[i].act==1)
//		{
//			for(dof=0; dof<2; dof++)
//			{
//				if(dof==0 && tts.cp[i].coor[dof]==0.)//constrain right boundary in x-direction
//				{
//					BCList.push_back(2*tts.paid[i]+dof);
//				}
//				else if(dof==1 && tts.cp[i].coor[dof]==0.)//constrain bottom boundary in y-direction
//				{
//					BCList.push_back(2*tts.paid[i]+dof);
//				}
//				//else if(tts.cp[i].coor[0]==-4. || tts.cp[i].coor[1]==4.)//apply disp on left and top boundaries
//				//{
//				//	BCList.push_back(2*tts.paid[i]+dof);
//				//}
//				npt++;
//			}
//		}
//	}
//	gh.clear();
//	gh.resize(npt,0.);
//	//for(i=0; i<tts.cp.size(); i++)
//	//{
//	//	if(tts.cp[i].act==1)
//	//	{
//	//		for(dof=0; dof<2; dof++)
//	//		{
//	//			if(tts.cp[i].coor[0]==-4. || tts.cp[i].coor[1]==4.)//apply disp on left and top boundaries
//	//			{
//	//				double uxy[2];
//	//				ExactDisp_PlateHole(tts.cp[i].coor[0],tts.cp[i].coor[1],uxy);
//	//				if(dof==0)
//	//					gh[2*tts.paid[i]+dof]=uxy[0];
//	//				else
//	//					gh[2*tts.paid[i]+dof]=uxy[1];
//	//			}
//	//		}
//	//	}
//	//}
//}
//
//void LinearElasticity::Quantityh_TSP(double u,double v,const BezierElement3D& bzel,double disp[2],double ss[3], double& detJ)
//{
//	vector<double> Nx;
//	vector<array<double,2>> dNdx;
//	BasisFunction_TSP(u,v,bzel,Nx,dNdx,detJ);
//	unsigned int i,j;
//	//displacement
//	disp[0]=0.; disp[1]=0.;
//	for(i=0; i<Nx.size(); i++)
//	{
//		disp[0]+=Nx[i]*uh[2*bzel.IEN[i]];
//		disp[1]+=Nx[i]*uh[2*bzel.IEN[i]+1];
//	}
//	//strain
//	double se[3];
//	se[0]=0.; se[1]=0.; se[2]=0.;
//	for(i=0; i<dNdx.size(); i++)
//	{
//		se[0]+=dNdx[i][0]*uh[2*bzel.IEN[i]];
//		se[1]+=dNdx[i][1]*uh[2*bzel.IEN[i]+1];
//		se[2]+=(dNdx[i][1]*uh[2*bzel.IEN[i]]+dNdx[i][0]*uh[2*bzel.IEN[i]+1]);
//	}
//	//stress
//	ss[0]=(lambda+2.*mu)*se[0]+lambda*se[1];
//	ss[1]=(lambda+2.*mu)*se[1]+lambda*se[0];
//	ss[2]=mu*se[2];
//}
//
//void LinearElasticity::ElementError(const BezierElement3D& bzel, double& L2)
//{
//	L2=0.;
//	for(unsigned int gid1=0;gid1<Gpt.size();gid1++)
//	{
//		for(unsigned int gid2=0;gid2<Gpt.size();gid2++)
//		{
//			double detJ, uh[2], ssh[3], u[2], ss[3];
//			array<double,3> x;
//			Quantityh_TSP(Gpt[gid1],Gpt[gid2],bzel,uh,ssh, detJ);
//			Para2Phys_TSP(Gpt[gid1],Gpt[gid2],bzel,x);
//			ExactStress_PlateHole(x[0],x[1],ss);
//			L2+=wght[gid1]*wght[gid2]*detJ*(ssh[0]-ss[0])*(ssh[0]-ss[0]);
//			//H1+=wght[gid1]*wght[gid2]*detJ*((val[1]-ux)*(val[1]-ux)+(val[2]-uy)*(val[2]-uy));
//		}
//	}
//}
//
//double LinearElasticity::OverallError(const vector<BezierElement3D>& bzmesh)
//{
//	double err(0.);
//	for(unsigned int i=0; i<bzmesh.size(); i++)
//	{ 
//		double etmp;
//		ElementError(bzmesh[i],etmp);
//		err+=etmp;
//	}
//	return sqrt(err);
//}
//
//void ExactDisp_PlateHole(double x, double y, double u[2])
//{
//	double Ey(1.e5),v(.3),Tx(10.),a(1.),ur[2];
//	double mu=Ey/(2*(1.+v));
//	double k=(3.-v)/(1.+v);
//	double r=sqrt(x*x+y*y);
//	double theta=acos(x/r);
//	//cout<<theta*180./3.1415926<<"\n";
//	//getchar();
//	ur[0]=Tx/(4.*mu)*(r*((k-1)/2.+cos(2*theta))+a*a/r*(1.+(1.+k)*cos(2.*theta))-a*a*a*a/(r*r*r)*cos(2.*theta));
//	ur[1]=Tx/(4.*mu)*((1.-k)*a*a/r-r-a*a*a*a/(r*r*r))*sin(2.*theta);
//	u[0]=ur[0]*cos(theta)-ur[1]*sin(theta);
//	u[1]=ur[0]*sin(theta)+ur[1]*cos(theta);
//}
//
//void ExactStress_PlateHole(double x, double y, double ss[3])
//{
//	double Ey(1.e5),v(.3),Tx(10.),a(1.),ur[2];
//	//double mu=Ey/(2*(1.+v));
//	//double k=(3.-v)/(1.+v);
//	double r=sqrt(x*x+y*y);
//	double theta=acos(x/r);
//	//ur[0]=Tx/(4.*mu)*(r*((k-1)/2.+cos(2*theta))+a*a/r*(1.+(1.+k)*cos(2*theta))-a*a*a*a/(r*r*r)*cos(2*theta));
//	//ur[1]=Tx/(4.*mu)*((1.-k)*a*a/r-r-a*a*a*a/(r*r*r))*sin(2*theta);
//	//u[0]=ur[0]*cos(theta)-ur[1]*sin(theta);
//	//u[1]=ur[0]*sin(theta)+ur[1]*cos(theta);
//	ss[0]=Tx*(1.-a*a/(r*r)*(1.5*cos(2.*theta)+cos(4.*theta))+1.5*a*a*a*a/(r*r*r*r)*cos(4.*theta));
//	ss[1]=-Tx*(a*a/(r*r)*(.5*cos(2.*theta)-cos(4.*theta))+1.5*a*a*a*a/(r*r*r*r)*cos(4.*theta));
//	ss[2]=-Tx*(a*a/(r*r)*(.5*sin(2.*theta)+sin(4.*theta))-1.5*a*a*a*a/(r*r*r*r)*sin(4.*theta));
//}