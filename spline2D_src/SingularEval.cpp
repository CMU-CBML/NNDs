#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SingularEval.h"
using namespace std;

#define IX(i,j,n) ((i)+(n)*(j))

SingularPatchEval::SingularPatchEval()
{
}

void SingularPatchEval::ReadEigenStruct()
{
	FILE * f;
    int i, N, K;
	if ( !(f = fopen ( "../src/ccdata50NT.dat", "rb" )) ) cerr<<"Cannot open!";
	else
	{
		fread ( &Nmax, sizeof(int), 1, f );
		eigen.resize(Nmax-2);
		for(i=0;i<Nmax-2;i++)
		{
			N = i+3;
			K = 2*N+8;
			eigen[i].L.resize(K);
			eigen[i].iV.resize(K*K);
			eigen[i].x.resize(3);
			eigen[i].x[0].resize(16*K);
			eigen[i].x[1].resize(16*K);
			eigen[i].x[2].resize(16*K);

			double* val=new double[K];
			double* vecI=new double[K*K];
			double* Phi=new double[16*K];
			fread ( val, sizeof(double), K, f );
			eigen[i].L.assign(val,val+K);
			fread ( vecI, sizeof(double), K*K, f );
			eigen[i].iV.assign(vecI,vecI+K*K);
			fread ( Phi, sizeof(double), K*16, f );
			eigen[i].x[0].assign(Phi,Phi+16*K);
			fread ( Phi, sizeof(double), K*16, f );
			eigen[i].x[1].assign(Phi,Phi+16*K);
			fread ( Phi, sizeof(double), K*16, f );
			eigen[i].x[2].assign(Phi,Phi+16*K);
			delete []val;
			delete []vecI;
			delete []Phi;
		}
		fclose(f);
	}
}

void SingularPatchEval::Print()
{
	int i, j, k, l, N, K;

	printf ( "Nmax = %d\n\n", Nmax );

   for ( i=0 ; i<Nmax-2 ; i++ )
   {
      N = i+3;
      K = 2*N+8;

      printf ( "N=%d\n\n", N );

      printf ( "Eigenvalues:\n\n" );
	  for ( j=0 ; j<K ; j++ ) printf ( "%6.3f ", eigen[i].L[j] );
      printf ( "\n\n" );

	  getchar();

      printf ( "Inverse of Eigenvectors:\n\n" );

      for ( j=0 ; j<K ; j++ )
      {
		  for ( k=0 ; k<K ; k++ ) printf ( "%6.3f ", eigen[i].iV[IX(j,k,K)] );
         printf ( "\n" );
      }
      printf ( "\n\n" );

	  getchar();

      printf ( "Coefficients of the Eigenbasis functions:\n" );

      for ( k=0 ; k<3 ; k++ )
      {
         printf ( "k=%d:\n", j );
         for ( j=0 ; j<K ; j++ )
         {
            for ( l=0 ; l<16 ; l++ ) printf("%6.3f ",eigen[i].x[k][IX(j,l,K)]);
            printf ( "\n" );
         }
         printf ( "\n" );

		 getchar();
      }
   }
}

void SingularPatchEval::ProjectPoints(vector<Vertex>& Cp,vector<Vertex>& C,int N)
{
	Cp.resize(2*N+8);
	if(Cp.size()!=C.size()) return;
	int num(0);
	for(int i=0;i<2*N+8;i++)
	{
		Cp[i].coor[0]=0.; Cp[i].coor[1]=0.; Cp[i].coor[2]=0.;
		for(int j=0;j<2*N+8;j++)
		{
			Cp[i]=Cp[i]+C[j]*eigen[N-3].iV[IX(i,j,2*N+8)];
		}
	}
}

void SingularPatchEval::ObtainIVcoef(int N, vector<double>& vec)
{
	vec.clear();
	vec.resize(2*N+1);
	for(unsigned int i=0; i<vec.size(); i++)
	{
		vec[i]=eigen[N-3].iV[IX(0,i,2*N+8)];
	}
}

void SingularPatchEval::ObtainIVMat(int N, vector<vector<double>>& mat)
{
	int dim(2*N+8);
	mat.resize(dim,vector<double>(dim,0.));
	for(int i=0; i<dim; i++)
	{
		for(int j=0; j<dim; j++)
		{
			mat[i][j]=eigen[N-3].iV[IX(i,j,2*N+8)];
		}
	}
}

void SingularPatchEval::EvalSpline(const vector<double>& coef,double u,double v,int ip,int N,double& val,double Dval[2])
{
	//double val(0.);
	val=0.;
	Dval[0]=0.; Dval[1]=0.;
	RegularPatchBasis Nu,Nv;
	Nu.Evaluate(u);
	Nv.Evaluate(v);
	//coef[IX(i,loc,2*N+8)]
	int loc(0);
	for(int j=0;j<4;j++)
	{
		for(int i=0;i<4;i++)
		{
			val+=coef[IX(ip,loc,2*N+8)]*Nu.val[i]*Nv.val[j];
			Dval[0]+=coef[IX(ip,loc,2*N+8)]*Nu.Dval[i]*Nv.val[j];
			Dval[1]+=coef[IX(ip,loc,2*N+8)]*Nu.val[i]*Nv.Dval[j];
			loc++;
		}
	}
	//return val;
}
//
//void SingularPatchEval::Evaluate(Vertex& P,double u,double v,vector<Vertex>& C,int N)
//{
//	vector<Vertex> Cp;
//	ProjectPoints(Cp,C,N);
//
//	double tol(0.00001);
//	if(u<tol && v<tol)
//	{
//		P=Cp[0];
//		//cout<<"C: "<<C[0].coor[0]<<" "<<C[0].coor[1]<<" "<<C[0].coor[2]<<'\n';
//		//cout<<"Cp: "<<Cp[0].coor[0]<<" "<<Cp[0].coor[1]<<" "<<Cp[0].coor[2]<<'\n';
//		return;
//	}
//	u=max(tol,u);
//	v=max(tol,v);
//	double tmp1=log(u)/log(2.),tmp2=log(v)/log(2.);
//	//int n=floor(min(-log2(u),-log2(v)))+1,k;
//	int n=floor(min(-tmp1,-tmp2))+1,k;
//	double pow2=pow(2.,double(n-1));
//	u*=pow2; v*=pow2;
//	if(v<0.5)
//	{
//		k=0; u=2.*u-1.; v=2.*v;
//	}
//	else if(u<0.5)
//	{
//		k=2; u=2.*u; v=2.*v-1.;
//	}
//	else
//	{
//		k=1; u=2.*u-1.; v=2.*v-1.;
//	}
//	P.coor[0]=0.; P.coor[1]=0.; P.coor[2]=0.;
//	for(int i=0;i<2*N+8;i++)
//	{
//		double val,Dval[2];
//		EvalSpline(eigen[N-3].x[k],u,v,i,N,val,Dval);
//		P=P+Cp[i]*(pow(eigen[N-3].L[i],double(n-1))*val);
//	}
//}
//
double SingularPatchEval::ExplicitBasis(int bsid,double u,double v,int N)
{
	double tol(0.00001);
	//if(u<tol && v<tol)
	//{
	//	P=Cp[0];
	//	//cout<<"C: "<<C[0].coor[0]<<" "<<C[0].coor[1]<<" "<<C[0].coor[2]<<'\n';
	//	//cout<<"Cp: "<<Cp[0].coor[0]<<" "<<Cp[0].coor[1]<<" "<<Cp[0].coor[2]<<'\n';
	//	return;
	//}
	u=max(tol,u);
	v=max(tol,v);
	double tmp1=log(u)/log(2.),tmp2=log(v)/log(2.);
	//int n=floor(min(-log2(u),-log2(v)))+1,k;
	int n=floor(min(-tmp1,-tmp2))+1,k;
	double pow2=pow(2.,double(n-1));
	u*=pow2; v*=pow2;
	if(v<0.5)
	{
		k=0; u=2.*u-1.; v=2.*v;
	}
	else if(u<0.5)
	{
		k=2; u=2.*u; v=2.*v-1.;
	}
	else
	{
		k=1; u=2.*u-1.; v=2.*v-1.;
	}

	vector<double> bs_all_1(2*N+8);
	vector<double> bs_all_2(2*N+8,0.);
	for(int i=0;i<2*N+8;i++)
	{
		double val,Dval[2];
		EvalSpline(eigen[N-3].x[k],u,v,i,N,val,Dval);
		bs_all_1[i]=pow(eigen[N-3].L[i],double(n-1))*val;
	}
	for(int i=0;i<2*N+8;i++)
	{
		for(int j=0;j<2*N+8;j++)
		{
			bs_all_2[i]+=eigen[N-3].iV[IX(j,i,2*N+8)]*bs_all_1[j];
		}
	}
	return bs_all_2[bsid];
}
//
//void SingularPatchEval::SingularBasis(double u,double v,int N,VectorXd& Nt,MatrixXd& dNdt)
//{
//	double tol(0.00001);
//	u=max(tol,u);
//	v=max(tol,v);
//	double tmp1=log(u)/log(2.),tmp2=log(v)/log(2.);
//	int n=floor(min(-tmp1,-tmp2))+1,k;
//	double pow2=pow(2.,double(n-1));
//	u*=pow2; v*=pow2;
//	if(v<0.5)
//	{
//		k=0; u=2.*u-1.; v=2.*v;
//	}
//	else if(u<0.5)
//	{
//		k=2; u=2.*u; v=2.*v-1.;
//	}
//	else
//	{
//		k=1; u=2.*u-1.; v=2.*v-1.;
//	}
//
//	vector<double> bs_all_1(2*N+8);
//	vector<double> Dbs_all_u(2*N+8);
//	vector<double> Dbs_all_v(2*N+8);
//	//vector<double> bs_all_2(2*N+8,0.);
//	Nt=VectorXd::Zero(2*N+8);
//	dNdt=MatrixXd::Zero(2*N+8,2);
//	for(int i=0;i<2*N+8;i++)
//	{
//		double val,Dval[2];
//		EvalSpline(eigen[N-3].x[k],u,v,i,N,val,Dval);
//		bs_all_1[i]=pow(eigen[N-3].L[i],double(n-1))*val;
//		Dbs_all_u[i]=pow(eigen[N-3].L[i],double(n-1))*Dval[0];
//		Dbs_all_v[i]=pow(eigen[N-3].L[i],double(n-1))*Dval[1];
//	}
//	for(int i=0;i<2*N+8;i++)
//	{
//		for(int j=0;j<2*N+8;j++)
//		{
//			Nt(i)+=eigen[N-3].iV[IX(j,i,2*N+8)]*bs_all_1[j];
//			dNdt(i,0)+=eigen[N-3].iV[IX(j,i,2*N+8)]*Dbs_all_u[j];
//			dNdt(i,1)+=eigen[N-3].iV[IX(j,i,2*N+8)]*Dbs_all_v[j];
//		}
//	}
//}
//
//void SingularPatchEval::EvalSpline2(const vector<double>& coef,double u,double v,int ip,int N,double& val,double Dval[2],double D2val[3])
//{
//	//double val(0.);
//	val=0.;
//	Dval[0]=0.; Dval[1]=0.;
//	D2val[0]=0.; D2val[1]=0.; D2val[2]=0.;
//	RegularPatchBasis Nu,Nv;
//	Nu.Evaluate(u);
//	Nv.Evaluate(v);
//	//coef[IX(i,loc,2*N+8)]
//	int loc(0);
//	for(int j=0;j<4;j++)
//	{
//		for(int i=0;i<4;i++)
//		{
//			val+=coef[IX(ip,loc,2*N+8)]*Nu.val[i]*Nv.val[j];
//			Dval[0]+=coef[IX(ip,loc,2*N+8)]*Nu.Dval[i]*Nv.val[j];
//			Dval[1]+=coef[IX(ip,loc,2*N+8)]*Nu.val[i]*Nv.Dval[j];
//			D2val[0]+=coef[IX(ip,loc,2*N+8)]*Nu.D2val[i]*Nv.val[j];
//			D2val[1]+=coef[IX(ip,loc,2*N+8)]*Nu.Dval[i]*Nv.Dval[j];
//			D2val[2]+=coef[IX(ip,loc,2*N+8)]*Nu.val[i]*Nv.D2val[j];
//			loc++;
//		}
//	}
//	//return val;
//}
//
//void SingularPatchEval::SingularBasis2(double u,double v,int N,VectorXd& Nt,MatrixXd& dNdt,MatrixXd& d2Ndt)
//{
//	double tol(0.00001);
//	u=max(tol,u);
//	v=max(tol,v);
//	double tmp1=log(u)/log(2.),tmp2=log(v)/log(2.);
//	int n=floor(min(-tmp1,-tmp2))+1,k;
//	double pow2=pow(2.,double(n-1));
//	u*=pow2; v*=pow2;
//	if(v<0.5)
//	{
//		k=0; u=2.*u-1.; v=2.*v;
//	}
//	else if(u<0.5)
//	{
//		k=2; u=2.*u; v=2.*v-1.;
//	}
//	else
//	{
//		k=1; u=2.*u-1.; v=2.*v-1.;
//	}
//
//	vector<double> bs_all_1(2*N+8);
//	vector<double> Dbs_all_u(2*N+8);
//	vector<double> Dbs_all_v(2*N+8);
//	vector<double> Dbs_all_uu(2*N+8);
//	vector<double> Dbs_all_uv(2*N+8);
//	vector<double> Dbs_all_vv(2*N+8);
//	//vector<double> bs_all_2(2*N+8,0.);
//	Nt=VectorXd::Zero(2*N+8);
//	dNdt=MatrixXd::Zero(2*N+8,2);
//	d2Ndt=MatrixXd::Zero(2*N+8,3);
//	for(int i=0;i<2*N+8;i++)
//	{
//		double val,Dval[2],D2val[3];
//		EvalSpline2(eigen[N-3].x[k],u,v,i,N,val,Dval,D2val);
//		bs_all_1[i]=pow(eigen[N-3].L[i],double(n-1))*val;
//		Dbs_all_u[i]=pow(eigen[N-3].L[i],double(n-1))*Dval[0];
//		Dbs_all_v[i]=pow(eigen[N-3].L[i],double(n-1))*Dval[1];
//		Dbs_all_uu[i]=pow(eigen[N-3].L[i],double(n-1))*D2val[0];
//		Dbs_all_uv[i]=pow(eigen[N-3].L[i],double(n-1))*D2val[1];
//		Dbs_all_vv[i]=pow(eigen[N-3].L[i],double(n-1))*D2val[2];
//	}
//	for(int i=0;i<2*N+8;i++)
//	{
//		for(int j=0;j<2*N+8;j++)
//		{
//			Nt(i)+=eigen[N-3].iV[IX(j,i,2*N+8)]*bs_all_1[j];
//			dNdt(i,0)+=eigen[N-3].iV[IX(j,i,2*N+8)]*Dbs_all_u[j];
//			dNdt(i,1)+=eigen[N-3].iV[IX(j,i,2*N+8)]*Dbs_all_v[j];
//			d2Ndt(i,0)+=eigen[N-3].iV[IX(j,i,2*N+8)]*Dbs_all_uu[j];
//			d2Ndt(i,1)+=eigen[N-3].iV[IX(j,i,2*N+8)]*Dbs_all_uv[j];
//			d2Ndt(i,2)+=eigen[N-3].iV[IX(j,i,2*N+8)]*Dbs_all_vv[j];
//		}
//	}
//}