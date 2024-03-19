#include "TTSP_2D.h"
#include "KnotInsertion.h"
#include "BSplineBasis.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>
//#include "Matlab_Solver_wap.h"
#include "SingularEval.h"

typedef unsigned int uint;

TruncatedTspline_2D::TruncatedTspline_2D()
{
	//cp.clear();
	//tmesh.clear();
}

double TruncatedTspline_2D::PartitionOfUnity(int eid,double u,double v)//later
{
	unsigned int i,j;
	int loc,loc1;
	double Ni,sum(0.);

	return sum;
}

void TruncatedTspline_2D::VisualizeControlMesh(string fn)
{
	string fname(fn+".vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if(fout.is_open())
	{
		fout<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout<<"POINTS "<<cp.size()<<" float\n";
		for(uint i=0;i<cp.size();i++)
		{
			fout<<cp[i].coor[0]<<" "<<cp[i].coor[1]<<" "<<cp[i].coor[2]<<"\n";
		}
		int nel_act(0);
		for(uint i=0; i<tmesh.size(); i++)
		{
			if(tmesh[i].act==1 && tmesh[i].type!=2 && tmesh[i].type!=3) nel_act++;
		}
		fout<<"\nCELLS "<<nel_act<<" "<<5*nel_act<<'\n';
		for(uint i=0; i<tmesh.size(); i++)
		{
			if(tmesh[i].act==1 && tmesh[i].type!=2 && tmesh[i].type!=3)
			{
			fout<<"4 "<<tmesh[i].cnct[0]<<' '<<tmesh[i].cnct[1]<<' '<<tmesh[i].cnct[2]<<' '<<tmesh[i].cnct[3]<<'\n';
			}
		}
		fout<<"\nCELL_TYPES "<<nel_act<<'\n';
		for(uint i=0; i<nel_act; i++)
		{
			fout<<"9\n";
		}
		//fout<<"\nCELLS "<<eleH[lev].size()<<" "<<5*eleH[lev].size()<<'\n';
		//for(uint i=0;i<eleH[lev].size();i++)
		//{
		//	fout<<"4 "<<eleH[lev][i].cnct[0]<<" "<<eleH[lev][i].cnct[1]<<" "<<eleH[lev][i].cnct[2]<<" "<<eleH[lev][i].cnct[3]<<'\n';
		//}
		//fout<<"\nCELL_TYPES "<<eleH[lev].size()<<'\n';
		//for(uint i=0;i<eleH[lev].size();i++)
		//{
		//	fout<<"9\n";
		//}
		//fout<<"\nCELL_DATA "<<eleH[lev].size()<<"\nSCALARS eact float 1\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<eleH[lev].size();i++)
		//{
		//	fout<<eleH[lev][i].act<<"\n";
		//	//fout<<eleH[lev][i].type<<"\n";
		//}
		//fout<<"POINT_DATA "<<cp.size()<<"\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<cp.size();i++)
		//{
		//	fout<<cp[i].act<<"\n";
		//}


		//fout<<"\nCELLS "<<cp.size()<<" "<<2*cp.size()<<'\n';
		//for(uint i=0;i<cp.size();i++)
		//{
		//	fout<<"1 "<<i<<'\n';
		//}
		//fout<<"\nCELL_TYPES "<<cp.size()<<'\n';
		//for(uint i=0;i<cp.size();i++)
		//{
		//	fout<<"1\n";
		//}
		//fout<<"POINT_DATA "<<cp.size()<<"\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<cp.size();i++)
		//{
		//	fout<<cp[i].trun<<"\n";
		//}

		fout.close();
	}
	else
	{
		cout<<"Cannot open "<<fname<<"!\n";
	}
}

void TruncatedTspline_2D::CollectActives()
{
	eaid.clear();
	for(uint i=0; i<tmesh.size(); i++)
	{
		if(tmesh[i].act==1 && tmesh[i].type!=2 && tmesh[i].type!=3)
		{
			eaid.push_back(i);
		}
	}
}

void TruncatedTspline_2D::StrongBalanceCheck(const vector<int>& rid, vector<int>& rid2)//later
{
	uint i,j,k;
	rid2.clear();
	for(i=0; i<rid.size(); i++)
	{
		if(tmesh[rid[i]].act==1)
		{
		}
	}
}

void TruncatedTspline_2D::FaceIntersectCheck(vector<int>& rid2)//later
{
	rid2.clear();
	for(uint i=0; i<tmesh.size(); i++)
	{
		//tmesh[i].ref=0;
		if(tmesh[i].act==1)
		{
			int ntjc(0);
			for(int j=0; j<4; j++)
			{
				if(tmedge[tmesh[i].edge[j]].act==0)
				{
					ntjc++;
				}
			}
			if(ntjc==1 && tmesh[i].type==2)//boundary element
			{
				//tmesh[i].ref=10;
				rid2.push_back(i);
			}
			else if(ntjc==2)
			{
				int pos(0);
				for(int j=0; j<4; j++)
				{
					if(tmedge[tmesh[i].edge[j]].act==0)
					{
						pos=j;
						break;
					}
				}
				if(tmedge[tmesh[i].edge[(pos+1)%4]].act==0)
				{
					//tmesh[i].ref=20;
				}
				else if(tmedge[tmesh[i].edge[(pos+2)%4]].act==0)
				{
					//tmesh[i].ref=21;
				}
				rid2.push_back(i);
			}
			else if(ntjc==3)
			{
				//tmesh[i].ref=3;
				rid2.push_back(i);
			}
			else if(ntjc==4)
			{
				//tmesh[i].ref=4;
				rid2.push_back(i);
			}
		}
	}
}

void TruncatedTspline_2D::StrongBalanceRefine(const vector<int>& ridsb)
{
	for(uint i=0; i<ridsb.size(); i++)
	{
		if(tmesh[ridsb[i]].act==1)
		{
			//if(tmesh[ridsb[i]].type==0)
			//{
			//	ElementRefine_Square_4(ridsb[i]);
			//}
			//else if(tmesh[ridsb[i]].type==1)
			//{
			//	ElementRefine_Rectangular(ridsb[i]);
			//}
			//else if(tmesh[ridsb[i]].type==2)
			//{
			//	ElementRefine_Boundary(ridsb[i]);
			//}
		}
	}
}

void TruncatedTspline_2D::TargetRefine(const vector<int>& rid)
{
	for(uint i=0; i<rid.size(); i++)
	{
		if(tmesh[rid[i]].act==1)
		{
			//if(tmesh[rid[i]].type==0)
			//{
			//	ElementRefine_Square_4(rid[i]);
			//}
			//else if(tmesh[rid[i]].type==1)
			//{
			//	ElementRefine_Rectangular(rid[i]);
			//}
		}
	}
}

void TruncatedTspline_2D::OneTjunctionRefine(const vector<int>& ridtp)//later
{
	
}

void TruncatedTspline_2D::TjuncExtentCheck()
{
	for(uint i=0; i<tmesh.size(); i++)
	{
		if(tmesh[i].act==1 && tmesh[i].type==0)
		{
			int pos(-1);
			for(int j=0; j<4; j++)
			{
				if(tmedge[tmesh[i].edge[j]].act==0)
				{
					pos=j;
					break;
				}
			}
			if(pos!=-1)
			{
				int ed[2]={(pos+1)%4,(pos+3)%4},ref(0);
				for(int k=0; k<2; k++)
				{
					int ednb(tmedge[tmesh[i].edge[ed[k]]].face[0]);
					if(ednb==i) ednb=tmedge[tmesh[i].edge[ed[k]]].face[1];
					int* it=find(tmesh[ednb].edge,tmesh[ednb].edge+4,ed[k]);
					int loc=it-tmesh[ednb].edge;
					loc=(loc+2)%4;
					if(tmedge[tmesh[ednb].edge[loc]].act==0)
						ref=1;
				}
				if(ref==1)
				{
					//ElementRefine_Square_2(i);
				}
			}
		}
	}
}

void TruncatedTspline_2D::TjuncExtentCheck_1(vector<int>& ridtjx)
{
	ridtjx.clear();
	for(uint i=0; i<tmesh.size(); i++)
	{
		if(tmesh[i].act==1)
		{
			int pos(-1);
			for(int j=0; j<4; j++)
			{
				if(tmedge[tmesh[i].edge[j]].act==0)
				{
					pos=j;
					break;
				}
			}
			if(pos!=-1)
			{
				int ed[2]={(pos+1)%4,(pos+3)%4},ref(0);
				for(int k=0; k<2; k++)
				{
					int ednb(tmedge[tmesh[i].edge[ed[k]]].face[0]);
					if(ednb==i) ednb=tmedge[tmesh[i].edge[ed[k]]].face[1];
					int* it=find(tmesh[ednb].edge,tmesh[ednb].edge+4,tmesh[i].edge[ed[k]]);
					if(it==tmesh[ednb].edge+4)
					{
						vector<int>::iterator it1=find(ridtjx.begin(),ridtjx.end(),ednb);
						if(it1==ridtjx.end())
							ridtjx.push_back(ednb);
					}
				}
			}
		}
	}
}

void TruncatedTspline_2D::TjuncExtentCheck_2(vector<int>& ridtjx)
{
	ridtjx.clear();
	for(uint i=0; i<tmesh.size(); i++)
	{
		if(tmesh[i].act==1 && tmesh[i].type==0)
		{
			int pos(-1);
			for(int j=0; j<4; j++)
			{
				if(tmedge[tmesh[i].edge[j]].act==0)
				{
					pos=j;
					break;
				}
			}
			if(pos!=-1)
			{
				int ed[2]={(pos+1)%4,(pos+3)%4},ref(0);
				for(int k=0; k<2; k++)
				{
					int ednb(tmedge[tmesh[i].edge[ed[k]]].face[0]);
					if(ednb==i) ednb=tmedge[tmesh[i].edge[ed[k]]].face[1];
					int* it=find(tmesh[ednb].edge,tmesh[ednb].edge+4,tmesh[i].edge[ed[k]]);
					int loc=it-tmesh[ednb].edge;
					loc=(loc+2)%4;
					if(tmedge[tmesh[ednb].edge[loc]].act==0)
					{
						vector<int>::iterator it=find(ridtjx.begin(),ridtjx.end(),i);
						if(it==ridtjx.end())
							ridtjx.push_back(i);
					}
				}
			}
		}
	}
}

void TruncatedTspline_2D::VisualizeTMesh(string fn)
{
	string fname(fn+".vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if(fout.is_open())
	{
		fout<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout<<"POINTS "<<cp.size()<<" float\n";
		for(uint i=0;i<cp.size();i++)
		{
			fout<<cp[i].coor[0]<<" "<<cp[i].coor[1]<<" "<<cp[i].coor[2]<<"\n";
		}
		int ned(0);
		for(uint i=0;i<tmedge.size();i++)
		{
			if(tmedge[i].act==1) ned++;
		}
		fout<<"\nCELLS "<<ned<<" "<<3*ned<<'\n';
		for(uint i=0;i<tmedge.size();i++)
		{
			if(tmedge[i].act==1)
				fout<<"2 "<<tmedge[i].pt[0]<<" "<<tmedge[i].pt[1]<<'\n';
		}
		fout<<"\nCELL_TYPES "<<ned<<'\n';
		for(uint i=0;i<tmedge.size();i++)
		{
			if(tmedge[i].act==1)
				fout<<"3\n";
		}
		fout<<"\nCELL_DATA "<<ned<<"\nSCALARS len float 1\nLOOKUP_TABLE default\n";
		for(uint i=0;i<tmedge.size();i++)
		{
			if(tmedge[i].act==1)
				fout<<tmedge[i].len<<"\n";
		}
		fout<<"POINT_DATA "<<cp.size()<<"\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		for(uint i=0;i<cp.size();i++)
		{
			fout<<cp[i].trun<<"\n";
		}


		//fout<<"\nCELLS "<<cp.size()<<" "<<2*cp.size()<<'\n';
		//for(uint i=0;i<cp.size();i++)
		//{
		//	fout<<"1 "<<i<<'\n';
		//}
		//fout<<"\nCELL_TYPES "<<cp.size()<<'\n';
		//for(uint i=0;i<cp.size();i++)
		//{
		//	fout<<"1\n";
		//}
		//fout<<"POINT_DATA "<<cp.size()<<"\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<cp.size();i++)
		//{
		//	fout<<cp[i].trun<<"\n";
		//}

		fout.close();
	}
	else
	{
		cout<<"Cannot open "<<fname<<"!\n";
	}
}

//void TruncatedTspline_2D::ShootRay_Edge(int edid, int pid, double kv[4])
//{
//	kv[2]=tmedge[edid].len;
//	int start(0),end(1);
//	if(tmedge[edid].pt[1]==pid)
//	{
//		start=1; end=0;
//	}
//	if(tmedge[edid].pn[end][0]==0) kv[3]=tmedge[tmedge[edid].pn[end][1]].len;
//	else if(tmedge[edid].pn[end][0]==1)
//	{
//		int fcid(tmedge[edid].pn[end][1]);
//		kv[3]=tmedge[tmesh[fcid].edge[0]].len;
//	}
//	else if(tmedge[edid].pn[end][0]==2)
//	{
//		int fcid(tmedge[edid].pn[end][1]);
//		kv[3]=tmedge[tmesh[fcid].edge[1]].len;
//	}
//	else
//		kv[3]=0.;
//
//	if(tmedge[edid].pn[start][0]==0)
//	{
//		kv[1]=tmedge[tmedge[edid].pn[start][1]].len;
//		int edpre(tmedge[edid].pn[start][1]);
//		int start1(0),end1(1);
//		if(tmedge[edpre].pt[0]==pid)
//		{
//			start1=1; end1=0;
//		}
//		if(tmedge[edpre].pn[start1][0]==0)
//			kv[0]=tmedge[tmedge[edpre].pn[start1][1]].len;
//		else if(tmedge[edpre].pn[start1][0]==1)
//		{
//			int fcid(tmedge[edpre].pn[start1][1]);
//			kv[0]=tmedge[tmesh[fcid].edge[0]].len;
//		}
//		else if(tmedge[edpre].pn[start1][0]==2)
//		{
//			int fcid(tmedge[edpre].pn[start1][1]);
//			kv[0]=tmedge[tmesh[fcid].edge[1]].len;
//		}
//		else
//			kv[0]=0.;
//	}
//	else if(tmedge[edid].pn[start][0]==1)
//	{
//		int fcid(tmedge[edid].pn[start][1]);
//		kv[1]=tmedge[tmesh[fcid].edge[0]].len;
//		int loc;
//		for(int j=0; j<4; j++)
//		{
//			if(tmesh[fcid].pn[j][0]==0 && tmesh[fcid].pn[j][1]==edid)
//			{
//				loc=j; break;
//			}
//		}
//		loc=(loc+2)%4;
//		if(tmesh[fcid].pn[loc][0]==0)
//		{
//			kv[0]=tmedge[tmesh[fcid].pn[loc][1]].len;
//		}
//		else if(tmesh[fcid].pn[loc][0]==1)
//		{
//			kv[0]=tmedge[tmesh[tmesh[fcid].pn[loc][0]].edge[0]].len;
//		}
//		else if(tmesh[fcid].pn[loc][0]==2)
//		{
//			kv[0]=tmedge[tmesh[tmesh[fcid].pn[loc][0]].edge[1]].len;
//		}
//		else
//		{
//			kv[0]=0.;
//		}
//	}
//	else if(tmedge[edid].pn[start][0]==2)
//	{
//		int fcid(tmedge[edid].pn[start][1]);
//		kv[1]=tmedge[tmesh[fcid].edge[1]].len;
//		int loc;
//		for(int j=0; j<4; j++)
//		{
//			if(tmesh[fcid].pn[j][0]==0 && tmesh[fcid].pn[j][1]==edid)
//			{
//				loc=j; break;
//			}
//		}
//		loc=(loc+2)%4;
//		if(tmesh[fcid].pn[loc][0]==0)
//		{
//			kv[0]=tmedge[tmesh[fcid].pn[loc][1]].len;
//		}
//		else if(tmesh[fcid].pn[loc][0]==1)
//		{
//			kv[0]=tmedge[tmesh[tmesh[fcid].pn[loc][0]].edge[0]].len;
//		}
//		else if(tmesh[fcid].pn[loc][0]==2)
//		{
//			kv[0]=tmedge[tmesh[tmesh[fcid].pn[loc][0]].edge[1]].len;
//		}
//		else
//		{
//			kv[0]=0.;
//		}
//	}
//	else
//	{
//		kv[1]=0.;
//		kv[0]=0.;
//	}
//}
//
//void TruncatedTspline_2D::ShootRay_Face(int fcid, int pid, double kv[4])
//{
//	int dir(0),loc;
//	for(int j=0; j<4; j++)
//	{
//		if(tmedge[tmesh[fcid].edge[j]].act==0 && tmedge[tmesh[fcid].edge[j]].midpt==pid)
//		{
//			loc=j; break;
//		}
//	}
//	if(loc==0 || loc==2) dir=1;
//	kv[2]=tmedge[tmesh[fcid].edge[dir]].len;
//	int loc1=(loc+2)%4;
//	if(tmesh[fcid].pn[loc1][0]==0)
//	{
//		kv[3]=tmedge[tmesh[fcid].pn[loc1][1]].len;
//	}
//	else if(tmesh[fcid].pn[loc1][0]==1)
//	{
//		int fcid1(tmesh[fcid].pn[loc1][1]);
//		kv[3]=tmedge[tmesh[fcid1].edge[0]].len;
//	}
//	else if(tmesh[fcid].pn[loc1][0]==2)
//	{
//		int fcid1(tmesh[fcid].pn[loc1][1]);
//		kv[3]=tmedge[tmesh[fcid1].edge[1]].len;
//	}
//	else
//	{
//		kv[3]=0.;
//	}
//
//	if(tmesh[fcid].pn[loc][0]==0)
//	{
//		int edid(tmesh[fcid].pn[loc][1]);
//		kv[1]=tmedge[edid].len;
//		int start(0),end(1);
//		if(tmedge[edid].pt[1]==pid)
//		{
//			start=1; end=0;
//		}
//		if(tmedge[edid].pn[end][0]==0)
//		{
//			kv[0]=tmedge[tmedge[edid].pn[end][1]].len;
//		}
//		else if(tmedge[edid].pn[end][0]==1)
//		{
//			int fcid1(tmedge[edid].pn[end][1]);
//			kv[0]=tmedge[tmesh[fcid1].edge[0]].len;
//		}
//		else if(tmedge[edid].pn[end][0]==2)
//		{
//			int fcid1(tmedge[edid].pn[end][1]);
//			kv[0]=tmedge[tmesh[fcid1].edge[0]].len;
//		}
//		else
//		{
//			kv[0]=0.;
//		}
//	}
//	else
//	{
//		cerr<<"Wrong edge connectivity!\n";
//	}
//}

bool TruncatedTspline_2D::CheckSubKnotVector(const double ku1[5], const double kv1[5], const double ku2[5], const double kv2[5])
{
	if(ku1[0]>=ku2[0] && ku1[4]<=ku2[4] && kv1[0]>=kv2[0] && kv1[4]<=kv2[4])
	{
		double min_u1(10.),min_v1(10.),min_u2(10.),min_v2(10.);
		for(int i=0; i<4; i++)
		{
			double tmp=ku1[i+1]-ku1[i];
			if(tmp!=0. && tmp<min_u1) min_u1=tmp;
			tmp=kv1[i+1]-kv1[i];
			if(tmp!=0. && tmp<min_v1) min_v1=tmp;
			tmp=ku2[i+1]-ku2[i];
			if(tmp!=0. && tmp<min_u2) min_u2=tmp;
			tmp=kv2[i+1]-kv2[i];
			if(tmp!=0. && tmp<min_v2) min_v2=tmp;
		}
		if(min_u1<=min_u2 && min_v1<=min_v2)
		{
			for(int i=0; i<5; i++)
			{
				const double* it1=find(ku2,ku2+5,ku1[i]);
				if(it1==ku2+5)
				{
					return true;
				}
				const double* it2=find(kv2,kv2+5,kv1[i]);
				if(it2==kv2+5)
				{
					return true;
				}
			}
			return false;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

bool TruncatedTspline_2D::CheckSubKnotVector(const array<double,5>& ku1, const array<double,5>& kv1, const array<double,5>& ku2, const array<double,5>& kv2)
{
	if(ku1[0]>=ku2[0] && ku1[4]<=ku2[4] && kv1[0]>=kv2[0] && kv1[4]<=kv2[4])
	{
		double min_u1(10.),min_v1(10.),min_u2(10.),min_v2(10.);
		for(int i=0; i<4; i++)
		{
			double tmp=ku1[i+1]-ku1[i];
			if(tmp!=0. && tmp<min_u1) min_u1=tmp;
			tmp=kv1[i+1]-kv1[i];
			if(tmp!=0. && tmp<min_v1) min_v1=tmp;
			tmp=ku2[i+1]-ku2[i];
			if(tmp!=0. && tmp<min_u2) min_u2=tmp;
			tmp=kv2[i+1]-kv2[i];
			if(tmp!=0. && tmp<min_v2) min_v2=tmp;
		}
		if(min_u1<=min_u2 && min_v1<=min_v2)
		{
			for(int i=0; i<5; i++)
			{
				array<double,5>::const_iterator it1=find(ku2.begin(),ku2.end(),ku1[i]);
				if(it1==ku2.end())
				{
					return true;
				}
				array<double,5>::const_iterator it2=find(kv2.begin(),kv2.end(),kv1[i]);
				if(it2==kv2.end())
				{
					return true;
				}
			}
			return false;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

void TruncatedTspline_2D::InitialConnect()
{
	uint i,j;
	tmedge.clear();
	for(i=0; i<tmesh.size(); i++)
	{
		tmesh[i].act=1;
		tmesh[i].type=0;
		for(j=0; j<4; j++)
		{
			cp[tmesh[i].cnct[j]].face.push_back(i);
		}
	}
	for(i=0; i<cp.size(); i++)
	{
		if(cp[i].face.size()==3 || cp[i].face.size()>4)//extraordinary point
		{
			cp[i].type=2;
		}
	}
	for(i=0; i<tmesh.size(); i++)
	{
		int flag(0),pos(0);
		for(j=0; j<4; j++)
		{
			if(cp[tmesh[i].cnct[j]].type==2)
			{
				flag=1;
				pos=j;
				break;
			}
		}
		if(flag==1)//reorder local node indices starting from EP
		{
			tmesh[i].type=4;
			int cnctnew[4]={tmesh[i].cnct[pos],tmesh[i].cnct[(pos+1)%4],tmesh[i].cnct[(pos+2)%4],tmesh[i].cnct[(pos+3)%4]};
			for(j=0; j<4; j++)
			{
				tmesh[i].cnct[j]=cnctnew[j];
			}
		}
	}
	for(i=0; i<tmesh.size(); i++)//construct edges
	{
		for(j=0; j<4; j++)
		{
			Edge edtmp;
			edtmp.act=1;
			edtmp.pt[0]=tmesh[i].cnct[j];
			edtmp.pt[1]=tmesh[i].cnct[(j+1)%4];
			edtmp.len=1.;
			vector<Edge>::iterator it=find(tmedge.begin(),tmedge.end(),edtmp);
			int edid(it-tmedge.begin());
			if(it==tmedge.end())
			{
				tmedge.push_back(edtmp);
			}
			tmesh[i].edge[j]=edid;
			tmedge[edid].face.push_back(i);
		}
	}
	for(i=0; i<tmedge.size(); i++)
	{
		cp[tmedge[i].pt[0]].edge.push_back(i);
		cp[tmedge[i].pt[1]].edge.push_back(i);
	}
	for(i=0; i<tmedge.size(); i++)
	{
		if(tmedge[i].face.size()==1)
		{
			int eid(tmedge[i].face[0]);
			if(tmesh[eid].type==0) tmesh[eid].type=2;
			else if(tmesh[eid].type==2) tmesh[eid].type=3;
			int* it=find(tmesh[eid].edge,tmesh[eid].edge+4,i);
			int loc(it-tmesh[eid].edge);
			tmedge[tmesh[eid].edge[(loc+1)%4]].len=0.;
			tmedge[tmesh[eid].edge[(loc+3)%4]].len=0.;
		}
	}

	//FindEdgeTopoDirec_1();
	//FindKnotInterval_1();
	//UpdateKnotInterval_1();
	//SetLocalCoorSystem();
	//FindIEN_1();
}

void TruncatedTspline_2D::UpdateConnect()
{
	uint i,j,k;
	for(i=0; i<cp.size(); i++)
	{
		cp[i].face.clear();
		cp[i].edge.clear();
	}
	for(i=0; i<tmedge.size(); i++)
	{
		tmedge[i].face.clear();
	}
	//loop all faces
	for(i=0; i<tmesh.size(); i++)
	{
		if(tmesh[i].act==1)
		{
			for(j=0; j<4; j++)
			{
				cp[tmesh[i].cnct[j]].face.push_back(i);
				if(tmedge[tmesh[i].edge[j]].act==1)
				{
					tmedge[tmesh[i].edge[j]].face.push_back(i);
				}
				else
				{
					cp[tmedge[tmesh[i].edge[j]].midpt].face.push_back(i);
					int chdid[2]={tmedge[tmesh[i].edge[j]].chd[0],tmedge[tmesh[i].edge[j]].chd[1]};
					if(tmedge[chdid[0]].act==1 && tmedge[chdid[1]].act==1)
					{
						tmedge[chdid[0]].face.push_back(i);
						tmedge[chdid[1]].face.push_back(i);
					}
					else
					{
						cerr<<"Configuration not recognized!\n";
					}
				}
			}
		}
	}
	//loop all edges
	for(i=0; i<tmedge.size(); i++)
	{
		if(tmedge[i].act==1)
		{
			cp[tmedge[i].pt[0]].edge.push_back(i);
			cp[tmedge[i].pt[1]].edge.push_back(i);
		}
	}

	for(i=0; i<cp.size(); i++)
	{
		if(cp[i].face.size()==3 && cp[i].type!=2)
		{
			cp[i].type=1;
		}
	}

	//FindEdgeTopoDirec_1();
	//FindKnotInterval_1();//find kitvtmp
	//UpdateKnotInterval_1();
	//SetLocalCoorSystem();
	//FindIEN_1();
}

void TruncatedTspline_2D::FindEdgeTopoDirec()
{
	for(uint i=0; i<tmedge.size(); i++)
	{
		tmedge[i].pn[0][0]=3; tmedge[i].pn[0][1]=-1;
		tmedge[i].pn[1][0]=3; tmedge[i].pn[1][1]=-1;//initialize as end
		if(tmedge[i].act==1)
		{
			for(int j=0; j<2; j++)
			{
				if(cp[tmedge[i].pt[j]].type==0)//regular
				{
					for(uint k=0; k<cp[tmedge[i].pt[j]].edge.size(); k++)
					{
						if(cp[tmedge[i].pt[j]].edge[k] != i)
						{
							int flag(0);
							for(uint i1=0; i1<tmedge[i].face.size(); i1++)
							{
								for(uint j1=0; j1<4; j1++)
								{
									int ed(tmesh[tmedge[i].face[i1]].edge[j1]);
									if(tmedge[ed].act==1 && ed==cp[tmedge[i].pt[j]].edge[k])
									{
										flag=1; break;
									}
									else if(tmedge[ed].act==0 && (tmedge[ed].chd[0]==cp[tmedge[i].pt[j]].edge[k] || tmedge[ed].chd[1]==cp[tmedge[i].pt[j]].edge[k]))
									{
										flag=1; break;
									}
								}
							}
							if(flag==0)
							{
								tmedge[i].pn[j][0]=0;
								tmedge[i].pn[j][1]=cp[tmedge[i].pt[j]].edge[k];
								break;
							}
						}
					}
				}
				else if(cp[tmedge[i].pt[j]].type==1)//T-junctions
				{
					int fid(-1);
					int loc(0);
					for(uint k=0; k<cp[tmedge[i].pt[j]].face.size(); k++)
					{
						int ftmp(cp[tmedge[i].pt[j]].face[k]);
						for(int k1=0; k1<4; k1++)
						{
							if(tmedge[tmesh[ftmp].edge[k1]].act==0 && tmedge[tmesh[ftmp].edge[k1]].midpt==tmedge[i].pt[j])
							{
								loc=k1; fid=ftmp; break;
							}
						}
						if(fid!=-1)
						{
							break;
						}
					}
					if(fid==-1)
					{
						cout<<"edge id: "<<i<<"\n";
						cout<<tmedge[i].pt[0]<<" "<<tmedge[i].pt[1]<<"\n";
						cout<<tmedge[i].prt<<"\n";
						//cout<<cp[tmedge[i].pt[0]].face.size()<<" "<<cp[tmedge[i].pt[1]].face.size()<<"\n";
						cout<<cp[tmedge[i].pt[0]].face[0]<<" "<<cp[tmedge[i].pt[0]].face[1]<<" "<<cp[tmedge[i].pt[0]].face[2]<<"\n";
						cerr<<"T-junction cannot be found in any neighboring elements!\n";
						getchar();
					}
					if(tmedge[tmesh[fid].edge[loc]].chd[0]==i)
					{
						tmedge[i].pn[j][0]=0;
						tmedge[i].pn[j][1]=tmedge[tmesh[fid].edge[loc]].chd[1];
					}
					else if(tmedge[tmesh[fid].edge[loc]].chd[1]==i)
					{
						tmedge[i].pn[j][0]=0;
						tmedge[i].pn[j][1]=tmedge[tmesh[fid].edge[loc]].chd[0];
					}
					else
					{
						tmedge[i].pn[j][0]=1;
						tmedge[i].pn[j][1]=fid;
					}
				}
				else if(cp[tmedge[i].pt[j]].type==2)//extraordinary
				{
					tmedge[i].pn[j][0]=2;
				}
			}
		}
	}
}

void TruncatedTspline_2D::FindKnotInterval()
{
	for(uint i=0; i<cp.size(); i++)
	{
		for(int j=0; j<4; j++)
		{
			cp[i].kitvUtmp[j]=1.; cp[i].kitvVtmp[j]=1.;
		}
		if(cp[i].type!=2)
		{
			int pos(0);
			cp[i].rfc=-1;
			for(uint j=0; j<cp[i].face.size(); j++)
			{
				int* it=find(tmesh[cp[i].face[j]].cnct,tmesh[cp[i].face[j]].cnct+4,i);
				if(it!=tmesh[cp[i].face[j]].cnct+4)
				{
					cp[i].rfc=cp[i].face[j];
					pos=it-tmesh[cp[i].face[j]].cnct;
					break;
				}
			}
			if(cp[i].rfc==-1)
			{
				cerr<<"Cannot find correct reference face!\n";
				getchar();
			}
			cp[i].uved[0]=tmesh[cp[i].rfc].edge[pos];
			cp[i].uved[1]=tmesh[cp[i].rfc].edge[(pos+3)%4];
			if(tmedge[cp[i].uved[0]].act==0)
			{
				int edtmp(tmedge[cp[i].uved[0]].chd[0]);
				if(tmedge[edtmp].pt[0]==i || tmedge[edtmp].pt[1]==i)
				{
					cp[i].uved[0]=tmedge[cp[i].uved[0]].chd[0];
				}
				else
				{
					cp[i].uved[0]=tmedge[cp[i].uved[0]].chd[1];
				}
			}
			if(tmedge[cp[i].uved[1]].act==0)
			{
				int edtmp(tmedge[cp[i].uved[1]].chd[0]);
				if(tmedge[edtmp].pt[0]==i || tmedge[edtmp].pt[1]==i)
				{
					cp[i].uved[1]=tmedge[cp[i].uved[1]].chd[0];
				}
				else
				{
					cp[i].uved[1]=tmedge[cp[i].uved[1]].chd[1];
				}
			}
			vector<int>::iterator it1=find(cp[i].edge.begin(),cp[i].edge.end(),cp[i].uved[0]);
			vector<int>::iterator it2=find(cp[i].edge.begin(),cp[i].edge.end(),cp[i].uved[1]);
			if(it1==cp[i].edge.end() || it2==cp[i].edge.end())
			{
				cerr<<"Cannot find correct uv edges!\n";
				getchar();
			}
			int trun_flag(0);
			ShootRay(i,cp[i].uved[0],cp[i].kitvUtmp,trun_flag);
			ShootRay(i,cp[i].uved[1],cp[i].kitvVtmp,trun_flag);
			//cp[i].trun=trun_flag;
		}
		else
		{
			for(int j=0; j<4; j++)
			{
				cp[i].kitvUtmp[j]=tmedge[cp[i].edge[0]].len; cp[i].kitvVtmp[j]=tmedge[cp[i].edge[0]].len;
			}
		}
	}
}

void TruncatedTspline_2D::UpdateKnotInterval()
{
	for(uint i=0; i<cp.size(); i++)
	{
		//if(cp[i].trun==0)
		//{
			for(int j=0; j<4; j++)
			{
				cp[i].kitvU[j]=cp[i].kitvUtmp[j];
				cp[i].kitvV[j]=cp[i].kitvVtmp[j];
			}
		//}
	}
}

void TruncatedTspline_2D::ShootRay(int pid, int edid, double kv[4], int& trun_flag)
{
	int loc0(0),loc1(1);
	int edge_lev[2]={-1,-1};
	trun_flag=0;
	int flag[2]={-1,-1};
	//positive direction
	if(pid!=tmedge[edid].pt[0])
	{
		loc0=1; loc1=0;
	}
	kv[2]=tmedge[edid].len;
	edge_lev[0]=tmedge[edid].lev;//first edge level
	if(tmedge[edid].pn[loc1][0]==0)//next is edge
	{
		kv[3]=tmedge[tmedge[edid].pn[loc1][1]].len;
		//find skip
		int edtmp(tmedge[edid].pn[loc1][1]);
		edge_lev[1]=tmedge[edtmp].lev;//second edge level
		if(edge_lev[1]>edge_lev[0] && tmedge[edtmp].face.size()==2)
		{
			int eid[2]={tmedge[edtmp].face[0],tmedge[edtmp].face[1]};
			if(tmesh[eid[0]].lev!=tmesh[eid[1]].lev)
			{
				int e_fe(eid[0]);
				if(tmesh[eid[0]].lev > tmesh[eid[1]].lev) e_fe=eid[1];
				int pos(0);
				for(int k=0; k<4; k++)
				{
					if(tmedge[tmesh[e_fe].edge[k]].act==0)
					{
						pos=k;
						break;
					}
				}
				if(cp[tmesh[e_fe].cnct[(pos+3)%4]].type==1)
				{
					flag[0]=1;
				}
			}
		}
	}
	else if(tmedge[edid].pn[loc1][0]==1)//next is face
	{
		int fid(tmedge[edid].pn[loc1][1]), pos(0);
		for(int i=0; i<4; i++)
		{
			if(tmedge[tmesh[fid].edge[i]].act==0 && tmedge[tmesh[fid].edge[i]].midpt==tmedge[edid].pt[loc1])
			{
				pos=i; break;
			}
		}
		kv[3]=tmedge[tmesh[fid].edge[(pos+1)%4]].len;
	}
	else if(tmedge[edid].pn[loc1][0]==2)//next is XP
	{
		kv[3]=kv[2];
	}
	else if(tmedge[edid].pn[loc1][0]==3)//end
	{
		kv[3]=0.;
	}
	if(flag[0]==1)//skip
	{
		kv[3]=2.*kv[3];
		trun_flag=1;
		//cout<<"trun!\n";
		//getchar();
	}
	//negative direction
	edge_lev[0]=-1; edge_lev[1]=-1; 
	if(tmedge[edid].pn[loc0][0]==0)//previous is edge
	{
		int ed0=tmedge[edid].pn[loc0][1];
		kv[1]=tmedge[ed0].len;
		edge_lev[0]=tmedge[ed0].lev;//first edge level
		int a0=0, a1=1;
		if(tmedge[ed0].pt[0]!=pid)
		{
			a0=1; a1=0;
		}
		if(tmedge[ed0].pn[a1][0]==0)
		{
			kv[0]=tmedge[tmedge[ed0].pn[a1][1]].len;
			//find skip
			int edtmp(tmedge[ed0].pn[a1][1]);
			edge_lev[1]=tmedge[edtmp].lev;//second edge level
			if(edge_lev[1]>edge_lev[0] && tmedge[edtmp].face.size()==2)
			{
				int eid[2]={tmedge[edtmp].face[0],tmedge[edtmp].face[1]};
				if(tmesh[eid[0]].lev!=tmesh[eid[1]].lev)
				{
					int e_fe(eid[0]);
					if(tmesh[eid[0]].lev > tmesh[eid[1]].lev) e_fe=eid[1];
					int pos(0);
					for(int k=0; k<4; k++)
					{
						if(tmedge[tmesh[e_fe].edge[k]].act==0)
						{
							pos=k;
							break;
						}
					}
					if(cp[tmesh[e_fe].cnct[(pos+3)%4]].type==1)
					{
						flag[1]=1;
					}
				}
			}
		}
		else if(tmedge[ed0].pn[a1][0]==1)
		{
			int pt0(tmedge[ed0].pt[a1]), fid(tmedge[ed0].pn[a1][1]), pos(0);
			for(int i=0; i<4; i++)
			{
				if(tmedge[tmesh[fid].edge[i]].act==0 && tmedge[tmesh[fid].edge[i]].midpt==pt0)
				{
					pos=i; break;
				}
			}
			kv[0]=tmedge[tmesh[fid].edge[(pos+1)%4]].len;
		}
		else if(tmedge[ed0].pn[a1][0]==2)
		{
			kv[0]=kv[1];
		}
		else if(tmedge[ed0].pn[a1][0]==3)
		{
			kv[0]=0.;
		}
		
	}
	else if(tmedge[edid].pn[loc0][0]==1)
	{
		int fid0(tmedge[edid].pn[loc0][1]), pos(0);
		for(int i=0; i<4; i++)
		{
			if(tmedge[tmesh[fid0].edge[i]].act==0 && tmedge[tmesh[fid0].edge[i]].midpt==pid)
			{
				pos=i; break;
			}
		}
		kv[1]=tmedge[tmesh[fid0].edge[(pos+1)%4]].len;
		int ed0(tmesh[fid0].edge[(pos+2)%4]);
		if(tmedge[ed0].act==1)
		{
			if(tmedge[ed0].face.size()==2)
			{
				int fid1(tmedge[ed0].face[0]);
				if(fid1==fid0) fid1=tmedge[ed0].face[1];
				int* it=find(tmesh[fid1].edge,tmesh[fid1].edge+4,ed0);
				int pos1(it-tmesh[fid1].edge);
				kv[0]=tmedge[tmesh[fid1].edge[(pos1+1)%4]].len;
			}
			else
			{
				kv[0]=0.;
			}
		}
		else
		{
			int pt0(tmedge[ed0].midpt), ed1;
			for(uint i=0; i<cp[pt0].edge.size(); i++)
			{
				if(cp[pt0].edge[i]!=tmedge[ed0].chd[0] && cp[pt0].edge[i]!=tmedge[ed0].chd[1])
				{
					ed1=cp[pt0].edge[i]; break;
				}
			}
			kv[0]=tmedge[ed1].len;
		}
	}
	else if(tmedge[edid].pn[loc0][0]==2)
	{
		kv[1]=kv[2]; kv[0]=kv[2];
	}
	else
	{
		kv[1]=0.; kv[0]=0.;
	}
	if(flag[1]==1)//skip
	{
		kv[0]=2.*kv[0];
		trun_flag=1;
		//cout<<"trun!\n";
		//getchar();
	}
}

void TruncatedTspline_2D::FindIEN_Unstruct()//IEN, not IENtmp
{
	for(uint eid=0; eid<tmesh.size(); eid++)
	{
		tmesh[eid].IEN.clear();
		tmesh[eid].patch_ku.clear();
		tmesh[eid].patch_kv.clear();
		if(tmesh[eid].act==1 && (tmesh[eid].type==0 || tmesh[eid].type==1))//find two ring neighorhood
		{
			array<double,2> urang={0.,tmedge[tmesh[eid].edge[0]].len};
			array<double,2> vrang={0.,tmedge[tmesh[eid].edge[3]].len};
			int count(0);
			//initial
			vector<int> pr0(tmesh[eid].node),er0(1,eid),pr1,er1,pr1_pref,pr1_eref;
			vector<int> rot_ref(tmesh[eid].node.size());
			vector<array<double,2>> uv_ref(tmesh[eid].node.size());
			for(uint i=0; i<tmesh[eid].node.size(); i++)
			{
				rot_ref[i]=tmesh[eid].lcs[i].rot;
				uv_ref[i][0]=tmesh[eid].lcs[i].u[0];
				uv_ref[i][1]=tmesh[eid].lcs[i].u[1];
			}
			for(uint i=0; i<tmesh[eid].node.size(); i++)
			{
				array<double,5> kui,kvi;
				FindLocalKnotVector(tmesh[eid].node[i],rot_ref[i],uv_ref[i],kui,kvi);
				if(CheckSupport(urang,vrang,kui,kvi))
				{
					tmesh[eid].IEN.push_back(tmesh[eid].node[i]);
					tmesh[eid].patch_ku.push_back(kui);
					tmesh[eid].patch_kv.push_back(kvi);
				}
			}
			while(count<2)
			{
				FindNextRing(pr0,er0,pr1,er1,pr1_pref,pr1_eref);
				vector<int> rot_tmp(pr1.size());
				vector<array<double,2>> uv_tmp(pr1.size());
				for(uint i=0; i<pr1.size(); i++)
				{
					array<double,5> kui,kvi;
					FindRotateAndUVCoor(pr0[pr1_pref[i]],rot_ref[pr1_pref[i]],uv_ref[pr1_pref[i]],pr1_eref[i],pr1[i],rot_tmp[i],uv_tmp[i]);
					FindLocalKnotVector(pr1[i],rot_tmp[i],uv_tmp[i],kui,kvi);
					if(CheckSupport(urang,vrang,kui,kvi))
					{
						tmesh[eid].IEN.push_back(pr1[i]);
						tmesh[eid].patch_ku.push_back(kui);
						tmesh[eid].patch_kv.push_back(kvi);
					}
				}
				pr0.clear();
				er0.clear();
				rot_ref.clear();
				uv_ref.clear();
				pr0=pr1;
				er0=er1;
				rot_ref=rot_tmp;
				uv_ref=uv_tmp;
				count++;
			}
		}
		else if(tmesh[eid].act==1 && tmesh[eid].type==4)
		{
			tmesh[eid].IEN.push_back(tmesh[eid].cnct[0]);
			int fc_pre(tmedge[tmesh[eid].edge[3]].face[0]);
			if(fc_pre==eid) fc_pre=tmedge[tmesh[eid].edge[3]].face[1];
			tmesh[eid].IEN.push_back(tmesh[fc_pre].cnct[3]);
			tmesh[eid].IEN.push_back(tmesh[fc_pre].cnct[2]);
			tmesh[eid].IEN.push_back(tmesh[eid].cnct[3]);
			tmesh[eid].IEN.push_back(tmesh[eid].cnct[2]);
			int fc_next(tmedge[tmesh[eid].edge[0]].face[0]);
			if(fc_next==eid) fc_next=tmedge[tmesh[eid].edge[0]].face[1];
			int fc_next0=fc_next;
			while(fc_next!=fc_pre)
			{
				tmesh[eid].IEN.push_back(tmesh[fc_next].cnct[3]);
				tmesh[eid].IEN.push_back(tmesh[fc_next].cnct[2]);
				int fc_nn(tmedge[tmesh[fc_next].edge[0]].face[0]);
				if(fc_nn==fc_next) fc_nn=tmedge[tmesh[fc_next].edge[0]].face[1];
				fc_next=fc_nn;
			}
			for(int j=1; j<4; j++)
			{
				for(uint k=0; k<cp[tmesh[eid].cnct[j]].face.size(); k++)
				{
					int fcid(cp[tmesh[eid].cnct[j]].face[k]);
					if(fcid!=eid && fcid!=fc_pre && fcid!=fc_next0)
					{
						for(uint k1=0; k1<tmesh[fcid].node.size(); k1++)
						{
							vector<int>::iterator it=find(tmesh[eid].IEN.begin(),tmesh[eid].IEN.end(),tmesh[fcid].node[k1]);
							if(it==tmesh[eid].IEN.end())
							{
								array<double,2> uv_ref={tmesh[eid].lcs[j].u[0],tmesh[eid].lcs[j].u[1]},uv;
								int rot;
								array<double,5> kui, kvi;
								FindRotateAndUVCoor(tmesh[eid].cnct[j],tmesh[eid].lcs[j].rot,uv_ref,fcid,tmesh[fcid].node[k1],rot,uv);
								FindLocalKnotVector(tmesh[fcid].node[k1],rot,uv,kui,kvi);
								tmesh[eid].IEN.push_back(tmesh[fcid].node[k1]);
								tmesh[eid].patch_ku.push_back(kui);
								tmesh[eid].patch_kv.push_back(kvi);
							}
						}
					}
				}
			}
		}
	}
}

void TruncatedTspline_2D::FindIEN_Invalid()
{
	for(uint eid=0; eid<tmesh.size(); eid++)
	{
		tmesh[eid].IENtmp.clear();
		tmesh[eid].patch_kutmp.clear();
		tmesh[eid].patch_kvtmp.clear();
		if(tmesh[eid].act==1 && (tmesh[eid].type==0 || tmesh[eid].type==1))//find two ring neighorhood
		{
			array<double,2> urang={0.,tmedge[tmesh[eid].edge[0]].len};
			array<double,2> vrang={0.,tmedge[tmesh[eid].edge[3]].len};
			int count(0);
			//initial
			vector<int> pr0(tmesh[eid].node),er0(1,eid),pr1,er1,pr1_pref,pr1_eref;
			vector<int> rot_ref(tmesh[eid].node.size());
			vector<array<double,2>> uv_ref(tmesh[eid].node.size());
			for(uint i=0; i<tmesh[eid].node.size(); i++)
			{
				rot_ref[i]=tmesh[eid].lcs[i].rot;
				uv_ref[i][0]=tmesh[eid].lcs[i].u[0];
				uv_ref[i][1]=tmesh[eid].lcs[i].u[1];
			}
			for(uint i=0; i<tmesh[eid].node.size(); i++)
			{
				array<double,5> kui,kvi;
				FindLocalKnotVector(tmesh[eid].node[i],rot_ref[i],uv_ref[i],kui,kvi);
				if(CheckSupport(urang,vrang,kui,kvi))
				{
					//if(tmesh[eid].node[i]>=npt_old) tmesh[eid].aff=1;
					tmesh[eid].IENtmp.push_back(tmesh[eid].node[i]);
					tmesh[eid].patch_kutmp.push_back(kui);
					tmesh[eid].patch_kvtmp.push_back(kvi);
				}
			}
			while(count<2)
			{
				FindNextRing(pr0,er0,pr1,er1,pr1_pref,pr1_eref);
				vector<int> rot_tmp(pr1.size());
				vector<array<double,2>> uv_tmp(pr1.size());
				for(uint i=0; i<pr1.size(); i++)
				{
					array<double,5> kui,kvi;
					FindRotateAndUVCoor(pr0[pr1_pref[i]],rot_ref[pr1_pref[i]],uv_ref[pr1_pref[i]],pr1_eref[i],pr1[i],rot_tmp[i],uv_tmp[i]);
					FindLocalKnotVector(pr1[i],rot_tmp[i],uv_tmp[i],kui,kvi);
					if(CheckSupport(urang,vrang,kui,kvi))
					{
						//if(pr1[i]>=npt_old) tmesh[eid].aff=1;
						tmesh[eid].IENtmp.push_back(pr1[i]);
						tmesh[eid].patch_kutmp.push_back(kui);
						tmesh[eid].patch_kvtmp.push_back(kvi);
					}
				}
				pr0.clear();
				er0.clear();
				rot_ref.clear();
				uv_ref.clear();
				pr0=pr1;
				er0=er1;
				rot_ref=rot_tmp;
				uv_ref=uv_tmp;
				count++;
			}
		}
		else if(tmesh[eid].act==1 && tmesh[eid].type==4)
		{
			tmesh[eid].IENtmp.push_back(tmesh[eid].cnct[0]);
			int fc_pre(tmedge[tmesh[eid].edge[3]].face[0]);
			if(fc_pre==eid) fc_pre=tmedge[tmesh[eid].edge[3]].face[1];
			int pt_pre[2]={tmesh[fc_pre].cnct[3],tmesh[fc_pre].cnct[2]};
			if(tmesh[fc_pre].type==5)
			{
				int* it1=find(tmesh[fc_pre].cnct,tmesh[fc_pre].cnct+4,tmesh[eid].cnct[0]);
				int loc1(it1-tmesh[fc_pre].cnct);
				pt_pre[0]=tmesh[fc_pre].cnct[(loc1+3)%4];
				pt_pre[1]=tmesh[fc_pre].cnct[(loc1+2)%4];
			}
			tmesh[eid].IENtmp.push_back(pt_pre[0]);
			tmesh[eid].IENtmp.push_back(pt_pre[1]);
			tmesh[eid].IENtmp.push_back(tmesh[eid].cnct[3]);
			tmesh[eid].IENtmp.push_back(tmesh[eid].cnct[2]);
			int fc_next(tmedge[tmesh[eid].edge[0]].face[0]);
			if(fc_next==eid) fc_next=tmedge[tmesh[eid].edge[0]].face[1];
			int fc_next0=fc_next;
			int count(0);
			while(fc_next!=fc_pre)
			{
				if(fc_next>=tmesh.size())
				{
					//cout<<eid<<" "<<fc_next<<" "<<tmesh.size()<<"\n";
					//cout<<tmesh[eid].cnct[0]<<" "<<tmesh[eid].cnct[1]<<" "<<tmesh[eid].cnct[2]<<" "<<tmesh[eid].cnct[3]<<"\n";
					cout<<tmesh[eid].type<<"\n";
					cout<<cp[tmesh[eid].cnct[0]].face.size()<<"\n";
					cout<<cp[tmesh[eid].cnct[1]].face.size()<<"\n";
					cout<<cp[tmesh[eid].cnct[2]].face.size()<<"\n";
					cout<<cp[tmesh[eid].cnct[3]].face.size()<<"\n";
					getchar();
				}
				int pt_next[2]={tmesh[fc_next].cnct[3],tmesh[fc_next].cnct[2]};
				int loc2(0);
				if(tmesh[fc_next].type==5)
				{
					int* it2=find(tmesh[fc_next].cnct,tmesh[fc_next].cnct+4,tmesh[eid].cnct[0]);
					loc2=it2-tmesh[fc_next].cnct;
					pt_next[0]=tmesh[fc_next].cnct[(loc2+3)%4];
					pt_next[1]=tmesh[fc_next].cnct[(loc2+2)%4];
				}
				tmesh[eid].IENtmp.push_back(pt_next[0]);
				tmesh[eid].IENtmp.push_back(pt_next[1]);
				int fc_nn(tmedge[tmesh[fc_next].edge[loc2]].face[0]);
				if(fc_nn==fc_next) fc_nn=tmedge[tmesh[fc_next].edge[loc2]].face[1];
				fc_next=fc_nn;
				count++;
				if(count>20)
				{
					cerr<<"Loop more than 20 times!\n";
					getchar();
					break;
				}
			}
			for(int j=1; j<4; j++)
			{
				for(uint k=0; k<cp[tmesh[eid].cnct[j]].face.size(); k++)
				{
					int fcid(cp[tmesh[eid].cnct[j]].face[k]);
					if(fcid!=eid && fcid!=fc_pre && fcid!=fc_next0)
					{
						for(uint k1=0; k1<tmesh[fcid].node.size(); k1++)
						{
							vector<int>::iterator it=find(tmesh[eid].IENtmp.begin(),tmesh[eid].IENtmp.end(),tmesh[fcid].node[k1]);
							if(it==tmesh[eid].IENtmp.end())
							{
								array<double,2> uv_ref={tmesh[eid].lcs[j].u[0],tmesh[eid].lcs[j].u[1]},uv;
								int rot;
								array<double,5> kui, kvi;
								FindRotateAndUVCoor(tmesh[eid].cnct[j],tmesh[eid].lcs[j].rot,uv_ref,fcid,tmesh[fcid].node[k1],rot,uv);
								FindLocalKnotVector(tmesh[fcid].node[k1],rot,uv,kui,kvi);
								tmesh[eid].IENtmp.push_back(tmesh[fcid].node[k1]);
								tmesh[eid].patch_kutmp.push_back(kui);
								tmesh[eid].patch_kvtmp.push_back(kvi);
							}
						}
					}
				}
			}
			//tmesh[eid].IEN.clear();
			//tmesh[eid].patch_ku.clear();
			//tmesh[eid].patch_kv.clear();
			//tmesh[eid].IEN=tmesh[eid].IENtmp;
			//tmesh[eid].patch_ku=tmesh[eid].patch_kutmp;
			//tmesh[eid].patch_kv=tmesh[eid].patch_kvtmp;
		}
	}
}

void TruncatedTspline_2D::FindNextRing(const vector<int>& pr0, const vector<int>& er0, vector<int>& pr1, vector<int>& er1, vector<int>& pr1_pref, vector<int>& pr1_eref)
{
	pr1.clear();
	er1.clear();
	pr1_pref.clear();
	pr1_eref.clear();
	for(uint i=0; i<pr0.size(); i++)
	{
		if(cp[pr0[i]].type!=2)
		{
			for(uint j=0; j<cp[pr0[i]].face.size(); j++)
			{
				int fc(cp[pr0[i]].face[j]);
				vector<int>::const_iterator it1=find(er0.begin(),er0.end(),fc);
				if(it1==er0.end())
				{
					er1.push_back(fc);
					for(uint k=0; k<tmesh[fc].node.size(); k++)
					{
						vector<int>::const_iterator it2=find(pr0.begin(),pr0.end(),tmesh[fc].node[k]);
						vector<int>::iterator it3=find(pr1.begin(),pr1.end(),tmesh[fc].node[k]);
						if(it2==pr0.end() && it3==pr1.end())
						{
							pr1.push_back(tmesh[fc].node[k]);
							pr1_pref.push_back(i);
							pr1_eref.push_back(fc);
						}
					}
				}
			}
		}
	}
}

void TruncatedTspline_2D::FindRotateAndUVCoor(int pref,int rot_ref,const array<double,2>& uv_ref,int eid,int pid,int& rot,array<double,2>& uv)//uv_ref is "global"
{
	vector<int>::iterator it0=find(tmesh[eid].node.begin(),tmesh[eid].node.end(),pref);
	vector<int>::iterator it1=find(tmesh[eid].node.begin(),tmesh[eid].node.end(),pid);
	int loc0=it0-tmesh[eid].node.begin();
	int loc1=it1-tmesh[eid].node.begin();
	int rot1=(tmesh[eid].lcs[loc1].rot+4-tmesh[eid].lcs[loc0].rot)%4;
	rot=(rot1+rot_ref)%4;
	double tmp[2]={tmesh[eid].lcs[loc1].u[0]-tmesh[eid].lcs[loc0].u[0],tmesh[eid].lcs[loc1].u[1]-tmesh[eid].lcs[loc0].u[1]};//direction vecton in element eid
	int rot2=(rot_ref-tmesh[eid].lcs[loc0].rot+4)%4;
	double tmp1[2]={tmp[0],tmp[1]};
	if(rot2==1)
	{
		tmp1[0]=-tmp[1]; tmp1[1]=tmp[0];
	}
	else if(rot2==2)
	{
		tmp1[0]=-tmp[0]; tmp1[1]=-tmp[1];
	}
	else if(rot2==3)
	{
		tmp1[0]=tmp[1]; tmp1[1]=-tmp[0];
	}
	uv[0]=uv_ref[0]+tmp1[0];
	uv[1]=uv_ref[1]+tmp1[1];
}

void TruncatedTspline_2D::FindLocalKnotVector(int id,int rot,const array<double,2>& uv,array<double,5>& ku,array<double,5>& kv)
{
	ku[2]=uv[0]; kv[2]=uv[1];
	if(rot==0)
	{
		for(int i=0; i<2; i++)
		{
			ku[i+3]=ku[i+2]+cp[id].kitvU[i+2];
			kv[i+3]=kv[i+2]+cp[id].kitvV[i+2];
			ku[1-i]=ku[2-i]-cp[id].kitvU[1-i];
			kv[1-i]=kv[2-i]-cp[id].kitvV[1-i];
		}
	}
	else if(rot==1)
	{
		for(int i=0; i<2; i++)
		{
			ku[i+3]=ku[i+2]+cp[id].kitvV[1-i];
			kv[i+3]=kv[i+2]+cp[id].kitvU[i+2];
			ku[1-i]=ku[2-i]-cp[id].kitvV[i+2];
			kv[1-i]=kv[2-i]-cp[id].kitvU[1-i];
		}
	}
	else if(rot==2)
	{
		for(int i=0; i<2; i++)
		{
			ku[i+3]=ku[i+2]+cp[id].kitvU[1-i];
			kv[i+3]=kv[i+2]+cp[id].kitvV[1-i];
			ku[1-i]=ku[2-i]-cp[id].kitvU[2+i];
			kv[1-i]=kv[2-i]-cp[id].kitvV[2+i];
		}
	}
	else if(rot==3)
	{
		for(int i=0; i<2; i++)
		{
			ku[i+3]=ku[i+2]+cp[id].kitvV[i+2];
			kv[i+3]=kv[i+2]+cp[id].kitvU[1-i];
			ku[1-i]=ku[2-i]-cp[id].kitvV[1-i];
			kv[1-i]=kv[2-i]-cp[id].kitvU[i+2];
		}
	}
	else
	{
		cerr<<"rot cannot be other numbers!\n";
		getchar();
	}
}

bool TruncatedTspline_2D::CheckSupport(const array<double,2>& u,const array<double,2>& v,const array<double,5>& ku,const array<double,5>& kv)
{
	if(ku[0]<u[1] && ku[4]>u[0] && kv[0]<v[1] && kv[4]>v[0])
	{
		return true;
	}
	else
	{
		return false;
	}
}

void TruncatedTspline_2D::SetLocalCoorSystem()
{
	for(uint eid=0; eid<tmesh.size(); eid++)
	{
		tmesh[eid].node.clear();
		tmesh[eid].lcs.clear();
		if(tmesh[eid].act==1)
		{
			double ul[2]={tmedge[tmesh[eid].edge[0]].len,tmedge[tmesh[eid].edge[1]].len};
			double uvcoor[8][2]={{0.,0.},{ul[0]/2.,0.},{ul[0],0.},{ul[0],ul[1]/2.},{ul[0],ul[1]},{ul[0]/2.,ul[1]},{0.,ul[1]},{0.,ul[1]/2.}};
			for(int i=0; i<4; i++)
			{
				int uved[2]={tmesh[eid].edge[i],tmesh[eid].edge[(i+3)%4]};
				for(int j=0; j<2; j++)
				{
					if(tmedge[uved[j]].act==0)
					{
						int edtmp(tmedge[uved[j]].chd[0]);
						if(tmedge[edtmp].pt[0]!= tmesh[eid].cnct[i] && tmedge[edtmp].pt[1]!= tmesh[eid].cnct[i])
						{
							edtmp=tmedge[uved[j]].chd[1];
						}
						uved[j]=edtmp;
					}
				}
				tmesh[eid].node.push_back(tmesh[eid].cnct[i]);
				ELCS2D tmp;
				tmp.u[0]=uvcoor[2*i][0]; tmp.u[1]=uvcoor[2*i][1];
				if(cp[tmesh[eid].cnct[i]].uved[0]==uved[0])
				{
					tmp.rot=i;
				}
				else if(cp[tmesh[eid].cnct[i]].uved[0]==uved[1])
				{
					tmp.rot=(i+1)%4;
				}
				else if(cp[tmesh[eid].cnct[i]].uved[1]==uved[0])
				{
					tmp.rot=(i+3)%4;
				}
				else
				{
					tmp.rot=(i+2)%4;
				}
				tmesh[eid].lcs.push_back(tmp);

				if(tmedge[tmesh[eid].edge[i]].act==0)
				{
					int ued(tmedge[tmesh[eid].edge[i]].chd[1]);
					if(tmedge[ued].pt[0]==tmesh[eid].cnct[i] || tmedge[ued].pt[1]==tmesh[eid].cnct[i])
					{
						ued=tmedge[tmesh[eid].edge[i]].chd[0];
					}
					tmesh[eid].node.push_back(tmedge[tmesh[eid].edge[i]].midpt);
					ELCS2D tmp1;
					tmp1.u[0]=uvcoor[2*i+1][0]; tmp1.u[1]=uvcoor[2*i+1][1];
					if(cp[tmedge[tmesh[eid].edge[i]].midpt].uved[1]==ued)
					{
						tmp1.rot=(i+3)%4;
					}
					else
					{
						tmp1.rot=(i+2)%4;
					}
					tmesh[eid].lcs.push_back(tmp1);
				}
			}
		}
	}
}

void TruncatedTspline_2D::UpdatePatchCP_Unstruct(int eid)
{
	vector<int> IEN_new_all;
	vector<array<double,5>> KU_new_all;
	vector<array<double,5>> KV_new_all;
	for(uint i=0; i<tmesh[eid].chd.size(); i++)
	{
		int cid(tmesh[eid].chd[i]);
		for(uint j=0; j<tmesh[cid].IENtmp.size(); j++)
		{
			vector<int>::iterator it=find(IEN_new_all.begin(),IEN_new_all.end(),tmesh[cid].IENtmp[j]);
			if(it==IEN_new_all.end())
			{
				IEN_new_all.push_back(tmesh[cid].IENtmp[j]);
				array<double,5> kutmp,kvtmp;
				for(int k=0; k<5; k++)
				{
					kutmp[k]=tmesh[cid].patch_kutmp[j][k]+tmesh[eid].chd_o[i][0];
					kvtmp[k]=tmesh[cid].patch_kvtmp[j][k]+tmesh[eid].chd_o[i][1];
				}
				KU_new_all.push_back(kutmp);
				KV_new_all.push_back(kvtmp);
			}
		}
	}

	vector<int> IEN_old(tmesh[eid].IEN);
	vector<vector<double>> cmat(IEN_new_all.size(),vector<double>(IEN_old.size(),0.));
	//vector<vector<double>> tmat(IEN_new_all.size(),vector<double>(IEN_old.size(),0.));
	for(uint i=0; i<IEN_old.size(); i++)
	{
		vector<int>::iterator it=find(IEN_new_all.begin(),IEN_new_all.end(),IEN_old[i]);
		if(it!=IEN_new_all.end())
		{
			int loc(it-IEN_new_all.begin());
			if(tmesh[eid].patch_ku[i]!=KU_new_all[loc] || tmesh[eid].patch_kv[i]!=KV_new_all[loc])
			{
				cp[IEN_old[i]].aff=1;
			}
			else
			{
				cmat[loc][i]=1.;
			}
		}
	}

	for(uint i=0; i<IEN_new_all.size(); i++)
	{
		for(uint j=0; j<IEN_old.size(); j++)
		{
			if(CheckSubKnotVector(KU_new_all[i],KV_new_all[i],tmesh[eid].patch_ku[j],tmesh[eid].patch_kv[j]))
			{
				vector<double> ku1(10),kv1(10);
				vector<vector<double>> Tu,Tv;
				vector<double> ku(tmesh[eid].patch_ku[j].begin(),tmesh[eid].patch_ku[j].end());
				vector<double> kv(tmesh[eid].patch_kv[j].begin(),tmesh[eid].patch_kv[j].end());
				vector<double>::iterator it1,it2;
				it1=set_union(tmesh[eid].patch_ku[j].begin(),tmesh[eid].patch_ku[j].end(),KU_new_all[i].begin(),KU_new_all[i].end(),ku1.begin());
				it2=set_union(tmesh[eid].patch_kv[j].begin(),tmesh[eid].patch_kv[j].end(),KV_new_all[i].begin(),KV_new_all[i].end(),kv1.begin());
				ku1.resize(it1-ku1.begin());
				kv1.resize(it2-kv1.begin());
				TMatrix(ku,ku1,3,Tu);
				TMatrix(kv,kv1,3,Tv);
				it1=search(ku1.begin(),ku1.end(),KU_new_all[i].begin(),KU_new_all[i].end());
				it2=search(kv1.begin(),kv1.end(),KV_new_all[i].begin(),KV_new_all[i].end());
				if(it1!=ku1.end() && it2!=kv1.end())
				{
					int loc1=it1-ku1.begin();
					int loc2=it2-kv1.begin();
					cmat[i][j]=Tu[loc1][0]*Tv[loc2][0];
				}
			}
		}
	}
	for(uint i=0; i<IEN_new_all.size(); i++)
	{
		//cp[IEN_new_all[i]].coortmp[0]=0.; cp[IEN_new_all[i]].coortmp[1]=0.; cp[IEN_new_all[i]].coortmp[2]=0.;
		if(cp[IEN_new_all[i]].update==0)
		{
			for(uint j=0; j<IEN_old.size(); j++)
			{
				if(cmat[i][j]!=0.)
				{
					cp[IEN_new_all[i]].coortmp[0]+=cmat[i][j]*cp[IEN_old[j]].coor[0];
					cp[IEN_new_all[i]].coortmp[1]+=cmat[i][j]*cp[IEN_old[j]].coor[1];
					cp[IEN_new_all[i]].coortmp[2]+=cmat[i][j]*cp[IEN_old[j]].coor[2];
				}
			}
			cp[IEN_new_all[i]].update=1;
		}
	}
	//check truncation
	/*double sp0[4][2]={{0.,0.},{tmedge[tmesh[eid].edge[0]].len,0.},{tmedge[tmesh[eid].edge[0]].len,tmedge[tmesh[eid].edge[1]].len},{0.,tmedge[tmesh[eid].edge[1]].len}};
	vector<array<double,2>> spt(9);
	for(int i=0; i<4; i++)
	{
		spt[i][0]=sp0[i][0]; spt[i][1]=sp0[i][1];
		spt[i+4][0]=(sp0[i][0]+sp0[(i+1)%4][0])/2.; spt[i+4][1]=(sp0[i][1]+sp0[(i+1)%4][1])/2.; 
	}
	spt[8][0]=spt[4][0]; spt[8][1]=spt[5][1];
	for(uint i=0; i<IEN_old.size(); i++)
	{
		vector<double> coef;
		for(uint j=0; j<IEN_new_all.size(); j++)
		{
			coef.push_back(cmat[j][i]);
		}
		if(!CheckFullChildren_1(spt,tmesh[eid].patch_ku[i],tmesh[eid].patch_kv[i],KU_new_all,KV_new_all,coef))
		{
			cp[IEN_old[i]].truntmp=1;
			for(uint j=0; j<IEN_new_all.size(); j++)
			{
				if(cmat[j][i]!=0. && IEN_old[i]!=IEN_new_all[j])
				{
					vector<int>::iterator it=find(cp[IEN_old[i]].tbftmp.begin(),cp[IEN_old[i]].tbftmp.end(),IEN_new_all[j]);
					if(it==cp[IEN_old[i]].tbftmp.end())
					{
						cp[IEN_old[i]].tbftmp.push_back(IEN_new_all[j]);
						cp[IEN_old[i]].tctmp.push_back(cmat[j][i]);
					}
				}
			}
		}
	}*/
}

void TruncatedTspline_2D::Truncation()
{
	for(uint eid=0; eid<tmesh.size(); eid++)
	{
		if(tmesh[eid].act==1 && (tmesh[eid].type==0 || tmesh[eid].type==1))
		{
			for(uint i=0; i<tmesh[eid].IEN.size(); i++)
			{
				for(uint j=i+1; j<tmesh[eid].IEN.size(); j++)
				{
					if(CheckSubKnotVector(tmesh[eid].patch_ku[j],tmesh[eid].patch_kv[j],tmesh[eid].patch_ku[i],tmesh[eid].patch_kv[i]))//i to be truncated, j to be discarded children
					{
						int Bt(tmesh[eid].IEN[i]), Bc(tmesh[eid].IEN[j]);
						vector<int>::iterator it=find(cp[Bt].tbf.begin(),cp[Bt].tbf.end(),Bc);
						if(it==cp[Bt].tbf.end())
						{
							vector<double> ku1(10),kv1(10);
							vector<vector<double>> Tu,Tv;
							vector<double> ku(tmesh[eid].patch_ku[i].begin(),tmesh[eid].patch_ku[i].end());
							vector<double> kv(tmesh[eid].patch_kv[i].begin(),tmesh[eid].patch_kv[i].end());
							vector<double>::iterator it1,it2;
							it1=set_union(tmesh[eid].patch_ku[i].begin(),tmesh[eid].patch_ku[i].end(),tmesh[eid].patch_ku[j].begin(),tmesh[eid].patch_ku[j].end(),ku1.begin());
							it2=set_union(tmesh[eid].patch_kv[i].begin(),tmesh[eid].patch_kv[i].end(),tmesh[eid].patch_kv[j].begin(),tmesh[eid].patch_kv[j].end(),kv1.begin());
							ku1.resize(it1-ku1.begin());
							kv1.resize(it2-kv1.begin());
							TMatrix(ku,ku1,3,Tu);
							TMatrix(kv,kv1,3,Tv);
							it1=search(ku1.begin(),ku1.end(),tmesh[eid].patch_ku[j].begin(),tmesh[eid].patch_ku[j].end());
							it2=search(kv1.begin(),kv1.end(),tmesh[eid].patch_kv[j].begin(),tmesh[eid].patch_kv[j].end());
							if(it1!=ku1.end() && it2!=kv1.end())
							{
								int loc1=it1-ku1.begin();
								int loc2=it2-kv1.begin();
								cp[Bt].trun=1;
								cp[Bt].tbf.push_back(Bc);
								cp[Bt].tc.push_back(Tu[loc1][0]*Tv[loc2][0]);
							}
						}
					}
					else if(CheckSubKnotVector(tmesh[eid].patch_ku[i],tmesh[eid].patch_kv[i],tmesh[eid].patch_ku[j],tmesh[eid].patch_kv[j]))//j to be truncated, i to be discarded children
					{
						int Bt(tmesh[eid].IEN[j]), Bc(tmesh[eid].IEN[i]);
						vector<int>::iterator it=find(cp[Bt].tbf.begin(),cp[Bt].tbf.end(),Bc);
						if(it==cp[Bt].tbf.end())
						{
							vector<double> ku1(10),kv1(10);
							vector<vector<double>> Tu,Tv;
							vector<double> ku(tmesh[eid].patch_ku[j].begin(),tmesh[eid].patch_ku[j].end());
							vector<double> kv(tmesh[eid].patch_kv[j].begin(),tmesh[eid].patch_kv[j].end());
							vector<double>::iterator it1,it2;
							it1=set_union(tmesh[eid].patch_ku[j].begin(),tmesh[eid].patch_ku[j].end(),tmesh[eid].patch_ku[i].begin(),tmesh[eid].patch_ku[i].end(),ku1.begin());
							it2=set_union(tmesh[eid].patch_kv[j].begin(),tmesh[eid].patch_kv[j].end(),tmesh[eid].patch_kv[i].begin(),tmesh[eid].patch_kv[i].end(),kv1.begin());
							ku1.resize(it1-ku1.begin());
							kv1.resize(it2-kv1.begin());
							TMatrix(ku,ku1,3,Tu);
							TMatrix(kv,kv1,3,Tv);
							it1=search(ku1.begin(),ku1.end(),tmesh[eid].patch_ku[i].begin(),tmesh[eid].patch_ku[i].end());
							it2=search(kv1.begin(),kv1.end(),tmesh[eid].patch_kv[i].begin(),tmesh[eid].patch_kv[i].end());
							if(it1!=ku1.end() && it2!=kv1.end())
							{
								int loc1=it1-ku1.begin();
								int loc2=it2-kv1.begin();
								cp[Bt].trun=1;
								cp[Bt].tbf.push_back(Bc);
								cp[Bt].tc.push_back(Tu[loc1][0]*Tv[loc2][0]);
							}
						}
					}
				}
			}
		}
	}
}

void TruncatedTspline_2D::ElementBasis(int eid, double u, double v, vector<double>& Nt, vector<array<double,2>>& dNdt)
{
	if(tmesh[eid].act==1)
	{
		if(tmesh[eid].type==0 || tmesh[eid].type==1)
		{
			ElementBasis_Regular(eid,u,v,Nt,dNdt);
		}
		else if(tmesh[eid].type==4)
		{
			ElementBasis_Irregular(eid,u,v,Nt,dNdt);
		}
	}
}

void TruncatedTspline_2D::ElementBasis_Regular(int eid, double u, double v, vector<double>& Nt, vector<array<double,2>>& dNdt)
{
	Nt.clear();
	dNdt.clear();
	Nt.resize(tmesh[eid].IEN.size());
	dNdt.resize(tmesh[eid].IEN.size());
	vector<double> Nt0(tmesh[eid].IEN.size());
	vector<array<double,2>> dNdt0(tmesh[eid].IEN.size());
	vector<double> ku(5,0.),kv(5,0.),uval,vval;
	BSplineBasis bu,bv;
	for(uint i=0;i<tmesh[eid].IEN.size();i++)
	{
		ku.assign(tmesh[eid].patch_ku[i].begin(),tmesh[eid].patch_ku[i].end());
		kv.assign(tmesh[eid].patch_kv[i].begin(),tmesh[eid].patch_kv[i].end());
		bu.Set(3,ku);
		bv.Set(3,kv);
		bu.BasisFunction(0,u,1,uval);
		bv.BasisFunction(0,v,1,vval);
		Nt0[i]=uval[0]*vval[0];
		dNdt0[i][0]=uval[1]*vval[0];
		dNdt0[i][1]=uval[0]*vval[1];
	}
	for(uint i=0;i<tmesh[eid].IEN.size();i++)
	{
		Nt[i]=Nt0[i];
		dNdt[i][0]=dNdt0[i][0];
		dNdt[i][1]=dNdt0[i][1];
		if(cp[tmesh[eid].IEN[i]].trun==1)
		{
			int pid(tmesh[eid].IEN[i]);
			for(uint j=0; j<cp[pid].tbf.size(); j++)
			{
				vector<int>::iterator it=find(tmesh[eid].IEN.begin(),tmesh[eid].IEN.end(),cp[pid].tbf[j]);
				if(it!=tmesh[eid].IEN.end())
				{
					int loc(it-tmesh[eid].IEN.begin());
					Nt[i]-=cp[pid].tc[j]*Nt0[loc];
					dNdt[i][0]-=cp[pid].tc[j]*dNdt0[loc][0];
					dNdt[i][1]-=cp[pid].tc[j]*dNdt0[loc][1];
				}
			}
		}
	}
}

void TruncatedTspline_2D::ElementBasis_Irregular(int eid, double u, double v, vector<double>& Nt, vector<array<double,2>>& dNdt)
{
	//need to change variabe for Bezier, (u,v) -> (u_b,v_b)
	Nt.clear();
	dNdt.clear();
	Nt.resize(tmesh[eid].IEN.size(),0.);
	dNdt.resize(tmesh[eid].IEN.size());
	vector<double> Nt1(tmesh[eid].IEN.size(),0.);
	vector<array<double,2>> dNdt1(tmesh[eid].IEN.size());
	uint nv(cp[tmesh[eid].cnct[0]].face.size());
	vector<vector<double>> bmat(tmesh[eid].bemat.size(),vector<double>(tmesh[eid].bemat[0].size()));
	for(uint i=0; i<tmesh[eid].bemat.size(); i++)
	{
		for(uint j=0; j<tmesh[eid].bemat[i].size(); j++)
		{
			bmat[i][j]=tmesh[eid].bemat[i][j];
		}
	}
	//SetBezier4TranMatOP(nv,bmat);
	BezierElement2D be(4);
	vector<double> Nt0;
	vector<array<double,2>> dNdt0;
	double u_b(u/tmedge[tmesh[eid].edge[0]].len), v_b(v/tmedge[tmesh[eid].edge[3]].len);
	be.Basis(u_b,v_b,Nt0,dNdt0);
	for(uint i=0; i<2*nv+1; i++)
	{
		dNdt1[i][0]=0.; dNdt1[i][1]=0.;
		for(int j=0; j<25; j++)
		{
			Nt1[i]+=bmat[i][j]*Nt0[j];
			dNdt1[i][0]+=bmat[i][j]*dNdt0[j][0];
			dNdt1[i][1]+=bmat[i][j]*dNdt0[j][1];
		}
		dNdt1[i][0]/=tmedge[tmesh[eid].edge[0]].len;
		dNdt1[i][1]/=tmedge[tmesh[eid].edge[3]].len;
	}
	vector<double> ku(5,0.),kv(5,0.),uval,vval;
	BSplineBasis bu,bv;
	for(uint i=2*nv+1; i<tmesh[eid].IEN.size(); i++)
	{
		ku.assign(tmesh[eid].patch_ku[i-(2*nv+1)].begin(),tmesh[eid].patch_ku[i-(2*nv+1)].end());
		kv.assign(tmesh[eid].patch_kv[i-(2*nv+1)].begin(),tmesh[eid].patch_kv[i-(2*nv+1)].end());
		bu.Set(3,ku);
		bv.Set(3,kv);
		bu.BasisFunction(0,u,1,uval);
		bv.BasisFunction(0,v,1,vval);
		Nt1[i]=uval[0]*vval[0];
		dNdt1[i][0]=uval[1]*vval[0];
		dNdt1[i][1]=uval[0]*vval[1];
	}
	for(uint i=0;i<tmesh[eid].IEN.size();i++)
	{
		Nt[i]=Nt1[i];
		dNdt[i][0]=dNdt1[i][0];
		dNdt[i][1]=dNdt1[i][1];
		if(cp[tmesh[eid].IEN[i]].trun==1)
		{
			int pid(tmesh[eid].IEN[i]);
			for(uint j=0; j<cp[pid].tbf.size(); j++)
			{
				vector<int>::iterator it=find(tmesh[eid].IEN.begin(),tmesh[eid].IEN.end(),cp[pid].tbf[j]);
				if(it!=tmesh[eid].IEN.end())
				{
					int loc(it-tmesh[eid].IEN.begin());
					Nt[i]-=cp[pid].tc[j]*Nt1[loc];
					dNdt[i][0]-=cp[pid].tc[j]*dNdt1[loc][0];
					dNdt[i][1]-=cp[pid].tc[j]*dNdt1[loc][1];
				}
			}
		}
	}
}

void TruncatedTspline_2D::SetBezierMatIrrPatch()
{
	for(uint eid=0; eid<tmesh.size(); eid++)
	{
		if(tmesh[eid].act==1 && tmesh[eid].type==4)
		{
			tmesh[eid].bemat.clear();
			uint nv(cp[tmesh[eid].cnct[0]].face.size());
			//vector<vector<double>> bmat;
			SetBezier4TranMatOP(nv,tmesh[eid].bemat);
			//SetBezier4TranMat(nv,tmesh[eid].bemat);
			vector<array<double,5>> patch_ku, patch_kv;
			FindPatchKnotVector_Irr(eid,patch_ku,patch_kv);
			vector<vector<double>> coef(patch_ku.size(),vector<double>(16,0.));
			vector<vector<double>> coef1(patch_ku.size(),vector<double>(25,0.));
			array<double,2> ktsU={0.,tmedge[tmesh[eid].edge[0]].len};
			array<double,2> ktsV={0.,tmedge[tmesh[eid].edge[3]].len};
			double bzku[6]={ktsU[0],ktsU[0],ktsU[0],ktsU[1],ktsU[1],ktsU[1]};
			double bzkv[6]={ktsV[0],ktsV[0],ktsV[0],ktsV[1],ktsV[1],ktsV[1]};
			for(uint i=0; i<patch_ku.size(); i++)
			{
				if(patch_ku[i][0]<ktsU[1] && patch_ku[i][4]>ktsU[0] && patch_kv[i][0]<ktsV[1] && patch_kv[i][4]>ktsV[0])
				{
					vector<double> ku(patch_ku[i].begin(),patch_ku[i].end());
					vector<double> kv(patch_kv[i].begin(),patch_kv[i].end());
					vector<double> ku1,kv1;
					vector<vector<double>> Tu,Tv;
					BezierInsertKnots(ku,ktsU,ku1);
					BezierInsertKnots(kv,ktsV,kv1);
					TMatrix(ku,ku1,3,Tu);
					TMatrix(kv,kv1,3,Tv);
					vector<double>::iterator it1=search(ku1.begin(),ku1.end(),bzku,bzku+6);
					vector<double>::iterator it2=search(kv1.begin(),kv1.end(),bzkv,bzkv+6);
					//if(it1!=ku1.end() && it2!=kv1.end())
					//{
					int loc1(it1-ku1.begin()-1),loc2(it2-kv1.begin()-1),count(0);
					for(int j=0;j<4;j++)
					{
						for(int k=0;k<4;k++)
						{
							coef[i][count]=Tu[loc1+k][0]*Tv[loc2+j][0];
							count++;
						}
					}
					//}
				}
			}
			vector<vector<double>> demat;
			DegreeElevate(demat);
			for(uint i=0; i<patch_ku.size(); i++)
			{
				for(int j=0; j<16; j++)
				{
					for(int k=0; k<25; k++)
					{
						coef1[i][k]+=coef[i][j]*demat[k][j];
					}
				}
			}
			int ploc[5]={2,3,4,5,6};
			int bzloc[16]={3,4,8,9,13,14,15,16,17,18,19,20,21,22,23,24};
			for(int i=0; i<5; i++)
			{
				if(cp[tmesh[eid].IEN[ploc[i]]].trun==0)
				{
					for(int j=0; j<16; j++)
					{
						tmesh[eid].bemat[ploc[i]][bzloc[j]]=coef1[i][bzloc[j]];
					}
				}
			}
		}
	}
}

void TruncatedTspline_2D::FindPatchKnotVector_Irr(int eid, vector<array<double,5>>& patch_ku, vector<array<double,5>>& patch_kv)
{
	patch_ku.clear();
	patch_kv.clear();
	patch_ku.resize(5);
	patch_kv.resize(5);
	//find knot vectors for node 2 to 6
	if(tmesh[eid].act==1 && tmesh[eid].type==4)
	{
		//first for node 2
		array<double,2> uv_ref0={tmesh[eid].lcs[3].u[0],tmesh[eid].lcs[3].u[1]},uv0;
		int fcid0(tmedge[tmesh[eid].edge[3]].face[0]), rot0;
		if(fcid0==eid) fcid0=tmedge[tmesh[eid].edge[3]].face[1];
		FindRotateAndUVCoor(tmesh[eid].cnct[3],tmesh[eid].lcs[3].rot,uv_ref0,fcid0,tmesh[fcid0].cnct[2],rot0,uv0);
		FindLocalKnotVector(tmesh[fcid0].cnct[2],rot0,uv0,patch_ku[0],patch_kv[0]);
		//node 3, 4, 5
		for(int i=3; i>0; i--)
		{
			array<double,2> uv1={tmesh[eid].lcs[i].u[0],tmesh[eid].lcs[i].u[1]};
			FindLocalKnotVector(tmesh[eid].cnct[i],tmesh[eid].lcs[i].rot,uv1,patch_ku[4-i],patch_kv[4-i]);
		}
		//node 6
		array<double,2> uv_ref2={tmesh[eid].lcs[1].u[0],tmesh[eid].lcs[1].u[1]},uv2;
		int fcid2(tmedge[tmesh[eid].edge[0]].face[0]), rot2;
		if(fcid2==eid) fcid2=tmedge[tmesh[eid].edge[0]].face[1];
		FindRotateAndUVCoor(tmesh[eid].cnct[1],tmesh[eid].lcs[1].rot,uv_ref2,fcid2,tmesh[fcid2].cnct[2],rot2,uv2);
		FindLocalKnotVector(tmesh[fcid2].cnct[2],rot2,uv2,patch_ku[4],patch_kv[4]);
	}
}

void TruncatedTspline_2D::SurfacePointMap(int eid, double u, double v, array<double,3>& pt, array<double,3>& norm)
{
	vector<double> Nt;
	vector<array<double,2>> dNdt;
	ElementBasis(eid,u,v,Nt,dNdt);
	pt[0]=0.; pt[1]=0.; pt[2]=0.;
	double nmtmp[2][3]={{0.,0.,0.},{0.,0.,0.}};
	for(uint i=0; i<tmesh[eid].IEN.size(); i++)
	{
		pt[0]+=cp[tmesh[eid].IEN[i]].coor[0]*Nt[i];
		pt[1]+=cp[tmesh[eid].IEN[i]].coor[1]*Nt[i];
		pt[2]+=cp[tmesh[eid].IEN[i]].coor[2]*Nt[i];
		nmtmp[0][0]+=cp[tmesh[eid].IEN[i]].coor[0]*dNdt[i][0];
		nmtmp[0][1]+=cp[tmesh[eid].IEN[i]].coor[1]*dNdt[i][0];
		nmtmp[0][2]+=cp[tmesh[eid].IEN[i]].coor[2]*dNdt[i][0];
		nmtmp[1][0]+=cp[tmesh[eid].IEN[i]].coor[0]*dNdt[i][1];
		nmtmp[1][1]+=cp[tmesh[eid].IEN[i]].coor[1]*dNdt[i][1];
		nmtmp[1][2]+=cp[tmesh[eid].IEN[i]].coor[2]*dNdt[i][1];
	}
	norm[0]=nmtmp[0][1]*nmtmp[1][2]-nmtmp[0][2]*nmtmp[1][1];
	norm[1]=nmtmp[0][2]*nmtmp[1][0]-nmtmp[0][0]*nmtmp[1][2];
	norm[2]=nmtmp[0][0]*nmtmp[1][1]-nmtmp[0][1]*nmtmp[1][0];
	double len=sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
	norm[0]/=len; norm[1]/=len; norm[2]/=len;
}

void TruncatedTspline_2D::ElementRefine_Unstruct_4(int eid)
{
	int pid[5];
	int edid[12];
	Vertex2D ptmp1;
	//ptmp1.coor[0]=(cp[tmesh[eid].cnct[0]].coor[0]+cp[tmesh[eid].cnct[1]].coor[0]+cp[tmesh[eid].cnct[2]].coor[0]+cp[tmesh[eid].cnct[3]].coor[0])/4.;
	//ptmp1.coor[1]=(cp[tmesh[eid].cnct[0]].coor[1]+cp[tmesh[eid].cnct[1]].coor[1]+cp[tmesh[eid].cnct[2]].coor[1]+cp[tmesh[eid].cnct[3]].coor[1])/4.;
	//ptmp1.coor[2]=(cp[tmesh[eid].cnct[0]].coor[2]+cp[tmesh[eid].cnct[1]].coor[2]+cp[tmesh[eid].cnct[2]].coor[2]+cp[tmesh[eid].cnct[3]].coor[2])/4.;
	cp.push_back(ptmp1);
	pid[0]=cp.size()-1;
	for(int j=0; j<4; j++)
	{
		int itmp[2]={tmesh[eid].cnct[j],tmesh[eid].cnct[(j+1)%4]};
		if(tmedge[tmesh[eid].edge[j]].act==1)
		{
			Vertex2D ptmp;
			//ptmp.coor[0]=(cp[itmp[0]].coor[0]+cp[itmp[1]].coor[0])/2.;
			//ptmp.coor[1]=(cp[itmp[0]].coor[1]+cp[itmp[1]].coor[1])/2.;
			//ptmp.coor[2]=(cp[itmp[0]].coor[2]+cp[itmp[1]].coor[2])/2.;
			cp.push_back(ptmp);
			pid[j+1]=cp.size()-1;
			tmedge[tmesh[eid].edge[j]].midpt=pid[j+1];
			Edge edtmp1,edtmp2,edtmp3;
			edtmp1.act=1;
			edtmp1.lev=tmedge[tmesh[eid].edge[j]].lev+1;
			edtmp1.pt[0]=itmp[0]; edtmp1.pt[1]=pid[j+1];
			edtmp1.len=tmedge[tmesh[eid].edge[j]].len/2.;
			edtmp1.prt=tmesh[eid].edge[j];
			edtmp2.act=1;
			edtmp2.lev=tmedge[tmesh[eid].edge[j]].lev+1;
			edtmp2.pt[0]=pid[j+1]; edtmp2.pt[1]=itmp[1];
			edtmp2.len=tmedge[tmesh[eid].edge[j]].len/2.;
			edtmp2.prt=tmesh[eid].edge[j];
			edtmp3.act=1;
			edtmp3.lev=tmedge[tmesh[eid].edge[(j+1)%4]].lev+1;
			edtmp3.pt[0]=pid[0]; edtmp3.pt[1]=pid[j+1];
			edtmp3.len=tmedge[tmesh[eid].edge[(j+1)%4]].len/2.;
			tmedge.push_back(edtmp1);
			edid[3*j]=tmedge.size()-1;
			tmedge[tmesh[eid].edge[j]].chd[0]=edid[3*j];
			tmedge.push_back(edtmp2);
			edid[3*j+1]=tmedge.size()-1;
			tmedge[tmesh[eid].edge[j]].chd[1]=edid[3*j+1];
			tmedge.push_back(edtmp3);
			edid[3*j+2]=tmedge.size()-1;
			tmedge[tmesh[eid].edge[j]].act=0;
		}
		else
		{
			pid[j+1]=tmedge[tmesh[eid].edge[j]].midpt;
			Edge edtmp3;
			edtmp3.act=1;
			edtmp3.lev=tmedge[tmesh[eid].edge[(j+1)%4]].lev+1;
			edtmp3.pt[0]=pid[0]; edtmp3.pt[1]=pid[j+1];
			edtmp3.len=tmedge[tmesh[eid].edge[(j+1)%4]].len/2.;
			tmedge.push_back(edtmp3);
			edid[3*j+2]=tmedge.size()-1;
			int ied(tmedge[tmesh[eid].edge[j]].chd[0]);
			if(tmedge[ied].pt[0]==itmp[0] || tmedge[ied].pt[1]==itmp[0])
			{
				edid[3*j]=ied;
				edid[3*j+1]=tmedge[tmesh[eid].edge[j]].chd[1];
			}
			else
			{
				edid[3*j]=tmedge[tmesh[eid].edge[j]].chd[1];
				edid[3*j+1]=ied;
			}
		}
	}

	int e_cnct[4][4]={{tmesh[eid].cnct[0],pid[1],pid[0],pid[4]},{pid[1],tmesh[eid].cnct[1],pid[2],pid[0]},{pid[0],pid[2],tmesh[eid].cnct[2],pid[3]},{pid[4],pid[0],pid[3],tmesh[eid].cnct[3]}};
	int e_edge[4][4]={{edid[0],edid[2],edid[11],edid[10]},{edid[1],edid[3],edid[5],edid[2]},{edid[5],edid[4],edid[6],edid[8]},{edid[11],edid[8],edid[7],edid[9]}};
	int enewid[4];
	double chd_org[4][2]={{0.,0.},{tmedge[tmesh[eid].edge[0]].len/2.,0.},{tmedge[tmesh[eid].edge[0]].len/2.,tmedge[tmesh[eid].edge[1]].len/2.},{0.,tmedge[tmesh[eid].edge[1]].len/2.}};
	vector<Element2D> etmp(4);
	for(int i=0; i<4; i++)
	{
		etmp[i].act=1;
		etmp[i].type=0;
		etmp[i].lev=tmesh[eid].lev+1;
		etmp[i].prt=eid;
		for(int j=0; j<4; j++)
		{
			etmp[i].cnct[j]=e_cnct[i][j];
			etmp[i].edge[j]=e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i]=tmesh.size()-1;
		tmesh[eid].chd.push_back(enewid[i]);
		array<double,2> chd_o_tmp={chd_org[i][0],chd_org[i][1]};
		tmesh[eid].chd_o.push_back(chd_o_tmp);
	}

	tmesh[eid].act=0;
}

void TruncatedTspline_2D::ElementRefine_Irregular_4(int eid)
{
	int pid[5];
	int edid[12];
	Vertex2D ptmp1;
	ptmp1.update=1;
	ptmp1.coortmp[0]=(cp[tmesh[eid].cnct[0]].coor[0]+cp[tmesh[eid].cnct[1]].coor[0]+cp[tmesh[eid].cnct[2]].coor[0]+cp[tmesh[eid].cnct[3]].coor[0])/4.;
	ptmp1.coortmp[1]=(cp[tmesh[eid].cnct[0]].coor[1]+cp[tmesh[eid].cnct[1]].coor[1]+cp[tmesh[eid].cnct[2]].coor[1]+cp[tmesh[eid].cnct[3]].coor[1])/4.;
	ptmp1.coortmp[2]=(cp[tmesh[eid].cnct[0]].coor[2]+cp[tmesh[eid].cnct[1]].coor[2]+cp[tmesh[eid].cnct[2]].coor[2]+cp[tmesh[eid].cnct[3]].coor[2])/4.;
	cp.push_back(ptmp1);
	pid[0]=cp.size()-1;
	for(int j=0; j<4; j++)
	{
		int itmp[2]={tmesh[eid].cnct[j],tmesh[eid].cnct[(j+1)%4]};
		if(tmedge[tmesh[eid].edge[j]].act==1)
		{
			Vertex2D ptmp;
			int id0[4]={tmesh[eid].cnct[j],tmesh[eid].cnct[(j+1)%4],tmesh[eid].cnct[(j+2)%4],tmesh[eid].cnct[(j+3)%4]};
			int fcnb(tmedge[tmesh[eid].edge[j]].face[0]);
			if(fcnb==eid) fcnb=tmedge[tmesh[eid].edge[j]].face[1];
			int* it1=find(tmesh[fcnb].cnct,tmesh[fcnb].cnct+4,tmesh[eid].cnct[j]);
			int loc1(it1-tmesh[fcnb].cnct);
			int id1[2]={tmesh[fcnb].cnct[(loc1+2)%4],tmesh[fcnb].cnct[(loc1+3)%4]};
			//if(j==0 || j==3)
			//{
				ptmp.coortmp[0]=(6.*cp[id0[0]].coor[0]+6.*cp[id0[1]].coor[0]+cp[id0[2]].coor[0]+cp[id0[3]].coor[0]+cp[id1[0]].coor[0]+cp[id1[1]].coor[0])/16.;
				ptmp.coortmp[1]=(6.*cp[id0[0]].coor[1]+6.*cp[id0[1]].coor[1]+cp[id0[2]].coor[1]+cp[id0[3]].coor[1]+cp[id1[0]].coor[1]+cp[id1[1]].coor[1])/16.;
				ptmp.coortmp[2]=(6.*cp[id0[0]].coor[2]+6.*cp[id0[1]].coor[2]+cp[id0[2]].coor[2]+cp[id0[3]].coor[2]+cp[id1[0]].coor[2]+cp[id1[1]].coor[2])/16.;
			//}
			//else
			//{
			//	ptmp.coor[0]=cp[tmesh[eid].cnct[0]].coor[0]/16.;
			//	ptmp.coor[1]=cp[tmesh[eid].cnct[0]].coor[1]/16.;
			//	ptmp.coor[2]=cp[tmesh[eid].cnct[0]].coor[2]/16.;
			//}
			cp.push_back(ptmp);
			pid[j+1]=cp.size()-1;
			tmedge[tmesh[eid].edge[j]].midpt=pid[j+1];
			Edge edtmp1,edtmp2,edtmp3;
			edtmp1.act=1;
			edtmp1.lev=tmedge[tmesh[eid].edge[j]].lev+1;
			edtmp1.pt[0]=itmp[0]; edtmp1.pt[1]=pid[j+1];
			edtmp1.len=tmedge[tmesh[eid].edge[j]].len/2.;
			edtmp1.prt=tmesh[eid].edge[j];
			edtmp2.act=1;
			edtmp2.lev=tmedge[tmesh[eid].edge[j]].lev+1;
			edtmp2.pt[0]=pid[j+1]; edtmp2.pt[1]=itmp[1];
			edtmp2.len=tmedge[tmesh[eid].edge[j]].len/2.;
			edtmp2.prt=tmesh[eid].edge[j];
			edtmp3.act=1;
			edtmp3.lev=tmedge[tmesh[eid].edge[(j+1)%4]].lev+1;
			edtmp3.pt[0]=pid[0]; edtmp3.pt[1]=pid[j+1];
			edtmp3.len=tmedge[tmesh[eid].edge[(j+1)%4]].len/2.;
			tmedge.push_back(edtmp1);
			edid[3*j]=tmedge.size()-1;
			tmedge[tmesh[eid].edge[j]].chd[0]=edid[3*j];
			tmedge.push_back(edtmp2);
			edid[3*j+1]=tmedge.size()-1;
			tmedge[tmesh[eid].edge[j]].chd[1]=edid[3*j+1];
			tmedge.push_back(edtmp3);
			edid[3*j+2]=tmedge.size()-1;
			tmedge[tmesh[eid].edge[j]].act=0;
		}
		else
		{
			pid[j+1]=tmedge[tmesh[eid].edge[j]].midpt;
			int id1[2]={tmesh[eid].cnct[(j+2)%4],tmesh[eid].cnct[(j+3)%4]};
			//cp[pid[j+1]].update=1;
			//if(j==0 || j==3)
			//{
			//	cp[pid[j+1]].coor[0]+=(cp[id1[0]].coor[0]+cp[id1[1]].coor[0])/16.;
			//	cp[pid[j+1]].coor[1]+=(cp[id1[0]].coor[1]+cp[id1[1]].coor[1])/16.;
			//	cp[pid[j+1]].coor[2]+=(cp[id1[0]].coor[2]+cp[id1[1]].coor[2])/16.;
			//}
			//else
			//{
			//	cp[pid[j+1]].coor[0]+=cp[tmesh[eid].cnct[0]].coor[0]/16.;
			//	cp[pid[j+1]].coor[1]+=cp[tmesh[eid].cnct[0]].coor[1]/16.;
			//	cp[pid[j+1]].coor[2]+=cp[tmesh[eid].cnct[0]].coor[2]/16.;
			//}
			Edge edtmp3;
			edtmp3.act=1;
			edtmp3.lev=tmedge[tmesh[eid].edge[(j+1)%4]].lev+1;
			edtmp3.pt[0]=pid[0]; edtmp3.pt[1]=pid[j+1];
			edtmp3.len=tmedge[tmesh[eid].edge[(j+1)%4]].len/2.;
			tmedge.push_back(edtmp3);
			edid[3*j+2]=tmedge.size()-1;
			int ied(tmedge[tmesh[eid].edge[j]].chd[0]);
			if(tmedge[ied].pt[0]==itmp[0] || tmedge[ied].pt[1]==itmp[0])
			{
				edid[3*j]=ied;
				edid[3*j+1]=tmedge[tmesh[eid].edge[j]].chd[1];
			}
			else
			{
				edid[3*j]=tmedge[tmesh[eid].edge[j]].chd[1];
				edid[3*j+1]=ied;
			}
		}
	}

	int nvl(cp[tmesh[eid].cnct[0]].face.size());
	if(cp[tmesh[eid].cnct[0]].update==0)
	{
		double ccf[3]={1.-7./(4.*nvl),3./(2.*nvl*nvl),1./(4.*nvl*nvl)};
		cp[tmesh[eid].cnct[0]].coortmp[0]=ccf[0]*cp[tmesh[eid].cnct[0]].coor[0]+ccf[1]*cp[tmesh[eid].cnct[1]].coor[0]+ccf[2]*cp[tmesh[eid].cnct[2]].coor[0];
		cp[tmesh[eid].cnct[0]].coortmp[1]=ccf[0]*cp[tmesh[eid].cnct[0]].coor[1]+ccf[1]*cp[tmesh[eid].cnct[1]].coor[1]+ccf[2]*cp[tmesh[eid].cnct[2]].coor[1];
		cp[tmesh[eid].cnct[0]].coortmp[2]=ccf[0]*cp[tmesh[eid].cnct[0]].coor[2]+ccf[1]*cp[tmesh[eid].cnct[1]].coor[2]+ccf[2]*cp[tmesh[eid].cnct[2]].coor[2];
		cp[tmesh[eid].cnct[0]].update=2;
	}
	else if(cp[tmesh[eid].cnct[0]].update==2)
	{
		double ccf[2]={3./(2.*nvl*nvl),1./(4.*nvl*nvl)};
		cp[tmesh[eid].cnct[0]].coortmp[0]+=ccf[0]*cp[tmesh[eid].cnct[1]].coor[0]+ccf[1]*cp[tmesh[eid].cnct[2]].coor[0];
		cp[tmesh[eid].cnct[0]].coortmp[1]+=ccf[0]*cp[tmesh[eid].cnct[1]].coor[1]+ccf[1]*cp[tmesh[eid].cnct[2]].coor[1];
		cp[tmesh[eid].cnct[0]].coortmp[2]+=ccf[0]*cp[tmesh[eid].cnct[1]].coor[2]+ccf[1]*cp[tmesh[eid].cnct[2]].coor[2];
	}
	else
	{
		cout<<"Not supported for other udpate types!\n";
		getchar();
	}

	int e_cnct[4][4]={{tmesh[eid].cnct[0],pid[1],pid[0],pid[4]},{pid[1],tmesh[eid].cnct[1],pid[2],pid[0]},{pid[0],pid[2],tmesh[eid].cnct[2],pid[3]},{pid[4],pid[0],pid[3],tmesh[eid].cnct[3]}};
	int e_edge[4][4]={{edid[0],edid[2],edid[11],edid[10]},{edid[1],edid[3],edid[5],edid[2]},{edid[5],edid[4],edid[6],edid[8]},{edid[11],edid[8],edid[7],edid[9]}};
	int enewid[4];
	int e_type[4]={4,0,0,0};
	vector<Element2D> etmp(4);
	for(int i=0; i<4; i++)
	{
		etmp[i].act=1;
		etmp[i].type=e_type[i];
		etmp[i].lev=tmesh[eid].lev+1;
		etmp[i].prt=eid;
		for(int j=0; j<4; j++)
		{
			etmp[i].cnct[j]=e_cnct[i][j];
			etmp[i].edge[j]=e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i]=tmesh.size()-1;
		tmesh[eid].chd.push_back(enewid[i]);
	}

	tmesh[eid].act=0;
}

void TruncatedTspline_2D::ElementRefine_Unstruct_2(int eid, int dir)//need modification
{
	int pid[2],edid[7],pos(dir);
	int cnid[2]={pos,(pos+2)%4};
	for(int j=0; j<2; j++)
	{
		int itmp[2]={tmesh[eid].cnct[cnid[j]],tmesh[eid].cnct[(cnid[j]+1)%4]};
		if(tmedge[tmesh[eid].edge[cnid[j]]].act==1)
		{
			Vertex2D ptmp;
			cp.push_back(ptmp);
			pid[j]=cp.size()-1;
			tmedge[tmesh[eid].edge[cnid[j]]].midpt=pid[j];
			Edge edtmp1,edtmp2;
			edtmp1.act=1;
			edtmp1.lev=tmedge[tmesh[eid].edge[j]].lev+1;
			edtmp1.pt[0]=itmp[0]; edtmp1.pt[1]=pid[j];
			edtmp1.len=tmedge[tmesh[eid].edge[cnid[j]]].len/2.;
			edtmp1.prt=tmesh[eid].edge[cnid[j]];
			edtmp2.act=1;
			edtmp2.lev=tmedge[tmesh[eid].edge[j]].lev+1;
			edtmp2.pt[0]=pid[j]; edtmp2.pt[1]=itmp[1];
			edtmp2.len=tmedge[tmesh[eid].edge[cnid[j]]].len/2.;
			edtmp2.prt=tmesh[eid].edge[cnid[j]];
			tmedge.push_back(edtmp1);
			edid[2*j]=tmedge.size()-1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[0]=edid[2*j];
			tmedge.push_back(edtmp2);
			edid[2*j+1]=tmedge.size()-1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[1]=edid[2*j+1];
			tmedge[tmesh[eid].edge[cnid[j]]].act=0;
		}
		else
		{
			pid[j]=tmedge[tmesh[eid].edge[cnid[j]]].midpt;
			int ied(tmedge[tmesh[eid].edge[cnid[j]]].chd[0]);
			if(tmedge[ied].pt[0]==itmp[0] || tmedge[ied].pt[1]==itmp[0])
			{
				edid[2*j]=ied;
				edid[2*j+1]=tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
			}
			else
			{
				edid[2*j]=tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
				edid[2*j+1]=ied;
			}
		}
	}
	Edge edtmp;
	edtmp.act=1;
	edtmp.lev=tmedge[tmesh[eid].edge[(pos+1)%4]].lev;
	edtmp.pt[0]=pid[0]; edtmp.pt[1]=pid[1];
	edtmp.len=tmedge[tmesh[eid].edge[(pos+1)%4]].len;
	tmedge.push_back(edtmp);
	edid[4]=tmedge.size()-1;
	edid[5]=tmesh[eid].edge[(pos+1)%4];
	edid[6]=tmesh[eid].edge[(pos+3)%4];

	int e_cnct[2][4]={{tmesh[eid].cnct[pos],pid[0],pid[1],tmesh[eid].cnct[(pos+3)%4]},{pid[0],tmesh[eid].cnct[(pos+1)%4],tmesh[eid].cnct[(pos+2)%4],pid[1]}};
	int e_edge[2][4]={{edid[0],edid[4],edid[3],edid[6]},{edid[1],edid[5],edid[2],edid[4]}};
	int enewid[2];
	double chd_org[2][2]={{0.,0.},{tmedge[tmesh[eid].edge[0]].len/2.,0.}};
	if(dir==1)
	{
		chd_org[1][0]=0.; chd_org[1][1]=tmedge[tmesh[eid].edge[1]].len/2.;
	}
	vector<Element2D> etmp(2);
	for(int i=0; i<2; i++)
	{
		etmp[i].act=1;
		etmp[i].type=1;//need to distinguish square and rectangular
		etmp[i].lev=tmesh[eid].lev;//modify level later
		etmp[i].prt=eid;
		for(int j=0; j<4; j++)
		{
			etmp[i].cnct[(pos+j)%4]=e_cnct[i][j];
			etmp[i].edge[(pos+j)%4]=e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i]=tmesh.size()-1;
		tmesh[eid].chd.push_back(enewid[i]);
		array<double,2> chd_o_tmp={chd_org[i][0],chd_org[i][1]};
		tmesh[eid].chd_o.push_back(chd_o_tmp);
	}

	tmesh[eid].act=0;
}

void TruncatedTspline_2D::ElementRefine_Unstruct_b(int eid)
{
	int pid[2],edid[7],pos(0);
	if(tmedge[tmesh[eid].edge[0]].len == 0.) pos=1;
	int cnid[2]={pos,(pos+2)%4};
	for(int j=0; j<2; j++)
	{
		int itmp[2]={tmesh[eid].cnct[cnid[j]],tmesh[eid].cnct[(cnid[j]+1)%4]};
		if(tmedge[tmesh[eid].edge[cnid[j]]].act==1)
		{
			Vertex2D ptmp;
			//ptmp.coor[0]=(cp[itmp[0]].coor[0]+cp[itmp[1]].coor[0])/2.;
			//ptmp.coor[1]=(cp[itmp[0]].coor[1]+cp[itmp[1]].coor[1])/2.;
			//ptmp.coor[2]=(cp[itmp[0]].coor[2]+cp[itmp[1]].coor[2])/2.;
			cp.push_back(ptmp);
			pid[j]=cp.size()-1;
			tmedge[tmesh[eid].edge[cnid[j]]].midpt=pid[j];
			Edge edtmp1,edtmp2;
			edtmp1.act=1;
			edtmp1.lev=tmedge[tmesh[eid].edge[j]].lev+1;
			edtmp1.pt[0]=itmp[0]; edtmp1.pt[1]=pid[j];
			edtmp1.len=tmedge[tmesh[eid].edge[cnid[j]]].len/2.;
			edtmp1.prt=tmesh[eid].edge[cnid[j]];
			edtmp2.act=1;
			edtmp2.lev=tmedge[tmesh[eid].edge[j]].lev+1;
			edtmp2.pt[0]=pid[j]; edtmp2.pt[1]=itmp[1];
			edtmp2.len=tmedge[tmesh[eid].edge[cnid[j]]].len/2.;
			edtmp2.prt=tmesh[eid].edge[cnid[j]];
			tmedge.push_back(edtmp1);
			edid[2*j]=tmedge.size()-1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[0]=edid[2*j];
			tmedge.push_back(edtmp2);
			edid[2*j+1]=tmedge.size()-1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[1]=edid[2*j+1];
			tmedge[tmesh[eid].edge[cnid[j]]].act=0;
		}
		else
		{
			pid[j]=tmedge[tmesh[eid].edge[cnid[j]]].midpt;
			int ied(tmedge[tmesh[eid].edge[cnid[j]]].chd[0]);
			if(tmedge[ied].pt[0]==itmp[0] || tmedge[ied].pt[1]==itmp[0])
			{
				edid[2*j]=ied;
				edid[2*j+1]=tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
			}
			else
			{
				edid[2*j]=tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
				edid[2*j+1]=ied;
			}
		}
	}
	Edge edtmp;
	edtmp.act=1;
	edtmp.lev=tmedge[tmesh[eid].edge[(pos+1)%4]].lev+1;
	edtmp.pt[0]=pid[0]; edtmp.pt[1]=pid[1];
	edtmp.len=tmedge[tmesh[eid].edge[(pos+1)%4]].len;
	tmedge.push_back(edtmp);
	edid[4]=tmedge.size()-1;
	edid[5]=tmesh[eid].edge[(pos+1)%4];
	edid[6]=tmesh[eid].edge[(pos+3)%4];

	int e_cnct[2][4]={{tmesh[eid].cnct[pos],pid[0],pid[1],tmesh[eid].cnct[(pos+3)%4]},{pid[0],tmesh[eid].cnct[(pos+1)%4],tmesh[eid].cnct[(pos+2)%4],pid[1]}};
	int e_edge[2][4]={{edid[0],edid[4],edid[3],edid[6]},{edid[1],edid[5],edid[2],edid[4]}};
	int enewid[2];
	double chd_org[2][2]={{0.,0.},{tmedge[tmesh[eid].edge[0]].len/2.,0.}};
	if(pos==1)
	{
		chd_org[1][0]=0.; chd_org[1][1]=tmedge[tmesh[eid].edge[1]].len/2.;
	}
	vector<Element2D> etmp(2);
	for(int i=0; i<2; i++)
	{
		etmp[i].act=1;
		etmp[i].type=2;
		etmp[i].lev=tmesh[eid].lev+1;
		etmp[i].prt=eid;
		for(int j=0; j<4; j++)
		{
			etmp[i].cnct[(pos+j)%4]=e_cnct[i][j];
			etmp[i].edge[(pos+j)%4]=e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i]=tmesh.size()-1;
		tmesh[eid].chd.push_back(enewid[i]);
		array<double,2> chd_o_tmp={chd_org[i][0],chd_org[i][1]};
		tmesh[eid].chd_o.push_back(chd_o_tmp);
	}

	tmesh[eid].act=0;
}

void TruncatedTspline_2D::ElementRefine_Invalid_4(int eid)
{
	int pid[5];
	int edid[12];
	Vertex2D ptmp1;
	ptmp1.update=1;
	ptmp1.coortmp[0]=(cp[tmesh[eid].cnct[0]].coor[0]+cp[tmesh[eid].cnct[1]].coor[0]+cp[tmesh[eid].cnct[2]].coor[0]+cp[tmesh[eid].cnct[3]].coor[0])/4.;
	ptmp1.coortmp[1]=(cp[tmesh[eid].cnct[0]].coor[1]+cp[tmesh[eid].cnct[1]].coor[1]+cp[tmesh[eid].cnct[2]].coor[1]+cp[tmesh[eid].cnct[3]].coor[1])/4.;
	ptmp1.coortmp[2]=(cp[tmesh[eid].cnct[0]].coor[2]+cp[tmesh[eid].cnct[1]].coor[2]+cp[tmesh[eid].cnct[2]].coor[2]+cp[tmesh[eid].cnct[3]].coor[2])/4.;
	cp.push_back(ptmp1);
	pid[0]=cp.size()-1;
	for(int j=0; j<4; j++)
	{
		int itmp[2]={tmesh[eid].cnct[j],tmesh[eid].cnct[(j+1)%4]};
		if(tmedge[tmesh[eid].edge[j]].act==1)
		{
			Vertex2D ptmp;
			int id0[4]={tmesh[eid].cnct[j],tmesh[eid].cnct[(j+1)%4],tmesh[eid].cnct[(j+2)%4],tmesh[eid].cnct[(j+3)%4]};
			int fcnb(tmedge[tmesh[eid].edge[j]].face[0]);
			if(fcnb==eid) fcnb=tmedge[tmesh[eid].edge[j]].face[1];
			int* it1=find(tmesh[fcnb].cnct,tmesh[fcnb].cnct+4,tmesh[eid].cnct[j]);
			int loc1(it1-tmesh[fcnb].cnct);
			int id1[2]={tmesh[fcnb].cnct[(loc1+2)%4],tmesh[fcnb].cnct[(loc1+3)%4]};
			ptmp.coortmp[0]=(6.*cp[id0[0]].coor[0]+6.*cp[id0[1]].coor[0]+cp[id0[2]].coor[0]+cp[id0[3]].coor[0]+cp[id1[0]].coor[0]+cp[id1[1]].coor[0])/16.;
			ptmp.coortmp[1]=(6.*cp[id0[0]].coor[1]+6.*cp[id0[1]].coor[1]+cp[id0[2]].coor[1]+cp[id0[3]].coor[1]+cp[id1[0]].coor[1]+cp[id1[1]].coor[1])/16.;
			ptmp.coortmp[2]=(6.*cp[id0[0]].coor[2]+6.*cp[id0[1]].coor[2]+cp[id0[2]].coor[2]+cp[id0[3]].coor[2]+cp[id1[0]].coor[2]+cp[id1[1]].coor[2])/16.;
			cp.push_back(ptmp);
			pid[j+1]=cp.size()-1;
			tmedge[tmesh[eid].edge[j]].midpt=pid[j+1];
			Edge edtmp1,edtmp2,edtmp3;
			edtmp1.act=1;
			edtmp1.lev=tmedge[tmesh[eid].edge[j]].lev+1;
			edtmp1.pt[0]=itmp[0]; edtmp1.pt[1]=pid[j+1];
			edtmp1.len=tmedge[tmesh[eid].edge[j]].len/2.;
			edtmp1.prt=tmesh[eid].edge[j];
			edtmp2.act=1;
			edtmp2.lev=tmedge[tmesh[eid].edge[j]].lev+1;
			edtmp2.pt[0]=pid[j+1]; edtmp2.pt[1]=itmp[1];
			edtmp2.len=tmedge[tmesh[eid].edge[j]].len/2.;
			edtmp2.prt=tmesh[eid].edge[j];
			edtmp3.act=1;
			edtmp3.lev=tmedge[tmesh[eid].edge[(j+1)%4]].lev+1;
			edtmp3.pt[0]=pid[0]; edtmp3.pt[1]=pid[j+1];
			edtmp3.len=tmedge[tmesh[eid].edge[(j+1)%4]].len/2.;
			tmedge.push_back(edtmp1);
			edid[3*j]=tmedge.size()-1;
			tmedge[tmesh[eid].edge[j]].chd[0]=edid[3*j];
			tmedge.push_back(edtmp2);
			edid[3*j+1]=tmedge.size()-1;
			tmedge[tmesh[eid].edge[j]].chd[1]=edid[3*j+1];
			tmedge.push_back(edtmp3);
			edid[3*j+2]=tmedge.size()-1;
			tmedge[tmesh[eid].edge[j]].act=0;
		}
		else
		{
			pid[j+1]=tmedge[tmesh[eid].edge[j]].midpt;
			int id1[2]={tmesh[eid].cnct[(j+2)%4],tmesh[eid].cnct[(j+3)%4]};
			Edge edtmp3;
			edtmp3.act=1;
			edtmp3.lev=tmedge[tmesh[eid].edge[(j+1)%4]].lev+1;
			edtmp3.pt[0]=pid[0]; edtmp3.pt[1]=pid[j+1];
			edtmp3.len=tmedge[tmesh[eid].edge[(j+1)%4]].len/2.;
			tmedge.push_back(edtmp3);
			edid[3*j+2]=tmedge.size()-1;
			int ied(tmedge[tmesh[eid].edge[j]].chd[0]);
			if(tmedge[ied].pt[0]==itmp[0] || tmedge[ied].pt[1]==itmp[0])
			{
				edid[3*j]=ied;
				edid[3*j+1]=tmedge[tmesh[eid].edge[j]].chd[1];
			}
			else
			{
				edid[3*j]=tmedge[tmesh[eid].edge[j]].chd[1];
				edid[3*j+1]=ied;
			}
		}
	}

	for(int j=0; j<4; j++)
	{
		int nvl(cp[tmesh[eid].cnct[j]].face.size());
		if(nvl==3 || nvl>4)
		{
			if(cp[tmesh[eid].cnct[j]].update==0)
			{
				double ccf[3]={1.-7./(4.*nvl),3./(2.*nvl*nvl),1./(4.*nvl*nvl)};
				cp[tmesh[eid].cnct[j]].coortmp[0]=ccf[0]*cp[tmesh[eid].cnct[j]].coor[0]+ccf[1]*cp[tmesh[eid].cnct[(j+1)%4]].coor[0]+ccf[2]*cp[tmesh[eid].cnct[(j+2)%4]].coor[0];
				cp[tmesh[eid].cnct[j]].coortmp[1]=ccf[0]*cp[tmesh[eid].cnct[j]].coor[1]+ccf[1]*cp[tmesh[eid].cnct[(j+1)%4]].coor[1]+ccf[2]*cp[tmesh[eid].cnct[(j+2)%4]].coor[1];
				cp[tmesh[eid].cnct[j]].coortmp[2]=ccf[0]*cp[tmesh[eid].cnct[j]].coor[2]+ccf[1]*cp[tmesh[eid].cnct[(j+1)%4]].coor[2]+ccf[2]*cp[tmesh[eid].cnct[(j+2)%4]].coor[2];
				cp[tmesh[eid].cnct[j]].update=2;
			}
			else if(cp[tmesh[eid].cnct[j]].update==2)
			{
				double ccf[2]={3./(2.*nvl*nvl),1./(4.*nvl*nvl)};
				cp[tmesh[eid].cnct[j]].coortmp[0]+=ccf[0]*cp[tmesh[eid].cnct[(j+1)%4]].coor[0]+ccf[1]*cp[tmesh[eid].cnct[(j+2)%4]].coor[0];
				cp[tmesh[eid].cnct[j]].coortmp[1]+=ccf[0]*cp[tmesh[eid].cnct[(j+1)%4]].coor[1]+ccf[1]*cp[tmesh[eid].cnct[(j+2)%4]].coor[1];
				cp[tmesh[eid].cnct[j]].coortmp[2]+=ccf[0]*cp[tmesh[eid].cnct[(j+1)%4]].coor[2]+ccf[1]*cp[tmesh[eid].cnct[(j+2)%4]].coor[2];
			}
			else
			{
				cout<<"Not supported for other udpate types!\n";
				getchar();
			}
		}
	}

	int e_cnct[4][4]={{tmesh[eid].cnct[0],pid[1],pid[0],pid[4]},{pid[1],tmesh[eid].cnct[1],pid[2],pid[0]},{pid[0],pid[2],tmesh[eid].cnct[2],pid[3]},{pid[4],pid[0],pid[3],tmesh[eid].cnct[3]}};
	int e_edge[4][4]={{edid[0],edid[2],edid[11],edid[10]},{edid[1],edid[3],edid[5],edid[2]},{edid[5],edid[4],edid[6],edid[8]},{edid[11],edid[8],edid[7],edid[9]}};
	int enewid[4];
	int e_type[4]={4,0,0,0};
	for(int i=1; i<4; i++)
	{
		if(cp[tmesh[eid].cnct[i]].face.size()==3 || cp[tmesh[eid].cnct[i]].face.size()>4)
		{
			e_type[i]=4;
			int cnct_tmp[4]={e_cnct[i][0],e_cnct[i][1],e_cnct[i][2],e_cnct[i][3]};
			int edge_tmp[4]={e_edge[i][0],e_edge[i][1],e_edge[i][2],e_edge[i][3]};
			for(int j=0; j<4; j++)
			{
				e_cnct[i][j]=cnct_tmp[(i+j)%4];
				e_edge[i][j]=edge_tmp[(i+j)%4];
			}
		}
	}
	vector<Element2D> etmp(4);
	for(int i=0; i<4; i++)
	{
		etmp[i].act=1;
		etmp[i].type=e_type[i];
		etmp[i].lev=tmesh[eid].lev+1;
		etmp[i].prt=eid;
		for(int j=0; j<4; j++)
		{
			etmp[i].cnct[j]=e_cnct[i][j];
			etmp[i].edge[j]=e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i]=tmesh.size()-1;
		tmesh[eid].chd.push_back(enewid[i]);
	}

	tmesh[eid].act=0;
}

void TruncatedTspline_2D::Topo_Refine_Unstruct(const vector<int>& rfid, const vector<int>& rftype, vector<int>& rfid_more, vector<int>& rftype_more)
{
	//initialize for refinement
	npt_old=cp.size();
	nel_old=tmesh.size();
	rfid_more.clear();
	rftype_more.clear();
	rfid_more=rfid;
	rftype_more=rftype;
	for(uint i=0; i<cp.size(); i++)
	{
		cp[i].update=0;
		cp[i].aff=0;
		cp[i].coortmp[0]=0.; cp[i].coortmp[1]=0.; cp[i].coortmp[2]=0.;
		for(int j=0; j<4; j++)
		{
			cp[i].kitvUtmp[j]=0.; cp[i].kitvVtmp[j]=0.;
		}
		//cp[i].truntmp=0;
		//vector<int>().swap(cp[i].tbftmp);
		//vector<double>().swap(cp[i].tctmp);
	}
	for(uint i=0; i<tmesh.size(); i++)
	{
		vector<int>().swap(tmesh[i].IENtmp);
		vector<array<double,5>>().swap(tmesh[i].patch_kutmp);
		vector<array<double,5>>().swap(tmesh[i].patch_kvtmp);
	}
	//refine
	for(uint i=0; i<rfid.size(); i++)
	{
		if(rftype[i]==0)
		{
			ElementRefine_Unstruct_4(rfid[i]);
		}
		else if(rftype[i]==1)
		{
			ElementRefine_Unstruct_2(rfid[i],rftype[i]-1);
		}
		else if(rftype[i]==2)
		{
			ElementRefine_Unstruct_2(rfid[i],rftype[i]-1);
		}
		else if(rftype[i]==3)
		{
			ElementRefine_Unstruct_b(rfid[i]);
		}
		else if(rftype[i]==4)
		{
			ElementRefine_Irregular_4(rfid[i]);
		}
		else if(rftype[i]==5)
		{
			ElementRefine_Invalid_4(rfid[i]);
		}
		else
		{
			cout<<"Other types of elements are not supported to be refined!\n";
			getchar();
		}
	}

	vector<int> rid_b;
	for(uint i=0; i<tmesh.size(); i++)
	{
		if(tmesh[i].act==1)
		{
			if(tmesh[i].type==2)
			{
				for(int j=0; j<4; j++)
				{
					if(tmedge[tmesh[i].edge[j]].act==0)
					{
						rid_b.push_back(i);
						rfid_more.push_back(i);
						rftype_more.push_back(3);
						break;
					}
				}
			}
		}
	}
	for(uint i=0; i<rid_b.size(); i++)
	{
		ElementRefine_Unstruct_b(rid_b[i]);
	}
}

void TruncatedTspline_2D::Geom_Refine_Unstruct(const vector<int>& rfid, const vector<int>& rftype)
{
	UpdateConnect();
	FindEdgeTopoDirec();
	FindKnotInterval();
	UpdateKnotInterval();//might not work
	SetLocalCoorSystem();
	FindIEN_Invalid();

	//calculate control points
	for(uint i=0; i<rfid.size(); i++)
	{
		if(tmesh[rfid[i]].type==0 || tmesh[rfid[i]].type==1 || tmesh[rfid[i]].type==2)
		{
			UpdatePatchCP_Unstruct(rfid[i]);
		}
	}
	//update basis functions
	for(uint i=0; i<cp.size(); i++)
	{
		if(i<npt_old)
		{
			if(cp[i].aff==1)
			{
				cp[i].coor[0]=cp[i].coortmp[0];
				cp[i].coor[1]=cp[i].coortmp[1];
				cp[i].coor[2]=cp[i].coortmp[2];
			}
		}
		else
		{
			cp[i].coor[0]=cp[i].coortmp[0];
			cp[i].coor[1]=cp[i].coortmp[1];
			cp[i].coor[2]=cp[i].coortmp[2];
		}
	}

	//update new elements
	for(uint i=0; i<rfid.size(); i++)
	{
		if(tmesh[rfid[i]].type==0 || tmesh[rfid[i]].type==1)
		{
			vector<int> trunIEN;
			vector<array<double,5>> trun_patch_ku;
			vector<array<double,5>> trun_patch_kv;
			for(uint j=0; j<tmesh[rfid[i]].IEN.size(); j++)
			{
				if(cp[tmesh[rfid[i]].IEN[j]].trun==1)
				{
					trunIEN.push_back(tmesh[rfid[i]].IEN[j]);
					trun_patch_ku.push_back(tmesh[rfid[i]].patch_ku[j]);
					trun_patch_kv.push_back(tmesh[rfid[i]].patch_kv[j]);
				}
			}
			for(uint j=0; j<tmesh[rfid[i]].chd.size(); j++)
			{
				int chdid(tmesh[rfid[i]].chd[j]);
				tmesh[chdid].IEN.clear();
				tmesh[chdid].patch_ku.clear();
				tmesh[chdid].patch_kv.clear();
				for(uint k=0; k<tmesh[chdid].IENtmp.size(); k++)
				{
					tmesh[chdid].IEN.push_back(tmesh[chdid].IENtmp[k]);
					if(cp[tmesh[chdid].IENtmp[k]].trun==0)
					{
						tmesh[chdid].patch_ku.push_back(tmesh[chdid].patch_kutmp[k]);
						tmesh[chdid].patch_kv.push_back(tmesh[chdid].patch_kvtmp[k]);
					}
					else
					{
						vector<int>::iterator it=find(tmesh[rfid[i]].IEN.begin(),tmesh[rfid[i]].IEN.end(),tmesh[chdid].IENtmp[k]);
						if(it!=tmesh[rfid[i]].IEN.end())
						{
							int loc(it-tmesh[rfid[i]].IEN.begin());
							array<double,5> kutmp,kvtmp;
							for(int a=0; a<5; a++)
							{
								kutmp[a]=tmesh[rfid[i]].patch_ku[loc][a]-tmesh[rfid[i]].chd_o[j][0];
								kvtmp[a]=tmesh[rfid[i]].patch_kv[loc][a]-tmesh[rfid[i]].chd_o[j][1];
							}
							tmesh[chdid].patch_ku.push_back(kutmp);
							tmesh[chdid].patch_kv.push_back(kvtmp);
						}
						else
						{
							//cout<<i<<" "<<rfid[i]<<" "<<tmesh[rfid[i]].IEN.size()<<"\n";
							cout<<"Cannot find truncated ID in the old set!\n";
							getchar();
						}
					}
				}
				for(uint k=0; k<trunIEN.size(); k++)
				{
					vector<int>::iterator it=find(tmesh[chdid].IEN.begin(),tmesh[chdid].IEN.end(),trunIEN[k]);
					if(it==tmesh[chdid].IEN.end())
					{
						tmesh[chdid].IEN.push_back(trunIEN[k]);
						array<double,5> kutmp,kvtmp;
						for(int a=0; a<5; a++)
						{
							kutmp[a]=trun_patch_ku[k][a]-tmesh[rfid[i]].chd_o[j][0];
							kvtmp[a]=trun_patch_kv[k][a]-tmesh[rfid[i]].chd_o[j][1];
						}
						tmesh[chdid].patch_ku.push_back(kutmp);
						tmesh[chdid].patch_kv.push_back(kvtmp);
					}
				}
			}
		}
		else if(tmesh[rfid[i]].type==4 || tmesh[rfid[i]].type==5)
		{
			for(uint j=0; j<tmesh[rfid[i]].chd.size(); j++)
			{
				int chdid(tmesh[rfid[i]].chd[j]);
				tmesh[chdid].IEN.clear();
				tmesh[chdid].patch_ku.clear();
				tmesh[chdid].patch_kv.clear();
				tmesh[chdid].IEN=tmesh[chdid].IENtmp;
				tmesh[chdid].patch_ku=tmesh[chdid].patch_kutmp;
				tmesh[chdid].patch_kv=tmesh[chdid].patch_kvtmp;
				vector<int>().swap(tmesh[chdid].IENtmp);
				vector<array<double,5>>().swap(tmesh[chdid].patch_kutmp);
				vector<array<double,5>>().swap(tmesh[chdid].patch_kvtmp);
			}
		}
	}

	//update old elements
	for(uint i=0; i<nel_old; i++)
	{
		if(tmesh[i].act==1 && (tmesh[i].type==0 || tmesh[i].type==1))
		{
			vector<int> ien_old(tmesh[i].IEN);
			vector<array<double,5>> ku_old(tmesh[i].patch_ku);
			vector<array<double,5>> kv_old(tmesh[i].patch_kv);
			tmesh[i].IEN.clear();
			tmesh[i].patch_ku.clear();
			tmesh[i].patch_kv.clear();
			for(uint j=0; j<tmesh[i].IENtmp.size(); j++)
			{
				tmesh[i].IEN.push_back(tmesh[i].IENtmp[j]);
				if(cp[tmesh[i].IENtmp[j]].trun==0)
				{
					tmesh[i].patch_ku.push_back(tmesh[i].patch_kutmp[j]);
					tmesh[i].patch_kv.push_back(tmesh[i].patch_kvtmp[j]);
				}
				else
				{
					vector<int>::iterator it=find(ien_old.begin(),ien_old.end(),tmesh[i].IENtmp[j]);
					int loc(it-ien_old.begin());
					tmesh[i].patch_ku.push_back(ku_old[loc]);
					tmesh[i].patch_kv.push_back(kv_old[loc]);
				}
			}
			vector<int>().swap(tmesh[i].IENtmp);
			vector<array<double,5>>().swap(tmesh[i].patch_kutmp);
			vector<array<double,5>>().swap(tmesh[i].patch_kvtmp);
		}
		else if(tmesh[i].act==1 && tmesh[i].type==4)
		{
			vector<int> ien_old(tmesh[i].IEN);
			vector<array<double,5>> ku_old(tmesh[i].patch_ku);
			vector<array<double,5>> kv_old(tmesh[i].patch_kv);
			tmesh[i].IEN.clear();
			tmesh[i].patch_ku.clear();
			tmesh[i].patch_kv.clear();
			uint n_1r(2*cp[tmesh[i].cnct[0]].face.size()+1);
			for(uint j=0; j<tmesh[i].IENtmp.size(); j++)
			{
				tmesh[i].IEN.push_back(tmesh[i].IENtmp[j]);
				if(j>=n_1r)
				{
					if(cp[tmesh[i].IENtmp[j]].trun==0)
					{
						tmesh[i].patch_ku.push_back(tmesh[i].patch_kutmp[j-n_1r]);
						tmesh[i].patch_kv.push_back(tmesh[i].patch_kvtmp[j-n_1r]);
					}
					else
					{
						vector<int>::iterator it=find(ien_old.begin(),ien_old.end(),tmesh[i].IENtmp[j]);
						int loc(it-ien_old.begin());
						tmesh[i].patch_ku.push_back(ku_old[loc-n_1r]);
						tmesh[i].patch_kv.push_back(kv_old[loc-n_1r]);
					}
				}
			}
			vector<int>().swap(tmesh[i].IENtmp);
			vector<array<double,5>>().swap(tmesh[i].patch_kutmp);
			vector<array<double,5>>().swap(tmesh[i].patch_kvtmp);
		}
	}
}

void TruncatedTspline_2D::BezierExtract_Unstruct(vector<BezierElement2D>& bzmesh)
{
	bzmesh.clear();
	for(uint eid=0;eid<tmesh.size();eid++)
	{
		if(tmesh[eid].act==1)
		{
			if(tmesh[eid].type==0 || tmesh[eid].type==1)
			{
				BezierElementExtract_Unstruct(eid,bzmesh);
			}
			else if(tmesh[eid].type==4)
			{
				BezierElementExtract_Unstruct_Irr(eid,bzmesh);
			}
		}
	}
}

void TruncatedTspline_2D::BezierElementExtract_Unstruct(int eid,vector<BezierElement2D>& bzmesh)
{
	vector<array<double,4>> be;
	BezierFinder_Unstruct(eid,be);
	uint nbzold=bzmesh.size();
	for(uint i=0;i<be.size();i++)
	{
		//BezierUnit_Unstruct(eid,be[i],bzmesh);
		BezierUnit_Unstruct_Trun(eid,be[i],bzmesh);
	}
	for(uint i=nbzold; i<bzmesh.size(); i++)
	{
		bzmesh[i].prt=eid;
	}
}

void TruncatedTspline_2D::BezierFinder_Unstruct(int eid,vector<array<double,4>>& be)
{
	be.clear();
	array<double,4> e0={0.,tmedge[tmesh[eid].edge[0]].len,0.,tmedge[tmesh[eid].edge[3]].len};
	vector<double> ukt,vkt;
	for(uint i=0; i<tmesh[eid].patch_ku.size(); i++)
	{
		for(int j=0; j<5; j++)
		{
			if(tmesh[eid].patch_ku[i][j]>e0[0] && tmesh[eid].patch_ku[i][j]<e0[1])
			{
				vector<double>::iterator it=find(ukt.begin(),ukt.end(),tmesh[eid].patch_ku[i][j]);
				if(it==ukt.end()) ukt.push_back(tmesh[eid].patch_ku[i][j]);
			}
			if(tmesh[eid].patch_kv[i][j]>e0[2] && tmesh[eid].patch_kv[i][j]<e0[3])
			{
				vector<double>::iterator it=find(vkt.begin(),vkt.end(),tmesh[eid].patch_kv[i][j]);
				if(it==vkt.end()) vkt.push_back(tmesh[eid].patch_kv[i][j]);
			}
		}
	}
	if(ukt.size()!=0 && vkt.size()==0)
	{
		sort(ukt.begin(),ukt.end());
		vector<double> ukt1;
		ukt1.push_back(e0[0]);
		for(uint i=0; i<ukt.size(); i++)
		{
			ukt1.push_back(ukt[i]);
		}
		ukt1.push_back(e0[1]);
		for(uint i=0; i<ukt1.size()-1; i++)
		{
			array<double,4> betmp={ukt1[i],ukt1[i+1],e0[2],e0[3]};
			be.push_back(betmp);
		}
	}
	else if(ukt.size()==0 && vkt.size()!=0)
	{
		sort(vkt.begin(),vkt.end());
		vector<double> vkt1;
		vkt1.push_back(e0[2]);
		for(uint i=0; i<vkt.size(); i++)
		{
			vkt1.push_back(vkt[i]);
		}
		vkt1.push_back(e0[3]);
		for(uint i=0; i<vkt1.size()-1; i++)
		{
			array<double,4> betmp={e0[0],e0[1],vkt1[i],vkt1[i+1]};
			be.push_back(betmp);
		}
	}
	else if(ukt.size()==0 && vkt.size()==0)
	{
		be.push_back(e0);
	}
	else
	{
		cerr<<"Not available for face intersection!\n";
		getchar();
	}
}

void TruncatedTspline_2D::BezierUnit_Unstruct(int eid,array<double,4>& kts,vector<BezierElement2D>& bzmesh)
{
	vector<int> pid=tmesh[eid].IEN;
	BezierElement2D bzel;
	bzel.cmat.resize(pid.size());
	for(uint i=0;i<pid.size();i++)
	{
		for(int j=0;j<16;j++)
		{
			bzel.cmat[i][j]=0.;
		}
	}
	vector<vector<double>> coef(pid.size(),vector<double>(16));
	array<double,2> ktsU={kts[0],kts[1]};
	array<double,2> ktsV={kts[2],kts[3]};
	double bzku[6]={ktsU[0],ktsU[0],ktsU[0],ktsU[1],ktsU[1],ktsU[1]};
	double bzkv[6]={ktsV[0],ktsV[0],ktsV[0],ktsV[1],ktsV[1],ktsV[1]};
	for(uint i=0;i<pid.size();i++)
	{
		if(tmesh[eid].patch_ku[i][0]<ktsU[1] && tmesh[eid].patch_ku[i][4]>ktsU[0] && tmesh[eid].patch_kv[i][0]<ktsV[1] && tmesh[eid].patch_kv[i][4]>ktsV[0])
		{
			vector<double> ku(tmesh[eid].patch_ku[i].begin(),tmesh[eid].patch_ku[i].end());
			vector<double> kv(tmesh[eid].patch_kv[i].begin(),tmesh[eid].patch_kv[i].end());
			vector<double> ku1,kv1;
			vector<vector<double>> Tu,Tv;
			BezierInsertKnots(ku,ktsU,ku1);
			BezierInsertKnots(kv,ktsV,kv1);
			TMatrix(ku,ku1,3,Tu);
			TMatrix(kv,kv1,3,Tv);
			vector<double>::iterator it1=search(ku1.begin(),ku1.end(),bzku,bzku+6);
			vector<double>::iterator it2=search(kv1.begin(),kv1.end(),bzkv,bzkv+6);
			int loc1(it1-ku1.begin()-1),loc2(it2-kv1.begin()-1),count(0);
			for(int j=0;j<4;j++)
			{
				for(int k=0;k<4;k++)
				{
					coef[i][count]=Tu[loc1+k][0]*Tv[loc2+j][0];
					count++;
				}
			}
		}
	}

	for(uint i=0;i<pid.size();i++)
	{
		if(tmesh[eid].patch_ku[i][0]<ktsU[1] && tmesh[eid].patch_ku[i][4]>ktsU[0] && tmesh[eid].patch_kv[i][0]<ktsV[1] && tmesh[eid].patch_kv[i][4]>ktsV[0])
		{
			for(int j=0;j<16;j++)
			{
				double tmp(coef[i][j]);
				bzel.pts[j][0]+=tmp*cp[pid[i]].coor[0];
				bzel.pts[j][1]+=tmp*cp[pid[i]].coor[1];
				bzel.pts[j][2]+=tmp*cp[pid[i]].coor[2];
				bzel.cmat[i][j]=tmp;
			}
		}
	}

	for(uint i=0;i<pid.size();i++)
	{
		bzel.IEN.push_back(pid[i]);
	}
	bzmesh.push_back(bzel);
}

void TruncatedTspline_2D::BezierUnit_Unstruct_Irr(int eid,array<double,4>& kts,vector<BezierElement2D>& bzmesh)
{
	vector<int> pid=tmesh[eid].IEN;
	BezierElement2D bzel(4);
	bzel.cmat.resize(pid.size(),vector<double>(25,0.));
	//for(uint i=0;i<pid.size();i++)
	//{
	//	for(int j=0;j<25;j++)
	//	{
	//		bzel.cmat[i][j]=0.;
	//	}
	//}

	//first 2N+1
	uint nv(cp[tmesh[eid].cnct[0]].face.size());
	vector<vector<double>> bmat;
	SetBezier4TranMatOP(nv,bmat);
	//SetBezier4TranMat(nv,bmat);
	//SetBezier3TranMat(nv,bmat3);
	int dir(-1);
	double insU(kts[0]),insV(kts[2]);
	//double urang[2]={0.,1.}, vrang[2]={0.,1.};
	double urang[2]={0.,tmedge[tmesh[eid].edge[0]].len}, vrang[2]={0.,tmedge[tmesh[eid].edge[3]].len};
	int pos(0);
	for(int i=0; i<2; i++)
	{
		if(kts[i]>urang[0] && kts[i]<urang[1])
		{
			if(i==0) pos=4;
			dir=0; insU=kts[i]; break;
		}
		if(kts[i+2]>vrang[0] && kts[i+2]<vrang[1])
		{
			if(i==0) pos=4;
			dir=1; insV=kts[i+2]; break;
		}
	}
	vector<double> knv0, knv1;
	if(dir==0)
	{
		double bzku_0[10]={urang[0],urang[0],urang[0],urang[0],urang[0],urang[1],urang[1],urang[1],urang[1],urang[1]};
		double bzku_1[14]={urang[0],urang[0],urang[0],urang[0],urang[0],insU,insU,insU,insU,urang[1],urang[1],urang[1],urang[1],urang[1]};
		knv0.assign(bzku_0,bzku_0+10);
		knv1.assign(bzku_1,bzku_1+14);
	}
	else if(dir==1)
	{
		double bzkv_0[10]={vrang[0],vrang[0],vrang[0],vrang[0],vrang[0],vrang[1],vrang[1],vrang[1],vrang[1],vrang[1]};
		double bzkv_1[14]={vrang[0],vrang[0],vrang[0],vrang[0],vrang[0],insV,insV,insV,insV,vrang[1],vrang[1],vrang[1],vrang[1],vrang[1]};
		knv0.assign(bzkv_0,bzkv_0+10);
		knv1.assign(bzkv_1,bzkv_1+14);
	}
	if(dir==0 || dir==1)
	{
		vector<vector<double>> tmat0,tmat_u,tmat_v;
		TMatrix(knv0,knv1,4,tmat0);
		tmat_u.resize(5,vector<double>(5,0.));
		tmat_v.resize(5,vector<double>(5,0.));
		if(dir==0)
		{
			for(uint i=0; i<5; i++)
			{
				tmat_v[i][i]=1.;
				for(uint j=0; j<5; j++)
				{
					tmat_u[i][j]=tmat0[pos+i][j];
				}
			}
		}
		else
		{
			for(uint i=0; i<5; i++)
			{
				tmat_u[i][i]=1.;
				for(uint j=0; j<5; j++)
				{
					tmat_v[i][j]=tmat0[pos+i][j];
				}
			}
		}
		//vector<vector<double>> Tmat(25,vector<double>(25));
		//int loc_old(0), loc_new(0);
		//for(int i=0; i<5; i++)//v-direction, old
		//{
		//	for(int j=0; j<5; j++)//u-direction, old
		//	{
		//		loc_new=0;
		//		for(int k=0; k<5; k++)//v-direction, new
		//		{
		//			for(int l=0; l<5; l++)//u-direction, new
		//			{
		//				Tmat[loc_old][loc_new]=tmat_u[l][j]*tmat_v[k][i];
		//				loc_new++;
		//			}
		//		}
		//		loc_old++;
		//	}
		//}
		for(uint i=0; i<2*nv+1; i++)
		{
			for(int j=0;j<25;j++)
			{
				for(int k=0; k<25; k++)
				{
					int iold(j/5), jold(j%5), knew(k/5), lnew(k%5);
					bzel.cmat[i][k]+=bmat[i][j]*tmat_u[lnew][jold]*tmat_v[knew][iold];
					//bzel.cmat4[i][k]+=bmat[i][j]*Tmat[j][k];
				}
			}
		}
	}
	else
	{
		for(uint i=0; i<2*nv+1; i++)
		{
			for(int j=0;j<25;j++)
			{
				bzel.cmat[i][j]=bmat[i][j];
			}
		}
		//for(uint i=0; i<2*nv+1; i++)
		//{
		//	for(int j=0;j<16;j++)
		//	{
		//		bzel.cmat[i][j]=bmat3[i][j];
		//	}
		//}
	}

	//the others
	vector<vector<double>> coef(tmesh[eid].patch_ku.size(),vector<double>(16,0.));
	array<double,2> ktsU={kts[0],kts[1]};
	array<double,2> ktsV={kts[2],kts[3]};
	double bzku[6]={ktsU[0],ktsU[0],ktsU[0],ktsU[1],ktsU[1],ktsU[1]};
	double bzkv[6]={ktsV[0],ktsV[0],ktsV[0],ktsV[1],ktsV[1],ktsV[1]};
	for(uint i=0;i<tmesh[eid].patch_ku.size();i++)
	{
		if(tmesh[eid].patch_ku[i][0]<ktsU[1] && tmesh[eid].patch_ku[i][4]>ktsU[0] && tmesh[eid].patch_kv[i][0]<ktsV[1] && tmesh[eid].patch_kv[i][4]>ktsV[0])
		{
			vector<double> ku(tmesh[eid].patch_ku[i].begin(),tmesh[eid].patch_ku[i].end());
			vector<double> kv(tmesh[eid].patch_kv[i].begin(),tmesh[eid].patch_kv[i].end());
			vector<double> ku1,kv1;
			vector<vector<double>> Tu,Tv;
			BezierInsertKnots(ku,ktsU,ku1);
			BezierInsertKnots(kv,ktsV,kv1);
			TMatrix(ku,ku1,3,Tu);
			TMatrix(kv,kv1,3,Tv);
			vector<double>::iterator it1=search(ku1.begin(),ku1.end(),bzku,bzku+6);
			vector<double>::iterator it2=search(kv1.begin(),kv1.end(),bzkv,bzkv+6);
			int loc1(it1-ku1.begin()-1),loc2(it2-kv1.begin()-1),count(0);
			for(int j=0;j<4;j++)
			{
				for(int k=0;k<4;k++)
				{
					coef[i][count]=Tu[loc1+k][0]*Tv[loc2+j][0];
					count++;
				}
			}
		}
	}
	//for(uint i=0; i<tmesh[eid].patch_ku.size(); i++)
	//{
	//	for(int j=0; j<16; j++)
	//	{
	//		bzel.cmat[i+(2*nv+1)][j]=coef[i][j];
	//	}
	//}
	vector<vector<double>> demat;
	DegreeElevate(demat);
	for(uint i=0; i<tmesh[eid].patch_ku.size(); i++)
	{
		for(int j=0; j<16; j++)
		{
			for(int k=0; k<25; k++)
			{
				bzel.cmat[i+(2*nv+1)][k]+=coef[i][j]*demat[k][j];
			}
		}
	}
	for(uint i=0;i<pid.size();i++)
	{
		for(int j=0;j<25;j++)
		{
			bzel.pts[j][0]+=bzel.cmat[i][j]*cp[pid[i]].coor[0];
			bzel.pts[j][1]+=bzel.cmat[i][j]*cp[pid[i]].coor[1];
			bzel.pts[j][2]+=bzel.cmat[i][j]*cp[pid[i]].coor[2];
		}
	}
	//for(uint i=0;i<pid.size();i++)
	//{
	//	for(int j=0;j<16;j++)
	//	{
	//		bzel.pts[j][0]+=bzel.cmat[i][j]*cp[pid[i]].coor[0];
	//		bzel.pts[j][1]+=bzel.cmat[i][j]*cp[pid[i]].coor[1];
	//		bzel.pts[j][2]+=bzel.cmat[i][j]*cp[pid[i]].coor[2];
	//	}
	//}
	for(uint i=0;i<pid.size();i++)
	{
		bzel.IEN.push_back(pid[i]);
	}
	bzmesh.push_back(bzel);
}

void TruncatedTspline_2D::BezierUnit_Unstruct_Trun(int eid,array<double,4>& kts,vector<BezierElement2D>& bzmesh)
{
	vector<int> pid=tmesh[eid].IEN;
	BezierElement2D bzel;
	bzel.cmat.resize(pid.size(),vector<double>(16,0.));
	for(uint i=0;i<pid.size();i++)
	{
		for(int j=0;j<16;j++)
		{
			bzel.cmat[i][j]=0.;
		}
	}
	vector<vector<double>> coef(pid.size(),vector<double>(16));
	array<double,2> ktsU={kts[0],kts[1]};
	array<double,2> ktsV={kts[2],kts[3]};
	double bzku[6]={ktsU[0],ktsU[0],ktsU[0],ktsU[1],ktsU[1],ktsU[1]};
	double bzkv[6]={ktsV[0],ktsV[0],ktsV[0],ktsV[1],ktsV[1],ktsV[1]};
	for(uint i=0;i<pid.size();i++)
	{
		if(tmesh[eid].patch_ku[i][0]<ktsU[1] && tmesh[eid].patch_ku[i][4]>ktsU[0] && tmesh[eid].patch_kv[i][0]<ktsV[1] && tmesh[eid].patch_kv[i][4]>ktsV[0])
		{
			vector<double> ku(tmesh[eid].patch_ku[i].begin(),tmesh[eid].patch_ku[i].end());
			vector<double> kv(tmesh[eid].patch_kv[i].begin(),tmesh[eid].patch_kv[i].end());
			vector<double> ku1,kv1;
			vector<vector<double>> Tu,Tv;
			BezierInsertKnots(ku,ktsU,ku1);
			BezierInsertKnots(kv,ktsV,kv1);
			TMatrix(ku,ku1,3,Tu);
			TMatrix(kv,kv1,3,Tv);
			vector<double>::iterator it1=search(ku1.begin(),ku1.end(),bzku,bzku+6);
			vector<double>::iterator it2=search(kv1.begin(),kv1.end(),bzkv,bzkv+6);
			int loc1(it1-ku1.begin()-1),loc2(it2-kv1.begin()-1),count(0);
			for(int j=0;j<4;j++)
			{
				for(int k=0;k<4;k++)
				{
					coef[i][count]=Tu[loc1+k][0]*Tv[loc2+j][0];
					count++;
				}
			}
		}
	}
	for(uint i=0;i<pid.size();i++)
	{
		if(tmesh[eid].patch_ku[i][0]<ktsU[1] && tmesh[eid].patch_ku[i][4]>ktsU[0] && tmesh[eid].patch_kv[i][0]<ktsV[1] && tmesh[eid].patch_kv[i][4]>ktsV[0])
		{
			for(int j=0;j<16;j++)
			{
				double tmp(coef[i][j]);
				for(uint k=0;k<cp[pid[i]].tbf.size();k++)
				{
					vector<int>::iterator it=find(pid.begin(),pid.end(),cp[pid[i]].tbf[k]);
					if(it!=pid.end())
					{
						int loc=it-pid.begin();
						tmp-=cp[pid[i]].tc[k]*coef[loc][j];
					}
				}
				bzel.pts[j][0]+=tmp*cp[pid[i]].coor[0];
				bzel.pts[j][1]+=tmp*cp[pid[i]].coor[1];
				bzel.pts[j][2]+=tmp*cp[pid[i]].coor[2];
				bzel.cmat[i][j]=tmp;
			}
		}
	}
	for(uint i=0;i<pid.size();i++)
	{
		bzel.IEN.push_back(pid[i]);
	}
	bzmesh.push_back(bzel);
}

void TruncatedTspline_2D::BezierUnit_Unstruct_Irr_Trun(int eid,array<double,4>& kts,vector<BezierElement2D>& bzmesh)
{
	vector<int> pid=tmesh[eid].IEN;
	BezierElement2D bzel(4);
	bzel.cmat.resize(pid.size(),vector<double>(25,0.));
	for(uint i=0;i<pid.size();i++)
	{
		for(int j=0;j<25;j++)
		{
			bzel.cmat[i][j]=0.;
		}
	}
	//for(uint i=0;i<pid.size();i++)
	//{
	//	for(int j=0;j<16;j++)
	//	{
	//		bzel.cmat[i][j]=0.;
	//	}
	//}

	//first 2N+1
	uint nv(cp[tmesh[eid].cnct[0]].face.size());
	//vector<vector<double>> bmat,bmat3;
	//SetBezier4TranMatOP(nv,bmat);
	//SetBezier4TranMat(nv,bmat);
	//SetBezier3TranMat(nv,bmat3);
	vector<vector<double>> cmat(pid.size(),vector<double>(25));
	int dir(-1);
	double insU(kts[0]),insV(kts[2]);
	//double urang[2]={0.,1.}, vrang[2]={0.,1.};
	double urang[2]={0.,tmedge[tmesh[eid].edge[0]].len}, vrang[2]={0.,tmedge[tmesh[eid].edge[3]].len};
	int pos(0);
	for(int i=0; i<2; i++)
	{
		if(kts[i]>urang[0] && kts[i]<urang[1])
		{
			if(i==0) pos=4;
			dir=0; insU=kts[i]; break;
		}
		if(kts[i+2]>vrang[0] && kts[i+2]<vrang[1])
		{
			if(i==0) pos=4;
			dir=1; insV=kts[i+2]; break;
		}
	}
	vector<double> knv0, knv1;
	if(dir==0)
	{
		double bzku_0[10]={urang[0],urang[0],urang[0],urang[0],urang[0],urang[1],urang[1],urang[1],urang[1],urang[1]};
		double bzku_1[14]={urang[0],urang[0],urang[0],urang[0],urang[0],insU,insU,insU,insU,urang[1],urang[1],urang[1],urang[1],urang[1]};
		knv0.assign(bzku_0,bzku_0+10);
		knv1.assign(bzku_1,bzku_1+14);
	}
	else if(dir==1)
	{
		double bzkv_0[10]={vrang[0],vrang[0],vrang[0],vrang[0],vrang[0],vrang[1],vrang[1],vrang[1],vrang[1],vrang[1]};
		double bzkv_1[14]={vrang[0],vrang[0],vrang[0],vrang[0],vrang[0],insV,insV,insV,insV,vrang[1],vrang[1],vrang[1],vrang[1],vrang[1]};
		knv0.assign(bzkv_0,bzkv_0+10);
		knv1.assign(bzkv_1,bzkv_1+14);
	}
	if(dir==0 || dir==1)
	{
		vector<vector<double>> tmat0,tmat_u,tmat_v;
		TMatrix(knv0,knv1,4,tmat0);
		tmat_u.resize(5,vector<double>(5,0.));
		tmat_v.resize(5,vector<double>(5,0.));
		if(dir==0)
		{
			for(uint i=0; i<5; i++)
			{
				tmat_v[i][i]=1.;
				for(uint j=0; j<5; j++)
				{
					tmat_u[i][j]=tmat0[pos+i][j];
				}
			}
		}
		else
		{
			for(uint i=0; i<5; i++)
			{
				tmat_u[i][i]=1.;
				for(uint j=0; j<5; j++)
				{
					tmat_v[i][j]=tmat0[pos+i][j];
				}
			}
		}
		for(uint i=0; i<2*nv+1; i++)
		{
			for(int j=0;j<25;j++)
			{
				for(int k=0; k<25; k++)
				{
					int iold(j/5), jold(j%5), knew(k/5), lnew(k%5);
					cmat[i][k]+=tmesh[eid].bemat[i][j]*tmat_u[lnew][jold]*tmat_v[knew][iold];
				}
			}
		}
	}
	else
	{
		for(uint i=0; i<2*nv+1; i++)
		{
			for(int j=0;j<25;j++)
			{
				cmat[i][j]=tmesh[eid].bemat[i][j];
			}
		}
		//for(uint i=0; i<2*nv+1; i++)
		//{
		//	for(int j=0;j<16;j++)
		//	{
		//		bzel.cmat[i][j]=bmat3[i][j];
		//	}
		//}
	}

	//the others
	vector<vector<double>> coef(tmesh[eid].patch_ku.size(),vector<double>(16,0.));
	array<double,2> ktsU={kts[0],kts[1]};
	array<double,2> ktsV={kts[2],kts[3]};
	double bzku[6]={ktsU[0],ktsU[0],ktsU[0],ktsU[1],ktsU[1],ktsU[1]};
	double bzkv[6]={ktsV[0],ktsV[0],ktsV[0],ktsV[1],ktsV[1],ktsV[1]};
	for(uint i=0;i<tmesh[eid].patch_ku.size();i++)
	{
		if(tmesh[eid].patch_ku[i][0]<ktsU[1] && tmesh[eid].patch_ku[i][4]>ktsU[0] && tmesh[eid].patch_kv[i][0]<ktsV[1] && tmesh[eid].patch_kv[i][4]>ktsV[0])
		{
			vector<double> ku(tmesh[eid].patch_ku[i].begin(),tmesh[eid].patch_ku[i].end());
			vector<double> kv(tmesh[eid].patch_kv[i].begin(),tmesh[eid].patch_kv[i].end());
			vector<double> ku1,kv1;
			vector<vector<double>> Tu,Tv;
			BezierInsertKnots(ku,ktsU,ku1);
			BezierInsertKnots(kv,ktsV,kv1);
			TMatrix(ku,ku1,3,Tu);
			TMatrix(kv,kv1,3,Tv);
			vector<double>::iterator it1=search(ku1.begin(),ku1.end(),bzku,bzku+6);
			vector<double>::iterator it2=search(kv1.begin(),kv1.end(),bzkv,bzkv+6);
			int loc1(it1-ku1.begin()-1),loc2(it2-kv1.begin()-1),count(0);
			for(int j=0;j<4;j++)
			{
				for(int k=0;k<4;k++)
				{
					coef[i][count]=Tu[loc1+k][0]*Tv[loc2+j][0];
					count++;
				}
			}
		}
	}
	//for(uint i=0; i<tmesh[eid].patch_ku.size(); i++)
	//{
	//	for(int j=0; j<16; j++)
	//	{
	//		bzel.cmat[i+(2*nv+1)][j]=coef[i][j];
	//	}
	//}
	vector<vector<double>> demat;
	DegreeElevate(demat);
	for(uint i=0; i<tmesh[eid].patch_ku.size(); i++)
	{
		for(int j=0; j<16; j++)
		{
			for(int k=0; k<25; k++)
			{
				cmat[i+(2*nv+1)][k]+=coef[i][j]*demat[k][j];
			}
		}
	}
	//truncation
	for(uint i=0;i<pid.size();i++)
	{
		for(int j=0;j<25;j++)
		{
			bzel.cmat[i][j]=cmat[i][j];
			if(cp[pid[i]].trun==1)
			{
				for(uint k=0;k<cp[pid[i]].tbf.size();k++)
				{
					vector<int>::iterator it=find(pid.begin(),pid.end(),cp[pid[i]].tbf[k]);
					if(it!=pid.end())
					{
						int loc=it-pid.begin();
						bzel.cmat[i][j]-=cp[pid[i]].tc[k]*cmat[loc][j];
					}
				}
			}
		}
	}
	for(uint i=0;i<pid.size();i++)
	{
		for(int j=0;j<25;j++)
		{
			bzel.pts[j][0]+=bzel.cmat[i][j]*cp[pid[i]].coor[0];
			bzel.pts[j][1]+=bzel.cmat[i][j]*cp[pid[i]].coor[1];
			bzel.pts[j][2]+=bzel.cmat[i][j]*cp[pid[i]].coor[2];
		}
	}
	//for(uint i=0;i<pid.size();i++)
	//{
	//	for(int j=0;j<16;j++)
	//	{
	//		bzel.pts[j][0]+=bzel.cmat[i][j]*cp[pid[i]].coor[0];
	//		bzel.pts[j][1]+=bzel.cmat[i][j]*cp[pid[i]].coor[1];
	//		bzel.pts[j][2]+=bzel.cmat[i][j]*cp[pid[i]].coor[2];
	//	}
	//}
	for(uint i=0;i<pid.size();i++)
	{
		bzel.IEN.push_back(pid[i]);
	}
	bzmesh.push_back(bzel);
}

void TruncatedTspline_2D::BezierElementExtract_Unstruct_Irr(int eid,vector<BezierElement2D>& bzmesh)
{
	vector<array<double,4>> be;
	BezierFinder_Unstruct(eid,be);
	uint nbzold=bzmesh.size();
	for(uint i=0;i<be.size();i++)
	{
		//BezierUnit_Unstruct_Irr(eid,be[i],bzmesh);
		BezierUnit_Unstruct_Irr_Trun(eid,be[i],bzmesh);
	}
	for(uint i=nbzold; i<bzmesh.size(); i++)
	{
		bzmesh[i].prt=eid;
	}
}

void TruncatedTspline_2D::BezierVTK_Unstruct(string fn,vector<BezierElement2D>& bzmesh)
{
	vector<array<double,3>> spt;
	vector<double> sval;
	vector<array<double,3>> norm;
	vector<array<int,4>> sele;
	vector<array<double,3>> lpt;//visulize parameter lines
	vector<array<int,2>> led;//line connectivity
	int ns(5),ecount(0);
	vector<double> su(ns),sv(ns);
	for(int i=0;i<ns;i++)
	{
		su[i]=double(i)/double(ns-1);
		sv[i]=double(i)/double(ns-1);
	}

	for(uint e=0;e<bzmesh.size();e++)
	{
		int loc(0);
		for(int a=0;a<ns;a++)
		{
			for(int b=0;b<ns;b++)
			{
				array<double,3> pt,nm;
				bzmesh[e].SurfPointNormal(su[b],sv[a],pt,nm);
				spt.push_back(pt);
				norm.push_back(nm);
				if(a==0||a==ns-1||b==0||b==ns-1)
				{
					lpt.push_back(pt);
				}
			}
		}
		for(int a=0;a<ns-1;a++)
		{
			for(int b=0;b<ns-1;b++)
			{
				array<int,4> el;
				el[0]=ecount*ns*ns+a*ns+b;
				el[1]=ecount*ns*ns+a*ns+b+1;
				el[2]=ecount*ns*ns+(a+1)*ns+b+1;
				el[3]=ecount*ns*ns+(a+1)*ns+b;
				sele.push_back(el);
			}
		}
		for(int a=0;a<ns-1;a++)
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
		led.push_back(lc1);
		ecount++;
	}

	string fname=fn+".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if(fout.is_open())
	{
		fout<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout<<"POINTS "<<spt.size()<<" float\n";
		for(uint i=0;i<spt.size();i++)
		{
			fout<<spt[i][0]<<" "<<spt[i][1]<<" "<<spt[i][2]<<"\n";
		}
		fout<<"\nCELLS "<<sele.size()<<" "<<5*sele.size()<<'\n';
		for(uint i=0;i<sele.size();i++)
		{
			fout<<"4 "<<sele[i][0]<<" "<<sele[i][1]<<" "<<sele[i][2]<<" "<<sele[i][3]<<'\n';
		}
		fout<<"\nCELL_TYPES "<<sele.size()<<'\n';
		for(uint i=0;i<sele.size();i++)
		{
			fout<<"9\n";
		}
		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i]<<"\n";
		//}
		fout<<"\nPOINT_DATA "<<norm.size()<<"\nNORMALS Normal FLOAT\n";
		for(uint i=0;i<norm.size();i++)
		{
			fout<<norm[i][0]<<" "<<norm[i][1]<<" "<<norm[i][2]<<"\n";
		}
		fout.close();
	}
	else
	{
		cout<<"Cannot open "<<fname<<"!\n";
	}

	string fname1(fn+"-lines.vtk");
	ofstream fout1;
	fout1.open(fname1.c_str());
	if(fout1.is_open())
	{
		fout1<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout1<<"POINTS "<<lpt.size()<<" float\n";
		for(uint i=0;i<lpt.size();i++)
		{
			fout1<<lpt[i][0]<<" "<<lpt[i][1]<<" "<<lpt[i][2]<<"\n";
		}
		fout1<<"\nCELLS "<<led.size()<<" "<<3*led.size()<<'\n';
		for(uint i=0;i<led.size();i++)
		{
			fout1<<"2 "<<led[i][0]<<" "<<led[i][1]<<'\n';
		}
		fout1<<"\nCELL_TYPES "<<led.size()<<'\n';
		for(uint i=0;i<led.size();i++)
		{
			fout1<<"3\n";
		}
		fout1.close();
	}
	else
	{
		cout<<"Cannot open "<<fname1<<"!\n";
	}
}

void TruncatedTspline_2D::BezierControlMesh_Unstruct(string fn,vector<BezierElement2D>& bzmesh)
{
	vector<array<double,3>> spt;
	vector<double> sval;
	vector<array<int,4>> sele;
	//int ns(2);
	int ecount(0);
	//vector<double> su(ns),sv(ns);
	//for(int i=0;i<ns;i++)
	//{
	//	su[i]=double(i)/double(ns-1);
	//	sv[i]=double(i)/double(ns-1);
	//}

	for(uint e=0;e<bzmesh.size();e++)
	{
		int loc(0);
		int ns=bzmesh[e].order+1;
		int nold=spt.size();
		for(int a=0;a<ns;a++)
		{
			for(int b=0;b<ns;b++)
			{
				spt.push_back(bzmesh[e].pts[loc]);
				loc++;
			}
		}
		for(int a=0;a<ns-1;a++)
		{
			for(int b=0;b<ns-1;b++)
			{
				array<int,4> el;
				el[0]=nold+a*ns+b;
				el[1]=nold+a*ns+b+1;
				el[2]=nold+(a+1)*ns+b+1;
				el[3]=nold+(a+1)*ns+b;
				sele.push_back(el);
			}
		}
	}

	string fname=fn+".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if(fout.is_open())
	{
		fout<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout<<"POINTS "<<spt.size()<<" float\n";
		for(uint i=0;i<spt.size();i++)
		{
			fout<<spt[i][0]<<" "<<spt[i][1]<<" "<<spt[i][2]<<"\n";
		}
		fout<<"\nCELLS "<<sele.size()<<" "<<5*sele.size()<<'\n';
		for(uint i=0;i<sele.size();i++)
		{
			fout<<"4 "<<sele[i][0]<<" "<<sele[i][1]<<" "<<sele[i][2]<<" "<<sele[i][3]<<'\n';
		}
		fout<<"\nCELL_TYPES "<<sele.size()<<'\n';
		for(uint i=0;i<sele.size();i++)
		{
			fout<<"9\n";
		}
		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i]<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout<<"Cannot open "<<fname<<"!\n";
	}
}

void TruncatedTspline_2D::run_surf_XP(string fn, vector<int> ids)
{
	// SetProblem_surf_XP(fn + "controlmesh");
	SetProblem_surf_XP(fn + "controlmesh_initial");
	Refine_Surf_Test_0(ids);
	SetBezierMatIrrPatch();
	Truncation();
	//CollectActives();
	// VisualizeTMesh(fn + "Tmesh_1");
	// VisualizeSurface(fn + "controlmesh_surf");
	// VisualizeControlMesh(fn + "controlmesh_locallyrefined");
	VisualizeControlMesh(fn + "controlmesh");

	vector<BezierElement2D> bzmesh;
	BezierExtract_Unstruct(bzmesh);
	// BezierControlMesh_Unstruct(fn + "bzmesh_2",bzmesh);
	// BezierVTK_Unstruct(fn + "bzmesh_3",bzmesh);

	OutputMesh(bzmesh, fn);
}

void TruncatedTspline_2D::OutputMesh(const vector<BezierElement2D>& bzmesh, string fn)
{
	// int cn[8] = { 0, 3, 15, 12, 48, 51, 63, 60 };
	int cn[4] = { 0, 3, 15, 12};
	string fname = fn + "bzmesh.vtk";
	ofstream fout;
	fout.open(fname.c_str(), std::ofstream::trunc);
	// fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nBezier mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << 4 * bzmesh.size() << " float\n";
		for (uint i = 0; i<bzmesh.size(); i++)
		{
			for (int j = 0; j < 4; j++)
			{
				fout << bzmesh[i].pts[cn[j]][0] << " " << bzmesh[i].pts[cn[j]][1] << " " << bzmesh[i].pts[cn[j]][2] << "\n";
			}
		}
		fout << "\nCELLS " << bzmesh.size() << " " << 5 * bzmesh.size() << '\n';
		for (uint i = 0; i<bzmesh.size(); i++)
		{
			fout << "4 " << 4 * i << " " << 4 * i + 1 << " " << 4 * i + 2 << " " << 4 * i + 3 << '\n';
		}
		fout << "\nCELL_TYPES " << bzmesh.size() << '\n';
		for (uint i = 0; i<bzmesh.size(); i++)
		{
			fout << "9\n";
		}
		// //fout << "POINT_DATA " << sdisp.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		// //for (uint i = 0; i<sdisp.size(); i++)
		// //{
		// //	fout << sdisp[i] << "\n";
		// //}
		// fout << "\nCELL_DATA " << bzmesh.size() << "\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		// for (uint i = 0; i < bzmesh.size(); i++)
		// {
		// 	// fout << bzmesh[i].type << "\n";
		// 	fout << "9" << "\n";
		// }
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
	string fname3(fn + "bzmeshinfo.txt");
	//ofstream fout;
	fout.open(fname3.c_str(), std::ofstream::trunc);
	// fout.open(fname3.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() << "\n";
		for (uint i = 0; i<bzmesh.size(); i++)
		{
			for (uint l = 0; l < bzmesh[i].IEN.size(); l++)
			{
				fout << bzmesh[i].IEN[l] + 1;
				if (l == bzmesh[i].IEN.size() - 1)
				{
					fout << "\n";
				}
				else
				{
					fout << " ";
				}
			}
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname3 << '\n';
	}

	string fname1(fn + "cmat.txt");
	//ofstream fout;
	fout.open(fname1.c_str(), std::ofstream::trunc);
	// fout.open(fname1.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() << "\n";
		for (uint i = 0; i<bzmesh.size(); i++)
		{
			// fout << i << " " << bzmesh[i].IEN.size() << " " << bzmesh[i].type << "\n";
			fout << i << " " << bzmesh[i].IEN.size() << " " << "9" << "\n";
			for (uint l = 0; l < bzmesh[i].IEN.size(); l++)
			{
				fout << bzmesh[i].IEN[l];
				if (l == bzmesh[i].IEN.size() - 1)
				{
					fout << "\n";
				}
				else
				{
					fout << " ";
				}
			}
			for (uint j = 0; j < bzmesh[i].cmat.size(); j++)
			{
				for (uint k = 0; k < bzmesh[i].cmat[j].size(); k++)
				{
					fout << bzmesh[i].cmat[j][k];
					if (k == bzmesh[i].cmat[j].size() - 1)
					{
						fout << "\n";
					}
					else
					{
						fout << " ";
					}
				}
			}

		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname1 << '\n';
	}

	string fname2(fn + "bzpt.txt");
	//ofstream fout;
	fout.open(fname2.c_str(), std::ofstream::trunc);
	// fout.open(fname2.c_str());
	if (fout.is_open())
	{
		// fout << bzmesh.size() * 64 << "\n";
		fout << bzmesh.size() * 16 << "\n";
		for (uint i = 0; i<bzmesh.size(); i++)
		{
			// for (int j = 0; j < 64; j++)
			for (int j = 0; j < 16; j++)
			{
				fout << bzmesh[i].pts[j][0] << " " << bzmesh[i].pts[j][1] << " " << bzmesh[i].pts[j][2] << "\n";
			}
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname2 << '\n';
	}
}

void TruncatedTspline_2D::runXP_Laplace()
{
	SetLshapeProblem_XP("Lshape_new/input_CM_XP");
	//InitialConnect();
	//int rf_id0[2]={25,26};
	//int rf_type0[2]={0,0};
	//vector<int> rf_id(rf_id0,rf_id0+1);
	//vector<int> rf_type(rf_type0,rf_type0+1);
	//Refine_Unstruct(rf_id,rf_type);
	//LshapeRefine_XP();
	SetBezierMatIrrPatch();
	CollectActives();
	cout<<"DOF: "<<cp.size()<<"\n";
	//VisualizeTMesh("Lshape_new/Tmesh_XP_1");
	//VisualizeControlMesh("Lshape_new/CM_XP_4");
	VisualizeSurface("Lshape_new/surf_XP_4");

	//vector<BezierElement> bzmesh;
	//BezierExtract_Unstruct(bzmesh);
	//BezierControlMesh_Unstruct("test8/Lshape_bzmesh_CM_2.vtk",bzmesh);
	//BezierVTK_Unstruct("test8/Lshape_bzmesh_2.vtk",bzmesh);
}

void TruncatedTspline_2D::runXP_complex(string fn)
{
	SetProblem_complex(fn);
	CollectActives();
	VisualizeControlMesh("../iotest/output");
}

void TruncatedTspline_2D::SetBezier3TranMat(int N, vector<vector<double>>& bmat)
{
	int nb=2*N+8;
	bmat.resize(nb,vector<double>(16,0.));
	double a(4./9.), bN(4./(9.*double(N))), cN(1./(9.*double(N))), b(2./9.), c(1./9.), d(1./18.), e(1./36);

	//1 extraordinary node coef
	bmat[0][0]=a;
	for(int i=0; i<N; i++)
	{
		bmat[2*i+1][0]=bN;
		bmat[2*i+2][0]=cN;
	}
	//3 regular corner coefs
	int qv[3]={4,16,13};
	int pv[3][9]={{6,5,2*N+3,2*N+4,2*N+5,7,8,1,4},{5,4,2*N+7,2*N+6,2*N+2,2*N+3,2*N+4,6,1},{4,1,2,3,2*N+8,2*N+7,2*N+6,5,6}};
	if(N==3) pv[0][6]=2;
	double dv[9]={a,c,e,c,e,c,e,c,e};
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<9; j++)
		{
			bmat[pv[i][j]-1][qv[i]-1]=dv[j];
		}
	}
	//4 face coefs
	int qf[4]={6,7,10,11};
	int pf[4][4]={{1,6,4,5},{6,5,1,4},{4,1,5,6},{5,4,6,1}};
	double df[4]={a,b,b,c};
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
		{
			bmat[pf[i][j]-1][qf[i]-1]=df[j];
		}
	}
	//8 edge coefs
	int qe[8]={2,3,5,9,8,12,14,15};
	int pe[8][6]={{1,6,4,8,5,7},{6,1,5,7,4,8},{1,4,2,6,3,5},{4,1,3,5,2,6},{6,5,1,2*N+4,4,2*N+3},{5,6,4,2*N+3,1,2*N+4},{4,5,1,2*N+7,6,2*N+6},{5,4,6,2*N+6,1,2*N+7}};
	if(N==3)
	{
		pe[0][3]=2; pe[1][5]=2;
	}
	double de[6]={a,b,c,c,d,d};
	for(int i=0; i<8; i++)
	{
		for(int j=0; j<6; j++)
		{
			bmat[pe[i][j]-1][qe[i]-1]=de[j];
		}
	}

	//update first point as Catmull-Clark limit point
	SingularPatchEval cceigen;
	cceigen.ReadEigenStruct();
	vector<double> cc_coef;
	cceigen.ObtainIVcoef(N,cc_coef);
	for(uint i=0; i<cc_coef.size(); i++)
	{
		bmat[i][0]=cc_coef[i];
	}

	//vector<double> cc_coef0;
	//cceigen.ObtainIVcoef(N,cc_coef0);
	/*vector<vector<double>> cc_coef1(4,vector<double>(2*N+1));
	cceigen.ObtainIVcoef(N,cc_coef1[0]);
	double uv[3][2]={{1.,0.},{1.,1.},{0.,1.}};
	for(uint i=1; i<cc_coef1.size(); i++)
	{
		for(uint j=0; j<cc_coef1[i].size(); j++)
		{
			cc_coef1[i][j]=cceigen.ExplicitBasis(j,uv[i-1][0],uv[i-1][1],N);
			if(cc_coef1[i][j]<1.e-10) cc_coef1[i][j]=0.;
		}
	}
	vector<vector<double>> cc_coef2(4,vector<double>(2*N+1,0.));
	int bzfid[4]={0,1,4,5};
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<2*N+1; j++)
		{
			for(int k=0; k<2*N+1; k++)
			{
				cc_coef2[i][j]+=bmat[j][bzfid[i]]*cc_coef1[i][k];
			}
		}
	}
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<2*N+1; j++)
		{
			bmat[j][bzfid[i]]=cc_coef2[i][j];
		}
	}*/

	//vector<vector<double>> cc_mat;
	//cceigen.ObtainIVMat(N,cc_mat);
	//for(uint i=0; i<cc_mat.size(); i++)
	//{
	//	for(uint j=0; j<cc_mat[i].size(); j++)
	//	{
	//		cout<<cc_mat[i][j]<<" ";
	//	}
	//	cout<<"\n";
	//}
	//getchar();
	//vector<vector<double>> bmat1(16,vector<double>(nb,0.));
	//for(uint i=0; i<bmat1.size(); i++)
	//{
	//	for(uint j=0; j<bmat1[i].size(); j++)
	//	{
	//		for(uint k=0; k<cc_mat[j].size(); k++)
	//		{
	//			bmat1[i][j]+=bmat[j][i]*cc_mat[j][k];
	//		}
	//	}
	//}
	//for(uint i=0; i<bmat.size(); i++)
	//{
	//	for(uint j=0; j<bmat[i].size(); j++)
	//	{
	//		bmat[i][j]=bmat1[j][i];
	//	}
	//}
}

void TruncatedTspline_2D::SetBezier4TranMat(int N, vector<vector<double>>& bmat)
{
	int p(3);
	vector<vector<double>> b3mat;
	SetBezier3TranMat(N,b3mat);
	//degree elevation
	vector<vector<double>> demat(25,vector<double>(16,0.));
	int loc1(0);
	for(int j=0; j<5; j++)
	{
		for(int i=0; i<5; i++)
		{
			double a(double(i)/double(p+1)), b(double(j)/double(p+1));
			double coef[4]={(1.-a)*(1.-b),a*(1.-b),(1.-a)*b,a*b};
			int loc0[4]={4*j+i,4*j+i-1,4*(j-1)+i,4*(j-1)+i-1};
			if(i==0)
			{
				loc0[1]=-1; loc0[3]=-1;
			}
			if(j==0)
			{
				loc0[2]=-1; loc0[3]=-1;
			}
			if(i==4)
			{
				loc0[0]=-1; loc0[2]=-1;
			}
			if(j==4)
			{
				loc0[0]=-1; loc0[1]=-1;
			}
			for(int k=0; k<4; k++)
			{
				if(loc0[k]!=-1)
				{
					demat[loc1][loc0[k]]=coef[k];
				}
			}
			loc1++;
		}
	}
	bmat.resize(2*N+8,vector<double>(25,0.));
	for(uint i=0; i<bmat.size(); i++)
	{
		for(uint j=0; j<bmat[i].size(); j++)
		{
			for(int k=0; k<16; k++)
			{
				bmat[i][j]+=b3mat[i][k]*demat[j][k];
			}
		}
	}
	//cout<<"b3mat:\n";
	//for(int i=0; i<16; i++)
	//{
	//	cout<<b3mat[1][i]<<" ";
	//	if(i%4==3) cout<<"\n";
	//}
	//cout<<"bmat:\n";
	//for(int i=0; i<25; i++)
	//{
	//	cout<<bmat[1][i]<<" ";
	//	if(i%5==4) cout<<"\n";
	//}
	//getchar();
}

void TruncatedTspline_2D::SetBezier4TranMatOP(int N, vector<vector<double>>& bmat)
{
	bmat.clear();
	ifstream fin;
	fin.open("../src/bmat1/bmat_"+to_string((long long)(N))+".txt");
	if(fin.is_open())
	{
		int dim1, dim2;
		fin>>dim1>>dim2;
		bmat.resize(dim1,vector<double>(dim2));
		for(int i=0; i<dim1; i++)
		{
			for(int j=0; j<dim2; j++)
			{
				fin>>bmat[i][j];
			}
		}
	}
	else
	{
		cerr<<"Cannot open bmat file!\n";
	}
}

void TruncatedTspline_2D::VisualizeSurface(string fn)
{
	vector<array<double,3>> spt;
	vector<array<double,3>> sval;
	vector<array<int,4>> sele;
	vector<array<double,3>> lpt;//visulize parameter lines
	vector<array<int,2>> led;//line connectivity
	int ns(5),ecount(0),loc0,loc1,loc2;
	//vector<double> su(ns),sv(ns);
	//for(int i=0; i<ns; i++)
	//{
	//	su[i]=i*1./(ns-1);
	//	sv[i]=i*1./(ns-1);
	//}

	for(uint e=0;e<tmesh.size();e++)
	{
		if(tmesh[e].act==1 && (tmesh[e].type==0 || tmesh[e].type==1 || tmesh[e].type==4))
		{
			int loc(0);
			vector<double> su(ns),sv(ns);
			for(int i=0; i<ns; i++)
			{
				su[i]=i*tmedge[tmesh[e].edge[0]].len/(ns-1);
				sv[i]=i*tmedge[tmesh[e].edge[3]].len/(ns-1);
			}

			for(int a=0;a<ns;a++)
			{
				for(int b=0;b<ns;b++)
				{
					array<double,3> pt;
					array<double,3> nm;
					SurfacePointMap(e,su[b],sv[a],pt,nm);
					spt.push_back(pt);
					sval.push_back(nm);
					if(a==0||a==ns-1||b==0||b==ns-1)
					{
						lpt.push_back(pt);
					}
				}
			}

			for(int a=0;a<ns-1;a++)
			{
				for(int b=0;b<ns-1;b++)
				{
					array<int,4> el;
					el[0]=ecount*ns*ns+a*ns+b;
					el[1]=ecount*ns*ns+a*ns+b+1;
					el[2]=ecount*ns*ns+(a+1)*ns+b+1;
					el[3]=ecount*ns*ns+(a+1)*ns+b;
					sele.push_back(el);
				}
			}
			for(int a=0;a<ns-1;a++)
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
			led.push_back(lc1);
			ecount++;
		}
	}

	string fname=fn+".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if(fout.is_open())
	{
		fout<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout<<"POINTS "<<spt.size()<<" float\n";
		for(uint i=0;i<spt.size();i++)
		{
			fout<<spt[i][0]<<" "<<spt[i][1]<<" "<<spt[i][2]<<"\n";
		}
		fout<<"\nCELLS "<<sele.size()<<" "<<5*sele.size()<<'\n';
		for(uint i=0;i<sele.size();i++)
		{
			fout<<"4 "<<sele[i][0]<<" "<<sele[i][1]<<" "<<sele[i][2]<<" "<<sele[i][3]<<'\n';
		}
		fout<<"\nCELL_TYPES "<<sele.size()<<'\n';
		for(uint i=0;i<sele.size();i++)
		{
			fout<<"9\n";
		}
		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i]<<"\n";
		//}
		fout<<"\nPOINT_DATA "<<sval.size()<<"\nNORMALS Normal FLOAT\n";
		for(uint i=0;i<sval.size();i++)
		{
			fout<<sval[i][0]<<" "<<sval[i][1]<<" "<<sval[i][2]<<"\n";
		}
		fout.close();
	}
	else
	{
		cout<<"Cannot open "<<fname<<"!\n";
	}

	string fname1(fn+"-lines.vtk");
	ofstream fout1;
	fout1.open(fname1.c_str());
	if(fout1.is_open())
	{
		fout1<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout1<<"POINTS "<<lpt.size()<<" float\n";
		for(uint i=0;i<lpt.size();i++)
		{
			fout1<<lpt[i][0]<<" "<<lpt[i][1]<<" "<<lpt[i][2]<<"\n";
		}
		fout1<<"\nCELLS "<<led.size()<<" "<<3*led.size()<<'\n';
		for(uint i=0;i<led.size();i++)
		{
			fout1<<"2 "<<led[i][0]<<" "<<led[i][1]<<'\n';
		}
		fout1<<"\nCELL_TYPES "<<led.size()<<'\n';
		for(uint i=0;i<led.size();i++)
		{
			fout1<<"3\n";
		}
		fout1.close();
	}
	else
	{
		cout<<"Cannot open "<<fname1<<"!\n";
	}
}

void TruncatedTspline_2D::CapIndex_Loc2Glb(int N, vector<vector<int>>& ICN)
{
	ICN.resize(N,vector<int>(25));
	for(int i=0; i<N-1; i++)
	{
		ICN[i][0]=0;
		int loc(1);
		for(int j=1; j<25; j++)
		{
			int a(j%5), b(j/5);
			if(a==0 && b!=0) ICN[i][j]=20*(i+1)+b;
			else
			{
				ICN[i][j]=20*i+loc;
				loc++;
			}
		}
	}
	ICN[N-1][0]=0;
	int loc(1);
	for(int j=1; j<25; j++)
	{
		int a(j%5), b(j/5);
		if(a==0 && b!=0) ICN[N-1][j]=b;
		else
		{
			ICN[N-1][j]=20*(N-1)+loc;
			loc++;
		}
	}
}

void TruncatedTspline_2D::BuildGeomMat(int N, const vector<int>& bloc, const vector<vector<double>>& bmat, const vector<vector<int>>& ICN, vector<vector<double>>& gmat, vector<double>& gvec)
{
	double PI(3.1415926535);
	double lamda(-2.*cos(2.*PI/N));
	int fix[14]={3,4,8,9,13,14,18,19,23,24,16,17,21,22};
	gmat.resize(17*N+2,vector<double>(20*N+1,0.));
	gvec.resize(17*N+2,0.);
	int rs(0);
	//gmat[rs][ICN[N-1][1]]=1.; gmat[rs][ICN[0][0]]=-2.-lamda; gmat[rs][ICN[0][1]]=lamda; gmat[rs][ICN[0][5]]=1.;
	//rs++;
	gmat[rs][ICN[0][0]]=-double(N);
	for(int i=0; i<N; i++)
	{
		gmat[rs][ICN[i][1]]=1.;
	}
	rs++;
	//gmat[rs][ICN[0][0]]=-double(N);
	//for(int i=0; i<N; i++)
	//{
	//	gmat[rs][ICN[i][6]]=1.;
	//}
	//rs++;
	gmat[rs][0]=1.;
	gvec[rs]=bmat[bloc[0]][0];
	rs++;
	for(int i=0; i<N; i++)//loop each face
	{
		for(int j=0; j<14; j++)// 14 fixed coefficients
		{
			gmat[rs][ICN[i][fix[j]]]=1.;
			gvec[rs]=bmat[bloc[i]][fix[j]];
			rs++;
		}
		//6 boundary equations
		int i0((i-1+N)%N);
		gmat[rs][ICN[i][0]]=1.; gmat[rs][ICN[i][1]]=-4.; gmat[rs][ICN[i][2]]=6.; gmat[rs][ICN[i][3]]=-4.; gmat[rs][ICN[i][4]]=1.;
		rs++;
		//gmat[rs][ICN[i0][1]]=1.; gmat[rs][ICN[i][0]]=-2.-lamda; gmat[rs][ICN[i][1]]=lamda; gmat[rs][ICN[i][5]]=1.;
		//rs++;
		gmat[rs][ICN[i0][6]]=4.; gmat[rs][ICN[i][1]]=-8.-2.*lamda; gmat[rs][ICN[i][0]]=lamda/2.; gmat[rs][ICN[i][3]]=2.*lamda; gmat[rs][ICN[i][4]]=-lamda/2.; gmat[rs][ICN[i][6]]=4.;
		rs++;
		gmat[rs][ICN[i0][11]]=4.; gmat[rs][ICN[i][2]]=-8.; gmat[rs][ICN[i][4]]=2.*lamda/3.; gmat[rs][ICN[i][3]]=-2.*lamda/3.; gmat[rs][ICN[i][7]]=4.;
		rs++;
		//gmat[rs][ICN[i0][16]]=1.; gmat[rs][ICN[i][3]]=-2.; gmat[rs][ICN[i][8]]=1.;
		//rs++;
		//gmat[rs][ICN[i0][21]]=1.; gmat[rs][ICN[i][4]]=-2.; gmat[rs][ICN[i][9]]=1.;
		//rs++;
		//gmat[rs][ICN[i][0]]=1.; gmat[rs][ICN[i][1]]=-4.; gmat[rs][ICN[i][2]]=6.; gmat[rs][ICN[i][3]]=-4.; gmat[rs][ICN[i][4]]=1.;
		//rs++;

		//gmat[rs][ICN[i0][1]]=1.; gmat[rs][ICN[i][0]]=-2.-lamda; gmat[rs][ICN[i][1]]=lamda; gmat[rs][ICN[i][5]]=1.;
		//rs++;
		//gmat[rs][ICN[i0][6]]=4.; gmat[rs][ICN[i][1]]=-8.-2.*lamda; gmat[rs][ICN[i][0]]=lamda/2.; /*gmat[rs][ICN[i][3]]=2.*lamda; gmat[rs][ICN[i][4]]=-lamda/2.;*/ gmat[rs][ICN[i][6]]=4.;
		//gvec[rs]=-2.*lamda*bmat[bloc[i]][3]+lamda/2.*bmat[bloc[i]][4];
		//rs++;
		//gmat[rs][ICN[i0][11]]=4.; gmat[rs][ICN[i][2]]=-8.; /*gmat[rs][ICN[i][4]]=2.*lamda/3.; gmat[rs][ICN[i][3]]=-2.*lamda/3.;*/ gmat[rs][ICN[i][7]]=4.;
		//gvec[rs]=-2.*lamda/3.*(bmat[bloc[i]][4]-bmat[bloc[i]][3]);
		//rs++;
		////gmat[rs][ICN[i0][16]]=1.; gmat[rs][ICN[i][3]]=-2.; gmat[rs][ICN[i][8]]=1.;
		////rs++;
		////gmat[rs][ICN[i0][21]]=1.; gmat[rs][ICN[i][4]]=-2.; gmat[rs][ICN[i][9]]=1.;
		////rs++;
		//gmat[rs][ICN[i][0]]=1.; gmat[rs][ICN[i][1]]=-4.; gmat[rs][ICN[i][2]]=6.; /*gmat[rs][ICN[i][3]]=-4.; gmat[rs][ICN[i][4]]=1.;*/
		//gvec[rs]=4.*bmat[bloc[i]][3]-bmat[bloc[i]][4];
		//rs++;
	}

	//int fix1[5]={1,2,6,7,11};
	//for(int i=0; i<N; i++)
	//{
	//	for(int j=0; j<5; j++)
	//	{
	//		if(bmat[bloc[i]][fix1[j]]==0.)
	//		{
	//			vector<double> gtmp(20*N+1,0.);
	//			gtmp[ICN[i][fix1[j]]]=1.;
	//			gmat.push_back(gtmp);
	//			gvec.push_back(0.);
	//		}
	//	}
	//}
}

void TruncatedTspline_2D::BuildFairMat(int N, const vector<int>& bloc, const vector<vector<double>>& bmat, const vector<vector<int>>& ICN, vector<vector<double>>& fmat, vector<double>& fvec)
{
	fmat.resize(40*N,vector<double>(20*N+1,0.));
	fvec.resize(40*N,0.);
	int count(0);
	for(int i=0; i<N; i++)//loop each face
	{
		for(int k=0; k<5; k++)
		{
			for(int j=0; j<4; j++)
			{
				int loc(5*k+j), loc1(5*k+j+1);
				fmat[count][ICN[i][loc]]=1.; fmat[count][ICN[i][loc1]]=-1.;
				fvec[count]=bmat[bloc[i]][loc]-bmat[bloc[i]][loc1];
				count++;
			}
		}
		for(int k=0; k<4; k++)
		{
			for(int j=0; j<5; j++)
			{
				int loc(5*k+j), loc1(5*(k+1)+j);
				fmat[count][ICN[i][loc]]=1.; fmat[count][ICN[i][loc1]]=-1.;
				fvec[count]=bmat[bloc[i]][loc]-bmat[bloc[i]][loc1];
				count++;
			}
		}
	}
}

void TruncatedTspline_2D::Convert2MatlabData(const vector<vector<double>>& mat, vector<double>& row_id, vector<double>& col_id, vector<double>& coef)
{
	for(uint j=0; j<mat[0].size(); j++)
	{
		for(uint i=0; i<mat.size(); i++)
		{
			if(mat[i][j]!=0.)
			{
				row_id.push_back(double(i+1));
				col_id.push_back(double(j+1));
				coef.push_back(mat[i][j]);
			}
		}
	}
}

void TruncatedTspline_2D::PreConditionCoef(int N, const vector<int>& bloc, const vector<vector<double>>& bmat, const vector<vector<int>>& ICN, vector<double>& coef)
{
	coef.resize(20*N+1,0.);
	for(uint i=0; i<ICN.size(); i++)
	{
		for(uint j=0; j<ICN[i].size(); j++)
		{
			coef[ICN[i][j]]=bmat[bloc[i]][j];
		}
	}
}

//void TruncatedTspline_2D::SolveOptimizeBezierMat(int N, vector<vector<double>>& bmatop)
//{
//	vector<vector<int>> ICN;
//	vector<vector<double>> bmat;
//	CapIndex_Loc2Glb(N,ICN);
//	SetBezier4TranMat(N,bmat);
//
//	bmatop.resize(bmat.size(),vector<double>(bmat[0].size()));
//	//ofstream fout;
//	//fout.open("bmat1/bmat_5_0.txt");
//	for(uint i=0; i<bmat.size(); i++)
//	{
//		for(uint j=0; j<bmat[i].size(); j++)
//		{
//			bmatop[i][j]=bmat[i][j];
//			//fout<<bmat[i][j]<<" ";
//		}
//		//fout<<"\n";
//	}
//	//fout.close();
//	//getchar();
//	vector<vector<int>> G2L(2*N+1,vector<int>(N));
//	vector<vector<int>> L2G(N,vector<int>(2*N+1));
//	for(int i=0; i<N; i++)
//	{
//		L2G[i][0]=0;
//		G2L[0][i]=0;
//		for(int j=1; j<2*N+1; j++)
//		{
//			L2G[i][j]=(-2*i+j-1+2*N)%(2*N)+1;
//			G2L[L2G[i][j]][i]=j;
//		}
//	}
//
//	for(int bfid=0; bfid<3; bfid++)//optimize for one-ring basis functions
//	{
//		cout<<"Basis function: "<<bfid<<"\n";
//		vector<vector<double>> gmat,fmat;
//		vector<double> gvec,fvec;
//		BuildGeomMat(N,G2L[bfid],bmat,ICN,gmat,gvec);
//		BuildFairMat(N,G2L[bfid],bmat,ICN,fmat,fvec);
//
//		vector<double> gid1,gid2,gc,fid1,fid2,fc,precf;
//		Convert2MatlabData(gmat,gid1,gid2,gc);
//		Convert2MatlabData(fmat,fid1,fid2,fc);
//		PreConditionCoef(N,G2L[bfid],bmat,ICN,precf);
//
//		MatlabSolver solver;
//		cout<<"solving\n";
//		vector<double> sol(20*N+1);
//		solver.Initilize2(gid1,gid2,gc,gvec,fid1,fid2,fc,fvec,precf);
//		solver.Solve3(sol.data());
//
//		for(int i=0; i<25; i++)
//		{
//			if(bfid==0)
//			{
//				bmatop[0][i]=sol[ICN[0][i]];
//			}
//			else if(bfid==1)
//			{
//				for(int j=0; j<N; j++)
//				{
//					int rowid=2*j+1;
//					bmatop[rowid][i]=sol[ICN[j][i]];
//					if(bmatop[rowid][i]<1.e-10) bmatop[rowid][i]=0.;
//				}
//			}
//			else if(bfid==2)
//			{
//				for(int j=0; j<N; j++)
//				{
//					int rowid=2*j+2;
//					bmatop[rowid][i]=sol[ICN[j][i]];
//					if(bmatop[rowid][i]<1.e-10) bmatop[rowid][i]=0.;
//				}
//			}
//		}
//	}
//
//	for(uint i=0; i<bmatop[0].size(); i++)
//	{
//		double sum(0.);
//		for(uint j=0; j<bmatop.size(); j++)
//		{
//			sum+=bmatop[j][i];
//		}
//		for(uint j=0; j<bmatop.size(); j++)
//		{
//			bmatop[j][i]/=sum;
//		}
//	}
//
//	ofstream fout;
//	fout.open("bmat1/bmat_"+to_string((long long)(N))+".txt");
//	if(fout.is_open())
//	{
//		fout<<bmatop.size()<<" "<<bmatop[0].size()<<"\n";
//		for(uint i=0; i<bmatop.size(); i++)
//		{
//			for(uint j=0; j<bmatop[i].size(); j++)
//			{
//				fout<<bmatop[i][j]<<" ";
//			}
//			fout<<"\n";
//		}
//		fout.close();
//	}
//	else
//	{
//		cerr<<"Cannot open bmat file!\n";
//	}
//}

void TruncatedTspline_2D::SetProblem_surf_XP(string fn)
{
	string fname(fn+".vtk"),stmp;
	int npts,neles,itmp;
	ifstream fin;
	fin.open(fname);
	if(fin.is_open())
	{
		for(int i=0;i<4;i++) getline(fin,stmp);//skip lines
		fin>>stmp>>npts>>stmp;
		cp.resize(npts);
		for(int i=0;i<npts;i++)
		{
			fin>>cp[i].coor[0]>>cp[i].coor[1]>>cp[i].coor[2];
		}
		getline(fin,stmp);
		fin>>stmp>>neles>>itmp;
		tmesh.resize(neles);
		for(int i=0;i<neles;i++)
		{
			fin>>itmp>>tmesh[i].cnct[0]>>tmesh[i].cnct[1]>>tmesh[i].cnct[2]>>tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr<<"Cannot open "<<fname<<"!\n";
	}

	InitialConnect();
	FindEdgeTopoDirec();
	FindKnotInterval();
	UpdateKnotInterval();
	SetLocalCoorSystem();
	FindIEN_Invalid();
	Update_IEN();
}

void TruncatedTspline_2D::SetLshapeProblem_XP(string fn)
{
	//read quad vtk
	string fname(fn+".vtk"),stmp;
	int npts,neles,itmp;
	ifstream fin;
	fin.open(fname);
	if(fin.is_open())
	{
		for(int i=0;i<4;i++) getline(fin,stmp);//skip lines
		fin>>stmp>>npts>>stmp;
		cp.resize(npts);
		for(int i=0;i<npts;i++)
		{
			fin>>cp[i].coor[0]>>cp[i].coor[1]>>cp[i].coor[2];
		}
		getline(fin,stmp);
		fin>>stmp>>neles>>itmp;
		tmesh.resize(neles);
		for(int i=0;i<neles;i++)
		{
			fin>>itmp>>tmesh[i].cnct[0]>>tmesh[i].cnct[1]>>tmesh[i].cnct[2]>>tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr<<"Cannot open "<<fname<<"!\n";
	}

	//set T-spline info
	InitialConnect();
	//VisualizeTMesh("test8/Lshape_Tmesh_1");
	int edge_zero[26]={15,17,19,26,28,53,55,62,64,157,159,166,168,212,214,221,243,247,249,255,275,279,349,355,375,379};
	//int el_type3[4]={5,74,96,175};
	for(int i=0; i<26; i++)
	{
		tmedge[edge_zero[i]].len=0.;
	}
	//for(int i=0; i<4; i++)
	//{
	//	tmesh[el_type3[i]].type=3;
	//}
	for(uint i=0; i<tmesh.size(); i++)
	{
		if(tmedge[tmesh[i].edge[0]].len==0. || tmedge[tmesh[i].edge[1]].len==0.)
		{
			tmesh[i].type=2;
		}
		if(tmedge[tmesh[i].edge[0]].len==0. && tmedge[tmesh[i].edge[1]].len==0.)
		{
			tmesh[i].type=3;
		}
	}
	FindEdgeTopoDirec();
	FindKnotInterval();
	UpdateKnotInterval();
	SetLocalCoorSystem();
	FindIEN_Unstruct();
}

void TruncatedTspline_2D::SetProblem_complex(string fn)
{
	//read quad vtk
	string fname(fn+".vtk"),stmp;
	int npts,neles,itmp;
	ifstream fin;
	fin.open(fname);
	if(fin.is_open())
	{
		for(int i=0;i<4;i++) getline(fin,stmp);//skip lines
		fin>>stmp>>npts>>stmp;
		cp.resize(npts);
		for(int i=0;i<npts;i++)
		{
			fin>>cp[i].coor[0]>>cp[i].coor[1]>>cp[i].coor[2];
		}
		getline(fin,stmp);
		fin>>stmp>>neles>>itmp;
		tmesh.resize(neles);
		for(int i=0;i<neles;i++)
		{
			fin>>itmp>>tmesh[i].cnct[0]>>tmesh[i].cnct[1]>>tmesh[i].cnct[2]>>tmesh[i].cnct[3];
		}
		fin.close();
		std::cout << "Read " << fname << std::endl;
	}
	else
	{
		cerr<<"Cannot open "<<fname<<"!\n";
	}

	InitialConnect();
	FindEdgeTopoDirec();
	FindKnotInterval();
	UpdateKnotInterval();
	SetLocalCoorSystem();
	FindIEN_Invalid();
	Update_IEN();
	//Validate_Tmesh();
	//CollectActives();
}

void TruncatedTspline_2D::SetBounaryPoints(vector<int>& pid, vector<double>& disp, int& ncp)//fertility
{
	ncp=cp.size();
	int nebc1(9);
	int ebcid1[]={293,294,295,441,442,274,277,279,292};//disp 100
	pid.clear();
	disp.clear();
	for(int i=0; i<nebc1; i++)
	{
		if(tmesh[ebcid1[i]].act==0)
		{
			for(int k=0; k<4; k++)
			{
				if(tmesh[ebcid1[i]].chd[k]!=-1)
				{
					int echd(tmesh[ebcid1[i]].chd[k]);
					for(int j=0; j<4; j++)
					{
						vector<int>::iterator it=find(pid.begin(),pid.end(),tmesh[echd].cnct[j]);
						if(it==pid.end())
						{
							pid.push_back(tmesh[echd].cnct[j]);
							disp.push_back(100.);
						}
					}
				}
			}
		}
		else
		{
			for(int j=0; j<4; j++)
			{
				vector<int>::iterator it=find(pid.begin(),pid.end(),tmesh[ebcid1[i]].cnct[j]);
				if(it==pid.end())
				{
					pid.push_back(tmesh[ebcid1[i]].cnct[j]);
					disp.push_back(100.);
				}
			}
		}
	}
	int nebc2(6);
	int ebcid2[]={341,343,344,345,346,347};//disp 0
	for(int i=0; i<nebc2; i++)
	{
		if(tmesh[ebcid2[i]].act==0)
		{
			for(int k=0; k<4; k++)
			{
				if(tmesh[ebcid2[i]].chd[k]!=-1)
				{
					int echd(tmesh[ebcid2[i]].chd[k]);
					for(int j=0; j<4; j++)
					{
						vector<int>::iterator it=find(pid.begin(),pid.end(),tmesh[echd].cnct[j]);
						if(it==pid.end())
						{
							pid.push_back(tmesh[echd].cnct[j]);
							disp.push_back(0.);
						}
					}
				}
			}
		}
		else
		{
			for(int j=0; j<4; j++)
			{
				vector<int>::iterator it=find(pid.begin(),pid.end(),tmesh[ebcid2[i]].cnct[j]);
				if(it==pid.end())
				{
					pid.push_back(tmesh[ebcid2[i]].cnct[j]);
					disp.push_back(0.);
				}
			}
		}
	}
}

void TruncatedTspline_2D::Identify_Invalid_Elements(vector<int>& rid)
{
	//find two ring elements of XPs
	vector<int> xpid;
	vector<vector<int>> e2r;
	vector<int> loc(cp.size(),-1);
	int count(0);
	for(uint i=0; i<cp.size(); i++)
	{
		if(cp[i].face.size()==3 || cp[i].face.size()>4)
		{
			loc[i]=count++;
			xpid.push_back(i);
			vector<int> etmp;
			vector<int> p1r;
			for(uint j=0; j<cp[i].face.size(); j++)
			{
				etmp.push_back(cp[i].face[j]);
				for(uint k=0; k<4; k++)
				{
					if(tmesh[cp[i].face[j]].cnct[k]!=i)
					{
						vector<int>::iterator it=find(p1r.begin(),p1r.end(),tmesh[cp[i].face[j]].cnct[k]);
						if(it==p1r.end())
						{
							p1r.push_back(tmesh[cp[i].face[j]].cnct[k]);
						}
					}
				}
			}
			for(uint j=0; j<p1r.size(); j++)
			{
				for(uint k=0; k<cp[p1r[j]].face.size(); k++)
				{
					vector<int>::iterator it=find(etmp.begin(),etmp.end(),cp[p1r[j]].face[k]);
					if(it==etmp.end())
					{
						etmp.push_back(cp[p1r[j]].face[k]);
					}
				}
			}
			e2r.push_back(etmp);
		}
	}
	//find invalid element and its neighborhood
	vector<int> xpid0;
	for(uint i=0; i<tmesh.size(); i++)
	{
		if(tmesh[i].type==5)
		{
			for(int j=0; j<4; j++)
			{
				if(cp[tmesh[i].cnct[j]].face.size()==3 || cp[tmesh[i].cnct[j]].face.size()>4)
				{
					vector<int>::iterator it=find(xpid0.begin(),xpid0.end(),loc[tmesh[i].cnct[j]]);
					if(it==xpid0.end())
					{
						xpid0.push_back(loc[tmesh[i].cnct[j]]);
					}
				}
			}
		}
	}
	vector<int> rid0;
	for(uint i=0; i<xpid0.size(); i++)
	{
		for(uint j=0; j<e2r[xpid0[i]].size(); j++)
		{
			vector<int>::iterator it=find(rid0.begin(),rid0.end(),e2r[xpid0[i]][j]);
			if(it==rid0.end())
			{
				rid0.push_back(e2r[xpid0[i]][j]);
			}
		}
	}
	for(uint i=0; i<rid0.size(); i++)
	{
		tmesh[rid0[i]].act=0;
	}
	for(int iter=0; iter<3; iter++)
	{
		for(uint i=0; i<e2r.size(); i++)
		{
			int flag(0);
			for(uint j=0; j<e2r[i].size(); j++)
			{
				if(tmesh[e2r[i][j]].act==0)
				{
					flag=1; break;
				}
			}
			if(flag==1)
			{
				for(uint j=0; j<e2r[i].size(); j++)
				{
					vector<int>::iterator it=find(rid0.begin(),rid0.end(),e2r[i][j]);
					if(it==rid0.end())
					{
						rid0.push_back(e2r[i][j]);
						tmesh[e2r[i][j]].act=0;
					}
				}
			}
		}
	}

	rid.clear();
	rid=rid0;
	rid0.clear();
}

void TruncatedTspline_2D::Identify_More_Elements(vector<int>& rid)
{
	//find two ring elements of XPs
	vector<int> xpid;
	vector<vector<int>> e2r;
	vector<int> loc(cp.size(),-1);
	int count(0);
	for(uint i=0; i<cp.size(); i++)
	{
		if(cp[i].face.size()==3 || cp[i].face.size()>4)
		{
			loc[i]=count++;
			xpid.push_back(i);
			vector<int> etmp;
			vector<int> p1r;
			for(uint j=0; j<cp[i].face.size(); j++)
			{
				etmp.push_back(cp[i].face[j]);
				for(uint k=0; k<4; k++)
				{
					if(tmesh[cp[i].face[j]].cnct[k]!=i)
					{
						vector<int>::iterator it=find(p1r.begin(),p1r.end(),tmesh[cp[i].face[j]].cnct[k]);
						if(it==p1r.end())
						{
							p1r.push_back(tmesh[cp[i].face[j]].cnct[k]);
						}
					}
				}
			}
			for(uint j=0; j<p1r.size(); j++)
			{
				for(uint k=0; k<cp[p1r[j]].face.size(); k++)
				{
					vector<int>::iterator it=find(etmp.begin(),etmp.end(),cp[p1r[j]].face[k]);
					if(it==etmp.end())
					{
						etmp.push_back(cp[p1r[j]].face[k]);
					}
				}
			}
			e2r.push_back(etmp);
		}
	}
	//find invalid element and its neighborhood
	vector<int> xpid0;
	for(uint i=0; i<tmesh.size(); i++)
	{
		if(tmesh[i].type==5)
		{
			for(int j=0; j<4; j++)
			{
				if(cp[tmesh[i].cnct[j]].face.size()==3 || cp[tmesh[i].cnct[j]].face.size()>4)
				{
					vector<int>::iterator it=find(xpid0.begin(),xpid0.end(),loc[tmesh[i].cnct[j]]);
					if(it==xpid0.end())
					{
						xpid0.push_back(loc[tmesh[i].cnct[j]]);
					}
				}
			}
		}
	}
	vector<int> rid0;
	for(uint i=0; i<xpid0.size(); i++)
	{
		for(uint j=0; j<e2r[xpid0[i]].size(); j++)
		{
			vector<int>::iterator it=find(rid0.begin(),rid0.end(),e2r[xpid0[i]][j]);
			if(it==rid0.end())
			{
				rid0.push_back(e2r[xpid0[i]][j]);
			}
		}
	}
	for(uint i=0; i<rid.size(); i++)
	{
		rid0.push_back(rid[i]);
	}
	for(uint i=0; i<rid0.size(); i++)
	{
		tmesh[rid0[i]].act=0;
	}
	for(int iter=0; iter<6; iter++)
	{
		for(uint i=0; i<e2r.size(); i++)
		{
			int flag(0);
			for(uint j=0; j<e2r[i].size(); j++)
			{
				if(tmesh[e2r[i][j]].act==0)
				{
					flag=1; break;
				}
			}
			if(flag==1)
			{
				for(uint j=0; j<e2r[i].size(); j++)
				{
					vector<int>::iterator it=find(rid0.begin(),rid0.end(),e2r[i][j]);
					if(it==rid0.end())
					{
						rid0.push_back(e2r[i][j]);
						tmesh[e2r[i][j]].act=0;
					}
				}
			}
		}
	}

	rid.clear();
	rid=rid0;
	rid0.clear();
}

void TruncatedTspline_2D::Update_IEN()
{
	for(uint eid=0; eid<tmesh.size(); eid++)
	{
		vector<int>().swap(tmesh[eid].IEN);
		vector<array<double,5>>().swap(tmesh[eid].patch_ku);
		vector<array<double,5>>().swap(tmesh[eid].patch_kv);
		tmesh[eid].IEN=tmesh[eid].IENtmp;
		tmesh[eid].patch_ku=tmesh[eid].patch_kutmp;
		tmesh[eid].patch_kv=tmesh[eid].patch_kvtmp;
		vector<int>().swap(tmesh[eid].IENtmp);
		vector<array<double,5>>().swap(tmesh[eid].patch_kutmp);
		vector<array<double,5>>().swap(tmesh[eid].patch_kvtmp);
	}
}

void TruncatedTspline_2D::Refine_Surf_Select_0(vector<int>& rfid, vector<int>& rftype)
{
	rfid.clear();
	rftype.clear();
	int ids[]={210};
	int types[1]={0};
	int nrid(1);

	rfid.assign(ids,ids+nrid);
	rftype.assign(types,types+nrid);
}

void TruncatedTspline_2D::Refine_Surf_Test_0(vector<int> ids)
{
	vector<int> rfid, rftype, rfid1, rftype1;
	rfid.clear();
	rftype.clear();
	
	// Refine_Surf_Select_0(rfid, rftype);

	// rfid = ids;
	// for (int i = 0; i<ids.size(); i++) {
	// 	rftype.push_back(1);
	// }
	for (int i = 0; i<ids.size()/2; i++) {
		rfid.push_back(ids[i]);
		rftype.push_back(ids[i+ids.size()/2]);
		// std::cout << rfid[i] << " " << rftype[i] << std::endl;
	}
 
	// std::cout << "Outputting rfid ##############" << std::endl;
	// for (int i = 0; i < rfid.size(); i++){
	// 	std::cout << rfid[i] << " ";
	// }
	// std::cout << std::endl;
	// std::cout << "End of rfid ##############" << std::endl;

	// std::cout << "Outputting rftype ##############" << std::endl;
	// for (int i = 0; i < rftype.size(); i++){
	// 	std::cout << rftype[i] << " ";
	// }
	// std::cout << std::endl;
	// std::cout << "End of rftype ##############" << std::endl;

	Topo_Refine_Unstruct(rfid,rftype,rfid1,rftype1);
	Geom_Refine_Unstruct(rfid1,rftype1);
}

//void TruncatedTspline_2D::ElementRefine_Unstruct_Topo_Geom_4(int eid)
//{
//	int pid[5];
//	int edid[12];
//	Vertex2D ptmp1;
//	//ptmp1.update=1;
//	ptmp1.coor[0]=(cp[tmesh[eid].cnct[0]].coor[0]+cp[tmesh[eid].cnct[1]].coor[0]+cp[tmesh[eid].cnct[2]].coor[0]+cp[tmesh[eid].cnct[3]].coor[0])/4.;
//	ptmp1.coor[1]=(cp[tmesh[eid].cnct[0]].coor[1]+cp[tmesh[eid].cnct[1]].coor[1]+cp[tmesh[eid].cnct[2]].coor[1]+cp[tmesh[eid].cnct[3]].coor[1])/4.;
//	ptmp1.coor[2]=(cp[tmesh[eid].cnct[0]].coor[2]+cp[tmesh[eid].cnct[1]].coor[2]+cp[tmesh[eid].cnct[2]].coor[2]+cp[tmesh[eid].cnct[3]].coor[2])/4.;
//	cp.push_back(ptmp1);
//	pid[0]=cp.size()-1;
//	for(int j=0; j<4; j++)
//	{
//		int itmp[2]={tmesh[eid].cnct[j],tmesh[eid].cnct[(j+1)%4]};
//		if(tmedge[tmesh[eid].edge[j]].act==1)
//		{
//			Vertex2D ptmp;
//			//ptmp.update=1;
//			int id0[4]={tmesh[eid].cnct[j],tmesh[eid].cnct[(j+1)%4],tmesh[eid].cnct[(j+2)%4],tmesh[eid].cnct[(j+3)%4]};
//			ptmp.coor[0]=(cp[id0[0]].coor[0]+cp[id0[1]].coor[0])/2.;
//			ptmp.coor[1]=(cp[id0[0]].coor[1]+cp[id0[1]].coor[1])/2.;
//			ptmp.coor[2]=(cp[id0[0]].coor[2]+cp[id0[1]].coor[2])/2.;
//			cp.push_back(ptmp);
//			pid[j+1]=cp.size()-1;
//			tmedge[tmesh[eid].edge[j]].midpt=pid[j+1];
//			Edge edtmp1,edtmp2,edtmp3;
//			edtmp1.act=1;
//			edtmp1.pt[0]=itmp[0]; edtmp1.pt[1]=pid[j+1];
//			edtmp1.len=tmedge[tmesh[eid].edge[j]].len/2.;
//			edtmp1.prt=tmesh[eid].edge[j];
//			edtmp2.act=1;
//			edtmp2.pt[0]=pid[j+1]; edtmp2.pt[1]=itmp[1];
//			edtmp2.len=tmedge[tmesh[eid].edge[j]].len/2.;
//			edtmp2.prt=tmesh[eid].edge[j];
//			edtmp3.act=1;
//			edtmp3.pt[0]=pid[0]; edtmp3.pt[1]=pid[j+1];
//			edtmp3.len=tmedge[tmesh[eid].edge[(j+1)%4]].len/2.;
//			tmedge.push_back(edtmp1);
//			edid[3*j]=tmedge.size()-1;
//			tmedge[tmesh[eid].edge[j]].chd[0]=edid[3*j];
//			tmedge.push_back(edtmp2);
//			edid[3*j+1]=tmedge.size()-1;
//			tmedge[tmesh[eid].edge[j]].chd[1]=edid[3*j+1];
//			tmedge.push_back(edtmp3);
//			edid[3*j+2]=tmedge.size()-1;
//			tmedge[tmesh[eid].edge[j]].act=0;
//		}
//		else
//		{
//			pid[j+1]=tmedge[tmesh[eid].edge[j]].midpt;
//			int id1[2]={tmesh[eid].cnct[(j+2)%4],tmesh[eid].cnct[(j+3)%4]};
//			//cp[pid[j+1]].update=1;
//			//if(j==0 || j==3)
//			//{
//				//cp[pid[j+1]].coortmp[0]+=(cp[id1[0]].coor[0]+cp[id1[1]].coor[0])/16.;
//				//cp[pid[j+1]].coortmp[1]+=(cp[id1[0]].coor[1]+cp[id1[1]].coor[1])/16.;
//				//cp[pid[j+1]].coortmp[2]+=(cp[id1[0]].coor[2]+cp[id1[1]].coor[2])/16.;
//			//}
//			//else
//			//{
//			//	cp[pid[j+1]].coor[0]+=cp[tmesh[eid].cnct[0]].coor[0]/16.;
//			//	cp[pid[j+1]].coor[1]+=cp[tmesh[eid].cnct[0]].coor[1]/16.;
//			//	cp[pid[j+1]].coor[2]+=cp[tmesh[eid].cnct[0]].coor[2]/16.;
//			//}
//			Edge edtmp3;
//			edtmp3.act=1;
//			edtmp3.pt[0]=pid[0]; edtmp3.pt[1]=pid[j+1];
//			edtmp3.len=tmedge[tmesh[eid].edge[(j+1)%4]].len/2.;
//			tmedge.push_back(edtmp3);
//			edid[3*j+2]=tmedge.size()-1;
//			int ied(tmedge[tmesh[eid].edge[j]].chd[0]);
//			if(tmedge[ied].pt[0]==itmp[0] || tmedge[ied].pt[1]==itmp[0])
//			{
//				edid[3*j]=ied;
//				edid[3*j+1]=tmedge[tmesh[eid].edge[j]].chd[1];
//			}
//			else
//			{
//				edid[3*j]=tmedge[tmesh[eid].edge[j]].chd[1];
//				edid[3*j+1]=ied;
//			}
//		}
//	}
//
//	//int nvl(cp[tmesh[eid].cnct[0]].face.size());
//	//if(cp[tmesh[eid].cnct[0]].update==0)
//	//{
//	//	double ccf[3]={1.-7./(4.*nvl),3./(2.*nvl*nvl),1./(4.*nvl*nvl)};
//	//	cp[tmesh[eid].cnct[0]].coortmp[0]=ccf[0]*cp[tmesh[eid].cnct[0]].coor[0]+ccf[1]*cp[tmesh[eid].cnct[1]].coor[0]+ccf[2]*cp[tmesh[eid].cnct[2]].coor[0];
//	//	cp[tmesh[eid].cnct[0]].coortmp[1]=ccf[0]*cp[tmesh[eid].cnct[0]].coor[1]+ccf[1]*cp[tmesh[eid].cnct[1]].coor[1]+ccf[2]*cp[tmesh[eid].cnct[2]].coor[1];
//	//	cp[tmesh[eid].cnct[0]].coortmp[2]=ccf[0]*cp[tmesh[eid].cnct[0]].coor[2]+ccf[1]*cp[tmesh[eid].cnct[1]].coor[2]+ccf[2]*cp[tmesh[eid].cnct[2]].coor[2];
//	//	cp[tmesh[eid].cnct[0]].update=2;
//	//}
//	//else if(cp[tmesh[eid].cnct[0]].update==2)
//	//{
//	//	double ccf[2]={3./(2.*nvl*nvl),1./(4.*nvl*nvl)};
//	//	cp[tmesh[eid].cnct[0]].coortmp[0]+=ccf[0]*cp[tmesh[eid].cnct[1]].coor[0]+ccf[1]*cp[tmesh[eid].cnct[2]].coor[0];
//	//	cp[tmesh[eid].cnct[0]].coortmp[1]+=ccf[0]*cp[tmesh[eid].cnct[1]].coor[1]+ccf[1]*cp[tmesh[eid].cnct[2]].coor[1];
//	//	cp[tmesh[eid].cnct[0]].coortmp[2]+=ccf[0]*cp[tmesh[eid].cnct[1]].coor[2]+ccf[1]*cp[tmesh[eid].cnct[2]].coor[2];
//	//}
//	//else
//	//{
//	//	cout<<"Not supported for other udpate types!\n";
//	//	getchar();
//	//}
//
//	//for(int j=0; j<4; j++)
//	//{
//	//	int nvl(cp[tmesh[eid].cnct[j]].face.size());
//	//	if(cp[tmesh[eid].cnct[j]].update==0)
//	//	{
//	//		cp[tmesh[eid].cnct[j]].update=1;
//	//		double ccf[3]={1.-7./(4.*nvl),3./(2.*nvl*nvl),1./(4.*nvl*nvl)};
//	//		cp[tmesh[eid].cnct[j]].coortmp[0]=ccf[0]*cp[tmesh[eid].cnct[j]].coor[0]+ccf[1]*cp[tmesh[eid].cnct[(j+1)%4]].coor[0]+ccf[2]*cp[tmesh[eid].cnct[(j+2)%4]].coor[0];
//	//		cp[tmesh[eid].cnct[j]].coortmp[1]=ccf[0]*cp[tmesh[eid].cnct[j]].coor[1]+ccf[1]*cp[tmesh[eid].cnct[(j+1)%4]].coor[1]+ccf[2]*cp[tmesh[eid].cnct[(j+2)%4]].coor[1];
//	//		cp[tmesh[eid].cnct[j]].coortmp[2]=ccf[0]*cp[tmesh[eid].cnct[j]].coor[2]+ccf[1]*cp[tmesh[eid].cnct[(j+1)%4]].coor[2]+ccf[2]*cp[tmesh[eid].cnct[(j+2)%4]].coor[2];
//	//	}
//	//	else
//	//	{
//	//		double ccf[2]={3./(2.*nvl*nvl),1./(4.*nvl*nvl)};
//	//		cp[tmesh[eid].cnct[j]].coortmp[0]+=ccf[0]*cp[tmesh[eid].cnct[(j+1)%4]].coor[0]+ccf[1]*cp[tmesh[eid].cnct[(j+2)%4]].coor[0];
//	//		cp[tmesh[eid].cnct[j]].coortmp[1]+=ccf[0]*cp[tmesh[eid].cnct[(j+1)%4]].coor[1]+ccf[1]*cp[tmesh[eid].cnct[(j+2)%4]].coor[1];
//	//		cp[tmesh[eid].cnct[j]].coortmp[2]+=ccf[0]*cp[tmesh[eid].cnct[(j+1)%4]].coor[2]+ccf[1]*cp[tmesh[eid].cnct[(j+2)%4]].coor[2];
//	//	}
//	//}
//
//	int e_cnct[4][4]={{tmesh[eid].cnct[0],pid[1],pid[0],pid[4]},{pid[1],tmesh[eid].cnct[1],pid[2],pid[0]},{pid[0],pid[2],tmesh[eid].cnct[2],pid[3]},{pid[4],pid[0],pid[3],tmesh[eid].cnct[3]}};
//	int e_edge[4][4]={{edid[0],edid[2],edid[11],edid[10]},{edid[1],edid[3],edid[5],edid[2]},{edid[5],edid[4],edid[6],edid[8]},{edid[11],edid[8],edid[7],edid[9]}};
//	int enewid[4];
//	int e_type[4]={4,0,0,0};
//	vector<Element> etmp(4);
//	for(int i=0; i<4; i++)
//	{
//		etmp[i].act=1;
//		etmp[i].type=e_type[i];
//		etmp[i].lv=tmesh[eid].lv+1.;
//		etmp[i].prt=eid;
//		for(int j=0; j<4; j++)
//		{
//			etmp[i].cnct[j]=e_cnct[i][j];
//			etmp[i].edge[j]=e_edge[i][j];
//		}
//		tmesh.push_back(etmp[i]);
//		enewid[i]=tmesh.size()-1;
//		tmesh[eid].chd[i]=enewid[i];
//	}
//
//	tmesh[eid].act=0;
//}
//
//void TruncatedTspline_2D::ElementRefine_Unstruct_Topo_Geom_2(int eid, int dir)
//{
//	int pid[2],edid[7],pos(dir);
//	int cnid[2]={pos,(pos+2)%4};
//	for(int j=0; j<2; j++)
//	{
//		int itmp[2]={tmesh[eid].cnct[cnid[j]],tmesh[eid].cnct[(cnid[j]+1)%4]};
//		if(tmedge[tmesh[eid].edge[cnid[j]]].act==1)
//		{
//			Vertex ptmp;
//			ptmp.coor[0]=(cp[itmp[0]].coor[0]+cp[itmp[1]].coor[0])/2.;
//			ptmp.coor[1]=(cp[itmp[0]].coor[1]+cp[itmp[1]].coor[1])/2.;
//			ptmp.coor[2]=(cp[itmp[0]].coor[2]+cp[itmp[1]].coor[2])/2.;
//			cp.push_back(ptmp);
//			pid[j]=cp.size()-1;
//			tmedge[tmesh[eid].edge[cnid[j]]].midpt=pid[j];
//			Edge edtmp1,edtmp2;
//			edtmp1.act=1;
//			edtmp1.pt[0]=itmp[0]; edtmp1.pt[1]=pid[j];
//			edtmp1.len=tmedge[tmesh[eid].edge[cnid[j]]].len/2.;
//			edtmp1.prt=tmesh[eid].edge[cnid[j]];
//			edtmp2.act=1;
//			edtmp2.pt[0]=pid[j]; edtmp2.pt[1]=itmp[1];
//			edtmp2.len=tmedge[tmesh[eid].edge[cnid[j]]].len/2.;
//			edtmp2.prt=tmesh[eid].edge[cnid[j]];
//			tmedge.push_back(edtmp1);
//			edid[2*j]=tmedge.size()-1;
//			tmedge[tmesh[eid].edge[cnid[j]]].chd[0]=edid[2*j];
//			tmedge.push_back(edtmp2);
//			edid[2*j+1]=tmedge.size()-1;
//			tmedge[tmesh[eid].edge[cnid[j]]].chd[1]=edid[2*j+1];
//			tmedge[tmesh[eid].edge[cnid[j]]].act=0;
//		}
//		else
//		{
//			pid[j]=tmedge[tmesh[eid].edge[cnid[j]]].midpt;
//			int ied(tmedge[tmesh[eid].edge[cnid[j]]].chd[0]);
//			if(tmedge[ied].pt[0]==itmp[0] || tmedge[ied].pt[1]==itmp[0])
//			{
//				edid[2*j]=ied;
//				edid[2*j+1]=tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
//			}
//			else
//			{
//				edid[2*j]=tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
//				edid[2*j+1]=ied;
//			}
//		}
//	}
//	Edge edtmp;
//	edtmp.act=1;
//	edtmp.pt[0]=pid[0]; edtmp.pt[1]=pid[1];
//	edtmp.len=tmedge[tmesh[eid].edge[(pos+1)%4]].len;
//	tmedge.push_back(edtmp);
//	edid[4]=tmedge.size()-1;
//	edid[5]=tmesh[eid].edge[(pos+1)%4];
//	edid[6]=tmesh[eid].edge[(pos+3)%4];
//
//	int e_cnct[2][4]={{tmesh[eid].cnct[pos],pid[0],pid[1],tmesh[eid].cnct[(pos+3)%4]},{pid[0],tmesh[eid].cnct[(pos+1)%4],tmesh[eid].cnct[(pos+2)%4],pid[1]}};
//	int e_edge[2][4]={{edid[0],edid[4],edid[3],edid[6]},{edid[1],edid[5],edid[2],edid[4]}};
//	int enewid[2];
//	double chd_org[2][2]={{0.,0.},{tmedge[tmesh[eid].edge[0]].len/2.,0.}};
//	if(dir==1)
//	{
//		chd_org[1][0]=0.; chd_org[1][1]=tmedge[tmesh[eid].edge[1]].len/2.;
//	}
//	vector<Element> etmp(2);
//	for(int i=0; i<2; i++)
//	{
//		etmp[i].act=1;
//		etmp[i].type=0;
//		etmp[i].lv=tmesh[eid].lv+0.5;
//		etmp[i].prt=eid;
//		for(int j=0; j<4; j++)
//		{
//			etmp[i].cnct[(pos+j)%4]=e_cnct[i][j];
//			etmp[i].edge[(pos+j)%4]=e_edge[i][j];
//		}
//		tmesh.push_back(etmp[i]);
//		enewid[i]=tmesh.size()-1;
//		tmesh[eid].chd[i]=enewid[i];
//		tmesh[eid].chd_o[i][0]=chd_org[i][0];
//		tmesh[eid].chd_o[i][1]=chd_org[i][1];
//	}
//
//	tmesh[eid].act=0;
//}



//void TruncatedTspline_2D::SurfRefine()
//{
//	CollectActives();
//	int rid_tmp1[18]={5,4,3,17,16,15,29,28,41,6,7,8,18,19,20,30,31,42};
//	int nsub1(18);
//	vector<int> rid1(nsub1);
//	for(int i=0; i<nsub1; i++)
//	{
//		rid1[i]=eaid[rid_tmp1[i]];
//	}
//	Refine_Target(rid1);
//
//	CollectActives();
//	int rid_tmp2[18]={51,50,55,52,53,56,63,62,64,86,87,90,89,88,93,98,99,101};
//	int nsub2(18);
//	vector<int> rid2(nsub2);
//	for(int i=0; i<nsub2; i++)
//	{
//		rid2[i]=eaid[rid_tmp2[i]];
//	}
//	Refine_Addition(rid2);
//	Refine_Target(rid2);
//
//	CollectActives();
//	int rid_tmp3[18]={109,108,113,110,111,114,121,120,122,144,145,148,147,146,151,156,157,159};
//	int nsub3(18);
//	vector<int> rid3(nsub3);
//	for(int i=0; i<nsub3; i++)
//	{
//		rid3[i]=eaid[rid_tmp3[i]];
//	}
//	Refine_Addition(rid3);
//	Refine_Target(rid3);
//
//	CollectActives();
//	int rid_tmp4[18]={167,166,171,168,169,172,179,178,180,202,203,206,205,204,209,214,215,217};
//	int nsub4(18);
//	vector<int> rid4(nsub4);
//	for(int i=0; i<nsub4; i++)
//	{
//		rid4[i]=eaid[rid_tmp4[i]];
//	}
//	Refine_Addition(rid4);
//	Refine_Target(rid4);
//
//	CollectActives();
//	int rid_tmp5[20]={225,224,229,226,227,230,237,236,238,260,261,264,263,262,267,272,273,275,38,49};
//	int nsub5(20);
//	vector<int> rid5(nsub5);
//	for(int i=0; i<nsub5; i++)
//	{
//		rid5[i]=eaid[rid_tmp5[i]];
//	}
//	Refine_Addition(rid5);
//	Refine_Target(rid5);
//
//	CollectActives();
//}
//
//void TruncatedTspline_2D::LshapeRefine()
//{
//	int rid_tmp1[26]={20,21,22,36,37,38,52,53,54,68,69,70,97,
//		25,26,27,41,42,43,57,58,59,73,74,75,110};
//	vector<int> rid1(rid_tmp1,rid_tmp1+26);
//	Refine_Addition(rid1);
//	Refine_Target(rid1);
//
//	int rid_tmp2[22]={133,136,137,134,139,138,145,148,149,151,150,
//		180,181,184,183,182,187,192,193,196,195,194};
//	vector<int> rid2(rid_tmp2,rid_tmp2+22);
//	Refine_Addition(rid2);
//	Refine_Target(rid2);
//
//	int rid_tmp3[14]={273,276,277,279,278,288,289,
//		312,313,316,315,314,324,325};
//	vector<int> rid3(rid_tmp3,rid_tmp3+14);
//	Refine_Addition(rid3);
//	Refine_Target(rid3);
//
//	int rid_tmp4[10]={396,397,399,398,405,416,417,419,418,428};
//	vector<int> rid4(rid_tmp4,rid_tmp4+10);
//	Refine_Addition(rid4);
//	Refine_Target(rid4);
//
//	CollectActives();
//
//	//cout<<tmesh[123].type<<'\n';
//	//int ntjc=0;
//	//for(int i=0; i<4; i++)
//	//{
//	//	if(tmedge[tmesh[123].edge[i]].act==0) 
//	//	{
//	//		cout<<i<<'\n';
//	//		ntjc++;
//	//	}
//	//}
//	//cout<<ntjc<<'\n';
//	//getchar();
//}
//
//void TruncatedTspline_2D::LshapeRefine_H1()
//{
//	CollectActives();
//	int rid_tmp1[18]={5,4,3,17,16,15,29,28,41,6,7,8,18,19,20,30,31,42};
//	int nsub1(18);
//	vector<int> rid1(nsub1);
//	for(int i=0; i<nsub1; i++)
//	{
//		rid1[i]=eaid[rid_tmp1[i]];
//	}
//	Refine_Addition(rid1);
//	Refine_Target(rid1);
//
//	CollectActives();
//	int rid_tmp2[18]={51,50,55,52,53,56,63,62,64,86,87,90,89,88,93,98,99,101};
//	int nsub2(18);
//	vector<int> rid2(nsub2);
//	for(int i=0; i<nsub2; i++)
//	{
//		rid2[i]=eaid[rid_tmp2[i]];
//	}
//	Refine_Addition(rid2);
//	Refine_Target(rid2);
//
//	CollectActives();
//	int rid_tmp3[18]={109,108,113,110,111,114,121,120,122,144,145,148,147,146,151,156,157,159};
//	int nsub3(18);
//	vector<int> rid3(nsub3);
//	for(int i=0; i<nsub3; i++)
//	{
//		rid3[i]=eaid[rid_tmp3[i]];
//	}
//	Refine_Addition(rid3);
//	Refine_Target(rid3);
//
//	CollectActives();
//	int rid_tmp4[18]={167,166,171,168,169,172,179,178,180,202,203,206,205,204,209,214,215,217};
//	int nsub4(18);
//	vector<int> rid4(nsub4);
//	for(int i=0; i<nsub4; i++)
//	{
//		rid4[i]=eaid[rid_tmp4[i]];
//	}
//	Refine_Addition(rid4);
//	Refine_Target(rid4);
//
//	CollectActives();
//	int rid_tmp5[20]={225,224,229,226,227,230,237,236,238,260,261,264,263,262,267,272,273,275,38,49};
//	int nsub5(20);
//	vector<int> rid5(nsub5);
//	for(int i=0; i<nsub5; i++)
//	{
//		rid5[i]=eaid[rid_tmp5[i]];
//	}
//	Refine_Addition(rid5);
//	Refine_Target(rid5);
//
//	CollectActives();
//}
//
//void TruncatedTspline_2D::LshapeRefine_XP()
//{
//	int rid_tmp1[42]={7,2,3,8,13,12,11,14,15,42,43,41,40,38,39,16,19,18,17,23,20,
//		106,107,110,105,104,109,102,103,98,142,143,141,140,130,131,117,118,116,119,113,114};
//	vector<int> rid1(rid_tmp1,rid_tmp1+42), rid1_more;
//	vector<int> rftype1(42,0), rftype1_more;
//	for(int i=0; i<42; i++)
//	{
//		if(tmesh[rid_tmp1[i]].type==4)
//		{
//			rftype1[i]=4;
//		}
//	}
//	Topo_Refine_Unstruct(rid1,rftype1,rid1_more,rftype1_more);
//	Geom_Refine_Unstruct(rid1_more,rftype1_more);
//
//	CollectActives();
//	int rid_tmp2[22]={79,78,83,80,81,84,91,90,95,92,93,194,195,190,197,196,193,182,183,178,185,184};//active element id, need to convert
//	int rid_tmp2_1[2]={96,181};
//	int rid_tmp2_dir[2]={238,311};
//	int nsub2(22),nsub2_1(2);
//	for(int i=0; i<nsub2_1; i++)
//	{
//		int* it=find(tmesh[eaid[rid_tmp2_1[i]]].cnct,tmesh[eaid[rid_tmp2_1[i]]].cnct+4,rid_tmp2_dir[i]);
//		int loc(it-tmesh[eaid[rid_tmp2_1[i]]].cnct);
//		if(loc==0 || loc==2) rid_tmp2_dir[i]=0;
//		else if(loc==1 || loc==3) rid_tmp2_dir[i]=1;
//		else
//		{
//			cout<<"Cannot find refinement direction!\n";
//			getchar();
//		}
//	}
//	vector<int> rid2(nsub2), rid2_more;
//	vector<int> rftype2(nsub2,0), rftype2_more;
//	for(int i=0; i<nsub2; i++)
//	{
//		rid2[i]=eaid[rid_tmp2[i]];
//		if(tmesh[eaid[rid_tmp2[i]]].type==4)
//		{
//			rftype2[i]=4;
//		}
//	}
//	for(int i=0; i<nsub2_1; i++)
//	{
//		rid2.push_back(eaid[rid_tmp2_1[i]]);
//		rftype2.push_back(rid_tmp2_dir[i]+1);
//	}
//	Topo_Refine_Unstruct(rid2,rftype2,rid2_more,rftype2_more);
//	Geom_Refine_Unstruct(rid2_more,rftype2_more);
//
//	CollectActives();
//	int rid_tmp3[22]={223,222,227,224,225,228,235,234,239,236,237,266,267,270,269,268,273,278,279,282,281,280};//active element id, need to convert
//	int rid_tmp3_1[2]={240,285};
//	int rid_tmp3_dir[2]={395,432};
//	int nsub3(22),nsub3_1(2);
//	for(int i=0; i<nsub3_1; i++)
//	{
//		int* it=find(tmesh[eaid[rid_tmp3_1[i]]].cnct,tmesh[eaid[rid_tmp3_1[i]]].cnct+4,rid_tmp3_dir[i]);
//		int loc(it-tmesh[eaid[rid_tmp3_1[i]]].cnct);
//		if(loc==0 || loc==2) rid_tmp3_dir[i]=0;
//		else if(loc==1 || loc==3) rid_tmp3_dir[i]=1;
//		else
//		{
//			cout<<"Cannot find refinement direction!\n";
//			getchar();
//		}
//	}
//	vector<int> rid3(nsub3), rid3_more;
//	vector<int> rftype3(nsub3,0), rftype3_more;
//	for(int i=0; i<nsub3; i++)
//	{
//		rid3[i]=eaid[rid_tmp3[i]];
//		if(tmesh[eaid[rid_tmp3[i]]].type==4)
//		{
//			rftype3[i]=4;
//		}
//	}
//	for(int i=0; i<nsub3_1; i++)
//	{
//		rid3.push_back(eaid[rid_tmp3_1[i]]);
//		rftype3.push_back(rid_tmp3_dir[i]+1);
//	}
//	Topo_Refine_Unstruct(rid3,rftype3,rid3_more,rftype3_more);
//	Geom_Refine_Unstruct(rid3_more,rftype3_more);
//
//	CollectActives();
//	int rid_tmp4[16]={291,290,295,292,293,296,303,302,334,335,338,337,336,341,346,347};//active element id, need to convert
//	int rid_tmp4_1[2]={307,350};
//	int rid_tmp4_dir[2]={378,524};
//	int nsub4(16),nsub4_1(2);
//	for(int i=0; i<nsub4_1; i++)
//	{
//		int* it=find(tmesh[eaid[rid_tmp4_1[i]]].cnct,tmesh[eaid[rid_tmp4_1[i]]].cnct+4,rid_tmp4_dir[i]);
//		int loc(it-tmesh[eaid[rid_tmp4_1[i]]].cnct);
//		if(loc==0 || loc==2) rid_tmp4_dir[i]=0;
//		else if(loc==1 || loc==3) rid_tmp4_dir[i]=1;
//		else
//		{
//			cout<<"Cannot find refinement direction!\n";
//			getchar();
//		}
//	}
//	vector<int> rid4(nsub4), rid4_more;
//	vector<int> rftype4(nsub4,0), rftype4_more;
//	for(int i=0; i<nsub4; i++)
//	{
//		rid4[i]=eaid[rid_tmp4[i]];
//		if(tmesh[eaid[rid_tmp4[i]]].type==4)
//		{
//			rftype4[i]=4;
//		}
//	}
//	for(int i=0; i<nsub4_1; i++)
//	{
//		rid4.push_back(eaid[rid_tmp4_1[i]]);
//		rftype4.push_back(rid_tmp4_dir[i]+1);
//	}
//	Topo_Refine_Unstruct(rid4,rftype4,rid4_more,rftype4_more);
//	Geom_Refine_Unstruct(rid4_more,rftype4_more);
//
//	//CollectActives();
//}