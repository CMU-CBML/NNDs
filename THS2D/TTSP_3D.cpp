#include "TTSP_3D.h"
#include "KnotInsertion.h"
#include "BSplineBasis.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>
//#include "Matlab_Solver_wap.h"
#include "SingularEval.h"
#include <sstream>
#include <iomanip>
#include "LeastSquare.h"

typedef unsigned int uint;

int fc_cnct[4][4] = { { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
int ed_cnct[4][4] = { { 0, 5, 8, 4 }, { 1, 6, 9, 5 }, { 2, 7, 10, 6 }, { 3, 4, 11, 7 } };
int ed_uvw[8][3] = { { 0, 3, 4 }, { 1, 0, 5 }, { 2, 1, 6 }, { 3, 2, 7 }, { 11, 8, 4 }, { 8, 9, 5 }, { 9, 10, 6 }, { 10, 11, 7 } };
int fc_ppd_ed[6] = { 4, 3, 0, 3, 0, 4 };
int fc_opst[6] = { 5, 3, 4, 1, 2, 0 };

int solid_fc[6][4] = { { 0, 1, 2, 3 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3,0, 4, 7 }, {4,5,6,7} };
int solid_ed[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 } };

/*
1) InitialRotate() - set rotation matrix for each hex w.r.t. each of its non-boundary faces, note that this rotation is element-element rotation, that is, the transformation
   between local coordinate systems of two neighboring elements; so for each shared reference vertex, we need to fix its LCS, note this LCS is different from
   that of element, it is the LCS of vertex!
2) setReference() - for each vertex, set a LCS, note that this is generally different from what we mentioned in InitialRotate(); each vertex can always be located as
   a corner of a certain hex, then we set...
*/

TruncatedTspline_3D::TruncatedTspline_3D()
{
	//cp.clear();
	//tmesh.clear();
}

void TruncatedTspline_3D::CreateUniformCube(string fn)
{
	int nsmp(5);
	int npt = (nsmp + 1)*(nsmp + 1)*(nsmp + 1);
	int nel = nsmp*nsmp*nsmp;
	vector<array<double, 3>> spt(npt);
	vector<array<int, 8>> sele(nel);
	int i, j, k, loc(0);
	for (i = 0; i < nsmp + 1; i++)
	{
		for (j = 0; j < nsmp + 1; j++)
		{
			for (k = 0; k < nsmp + 1; k++)
			{
				spt[loc][0] = k;
				spt[loc][1] = j;
				spt[loc][2] = i;
				loc++;
			}
		}
	}
	loc = 0;
	for (i = 0; i < nsmp; i++)
	{
		for (j = 0; j < nsmp; j++)
		{
			for (k = 0; k < nsmp; k++)
			{
				sele[loc][0] = i*(nsmp + 1)*(nsmp + 1) + j*(nsmp + 1) + k;
				sele[loc][1] = i*(nsmp + 1)*(nsmp + 1) + j*(nsmp + 1) + k+1;
				sele[loc][2] = i*(nsmp + 1)*(nsmp + 1) + (j+1)*(nsmp + 1) + k + 1;
				sele[loc][3] = i*(nsmp + 1)*(nsmp + 1) + (j + 1)*(nsmp + 1) + k;
				sele[loc][4] = (i+1)*(nsmp + 1)*(nsmp + 1) + j*(nsmp + 1) + k;
				sele[loc][5] = (i + 1)*(nsmp + 1)*(nsmp + 1) + j*(nsmp + 1) + k + 1;
				sele[loc][6] = (i + 1)*(nsmp + 1)*(nsmp + 1) + (j + 1)*(nsmp + 1) + k + 1;
				sele[loc][7] = (i + 1)*(nsmp + 1)*(nsmp + 1) + (j + 1)*(nsmp + 1) + k;
				loc++;
			}
		}
	}

	string fname = fn + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nCube hex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (uint i = 0; i<sele.size(); i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i<sele.size(); i++)
		{
			fout << "12\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

double TruncatedTspline_3D::PartionOfUnity(int eid, const array<double, 3>& u)
{
	double sum(0.);
	vector<double> Nt;
	vector<array<double, 3>> dNdt;
	ElementBasis(eid,u,Nt,dNdt);
	for (uint i = 0; i < Nt.size(); i++) sum += Nt[i];
	return sum;
}

void TruncatedTspline_3D::VisualizeControlMesh(string fn)
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
			if (tmesh[i].act == 1 /*&& tmesh[i].type == 2*/ /*&& tmesh[i].type != 2 && tmesh[i].type != 3*/) nel_act++;
		}
		fout<<"\nCELLS "<<nel_act<<" "<<9*nel_act<<'\n';
		for(uint i=0; i<tmesh.size(); i++)
		{
			if (tmesh[i].act == 1 /*&& tmesh[i].type == 2*/ /*&& tmesh[i].type != 2 && tmesh[i].type != 3*/)
			{
				fout<<"8 ";
				for(int j=0; j<8; j++)
				{
					fout<<tmesh[i].cnct[j]<<' ';
				}
				fout<<'\n';
			}
		}
		fout<<"\nCELL_TYPES "<<nel_act<<'\n';
		for(uint i=0; i<nel_act; i++)
		{
			fout<<"12\n";
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
		fout<<"POINT_DATA "<<cp.size()<<"\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		for(uint i=0;i<cp.size();i++)
		{
			fout<<cp[i].trun<<"\n";
		}

		fout.close();
	}
	else
	{
		cout<<"Cannot open "<<fname<<"!\n";
	}
}

void TruncatedTspline_3D::CollectActives()
{
	haid.clear();
	faid.clear();
	eaid.clear();
	int loc(0);
	for(uint i=0; i<tmesh.size(); i++)
	{
		tmesh[i].id_act = -1;
		if(tmesh[i].act==1 /*&& tmesh[i].type!=2 && tmesh[i].type!=3*/)
		{
			haid.push_back(i);
			tmesh[i].id_act = loc++;
		}
	}
	loc = 0;
	for (uint i = 0; i<tmface.size(); i++)
	{
		tmface[i].id_act = -1;
		if (tmface[i].act == 1)
		{
			faid.push_back(i);
			tmface[i].id_act = loc++;
		}
	}
	loc = 0;
	for (uint i = 0; i<tmedge.size(); i++)
	{
		tmedge[i].id_act = -1;
		if (tmedge[i].act == 1)
		{
			eaid.push_back(i);
			tmedge[i].id_act = loc++;
		}
	}
}

//void TruncatedTspline_3D::StrongBalanceCheck(const vector<int>& rid, vector<int>& rid2)//later
//{
//	uint i,j,k;
//	rid2.clear();
//	for(i=0; i<rid.size(); i++)
//	{
//		if(tmesh[rid[i]].act==1)
//		{
//		}
//	}
//}
//
//void TruncatedTspline_3D::FaceIntersectCheck(vector<int>& rid2)//later
//{
//	rid2.clear();
//	for(uint i=0; i<tmesh.size(); i++)
//	{
//		//tmesh[i].ref=0;
//		if(tmesh[i].act==1)
//		{
//			int ntjc(0);
//			for(int j=0; j<4; j++)
//			{
//				if(tmedge[tmesh[i].edge[j]].act==0)
//				{
//					ntjc++;
//				}
//			}
//			if(ntjc==1 && tmesh[i].type==2)//boundary element
//			{
//				//tmesh[i].ref=10;
//				rid2.push_back(i);
//			}
//			else if(ntjc==2)
//			{
//				int pos(0);
//				for(int j=0; j<4; j++)
//				{
//					if(tmedge[tmesh[i].edge[j]].act==0)
//					{
//						pos=j;
//						break;
//					}
//				}
//				if(tmedge[tmesh[i].edge[(pos+1)%4]].act==0)
//				{
//					//tmesh[i].ref=20;
//				}
//				else if(tmedge[tmesh[i].edge[(pos+2)%4]].act==0)
//				{
//					//tmesh[i].ref=21;
//				}
//				rid2.push_back(i);
//			}
//			else if(ntjc==3)
//			{
//				//tmesh[i].ref=3;
//				rid2.push_back(i);
//			}
//			else if(ntjc==4)
//			{
//				//tmesh[i].ref=4;
//				rid2.push_back(i);
//			}
//		}
//	}
//}
//
//void TruncatedTspline_3D::StrongBalanceRefine(const vector<int>& ridsb)
//{
//	for(uint i=0; i<ridsb.size(); i++)
//	{
//		if(tmesh[ridsb[i]].act==1)
//		{
//			//if(tmesh[ridsb[i]].type==0)
//			//{
//			//	ElementRefine_Square_4(ridsb[i]);
//			//}
//			//else if(tmesh[ridsb[i]].type==1)
//			//{
//			//	ElementRefine_Rectangular(ridsb[i]);
//			//}
//			//else if(tmesh[ridsb[i]].type==2)
//			//{
//			//	ElementRefine_Boundary(ridsb[i]);
//			//}
//		}
//	}
//}
//
//void TruncatedTspline_3D::TargetRefine(const vector<int>& rid)
//{
//	for(uint i=0; i<rid.size(); i++)
//	{
//		if(tmesh[rid[i]].act==1)
//		{
//			//if(tmesh[rid[i]].type==0)
//			//{
//			//	ElementRefine_Square_4(rid[i]);
//			//}
//			//else if(tmesh[rid[i]].type==1)
//			//{
//			//	ElementRefine_Rectangular(rid[i]);
//			//}
//		}
//	}
//}
//
//void TruncatedTspline_3D::OneTjunctionRefine(const vector<int>& ridtp)//later
//{
//	
//}
//
//void TruncatedTspline_3D::TjuncExtentCheck()
//{
//	for(uint i=0; i<tmesh.size(); i++)
//	{
//		if(tmesh[i].act==1 && tmesh[i].type==0)
//		{
//			int pos(-1);
//			for(int j=0; j<4; j++)
//			{
//				if(tmedge[tmesh[i].edge[j]].act==0)
//				{
//					pos=j;
//					break;
//				}
//			}
//			if(pos!=-1)
//			{
//				int ed[2]={(pos+1)%4,(pos+3)%4},ref(0);
//				for(int k=0; k<2; k++)
//				{
//					int ednb(tmedge[tmesh[i].edge[ed[k]]].face[0]);
//					if(ednb==i) ednb=tmedge[tmesh[i].edge[ed[k]]].face[1];
//					int* it=find(tmesh[ednb].edge,tmesh[ednb].edge+4,ed[k]);
//					int loc=it-tmesh[ednb].edge;
//					loc=(loc+2)%4;
//					if(tmedge[tmesh[ednb].edge[loc]].act==0)
//						ref=1;
//				}
//				if(ref==1)
//				{
//					//ElementRefine_Square_2(i);
//				}
//			}
//		}
//	}
//}
//
//void TruncatedTspline_3D::TjuncExtentCheck_1(vector<int>& ridtjx)
//{
//	ridtjx.clear();
//	for(uint i=0; i<tmesh.size(); i++)
//	{
//		if(tmesh[i].act==1)
//		{
//			int pos(-1);
//			for(int j=0; j<4; j++)
//			{
//				if(tmedge[tmesh[i].edge[j]].act==0)
//				{
//					pos=j;
//					break;
//				}
//			}
//			if(pos!=-1)
//			{
//				int ed[2]={(pos+1)%4,(pos+3)%4},ref(0);
//				for(int k=0; k<2; k++)
//				{
//					int ednb(tmedge[tmesh[i].edge[ed[k]]].face[0]);
//					if(ednb==i) ednb=tmedge[tmesh[i].edge[ed[k]]].face[1];
//					int* it=find(tmesh[ednb].edge,tmesh[ednb].edge+4,tmesh[i].edge[ed[k]]);
//					if(it==tmesh[ednb].edge+4)
//					{
//						vector<int>::iterator it1=find(ridtjx.begin(),ridtjx.end(),ednb);
//						if(it1==ridtjx.end())
//							ridtjx.push_back(ednb);
//					}
//				}
//			}
//		}
//	}
//}
//
//void TruncatedTspline_3D::TjuncExtentCheck_2(vector<int>& ridtjx)
//{
//	ridtjx.clear();
//	for(uint i=0; i<tmesh.size(); i++)
//	{
//		if(tmesh[i].act==1 && tmesh[i].type==0)
//		{
//			int pos(-1);
//			for(int j=0; j<4; j++)
//			{
//				if(tmedge[tmesh[i].edge[j]].act==0)
//				{
//					pos=j;
//					break;
//				}
//			}
//			if(pos!=-1)
//			{
//				int ed[2]={(pos+1)%4,(pos+3)%4},ref(0);
//				for(int k=0; k<2; k++)
//				{
//					int ednb(tmedge[tmesh[i].edge[ed[k]]].face[0]);
//					if(ednb==i) ednb=tmedge[tmesh[i].edge[ed[k]]].face[1];
//					int* it=find(tmesh[ednb].edge,tmesh[ednb].edge+4,tmesh[i].edge[ed[k]]);
//					int loc=it-tmesh[ednb].edge;
//					loc=(loc+2)%4;
//					if(tmedge[tmesh[ednb].edge[loc]].act==0)
//					{
//						vector<int>::iterator it=find(ridtjx.begin(),ridtjx.end(),i);
//						if(it==ridtjx.end())
//							ridtjx.push_back(i);
//					}
//				}
//			}
//		}
//	}
//}

void TruncatedTspline_3D::VisualizeTMesh(string fn)
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
				//fout << tmedge[i].type << "\n";
		}
		//fout<<"POINT_DATA "<<cp.size()<<"\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<cp.size();i++)
		//{
		//	fout<<cp[i].trun<<"\n";
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

void TruncatedTspline_3D::VisualizeFaceMesh(string fn)
{
	string fname(fn + ".vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << cp.size() << " float\n";
		for (uint i = 0; i<cp.size(); i++)
		{
			fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
		}
		int nfc(0);
		for (uint i = 0; i<tmface.size(); i++)
		{
			if (tmface[i].act == 1 /*&& tmface[i].type==1*/) nfc++;
		}
		fout << "\nCELLS " << nfc << " " << 5 * nfc << '\n';
		for (uint i = 0; i<tmface.size(); i++)
		{
			if (tmface[i].act == 1 /*&& tmface[i].type == 1*/)
				fout << "4 " << tmface[i].cnct[0] << " " << tmface[i].cnct[1] << " " << tmface[i].cnct[2] << " " << tmface[i].cnct[3] << '\n';
		}
		fout << "\nCELL_TYPES " << nfc << '\n';
		for (uint i = 0; i<tmface.size(); i++)
		{
			if (tmface[i].act == 1 /*&& tmface[i].type == 1*/)
				fout << "9\n";
		}
		//fout << "\nCELL_DATA " << nfc << "\nSCALARS len float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<tmface.size(); i++)
		//{
		//	if (tmface[i].act == 1)
		//		fout << tmedge[i].len << "\n";
		//}
		//fout<<"POINT_DATA "<<cp.size()<<"\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<cp.size();i++)
		//{
		//	fout<<cp[i].trun<<"\n";
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
		cout << "Cannot open " << fname << "!\n";
	}
}

//void TruncatedTspline_3D::ShootRay_Edge(int edid, int pid, double kv[4])
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
//void TruncatedTspline_3D::ShootRay_Face(int fcid, int pid, double kv[4])
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

bool TruncatedTspline_3D::CheckSubKnotVector(const double ku1[5], const double kv1[5], const double ku2[5], const double kv2[5])
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

bool TruncatedTspline_3D::CheckSubKnotVector(const array<double,5>& ku1, const array<double,5>& kv1, const array<double,5>& ku2, const array<double,5>& kv2)
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

void TruncatedTspline_3D::InitialConnect()
{
	uint i,j;
	tmedge.clear();
	tmface.clear();
	//construct edges
	for(i=0; i<tmesh.size(); i++)
	{
		for(j=0; j<4; j++)
		{
			Edge3D edtmp;
			edtmp.pt[0]=tmesh[i].cnct[j];
			edtmp.pt[1]=tmesh[i].cnct[(j+1)%4];
			vector<Edge3D>::iterator it=find(tmedge.begin(),tmedge.end(),edtmp);
			int edid(it-tmedge.begin());
			if(it==tmedge.end())
			{
				tmedge.push_back(edtmp);
			}
			tmesh[i].edge[j]=edid;
			//tmedge[edid].hex.push_back(i);
		}
		for(j=0; j<4; j++)
		{
			Edge3D edtmp;
			edtmp.pt[0]=tmesh[i].cnct[j];
			edtmp.pt[1]=tmesh[i].cnct[j+4];
			vector<Edge3D>::iterator it=find(tmedge.begin(),tmedge.end(),edtmp);
			int edid(it-tmedge.begin());
			if(it==tmedge.end())
			{
				tmedge.push_back(edtmp);
			}
			tmesh[i].edge[j+4]=edid;
			//tmedge[edid].hex.push_back(i);
		}
		for(j=0; j<4; j++)
		{
			Edge3D edtmp;
			edtmp.pt[0]=tmesh[i].cnct[j+4];
			edtmp.pt[1]=tmesh[i].cnct[(j+1)%4+4];
			vector<Edge3D>::iterator it=find(tmedge.begin(),tmedge.end(),edtmp);
			int edid(it-tmedge.begin());
			if(it==tmedge.end())
			{
				tmedge.push_back(edtmp);
			}
			tmesh[i].edge[j+8]=edid;
			//tmedge[edid].hex.push_back(i);
		}
	}
	//construct faces
	for(i=0; i<tmesh.size(); i++)
	{
		//one bottom face
		Face3D fc1;
		for(j=0; j<4; j++)
		{
			fc1.cnct[j]=tmesh[i].cnct[j];
			fc1.edge[j]=tmesh[i].edge[j];
		}
		vector<Face3D>::iterator it1=find(tmface.begin(),tmface.end(),fc1);
		int fc1id(it1-tmface.begin());
		if(it1==tmface.end())
		{
			tmface.push_back(fc1);
		}
		tmesh[i].face[0]=fc1id;
		//tmface[fc1id].hex.push_back(i);
		//4 side faces
		for(j=0; j<4; j++)
		{
			Face3D fc;
			for(int k=0; k<4; k++)
			{
				fc.cnct[k]=tmesh[i].cnct[fc_cnct[j][k]];
				fc.edge[k]=tmesh[i].edge[ed_cnct[j][k]];
			}
			vector<Face3D>::iterator it=find(tmface.begin(),tmface.end(),fc);
			int fcid(it-tmface.begin());
			if(it==tmface.end())
			{
				tmface.push_back(fc);
			}
			tmesh[i].face[j+1]=fcid;
			//tmface[fcid].hex.push_back(i);
		}
		//one top face
		Face3D fc2;
		for(j=0; j<4; j++)
		{
			fc2.cnct[j]=tmesh[i].cnct[j+4];
			fc2.edge[j]=tmesh[i].edge[j+8];
		}
		vector<Face3D>::iterator it2=find(tmface.begin(),tmface.end(),fc2);
		int fc2id(it2-tmface.begin());
		if(it2==tmface.end())
		{
			tmface.push_back(fc2);
		}
		tmesh[i].face[5]=fc2id;
		//tmface[fc2id].hex.push_back(i);
	}
	//vertex-to-hex, edge-to-hex, face-to-hex
	for(i=0; i<tmesh.size(); i++)
	{
		for(j=0; j<8; j++)
		{
			cp[tmesh[i].cnct[j]].hex.push_back(i);
		}
		for(j=0; j<12; j++)
		{
			tmedge[tmesh[i].edge[j]].hex.push_back(i);
		}
		for(j=0; j<6; j++)
		{
			tmface[tmesh[i].face[j]].hex.push_back(i);
		}
	}
	//vertex-to-face, edge-to-face
	for(i=0; i<tmface.size(); i++)
	{
		for(j=0; j<4; j++)
		{
			cp[tmface[i].cnct[j]].face.push_back(i);
			tmedge[tmface[i].edge[j]].face.push_back(i);
		}
	}
	//vertex-to-edge
	for(i=0; i<tmedge.size(); i++)
	{
		for(j=0; j<2; j++)
		{
			cp[tmedge[i].pt[j]].edge.push_back(i);
		}
	}
	//find BC face, edge, vertex
	int ed0[6][4] = { { 4, 5, 6, 7 }, { 1, 3, 9, 11 }, { 0, 2, 8, 10 }, { 1, 3, 9, 11 }, { 0, 2, 8, 10 }, { 4, 5, 6, 7 } };//order could be wrong, but doesn't matter
	for(i=0; i<tmface.size(); i++)
	{
		if(tmface[i].hex.size()==1)
		{
			tmface[i].type = 1;
			tmesh[tmface[i].hex[0]].type=1;
			for(j=0; j<4; j++)
			{
				cp[tmface[i].cnct[j]].type=1;
				tmedge[tmface[i].edge[j]].type=1;
			}
			//set zero length edges
			//int hexid(tmface[i].hex[0]);
			//int* it=find(tmesh[hexid].face,tmesh[hexid].face+6,i);
			//int fc_loc(it-tmesh[hexid].face);
			//for(j=0; j<4; j++)
			//{
			//	//tmedge[tmesh[hexid].edge[ed0[fc_loc][j]]].len=0.;
			//}
		}
	}
	//find extraordinary edges and vertices
	for(i=0; i<tmedge.size(); i++)
	{
		if(tmedge[i].type!=1 && tmedge[i].hex.size()!=4)
		{
			tmedge[i].type=2;
			if(cp[tmedge[i].pt[0]].type!=1)
				cp[tmedge[i].pt[0]].type=3;
			if(cp[tmedge[i].pt[1]].type!=1)
				cp[tmedge[i].pt[1]].type=3;
		}
	}
	//boundary
	for(i=0; i<cp.size(); i++)
	{
		if(cp[i].type==1)
		{
			//int val(0);
			//for(j=0; j<cp[i].face.size(); j++)
			//{
			//	if(tmface[cp[i].face[j]].type==1) val++;
			//}
			//if(val==3 || val>4) cp[i].type=13;
		}
	}
	//find irregular elements
	for(i=0; i<tmesh.size(); i++)
	{
		if(tmesh[i].type!=1)
		{
			for(j=0; j<12; j++)
			{
				if(tmedge[tmesh[i].edge[j]].type==2)
				{
					tmesh[i].type=2;
					break;
				}
			}
		}
	}
	//additional boundary elements
	for (i = 0; i<tmesh.size(); i++)
	{
		if (tmesh[i].type != 1)
		{
			for (j = 0; j<8; j++)
			{
				if (cp[tmesh[i].cnct[j]].type == 1)
				{
					tmesh[i].type = 1;
					break;
				}
			}
		}
	}
	//boundry extraordinary points
	for (i = 0; i<cp.size(); i++)
	{
		if (cp[i].type == 1)
		{
			int count(0);
			for (j = 0; j < cp[i].edge.size(); j++)
			{
				if (tmedge[cp[i].edge[j]].type == 2) count++;
			}
			if (count==1) cp[i].bcxp = 1;
			else if (count>1) cp[i].bcxp = 2;
		}
	}
}

void TruncatedTspline_3D::UpdateConnect()
{
	uint i, j, k;
	//initialization
	for (i = 0; i < cp.size(); i++)
	{
		cp[i].hex.clear();
		cp[i].face.clear();
		cp[i].edge.clear();
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		tmedge[i].face.clear();
		tmedge[i].hex.clear();
	}
	for (i = 0; i < tmface.size(); i++)
	{
		tmface[i].hex.clear();
	}
	//update
	for (i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			for (j = 0; j < 8; j++)
			{
				cp[tmesh[i].cnct[j]].hex.push_back(i);
			}
			for (j = 0; j < 12; j++)
			{
				if (tmedge[tmesh[i].edge[j]].act == 1)
				{
					tmedge[tmesh[i].edge[j]].hex.push_back(i);
				}
				else
				{
					cp[tmedge[tmesh[i].edge[j]].midpt].hex.push_back(i);
					int edchd[2] = { tmedge[tmesh[i].edge[j]].chd[0], tmedge[tmesh[i].edge[j]].chd[1] };
					tmedge[edchd[0]].hex.push_back(i);
					tmedge[edchd[1]].hex.push_back(i);
				}
			}
			for (j = 0; j < 6; j++)
			{
				if (tmface[tmesh[i].face[j]].act == 1)
				{
					tmface[tmesh[i].face[j]].hex.push_back(i);
				}
				else
				{
					if (tmface[tmesh[i].face[j]].ctpt != -1)
					{
						cp[tmface[tmesh[i].face[j]].ctpt].hex.push_back(i);
					}
					for (k = 0; k < tmface[tmesh[i].face[j]].Tedge.size(); k++)
					{
						if (tmedge[tmface[tmesh[i].face[j]].Tedge[k]].act==1)
							tmedge[tmface[tmesh[i].face[j]].Tedge[k]].hex.push_back(i);
						else
						{
							tmedge[tmedge[tmface[tmesh[i].face[j]].Tedge[k]].chd[0]].hex.push_back(i);
							tmedge[tmedge[tmface[tmesh[i].face[j]].Tedge[k]].chd[1]].hex.push_back(i);
						}
					}
					for (k = 0; k < tmface[tmesh[i].face[j]].chd.size(); k++)
					{
						tmface[tmface[tmesh[i].face[j]].chd[k]].hex.push_back(i);
					}
				}
			}
		}
	}
	for (i = 0; i < tmface.size(); i++)
	{
		if (tmface[i].act == 1)
		{
			for (j = 0; j < 4; j++)
			{
				cp[tmface[i].cnct[j]].face.push_back(i);
				if (tmedge[tmface[i].edge[j]].act == 1)
				{
					tmedge[tmface[i].edge[j]].face.push_back(i);
				}
				else
				{
					cp[tmedge[tmface[i].edge[j]].midpt].face.push_back(i);
					int edchd[2] = { tmedge[tmface[i].edge[j]].chd[0], tmedge[tmface[i].edge[j]].chd[1] };
					tmedge[edchd[0]].face.push_back(i);
					tmedge[edchd[1]].face.push_back(i);
				}
			}
		}
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].act == 1)
		{
			cp[tmedge[i].pt[0]].edge.push_back(i);
			cp[tmedge[i].pt[1]].edge.push_back(i);
		}
	}
}

void TruncatedTspline_3D::InitialRotate()
{
	int pref_loc[6] = {0,0,1,3,0,4};
	//int ed_uvw[8][3] = { { 0, 3, 4 }, { 1, 0, 5 }, { 2, 1, 6 }, { 3, 2, 7 }, { 11, 8, 4 }, { 8, 9, 5 }, { 9, 10, 6 }, { 10, 11, 7 } };
	int edref_loc[6][3] = { { 0, 3, -2 }, { 0, -2, 4 }, { -1, 1, 5 }, { 2, -1, 7 }, {-2,3,4}, { 8, 11, -1 }};
	vector<Matrix3d> mat_all;
	getElementRotate(mat_all);
	for (uint i = 0; i < mat_all.size(); i++)
	{
		Matrix3d tmp = mat_all[i].transpose();
		mat_all[i] = tmp;
	}
	//set 6 rotation matrix correpsonding to each face
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		for (int i = 0; i < 6; i++)
		{
			if (tmface[tmesh[eid].face[i]].type != 1)//not boundary
			{
				int fcid(tmesh[eid].face[i]);
				int hxid(tmface[tmesh[eid].face[i]].hex[0]);
				if (hxid == eid) hxid = tmface[tmesh[eid].face[i]].hex[1];
				int ptref(tmesh[eid].cnct[pref_loc[i]]);
				int* it = find(tmesh[hxid].cnct, tmesh[hxid].cnct+8, ptref);
				int loc(it - tmesh[hxid].cnct);
				int ed_ref[3];
				for (int j = 0; j < 3; j++)
				{
					if (edref_loc[i][j] < 0) ed_ref[j] = edref_loc[i][j];
					else ed_ref[j] = tmesh[eid].edge[edref_loc[i][j]];
				}
				Matrix3d mat_loc;
				getCornerRotate(hxid,ed_uvw[loc],ed_ref,mat_loc);
				tmesh[eid].nbrot[i] = mat_loc*mat_all[loc];
				//cout << "face nb: " << hxid << "\n";
				//cout << tmesh[eid].nbrot[i] << "\n";
				//if (i == 4)
				//{
				//	cout << "local: "<<loc<<"\n";
				//	cout << mat_loc << "\n";
				//	cout << mat_all[loc]<<"\n";
				//}
				//getchar();
			}
		}
	}
}

void TruncatedTspline_3D::getElementRotate(vector<Matrix3d>& mat)//rotatoin matrix for each corner vertex
{
	mat.clear();
	mat.resize(8);
	mat[0] << 1, 0, 0,
			  0, 1, 0, 
			  0, 0, 1;
	mat[1] << 0, -1, 0, 
		      1, 0, 0, 
			  0, 0, 1;
	mat[2] << -1, 0, 0, 
		      0, -1, 0, 
			  0, 0, 1;
	mat[3] << 0, 1, 0, 
		      -1, 0, 0, 
			  0, 0, 1;
	mat[4] << 0, 1, 0, 
		      1, 0, 0, 
			  0, 0, -1;
	mat[5] << -1, 0, 0, 
		      0, 1, 0, 
			  0, 0, -1;
	mat[6] << 0, -1, 0, 
		      -1, 0, 0, 
			  0, 0, -1;
	mat[7] << 1, 0, 0, 
		      0, -1, 0, 
			  0, 0, -1;
}

void TruncatedTspline_3D::getCornerRotate(int hxid, int uvw_loc[3], int uvw_ref[3], Matrix3d& mat)
{
	int pos(0), loc(0);
	Matrix3d tmp1 = Matrix3d::Identity();
	Matrix3d tmp2 = Matrix3d::Zero();
	tmp2(0, 0) = -1; tmp2(1, 1) = -1; tmp2(2, 2) = -1;
	for (int i = 0; i < 3; i++)
	{
		if (tmesh[hxid].edge[uvw_loc[i]] == uvw_ref[0])
		{
			//mat(i, 0) = 1.; mat(i, 1) = 0.; mat(i, 2) = 0.;
			mat(0, i) = 1.; mat(1, i) = 0.; mat(2, i) = 0.;
		}
		else if (tmesh[hxid].edge[uvw_loc[i]] == uvw_ref[1])
		{
			//mat(i, 0) = 0.; mat(i, 1) = 1.; mat(i, 2) = 0.;
			mat(0, i) = 0.; mat(1, i) = 1.; mat(2, i) = 0.;
		}
		else if (tmesh[hxid].edge[uvw_loc[i]] == uvw_ref[2])
		{
			//mat(i, 0) = 0.; mat(i, 1) = 0.; mat(i, 2) = 1.;
			mat(0, i) = 0.; mat(1, i) = 0.; mat(2, i) = 1.;
		}
		else
		{
			pos = i;
		}
		if (uvw_ref[i] < 0) loc = i;
	}
	if (uvw_ref[loc] == -1)//postive
	{
		//mat(pos, 0) = tmp1(loc, 0); mat(pos, 1) = tmp1(loc, 1); mat(pos, 2) = tmp1(loc, 2);
		mat(0, pos) = tmp1(0, loc); mat(1, pos) = tmp1(1, loc); mat(2, pos) = tmp1(2, loc);
	}
	else if (uvw_ref[loc] == -2)//negative
	{
		//mat(pos, 0) = tmp2(loc, 0); mat(pos, 1) = tmp2(loc, 1); mat(pos, 2) = tmp2(loc, 2);
		mat(0, pos) = tmp2(0, loc); mat(1, pos) = tmp2(1, loc); mat(2, pos) = tmp2(2, loc);
	}
}

//void TruncatedTspline_3D::UpdateConnect()
//{
//	uint i,j,k;
//	for(i=0; i<cp.size(); i++)
//	{
//		cp[i].face.clear();
//		cp[i].edge.clear();
//	}
//	for(i=0; i<tmedge.size(); i++)
//	{
//		tmedge[i].face.clear();
//	}
//	//loop all faces
//	for(i=0; i<tmesh.size(); i++)
//	{
//		if(tmesh[i].act==1)
//		{
//			for(j=0; j<4; j++)
//			{
//				cp[tmesh[i].cnct[j]].face.push_back(i);
//				if(tmedge[tmesh[i].edge[j]].act==1)
//				{
//					tmedge[tmesh[i].edge[j]].face.push_back(i);
//				}
//				else
//				{
//					cp[tmedge[tmesh[i].edge[j]].midpt].face.push_back(i);
//					int chdid[2]={tmedge[tmesh[i].edge[j]].chd[0],tmedge[tmesh[i].edge[j]].chd[1]};
//					if(tmedge[chdid[0]].act==1 && tmedge[chdid[1]].act==1)
//					{
//						tmedge[chdid[0]].face.push_back(i);
//						tmedge[chdid[1]].face.push_back(i);
//					}
//					else
//					{
//						cerr<<"Configuration not recognized!\n";
//					}
//				}
//			}
//		}
//	}
//	//loop all edges
//	for(i=0; i<tmedge.size(); i++)
//	{
//		if(tmedge[i].act==1)
//		{
//			cp[tmedge[i].pt[0]].edge.push_back(i);
//			cp[tmedge[i].pt[1]].edge.push_back(i);
//		}
//	}
//
//	for(i=0; i<cp.size(); i++)
//	{
//		if(cp[i].face.size()==3 && cp[i].type!=2)
//		{
//			cp[i].type=1;
//		}
//	}
//
//	//FindEdgeTopoDirec_1();
//	//FindKnotInterval_1();//find kitvtmp
//	//UpdateKnotInterval_1();
//	//SetLocalCoorSystem();
//	//FindIEN_1();
//}

void TruncatedTspline_3D::FindEdgeTopoDirec()
{
	for(uint i=0; i<tmedge.size(); i++)
	{
		tmedge[i].pn[0][0]=4; tmedge[i].pn[0][1]=-1;
		tmedge[i].pn[1][0]=4; tmedge[i].pn[1][1]=-1;//initialize as end
		if(tmedge[i].act==1)
		{
			for(int j=0; j<2; j++)
			{
				int type, next;
				EdgeConnect(i,tmedge[i].pt[j],type,next);//T-junction connectivity comes from refinement later
				tmedge[i].pn[j][0]=type;
				tmedge[i].pn[j][1]=next;

				//if (i == 347 && tmedge[i].pt[j] == 142)
				//{
				//	cout << tmedge[i].pn[j][0] << " " << tmedge[i].pn[j][1]<<"\n";
				//	getchar();
				//}

				//else if(cp[tmedge[i].pt[j]].type==1)//T-junctions
				//{
				//	int fid(-1);
				//	int loc(0);
				//	for(uint k=0; k<cp[tmedge[i].pt[j]].face.size(); k++)
				//	{
				//		int ftmp(cp[tmedge[i].pt[j]].face[k]);
				//		for(int k1=0; k1<4; k1++)
				//		{
				//			if(tmedge[tmesh[ftmp].edge[k1]].act==0 && tmedge[tmesh[ftmp].edge[k1]].midpt==tmedge[i].pt[j])
				//			{
				//				loc=k1; fid=ftmp; break;
				//			}
				//		}
				//		if(fid!=-1)
				//		{
				//			break;
				//		}
				//	}
				//	if(fid==-1)
				//	{
				//		cout<<"edge id: "<<i<<"\n";
				//		cout<<tmedge[i].pt[0]<<" "<<tmedge[i].pt[1]<<"\n";
				//		cout<<tmedge[i].prt<<"\n";
				//		//cout<<cp[tmedge[i].pt[0]].face.size()<<" "<<cp[tmedge[i].pt[1]].face.size()<<"\n";
				//		cout<<cp[tmedge[i].pt[0]].face[0]<<" "<<cp[tmedge[i].pt[0]].face[1]<<" "<<cp[tmedge[i].pt[0]].face[2]<<"\n";
				//		cerr<<"T-junction cannot be found in any neighboring elements!\n";
				//		getchar();
				//	}
				//	if(tmedge[tmesh[fid].edge[loc]].chd[0]==i)
				//	{
				//		tmedge[i].pn[j][0]=0;
				//		tmedge[i].pn[j][1]=tmedge[tmesh[fid].edge[loc]].chd[1];
				//	}
				//	else if(tmedge[tmesh[fid].edge[loc]].chd[1]==i)
				//	{
				//		tmedge[i].pn[j][0]=0;
				//		tmedge[i].pn[j][1]=tmedge[tmesh[fid].edge[loc]].chd[0];
				//	}
				//	else
				//	{
				//		tmedge[i].pn[j][0]=1;
				//		tmedge[i].pn[j][1]=fid;
				//	}
				//}
			}
		}
	}
}

//void TruncatedTspline_3D::EdgeConnect(int ed0, int pid0, int& type1, int& ed1)//ed0 and pid0 is the current edge and point, type1 and ed1 for next information
//{
//	type1=4;//end
//	ed1=-1;
//	if(tmedge[ed0].type!=1 && (cp[pid0].type==1 || cp[pid0].type==12 || cp[pid0].type==13)) return;
//	if(tmedge[ed0].type==1 && cp[pid0].type==13) {type1=3; return;}
//	vector<int> edtmp;
//	for(uint i=0; i<cp[pid0].edge.size(); i++)//find the edge that not belong to hexes
//	{
//		int edid(cp[pid0].edge[i]);//edge connect to the end point pid0
//		if(edid != ed0)
//		{
//			int flag(0);
//			for(uint j=0; j<tmedge[ed0].hex.size(); j++)
//			{
//				int hexid(tmedge[ed0].hex[j]);//hex sharing this edge
//				for(uint k=0; k<12; k++)
//				{
//					int ed(tmesh[hexid].edge[k]);
//					if(tmedge[ed].act==1 && ed==edid)
//					{
//						flag=1; break;
//					}
//					else if(tmedge[ed].act==0 && (tmedge[ed].chd[0]==edid || tmedge[ed].chd[1]==edid))
//					{
//						flag=1; break;
//					}
//				}
//				for(uint k=0; k<6; k++)
//				{
//					int fc(tmesh[hexid].face[k]);
//					vector<int>::iterator it=find(tmface[fc].Tedge.begin(),tmface[fc].Tedge.end(),edid);
//					if(it!=tmface[fc].Tedge.end())
//					{
//						flag=1; break;
//					}
//				}
//			}
//			if(flag==0)
//			{
//				edtmp.push_back(edid);
//			}
//		}
//	}
//	if(edtmp.size()==1)
//	{
//		type1=0;
//		ed1=edtmp[0];
//	}
//	else if(edtmp.size()==0)//4 possible cases
//	{
//		type1=3;//case 1 - XP
//		for(uint i=0; i<cp[pid0].hex.size(); i++)
//		{
//			for(int j=0; j<6; j++)
//			{
//				int fcid(tmesh[cp[pid0].hex[i]].face[j]);
//				if(tmface[fcid].act==0 && tmface[fcid].ctpt==pid0)
//				{
//					type1=2;
//					ed1=cp[pid0].hex[i];
//					break;
//				}
//			}
//			if(type1==2) break;//case 2 - hex
//		}
//		if(type1==3)
//		{
//			vector<int> fctmp;
//			for(uint i=0; i<cp[pid0].face.size(); i++)
//			{
//				int fcid(cp[pid0].face[i]), flag(0);
//				for(uint j=0; j<tmedge[ed0].hex.size(); j++)
//				{
//					int hxid(tmedge[ed0].hex[j]);
//					for(int k=0; k<6; k++)
//					{
//						if(tmface[tmesh[hxid].face[k]].act==1 && tmesh[hxid].face[k]==i) flag=1;
//						if(tmface[tmesh[hxid].face[k]].act==0)
//						{
//							int fcid1(tmesh[hxid].face[k]);
//							vector<int>::iterator it=find(tmface[fcid1].chd.begin(),tmface[fcid1].chd.end(),i);
//							if(it!=tmface[fcid1].chd.end()) flag=1;
//						}
//						if(flag==1) break;
//					}
//					if(flag==1) break;
//				}
//				if(flag==0)
//				{
//					fctmp.push_back(fcid);
//					//type1=1; ed1=fcid; break;//case 3 - face
//				}
//			}
//			if(fctmp.size()==1)
//			{
//				type1=1; ed1=fctmp[0];
//			}
//			else if(fctmp.size()==0)
//			{
//				if(tmedge[ed0].prt!=-1)
//				{
//					for(uint i=0; i<cp[pid0].edge.size(); i++)
//					{
//						if(tmedge[cp[pid0].edge[i]].prt==tmedge[ed0].prt)
//						{
//							type1=0; ed1=cp[pid0].edge[i]; break;
//						}
//					}
//				}
//			}
//		}
//	}
//	else
//	{
//		type1=3;
//	}
//}

void TruncatedTspline_3D::EdgeConnect(int ed0, int pid0, int& type1, int& ed1)//ed0 and pid0 is the current edge and point, type1 and ed1 for next information
{
	type1 = 4;//end
	ed1 = -1;
	if (tmedge[ed0].type != 1 && (cp[pid0].type == 1 || cp[pid0].type == 12 || cp[pid0].type == 13)) return;
	if (tmedge[ed0].type == 1 && cp[pid0].type == 13) { type1 = 3; return; }
	if (tmedge[ed0].prt != -1)
	{
		if (tmedge[tmedge[ed0].prt].midpt == pid0)
		{
			type1 = 0; ed1 = tmedge[tmedge[ed0].prt].chd[0];
			if (ed1 == ed0) ed1 = tmedge[tmedge[ed0].prt].chd[1];
			return;
		}
	}
	vector<int> edtmp;
	for (uint i = 0; i<cp[pid0].edge.size(); i++)//find the edge that not belong to faces
	{
		int edid(cp[pid0].edge[i]);//edge connect to the end point pid0
		if (edid != ed0)
		{
			int flag(0);
			for (uint j = 0; j<tmedge[ed0].face.size(); j++)
			{
				int fcid(tmedge[ed0].face[j]);//active face sharing input edge
				for (uint k = 0; k<4; k++)
				{
					int ed(tmface[fcid].edge[k]);
					if (tmedge[ed].act == 1 && ed == edid)
					{
						flag = 1; break;
					}
					else if (tmedge[ed].act == 0 && (tmedge[ed].chd[0] == edid || tmedge[ed].chd[1] == edid))
					{
						flag = 1; break;
					}
				}
			}
			if (flag == 0)
			{
				edtmp.push_back(edid);
			}
		}
	}
	if (edtmp.size() == 1)
	{
		type1 = 0;
		ed1 = edtmp[0];
	}
	else if (edtmp.size() == 0)//4 possible cases
	{
		type1 = 3;//initialized as XP
		for (uint i = 0; i<cp[pid0].hex.size(); i++)//pid0 is a face T-junction
		{
			for (int j = 0; j<6; j++)
			{
				int fcid(tmesh[cp[pid0].hex[i]].face[j]);
				if (tmface[fcid].act == 0 && tmface[fcid].ctpt == pid0)
				{
					type1 = 2;
					ed1 = cp[pid0].hex[i];//next is a solid
					return;
				}
			}
		}
		vector<int> fctmp;
		for (uint i = 0; i<cp[pid0].face.size(); i++)
		{
			int fcid(cp[pid0].face[i]), flag(0);
			for (uint j = 0; j<tmedge[ed0].hex.size(); j++)
			{
				int hxid(tmedge[ed0].hex[j]);
				for (int k = 0; k<6; k++)
				{
					if (tmface[tmesh[hxid].face[k]].act == 1 && tmesh[hxid].face[k] == fcid) flag = 1;
					if (tmface[tmesh[hxid].face[k]].act == 0)
					{
						int fcid1(tmesh[hxid].face[k]);
						vector<int>::iterator it = find(tmface[fcid1].chd.begin(), tmface[fcid1].chd.end(), fcid);
						if (it != tmface[fcid1].chd.end()) flag = 1;
					}
					if (flag == 1) break;
				}
				if (flag == 1) break;
			}
			if (flag == 0)
			{
				fctmp.push_back(fcid);
				//if (pid0 == 159 && ed0 == 403)
				//{
				//	cout << fcid <<" "<< tmface[fcid].id_act << "\n";
				//	getchar();
				//}
			}
			//if (pid0 == 142 && ed0 == 347)
			//{
			//	//cout << "here " << fctmp.size() << "\n";
			//	cout << fcid <<" "<<tmface[fcid].act << "\n";
			//	getchar();
			//}
		}
		if (fctmp.size() == 1)
		{
			type1 = 1; ed1 = fctmp[0];
			//if (pid0 == 142 && ed0 == 347)
			//{
			//	cout << "here " << fctmp.size() << "\n";
			//	getchar();
			//}
			return;
		}
		//else if (fctmp.size() == 0)
		//{
		//	if (tmedge[ed0].prt != -1)
		//	{
		//		for (uint i = 0; i<cp[pid0].edge.size(); i++)
		//		{
		//			if (tmedge[cp[pid0].edge[i]].prt == tmedge[ed0].prt)
		//			{
		//				type1 = 0; ed1 = cp[pid0].edge[i]; break;
		//			}
		//		}
		//	}
		//}
	}
	else
	{
		type1 = 3;
	}
}

//void TruncatedTspline_3D::FindFaceTopoDirec()
//{
//	for(uint i=0; i<tmface.size(); i++)
//	{
//		for(int j=0; j<4; j++)
//		{
//			tmface[i].pn[j][0]=4; tmface[i].pn[j][1]=-1; 
//		}
//		if(tmface[i].act==1)
//		{
//			for(int j=0; j<4; j++)
//			{
//				int type, next;
//				FaceConnect(i,tmface[i].edge[j],type,next);
//				tmface[i].pn[j][0]=type;
//				tmface[i].pn[j][1]=next;
//			}
//		}
//	}
//}

void TruncatedTspline_3D::FaceConnect(int fc0, int ed0, int& type1, int& fc1)
{
	type1=4; fc1=-1;
	if(tmface[fc0].type!=1 && tmedge[ed0].type==1) return;
	vector<int> fctmp;
	for(uint i=0; i<tmedge[ed0].face.size(); i++)//find the edge that not belong to hexes
	{
		int fcid(tmedge[ed0].face[i]);
		if(fcid != fc0)
		{
			int flag(0);
			for(uint j=0; j<tmface[fc0].hex.size(); j++)
			{
				int hexid(tmface[fc0].hex[j]);
				for(uint k=0; k<6; k++)
				{
					int fc(tmesh[hexid].face[k]);
					if(tmface[fc].act==1 && fc==fcid)
					{
						flag=1; break;
					}
					//else if(tmedge[ed].act==0 && (tmedge[ed].chd[0]==edid || tmedge[ed].chd[1]==edid))
					//{
					//	flag=1; break;
					//}
				}
			}
			if(flag==0)
			{
				fctmp.push_back(fcid);
			}
		}
	}
	if(fctmp.size()==1)
	{
		type1=1;//next is face
		fc1=fctmp[0];
	}
	else
	{
		type1=3;//next is extraordinary edge
	}
}

void TruncatedTspline_3D::setReference()//reference element and reference edges
{
	for(uint i=0; i<cp.size(); i++)
	{
		for(int j=0; j<4; j++)
		{
			cp[i].kitvU[j]=1.; cp[i].kitvV[j]=1.; cp[i].kitvW[j]=1.;
		}
		//if(cp[i].type!=3 && cp[i].type!=13)// not extraordinaty
		{
			int pos(0);
			cp[i].rhx=-1;
			for(uint j=0; j<cp[i].hex.size(); j++)
			{
				int* it=find(tmesh[cp[i].hex[j]].cnct,tmesh[cp[i].hex[j]].cnct+8,i);
				if(it!=tmesh[cp[i].hex[j]].cnct+8)
				{
					cp[i].rhx=cp[i].hex[j];
					pos=it-tmesh[cp[i].hex[j]].cnct;
					break;
				}
			}
			if(cp[i].rhx==-1)
			{
				cerr<<"Cannot find correct reference hex!\n";
				getchar();
			}
			getElementRotate_Unit(pos,cp[i].rot_ref);
			for(int j=0; j<3; j++)
			{
				cp[i].uved[j]=tmesh[cp[i].rhx].edge[ed_uvw[pos][j]];
				if(tmedge[cp[i].uved[j]].act==0)
				{
					int edtmp(tmedge[cp[i].uved[j]].chd[0]);
					if(tmedge[edtmp].pt[0]!=i && tmedge[edtmp].pt[1]!=i)
					{
						edtmp = tmedge[cp[i].uved[j]].chd[1];
					}
					cp[i].uved[j] = edtmp;
				}
				vector<int>::iterator it=find(cp[i].edge.begin(),cp[i].edge.end(),cp[i].uved[j]);
				if(it==cp[i].edge.end())
				{
					cerr<<"Cannot find correct uv edges!\n";
					//cout << "pid: " << i << "\n";
					//cout << "edge act: " << tmedge[tmesh[cp[i].rhx].edge[ed_uvw[pos][j]]].act << "\n";
					cout << "edge act: " << tmedge[tmesh[cp[i].rhx].edge[ed_uvw[pos][j]]].pt[0] << " " << tmedge[tmesh[cp[i].rhx].edge[ed_uvw[pos][j]]].pt[1] << "\n";
					//cout << "edge chd: " << tmedge[tmedge[tmesh[cp[i].rhx].edge[ed_uvw[pos][j]]].chd[0]].id_act << " " << tmedge[tmedge[tmesh[cp[i].rhx].edge[ed_uvw[pos][j]]].chd[1]].id_act << "\n";
					//cout << "edge: " << tmedge[cp[i].uved[j]].id_act << "\n";
					getchar();
				}
			}
		}
	}
}

void TruncatedTspline_3D::FindKnotInterval()
{
	setReference();
	//int fc_ppd_ed[6]={4,3,0,3,0,4};
	int trun_flag;
	for(uint i=0; i<cp.size(); i++)
	{
		//if(cp[i].type!=1 /*&& cp[i].type!=13*/)
		{
			ShootRay(i,cp[i].uved[0],cp[i].kitvU,trun_flag);
			ShootRay(i,cp[i].uved[1],cp[i].kitvV,trun_flag);
			ShootRay(i,cp[i].uved[2],cp[i].kitvW,trun_flag);
			//cp[i].trun=trun_flag;
		}
		//else if(cp[i].type==3)
		//{
		//	for(int j=0; j<4; j++)
		//	{
		//		cp[i].kitvU[j]=tmedge[cp[i].edge[0]].len;
		//		cp[i].kitvV[j]=tmedge[cp[i].edge[0]].len;
		//		cp[i].kitvW[j]=tmedge[cp[i].edge[0]].len;
		//	}
		//}
		//else if(cp[i].type==13)
		//{
		//}
	}
}

void TruncatedTspline_3D::ShootRay(int pid, int edid, double kv[4], int& trun_flag)
{
	int loc0(0),loc1(1);
	//int edge_lev[2]={-1,-1};
	//trun_flag=0;
	//int flag[2]={-1,-1};
	//positive direction
	if(pid!=tmedge[edid].pt[0])
	{
		loc0=1; loc1=0;
	}
	kv[2]=tmedge[edid].len;
	//edge_lev[0]=tmedge[edid].lev;//first edge level
	if(tmedge[edid].pn[loc1][0]==0)//next is edge
	{
		kv[3]=tmedge[tmedge[edid].pn[loc1][1]].len;
		////find skip
		//int edtmp(tmedge[edid].pn[loc1][1]);
		//edge_lev[1]=tmedge[edtmp].lev;//second edge level
		//if(edge_lev[1]>edge_lev[0] && tmedge[edtmp].face.size()==2)
		//{
		//	int eid[2]={tmedge[edtmp].face[0],tmedge[edtmp].face[1]};
		//	if(tmesh[eid[0]].lev!=tmesh[eid[1]].lev)
		//	{
		//		int e_fe(eid[0]);
		//		if(tmesh[eid[0]].lev > tmesh[eid[1]].lev) e_fe=eid[1];
		//		int pos(0);
		//		for(int k=0; k<4; k++)
		//		{
		//			if(tmedge[tmesh[e_fe].edge[k]].act==0)
		//			{
		//				pos=k;
		//				break;
		//			}
		//		}
		//		if(cp[tmesh[e_fe].cnct[(pos+3)%4]].type==1)
		//		{
		//			flag[0]=1;
		//		}
		//	}
		//}
	}
	else if(tmedge[edid].pn[loc1][0]==1)//next is face
	{
		int fid(tmedge[edid].pn[loc1][1]), pos(0);
		for(int i=0; i<4; i++)
		{
			if(tmedge[tmface[fid].edge[i]].act==0 && tmedge[tmface[fid].edge[i]].midpt==tmedge[edid].pt[loc1])
			{
				pos=i; break;
			}
		}
		kv[3]=tmedge[tmface[fid].edge[(pos+1)%4]].len;
	}
	else if(tmedge[edid].pn[loc1][0]==2)//next is hex
	{
		int hxid(tmedge[edid].pn[loc1][1]), pos(0);
		for(int i=0; i<6; i++)
		{
			if(tmface[tmesh[hxid].face[i]].act==0 && tmface[tmesh[hxid].face[i]].ctpt==tmedge[edid].pt[loc1])
			{
				pos=i; break;
			}
		}
		kv[3]=tmedge[tmesh[hxid].edge[fc_ppd_ed[pos]]].len;
	}
	else if(tmedge[edid].pn[loc1][0]==3)//next is XP
	{
		kv[3]=kv[2];
	}
	else if(tmedge[edid].pn[loc1][0]==4)//end
	{
		//kv[3]=0.;
		kv[3] = 1.;
	}
	//if(flag[0]==1)//skip
	//{
	//	kv[3]=2.*kv[3];
	//	trun_flag=1;
	//	//cout<<"trun!\n";
	//	//getchar();
	//}
	//negative direction
	//edge_lev[0]=-1; edge_lev[1]=-1; 
	if(tmedge[edid].pn[loc0][0]==0)//previous is edge
	{
		int ed0=tmedge[edid].pn[loc0][1];
		kv[1]=tmedge[ed0].len;
		//edge_lev[0]=tmedge[ed0].lev;//first edge level
		int a0=0, a1=1;
		if(tmedge[ed0].pt[0]!=pid)
		{
			a0=1; a1=0;
		}
		if(tmedge[ed0].pn[a1][0]==0)//previous previous is edge
		{
			kv[0]=tmedge[tmedge[ed0].pn[a1][1]].len;
			////find skip
			//int edtmp(tmedge[ed0].pn[a1][1]);
			//edge_lev[1]=tmedge[edtmp].lev;//second edge level
			//if(edge_lev[1]>edge_lev[0] && tmedge[edtmp].face.size()==2)
			//{
			//	int eid[2]={tmedge[edtmp].face[0],tmedge[edtmp].face[1]};
			//	if(tmesh[eid[0]].lev!=tmesh[eid[1]].lev)
			//	{
			//		int e_fe(eid[0]);
			//		if(tmesh[eid[0]].lev > tmesh[eid[1]].lev) e_fe=eid[1];
			//		int pos(0);
			//		for(int k=0; k<4; k++)
			//		{
			//			if(tmedge[tmesh[e_fe].edge[k]].act==0)
			//			{
			//				pos=k;
			//				break;
			//			}
			//		}
			//		if(cp[tmesh[e_fe].cnct[(pos+3)%4]].type==1)
			//		{
			//			flag[1]=1;
			//		}
			//	}
			//}
		}
		else if(tmedge[ed0].pn[a1][0]==1)//previous previous is face
		{
			int pt0(tmedge[ed0].pt[a1]), fid(tmedge[ed0].pn[a1][1]), pos(0);
			for(int i=0; i<4; i++)
			{
				if(tmedge[tmface[fid].edge[i]].act==0 && tmedge[tmface[fid].edge[i]].midpt==tmedge[ed0].pt[a1])
				{
					pos=i; break;
				}
			}
			kv[0]=tmedge[tmface[fid].edge[(pos+1)%4]].len;
		}
		else if(tmedge[ed0].pn[a1][0]==2)//previous previous is hex
		{
			int hxid(tmedge[ed0].pn[a1][1]), pos(0);
			for(int i=0; i<6; i++)
			{
				if(tmface[tmesh[hxid].face[i]].act==0 && tmface[tmesh[hxid].face[i]].ctpt==tmedge[ed0].pt[a1])
				{
					pos=i; break;
				}
			}
			kv[0]=tmedge[tmesh[hxid].edge[fc_ppd_ed[pos]]].len;
		}
		else if(tmedge[ed0].pn[a1][0]==3)
		{
			kv[0]=kv[1];
		}
		else if(tmedge[ed0].pn[a1][0]==4)
		{
			//kv[0]=0.;
			kv[0] = 1.;
		}
		
	}
	else if(tmedge[edid].pn[loc0][0]==1)//previous is face
	{
		int fid0(tmedge[edid].pn[loc0][1]), pos(0);
		for(int i=0; i<4; i++)
		{
			if(tmedge[tmface[fid0].edge[i]].act==0 && tmedge[tmface[fid0].edge[i]].midpt==pid)
			{
				pos=i; break;
			}
		}
		kv[1]=tmedge[tmface[fid0].edge[(pos+1)%4]].len;

		//if (pid == 142)
		//{
		//	cout << "here "<< kv[1] << "\n";
		//	getchar();
		//}

		int ed0(tmface[fid0].edge[(pos+2)%4]);//find opposite edge, next connetivity, can be either face or hex if ed0 is regular
		if(tmface[fid0].type!=1 && tmedge[ed0].type==1)//boundary
		{
			//kv[0]=0.;
			kv[0] = 1.;
		}
		else if(tmedge[ed0].type==2)//extraordinary
		{
			kv[0]=kv[1];
		}
		else if (tmface[fid0].prt != -1)//only consider two children first
		{
			int fid1(tmface[tmface[fid0].prt].chd[0]);
			if (fid1 == fid0) fid1 = tmface[tmface[fid0].prt].chd[1];
			int* it1 = find(tmface[fid1].edge, tmface[fid1].edge+4, ed0);
			if (it1 != tmface[fid1].edge + 4)
			{
				int pos1(it1 - tmface[fid1].edge);
				kv[0] = tmedge[tmface[fid1].edge[(pos+1)%4]].len;
			}
			else
			{
				cerr << "Can't find shared edge!\n";
				getchar();
			}
		}
		else// if(tmedge[ed0].type==0 || tmedge[ed0].type==3)//T-edge, not consider T-face
		{
			vector<int> fctmp;
			for(uint i=0; i<tmedge[ed0].face.size(); i++)
			{
				int fcid(tmedge[ed0].face[i]), flag(0);
				for(uint j=0; j<tmface[fcid].hex.size(); j++)
				{
					int hxid(tmface[fcid].hex[j]);
					for(int k=0; k<6; k++)
					{
						if(tmface[tmesh[hxid].face[k]].act==1 && tmesh[hxid].face[k]==fcid)
						{
							flag=1; break;
						}
						else if(tmface[tmesh[hxid].face[k]].act==0)
						{
							vector<int>::iterator it=find(tmface[tmesh[hxid].face[k]].chd.begin(),tmface[tmesh[hxid].face[k]].chd.end(),fcid);
							if(it!=tmface[tmesh[hxid].face[k]].chd.end())
							{
								flag=1; break;
							}
						}
					}
					if(flag==1) break;
				}
				if(flag==0) fctmp.push_back(fcid);
			}
			if(fctmp.size()==1)
			{
				int* it=find(tmface[fctmp[0]].edge,tmface[fctmp[0]].edge+4,ed0);
				if(it!=tmface[fctmp[0]].edge+4)
				{
					int pos=it-tmface[fctmp[0]].edge;
					kv[0]=tmedge[tmface[fctmp[0]].edge[(pos+1)%4]].len;
				}
				else
				{
					cerr<<"Cannot find next face that contains this edge!\n";
					getchar();
				}
			}
			else if(fctmp.size()==0)
			{
				kv[0]=kv[1];
				int hxid(-1), pos(-1);
				for(uint i=0; i<tmedge[ed0].hex.size(); i++)
				{
					for(int j=0; j<6; j++)
					{
						if(tmface[tmesh[tmedge[ed0].hex[i]].face[j]].act==0)
						{
							int fcid(tmesh[tmedge[ed0].hex[i]].face[j]);
							vector<int>::iterator it=find(tmface[fcid].Tedge.begin(),tmface[fcid].Tedge.end(),ed0);
							if(it!=tmface[fcid].Tedge.end())
							{
								hxid=tmedge[ed0].hex[i]; pos=j; break;
							}
						}
					}
					if(hxid!=-1) break;
				}
				if(hxid!=-1)
				{
					kv[0]=tmedge[tmesh[hxid].edge[fc_ppd_ed[pos]]].len;
				}
			}
		}

		//if(tmface[fid0].pn[ed0][0]==0)//previous previous is edge
		//{
		//	//do not allow currently
		//	kv[0]=tmedge[tmface[fid0].pn[ed0][1]].len;
		//}
		//else if(tmface[fid0].pn[ed0][0]==1)//previous previous is face
		//{
		//	//int fid1(tmface[fid0].pn[ed0][1]);
		//	//int* it=find(tmface[fid1].edge,tmface[fid1].edge+4,ed0);
		//	//pos=it-tmface[fid1].edge;
		//	//kv[0]=tmedge[tmface[fid1].edge[(pos+1)%4]].len;
		//}
		//else if(tmface[fid0].pn[ed0][0]==2)//previous previous is hex
		//{
		//	if(tmedge[tmface[fid0].edge[ed0]].act==1)
		//	{
		//	}
		//	else
		//	{
		//	}
		//}
		//else if(tmface[fid0].pn[ed0][0]==3)//previous previous is XP
		//{
		//	kv[0]=kv[1];
		//}
		//else if(tmface[fid0].pn[ed0][0]==4)//previous previous is end
		//{
		//	kv[0]=0.;
		//}
	}
	else if(tmedge[edid].pn[loc0][0]==2)//previous is hex
	{
		int hxid(tmedge[edid].pn[loc0][1]), pos(0);
		for(int i=0; i<6; i++)
		{
			if(tmface[tmesh[hxid].face[i]].act==0 && tmface[tmesh[hxid].face[i]].ctpt==tmedge[edid].pt[loc0])
			{
				pos=i; break;
			}
		}
		kv[1]=tmedge[tmesh[hxid].edge[fc_ppd_ed[pos]]].len;
		int fc0(tmesh[hxid].face[fc_opst[pos]]);
		if(tmface[fc0].act==1)
		{
			if(tmface[fc0].hex.size()==2)
			{
				int hxid1(tmface[fc0].hex[0]);
				if(hxid1==hxid) hxid1=tmface[fc0].hex[1];
				int* it=find(tmesh[hxid1].face,tmesh[hxid1].face+6,fc0);
				pos=it-tmesh[hxid1].face;
				kv[0]=tmedge[tmesh[hxid1].edge[fc_ppd_ed[pos]]].len;
			}
			else if(tmface[fc0].hex.size()==1)
			{
				//kv[0]=0.;
				kv[0] = 1.;
			}
		}
		else
		{
			//next could be edge, face, hex, end
		}
	}
	else if(tmedge[edid].pn[loc0][0]==3)//previous is XP
	{
		kv[1]=kv[2]; kv[0]=kv[2];
	}
	else
	{
		//kv[1]=0.; kv[0]=0.;
		kv[1] = 1.; kv[0] = 1.;
	}
	//if(flag[1]==1)//skip
	//{
	//	kv[0]=2.*kv[0];
	//	trun_flag=1;
	//	//cout<<"trun!\n";
	//	//getchar();
	//}
}

void TruncatedTspline_3D::SetLocalCoorSystem()
{
	for(uint eid=0; eid<tmesh.size(); eid++)
	{
		vector<int>().swap(tmesh[eid].node);
		vector<Matrix4d>().swap(tmesh[eid].lcs);
		if(tmesh[eid].act==1 /*&& tmesh[eid].type!=1*/)//allow boundary and try
		{
			double ul0[3] = { tmedge[tmesh[eid].edge[0]].len, tmedge[tmesh[eid].edge[1]].len, tmedge[tmesh[eid].edge[4]].len};
			double ul1[3] = { tmedge[tmesh[eid].edge[0]].len/2., tmedge[tmesh[eid].edge[1]].len/2., tmedge[tmesh[eid].edge[4]].len/2.};
			double uvcn[8][3] = { { 0., 0., 0. }, { ul0[0], 0., 0. }, { ul0[0], ul0[1], 0. }, { 0., ul0[1], 0. }, 
			{ 0., 0., ul0[2] }, { ul0[0], 0., ul0[2] }, { ul0[0], ul0[1], ul0[2] }, { 0., ul0[1], ul0[2] } };//corner point
			double uved[12][3] = { { ul1[0], 0., 0. }, { ul0[0], ul1[1], 0. }, { ul1[0], ul0[1], 0. }, {0.,ul1[1],0.},
			{ 0., 0., ul1[2] }, { ul0[0], 0., ul1[2] }, { ul0[0], ul0[1], ul1[2] }, { 0., ul0[1], ul1[2] },
			{ ul1[0], 0., ul0[2] }, { ul0[0], ul1[1], ul0[2] }, { ul1[0], ul0[1], ul0[2] }, { 0., ul1[1], ul0[2] } };//eddge point
			double uvfc[6][3] = { { ul1[0], ul1[1], 0. }, { ul1[0], 0., ul1[2] }, { ul0[0], ul1[1], ul1[2] }, 
			{ ul1[0], ul0[1], ul1[2] }, { 0., ul1[1], ul1[2] }, { ul1[0], ul1[1], ul0[2] } };//face point
			for (int i = 0; i < 8; i++)
			{
				Matrix3d tmp1;
				Matrix4d tmp2=Matrix4d::Zero();
				Find_Neighbor_Rot(eid,tmesh[eid].cnct[i],tmp1);
				for (int j = 0; j < 3; j++)
				{
					for (int k = 0; k < 3; k++)	tmp2(j, k) = tmp1(j, k);
				}
				tmp2(0, 3) = uvcn[i][0]; tmp2(1, 3) = uvcn[i][1]; tmp2(2, 3) = uvcn[i][2]; tmp2(3, 3) = 1.;
				tmesh[eid].node.push_back(tmesh[eid].cnct[i]);
				tmesh[eid].lcs.push_back(tmp2);
				//cout << "pid: " << tmesh[eid].cnct[i] << "\n";
				//cout << "point ref: " << cp[tmesh[eid].cnct[i]].rhx << "\n";
				//cout << tmp2 << "\n";
				//getchar();
			}
			for (int i = 0; i < 12; i++)
			{
				if (tmedge[tmesh[eid].edge[i]].act==0)
				{
					Matrix3d tmp1;
					Matrix4d tmp2 = Matrix4d::Zero();
					Find_Neighbor_Rot(eid, tmedge[tmesh[eid].edge[i]].midpt, tmp1);
					for (int j = 0; j < 3; j++)
					{
						for (int k = 0; k < 3; k++)	tmp2(j, k) = tmp1(j, k);
					}
					tmp2(0, 3) = uved[i][0]; tmp2(1, 3) = uved[i][1]; tmp2(2, 3) = uved[i][2]; tmp2(3, 3) = 1.;
					tmesh[eid].node.push_back(tmedge[tmesh[eid].edge[i]].midpt);
					tmesh[eid].lcs.push_back(tmp2);
				}
			}
			for (int i = 0; i < 6; i++)
			{
				if (tmface[tmesh[eid].face[i]].act == 0 && tmface[tmesh[eid].face[i]].ctpt!=-1)
				{
					Matrix3d tmp1;
					Matrix4d tmp2 = Matrix4d::Zero();
					Find_Neighbor_Rot(eid, tmface[tmesh[eid].face[i]].ctpt, tmp1);
					for (int j = 0; j < 3; j++)
					{
						for (int k = 0; k < 3; k++)	tmp2(j, k) = tmp1(j, k);
					}
					tmp2(0, 3) = uvfc[i][0]; tmp2(1, 3) = uvfc[i][1]; tmp2(2, 3) = uvfc[i][2]; tmp2(3, 3) = 1.;
					tmesh[eid].node.push_back(tmface[tmesh[eid].face[i]].ctpt);
					tmesh[eid].lcs.push_back(tmp2);
				}
			}
		}
	}
}

void TruncatedTspline_3D::Find_Neighbor_Rot(int hxid, int pid, Matrix3d& rot)
{
	if (cp[pid].rhx == hxid)//pid must be a corner point
	{
		rot = cp[pid].rot_ref;
		return;
	}
	//rhx share a face with hxid
	for (int i = 0; i < 6; i++)
	{
		if (tmface[tmesh[hxid].face[i]].act == 1)
		{
			int hxnb(tmface[tmesh[hxid].face[i]].hex[0]);
			if (hxnb == hxid) hxnb = tmface[tmesh[hxid].face[i]].hex[1];
			if (hxnb == cp[pid].rhx)
			{
				rot = tmesh[hxid].nbrot[i]*cp[pid].rot_ref;
				return;
			}
		}
		else
		{
			for (uint j = 0; j < tmface[tmesh[hxid].face[i]].chd.size(); j++)
			{
				int fcid(tmface[tmesh[hxid].face[i]].chd[j]);
				int hxnb(tmface[fcid].hex[0]);
				if (hxnb == hxid) hxnb = tmface[fcid].hex[1];
				if (hxnb == cp[pid].rhx)
				{
					rot = tmesh[hxid].nbrot[i] * cp[pid].rot_ref;
					return;
				}
			}
		}
	}
	//rhx does not share face with hxid, find intermediate hex first
	vector<int> remain;
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		if (cp[pid].hex[i] != hxid) remain.push_back(cp[pid].hex[i]);
	}
	int found(0);
	int cur(hxid), next, fc_loc;
	Matrix3d tmp = Matrix3d::Identity();
	Matrix3d tmp1;
	while (remain.size()!=0 && found==0)
	{
		next = -1; fc_loc = 0;
		for (int i = 0; i < 6; i++)
		{
			if (tmface[tmesh[cur].face[i]].type==0)//non-boundary
			{
				if (tmface[tmesh[cur].face[i]].act == 1)
				{
					int nb(tmface[tmesh[cur].face[i]].hex[0]);
					if (nb == cur) nb = tmface[tmesh[cur].face[i]].hex[1];
					vector<int>::iterator it = find(remain.begin(), remain.end(), nb);
					if (it != remain.end())
					{
						next = nb; fc_loc = i;
					}
				}
				else
				{
					for (uint j = 0; j < tmface[tmesh[cur].face[i]].chd.size(); j++)
					{
						int fcid(tmface[tmesh[cur].face[i]].chd[j]);
						int nb(tmface[fcid].hex[0]);
						if (nb == cur) nb = tmface[fcid].hex[1];
						vector<int>::iterator it = find(remain.begin(), remain.end(), nb);
						if (it != remain.end())
						{
							next = nb; fc_loc = i; break;
						}
					}
				}
				if (next != -1) break;
			}
		}
		tmp1 = tmp*tmesh[cur].nbrot[fc_loc];
		tmp = tmp1;
		if (next == cp[pid].rhx)
		{
			found = 1; break;
		}
		vector<int> rtmp(remain);
		remain.clear();
		for (uint i = 0; i < rtmp.size(); i++)
		{
			if (rtmp[i] != cur) remain.push_back(rtmp[i]);
		}
		cur = next;
	}
	rot = tmp*cp[pid].rot_ref;//
	if (found == 0)
	{
		cout << "Neighbor not found!\n";
		getchar();
	}
}

void TruncatedTspline_3D::getElementRotate_Unit(int loc, Matrix3d& rot)
{
	if (loc == 0) rot << 1, 0, 0, 0, 1, 0, 0, 0, 1;
	else if (loc == 1) rot << 0, -1, 0, 1, 0, 0, 0, 0, 1;
	else if (loc == 2) rot << -1, 0, 0, 0, -1, 0, 0, 0, 1;
	else if (loc == 3) rot << 0, 1, 0, -1, 0, 0, 0, 0, 1;
	else if (loc == 4) rot << 0, 1, 0, 1, 0, 0, 0, 0, -1;
	else if (loc == 5) rot << -1, 0, 0, 0, 1, 0, 0, 0, -1;
	else if (loc == 6) rot << 0, -1, 0, -1, 0, 0, 0, 0, -1;
	else rot << 1, 0, 0, 0, -1, 0, 0, 0, -1;
}

//void TruncatedTspline_3D::FindIEN_Unstruct()//IEN, not IENtmp
//{
//	for(uint eid=0; eid<tmesh.size(); eid++)
//	{
//		tmesh[eid].IEN.clear();
//		tmesh[eid].patch_ku.clear();
//		tmesh[eid].patch_kv.clear();
//		tmesh[eid].patch_kw.clear();
//		if(tmesh[eid].act==1 && tmesh[eid].type==0)//find two ring neighorhood
//		{
//			array<double,2> urang={0.,tmedge[tmesh[eid].edge[0]].len};
//			array<double,2> vrang={0.,tmedge[tmesh[eid].edge[3]].len};
//			array<double,2> wrang={0.,tmedge[tmesh[eid].edge[4]].len};
//			int count(0);
//			//initial
//			vector<int> pr0(tmesh[eid].node),er0(1,eid),pr1,er1,pr1_pref,pr1_eref;
//			vector<int> rot_ref(tmesh[eid].node.size());
//			vector<array<double,2>> uv_ref(tmesh[eid].node.size());
//			for(uint i=0; i<tmesh[eid].node.size(); i++)
//			{
//				rot_ref[i]=tmesh[eid].lcs[i].rot;
//				uv_ref[i][0]=tmesh[eid].lcs[i].u[0];
//				uv_ref[i][1]=tmesh[eid].lcs[i].u[1];
//			}
//			for(uint i=0; i<tmesh[eid].node.size(); i++)
//			{
//				array<double,5> kui,kvi;
//				FindLocalKnotVector(tmesh[eid].node[i],rot_ref[i],uv_ref[i],kui,kvi);
//				if(CheckSupport(urang,vrang,kui,kvi))
//				{
//					tmesh[eid].IEN.push_back(tmesh[eid].node[i]);
//					tmesh[eid].patch_ku.push_back(kui);
//					tmesh[eid].patch_kv.push_back(kvi);
//				}
//			}
//			while(count<2)
//			{
//				FindNextRing(pr0,er0,pr1,er1,pr1_pref,pr1_eref);
//				vector<int> rot_tmp(pr1.size());
//				vector<array<double,2>> uv_tmp(pr1.size());
//				for(uint i=0; i<pr1.size(); i++)
//				{
//					array<double,5> kui,kvi;
//					FindRotateAndUVCoor(pr0[pr1_pref[i]],rot_ref[pr1_pref[i]],uv_ref[pr1_pref[i]],pr1_eref[i],pr1[i],rot_tmp[i],uv_tmp[i]);
//					FindLocalKnotVector(pr1[i],rot_tmp[i],uv_tmp[i],kui,kvi);
//					if(CheckSupport(urang,vrang,kui,kvi))
//					{
//						tmesh[eid].IEN.push_back(pr1[i]);
//						tmesh[eid].patch_ku.push_back(kui);
//						tmesh[eid].patch_kv.push_back(kvi);
//					}
//				}
//				pr0.clear();
//				er0.clear();
//				rot_ref.clear();
//				uv_ref.clear();
//				pr0=pr1;
//				er0=er1;
//				rot_ref=rot_tmp;
//				uv_ref=uv_tmp;
//				count++;
//			}
//		}
//		else if(tmesh[eid].act==1 && tmesh[eid].type==4)
//		{
//			tmesh[eid].IEN.push_back(tmesh[eid].cnct[0]);
//			int fc_pre(tmedge[tmesh[eid].edge[3]].face[0]);
//			if(fc_pre==eid) fc_pre=tmedge[tmesh[eid].edge[3]].face[1];
//			tmesh[eid].IEN.push_back(tmesh[fc_pre].cnct[3]);
//			tmesh[eid].IEN.push_back(tmesh[fc_pre].cnct[2]);
//			tmesh[eid].IEN.push_back(tmesh[eid].cnct[3]);
//			tmesh[eid].IEN.push_back(tmesh[eid].cnct[2]);
//			int fc_next(tmedge[tmesh[eid].edge[0]].face[0]);
//			if(fc_next==eid) fc_next=tmedge[tmesh[eid].edge[0]].face[1];
//			int fc_next0=fc_next;
//			while(fc_next!=fc_pre)
//			{
//				tmesh[eid].IEN.push_back(tmesh[fc_next].cnct[3]);
//				tmesh[eid].IEN.push_back(tmesh[fc_next].cnct[2]);
//				int fc_nn(tmedge[tmesh[fc_next].edge[0]].face[0]);
//				if(fc_nn==fc_next) fc_nn=tmedge[tmesh[fc_next].edge[0]].face[1];
//				fc_next=fc_nn;
//			}
//			for(int j=1; j<4; j++)
//			{
//				for(uint k=0; k<cp[tmesh[eid].cnct[j]].face.size(); k++)
//				{
//					int fcid(cp[tmesh[eid].cnct[j]].face[k]);
//					if(fcid!=eid && fcid!=fc_pre && fcid!=fc_next0)
//					{
//						for(uint k1=0; k1<tmesh[fcid].node.size(); k1++)
//						{
//							vector<int>::iterator it=find(tmesh[eid].IEN.begin(),tmesh[eid].IEN.end(),tmesh[fcid].node[k1]);
//							if(it==tmesh[eid].IEN.end())
//							{
//								array<double,2> uv_ref={tmesh[eid].lcs[j].u[0],tmesh[eid].lcs[j].u[1]},uv;
//								int rot;
//								array<double,5> kui, kvi;
//								FindRotateAndUVCoor(tmesh[eid].cnct[j],tmesh[eid].lcs[j].rot,uv_ref,fcid,tmesh[fcid].node[k1],rot,uv);
//								FindLocalKnotVector(tmesh[fcid].node[k1],rot,uv,kui,kvi);
//								tmesh[eid].IEN.push_back(tmesh[fcid].node[k1]);
//								tmesh[eid].patch_ku.push_back(kui);
//								tmesh[eid].patch_kv.push_back(kvi);
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//}

void TruncatedTspline_3D::FindIEN_PatchKV()
{
	for(uint eid=0; eid<tmesh.size(); eid++)
	{
		vector<int>().swap(tmesh[eid].IENtmp);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kutmp);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kvtmp);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kwtmp);
		if (tmesh[eid].act == 1 && (tmesh[eid].type == 0))//not for boundary, find two ring neighorhood
		{
			array<double,2> urang={0.,tmedge[tmesh[eid].edge[0]].len};
			array<double,2> vrang={0.,tmedge[tmesh[eid].edge[3]].len};
			array<double, 2> wrang = { 0., tmedge[tmesh[eid].edge[4]].len };
			int count(0);
			//initial
			vector<int> pr0(tmesh[eid].node),er0(1,eid),pr1,er1,pr1_pref,pr1_eref;
			vector<Matrix4d> lcs0(tmesh[eid].lcs);
			for(uint i=0; i<pr0.size(); i++)
			{
				array<double,5> kui,kvi,kwi;
				FindLocalKnotVector(pr0[i], lcs0[i], kui, kvi, kwi);
				//cout <<"element id: "<< eid << "\n";
				////cout << lcs0[i] << "\n";
				//cout << pr0[i] << "\n";
				//cout << "referent element: " << cp[pr0[i]].rhx << "\n";
				//cout << "after\n";
				//cout << kui[0] << " " << kui[1] << " " << kui[2] << " " << kui[3] << " " << kui[4] << "\n";
				//cout << kvi[0] << " " << kvi[1] << " " << kvi[2] << " " << kvi[3] << " " << kvi[4] << "\n";
				//cout << kwi[0] << " " << kwi[1] << " " << kwi[2] << " " << kwi[3] << " " << kwi[4] << "\n";
				//getchar();
				if(CheckSupport(urang,vrang,wrang,kui,kvi,kwi))
				{
					//if(tmesh[eid].node[i]>=npt_old) tmesh[eid].aff=1;
					tmesh[eid].IENtmp.push_back(tmesh[eid].node[i]);
					tmesh[eid].patch_kutmp.push_back(kui);
					tmesh[eid].patch_kvtmp.push_back(kvi);
					tmesh[eid].patch_kwtmp.push_back(kwi);
				}
			}
			while(count<2)
			{
				FindNextRing(pr0,er0,pr1,er1,pr1_pref,pr1_eref);
				vector<Matrix4d> lcs1(pr1.size());
				//cout <<"pr1 size: "<< pr1.size()<<"\n";
				//cout << "er1 size: " << er1.size() << "\n";
				//getchar();
				for(uint i=0; i<pr1.size(); i++)
				{
					array<double,5> kui,kvi,kwi;
					TranslateLCS(pr0[pr1_pref[i]], lcs0[pr1_pref[i]], pr1_eref[i], pr1[i], lcs1[i]);
					FindLocalKnotVector(pr1[i],lcs1[i],kui,kvi,kwi);
					//cout << pr1[i] << "\n";
					//cout << "ref eid: "<<cp[pr1[i]].rhx<<"\n";
					//cout << "lcs1:\n";
					//cout << lcs1[i] << "\n\n";
					//cout << kui[0] << " " << kui[1] << " " << kui[2] << " " << kui[3] << " " << kui[4] << "\n";
					//cout << kvi[0] << " " << kvi[1] << " " << kvi[2] << " " << kvi[3] << " " << kvi[4] << "\n";
					//cout << kwi[0] << " " << kwi[1] << " " << kwi[2] << " " << kwi[3] << " " << kwi[4] << "\n";
					//getchar();
					if (CheckSupport(urang, vrang, wrang,kui, kvi,kwi))
					{
						//if(pr1[i]>=npt_old) tmesh[eid].aff=1;
						tmesh[eid].IENtmp.push_back(pr1[i]);
						tmesh[eid].patch_kutmp.push_back(kui);
						tmesh[eid].patch_kvtmp.push_back(kvi);
						tmesh[eid].patch_kwtmp.push_back(kwi);
					}
				}
				pr0.clear();
				er0.clear();
				vector<Matrix4d>().swap(lcs0);
				pr0=pr1;
				er0=er1;
				lcs0 = lcs1;
				vector<Matrix4d>().swap(lcs1);
				count++;
			}
		}
		else if (tmesh[eid].act == 1 && tmesh[eid].type == 2)//irregular element with extraordinary points
		{
			SetBezierMatIrrPatch(eid);
		}
		//else if (tmesh[eid].act == 1 && tmesh[eid].type == 1)//boundart element
		//{
		//	//will later be updated by ConstructBezier_Boundary
		//}
	}
}

void TruncatedTspline_3D::FindNextRing(const vector<int>& pr0, const vector<int>& er0, vector<int>& pr1, vector<int>& er1, vector<int>& pr1_pref, vector<int>& pr1_eref)
{
	pr1.clear();//next ring point
	er1.clear();//next ring hex
	pr1_pref.clear();//store certain local index of pr0
	pr1_eref.clear();//store the first neighbor hex containing pr1_pref
	for(uint i=0; i<pr0.size(); i++)
	{
		if (cp[pr0[i]].type != 3 && cp[pr0[i]].type != 13)
		{
			for(uint j=0; j<cp[pr0[i]].hex.size(); j++)
			{
				int hxid(cp[pr0[i]].hex[j]);
				vector<int>::const_iterator it1=find(er0.begin(),er0.end(),hxid);
				vector<int>::iterator it4 = find(er1.begin(), er1.end(), hxid);
				if (it1 == er0.end() && it4 == er1.end())
				{
					er1.push_back(hxid);
					for(uint k=0; k<tmesh[hxid].node.size(); k++)
					{
						vector<int>::const_iterator it2=find(pr0.begin(),pr0.end(),tmesh[hxid].node[k]);
						vector<int>::iterator it3=find(pr1.begin(),pr1.end(),tmesh[hxid].node[k]);
						if(it2==pr0.end() && it3==pr1.end())
						{
							pr1.push_back(tmesh[hxid].node[k]);
							pr1_pref.push_back(i);
							pr1_eref.push_back(hxid);
						}
					}
				}
			}
		}
	}
}

void TruncatedTspline_3D::TranslateLCS(int pref, Matrix4d& lcs_ref, int eid, int pid, Matrix4d& lcs)
{
	vector<int>::iterator it0=find(tmesh[eid].node.begin(),tmesh[eid].node.end(),pref);
	vector<int>::iterator it1=find(tmesh[eid].node.begin(),tmesh[eid].node.end(),pid);
	int loc0=it0-tmesh[eid].node.begin();
	int loc1=it1-tmesh[eid].node.begin();
	//find inverse of pref local eid
	Matrix4d ref_inv,mat1;
	getLCS_inverse(tmesh[eid].lcs[loc0],ref_inv);
	mat1 = ref_inv*tmesh[eid].lcs[loc1];
	lcs = lcs_ref*mat1;
}

void TruncatedTspline_3D::getLCS_inverse(const Matrix4d& mat_in, Matrix4d& mat_out)
{
	mat_out = mat_in.transpose();
	mat_out(0, 3) = -mat_in(0, 0)*mat_in(0, 3) - mat_in(1, 0)*mat_in(1, 3) - mat_in(2, 0)*mat_in(2, 3);
	mat_out(1, 3) = -mat_in(0, 1)*mat_in(0, 3) - mat_in(1, 1)*mat_in(1, 3) - mat_in(2, 1)*mat_in(2, 3);
	mat_out(2, 3) = -mat_in(0, 2)*mat_in(0, 3) - mat_in(1, 2)*mat_in(1, 3) - mat_in(2, 2)*mat_in(2, 3);
	mat_out(3, 0) = 0.; mat_out(3, 1) = 0.; mat_out(3, 2) = 0.;
}

void TruncatedTspline_3D::FindLocalKnotVector(int pid, Matrix4d lcs, array<double, 5>& ku, array<double, 5>& kv, array<double, 5>& kw)
{
	//array<double, 5> kutmp = { -cp[pid].kitvU[0] - cp[pid].kitvU[1], -cp[pid].kitvU[1], 0., cp[pid].kitvU[2], cp[pid].kitvU[2] + cp[pid].kitvU[3] };
	//array<double, 5> kvtmp = { -cp[pid].kitvV[0] - cp[pid].kitvV[1], -cp[pid].kitvV[1], 0., cp[pid].kitvV[2], cp[pid].kitvV[2] + cp[pid].kitvV[3] };
	//array<double, 5> kwtmp = { -cp[pid].kitvW[0] - cp[pid].kitvW[1], -cp[pid].kitvW[1], 0., cp[pid].kitvW[2], cp[pid].kitvW[2] + cp[pid].kitvW[3] };
	//cout << "before\n";
	//cout << kutmp[0] << " " << kutmp[1] << " " << kutmp[2] << " " << kutmp[3] << " " << kutmp[4] << "\n";
	//cout << kvtmp[0] << " " << kvtmp[1] << " " << kvtmp[2] << " " << kvtmp[3] << " " << kvtmp[4] << "\n";
	//cout << kwtmp[0] << " " << kwtmp[1] << " " << kwtmp[2] << " " << kwtmp[3] << " " << kwtmp[4] << "\n";
	//for (int i = 0; i < 5; i++)
	//{
	//	Vector4d tmp1(kutmp[i], kvtmp[i], kwtmp[i], 1.);
	//	Vector4d tmp2 = lcs*tmp1;
	//	ku[i] = tmp2(0); kv[i] = tmp2(1); kw[i] = tmp2(2);
	//}
	double kvit_pos[3][2] = { { cp[pid].kitvU[2], cp[pid].kitvU[2] + cp[pid].kitvU[3] }, { cp[pid].kitvV[2], cp[pid].kitvV[2] + cp[pid].kitvV[3] }, 
	{ cp[pid].kitvW[2], cp[pid].kitvW[2] + cp[pid].kitvW[3] } };
	double kvit_neg[3][2] = { { cp[pid].kitvU[1], cp[pid].kitvU[0] + cp[pid].kitvU[1] }, { cp[pid].kitvV[1], cp[pid].kitvV[0] + cp[pid].kitvV[1] },
	{ cp[pid].kitvW[1], cp[pid].kitvW[0] + cp[pid].kitvW[1] } };
	vector<Vector3d> v0(3);//local
	vector<Vector3d> v1(3);//global, correspond to ku, kv, kw
	v0[0] << lcs(0, 0), lcs(1, 0), lcs(2, 0);
	v0[1] << lcs(0, 1), lcs(1, 1), lcs(2, 1);
	v0[2] << lcs(0, 2), lcs(1, 2), lcs(2, 2);
	v1[0] << 1, 0, 0;
	v1[1] << 0, 1, 0;
	v1[2] << 0, 0, 1;
	int j;
	for (j = 0; j < 3; j++)
	{
		if (v1[0].dot(v0[j]) == 1.)
		{
			ku[0] = lcs(0, 3) - kvit_neg[j][1];
			ku[1] = lcs(0, 3) - kvit_neg[j][0];
			ku[2] = lcs(0, 3);
			ku[3] = lcs(0, 3) + kvit_pos[j][0];
			ku[4] = lcs(0, 3) + kvit_pos[j][1];
		}
		else if (v1[0].dot(v0[j]) == -1.)
		{
			ku[0] = lcs(0, 3) - kvit_pos[j][1];
			ku[1] = lcs(0, 3) - kvit_pos[j][0];
			ku[2] = lcs(0, 3);
			ku[3] = lcs(0, 3) + kvit_neg[j][0];
			ku[4] = lcs(0, 3) + kvit_neg[j][1];
		}
	}
	for (j = 0; j < 3; j++)
	{
		if (v1[1].dot(v0[j]) == 1.)
		{
			kv[0] = lcs(1, 3) - kvit_neg[j][1];
			kv[1] = lcs(1, 3) - kvit_neg[j][0];
			kv[2] = lcs(1, 3);
			kv[3] = lcs(1, 3) + kvit_pos[j][0];
			kv[4] = lcs(1, 3) + kvit_pos[j][1];
		}
		else if (v1[1].dot(v0[j]) == -1.)
		{
			kv[0] = lcs(1, 3) - kvit_pos[j][1];
			kv[1] = lcs(1, 3) - kvit_pos[j][0];
			kv[2] = lcs(1, 3);
			kv[3] = lcs(1, 3) + kvit_neg[j][0];
			kv[4] = lcs(1, 3) + kvit_neg[j][1];
		}
	}
	for (j = 0; j < 3; j++)
	{
		if (v1[2].dot(v0[j]) == 1.)
		{
			kw[0] = lcs(2, 3) - kvit_neg[j][1];
			kw[1] = lcs(2, 3) - kvit_neg[j][0];
			kw[2] = lcs(2, 3);
			kw[3] = lcs(2, 3) + kvit_pos[j][0];
			kw[4] = lcs(2, 3) + kvit_pos[j][1];
		}
		else if (v1[2].dot(v0[j]) == -1.)
		{
			kw[0] = lcs(2, 3) - kvit_pos[j][1];
			kw[1] = lcs(2, 3) - kvit_pos[j][0];
			kw[2] = lcs(2, 3);
			kw[3] = lcs(2, 3) + kvit_neg[j][0];
			kw[4] = lcs(2, 3) + kvit_neg[j][1];
		}
	}


	//for (int i = 0; i < 5; i++)
	//{
	//	Vector4d tmp1(kutmp[i],0.,0.,1.);
	//	Vector4d tmp2 = lcs*tmp1;
	//	ku[i] = tmp2(0);
	//}
	//for (int i = 0; i < 5; i++)
	//{
	//	Vector4d tmp1(0., kvtmp[i], 0., 1.);
	//	Vector4d tmp2 = lcs*tmp1;
	//	kv[i] = tmp2(1);
	//}
	//for (int i = 0; i < 5; i++)
	//{
	//	Vector4d tmp1(0., 0., kwtmp[i], 1.);
	//	Vector4d tmp2 = lcs*tmp1;
	//	kw[i] = tmp2(2);
	//}
}

bool TruncatedTspline_3D::CheckSupport(const array<double, 2>& u, const array<double, 2>& v, const array<double, 2>& w, const array<double, 5>& ku, const array<double, 5>& kv, const array<double, 5>& kw)
{
	if (ku[0]<u[1] && ku[4]>u[0] && kv[0]<v[1] && kv[4]>v[0] && kw[0]<w[1] && kw[4]>w[0])
	{
		return true;
	}
	else
	{
		return false;
	}
}

void TruncatedTspline_3D::Update_IEN()
{
	for (uint eid = 0; eid<tmesh.size(); eid++)
	{
		vector<int>().swap(tmesh[eid].IEN);
		vector<array<double, 5>>().swap(tmesh[eid].patch_ku);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kv);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kw);
		tmesh[eid].IEN = tmesh[eid].IENtmp;
		tmesh[eid].patch_ku = tmesh[eid].patch_kutmp;
		tmesh[eid].patch_kv = tmesh[eid].patch_kvtmp;
		tmesh[eid].patch_kw = tmesh[eid].patch_kwtmp;
		vector<int>().swap(tmesh[eid].IENtmp);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kutmp);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kvtmp);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kwtmp);
		//if (tmesh[eid].type == 0)
		//{
		//	cout << "element id: " << eid << "\n";
		//	cout << "#IEN: " << tmesh[eid].IEN.size() << "\n";
		//	cout << "#patch_ku: " << tmesh[eid].patch_ku.size() << "\n";
		//	cout << "#patch_kv: " << tmesh[eid].patch_kv.size() << "\n";
		//	cout << "#patch_kw: " << tmesh[eid].patch_kw.size() << "\n";
		//	//for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
		//	//{
		//	//	cout << "point id: " << tmesh[eid].IEN[i] << "\n";
		//	//	for (uint j = 0; j < tmesh[eid].patch_ku[i].size(); j++)
		//	//	{
		//	//		cout << tmesh[eid].patch_ku[i][j] << " ";
		//	//	}
		//	//	cout << "\n\n";
		//	//}
		//	getchar();
		//}
	}
}

//void TruncatedTspline_3D::UpdatePatchCP_Unstruct(int eid)
//{
//	vector<int> IEN_new_all;
//	vector<array<double,5>> KU_new_all;
//	vector<array<double,5>> KV_new_all;
//	for(uint i=0; i<tmesh[eid].chd.size(); i++)
//	{
//		int cid(tmesh[eid].chd[i]);
//		for(uint j=0; j<tmesh[cid].IENtmp.size(); j++)
//		{
//			vector<int>::iterator it=find(IEN_new_all.begin(),IEN_new_all.end(),tmesh[cid].IENtmp[j]);
//			if(it==IEN_new_all.end())
//			{
//				IEN_new_all.push_back(tmesh[cid].IENtmp[j]);
//				array<double,5> kutmp,kvtmp;
//				for(int k=0; k<5; k++)
//				{
//					kutmp[k]=tmesh[cid].patch_kutmp[j][k]+tmesh[eid].chd_o[i][0];
//					kvtmp[k]=tmesh[cid].patch_kvtmp[j][k]+tmesh[eid].chd_o[i][1];
//				}
//				KU_new_all.push_back(kutmp);
//				KV_new_all.push_back(kvtmp);
//			}
//		}
//	}
//
//	vector<int> IEN_old(tmesh[eid].IEN);
//	vector<vector<double>> cmat(IEN_new_all.size(),vector<double>(IEN_old.size(),0.));
//	//vector<vector<double>> tmat(IEN_new_all.size(),vector<double>(IEN_old.size(),0.));
//	for(uint i=0; i<IEN_old.size(); i++)
//	{
//		vector<int>::iterator it=find(IEN_new_all.begin(),IEN_new_all.end(),IEN_old[i]);
//		if(it!=IEN_new_all.end())
//		{
//			int loc(it-IEN_new_all.begin());
//			if(tmesh[eid].patch_ku[i]!=KU_new_all[loc] || tmesh[eid].patch_kv[i]!=KV_new_all[loc])
//			{
//				cp[IEN_old[i]].aff=1;
//			}
//			else
//			{
//				cmat[loc][i]=1.;
//			}
//		}
//	}
//
//	for(uint i=0; i<IEN_new_all.size(); i++)
//	{
//		for(uint j=0; j<IEN_old.size(); j++)
//		{
//			if(CheckSubKnotVector(KU_new_all[i],KV_new_all[i],tmesh[eid].patch_ku[j],tmesh[eid].patch_kv[j]))
//			{
//				vector<double> ku1(10),kv1(10);
//				vector<vector<double>> Tu,Tv;
//				vector<double> ku(tmesh[eid].patch_ku[j].begin(),tmesh[eid].patch_ku[j].end());
//				vector<double> kv(tmesh[eid].patch_kv[j].begin(),tmesh[eid].patch_kv[j].end());
//				vector<double>::iterator it1,it2;
//				it1=set_union(tmesh[eid].patch_ku[j].begin(),tmesh[eid].patch_ku[j].end(),KU_new_all[i].begin(),KU_new_all[i].end(),ku1.begin());
//				it2=set_union(tmesh[eid].patch_kv[j].begin(),tmesh[eid].patch_kv[j].end(),KV_new_all[i].begin(),KV_new_all[i].end(),kv1.begin());
//				ku1.resize(it1-ku1.begin());
//				kv1.resize(it2-kv1.begin());
//				TMatrix(ku,ku1,3,Tu);
//				TMatrix(kv,kv1,3,Tv);
//				it1=search(ku1.begin(),ku1.end(),KU_new_all[i].begin(),KU_new_all[i].end());
//				it2=search(kv1.begin(),kv1.end(),KV_new_all[i].begin(),KV_new_all[i].end());
//				if(it1!=ku1.end() && it2!=kv1.end())
//				{
//					int loc1=it1-ku1.begin();
//					int loc2=it2-kv1.begin();
//					cmat[i][j]=Tu[loc1][0]*Tv[loc2][0];
//				}
//			}
//		}
//	}
//	for(uint i=0; i<IEN_new_all.size(); i++)
//	{
//		//cp[IEN_new_all[i]].coortmp[0]=0.; cp[IEN_new_all[i]].coortmp[1]=0.; cp[IEN_new_all[i]].coortmp[2]=0.;
//		if(cp[IEN_new_all[i]].update==0)
//		{
//			for(uint j=0; j<IEN_old.size(); j++)
//			{
//				if(cmat[i][j]!=0.)
//				{
//					cp[IEN_new_all[i]].coortmp[0]+=cmat[i][j]*cp[IEN_old[j]].coor[0];
//					cp[IEN_new_all[i]].coortmp[1]+=cmat[i][j]*cp[IEN_old[j]].coor[1];
//					cp[IEN_new_all[i]].coortmp[2]+=cmat[i][j]*cp[IEN_old[j]].coor[2];
//				}
//			}
//			cp[IEN_new_all[i]].update=1;
//		}
//	}
//	//check truncation
//	/*double sp0[4][2]={{0.,0.},{tmedge[tmesh[eid].edge[0]].len,0.},{tmedge[tmesh[eid].edge[0]].len,tmedge[tmesh[eid].edge[1]].len},{0.,tmedge[tmesh[eid].edge[1]].len}};
//	vector<array<double,2>> spt(9);
//	for(int i=0; i<4; i++)
//	{
//		spt[i][0]=sp0[i][0]; spt[i][1]=sp0[i][1];
//		spt[i+4][0]=(sp0[i][0]+sp0[(i+1)%4][0])/2.; spt[i+4][1]=(sp0[i][1]+sp0[(i+1)%4][1])/2.; 
//	}
//	spt[8][0]=spt[4][0]; spt[8][1]=spt[5][1];
//	for(uint i=0; i<IEN_old.size(); i++)
//	{
//		vector<double> coef;
//		for(uint j=0; j<IEN_new_all.size(); j++)
//		{
//			coef.push_back(cmat[j][i]);
//		}
//		if(!CheckFullChildren_1(spt,tmesh[eid].patch_ku[i],tmesh[eid].patch_kv[i],KU_new_all,KV_new_all,coef))
//		{
//			cp[IEN_old[i]].truntmp=1;
//			for(uint j=0; j<IEN_new_all.size(); j++)
//			{
//				if(cmat[j][i]!=0. && IEN_old[i]!=IEN_new_all[j])
//				{
//					vector<int>::iterator it=find(cp[IEN_old[i]].tbftmp.begin(),cp[IEN_old[i]].tbftmp.end(),IEN_new_all[j]);
//					if(it==cp[IEN_old[i]].tbftmp.end())
//					{
//						cp[IEN_old[i]].tbftmp.push_back(IEN_new_all[j]);
//						cp[IEN_old[i]].tctmp.push_back(cmat[j][i]);
//					}
//				}
//			}
//		}
//	}*/
//}

//void TruncatedTspline_3D::Truncation()
//{
//	for(uint eid=0; eid<tmesh.size(); eid++)
//	{
//		if(tmesh[eid].act==1 && (tmesh[eid].type==0 || tmesh[eid].type==1))
//		{
//			for(uint i=0; i<tmesh[eid].IEN.size(); i++)
//			{
//				for(uint j=i+1; j<tmesh[eid].IEN.size(); j++)
//				{
//					if(CheckSubKnotVector(tmesh[eid].patch_ku[j],tmesh[eid].patch_kv[j],tmesh[eid].patch_ku[i],tmesh[eid].patch_kv[i]))//i to be truncated, j to be discarded children
//					{
//						int Bt(tmesh[eid].IEN[i]), Bc(tmesh[eid].IEN[j]);
//						vector<int>::iterator it=find(cp[Bt].tbf.begin(),cp[Bt].tbf.end(),Bc);
//						if(it==cp[Bt].tbf.end())
//						{
//							vector<double> ku1(10),kv1(10);
//							vector<vector<double>> Tu,Tv;
//							vector<double> ku(tmesh[eid].patch_ku[i].begin(),tmesh[eid].patch_ku[i].end());
//							vector<double> kv(tmesh[eid].patch_kv[i].begin(),tmesh[eid].patch_kv[i].end());
//							vector<double>::iterator it1,it2;
//							it1=set_union(tmesh[eid].patch_ku[i].begin(),tmesh[eid].patch_ku[i].end(),tmesh[eid].patch_ku[j].begin(),tmesh[eid].patch_ku[j].end(),ku1.begin());
//							it2=set_union(tmesh[eid].patch_kv[i].begin(),tmesh[eid].patch_kv[i].end(),tmesh[eid].patch_kv[j].begin(),tmesh[eid].patch_kv[j].end(),kv1.begin());
//							ku1.resize(it1-ku1.begin());
//							kv1.resize(it2-kv1.begin());
//							TMatrix(ku,ku1,3,Tu);
//							TMatrix(kv,kv1,3,Tv);
//							it1=search(ku1.begin(),ku1.end(),tmesh[eid].patch_ku[j].begin(),tmesh[eid].patch_ku[j].end());
//							it2=search(kv1.begin(),kv1.end(),tmesh[eid].patch_kv[j].begin(),tmesh[eid].patch_kv[j].end());
//							if(it1!=ku1.end() && it2!=kv1.end())
//							{
//								int loc1=it1-ku1.begin();
//								int loc2=it2-kv1.begin();
//								cp[Bt].trun=1;
//								cp[Bt].tbf.push_back(Bc);
//								cp[Bt].tc.push_back(Tu[loc1][0]*Tv[loc2][0]);
//							}
//						}
//					}
//					else if(CheckSubKnotVector(tmesh[eid].patch_ku[i],tmesh[eid].patch_kv[i],tmesh[eid].patch_ku[j],tmesh[eid].patch_kv[j]))//j to be truncated, i to be discarded children
//					{
//						int Bt(tmesh[eid].IEN[j]), Bc(tmesh[eid].IEN[i]);
//						vector<int>::iterator it=find(cp[Bt].tbf.begin(),cp[Bt].tbf.end(),Bc);
//						if(it==cp[Bt].tbf.end())
//						{
//							vector<double> ku1(10),kv1(10);
//							vector<vector<double>> Tu,Tv;
//							vector<double> ku(tmesh[eid].patch_ku[j].begin(),tmesh[eid].patch_ku[j].end());
//							vector<double> kv(tmesh[eid].patch_kv[j].begin(),tmesh[eid].patch_kv[j].end());
//							vector<double>::iterator it1,it2;
//							it1=set_union(tmesh[eid].patch_ku[j].begin(),tmesh[eid].patch_ku[j].end(),tmesh[eid].patch_ku[i].begin(),tmesh[eid].patch_ku[i].end(),ku1.begin());
//							it2=set_union(tmesh[eid].patch_kv[j].begin(),tmesh[eid].patch_kv[j].end(),tmesh[eid].patch_kv[i].begin(),tmesh[eid].patch_kv[i].end(),kv1.begin());
//							ku1.resize(it1-ku1.begin());
//							kv1.resize(it2-kv1.begin());
//							TMatrix(ku,ku1,3,Tu);
//							TMatrix(kv,kv1,3,Tv);
//							it1=search(ku1.begin(),ku1.end(),tmesh[eid].patch_ku[i].begin(),tmesh[eid].patch_ku[i].end());
//							it2=search(kv1.begin(),kv1.end(),tmesh[eid].patch_kv[i].begin(),tmesh[eid].patch_kv[i].end());
//							if(it1!=ku1.end() && it2!=kv1.end())
//							{
//								int loc1=it1-ku1.begin();
//								int loc2=it2-kv1.begin();
//								cp[Bt].trun=1;
//								cp[Bt].tbf.push_back(Bc);
//								cp[Bt].tc.push_back(Tu[loc1][0]*Tv[loc2][0]);
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//}
//
//void TruncatedTspline_3D::ElementBasis(int eid, double u, double v, vector<double>& Nt, vector<array<double,2>>& dNdt)
//{
//	if(tmesh[eid].act==1)
//	{
//		if(tmesh[eid].type==0 || tmesh[eid].type==1)
//		{
//			ElementBasis_Regular(eid,u,v,Nt,dNdt);
//		}
//		else if(tmesh[eid].type==4)
//		{
//			ElementBasis_Irregular(eid,u,v,Nt,dNdt);
//		}
//	}
//}
//
//void TruncatedTspline_3D::ElementBasis_Regular(int eid, double u, double v, vector<double>& Nt, vector<array<double,2>>& dNdt)
//{
//	Nt.clear();
//	dNdt.clear();
//	Nt.resize(tmesh[eid].IEN.size());
//	dNdt.resize(tmesh[eid].IEN.size());
//	vector<double> Nt0(tmesh[eid].IEN.size());
//	vector<array<double,2>> dNdt0(tmesh[eid].IEN.size());
//	vector<double> ku(5,0.),kv(5,0.),uval,vval;
//	BSplineBasis bu,bv;
//	for(uint i=0;i<tmesh[eid].IEN.size();i++)
//	{
//		ku.assign(tmesh[eid].patch_ku[i].begin(),tmesh[eid].patch_ku[i].end());
//		kv.assign(tmesh[eid].patch_kv[i].begin(),tmesh[eid].patch_kv[i].end());
//		bu.Set(3,ku);
//		bv.Set(3,kv);
//		bu.BasisFunction(0,u,1,uval);
//		bv.BasisFunction(0,v,1,vval);
//		Nt0[i]=uval[0]*vval[0];
//		dNdt0[i][0]=uval[1]*vval[0];
//		dNdt0[i][1]=uval[0]*vval[1];
//	}
//	for(uint i=0;i<tmesh[eid].IEN.size();i++)
//	{
//		Nt[i]=Nt0[i];
//		dNdt[i][0]=dNdt0[i][0];
//		dNdt[i][1]=dNdt0[i][1];
//		if(cp[tmesh[eid].IEN[i]].trun==1)
//		{
//			int pid(tmesh[eid].IEN[i]);
//			for(uint j=0; j<cp[pid].tbf.size(); j++)
//			{
//				vector<int>::iterator it=find(tmesh[eid].IEN.begin(),tmesh[eid].IEN.end(),cp[pid].tbf[j]);
//				if(it!=tmesh[eid].IEN.end())
//				{
//					int loc(it-tmesh[eid].IEN.begin());
//					Nt[i]-=cp[pid].tc[j]*Nt0[loc];
//					dNdt[i][0]-=cp[pid].tc[j]*dNdt0[loc][0];
//					dNdt[i][1]-=cp[pid].tc[j]*dNdt0[loc][1];
//				}
//			}
//		}
//	}
//}
//
//void TruncatedTspline_3D::ElementBasis_Irregular(int eid, double u, double v, vector<double>& Nt, vector<array<double,2>>& dNdt)
//{
//	//need to change variabe for Bezier, (u,v) -> (u_b,v_b)
//	Nt.clear();
//	dNdt.clear();
//	Nt.resize(tmesh[eid].IEN.size(),0.);
//	dNdt.resize(tmesh[eid].IEN.size());
//	vector<double> Nt1(tmesh[eid].IEN.size(),0.);
//	vector<array<double,2>> dNdt1(tmesh[eid].IEN.size());
//	uint nv(cp[tmesh[eid].cnct[0]].face.size());
//	vector<vector<double>> bmat(tmesh[eid].bemat.size(),vector<double>(tmesh[eid].bemat[0].size()));
//	for(uint i=0; i<tmesh[eid].bemat.size(); i++)
//	{
//		for(uint j=0; j<tmesh[eid].bemat[i].size(); j++)
//		{
//			bmat[i][j]=tmesh[eid].bemat[i][j];
//		}
//	}
//	//SetBezier4TranMatOP(nv,bmat);
//	BezierElement2D be(4);
//	vector<double> Nt0;
//	vector<array<double,2>> dNdt0;
//	double u_b(u/tmedge[tmesh[eid].edge[0]].len), v_b(v/tmedge[tmesh[eid].edge[3]].len);
//	be.Basis(u_b,v_b,Nt0,dNdt0);
//	for(uint i=0; i<2*nv+1; i++)
//	{
//		dNdt1[i][0]=0.; dNdt1[i][1]=0.;
//		for(int j=0; j<25; j++)
//		{
//			Nt1[i]+=bmat[i][j]*Nt0[j];
//			dNdt1[i][0]+=bmat[i][j]*dNdt0[j][0];
//			dNdt1[i][1]+=bmat[i][j]*dNdt0[j][1];
//		}
//		dNdt1[i][0]/=tmedge[tmesh[eid].edge[0]].len;
//		dNdt1[i][1]/=tmedge[tmesh[eid].edge[3]].len;
//	}
//	vector<double> ku(5,0.),kv(5,0.),uval,vval;
//	BSplineBasis bu,bv;
//	for(uint i=2*nv+1; i<tmesh[eid].IEN.size(); i++)
//	{
//		ku.assign(tmesh[eid].patch_ku[i-(2*nv+1)].begin(),tmesh[eid].patch_ku[i-(2*nv+1)].end());
//		kv.assign(tmesh[eid].patch_kv[i-(2*nv+1)].begin(),tmesh[eid].patch_kv[i-(2*nv+1)].end());
//		bu.Set(3,ku);
//		bv.Set(3,kv);
//		bu.BasisFunction(0,u,1,uval);
//		bv.BasisFunction(0,v,1,vval);
//		Nt1[i]=uval[0]*vval[0];
//		dNdt1[i][0]=uval[1]*vval[0];
//		dNdt1[i][1]=uval[0]*vval[1];
//	}
//	for(uint i=0;i<tmesh[eid].IEN.size();i++)
//	{
//		Nt[i]=Nt1[i];
//		dNdt[i][0]=dNdt1[i][0];
//		dNdt[i][1]=dNdt1[i][1];
//		if(cp[tmesh[eid].IEN[i]].trun==1)
//		{
//			int pid(tmesh[eid].IEN[i]);
//			for(uint j=0; j<cp[pid].tbf.size(); j++)
//			{
//				vector<int>::iterator it=find(tmesh[eid].IEN.begin(),tmesh[eid].IEN.end(),cp[pid].tbf[j]);
//				if(it!=tmesh[eid].IEN.end())
//				{
//					int loc(it-tmesh[eid].IEN.begin());
//					Nt[i]-=cp[pid].tc[j]*Nt1[loc];
//					dNdt[i][0]-=cp[pid].tc[j]*dNdt1[loc][0];
//					dNdt[i][1]-=cp[pid].tc[j]*dNdt1[loc][1];
//				}
//			}
//		}
//	}
//}

void TruncatedTspline_3D::SetBezierMatIrrPatch(int eid)
{
	//find IENtmp first
	tmesh[eid].IENtmp.clear();
	uint i, j, k, hxid;
	vector<int> loc(cp.size(),-1);
	vector<int> hx1r(1, eid);
	for (i = 0; i < 8; i++)
	{
		loc[tmesh[eid].cnct[i]] = tmesh[eid].IENtmp.size();
		tmesh[eid].IENtmp.push_back(tmesh[eid].cnct[i]);
	}
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < cp[tmesh[eid].cnct[i]].hex.size(); j++)
		{
			hxid = cp[tmesh[eid].cnct[i]].hex[j];
			vector<int>::iterator it1 = find(hx1r.begin(),hx1r.end(),hxid);
			if (it1 == hx1r.end())
			{
				hx1r.push_back(hxid);
				for (k = 0; k < 8; k++)
				{
					vector<int>::iterator it = find(tmesh[eid].IENtmp.begin(), tmesh[eid].IENtmp.end(), tmesh[hxid].cnct[k]);
					if (it == tmesh[eid].IENtmp.end())
					{
						loc[tmesh[hxid].cnct[k]] = tmesh[eid].IENtmp.size();
						tmesh[eid].IENtmp.push_back(tmesh[hxid].cnct[k]);
					}
				}
			}
		}
	}
	for (i = 0; i < tmesh[eid].bemat.size(); i++)
	{
		tmesh[eid].bemat[i].clear();
	}
	tmesh[eid].bemat.clear();
	tmesh[eid].bemat.resize(tmesh[eid].IENtmp.size(),vector<double>(64,0.));
	//8 body points, not consider boundary yet
	double w[2] = {2./3.,1./3.};
	double a[8] = { w[0] * w[0] * w[0], w[1] * w[0] * w[0], w[1] * w[1] * w[0], w[1] * w[0] * w[0],
		w[1] * w[0] * w[0], w[1] * w[1] * w[0], w[1] * w[1] * w[1], w[1] * w[1] * w[0] };
	int bpi[8][8] = { { 0, 1, 2, 3, 4, 5, 6, 7 }, { 1, 2, 3, 0, 5, 6, 7, 4 }, { 2, 3, 0, 1, 6, 7, 4, 5 }, { 3, 0, 1, 2, 7, 4, 5,6 },
	{ 4, 5, 6, 7, 0, 1, 2, 3 }, { 5, 6, 7, 4, 1, 2, 3, 0 }, { 6, 7, 4, 5, 2, 3, 0, 1 }, { 7, 4, 5, 6, 3, 0, 1, 2 } };
	vector<array<array<double, 8>, 8>> bpm(hx1r.size());
	vector<array<array<int, 8>, 8>> bpmap(hx1r.size());
	for (i = 0; i < hx1r.size(); i++)//which element
	{
		for (j = 0; j < 8; j++)//which body point, bezier
		{
			for (k = 0; k < 8; k++)//which local corner point, b-splines
			{
				bpm[i][j][k] = a[k];
				bpmap[i][j][k] = loc[tmesh[hx1r[i]].cnct[bpi[j][k]]];
			}
		}
	}
	int layer[4] = {0,16,32,48};
	int bpbz[8] = { 5 + layer[1], 6 + layer[1], 10 + layer[1], 9 + layer[1], 5 + layer[2], 6 + layer[2], 10 + layer[2], 9 + layer[2] };
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 8; j++)
		{
			tmesh[eid].bemat[bpmap[0][i][j]][bpbz[i]] = bpm[0][i][j];
		}
	}
	//2*12 edge points
	int edi[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, {7,4} };
	int edbz[12][2] = { { 1, 2 }, { 7, 11 }, { 14, 13 }, { 8, 4 }, { 0 + layer[1], 0 + layer[2] }, { 3 + layer[1], 3 + layer[2] }, 
	{ 15 + layer[1], 15 + layer[2] }, { 12 + layer[1], 12 + layer[2] }, { 1 + layer[3], 2 + layer[3] }, { 7 + layer[3], 11 + layer[3] }, { 14 + layer[3], 13 + layer[3] }, { 8 + layer[3], 4 + layer[3] } };
	int pos1,pos2;
	for (i = 0; i < 12; i++)
	{
		uint nhex = tmedge[tmesh[eid].edge[i]].hex.size();
		for (j = 0; j<tmedge[tmesh[eid].edge[i]].hex.size(); j++)
		{
			hxid = tmedge[tmesh[eid].edge[i]].hex[j];
			vector<int>::iterator it = find(hx1r.begin(),hx1r.end(),hxid);
			pos1 = it - hx1r.begin();
			int* it1 = find(tmesh[hxid].cnct, tmesh[hxid].cnct+8, tmesh[eid].cnct[edi[i][0]]);
			pos2 = it1 - tmesh[hxid].cnct;
			for (k = 0; k < 8; k++)
			{
				tmesh[eid].bemat[bpmap[pos1][pos2][k]][edbz[i][0]] += bpm[pos1][pos2][k]/nhex;
			}
			int* it2 = find(tmesh[hxid].cnct, tmesh[hxid].cnct + 8, tmesh[eid].cnct[edi[i][1]]);
			pos2 = it2 - tmesh[hxid].cnct;
			for (k = 0; k < 8; k++)
			{
				tmesh[eid].bemat[bpmap[pos1][pos2][k]][edbz[i][1]] += bpm[pos1][pos2][k] / nhex;
			}
		}
	}
	//8 corner points
	int cnbz[8] = { 0, 3, 15, 12, 0 + layer[3], 3 + layer[3], 15 + layer[3], 12 + layer[3]};
	for (i = 0; i < 8; i++)
	{
		uint nhex = cp[tmesh[eid].cnct[i]].hex.size();
		for (j = 0; j<nhex; j++)
		{
			hxid = cp[tmesh[eid].cnct[i]].hex[j];
			vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
			pos1 = it - hx1r.begin();
			int* it1 = find(tmesh[hxid].cnct, tmesh[hxid].cnct + 8, tmesh[eid].cnct[i]);
			pos2 = it1 - tmesh[hxid].cnct;
			for (k = 0; k < 8; k++)
			{
				tmesh[eid].bemat[bpmap[pos1][pos2][k]][cnbz[i]] += bpm[pos1][pos2][k] / nhex;
			}
		}
	}
	//4*6 face points
	int fci[6][4] = { { 0, 1, 3, 2 }, { 0, 1, 4, 5 }, { 1, 2, 5, 6 }, { 3, 2, 7, 6 }, { 0, 3, 4, 7 }, {4,5,7,6} };
	int fcbz[6][4] = { { 5, 6, 9, 10 }, { 1 + layer[1], 2 + layer[1], 1 + layer[2], 2 + layer[2] }, { 7 + layer[1], 11 + layer[1], 7 + layer[2], 11 + layer[2] },
	{ 13 + layer[1], 14 + layer[1], 13 + layer[2], 14 + layer[2] }, { 4 + layer[1], 8 + layer[1], 4 + layer[2], 8 + layer[2] }, { 5 + layer[3], 6 + layer[3], 9 + layer[3], 10 + layer[3] } };
	for (i = 0; i < 6; i++)
	{
		uint nhex = tmface[tmesh[eid].face[i]].hex.size();
		for (j = 0; j < nhex; j++)
		{
			hxid = tmface[tmesh[eid].face[i]].hex[j];
			vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
			pos1 = it - hx1r.begin();
			for (int j1 = 0; j1 < 4; j1++)
			{
				int* it1 = find(tmesh[hxid].cnct, tmesh[hxid].cnct + 8, tmesh[eid].cnct[fci[i][j1]]);
				pos2 = it1 - tmesh[hxid].cnct;
				for (k = 0; k < 8; k++)
				{
					tmesh[eid].bemat[bpmap[pos1][pos2][k]][fcbz[i][j1]] += bpm[pos1][pos2][k] / nhex;
				}
			}
		}
	}

	//cout << "bemat of eid: " << eid << "\n";
	//vector<double> col_sum(tmesh[eid].bemat[0].size(), 0.);
	//for (i = 0; i < tmesh[eid].bemat.size(); i++)
	//{
	//	for (j = 0; j < tmesh[eid].bemat[i].size(); j++)
	//	{
	//		col_sum[j] += tmesh[eid].bemat[i][j];
	//		//cout << tmesh[eid].bemat[i][j] << " ";
	//		//if (j % 16 == 0) cout << "\n";
	//	}
	//	//cout << "\n";
	//}
	//for (i = 0; i < col_sum.size(); i++)
	//{
	//	cout << col_sum[i] << "\n";
	//}
	//getchar();
}

void TruncatedTspline_3D::AllBezierPatch()//test purpose
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].type != 1)
		{
			SetBezierMatIrrPatch(eid);
		}
	}
}

void TruncatedTspline_3D::ElementBasis(int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt)
{
	if (tmesh[eid].act == 1)
	{
		if (tmesh[eid].type == 0)
		{
			ElementBasis_Regular(eid, u, Nt, dNdt);
		}
		else if (tmesh[eid].type == 2)
		{
			ElementBasis_Irregular(eid, u, Nt, dNdt);
		}
	}
}

void TruncatedTspline_3D::ElementBasis_Regular(int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt)
{
	Nt.clear();
	dNdt.clear();
	Nt.resize(tmesh[eid].IEN.size());
	dNdt.resize(tmesh[eid].IEN.size());
	vector<double> Nt0(tmesh[eid].IEN.size());
	vector<array<double, 3>> dNdt0(tmesh[eid].IEN.size());
	vector<double> ku(5, 0.), kv(5, 0.), kw(5,0.), uval, vval, wval;
	BSplineBasis bu, bv, bw;
	for (uint i = 0; i<tmesh[eid].IEN.size(); i++)
	{
		ku.assign(tmesh[eid].patch_ku[i].begin(), tmesh[eid].patch_ku[i].end());
		kv.assign(tmesh[eid].patch_kv[i].begin(), tmesh[eid].patch_kv[i].end());
		kw.assign(tmesh[eid].patch_kw[i].begin(), tmesh[eid].patch_kw[i].end());
		bu.Set(3, ku);
		bv.Set(3, kv);
		bw.Set(3, kw);
		bu.BasisFunction(0, u[0], 1, uval);
		bv.BasisFunction(0, u[1], 1, vval);
		bw.BasisFunction(0, u[2], 1, wval);
		Nt0[i] = uval[0] * vval[0] * wval[0];
		dNdt0[i][0] = uval[1] * vval[0] * wval[0];
		dNdt0[i][1] = uval[0] * vval[1] * wval[0];
		dNdt0[i][2] = uval[0] * vval[0] * wval[1];
	}
	for (uint i = 0; i<tmesh[eid].IEN.size(); i++)
	{
		Nt[i] = Nt0[i];
		dNdt[i][0] = dNdt0[i][0];
		dNdt[i][1] = dNdt0[i][1];
		dNdt[i][2] = dNdt0[i][2];
		if (cp[tmesh[eid].IEN[i]].trun == 1)//later improve by setPatchTruncation()
		{
			int pid(tmesh[eid].IEN[i]);
			for (uint j = 0; j<cp[pid].tbf.size(); j++)
			{
				vector<int>::iterator it = find(tmesh[eid].IEN.begin(), tmesh[eid].IEN.end(), cp[pid].tbf[j]);
				if (it != tmesh[eid].IEN.end())
				{
					int loc(it - tmesh[eid].IEN.begin());
					Nt[i] -= cp[pid].tc[j] * Nt0[loc];
					dNdt[i][0] -= cp[pid].tc[j] * dNdt0[loc][0];
					dNdt[i][1] -= cp[pid].tc[j] * dNdt0[loc][1];
					dNdt[i][2] -= cp[pid].tc[j] * dNdt0[loc][2];
				}
			}
		}
	}
}

void TruncatedTspline_3D::ElementBasis_Irregular(int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt)
{
	Nt.clear();
	dNdt.clear();
	Nt.resize(tmesh[eid].IEN.size());
	dNdt.resize(tmesh[eid].IEN.size());
	vector<double> Nt0(tmesh[eid].IEN.size());
	vector<array<double, 3>> dNdt0(tmesh[eid].IEN.size());
	BezierElement3D bzel;
	vector<double> Bt;
	vector<array<double, 3>> dBdt;
	bzel.Basis(u[0],u[1],u[2],Bt,dBdt);
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		Nt0[i] = 0.;
		dNdt0[i][0] = 0.; dNdt0[i][1] = 0.; dNdt0[i][2] = 0.;
		for (uint j = 0; j < 64; j++)
		{
			Nt0[i] += tmesh[eid].bemat[i][j] * Bt[j];
			dNdt0[i][0] += tmesh[eid].bemat[i][j] * dBdt[j][0];
			dNdt0[i][1] += tmesh[eid].bemat[i][j] * dBdt[j][1];
			dNdt0[i][2] += tmesh[eid].bemat[i][j] * dBdt[j][2];
		}
	}
	for (uint i = 0; i<tmesh[eid].IEN.size(); i++)
	{
		Nt[i] = Nt0[i];
		dNdt[i][0] = dNdt0[i][0];
		dNdt[i][1] = dNdt0[i][1];
		dNdt[i][2] = dNdt0[i][2];
		if (cp[tmesh[eid].IEN[i]].trun == 1)//later improve by setPatchTruncation()
		{
			int pid(tmesh[eid].IEN[i]);
			for (uint j = 0; j<cp[pid].tbf.size(); j++)
			{
				vector<int>::iterator it = find(tmesh[eid].IEN.begin(), tmesh[eid].IEN.end(), cp[pid].tbf[j]);
				if (it != tmesh[eid].IEN.end())
				{
					int loc(it - tmesh[eid].IEN.begin());
					Nt[i] -= cp[pid].tc[j] * Nt0[loc];
					dNdt[i][0] -= cp[pid].tc[j] * dNdt0[loc][0];
					dNdt[i][1] -= cp[pid].tc[j] * dNdt0[loc][1];
					dNdt[i][2] -= cp[pid].tc[j] * dNdt0[loc][2];
				}
			}
		}
	}
}

void TruncatedTspline_3D::Para2Physical(int eid, const array<double, 3>& u, array<double, 3>& pt)
{
	vector<double> Nt;
	vector<array<double, 3>> dNdt;
	ElementBasis(eid,u,Nt,dNdt);
	pt[0] = 0.; pt[1] = 0.; pt[2] = 0.;
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		pt[0] += Nt[i] * cp[tmesh[eid].IEN[i]].coor[0];
		pt[1] += Nt[i] * cp[tmesh[eid].IEN[i]].coor[1];
		pt[2] += Nt[i] * cp[tmesh[eid].IEN[i]].coor[2];
	}
}

void TruncatedTspline_3D::VisualizeSolidVTK(string fn)
{
	vector<array<double,3>> spt;
	vector<array<double,3>> sval;
	vector<double> ssum;
	vector<array<int,8>> sele;
	vector<array<double,3>> lpt;//visulize parameter lines
	vector<array<int,2>> led;//line connectivity
	int ns(2),ne_ref(0),loc0,loc1,loc2;
	//vector<double> su(ns),sv(ns);
	//for(int i=0; i<ns; i++)
	//{
	//	su[i]=i*1./(ns-1);
	//	sv[i]=i*1./(ns-1);
	//}

	for(uint e=0;e<tmesh.size();e++)
	{
		if(tmesh[e].act==1 && (tmesh[e].type==0 || tmesh[e].type==2))
		{
			int loc(0);
			vector<double> su(ns),sv(ns),sw(ns);
			for(int i=0; i<ns; i++)
			{
				su[i]=i*tmedge[tmesh[e].edge[0]].len/(ns-1);
				sv[i]=i*tmedge[tmesh[e].edge[3]].len/(ns-1);
				sw[i] = i*tmedge[tmesh[e].edge[4]].len / (ns - 1);
			}

			for(int a=0;a<ns;a++)
			{
				for(int b=0;b<ns;b++)
				{
					for (int c = 0; c < ns; c++)
					{
						array<double, 3> pt;
						array<double, 3> uval = { su[c], sv[b], sw[a] };
						Para2Physical(e, uval, pt);
						double sumtmp = PartionOfUnity(e, uval);
						spt.push_back(pt);
						ssum.push_back(sumtmp);
						//if(a==0||a==ns-1||b==0||b==ns-1)
						//{
						//	lpt.push_back(pt);
						//}
					}
				}
			}

			for(int a=0;a<ns-1;a++)
			{
				for(int b=0;b<ns-1;b++)
				{
					for (int c = 0; c < ns-1; c++)
					{
						array<int, 8> el;
						el[0] = ne_ref + a*ns*ns + b*ns+c;
						el[1] = ne_ref + a*ns*ns + b*ns + c+1;
						el[2] = ne_ref + a*ns*ns + (b+1)*ns + c+1;
						el[3] = ne_ref + a*ns*ns + (b+1)*ns + c;
						el[4] = ne_ref + (a+1)*ns*ns + b*ns + c;
						el[5] = ne_ref + (a + 1)*ns*ns + b*ns + c + 1;
						el[6] = ne_ref + (a + 1)*ns*ns + (b + 1)*ns + c + 1;
						el[7] = ne_ref + (a + 1)*ns*ns + (b + 1)*ns + c;
						sele.push_back(el);
					}
				}
			}
			//for(int a=0;a<ns-1;a++)
			//{
			//	array<int,2> lc;
			//	lc[0]=ecount*4*(ns-1)+a;
			//	lc[1]=ecount*4*(ns-1)+a+1;
			//	led.push_back(lc);
			//	lc[0]=ecount*4*(ns-1)+3*ns-4+a;
			//	lc[1]=ecount*4*(ns-1)+3*ns-4+a+1;
			//	led.push_back(lc);
			//}
			//for(int a=0;a<ns-2;a++)
			//{
			//	array<int,2> lc;
			//	lc[0]=ecount*4*(ns-1)+ns+2*a;
			//	lc[1]=ecount*4*(ns-1)+ns+2*a+2;
			//	led.push_back(lc);
			//	lc[0]=ecount*4*(ns-1)+ns+2*a-1;
			//	lc[1]=ecount*4*(ns-1)+ns+2*a+1;
			//	led.push_back(lc);
			//}
			//array<int,2> lc1;
			//lc1[0]=ecount*4*(ns-1);
			//lc1[1]=ecount*4*(ns-1)+ns;
			//led.push_back(lc1);
			//lc1[0]=ecount*4*(ns-1)+3*ns-5;
			//lc1[1]=ecount*4*(ns-1)+4*ns-5;
			//led.push_back(lc1);
			
			ne_ref += ns*ns*ns;
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
		fout<<"\nCELLS "<<sele.size()<<" "<<9*sele.size()<<'\n';
		for(uint i=0;i<sele.size();i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout<<"\nCELL_TYPES "<<sele.size()<<'\n';
		for(uint i=0;i<sele.size();i++)
		{
			fout<<"12\n";
		}
		fout<<"\nPOINT_DATA "<<ssum.size()<<"\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		for(uint i=0;i<ssum.size();i++)
		{
			fout<<ssum[i]<<"\n";
		}
		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nNORMALS Normal FLOAT\n";
		//for(uint i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i][0]<<" "<<sval[i][1]<<" "<<sval[i][2]<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout<<"Cannot open "<<fname<<"!\n";
	}

	/*string fname1(fn+"-lines.vtk");
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
	}*/
}

//void TruncatedTspline_3D::FindPatchKnotVector_Irr(int eid, vector<array<double,5>>& patch_ku, vector<array<double,5>>& patch_kv)
//{
//	patch_ku.clear();
//	patch_kv.clear();
//	patch_ku.resize(5);
//	patch_kv.resize(5);
//	//find knot vectors for node 2 to 6
//	if(tmesh[eid].act==1 && tmesh[eid].type==4)
//	{
//		//first for node 2
//		array<double,2> uv_ref0={tmesh[eid].lcs[3].u[0],tmesh[eid].lcs[3].u[1]},uv0;
//		int fcid0(tmedge[tmesh[eid].edge[3]].face[0]), rot0;
//		if(fcid0==eid) fcid0=tmedge[tmesh[eid].edge[3]].face[1];
//		FindRotateAndUVCoor(tmesh[eid].cnct[3],tmesh[eid].lcs[3].rot,uv_ref0,fcid0,tmesh[fcid0].cnct[2],rot0,uv0);
//		FindLocalKnotVector(tmesh[fcid0].cnct[2],rot0,uv0,patch_ku[0],patch_kv[0]);
//		//node 3, 4, 5
//		for(int i=3; i>0; i--)
//		{
//			array<double,2> uv1={tmesh[eid].lcs[i].u[0],tmesh[eid].lcs[i].u[1]};
//			FindLocalKnotVector(tmesh[eid].cnct[i],tmesh[eid].lcs[i].rot,uv1,patch_ku[4-i],patch_kv[4-i]);
//		}
//		//node 6
//		array<double,2> uv_ref2={tmesh[eid].lcs[1].u[0],tmesh[eid].lcs[1].u[1]},uv2;
//		int fcid2(tmedge[tmesh[eid].edge[0]].face[0]), rot2;
//		if(fcid2==eid) fcid2=tmedge[tmesh[eid].edge[0]].face[1];
//		FindRotateAndUVCoor(tmesh[eid].cnct[1],tmesh[eid].lcs[1].rot,uv_ref2,fcid2,tmesh[fcid2].cnct[2],rot2,uv2);
//		FindLocalKnotVector(tmesh[fcid2].cnct[2],rot2,uv2,patch_ku[4],patch_kv[4]);
//	}
//}

//void TruncatedTspline_3D::SurfacePointMap(int eid, double u, double v, array<double,3>& pt, array<double,3>& norm)
//{
//	vector<double> Nt;
//	vector<array<double,2>> dNdt;
//	ElementBasis(eid,u,v,Nt,dNdt);
//	pt[0]=0.; pt[1]=0.; pt[2]=0.;
//	double nmtmp[2][3]={{0.,0.,0.},{0.,0.,0.}};
//	for(uint i=0; i<tmesh[eid].IEN.size(); i++)
//	{
//		pt[0]+=cp[tmesh[eid].IEN[i]].coor[0]*Nt[i];
//		pt[1]+=cp[tmesh[eid].IEN[i]].coor[1]*Nt[i];
//		pt[2]+=cp[tmesh[eid].IEN[i]].coor[2]*Nt[i];
//		nmtmp[0][0]+=cp[tmesh[eid].IEN[i]].coor[0]*dNdt[i][0];
//		nmtmp[0][1]+=cp[tmesh[eid].IEN[i]].coor[1]*dNdt[i][0];
//		nmtmp[0][2]+=cp[tmesh[eid].IEN[i]].coor[2]*dNdt[i][0];
//		nmtmp[1][0]+=cp[tmesh[eid].IEN[i]].coor[0]*dNdt[i][1];
//		nmtmp[1][1]+=cp[tmesh[eid].IEN[i]].coor[1]*dNdt[i][1];
//		nmtmp[1][2]+=cp[tmesh[eid].IEN[i]].coor[2]*dNdt[i][1];
//	}
//	norm[0]=nmtmp[0][1]*nmtmp[1][2]-nmtmp[0][2]*nmtmp[1][1];
//	norm[1]=nmtmp[0][2]*nmtmp[1][0]-nmtmp[0][0]*nmtmp[1][2];
//	norm[2]=nmtmp[0][0]*nmtmp[1][1]-nmtmp[0][1]*nmtmp[1][0];
//	double len=sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
//	norm[0]/=len; norm[1]/=len; norm[2]/=len;
//}

void TruncatedTspline_3D::ElementSubdivide_8(int eid)
{
	int pid[19];
	int edid[54];
	int fcid[36];
	int i, j, loc(0), edloc(0), fcloc(0);
	Vertex3D ptmp1;//body point
	for (i = 0; i < 8; i++)
	{
		ptmp1.coor[0] += .125*cp[tmesh[eid].cnct[i]].coor[0];
		ptmp1.coor[1] += .125*cp[tmesh[eid].cnct[i]].coor[1];
		ptmp1.coor[2] += .125*cp[tmesh[eid].cnct[i]].coor[2];
	}
	cp.push_back(ptmp1);
	pid[loc] = cp.size() - 1; loc++;
	//6 faces
	for (i = 0; i < 6; i++)
	{
		int itmp(tmesh[eid].face[i]);
		int dir(0), pos(0);
		SolidFaceDirection(eid,i,dir,pos);
		if (tmface[itmp].act == 1)
		{
			FaceSubdivision_4(itmp);
		}
		else if (tmface[itmp].act == 0 && tmface[itmp].Tedge.size() == 1)//if the original face is split into two, we need to further split
		{
			int spl_dir(0);
			if (tmedge[tmface[itmp].Tedge[0]].pt[0] == tmedge[tmface[itmp].edge[0]].midpt || tmedge[tmface[itmp].Tedge[0]].pt[1] == tmedge[tmface[itmp].edge[0]].midpt)
				spl_dir = 1;
			FaceSubdivision_24(itmp,spl_dir);
		}
		pid[loc] = tmface[itmp].ctpt; loc++;
		if (dir == 0)//same direction
		{
			for (j = 0; j < 4; j++)
			{
				edid[edloc] = tmface[itmp].Tedge[(pos+j)%4]; edloc++;
				fcid[fcloc] = tmface[itmp].chd[(pos+j)%4]; fcloc++;
			}
		}
		else//opposite direction
		{
			for (j = 0; j < 4; j++)
			{
				edid[edloc] = tmface[itmp].Tedge[(pos+4-j)%4]; edloc++;
				fcid[fcloc] = tmface[itmp].chd[(pos+4-j)%4]; fcloc++;
			}
		}
	}
	//12 edges
	for (i = 0; i < 12; i++)
	{
		int itmp(tmesh[eid].edge[i]);
		int dir(0);
		if (tmedge[itmp].pt[0] == tmesh[eid].cnct[solid_ed[i][1]]) dir = 1;

		pid[loc] = tmedge[itmp].midpt; loc++;
		if (dir == 0)
		{
			edid[edloc] = tmedge[itmp].chd[0]; edloc++;
			edid[edloc] = tmedge[itmp].chd[1]; edloc++;
		}
		else
		{
			edid[edloc] = tmedge[itmp].chd[1]; edloc++;
			edid[edloc] = tmedge[itmp].chd[0]; edloc++;
		}
	}
	//construct 6 new edges
	for (i = 0; i < 6; i++)
	{
		Edge3D edtmp;
		edtmp.pt[0] = pid[0]; edtmp.pt[1] = tmface[tmesh[eid].face[i]].ctpt;
		edtmp.len = tmedge[tmesh[eid].edge[fc_ppd_ed[i]]].len/2.;
		edtmp.lev = tmedge[tmesh[eid].edge[fc_ppd_ed[i]]].lev+1;
		tmedge.push_back(edtmp);
		edid[edloc] = tmedge.size() - 1; edloc++;
	}
	//set new edge ids
	vector<int> ednew_id(edid,edid+54);
	//construct 12 new faces
	int fcnew_cnct[12][4] = { { 11, 2, 0, 5 }, { 2, 12, 3, 0 }, { 0, 3, 13, 4 }, { 5, 0, 4, 14 }, //xy-plane
	{ 7, 1, 0, 2 }, { 1, 9, 4, 0 }, { 0, 4, 17, 6 }, { 2, 0, 6, 15 },//yz
	{ 10, 1, 0, 5 }, { 1, 8, 3, 0 }, { 0, 3, 16, 6 }, { 5, 0, 6, 18 } };//xz
	//int fcnew_edge[12][4] = { { 7, 49, 52, 17 }, { 5, 11, 50, 49 }, { 50, 9, 15, 51 }, { 52, 51, 13, 19 }, //xy-plane
	//{ 3, 48, 49, 4 }, { 1, 12, 51, 48 }, { 51, 14, 22, 53 }, { 49, 53, 20, 6 },//yz
	//{ 0, 48, 52, 16 }, { 2, 8, 50, 48 }, { 50, 10, 21, 53 }, { 52, 53, 23, 18 } };//xz
	int fc_lev_ref[12] = {0,0,0,0,2,2,2,2,1,1,1,1};
	for (i = 0; i < 12; i++)
	{
		Face3D fctmp;
		for (j = 0; j < 4; j++)
		{
			fctmp.cnct[j] = pid[fcnew_cnct[i][j]];
		}
		EdgeIndex_in_Face(fctmp,ednew_id);
		fctmp.lev = tmface[tmesh[eid].face[fc_lev_ref[i]]].lev + 2;
		tmface.push_back(fctmp);
		fcid[fcloc] = tmface.size() - 1; fcloc++;
	}
	vector<int> fcnew_id(fcid, fcid + 36);
	//construct 8 new solids
	int cnid[8];
	for (i = 0; i < 8; i++) cnid[i] = tmesh[eid].cnct[i];
	int slnew_cnct[8][8] = { { cnid[0], pid[7], pid[1], pid[10], pid[11], pid[2], pid[0], pid[5] }, { pid[7], cnid[1], pid[8], pid[1], pid[2], pid[12], pid[3], pid[0] },
	{ pid[1], pid[8], cnid[2], pid[9], pid[0], pid[3], pid[13], pid[4] }, { pid[10], pid[1], pid[9], cnid[3], pid[5], pid[0], pid[4], pid[14] },
	{ pid[11], pid[2], pid[0], pid[5], cnid[4], pid[15], pid[6], pid[18] }, { pid[2], pid[12], pid[3], pid[0], pid[15], cnid[5], pid[16], pid[6] },
	{ pid[0], pid[3], pid[13], pid[4], pid[6], pid[16], cnid[6], pid[17] }, { pid[5], pid[0], pid[4], pid[14], pid[18], pid[6], pid[17], cnid[7] },
	};
	int fc_rot[8][3] = { { 0, 1, 4 }, { 0, 1, 2 }, { 0, 2, 3 }, { 0, 3, 4 }, { 1, 4, 5 }, { 1, 2, 5 }, { 2, 3, 5 }, { 3, 4, 5 } };
	double uvw[3][2] = { { 0., tmedge[tmesh[eid].edge[0]].len / 2.}, { 0., tmedge[tmesh[eid].edge[3]].len / 2.}, { 0., tmedge[tmesh[eid].edge[4]].len / 2.}};
	double chd_o[8][3] = { { uvw[0][0], uvw[1][0], uvw[2][0] }, { uvw[0][1], uvw[1][0], uvw[2][0] }, { uvw[0][1], uvw[1][1], uvw[2][0] }, { uvw[0][0], uvw[1][1], uvw[2][0] },
	{ uvw[0][0], uvw[1][0], uvw[2][1] }, { uvw[0][1], uvw[1][0], uvw[2][1] }, { uvw[0][1], uvw[1][1], uvw[2][1] }, { uvw[0][0], uvw[1][1], uvw[2][1] } };
	tmesh[eid].chd.clear();
	tmesh[eid].chd_o.clear();
	for (i = 0; i < 8; i++)
	{
		Element3D hxtmp;
		for (j = 0; j < 8; j++)
		{
			hxtmp.cnct[j] = slnew_cnct[i][j];
		}
		hxtmp.lev = tmesh[eid].lev + 3;
		hxtmp.prt = eid;
		EdgeFaceIndex_in_Solid(hxtmp,ednew_id,fcnew_id);
		for (j = 0; j < 3; j++)
		{
			hxtmp.nbrot[fc_rot[i][j]] = tmesh[eid].nbrot[fc_rot[i][j]];
		}
		tmesh.push_back(hxtmp);
		tmesh[eid].chd.push_back(tmesh.size()-1);
		array<double, 3> otmp;
		for (j = 0; j < 3; j++)
		{
			otmp[j] = chd_o[i][j];
		}
		tmesh[eid].chd_o.push_back(otmp);
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline_3D::ElementSubdivide_5(int eid, int sl_rfdir)//dir=x,y,z
{
	int pid[10];
	vector<int> edid(33);
	vector<int> fcid(20);
	int i, j, loc(0), edloc(0), fcloc(0);

	int fc_bsc[3][4] = { { 0, 1, 5, 3 }, { 0, 2, 5, 4 }, { 1, 2, 3, 4 } };
	int fc_dir[3][4] = { { 1, 1, 1, 1 }, { 0, 1, 0, 1 }, { 0, 0, 0, 0 } };
	int fc_sbd[3][2] = { { 2, 4 }, { 1, 3 }, { 0, 5 } };
	int ed_tgt[3][8] = { { 1, 3, 4, 5, 6, 7, 9, 11 }, { 0, 2, 4, 5, 6, 7, 8, 10 }, { 0, 1, 2, 3, 8, 9, 10, 11 } };

	int ed_0[3][4] = { { 0, 2, 8, 10 }, { 1, 3, 9, 11 }, { 4, 5, 6, 7 } };//original (base) edges from eid

	//4 faces to be bisected
	for (i = 0; i < 4; i++)
	{
		int itmp(tmesh[eid].face[fc_bsc[sl_rfdir][i]]);
		int dir(0), pos(0);
		SolidFaceDirection(eid, fc_bsc[sl_rfdir][i], dir, pos);
		int fc_rfdir = (fc_dir[sl_rfdir][i] + pos) % 4;
		if (dir == 1) fc_rfdir = (fc_dir[sl_rfdir][i] + pos + 1) % 4;
		if (fc_rfdir == 2) fc_rfdir = 0;
		if (fc_rfdir == 3) fc_rfdir = 1;
		if (tmface[itmp].act == 1)
		{
			FaceSubdivision_2(itmp, fc_rfdir);
			edid[edloc] = tmface[itmp].Tedge[0]; edloc++;
			fcid[fcloc] = tmface[itmp].chd[0]; fcloc++;
			fcid[fcloc] = tmface[itmp].chd[1]; fcloc++;
		}
		else
		{
			//if the original face is split into two, we need to further split
			if (tmface[itmp].ctpt != -1)//subdivided
			{
				int edbs, fcbs[2];
				Construct_BaseFace_Subdv(itmp, fc_rfdir, edbs, fcbs);
				edid[edloc] = edbs; edloc++;
				fcid[fcloc] = fcbs[0]; fcloc++;
				fcid[fcloc] = fcbs[1]; fcloc++;
			}
			else if (tmface[itmp].Tedge.size() == 1)//bisected
			{
				if (tmedge[tmface[itmp].Tedge[0]].pt[0] == tmedge[tmface[itmp].edge[fc_rfdir]].midpt
					|| tmedge[tmface[itmp].Tedge[0]].pt[1] == tmedge[tmface[itmp].edge[fc_rfdir]].midpt)//aligned
				{
					edid[edloc] = tmface[itmp].Tedge[0]; edloc++;
					fcid[fcloc] = tmface[itmp].chd[0]; fcloc++;
					fcid[fcloc] = tmface[itmp].chd[1]; fcloc++;
				}
				else//crossed
				{
					int edbs, fcbs[2];
					FaceSubdivision_24(itmp, fc_rfdir);
					Construct_BaseFace_Subdv(itmp, fc_rfdir, edbs, fcbs);
					edid[edloc] = edbs; edloc++;
					fcid[fcloc] = fcbs[0]; fcloc++;
					fcid[fcloc] = fcbs[1]; fcloc++;
				}

			}
		}
	}
	//2 faces to be subdivided
	for (i = 0; i < 2; i++)
	{
		int itmp(tmesh[eid].face[fc_sbd[sl_rfdir][i]]);
		int dir(0), pos(0);
		SolidFaceDirection(eid, fc_sbd[sl_rfdir][i], dir, pos);
		if (tmface[itmp].act == 1)
		{
			FaceSubdivision_4(itmp);
		}
		else
		{
			//if the original face is split into two, we need to further split
			if (tmface[itmp].act == 0 && tmface[itmp].Tedge.size() == 1)
			{
				int spl_dir(1);
				if (tmedge[tmface[itmp].Tedge[0]].pt[0] == tmedge[tmface[itmp].edge[0]].midpt || tmedge[tmface[itmp].Tedge[0]].pt[1] == tmedge[tmface[itmp].edge[0]].midpt)
					spl_dir = 0;
				FaceSubdivision_24(itmp, spl_dir);
			}
		}
		pid[loc] = tmface[itmp].ctpt; loc++;
		for (j = 0; j < 4; j++)
		{
			edid[edloc] = tmface[itmp].Tedge[j]; edloc++;
			fcid[fcloc] = tmface[itmp].chd[j]; fcloc++;
		}
	}
	//8 edges
	for (i = 0; i < 8; i++)
	{
		int itmp(tmesh[eid].edge[ed_tgt[sl_rfdir][i]]);
		pid[loc] = tmedge[itmp].midpt; loc++;
		edid[edloc] = tmedge[itmp].chd[0]; edloc++;
		edid[edloc] = tmedge[itmp].chd[1]; edloc++;
	}
	//4 untouched base edegs
	for (i = 0; i < 4; i++)
	{
		edid[edloc] = tmesh[eid].edge[ed_0[sl_rfdir][i]]; edloc++;
	}
	Edge3D edtmp;
	edtmp.pt[0] = tmface[tmesh[eid].face[fc_sbd[sl_rfdir][0]]].ctpt; edtmp.pt[1] = tmface[tmesh[eid].face[fc_sbd[sl_rfdir][1]]].ctpt;
	edtmp.len = tmedge[tmesh[eid].edge[ed_0[sl_rfdir][0]]].len;
	edtmp.lev = tmedge[tmesh[eid].edge[ed_0[sl_rfdir][0]]].lev;
	tmedge.push_back(edtmp);
	edid[edloc] = tmedge.size() - 1; edloc++;
	//construct 4 new faces
	int fcnew_cnct[4][4] = { { 4, 5, 0, 1 }, { 1, 0, 6, 7 }, { 3, 2, 0, 1 }, { 1, 0, 8, 9 } };
	int fc_lev_ref[4] = { 0, 0, 1, 1 };
	if (sl_rfdir == 1)
	{
		int fc_tmp[4][4] = { { 4, 0, 1, 7 }, { 0, 5, 6, 1 }, { 2, 3, 1, 0 }, { 0, 1, 9, 8 } };
		int ref_tmp[4] = { 0, 0, 2, 2 };
		for (i = 0; i < 4; i++)
		{
			fc_lev_ref[i] = ref_tmp[i];
			for (j = 0; j < 4; j++)
			{
				fcnew_cnct[i][j] = fc_tmp[i][j];
			}
		}
	}
	if (sl_rfdir == 2)
	{
		int fc_tmp[4][4] = { { 2, 0, 1, 6 }, { 0, 4, 8, 1 }, { 5, 0, 1, 9 }, { 0, 3, 7, 1 } };
		int ref_tmp[4] = { 2, 2, 1, 1 };
		for (i = 0; i < 4; i++)
		{
			fc_lev_ref[i] = ref_tmp[i];
			for (j = 0; j < 4; j++)
			{
				fcnew_cnct[i][j] = fc_tmp[i][j];
			}
		}
	}
	for (i = 0; i < 4; i++)
	{
		Face3D fctmp;
		for (j = 0; j < 4; j++)
		{
			fctmp.cnct[j] = pid[fcnew_cnct[i][j]];
		}
		EdgeIndex_in_Face(fctmp, edid);
		fctmp.lev = tmface[tmesh[eid].face[fc_lev_ref[i]]].lev + 1;
		tmface.push_back(fctmp);
		fcid[fcloc] = tmface.size() - 1; fcloc++;
	}
	//construct 4 new solids
	int cnid[8];
	for (i = 0; i < 8; i++) cnid[i] = tmesh[eid].cnct[i];
	int slnew_cnct[4][8] = { { cnid[0], cnid[1], pid[2], pid[3], pid[4], pid[5], pid[0], pid[1] }, { pid[3], pid[2], cnid[2], cnid[3], pid[1], pid[0], pid[6], pid[7] },
	{ pid[4], pid[5], pid[0], pid[1], cnid[4], cnid[5], pid[8], pid[9] }, { pid[1], pid[0], pid[6], pid[7], pid[9], pid[8], cnid[6], cnid[7] } };
	int fc_rot[4][4] = { { 0, 1, 2, 4 }, { 0, 2, 3, 4 }, { 1, 2, 4, 5 }, { 2, 3, 4, 5 } };
	double uvw[3][2] = { { 0., tmedge[tmesh[eid].edge[0]].len / 2. }, { 0., tmedge[tmesh[eid].edge[3]].len / 2. }, { 0., tmedge[tmesh[eid].edge[4]].len / 2. } };
	double chd_o[4][3] = { { uvw[0][0], uvw[1][0], uvw[2][0] }, { uvw[0][0], uvw[1][1], uvw[2][0] }, { uvw[0][0], uvw[1][0], uvw[2][1] }, { uvw[0][0], uvw[1][1], uvw[2][1] } };
	if (sl_rfdir == 1)
	{
		int cnct_tmp[4][8] = { { cnid[0], pid[2], pid[3], cnid[3], pid[4], pid[0], pid[1], pid[7] }, { pid[2], cnid[1], cnid[2], pid[3], pid[0], pid[5], pid[6], pid[1] },
		{ pid[4], pid[0], pid[1], pid[7], cnid[4], pid[8], pid[9], cnid[7] }, { pid[0], pid[5], pid[6], pid[1], pid[8], cnid[5], cnid[6], pid[9] } };
		int fcrot_tmp[4][4] = { { 0, 1, 3, 4 }, { 0, 1, 2, 3 }, { 1, 3, 4, 5 }, { 1, 2, 3, 5 } };
		double otmp[4][3] = { { uvw[0][0], uvw[1][0], uvw[2][0] }, { uvw[0][1], uvw[1][0], uvw[2][0] }, { uvw[0][0], uvw[1][0], uvw[2][1] }, { uvw[0][0], uvw[1][1], uvw[2][1] } };
		for (i = 0; i < 4; i++)
		{
			for (j = 0; j < 8; j++) slnew_cnct[i][j] = cnct_tmp[i][j];
			for (j = 0; j < 4; j++) fc_rot[i][j] = fcrot_tmp[i][j];
			for (j = 0; j < 3; j++) chd_o[i][j] = otmp[i][j];
		}
	}
	else if (sl_rfdir == 2)
	{
		int cnct_tmp[4][8] = { { cnid[0], pid[2], pid[0], pid[5], cnid[4], pid[6], pid[1], pid[9] }, { pid[2], cnid[1], pid[3], pid[0], pid[6], cnid[5], pid[7], pid[1] },
		{ pid[5], pid[0], pid[4], cnid[3], pid[9], pid[1], pid[8], cnid[7] }, { pid[0], pid[3], cnid[2], pid[4], pid[1], pid[7], cnid[6], pid[8] } };
		int fcrot_tmp[4][4] = { { 0, 1, 4, 5 }, { 0, 1, 2, 5 }, { 0, 3, 4, 5 }, { 0, 2, 3, 5 } };
		double otmp[4][3] = { { uvw[0][0], uvw[1][0], uvw[2][0] }, { uvw[0][1], uvw[1][0], uvw[2][0] }, { uvw[0][0], uvw[1][1], uvw[2][0] }, { uvw[0][1], uvw[1][1], uvw[2][0] } };
		for (i = 0; i < 4; i++)
		{
			for (j = 0; j < 8; j++) slnew_cnct[i][j] = cnct_tmp[i][j];
			for (j = 0; j < 4; j++) fc_rot[i][j] = fcrot_tmp[i][j];
			for (j = 0; j < 3; j++) chd_o[i][j] = otmp[i][j];
		}
	}
	tmesh[eid].chd.clear();
	tmesh[eid].chd_o.clear();
	for (i = 0; i < 4; i++)
	{
		Element3D hxtmp;
		for (j = 0; j < 8; j++)
		{
			hxtmp.cnct[j] = slnew_cnct[i][j];
		}
		hxtmp.lev = tmesh[eid].lev + 1;
		hxtmp.prt = eid;
		EdgeFaceIndex_in_Solid(hxtmp, edid, fcid);
		for (j = 0; j < 4; j++)
		{
			hxtmp.nbrot[fc_rot[i][j]] = tmesh[eid].nbrot[fc_rot[i][j]];
		}
		tmesh.push_back(hxtmp);
		tmesh[eid].chd.push_back(tmesh.size() - 1);
		array<double, 3> otmp;
		for (j = 0; j < 3; j++)
		{
			otmp[j] = chd_o[i][j];
		}
		tmesh[eid].chd_o.push_back(otmp);
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline_3D::ElementSubdivide_4(int eid, int sl_rfdir)//dir=x,y,z
{
	int pid[10];
	vector<int> edid(33);
	vector<int> fcid(20);
	int i, j, loc(0), edloc(0), fcloc(0);

	int fc_bsc[3][4] = { { 0, 1, 5, 3 }, { 0, 2, 5, 4 }, { 1, 2, 3, 4 } };
	int fc_dir[3][4] = { { 1, 1, 1, 1 }, { 0, 1, 0, 1 }, { 0, 0, 0, 0 } };
	int fc_sbd[3][2] = { { 2, 4 }, { 1, 3 }, {0, 5} };
	int ed_tgt[3][8] = { { 1, 3, 4, 5, 6, 7, 9, 11 }, { 0, 2, 4, 5, 6, 7, 8, 10 }, { 0,1,2,3,8,9,10,11 } };

	int ed_0[3][4] = { { 0, 2, 8, 10}, { 1, 3, 9, 11}, { 4, 5, 6, 7} };//original (base) edges from eid

	//4 faces to be bisected
	for (i = 0; i < 4; i++)
	{
		int itmp(tmesh[eid].face[fc_bsc[sl_rfdir][i]]);
		int dir(0), pos(0);
		SolidFaceDirection(eid, fc_bsc[sl_rfdir][i], dir, pos);
		int fc_rfdir = (fc_dir[sl_rfdir][i] + pos) % 4;
		if (dir == 1) fc_rfdir = (fc_dir[sl_rfdir][i] + pos + 1) % 4;
		if (fc_rfdir == 2) fc_rfdir = 0;
		if (fc_rfdir == 3) fc_rfdir = 1;
		if (tmface[itmp].act == 1)
		{
			FaceSubdivision_2(itmp, fc_rfdir);
			edid[edloc] = tmface[itmp].Tedge[0]; edloc++;
			fcid[fcloc] = tmface[itmp].chd[0]; fcloc++;
			fcid[fcloc] = tmface[itmp].chd[1]; fcloc++;
		}
		else
		{
			//if the original face is split into two, we need to further split
			if (tmface[itmp].ctpt != -1)//subdivided
			{
				int edbs, fcbs[2];
				Construct_BaseFace_Subdv(itmp, fc_rfdir, edbs, fcbs);
				edid[edloc] = edbs; edloc++;
				fcid[fcloc] = fcbs[0]; fcloc++;
				fcid[fcloc] = fcbs[1]; fcloc++;
			}
			else if (tmface[itmp].Tedge.size() == 1)//bisected
			{
				if (tmedge[tmface[itmp].Tedge[0]].pt[0] == tmedge[tmface[itmp].edge[fc_rfdir]].midpt
					|| tmedge[tmface[itmp].Tedge[0]].pt[1] == tmedge[tmface[itmp].edge[fc_rfdir]].midpt)//aligned
				{
					edid[edloc] = tmface[itmp].Tedge[0]; edloc++;
					fcid[fcloc] = tmface[itmp].chd[0]; fcloc++;
					fcid[fcloc] = tmface[itmp].chd[1]; fcloc++;
				}
				else//crossed
				{
					int edbs, fcbs[2];
					FaceSubdivision_24(itmp, fc_rfdir);
					Construct_BaseFace_Subdv(itmp, fc_rfdir, edbs, fcbs);
					edid[edloc] = edbs; edloc++;
					fcid[fcloc] = fcbs[0]; fcloc++;
					fcid[fcloc] = fcbs[1]; fcloc++;
				}

			}
		}
	}
	//2 faces to be subdivided
	for (i = 0; i < 2; i++)
	{
		int itmp(tmesh[eid].face[fc_sbd[sl_rfdir][i]]);
		int dir(0), pos(0);
		SolidFaceDirection(eid, fc_sbd[sl_rfdir][i], dir, pos);
		if (tmface[itmp].act == 1)
		{
			FaceSubdivision_4(itmp);
		}
		else
		{
			//if the original face is split into two, we need to further split
			if (tmface[itmp].act == 0 && tmface[itmp].Tedge.size() == 1)
			{
				int spl_dir(1);
				if (tmedge[tmface[itmp].Tedge[0]].pt[0] == tmedge[tmface[itmp].edge[0]].midpt || tmedge[tmface[itmp].Tedge[0]].pt[1] == tmedge[tmface[itmp].edge[0]].midpt)
					spl_dir = 0;
				FaceSubdivision_24(itmp, spl_dir);
			}
		}
		pid[loc] = tmface[itmp].ctpt; loc++;
		for (j = 0; j < 4; j++)
		{
			edid[edloc] = tmface[itmp].Tedge[j]; edloc++;
			fcid[fcloc] = tmface[itmp].chd[j]; fcloc++;
		}
	}
	//8 edges
	for (i = 0; i < 8; i++)
	{
		int itmp(tmesh[eid].edge[ed_tgt[sl_rfdir][i]]);
		pid[loc] = tmedge[itmp].midpt; loc++;
		edid[edloc] = tmedge[itmp].chd[0]; edloc++;
		edid[edloc] = tmedge[itmp].chd[1]; edloc++;
	}
	//4 untouched base edegs
	for (i = 0; i < 4; i++)
	{
		edid[edloc] = tmesh[eid].edge[ed_0[sl_rfdir][i]]; edloc++;
	}
	Edge3D edtmp;
	edtmp.pt[0] = tmface[tmesh[eid].face[fc_sbd[sl_rfdir][0]]].ctpt; edtmp.pt[1] = tmface[tmesh[eid].face[fc_sbd[sl_rfdir][1]]].ctpt;
	edtmp.len = tmedge[tmesh[eid].edge[ed_0[sl_rfdir][0]]].len;
	edtmp.lev = tmedge[tmesh[eid].edge[ed_0[sl_rfdir][0]]].lev;
	tmedge.push_back(edtmp);
	edid[edloc] = tmedge.size() - 1; edloc++;
	//construct 4 new faces
	int fcnew_cnct[4][4] = { { 4, 5, 0, 1 }, { 1, 0, 6, 7 }, { 3, 2, 0, 1 }, {1,0,8,9} };
	int fc_lev_ref[4] = {0,0,1,1};
	if (sl_rfdir == 1)
	{
		int fc_tmp[4][4] = { { 4, 0, 1, 7 }, { 0, 5, 6, 1 }, { 2, 3, 1, 0 }, {0,1,9,8} };
		int ref_tmp[4] = {0,0,2,2};
		for (i = 0; i < 4; i++)
		{
			fc_lev_ref[i] = ref_tmp[i];
			for (j = 0; j < 4; j++)
			{
				fcnew_cnct[i][j] = fc_tmp[i][j];
			}
		}
	}
	if (sl_rfdir == 2)
	{
		int fc_tmp[4][4] = { { 2, 0, 1, 6 }, { 0, 4, 8, 1 }, { 5, 0, 1, 9 }, { 0, 3, 7, 1 } };
		int ref_tmp[4] = { 2, 2, 1, 1 };
		for (i = 0; i < 4; i++)
		{
			fc_lev_ref[i] = ref_tmp[i];
			for (j = 0; j < 4; j++)
			{
				fcnew_cnct[i][j] = fc_tmp[i][j];
			}
		}
	}
	for (i = 0; i < 4; i++)
	{
		Face3D fctmp;
		for (j = 0; j < 4; j++)
		{
			fctmp.cnct[j] = pid[fcnew_cnct[i][j]];
		}
		EdgeIndex_in_Face(fctmp, edid);
		fctmp.lev = tmface[tmesh[eid].face[fc_lev_ref[i]]].lev+1;
		tmface.push_back(fctmp);
		fcid[fcloc] = tmface.size() - 1; fcloc++;
	}
	//construct 4 new solids
	int cnid[8];
	for (i = 0; i < 8; i++) cnid[i] = tmesh[eid].cnct[i];
	int slnew_cnct[4][8] = { { cnid[0], cnid[1], pid[2], pid[3], pid[4], pid[5], pid[0], pid[1] }, { pid[3], pid[2], cnid[2], cnid[3], pid[1], pid[0], pid[6], pid[7] },
	{ pid[4], pid[5], pid[0], pid[1], cnid[4], cnid[5], pid[8], pid[9] }, { pid[1], pid[0], pid[6], pid[7], pid[9], pid[8], cnid[6], cnid[7] } };
	int fc_rot[4][4] = { { 0, 1, 2, 4 }, { 0, 2, 3, 4 }, { 1, 2, 4, 5 }, {2,3,4,5} };
	double uvw[3][2] = { { 0., tmedge[tmesh[eid].edge[0]].len / 2. }, { 0., tmedge[tmesh[eid].edge[3]].len / 2. }, { 0., tmedge[tmesh[eid].edge[4]].len / 2. } };
	double chd_o[4][3] = { { uvw[0][0], uvw[1][0], uvw[2][0] }, { uvw[0][0], uvw[1][1], uvw[2][0] }, { uvw[0][0], uvw[1][0], uvw[2][1] }, { uvw[0][0], uvw[1][1], uvw[2][1] } };
	if (sl_rfdir == 1)
	{
		int cnct_tmp[4][8] = { { cnid[0], pid[2], pid[3], cnid[3], pid[4], pid[0], pid[1], pid[7] }, { pid[2], cnid[1], cnid[2], pid[3], pid[0], pid[5], pid[6], pid[1] },
		{ pid[4], pid[0], pid[1], pid[7], cnid[4], pid[8], pid[9], cnid[7] }, { pid[0], pid[5], pid[6], pid[1], pid[8], cnid[5], cnid[6], pid[9] } };
		int fcrot_tmp[4][4] = { { 0, 1, 3, 4 }, { 0, 1, 2, 3 }, { 1, 3, 4, 5 }, { 1, 2, 3, 5 } };
		double otmp[4][3] = { { uvw[0][0], uvw[1][0], uvw[2][0] }, { uvw[0][1], uvw[1][0], uvw[2][0] }, { uvw[0][0], uvw[1][0], uvw[2][1] }, { uvw[0][0], uvw[1][1], uvw[2][1] } };
		for (i = 0; i < 4; i++)
		{
			for (j = 0; j < 8; j++) slnew_cnct[i][j] = cnct_tmp[i][j];
			for (j = 0; j < 4; j++) fc_rot[i][j] = fcrot_tmp[i][j];
			for (j = 0; j < 3; j++) chd_o[i][j] = otmp[i][j];
		}
	}
	else if (sl_rfdir == 2)
	{
		int cnct_tmp[4][8] = { { cnid[0], pid[2], pid[0], pid[5], cnid[4], pid[6], pid[1], pid[9] }, { pid[2], cnid[1], pid[3], pid[0], pid[6], cnid[5], pid[7], pid[1] },
		{ pid[5], pid[0], pid[4], cnid[3], pid[9], pid[1], pid[8], cnid[7] }, { pid[0], pid[3], cnid[2], pid[4], pid[1], pid[7], cnid[6], pid[8] } };
		int fcrot_tmp[4][4] = { { 0, 1, 4, 5 }, { 0, 1, 2, 5 }, { 0, 3, 4, 5 }, { 0, 2, 3, 5 } };
		double otmp[4][3] = { { uvw[0][0], uvw[1][0], uvw[2][0] }, { uvw[0][1], uvw[1][0], uvw[2][0] }, { uvw[0][0], uvw[1][1], uvw[2][0] }, { uvw[0][1], uvw[1][1], uvw[2][0] } };
		for (i = 0; i < 4; i++)
		{
			for (j = 0; j < 8; j++) slnew_cnct[i][j] = cnct_tmp[i][j];
			for (j = 0; j < 4; j++) fc_rot[i][j] = fcrot_tmp[i][j];
			for (j = 0; j < 3; j++) chd_o[i][j] = otmp[i][j];
		}
	}
	tmesh[eid].chd.clear();
	tmesh[eid].chd_o.clear();
	for (i = 0; i < 4; i++)
	{
		Element3D hxtmp;
		for (j = 0; j < 8; j++)
		{
			hxtmp.cnct[j] = slnew_cnct[i][j];
		}
		hxtmp.lev = tmesh[eid].lev + 1;
		hxtmp.prt = eid;
		EdgeFaceIndex_in_Solid(hxtmp, edid, fcid);
		for (j = 0; j < 4; j++)
		{
			hxtmp.nbrot[fc_rot[i][j]] = tmesh[eid].nbrot[fc_rot[i][j]];
		}
		tmesh.push_back(hxtmp);
		tmesh[eid].chd.push_back(tmesh.size() - 1);
		array<double, 3> otmp;
		for (j = 0; j < 3; j++)
		{
			otmp[j] = chd_o[i][j];
		}
		tmesh[eid].chd_o.push_back(otmp);
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline_3D::ElementSubdivide_2(int eid, int sl_rfdir)//dir=x,y,z
{
	int pid[4];
	vector<int> edid(20);
	vector<int> fcid(11);
	int i, j, loc(0), edloc(0), fcloc(0);

	int fc_tgt[3][4] = { { 0, 1, 5, 3 }, { 0, 2, 5, 4 }, {1,2,3,4} };
	int fc_dir[3][4] = { { 0, 0, 0, 0 }, { 1, 0, 1, 0 }, {1,1,1,1} };
	int ed_tgt[3][4] = { { 0, 2, 10, 8 }, { 3, 1, 9, 11 }, {4,5,6,7} };

	int fc_0[3][2] = { { 2, 4 }, { 1, 3 }, { 0, 5 } };//original (base) faces from eid
	int ed_0[3][8] = { { 1,3,4,5,6,7,9,11 }, { 0,2,4,5,6,7,8,10 }, { 0,1,2,3,8,9,10,11 } };//original (base) edges from eid

	//4 faces
	for (i = 0; i < 4; i++)
	{
		int itmp(tmesh[eid].face[fc_tgt[sl_rfdir][i]]);
		int dir(0), pos(0);
		SolidFaceDirection(eid, fc_tgt[sl_rfdir][i], dir, pos);
		int fc_rfdir = (fc_dir[sl_rfdir][i] + pos) % 4;
		if (dir == 1) fc_rfdir = (fc_dir[sl_rfdir][i] + pos+1) % 4;
		if (fc_rfdir == 2) fc_rfdir = 0;
		if (fc_rfdir == 3) fc_rfdir = 1;
		if (tmface[itmp].act == 1)
		{
			FaceSubdivision_2(itmp, fc_rfdir);
			edid[edloc] = tmface[itmp].Tedge[0]; edloc++;
			fcid[fcloc] = tmface[itmp].chd[0]; fcloc++;
			fcid[fcloc] = tmface[itmp].chd[1]; fcloc++;
		}
		else
		{
			//if the original face is split into two, we need to further split
			if (tmface[itmp].ctpt!=-1)//subdivided
			{
				int edbs, fcbs[2];
				Construct_BaseFace_Subdv(itmp, fc_rfdir, edbs, fcbs);
				edid[edloc] = edbs; edloc++;
				fcid[fcloc] = fcbs[0]; fcloc++;
				fcid[fcloc] = fcbs[1]; fcloc++;
			}
			else if (tmface[itmp].Tedge.size() == 1)//bisected
			{
				if (tmedge[tmface[itmp].Tedge[0]].pt[0] == tmedge[tmface[itmp].edge[fc_rfdir]].midpt 
					|| tmedge[tmface[itmp].Tedge[0]].pt[1] == tmedge[tmface[itmp].edge[fc_rfdir]].midpt)//aligned
				{
					edid[edloc] = tmface[itmp].Tedge[0]; edloc++;
					fcid[fcloc] = tmface[itmp].chd[0]; fcloc++;
					fcid[fcloc] = tmface[itmp].chd[1]; fcloc++;
				}
				else//crossed
				{
					int edbs, fcbs[2];
					FaceSubdivision_24(itmp, fc_rfdir);
					Construct_BaseFace_Subdv(itmp, fc_rfdir, edbs, fcbs);
					edid[edloc] = edbs; edloc++;
					fcid[fcloc] = fcbs[0]; fcloc++;
					fcid[fcloc] = fcbs[1]; fcloc++;
				}

			}
		}
	}
	//2 untouched base faces
	for (i = 0; i < 2; i++)
	{
		fcid[fcloc] = tmesh[eid].face[fc_0[sl_rfdir][i]]; fcloc++;
	}
	//4 edges
	for (i = 0; i < 4; i++)
	{
		int itmp(tmesh[eid].edge[ed_tgt[sl_rfdir][i]]);
		int dir(0);
		if (tmedge[itmp].pt[0] == tmesh[eid].cnct[solid_ed[ed_tgt[sl_rfdir][i]][1]]) dir = 1;

		pid[loc] = tmedge[itmp].midpt; loc++;
		//if (dir == 0)
		//{
			edid[edloc] = tmedge[itmp].chd[0]; edloc++;
			edid[edloc] = tmedge[itmp].chd[1]; edloc++;
		//}
		//else
		//{
		//	edid[edloc] = tmedge[itmp].chd[1]; edloc++;
		//	edid[edloc] = tmedge[itmp].chd[0]; edloc++;
		//}
	}
	//8 untouched base edegs
	for (i = 0; i < 8; i++)
	{
		edid[edloc] = tmesh[eid].edge[ed_0[sl_rfdir][i]]; edloc++;
	}
	//construct 1 new face
	int fcnew_cnct[4] = { 0,1,2,3};
	int fc_lev_ref(2);
	if (sl_rfdir == 1) fc_lev_ref = 1;
	if (sl_rfdir == 2) fc_lev_ref = 0;
	Face3D fctmp;
	for (j = 0; j < 4; j++)
	{
		fctmp.cnct[j] = pid[fcnew_cnct[j]];
	}
	EdgeIndex_in_Face(fctmp, edid);
	fctmp.lev = tmface[tmesh[eid].face[fc_lev_ref]].lev;
	tmface.push_back(fctmp);
	fcid[fcloc] = tmface.size() - 1; fcloc++;
	//construct 2 new solids
	int cnid[8];
	for (i = 0; i < 8; i++) cnid[i] = tmesh[eid].cnct[i];
	int slnew_cnct[2][8] = { { cnid[0], pid[0], pid[1], cnid[3], cnid[4], pid[3], pid[2], cnid[7] }, { pid[0], cnid[1], cnid[2], pid[1], pid[3], cnid[5], cnid[6], pid[2] }};
	int fc_rot[2][5] = { { 0, 1, 3, 4, 5 }, { 0, 1, 2, 3, 5 } };
	double uvw[3][2] = { { 0., tmedge[tmesh[eid].edge[0]].len / 2. }, { 0., tmedge[tmesh[eid].edge[3]].len / 2. }, { 0., tmedge[tmesh[eid].edge[4]].len / 2. } };
	double chd_o[2][3] = { { uvw[0][0], uvw[1][0], uvw[2][0] }, { uvw[0][1], uvw[1][0], uvw[2][0] }};
	if (sl_rfdir == 1)
	{
		int cnct_tmp[2][8] = { { cnid[0], cnid[1], pid[1], pid[0], cnid[4], cnid[5], pid[2], pid[3] }, { pid[0], pid[1], cnid[2], cnid[3], pid[3], pid[2], cnid[6], cnid[7] } };
		int fcrot_tmp[2][5] = { { 0, 1, 2, 4, 5 }, { 0, 2, 3, 4, 5 } };
		double otmp[2][3] = { { uvw[0][0], uvw[1][0], uvw[2][0] }, { uvw[0][0], uvw[1][1], uvw[2][0] } };
		for (i = 0; i < 2; i++)
		{
			for (j = 0; j < 8; j++) slnew_cnct[i][j] = cnct_tmp[i][j];
			for (j = 0; j < 5; j++) fc_rot[i][j] = fcrot_tmp[i][j];
			for (j = 0; j < 3; j++) chd_o[i][j] = otmp[i][j];
		}
	}
	else if (sl_rfdir == 2)
	{
		int cnct_tmp[2][8] = { { cnid[0], cnid[1], cnid[2], cnid[3], pid[0], pid[1], pid[2], pid[3] }, { pid[0], pid[1], pid[2], pid[3], cnid[4], cnid[5], cnid[6], cnid[7] } };
		int fcrot_tmp[2][5] = { { 0, 1, 2, 3, 4 }, { 1, 2, 3, 4, 5 } };
		double otmp[2][3] = { { uvw[0][0], uvw[1][0], uvw[2][0] }, { uvw[0][0], uvw[1][0], uvw[2][1] } };
		for (i = 0; i < 2; i++)
		{
			for (j = 0; j < 8; j++) slnew_cnct[i][j] = cnct_tmp[i][j];
			for (j = 0; j < 5; j++) fc_rot[i][j] = fcrot_tmp[i][j];
			for (j = 0; j < 3; j++) chd_o[i][j] = otmp[i][j];
		}
	}
	tmesh[eid].chd.clear();
	for (i = 0; i < 2; i++)
	{
		Element3D hxtmp;
		for (j = 0; j < 8; j++)
		{
			hxtmp.cnct[j] = slnew_cnct[i][j];
		}
		hxtmp.lev = tmesh[eid].lev + 1;
		hxtmp.prt = eid;
		//cout << "Assigning index...\n";
		EdgeFaceIndex_in_Solid(hxtmp, edid, fcid);
		for (j = 0; j < 5; j++)
		{
			hxtmp.nbrot[fc_rot[i][j]] = tmesh[eid].nbrot[fc_rot[i][j]];
		}
		tmesh.push_back(hxtmp);
		tmesh[eid].chd.push_back(tmesh.size() - 1);
		//cout << "Finish assigning " << tmesh.size() - 1<<"\n";
		array<double, 3> otmp;
		for (j = 0; j < 3; j++)
		{
			otmp[j] = chd_o[i][j];
		}
		tmesh[eid].chd_o.push_back(otmp);
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline_3D::FaceSubdivision_4(int fcid)
{
	int pid[5];
	int edid[12];
	int i, j, loc(0), edloc(0);
	Vertex3D ptmp1;
	for (i = 0; i < 4; i++)
	{
		ptmp1.coor[0] += 0.25*cp[tmface[fcid].cnct[i]].coor[0];
		ptmp1.coor[1] += 0.25*cp[tmface[fcid].cnct[i]].coor[1];
		ptmp1.coor[2] += 0.25*cp[tmface[fcid].cnct[i]].coor[2];
	}
	cp.push_back(ptmp1); 
	pid[loc] = cp.size() - 1; loc++;
	for (i = 0; i<4; i++)//4 edges
	{
		int itmp(tmface[fcid].edge[i]);
		int dir(0);
		if (tmedge[itmp].pt[1] == tmface[fcid].cnct[i]) dir = 1;

		if (tmedge[itmp].act == 1)
		{
			EdgeSubdivision_2(itmp);//add one point and two edges
			pid[loc] = cp.size()-1; loc++;
			if (dir == 0)
			{
				edid[edloc] = tmedge.size() - 2; edloc++;
				edid[edloc] = tmedge.size() - 1; edloc++;
			}
			else
			{
				edid[edloc] = tmedge.size() - 1; edloc++;
				edid[edloc] = tmedge.size() - 2; edloc++;
			}
		}
		else
		{
			pid[loc] = tmedge[itmp].midpt; loc++;
			int ied(tmedge[itmp].chd[0]);
			if (tmedge[ied].pt[0] == tmface[fcid].cnct[i] || tmedge[ied].pt[1] == tmface[fcid].cnct[i])
			{
				edid[edloc] = ied; edloc++;
				edid[edloc] = tmedge[itmp].chd[1]; edloc++;
			}
			else
			{
				edid[edloc] = tmedge[itmp].chd[1]; edloc++;
				edid[edloc] = ied; edloc++;
			}
		}
		//construct a new edge
		Edge3D edtmp;
		edtmp.pt[0] = pid[0]; edtmp.pt[1] = tmedge[itmp].midpt;
		edtmp.lev = tmedge[tmface[fcid].edge[(i + 1) % 4]].lev+1;
		edtmp.len = tmedge[tmface[fcid].edge[(i + 1) % 4]].len/2.;
		tmedge.push_back(edtmp);
		edid[edloc] = tmedge.size() - 1; edloc++;
	}

	int e_cnct[4][4] = { { tmface[fcid].cnct[0], pid[1], pid[0], pid[4] }, { pid[1], tmface[fcid].cnct[1], pid[2], pid[0] }, 
	{ pid[0], pid[2], tmface[fcid].cnct[2], pid[3] }, { pid[4], pid[0], pid[3], tmface[fcid].cnct[3] } };
	int e_edge[4][4] = { { edid[0], edid[2], edid[11], edid[10] }, { edid[1], edid[3], edid[5], edid[2] }, 
	{ edid[5], edid[4], edid[6], edid[8] }, { edid[11], edid[8], edid[7], edid[9] } };
	int enewid[4];
	double chd_org[4][2] = { { 0., 0. }, { tmedge[tmface[fcid].edge[0]].len / 2., 0. }, 
	{ tmedge[tmface[fcid].edge[0]].len / 2., tmedge[tmface[fcid].edge[1]].len / 2. }, { 0., tmedge[tmface[fcid].edge[1]].len / 2. } };
	int Tedge_id[4] = {2,5,8,11};
	vector<Face3D> etmp(4);
	for (i = 0; i<4; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 0;
		etmp[i].lev = tmface[fcid].lev + 2;
		etmp[i].prt = fcid;
		for (int j = 0; j<4; j++)
		{
			etmp[i].cnct[j] = e_cnct[i][j];
			etmp[i].edge[j] = e_edge[i][j];
		}
		tmface.push_back(etmp[i]);
		enewid[i] = tmface.size() - 1;
		tmface[fcid].chd.push_back(enewid[i]);
		//array<double, 2> chd_o_tmp = { chd_org[i][0], chd_org[i][1] };
		//tmface[fcid].chd_o.push_back(chd_o_tmp);
		tmface[fcid].Tedge.push_back(edid[Tedge_id[i]]);
	}

	tmface[fcid].ctpt = pid[0];
	tmface[fcid].act = 0;
}

void TruncatedTspline_3D::FaceSubdivision_2(int fcid, int fc_dir)//dir (0 or 1) indicates the edge that is split
{
	int pid[2];
	int edid[5];
	int i, j, loc(0), edloc(0);
	int iloc[2] = {fc_dir, fc_dir+2};

	for (i = 0; i<2; i++)//2 edges
	{
		int itmp(tmface[fcid].edge[iloc[i]]);
		int dir(0);
		if (tmedge[itmp].pt[1] == tmface[fcid].cnct[iloc[i]]) dir = 1;

		if (tmedge[itmp].act == 1)
		{
			EdgeSubdivision_2(itmp);//add one point and two edges
			pid[loc] = cp.size() - 1; loc++;
		}
		else
		{
			pid[loc] = tmedge[itmp].midpt; loc++;
		}
		if (dir==0)
		{
			edid[edloc] = tmedge[itmp].chd[0]; edloc++;
			edid[edloc] = tmedge[itmp].chd[1]; edloc++;
		}
		else
		{
			edid[edloc] = tmedge[itmp].chd[1]; edloc++;
			edid[edloc] = tmedge[itmp].chd[0]; edloc++;
		}
	}
	//construct a new edge
	Edge3D edtmp;
	edtmp.pt[0] = tmedge[tmface[fcid].edge[iloc[0]]].midpt; edtmp.pt[1] = tmedge[tmface[fcid].edge[iloc[1]]].midpt;
	edtmp.lev = tmedge[tmface[fcid].edge[fc_dir+1]].lev;
	edtmp.len = tmedge[tmface[fcid].edge[fc_dir + 1]].len;
	tmedge.push_back(edtmp);
	edid[edloc] = tmedge.size() - 1; edloc++;

	int e_cnct[2][4] = { { tmface[fcid].cnct[fc_dir], pid[0], pid[1], tmface[fcid].cnct[(fc_dir + 3) % 4] }, { pid[0], tmface[fcid].cnct[fc_dir + 1], tmface[fcid].cnct[fc_dir + 2], pid[1] } };
	int e_edge[2][4] = { { edid[0], edid[4], edid[3], tmface[fcid].edge[(fc_dir + 3) % 4] }, { edid[1], tmface[fcid].edge[fc_dir+1], edid[2], edid[4] } };
	int enewid[2];
	//double chd_org[4][2] = { { 0., 0. }, { tmedge[tmface[fcid].edge[0]].len / 2., 0. },
	//{ tmedge[tmface[fcid].edge[0]].len / 2., tmedge[tmface[fcid].edge[1]].len / 2. }, { 0., tmedge[tmface[fcid].edge[1]].len / 2. } };
	vector<Face3D> etmp(2);
	for (i = 0; i<2; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 0;
		etmp[i].lev = tmface[fcid].lev + 1;
		etmp[i].prt = fcid;
		for (int j = 0; j<4; j++)
		{
			etmp[i].cnct[j] = e_cnct[i][j];
			etmp[i].edge[j] = e_edge[i][j];
		}
		tmface.push_back(etmp[i]);
		enewid[i] = tmface.size() - 1;
		tmface[fcid].chd.push_back(enewid[i]);
		//array<double, 2> chd_o_tmp = { chd_org[i][0], chd_org[i][1] };
		//tmface[fcid].chd_o.push_back(chd_o_tmp);
	}

	tmface[fcid].Tedge.push_back(edid[4]);
	tmface[fcid].act = 0;
}

void TruncatedTspline_3D::FaceSubdivision_24(int fcid, int fc_dir)//dir (0 or 1) indicates the edge that is to be split
{
	int pid[5];
	int edid[12];
	int i, j, loc(0), edloc(0);
	Vertex3D ptmp1;
	for (i = 0; i < 4; i++)
	{
		ptmp1.coor[0] += 0.25*cp[tmface[fcid].cnct[i]].coor[0];
		ptmp1.coor[1] += 0.25*cp[tmface[fcid].cnct[i]].coor[1];
		ptmp1.coor[2] += 0.25*cp[tmface[fcid].cnct[i]].coor[2];
	}
	cp.push_back(ptmp1);
	pid[loc] = cp.size() - 1; loc++;
	for (i = 0; i<4; i++)//4 edges
	{
		int itmp(tmface[fcid].edge[i]);
		int dir(0);
		if (tmedge[itmp].pt[1] == tmface[fcid].cnct[i]) dir = 1;

		if (tmedge[itmp].act == 1)
		{
			EdgeSubdivision_2(itmp);//add one point and two edges
			pid[loc] = cp.size() - 1; loc++;
			if (dir == 0)
			{
				edid[edloc] = tmedge.size() - 2; edloc++;
				edid[edloc] = tmedge.size() - 1; edloc++;
			}
			else
			{
				edid[edloc] = tmedge.size() - 1; edloc++;
				edid[edloc] = tmedge.size() - 2; edloc++;
			}
		}
		else
		{
			pid[loc] = tmedge[itmp].midpt; loc++;
			int ied(tmedge[itmp].chd[0]);
			if (tmedge[ied].pt[0] == tmface[fcid].cnct[i] || tmedge[ied].pt[1] == tmface[fcid].cnct[i])
			{
				edid[edloc] = ied; edloc++;
				edid[edloc] = tmedge[itmp].chd[1]; edloc++;
			}
			else
			{
				edid[edloc] = tmedge[itmp].chd[1]; edloc++;
				edid[edloc] = ied; edloc++;
			}
		}
		//construct a new edge
		Edge3D edtmp;
		edtmp.pt[0] = pid[0]; edtmp.pt[1] = tmedge[itmp].midpt;
		edtmp.lev = tmedge[tmface[fcid].edge[(i + 1) % 4]].lev + 1;
		edtmp.len = tmedge[tmface[fcid].edge[(i + 1) % 4]].len / 2.;
		tmedge.push_back(edtmp);
		edid[edloc] = tmedge.size() - 1; edloc++;
	}

	int e_cnct[4][4] = { { tmface[fcid].cnct[0], pid[1], pid[0], pid[4] }, { pid[1], tmface[fcid].cnct[1], pid[2], pid[0] },
	{ pid[0], pid[2], tmface[fcid].cnct[2], pid[3] }, { pid[4], pid[0], pid[3], tmface[fcid].cnct[3] } };
	int e_edge[4][4] = { { edid[0], edid[2], edid[11], edid[10] }, { edid[1], edid[3], edid[5], edid[2] },
	{ edid[5], edid[4], edid[6], edid[8] }, { edid[11], edid[8], edid[7], edid[9] } };
	int enewid[4];
	double chd_org[4][2] = { { 0., 0. }, { tmedge[tmface[fcid].edge[0]].len / 2., 0. },
	{ tmedge[tmface[fcid].edge[0]].len / 2., tmedge[tmface[fcid].edge[1]].len / 2. }, { 0., tmedge[tmface[fcid].edge[1]].len / 2. } };
	int Tedge_id[4] = { 2, 5, 8, 11 };
	vector<Face3D> etmp(4);
	for (i = 0; i<4; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 0;
		etmp[i].lev = tmface[fcid].lev + 2;
		etmp[i].prt = fcid;
		for (int j = 0; j<4; j++)
		{
			etmp[i].cnct[j] = e_cnct[i][j];
			etmp[i].edge[j] = e_edge[i][j];
		}
		tmface.push_back(etmp[i]);
		enewid[i] = tmface.size() - 1;
		//tmface[fcid].chd.push_back(enewid[i]);
		//tmface[fcid].Tedge.push_back(edid[Tedge_id[i]]);
	}

	int edchd[2] = {5,11};
	int fcchd[2][2] = { {0,1}, {3,2} };
	if (fc_dir == 1)
	{
		edchd[0] = 2; edchd[1] = 8;
		fcchd[0][0] = 0; fcchd[0][1] = 3;
		fcchd[1][0] = 1; fcchd[1][1] = 2;
	}
	if (tmface[fcid].Tedge.size() == 1)
	{
		tmedge[tmface[fcid].Tedge[0]].act = 0;
		tmedge[tmface[fcid].Tedge[0]].midpt = pid[0];
		for (i = 0; i < 2; i++)
		{
			tmedge[tmface[fcid].Tedge[0]].chd[i] = edid[edchd[i]];
			tmedge[edid[edchd[i]]].prt = tmface[fcid].Tedge[0];
		}
	}
	else
	{
		cerr << "Configuration of face subdivisiono is wrong!\n";
	}
	tmface[fcid].Tedge.clear();

	int Tedge_id1[2] = {2,8};
	if (fc_dir == 1)
	{
		Tedge_id1[0] = 11; Tedge_id1[1] = 5;
	}
	if (tmface[fcid].chd.size() == 2)
	{
		for (i = 0; i < 2; i++)
		{
			tmface[tmface[fcid].chd[i]].act = 0;
			tmface[tmface[fcid].chd[i]].chd.push_back(enewid[fcchd[i][0]]);
			tmface[tmface[fcid].chd[i]].chd.push_back(enewid[fcchd[i][1]]);
			tmface[tmface[fcid].chd[i]].Tedge.clear();
			tmface[tmface[fcid].chd[i]].Tedge.push_back(edid[Tedge_id1[i]]);

			//if (edid[Tedge_id1[i]] == 347)
			//{
			//	cout << "here\n"; getchar();
			//}

			//tmface[enewid[fcchd[i][0]]].prt = tmface[fcid].chd[i];
			//tmface[enewid[fcchd[i][0]]].prt = tmface[fcid].chd[i];
		}
	}
	else
	{
		cerr << "Configuration of face subdivisiono is wrong!\n";
	}
	tmface[fcid].chd.clear();
	for (i = 0; i<4; i++)
	{
		tmface[fcid].chd.push_back(enewid[i]);
		tmface[fcid].Tedge.push_back(edid[Tedge_id[i]]);
	}

	tmface[fcid].ctpt = pid[0];
	tmface[fcid].act = 0;
}

void TruncatedTspline_3D::EdgeSubdivision_2(int edid)
{
	Vertex3D ptmp;
	ptmp.coor[0] = (cp[tmedge[edid].pt[0]].coor[0] + cp[tmedge[edid].pt[1]].coor[0]) / 2.;
	ptmp.coor[1] = (cp[tmedge[edid].pt[0]].coor[1] + cp[tmedge[edid].pt[1]].coor[1]) / 2.;
	ptmp.coor[2] = (cp[tmedge[edid].pt[0]].coor[2] + cp[tmedge[edid].pt[1]].coor[2]) / 2.;
	int mid(cp.size());
	cp.push_back(ptmp);
	tmedge[edid].midpt = mid;
	Edge3D edtmp1, edtmp2;
	int chdid[2];
	//1st child edge
	edtmp1.act = 1;
	edtmp1.lev = tmedge[edid].lev + 1;
	edtmp1.pt[0] = tmedge[edid].pt[0]; edtmp1.pt[1] = mid;
	edtmp1.len = tmedge[edid].len / 2.;
	edtmp1.prt = edid;
	chdid[0] = tmedge.size();
	tmedge[edid].chd[0] = chdid[0];
	tmedge.push_back(edtmp1);
	//2nd child edge
	edtmp2.act = 1;
	edtmp2.lev = tmedge[edid].lev + 1;
	edtmp2.pt[1] = tmedge[edid].pt[1]; edtmp2.pt[0] = mid;
	edtmp2.len = tmedge[edid].len / 2.;
	edtmp2.prt = edid;
	chdid[1] = tmedge.size();
	tmedge[edid].chd[1] = chdid[1];
	tmedge.push_back(edtmp2);

	tmedge[edid].act = 0;
}

void TruncatedTspline_3D::SolidFaceDirection(int eid, int fcid, int& dir, int& pos)//fcid is the local index of the face in the solid
{
	dir = 0; pos = 0;
	int itmp(tmesh[eid].face[fcid]);
	int* it = find(tmface[itmp].cnct, tmface[itmp].cnct + 4, tmesh[eid].cnct[solid_fc[fcid][0]]);
	pos=it - tmface[itmp].cnct;
	if (it == tmface[itmp].cnct + 4)
	{
		cerr << "Cannot find point in the face!\n";
		getchar();
	}
	else
	{
		if (tmesh[eid].cnct[solid_fc[fcid][1]] != tmface[itmp].cnct[(pos + 1) % 4])
		{
			dir = 1;
		}
	}
}

void TruncatedTspline_3D::EdgeIndex_in_Face(Face3D& fc, const vector<int>& edid)
{
	for (int i = 0; i < 4; i++)
	{
		int loc(-1);
		for (uint j = 0; j < edid.size(); j++)
		{
			if ((fc.cnct[i] == tmedge[edid[j]].pt[0] && fc.cnct[(i + 1) % 4] == tmedge[edid[j]].pt[1]) ||
				(fc.cnct[i] == tmedge[edid[j]].pt[1] && fc.cnct[(i + 1) % 4] == tmedge[edid[j]].pt[0]))
			{
				loc = j; break;
			}
		}
		if (loc == -1)
		{
			loc = 0;
			cout << "Cant't find edge in face\n";
			getchar();
		}
		fc.edge[i] = edid[loc];
	}
}

void TruncatedTspline_3D::EdgeFaceIndex_in_Solid(Element3D& hx, const vector<int>& edid, const vector<int>& fcid)
{
	for (int i = 0; i < 12; i++)
	{
		int loc(-1);
		for (uint j = 0; j < edid.size(); j++)
		{
			if ((hx.cnct[solid_ed[i][0]] == tmedge[edid[j]].pt[0] && hx.cnct[solid_ed[i][1]] == tmedge[edid[j]].pt[1]) ||
				(hx.cnct[solid_ed[i][0]] == tmedge[edid[j]].pt[1] && hx.cnct[solid_ed[i][1]] == tmedge[edid[j]].pt[0]))
			{
				loc = j; break;
			}
		}
		if (loc == -1)
		{
			loc = 0;
			cout << "Can't find edge in a solid\n";
			getchar();
		}
		hx.edge[i] = edid[loc];
	}
	for (int i = 0; i < 6; i++)
	{
		int loc(-1);
		array<int, 4> cnct1 = { hx.cnct[solid_fc[i][0]], hx.cnct[solid_fc[i][1]], hx.cnct[solid_fc[i][2]], hx.cnct[solid_fc[i][3]]};
		sort(cnct1.begin(),cnct1.end());
		for (uint j = 0; j < fcid.size(); j++)
		{
			array<int, 4> cnct2 = { tmface[fcid[j]].cnct[0], tmface[fcid[j]].cnct[1], tmface[fcid[j]].cnct[2], tmface[fcid[j]].cnct[3]};
			sort(cnct2.begin(), cnct2.end());
			if (cnct1==cnct2)
			{
				loc = j; break;
			}
		}
		if (loc == -1)
		{
			loc = 0;
			cout << "Can't find face in a solid\n";
			cout << hx.cnct[solid_fc[i][0]] << " " << hx.cnct[solid_fc[i][1]] << " " << hx.cnct[solid_fc[i][2]] << " " << hx.cnct[solid_fc[i][3]] << "\n";
			getchar();
		}
		hx.face[i] = fcid[loc];
	}
}

void TruncatedTspline_3D::Index_Direction(int eid, int vloc[8], int edloc[12], int fcloc[6], int dir)
{
	int i;
	if (dir == 0)//x direction
	{
		for (i = 0; i < 8; i++) vloc[i] = tmesh[eid].cnct[i];
		for (i = 0; i < 12; i++) vloc[i] = tmesh[eid].edge[i];
		for (i = 0; i < 6; i++) vloc[i] = tmesh[eid].face[i];
	}
	else if (dir == 1)//y
	{
	}
	else//z
	{
	}
}

void TruncatedTspline_3D::Construct_NewEdge(int pt[2], int lev, double len)
{
	Edge3D edtmp;
	edtmp.pt[0] = pt[0]; edtmp.pt[1] = pt[1];
	edtmp.lev = lev;
	edtmp.len = len;
	tmedge.push_back(edtmp);
}

void TruncatedTspline_3D::Construct_BaseFace_Subdv(int fcid, int dir, int& ed_base, int fc_base[2])
{
	int edid[7];
	int iloc[2] = {dir,dir+2};
	for (int i = 0; i < 2; i++)
	{
		edid[2 * i] = tmedge[tmface[fcid].edge[iloc[i]]].chd[0];
		edid[2 * i+1] = tmedge[tmface[fcid].edge[iloc[i]]].chd[1];
		if (tmedge[edid[2 * i]].pt[0] != tmface[fcid].cnct[iloc[i]] && tmedge[edid[2 * i]].pt[1] != tmface[fcid].cnct[iloc[i]])
		{
			edid[2 * i] = tmedge[tmface[fcid].edge[iloc[i]]].chd[1];
			edid[2 * i + 1] = tmedge[tmface[fcid].edge[iloc[i]]].chd[0];
		}
	}
	edid[4] = tmface[fcid].edge[dir+1];
	edid[5] = tmface[fcid].edge[(dir + 3)%4];
	//construct base edge
	Edge3D etmp;
	etmp.act = 0;
	etmp.pt[0] = tmedge[tmface[fcid].edge[dir]].midpt;
	etmp.pt[1] = tmedge[tmface[fcid].edge[dir + 2]].midpt;
	etmp.lev = tmedge[tmface[fcid].edge[dir+1]].lev;
	etmp.len = tmedge[tmface[fcid].edge[dir + 1]].len;
	etmp.midpt = tmface[fcid].ctpt;
	etmp.chd[0] = tmface[fcid].Tedge[dir];
	etmp.chd[1] = tmface[fcid].Tedge[dir + 2];
	tmedge.push_back(etmp);
	edid[6] = tmedge.size() - 1;
	tmedge[tmface[fcid].Tedge[dir]].prt = edid[6];
	tmedge[tmface[fcid].Tedge[dir + 2]].prt = edid[6];

	//construct base face
	int fcnew_cnct[2][4] = { { tmface[fcid].cnct[dir], etmp.pt[0], etmp.pt[1], tmface[fcid].cnct[(dir + 3) % 4] }, 
	{ etmp.pt[0], tmface[fcid].cnct[dir + 1], tmface[fcid].cnct[dir + 2], etmp.pt[1] } };
	int fcnew_edge[2][4] = { { edid[0], edid[6], edid[3], edid[5] }, { edid[1], edid[4], edid[2], edid[6] }};
	int fc_chd[2][2] = { { tmface[fcid].chd[dir], tmface[fcid].chd[(dir + 3) % 4] }, { tmface[fcid].chd[dir + 1], tmface[fcid].chd[dir + 2] } };
	int fc_Ted[2] = { tmface[fcid].Tedge[(dir+3)%4], tmface[fcid].Tedge[dir+1] };
	for (int i = 0; i < 2; i++)
	{
		Face3D fctmp;
		fctmp.act = 0;
		fctmp.lev = tmface[fcid].lev - 1;
		for (int j = 0; j < 4; j++)
		{
			fctmp.cnct[j] = fcnew_cnct[i][j];
			fctmp.edge[j] = fcnew_edge[i][j];
		}
		fctmp.chd.push_back(fc_chd[i][0]);
		fctmp.chd.push_back(fc_chd[i][1]);
		fctmp.Tedge.push_back(fc_Ted[i]);
		tmface.push_back(fctmp);
		tmface[fc_chd[i][0]].prt = tmface.size() - 1;
		tmface[fc_chd[i][1]].prt = tmface.size() - 1;

		//if (fc_Ted[i] == 347)
		//{
		//	cout << "here\n"; getchar();
		//}
	}

	ed_base = edid[6];
	fc_base[0] = tmface.size() - 2;
	fc_base[1] = tmface.size() - 1;
}

//void TruncatedTspline_3D::Truncation()
//{
//	int diag[8] = { 6, 7, 4, 5, 2, 3, 0, 1 };
//	for (uint i = 0; i < tmesh.size(); i++)
//	{
//		if (tmesh[i].act == 1 && tmesh[i].type == 0)
//		{
//			vector<int> fcnb;
//			for (uint j = 0; j < 6; j++)
//			{
//				if (tmface[tmesh[i].face[j]].act == 1)
//				{
//					int tmp(tmface[tmesh[i].face[j]].hex[0]);
//					if (tmp == i) tmp = tmface[tmesh[i].face[j]].hex[1];
//					fcnb.push_back(tmp);
//				}
//				else
//				{
//					for (uint k = 0; k < tmface[tmesh[i].face[j]].chd.size(); k++)
//					{
//						int tmp(tmface[tmface[tmesh[i].face[j]].chd[k]].hex[0]);
//						if (tmp == i) tmp = tmface[tmface[tmesh[i].face[j]].chd[k]].hex[1];
//						fcnb.push_back(tmp);
//					}
//				}
//			}
//			int dc(-1), crnb(-1);
//			for (int j = 0; j < 8; j++)
//			{
//				for (uint k = 0; k < cp[tmesh[i].cnct[j]].hex.size(); k++)
//				{
//					if (tmesh[cp[tmesh[i].cnct[j]].hex[k]].lev == tmesh[i].lev - 3)
//					{
//						vector<int>::iterator it = find(fcnb.begin(), fcnb.end(), cp[tmesh[i].cnct[j]].hex[k]);
//						int hxid(cp[tmesh[i].cnct[j]].hex[k]);
//						int* it1 = find(tmesh[hxid].cnct, tmesh[hxid].cnct + 8, tmesh[i].cnct[j]);
//						if (it == fcnb.end() && it1 != tmesh[hxid].cnct + 8)
//						{
//							crnb = cp[tmesh[i].cnct[j]].hex[k];
//							dc = tmesh[i].cnct[j];
//							break;
//						}
//					}
//				}
//				if (dc != -1) break;
//			}
//			if (dc!=-1)
//			{
//				int* it = find(tmesh[crnb].cnct, tmesh[crnb].cnct+8, dc);
//				if (it != tmesh[crnb].cnct + 8)
//				{
//					int pos(it - tmesh[crnb].cnct);
//					int trid(tmesh[crnb].cnct[diag[pos]]);
//					cout << "trid: " << trid << "\n";
//					cout << "dc: " << dc << "\n";
//					getchar();
//					vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), dc);
//					int pos1(it1 - tmesh[i].IEN.begin());
//					vector<int>::iterator it2 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), trid);
//					int pos2(it2 - tmesh[i].IEN.begin());
//					if (CheckSubKnotVector(tmesh[i].patch_ku[pos1], tmesh[i].patch_kv[pos1], tmesh[i].patch_kw[pos1], tmesh[i].patch_ku[pos2], tmesh[i].patch_kv[pos2], tmesh[i].patch_kw[pos2]))
//					{
//					vector<double> ku1(10), kv1(10), kw1(10);
//					vector<vector<double>> Tu, Tv, Tw;
//					vector<double> ku(tmesh[i].patch_ku[pos2].begin(), tmesh[i].patch_ku[pos2].end());
//					vector<double> kv(tmesh[i].patch_kv[pos2].begin(), tmesh[i].patch_kv[pos2].end());
//					vector<double> kw(tmesh[i].patch_kw[pos2].begin(), tmesh[i].patch_kw[pos2].end());
//					vector<double>::iterator it1, it2, it3;
//					it1 = set_union(tmesh[i].patch_ku[pos2].begin(), tmesh[i].patch_ku[pos2].end(), tmesh[i].patch_ku[pos1].begin(), tmesh[i].patch_ku[pos1].end(), ku1.begin());
//					it2 = set_union(tmesh[i].patch_kv[pos2].begin(), tmesh[i].patch_kv[pos2].end(), tmesh[i].patch_kv[pos1].begin(), tmesh[i].patch_kv[pos1].end(), kv1.begin());
//					it3 = set_union(tmesh[i].patch_kw[pos2].begin(), tmesh[i].patch_kw[pos2].end(), tmesh[i].patch_kw[pos1].begin(), tmesh[i].patch_kw[pos1].end(), kw1.begin());
//					ku1.resize(it1 - ku1.begin());
//					kv1.resize(it2 - kv1.begin());
//					kw1.resize(it3 - kw1.begin());
//					TMatrix(ku, ku1, 3, Tu);
//					TMatrix(kv, kv1, 3, Tv);
//					TMatrix(kw, kw1, 3, Tw);
//					it1 = search(ku1.begin(), ku1.end(), tmesh[i].patch_ku[pos1].begin(), tmesh[i].patch_ku[pos1].end());
//					it2 = search(kv1.begin(), kv1.end(), tmesh[i].patch_kv[pos1].begin(), tmesh[i].patch_kv[pos1].end());
//					it3 = search(kw1.begin(), kw1.end(), tmesh[i].patch_kw[pos1].begin(), tmesh[i].patch_kw[pos1].end());
//					if (it1 != ku1.end() && it2 != kv1.end() && it3 != kw1.end())
//					{
//					int loc1 = it1 - ku1.begin();
//					int loc2 = it2 - kv1.begin();
//					int loc3 = it3 - kw1.begin();
//					double tc = Tu[loc1][0] * Tv[loc2][0] * Tw[loc3][0];
//					cp[trid].trun = 1;
//					cp[trid].tbf.push_back(dc);
//					cp[trid].tc.push_back(tc);
//					}
//					}
//				}
//			}
//		}
//	}
//}

void TruncatedTspline_3D::Truncation()
{
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].type == 0)
		{
			for (uint j = 0; j < tmesh[i].IEN.size(); j++)//trun
			{
				for (uint k = 0; k < tmesh[i].IEN.size(); k++)//child
				{
					if (j != k)
					{
						int pos1(k), pos2(j);
						if (CheckSubKnotVector(tmesh[i].patch_ku[pos1], tmesh[i].patch_kv[pos1], tmesh[i].patch_kw[pos1], tmesh[i].patch_ku[pos2], tmesh[i].patch_kv[pos2], tmesh[i].patch_kw[pos2]))
						{
							vector<double> ku1(10), kv1(10), kw1(10);
							vector<vector<double>> Tu, Tv, Tw;
							vector<double> ku(tmesh[i].patch_ku[pos2].begin(), tmesh[i].patch_ku[pos2].end());
							vector<double> kv(tmesh[i].patch_kv[pos2].begin(), tmesh[i].patch_kv[pos2].end());
							vector<double> kw(tmesh[i].patch_kw[pos2].begin(), tmesh[i].patch_kw[pos2].end());
							vector<double>::iterator it1, it2, it3;
							it1 = set_union(tmesh[i].patch_ku[pos2].begin(), tmesh[i].patch_ku[pos2].end(), tmesh[i].patch_ku[pos1].begin(), tmesh[i].patch_ku[pos1].end(), ku1.begin());
							it2 = set_union(tmesh[i].patch_kv[pos2].begin(), tmesh[i].patch_kv[pos2].end(), tmesh[i].patch_kv[pos1].begin(), tmesh[i].patch_kv[pos1].end(), kv1.begin());
							it3 = set_union(tmesh[i].patch_kw[pos2].begin(), tmesh[i].patch_kw[pos2].end(), tmesh[i].patch_kw[pos1].begin(), tmesh[i].patch_kw[pos1].end(), kw1.begin());
							ku1.resize(it1 - ku1.begin());
							kv1.resize(it2 - kv1.begin());
							kw1.resize(it3 - kw1.begin());
							TMatrix(ku, ku1, 3, Tu);
							TMatrix(kv, kv1, 3, Tv);
							TMatrix(kw, kw1, 3, Tw);
							it1 = search(ku1.begin(), ku1.end(), tmesh[i].patch_ku[pos1].begin(), tmesh[i].patch_ku[pos1].end());
							it2 = search(kv1.begin(), kv1.end(), tmesh[i].patch_kv[pos1].begin(), tmesh[i].patch_kv[pos1].end());
							it3 = search(kw1.begin(), kw1.end(), tmesh[i].patch_kw[pos1].begin(), tmesh[i].patch_kw[pos1].end());
							if (it1 != ku1.end() && it2 != kv1.end() && it3 != kw1.end())
							{
								int loc1 = it1 - ku1.begin();
								int loc2 = it2 - kv1.begin();
								int loc3 = it3 - kw1.begin();
								double tc = Tu[loc1][0] * Tv[loc2][0] * Tw[loc3][0];
								if (tc != 0.)
								{
									cp[tmesh[i].IEN[j]].trun = 1;
									vector<int>::iterator it = find(cp[tmesh[i].IEN[j]].tbf.begin(), cp[tmesh[i].IEN[j]].tbf.end(), tmesh[i].IEN[k]);
									if (it == cp[tmesh[i].IEN[j]].tbf.end())
									{
										cp[tmesh[i].IEN[j]].tbf.push_back(tmesh[i].IEN[k]);
										cp[tmesh[i].IEN[j]].tc.push_back(tc);
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void TruncatedTspline_3D::IdentifyTarget(vector<int>& target)
{
	target.clear();
	//int rfid(-1), pid(2 * 25 + 2 * 5 + 2);//pid(124), pid(2*25+2*5+2)
	//for (uint i = 0; i < tmesh.size(); i++)
	//{
	//	for (int j = 0; j < 8; j++)
	//	{
	//		if (tmesh[i].cnct[j] == pid)
	//		{
	//			rfid = i; break;
	//		}
	//		if (rfid != -1) break;
	//	}
	//}
	//target.push_back(rfid);
	//for (uint i = 0; i < cp[pid].hex.size(); i++)
	//{
	//	int flag(0);
	//	for (uint j = 0; j < 12; j++)
	//	{
	//		for (uint k = 0; k < tmedge[tmesh[cp[pid].hex[i]].edge[j]].hex.size(); k++)
	//		{
	//			if (tmedge[tmesh[cp[pid].hex[i]].edge[j]].hex[k] == rfid)//share an edge with rfid
	//			{
	//				flag = 1; break;
	//			}
	//		}
	//	}
	//	if (flag == 0)
	//	{
	//		target.push_back(cp[pid].hex[i]); break;
	//	}
	//}

	//for (uint i = 0; i < target.size(); i++)
	//{
	//	//ElementSubdivide_8(target[i]);
	//	cout << target[i] << "\n";
	//}
	//getchar();

	//cube5-targets
	int nh(9);
	for (int i = 1; i < nh-1; i++)
	{
		target.push_back(i*nh*nh+i*nh+i);
	}
}

void TruncatedTspline_3D::IdentifyAddition(vector<int>& target)
{
	uint i, j, k;
	for (i = 0; i < tmesh.size(); i++) tmesh[i].ref_flag = 0;
	for (i = 0; i < target.size(); i++) tmesh[target[i]].ref_flag = 1;//subdv into 8
	for (i = 0; i < cp.size(); i++)
	{
		if (cp[i].type == 0)//interior regular node
		{
			int nref(0);
			for (j = 0; j < cp[i].hex.size(); j++)
			{
				if (tmesh[cp[i].hex[j]].ref_flag == 1) nref++;
			}
			if (nref>0 && nref < 8)
			{
				int ref3d(0);
				for (j = 0; j < cp[i].edge.size(); j++)
				{
					int nref_ed(0);
					for (k = 0; k < tmedge[cp[i].edge[j]].hex.size(); k++)
					{
						if (tmesh[tmedge[cp[i].edge[j]].hex[k]].ref_flag == 1) nref_ed++;
					}
					if (nref_ed>0 && nref_ed < nref)//3D case
					{
						ref3d = 1;
					}
				}
				if (ref3d == 1)
				{
					for (j = 0; j < cp[i].hex.size(); j++)
					{
						if (tmesh[cp[i].hex[j]].ref_flag != 1) tmesh[cp[i].hex[j]].ref_flag = 21;//tmp
					}
				}
			}
		}
	}

	//int sub4[4] = {22,25,38,41};
	//int sub2[2] = { 26,37 };
	//for (i = 0; i < 4; i++) tmesh[sub4[i]].ref_flag = 23;
	//for (i = 0; i < 2; i++) tmesh[sub2[i]].ref_flag = 31;

	//int sub8[2] = { 1 * 5 * 5 + 1 * 5 + 1, 3 * 5 * 5 + 3 * 5 + 3};

	//tmesh[sub8[0]].ref_flag = 1;
	//tmesh[sub8[1]].ref_flag = 1;
}

bool TruncatedTspline_3D::IdentifyAddition()
{
	uint i, j, k;
	vector<int> r1;
	int flag(0);
	for (i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].ref_flag==1)
		{
			for (j = 0; j < 8; j++)
			{
				for (k = 0; k < cp[tmesh[i].cnct[j]].hex.size(); k++)
				{
					int hxid(cp[tmesh[i].cnct[j]].hex[k]);
					if (hxid != i && tmesh[hxid].ref_flag==0)
					{
						vector<int>::iterator it2 = find(r1.begin(), r1.end(), hxid);
						if (it2 == r1.end())
						{
							r1.push_back(hxid);
						}
					}
				}
			}
		}
	}

	for (i = 0; i < r1.size(); i++)
	{
		int flag_fc(0), flag_ed(0), flag_vt(0);
		vector<int> fcnb, ednb, vtnb;
		for (j = 0; j < 6; j++)
		{
			int fcid(tmesh[r1[i]].face[j]);
			if (tmface[fcid].hex.size() == 2)
			{
				int tmp1(tmface[fcid].hex[0]);
				if (tmp1 == r1[i]) tmp1 = tmface[fcid].hex[1];
				fcnb.push_back(tmp1);
			}
		}
		for (j = 0; j < 12; j++)
		{
			int edid(tmesh[r1[i]].edge[j]);
			if (tmedge[edid].hex.size() == 4)
			{
				for (k = 0; k < tmedge[edid].hex.size(); k++)
				{
					if (tmedge[edid].hex[k] != r1[i])
					{
						vector<int>::iterator it = find(fcnb.begin(), fcnb.end(), tmedge[edid].hex[k]);
						if (it == fcnb.end())
						{
							ednb.push_back(tmedge[edid].hex[k]);
							break;
						}
					}
				}
			}
		}
		for (j = 0; j < 8; j++)
		{
			int pid(tmesh[r1[i]].cnct[j]);
			if (cp[pid].hex.size() == 8)
			{
				for (k = 0; k < cp[pid].hex.size(); k++)
				{
					if (cp[pid].hex[k] != r1[i])
					{
						vector<int>::iterator it1 = find(fcnb.begin(), fcnb.end(), cp[pid].hex[k]);
						vector<int>::iterator it2 = find(ednb.begin(), ednb.end(), cp[pid].hex[k]);
						if (it1 == fcnb.end() && it2 == ednb.end())
						{
							vtnb.push_back(cp[pid].hex[k]);
							break;
						}
					}
				}
			}
		}

		for (j = 0; j < fcnb.size(); j++)
		{
			if (tmesh[fcnb[j]].ref_flag == 1) flag_fc++;
		}
		for (j = 0; j < ednb.size(); j++)
		{
			if (tmesh[ednb[j]].ref_flag == 1) flag_ed++;
		}
		for (j = 0; j < vtnb.size(); j++)
		{
			if (tmesh[vtnb[j]].ref_flag == 1) flag_vt++;
		}

		if (flag_fc == 1 && flag_ed==0 && flag_vt != 0)
		{
			//no finished
			tmesh[r1[i]].ref_flag = 1;
			flag = 1;
		}
		else if (flag_fc == 0 && (flag_vt != 0 || flag_ed!=0))
		{
			tmesh[r1[i]].ref_flag = 1;
			flag = 1;
		}
	}

	if (flag == 1) return true;
	else return false;
}

bool TruncatedTspline_3D::Identify_FaceFace(int eid)//apply this one first until no such case
{
	int nfc(0);
	for (int i = 0; i < 6; i++)
	{
		if (tmface[tmesh[eid].face[i]].ctpt != -1) nfc++;
	}
	if (nfc < 2) return false;
	else if (nfc == 2)
	{
		if (tmface[tmesh[eid].face[0]].act == 0 && tmface[tmesh[eid].face[5]].act == 0)//opposite pair, split into 4
		{
			tmesh[eid].ref_flag = 23;//w direction
		}
		else if (tmface[tmesh[eid].face[1]].act == 0 && tmface[tmesh[eid].face[3]].act == 0)//opposite pair, split into 4
		{
			tmesh[eid].ref_flag = 22;//v direction
		}
		else if (tmface[tmesh[eid].face[2]].act == 0 && tmface[tmesh[eid].face[4]].act == 0)//opposite pair, split into 4
		{
			tmesh[eid].ref_flag = 21;//u direction
		}
		else if ((tmface[tmesh[eid].face[0]].act == 0 || tmface[tmesh[eid].face[5]].act == 0) && (tmface[tmesh[eid].face[1]].act == 0 || tmface[tmesh[eid].face[3]].act == 0))
		{
			tmesh[eid].ref_flag = 32;//v direction
		}
		else if ((tmface[tmesh[eid].face[0]].act == 0 || tmface[tmesh[eid].face[5]].act == 0) && (tmface[tmesh[eid].face[2]].act == 0 || tmface[tmesh[eid].face[4]].act == 0))
		{
			tmesh[eid].ref_flag = 31;//u direction
		}
		else if ((tmface[tmesh[eid].face[1]].act == 0 || tmface[tmesh[eid].face[3]].act == 0) && (tmface[tmesh[eid].face[2]].act == 0 || tmface[tmesh[eid].face[4]].act == 0))
		{
			tmesh[eid].ref_flag = 33;//v direction
		}
		return true;
	}
	else if (nfc == 3 || nfc == 4)
	{
		if (tmface[tmesh[eid].face[0]].act == 0 && tmface[tmesh[eid].face[5]].act == 0)//opposite pair, split into 4
		{
			tmesh[eid].ref_flag = 23;//w direction
		}
		else if (tmface[tmesh[eid].face[1]].act == 0 && tmface[tmesh[eid].face[3]].act == 0)//opposite pair, split into 4
		{
			tmesh[eid].ref_flag = 22;//v direction
		}
		else if (tmface[tmesh[eid].face[2]].act == 0 && tmface[tmesh[eid].face[4]].act == 0)//opposite pair, split into 4
		{
			tmesh[eid].ref_flag = 21;//u direction
		}
		else
		{
			tmesh[eid].ref_flag = 21;//u direction
		}
		return true;
	}
	else if (nfc > 4)
	{
		tmesh[eid].ref_flag = 1;//subdivide into 8
		return true;
	}
}

bool TruncatedTspline_3D::Identify_FaceEdge(int eid)
{
	int nfc(0), pos(0);
	for (int i = 0; i < 6; i++)
	{
		if (tmface[tmesh[eid].face[i]].ctpt != -1)
		{
			pos = i; nfc++;
		}
	}
	if (nfc != 1) return false;

	int edppd[6][4] = { { 4, 5, 6, 7 }, { 1, 3, 9, 11 }, { 0, 2, 8, 10 }, { 1, 3, 9, 11 }, { 0, 2, 8, 10 }, { 4, 5, 6, 7 } };
	int edopo[6][4] = { { 8, 9, 10, 11 }, { 2, 6, 10, 7 }, { 3, 7, 11, 4 }, { 0, 5, 8, 4 }, { 1, 6, 9, 5 }, {0,1,2,3} };
	int opo_dir[6][4] = { { 1, 2, 1, 2 }, { 1, 3, 1, 3 }, { 2, 3, 2, 3 }, { 1, 3, 1, 3 }, { 2, 3, 2, 3 }, {1,2,1,2} };
	int ref_dir[3] = {0,0,0};
	int flag_ppd(0), flag_opo(0);
	for (int i = 0; i < 4; i++)
	{
		if (tmedge[tmesh[eid].edge[edppd[pos][i]]].act == 0) flag_ppd++;
		if (tmedge[tmesh[eid].edge[edopo[pos][i]]].act == 0) flag_opo++;
	}

	if (flag_ppd != 0) 
	{
		if (pos == 0 || pos == 5) ref_dir[2] = 1;
		else if (pos == 1 || pos == 3) ref_dir[1] = 1;
		else ref_dir[0] = 1;
	}
	if (flag_opo == 1)
	{
		if (tmedge[tmesh[eid].edge[edopo[pos][0]]].act == 0 || tmedge[tmesh[eid].edge[edopo[pos][2]]].act == 0)
		{
			ref_dir[opo_dir[pos][0]-1] = 1;
		}
		else
		{
			ref_dir[opo_dir[pos][1] - 1] = 1;
		}
	}
	else if (flag_opo == 2)
	{
		if (tmedge[tmesh[eid].edge[edopo[pos][0]]].act == 0 && tmedge[tmesh[eid].edge[edopo[pos][2]]].act == 0)
		{
			ref_dir[opo_dir[pos][0] - 1] = 1;
		}
		else if (tmedge[tmesh[eid].edge[edopo[pos][1]]].act == 0 && tmedge[tmesh[eid].edge[edopo[pos][3]]].act == 0)
		{
			ref_dir[opo_dir[pos][1] - 1] = 1;
		}
		else
		{
			ref_dir[opo_dir[pos][0] - 1] = 1; ref_dir[opo_dir[pos][1] - 1] = 1;
		}
	}
	else if (flag_opo == 3 || flag_opo == 4)
	{
		ref_dir[opo_dir[pos][0] - 1] = 1; ref_dir[opo_dir[pos][1] - 1] = 1;
	}

	if (ref_dir[0] == 1 && ref_dir[1] == 1 && ref_dir[2] == 1)
	{
		tmesh[eid].ref_flag = 1;
		return true;
	}
	else if (ref_dir[0] == 0 && ref_dir[1] == 1 && ref_dir[2] == 1)
	{
		tmesh[eid].ref_flag = 21;
		return true;
	}
	else if (ref_dir[0] == 1 && ref_dir[1] == 0 && ref_dir[2] == 1)
	{
		tmesh[eid].ref_flag = 22;
		return true;
	}
	else if (ref_dir[0] == 1 && ref_dir[1] == 1 && ref_dir[2] == 0)
	{
		tmesh[eid].ref_flag = 23;
		return true;
	}
	else if (ref_dir[0] == 1 && ref_dir[1] == 0 && ref_dir[2] == 0)
	{
		tmesh[eid].ref_flag = 31;
		return true;
	}
	else if (ref_dir[0] == 0 && ref_dir[1] == 1 && ref_dir[2] == 0)
	{
		tmesh[eid].ref_flag = 32;
		return true;
	}
	else if (ref_dir[0] == 0 && ref_dir[1] == 0 && ref_dir[2] == 1)
	{
		tmesh[eid].ref_flag = 33;
		return true;
	}
	else
	{
		return false;
	}

	//if (flag_ppd != 0 && flag_opo == 0)
	//{
	//	if (pos == 0 || pos == 5) tmesh[eid].ref_flag = 33;//w direction
	//	else if (pos == 1 || pos == 3) tmesh[eid].ref_flag = 32;//v direction
	//	else if (pos == 2 || pos == 4) tmesh[eid].ref_flag = 31;//u direction
	//	return true;
	//}
	//else if (flag_ppd == 0 && flag_opo != 0)
	//{

	//}
	//else if (flag_ppd != 0 && flag_opo != 0)
	//{
	//}
	//else
	//{
	//	return false;
	//}

	//int fcid(tmesh[eid].face[pos]);
	//for (int i = 0; i < 4; i++)
	//{
	//	for (int j = 0; j < 4; j++)
	//	{
	//		int edid(tmface[tmesh[eid].face[edppd[pos][i]]].edge[j]);
	//		int* it = find(tmface[fcid].edge, tmface[fcid].edge+4, edid);
	//		int edpos = it - tmface[fcid].edge;

	//	}
	//}
}

bool TruncatedTspline_3D::Template_FaceFace()
{
	bool flag(false);
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].ref_flag == 0 && Identify_FaceFace(i))
		{
			flag = true;
			if (tmesh[i].ref_flag == 1) ElementSubdivide_8(i);
			else if (tmesh[i].ref_flag == 21) ElementSubdivide_4(i, 0);
			else if (tmesh[i].ref_flag == 22) ElementSubdivide_4(i, 1);
			else if (tmesh[i].ref_flag == 23) ElementSubdivide_4(i, 2);
			else if (tmesh[i].ref_flag == 31) ElementSubdivide_2(i, 0);
			else if (tmesh[i].ref_flag == 32) ElementSubdivide_2(i, 1);
			else if (tmesh[i].ref_flag == 33) ElementSubdivide_2(i, 2);
		}
	}
	return flag;
}

bool TruncatedTspline_3D::Template_FaceEdge()
{
	bool flag(false);
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].ref_flag == 0 && Identify_FaceEdge(i))
		{
			flag = true;
			if (tmesh[i].ref_flag == 1) ElementSubdivide_8(i);
			else if (tmesh[i].ref_flag == 21) ElementSubdivide_4(i, 0);
			else if (tmesh[i].ref_flag == 22) ElementSubdivide_4(i, 1);
			else if (tmesh[i].ref_flag == 23) ElementSubdivide_4(i, 2);
			else if (tmesh[i].ref_flag == 31) ElementSubdivide_2(i, 0);
			else if (tmesh[i].ref_flag == 32) ElementSubdivide_2(i, 1);
			else if (tmesh[i].ref_flag == 33) ElementSubdivide_2(i, 2);
		}
	}
	return flag;
}

void TruncatedTspline_3D::RefineTopology()
{
	npt_old = cp.size();
	nel_old = tmesh.size();
	vector<int> target;
	IdentifyTarget(target);
	IdentifyAddition(target);
	//int nloop(0), loop_max(50);
	//while (Template_FaceFace() && nloop<loop_max)
	//{
	//	nloop++;
	//}
	//cout << "# of loops for face-face template: " << nloop << "\n";
	////getchar();
	//nloop = 0;
	//while (Template_FaceEdge() && nloop<loop_max)
	//{
	//	nloop++;
	//}
	//cout << "# of loops for face-face template: " << nloop << "\n";
	////getchar();

	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].ref_flag == 1) ElementSubdivide_8(i);
		else if (tmesh[i].ref_flag == 21) ElementSubdivide_4(i, 0);
		else if (tmesh[i].ref_flag == 22) ElementSubdivide_4(i, 1);
		else if (tmesh[i].ref_flag == 23) ElementSubdivide_4(i, 2);
		else if (tmesh[i].ref_flag == 31) ElementSubdivide_2(i,0);
		else if (tmesh[i].ref_flag == 32) ElementSubdivide_2(i, 1);
		else if (tmesh[i].ref_flag == 33) ElementSubdivide_2(i, 2);
	}
}

void TruncatedTspline_3D::RefineGeometry()
{
	CollectActives();
	UpdateConnect();
	FindEdgeTopoDirec();
	FindKnotInterval();

	//int kv_mod[2] = {42,82};
	//for (int i = 0; i < 2; i++)
	//{
	//	for (int j = 0; j < 4; j++)
	//	{
	//		cp[kv_mod[i]].kitvU[j] = 1.;
	//		cp[kv_mod[i]].kitvV[j] = 1.;
	//		cp[kv_mod[i]].kitvW[j] = 1.;
	//	}
	//}

	//int edid(403);
	//cout << "neighbor faces: ";
	//for (uint i = 0; i < tmedge[edid].face.size(); i++)
	//{
	//	cout << tmface[tmedge[edid].face[i]].id_act << " ";
	//}
	//cout << "\n";
	//getchar();
	//cout << "neighbor hex:\n";
	//for (uint i = 0; i < tmedge[edid].hex.size(); i++)
	//{
	//	cout << "hex id: "<< tmesh[tmedge[edid].hex[i]].id_act << "\n";
	//	for (int j = 0; j < 6; j++)
	//	{
	//		cout << tmface[tmesh[tmedge[edid].hex[i]].face[j]].id_act << "(" << tmface[tmesh[tmedge[edid].hex[i]].face[j]].act << ")  ";
	//	}
	//	cout << "\n";
	//}
	//cout << "\n";
	//getchar();
	//cout << tmesh[75].id_act << " " << tmesh[78].id_act << "\n";
	//getchar();

	//for (uint i = 0; i < cp.size(); i++)
	//{
	//	if (i==159)
	//	{
	//		cout << "pid: " << i << "\n";
	//		cout << "uv edges: " << cp[i].uved[0] << " " << cp[i].uved[1] << " " << cp[i].uved[2] << "\n";
	//		cout << "uv edges: " << tmedge[cp[i].uved[0]].id_act << " " << tmedge[cp[i].uved[1]].id_act << " " << tmedge[cp[i].uved[2]].id_act << "\n";
	//		cout << "pt[0] "<< tmedge[cp[i].uved[0]].pt[0] << ": " << tmedge[cp[i].uved[0]].pn[0][0] << "(type), " << tmedge[cp[i].uved[0]].pn[0][1] << "(ID)\n";
	//		cout << "pt[1] " << tmedge[cp[i].uved[0]].pt[1] << ": " << tmedge[cp[i].uved[0]].pn[1][0] << "(type), " << tmedge[cp[i].uved[0]].pn[1][1] << "(ID)\n";
	//		cout << "kitU: ";
	//		for (int j = 0; j < 4; j++) cout << cp[i].kitvU[j] << " ";
	//		cout << "\n";
	//		cout << "kitV: ";
	//		for (int j = 0; j < 4; j++) cout << cp[i].kitvV[j] << " ";
	//		cout << "\n";
	//		cout << "kitW: ";
	//		for (int j = 0; j < 4; j++) cout << cp[i].kitvW[j] << " ";
	//		cout << "\n";
	//		getchar();
	//	}
	//}

	SetLocalCoorSystem();
	FindIEN_PatchKV();
	//Update_IEN();

	Truncation();

	//int pid(33);
	//cout << "pid " << pid << " trun: " << cp[pid].trun << "\n";
	//for (uint i = 0; i < cp[pid].tbf.size(); i++)
	//{
	//	cout << cp[pid].tbf[i] << " ";
	//}
	//cout << "\n";
	//getchar();
}

void TruncatedTspline_3D::RefineGeometry_v1()
{
	CollectActives();
	UpdateConnect();
	FindEdgeTopoDirec();
	FindKnotInterval();

	SetLocalCoorSystem();
	FindIEN_PatchKV();

	for (uint i = 0; i < cp.size(); i++)
	{
		cp[i].coortmp[0] = 0.; cp[i].coortmp[1] = 0.; cp[i].coortmp[2] = 0.;
		if (i >= npt_old)
		{
			cp[i].coor[0] = 0.; cp[i].coor[1] = 0.; cp[i].coor[2] = 0.;
		}
	}
	vector<int> rfid;
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 0 && tmesh[i].ref_flag != 0)
		{
			rfid.push_back(i);
			CalPatchCP_Regular(i);
		}
	}

	UpdateGeometry(rfid);
}

void TruncatedTspline_3D::CalPatchCP_Regular(int eid)
{
	vector<int> IEN_new_all;
	vector<array<double, 5>> KU_new_all;
	vector<array<double, 5>> KV_new_all;
	vector<array<double, 5>> KW_new_all;
	for (uint i = 0; i<tmesh[eid].chd.size(); i++)
	{
		int cid(tmesh[eid].chd[i]);
		for (uint j = 0; j<tmesh[cid].IENtmp.size(); j++)
		{
			vector<int>::iterator it = find(IEN_new_all.begin(), IEN_new_all.end(), tmesh[cid].IENtmp[j]);
			if (it == IEN_new_all.end())
			{
				IEN_new_all.push_back(tmesh[cid].IENtmp[j]);
				array<double, 5> kutmp, kvtmp, kwtmp;
				for (int k = 0; k<5; k++)
				{
					kutmp[k] = tmesh[cid].patch_kutmp[j][k] + tmesh[eid].chd_o[i][0];
					kvtmp[k] = tmesh[cid].patch_kvtmp[j][k] + tmesh[eid].chd_o[i][1];
					kwtmp[k] = tmesh[cid].patch_kwtmp[j][k] + tmesh[eid].chd_o[i][2];
				}
				KU_new_all.push_back(kutmp);
				KV_new_all.push_back(kvtmp);
				KW_new_all.push_back(kwtmp);
			}
		}
	}

	vector<int> IEN_old(tmesh[eid].IEN);
	vector<vector<double>> cmat(IEN_new_all.size(), vector<double>(IEN_old.size(), 0.));
	for (uint i = 0; i<IEN_old.size(); i++)
	{
		vector<int>::iterator it = find(IEN_new_all.begin(), IEN_new_all.end(), IEN_old[i]);
		if (it != IEN_new_all.end())
		{
			int loc(it - IEN_new_all.begin());
			if (tmesh[eid].patch_ku[i] != KU_new_all[loc] || tmesh[eid].patch_kv[i] != KV_new_all[loc] || tmesh[eid].patch_kw[i] != KW_new_all[loc])
			{
				cp[IEN_old[i]].aff = 1;
			}
			else
			{
				cmat[loc][i] = 1.;
			}
		}
	}

	for (uint i = 0; i<IEN_new_all.size(); i++)
	{
		for (uint j = 0; j<IEN_old.size(); j++)
		{
			if (CheckSubKnotVector(KU_new_all[i], KV_new_all[i], KW_new_all[i], tmesh[eid].patch_ku[j], tmesh[eid].patch_kv[j], tmesh[eid].patch_kw[j]))
			{
				vector<double> ku1(10), kv1(10), kw1(10);
				vector<vector<double>> Tu, Tv, Tw;
				vector<double> ku(tmesh[eid].patch_ku[j].begin(), tmesh[eid].patch_ku[j].end());
				vector<double> kv(tmesh[eid].patch_kv[j].begin(), tmesh[eid].patch_kv[j].end());
				vector<double> kw(tmesh[eid].patch_kw[j].begin(), tmesh[eid].patch_kw[j].end());
				vector<double>::iterator it1, it2, it3;
				it1 = set_union(tmesh[eid].patch_ku[j].begin(), tmesh[eid].patch_ku[j].end(), KU_new_all[i].begin(), KU_new_all[i].end(), ku1.begin());
				it2 = set_union(tmesh[eid].patch_kv[j].begin(), tmesh[eid].patch_kv[j].end(), KV_new_all[i].begin(), KV_new_all[i].end(), kv1.begin());
				it3 = set_union(tmesh[eid].patch_kw[j].begin(), tmesh[eid].patch_kw[j].end(), KW_new_all[i].begin(), KW_new_all[i].end(), kw1.begin());
				ku1.resize(it1 - ku1.begin());
				kv1.resize(it2 - kv1.begin());
				kw1.resize(it3 - kw1.begin());
				TMatrix(ku, ku1, 3, Tu);
				TMatrix(kv, kv1, 3, Tv);
				TMatrix(kw, kw1, 3, Tw);
				it1 = search(ku1.begin(), ku1.end(), KU_new_all[i].begin(), KU_new_all[i].end());
				it2 = search(kv1.begin(), kv1.end(), KV_new_all[i].begin(), KV_new_all[i].end());
				it3 = search(kw1.begin(), kw1.end(), KW_new_all[i].begin(), KW_new_all[i].end());
				if (it1 != ku1.end() && it2 != kv1.end() && it3 != kw1.end())
				{
					int loc1 = it1 - ku1.begin();
					int loc2 = it2 - kv1.begin();
					int loc3 = it3 - kw1.begin();
					cmat[i][j] = Tu[loc1][0] * Tv[loc2][0] * Tw[loc3][0];
				}
			}
		}
	}
	for (uint i = 0; i<IEN_new_all.size(); i++)
	{
		if (cp[IEN_new_all[i]].update == 0)
		{
			for (uint j = 0; j<IEN_old.size(); j++)
			{
				if (cmat[i][j] != 0.)
				{
					cp[IEN_new_all[i]].coortmp[0] += cmat[i][j] * cp[IEN_old[j]].coor[0];
					cp[IEN_new_all[i]].coortmp[1] += cmat[i][j] * cp[IEN_old[j]].coor[1];
					cp[IEN_new_all[i]].coortmp[2] += cmat[i][j] * cp[IEN_old[j]].coor[2];
				}
			}
			cp[IEN_new_all[i]].update = 1;
		}
	}

	//check truncation
	double uvw[3][3] = { { 0., tmedge[tmesh[eid].edge[0]].len / 2., tmedge[tmesh[eid].edge[0]].len }, { 0., tmedge[tmesh[eid].edge[3]].len / 2., tmedge[tmesh[eid].edge[3]].len },
	{ 0., tmedge[tmesh[eid].edge[4]].len / 2., tmedge[tmesh[eid].edge[4]].len } };
	vector<array<double, 3>> spt(27);
	int loc(0);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				spt[loc][0] = uvw[0][k];
				spt[loc][1] = uvw[1][j];
				spt[loc][2] = uvw[2][i];
				loc++;
			}
		}
	}
	for (uint i = 0; i<IEN_old.size(); i++)
	{
		vector<double> coef;
		for (uint j = 0; j<IEN_new_all.size(); j++)
		{
			coef.push_back(cmat[j][i]);
		}
		if (!CheckFullRefine(spt, tmesh[eid].patch_ku[i], tmesh[eid].patch_kv[i], tmesh[eid].patch_kw[i], KU_new_all, KV_new_all, KW_new_all, coef))
		{
			cp[IEN_old[i]].truntmp = 1;
			for (uint j = 0; j<IEN_new_all.size(); j++)
			{
				if (cmat[j][i] != 0. && IEN_old[i] != IEN_new_all[j])
				{
					vector<int>::iterator it = find(cp[IEN_old[i]].tbftmp.begin(), cp[IEN_old[i]].tbftmp.end(), IEN_new_all[j]);
					if (it == cp[IEN_old[i]].tbftmp.end())
					{
						cp[IEN_old[i]].tbftmp.push_back(IEN_new_all[j]);
						cp[IEN_old[i]].tctmp.push_back(cmat[j][i]);
					}
				}
			}
		}
	}
}

void TruncatedTspline_3D::UpdateGeometry(vector<int>& rfid)
{
	for (uint i = 0; i<cp.size(); i++)
	{
		if (i<npt_old)
		{
			if (cp[i].truntmp == 0 && cp[i].aff == 1)
			{
				//cout << "update coor for old points...\n";
				//getchar();
				cp[i].coor[0] = cp[i].coortmp[0];
				cp[i].coor[1] = cp[i].coortmp[1];
				cp[i].coor[2] = cp[i].coortmp[2];
				//cp[i].w = cp[i].wtmp;
			}
			else if (cp[i].truntmp == 1)
			{
				//cout << "update trun points...\n";
				//getchar();
				cp[i].trun = 1;
				cp[i].tbf = cp[i].tbftmp;
				cp[i].tc = cp[i].tctmp;
				vector<int>().swap(cp[i].tbftmp);
				vector<double>().swap(cp[i].tctmp);
			}
		}
		else
		{
			cp[i].coor[0] = cp[i].coortmp[0];
			cp[i].coor[1] = cp[i].coortmp[1];
			cp[i].coor[2] = cp[i].coortmp[2];
			//cp[i].w = cp[i].wtmp;
		}
	}

	//update new elements
	for (uint i = 0; i<rfid.size(); i++)
	{
		if (tmesh[rfid[i]].type == 0 /*|| tmesh[rfid[i]].type == 1 || tmesh[rfid[i]].type == 2*/)
		{
			vector<int> trunIEN;
			vector<array<double, 5>> trun_patch_ku;
			vector<array<double, 5>> trun_patch_kv;
			vector<array<double, 5>> trun_patch_kw;
			for (uint j = 0; j<tmesh[rfid[i]].IEN.size(); j++)
			{
				if (cp[tmesh[rfid[i]].IEN[j]].trun == 1)
				{
					trunIEN.push_back(tmesh[rfid[i]].IEN[j]);
					trun_patch_ku.push_back(tmesh[rfid[i]].patch_ku[j]);
					trun_patch_kv.push_back(tmesh[rfid[i]].patch_kv[j]);
					trun_patch_kw.push_back(tmesh[rfid[i]].patch_kw[j]);
				}
			}
			for (uint j = 0; j<tmesh[rfid[i]].chd.size(); j++)
			{
				int chdid(tmesh[rfid[i]].chd[j]);
				tmesh[chdid].IEN.clear();
				tmesh[chdid].patch_ku.clear();
				tmesh[chdid].patch_kv.clear();
				tmesh[chdid].patch_kw.clear();
				for (uint k = 0; k<tmesh[chdid].IENtmp.size(); k++)
				{
					tmesh[chdid].IEN.push_back(tmesh[chdid].IENtmp[k]);
					if (cp[tmesh[chdid].IENtmp[k]].trun == 0)
					{
						tmesh[chdid].patch_ku.push_back(tmesh[chdid].patch_kutmp[k]);
						tmesh[chdid].patch_kv.push_back(tmesh[chdid].patch_kvtmp[k]);
						tmesh[chdid].patch_kw.push_back(tmesh[chdid].patch_kwtmp[k]);
					}
					else
					{
						vector<int>::iterator it = find(tmesh[rfid[i]].IEN.begin(), tmesh[rfid[i]].IEN.end(), tmesh[chdid].IENtmp[k]);
						if (it != tmesh[rfid[i]].IEN.end())
						{
							int loc(it - tmesh[rfid[i]].IEN.begin());
							array<double, 5> kutmp, kvtmp, kwtmp;
							for (int a = 0; a<5; a++)
							{
								kutmp[a] = tmesh[rfid[i]].patch_ku[loc][a] - tmesh[rfid[i]].chd_o[j][0];
								kvtmp[a] = tmesh[rfid[i]].patch_kv[loc][a] - tmesh[rfid[i]].chd_o[j][1];
								kwtmp[a] = tmesh[rfid[i]].patch_kw[loc][a] - tmesh[rfid[i]].chd_o[j][2];
							}
							tmesh[chdid].patch_ku.push_back(kutmp);
							tmesh[chdid].patch_kv.push_back(kvtmp);
							tmesh[chdid].patch_kw.push_back(kwtmp);
						}
						else
						{
							//cout<<i<<" "<<rfid[i]<<" "<<tmesh[rfid[i]].IEN.size()<<"\n";
							cout << "Cannot find truncated ID in the old set!\n";
							getchar();
						}
					}
				}
				for (uint k = 0; k<trunIEN.size(); k++)
				{
					vector<int>::iterator it = find(tmesh[chdid].IEN.begin(), tmesh[chdid].IEN.end(), trunIEN[k]);
					if (it == tmesh[chdid].IEN.end())
					{
						tmesh[chdid].IEN.push_back(trunIEN[k]);
						array<double, 5> kutmp, kvtmp, kwtmp;
						for (int a = 0; a<5; a++)
						{
							kutmp[a] = trun_patch_ku[k][a] - tmesh[rfid[i]].chd_o[j][0];
							kvtmp[a] = trun_patch_kv[k][a] - tmesh[rfid[i]].chd_o[j][1];
							kwtmp[a] = trun_patch_kw[k][a] - tmesh[rfid[i]].chd_o[j][2];
						}
						tmesh[chdid].patch_ku.push_back(kutmp);
						tmesh[chdid].patch_kv.push_back(kvtmp);
						tmesh[chdid].patch_kw.push_back(kwtmp);
					}
				}
			}
		}
		//else if (tmesh[rfid[i]].type == 4 || tmesh[rfid[i]].type == 5)
		//{
		//	for (int j = 0; j<4; j++)
		//	{
		//		int chdid(tmesh[rfid[i]].chd[j]);
		//		if (chdid != -1)
		//		{
		//			tmesh[chdid].IEN.clear();
		//			tmesh[chdid].patch_ku.clear();
		//			tmesh[chdid].patch_kv.clear();
		//			tmesh[chdid].IEN = tmesh[chdid].IENtmp;
		//			tmesh[chdid].patch_ku = tmesh[chdid].patch_kutmp;
		//			tmesh[chdid].patch_kv = tmesh[chdid].patch_kvtmp;
		//			vector<int>().swap(tmesh[chdid].IENtmp);
		//			vector<array<double, 5>>().swap(tmesh[chdid].patch_kutmp);
		//			vector<array<double, 5>>().swap(tmesh[chdid].patch_kvtmp);
		//		}
		//	}
		//}
	}

	//update old elements
	for (uint i = 0; i<nel_old; i++)
	{
		if (tmesh[i].act == 1 && (tmesh[i].type == 0 /*|| tmesh[i].type == 1 || tmesh[i].type == 2*/))
		{
			vector<int> ien_old(tmesh[i].IEN);
			vector<array<double, 5>> ku_old(tmesh[i].patch_ku);
			vector<array<double, 5>> kv_old(tmesh[i].patch_kv);
			vector<array<double, 5>> kw_old(tmesh[i].patch_kw);
			tmesh[i].IEN.clear();
			tmesh[i].patch_ku.clear();
			tmesh[i].patch_kv.clear();
			tmesh[i].patch_kw.clear();
			//
			for (uint j = 0; j<ien_old.size(); j++)
			{
				if (cp[ien_old[j]].trun == 1)
				{
					tmesh[i].IEN.push_back(ien_old[j]);
					tmesh[i].patch_ku.push_back(ku_old[j]);
					tmesh[i].patch_kv.push_back(kv_old[j]);
					tmesh[i].patch_kw.push_back(kw_old[j]);
				}
			}
			//
			for (uint j = 0; j<tmesh[i].IENtmp.size(); j++)
			{
				//tmesh[i].IEN.push_back(tmesh[i].IENtmp[j]);
				if (cp[tmesh[i].IENtmp[j]].trun == 0)
				{
					tmesh[i].IEN.push_back(tmesh[i].IENtmp[j]);
					tmesh[i].patch_ku.push_back(tmesh[i].patch_kutmp[j]);
					tmesh[i].patch_kv.push_back(tmesh[i].patch_kvtmp[j]);
					tmesh[i].patch_kw.push_back(tmesh[i].patch_kwtmp[j]);
				}
				//else
				//{
				//	vector<int>::iterator it=find(ien_old.begin(),ien_old.end(),tmesh[i].IENtmp[j]);
				//	int loc(it-ien_old.begin());
				//	tmesh[i].patch_ku.push_back(ku_old[loc]);
				//	tmesh[i].patch_kv.push_back(kv_old[loc]);
				//}
			}
			vector<int>().swap(tmesh[i].IENtmp);
			vector<array<double, 5>>().swap(tmesh[i].patch_kutmp);
			vector<array<double, 5>>().swap(tmesh[i].patch_kvtmp);
			vector<array<double, 5>>().swap(tmesh[i].patch_kwtmp);
		}
		/*else if (tmesh[i].act == 1 && tmesh[i].type == 4)
		{
			vector<int> ien_old(tmesh[i].IEN);
			vector<array<double, 5>> ku_old(tmesh[i].patch_ku);
			vector<array<double, 5>> kv_old(tmesh[i].patch_kv);
			tmesh[i].IEN.clear();
			tmesh[i].patch_ku.clear();
			tmesh[i].patch_kv.clear();
			uint n_1r(2 * cp[tmesh[i].cnct[0]].face.size() + 1);
			for (uint j = 0; j<tmesh[i].IENtmp.size(); j++)
			{
				tmesh[i].IEN.push_back(tmesh[i].IENtmp[j]);
				if (j >= n_1r)
				{
					if (cp[tmesh[i].IENtmp[j]].trun == 0)
					{
						tmesh[i].patch_ku.push_back(tmesh[i].patch_kutmp[j - n_1r]);
						tmesh[i].patch_kv.push_back(tmesh[i].patch_kvtmp[j - n_1r]);
					}
					else
					{
						vector<int>::iterator it = find(ien_old.begin(), ien_old.end(), tmesh[i].IENtmp[j]);
						int loc(it - ien_old.begin());
						tmesh[i].patch_ku.push_back(ku_old[loc - n_1r]);
						tmesh[i].patch_kv.push_back(kv_old[loc - n_1r]);
					}
				}
			}
			vector<int>().swap(tmesh[i].IENtmp);
			vector<array<double, 5>>().swap(tmesh[i].patch_kutmp);
			vector<array<double, 5>>().swap(tmesh[i].patch_kvtmp);
		}*/
	}
}

bool TruncatedTspline_3D::CheckFullRefine(const vector<array<double, 3>>& spt, const array<double, 5>& motu, const array<double, 5>& motv, const array<double, 5>& motw,
	const vector<array<double, 5>>& chdu, const vector<array<double, 5>>& chdv, const vector<array<double, 5>>& chdw, const vector<double>& coef)
{
	double tol(1.e-8);
	double sum(0.);
	for (uint i = 0; i<spt.size(); i++)
	{
		vector<double> ku(motu.begin(), motu.end());
		vector<double> kv(motv.begin(), motv.end());
		vector<double> kw(motw.begin(), motw.end());
		vector<double> uval, vval, wval;
		BSplineBasis bu, bv, bw;
		bu.Set(3, ku);
		bv.Set(3, kv);
		bw.Set(3, kw);
		bu.BasisFunction(0, spt[i][0], 0, uval);
		bv.BasisFunction(0, spt[i][1], 0, vval);
		bw.BasisFunction(0, spt[i][2], 0, wval);
		double Ni = uval[0] * vval[0] * wval[0];
		for (uint j = 0; j<chdu.size(); j++)
		{
			if (coef[j] != 0.)
			{
				ku.assign(chdu[j].begin(), chdu[j].end());
				kv.assign(chdv[j].begin(), chdv[j].end());
				kw.assign(chdw[j].begin(), chdw[j].end());
				bu.Set(3, ku);
				bv.Set(3, kv);
				bw.Set(3, kw);
				bu.BasisFunction(0, spt[i][0], 0, uval);
				bv.BasisFunction(0, spt[i][1], 0, vval);
				bw.BasisFunction(0, spt[i][2], 0, wval);
				Ni -= coef[j] * uval[0] * vval[0] * wval[0];
			}
		}
		if (Ni<0.) Ni = -Ni;
		sum += Ni;
	}
	if (sum / spt.size() < tol)
		return true;
	else
		return false;
}

bool TruncatedTspline_3D::CheckSubKnotVector(const array<double, 5>& ku1, const array<double, 5>& kv1, const array<double, 5>& kw1, 
	const array<double, 5>& ku2, const array<double, 5>& kv2, const array<double, 5>& kw2)// 1 is a child of 2
{
	if (ku1[0] >= ku2[0] && ku1[4] <= ku2[4] && kv1[0] >= kv2[0] && kv1[4] <= kv2[4] && kw1[0] >= kw2[0] && kw1[4] <= kw2[4])
	{
		double min_u1(10.), min_v1(10.), min_w1(10.), min_u2(10.), min_v2(10.), min_w2(10.);
		for (int i = 0; i<4; i++)
		{
			double tmp = ku1[i + 1] - ku1[i];
			if (tmp != 0. && tmp<min_u1) min_u1 = tmp;
			tmp = kv1[i + 1] - kv1[i];
			if (tmp != 0. && tmp<min_v1) min_v1 = tmp;
			tmp = kw1[i + 1] - kw1[i];
			if (tmp != 0. && tmp<min_w1) min_w1 = tmp;
			tmp = ku2[i + 1] - ku2[i];
			if (tmp != 0. && tmp<min_u2) min_u2 = tmp;
			tmp = kv2[i + 1] - kv2[i];
			if (tmp != 0. && tmp<min_v2) min_v2 = tmp;
			tmp = kw2[i + 1] - kw2[i];
			if (tmp != 0. && tmp<min_w2) min_w2 = tmp;
		}
		if (min_u1 <= min_u2 && min_v1 <= min_v2 && min_w1 <= min_w2)
		{
			for (int i = 0; i<5; i++)
			{
				array<double, 5>::const_iterator it1 = find(ku2.begin(), ku2.end(), ku1[i]);
				if (it1 == ku2.end())
				{
					return true;
				}
				array<double, 5>::const_iterator it2 = find(kv2.begin(), kv2.end(), kv1[i]);
				if (it2 == kv2.end())
				{
					return true;
				}
				array<double, 5>::const_iterator it3 = find(kw2.begin(), kw2.end(), kw1[i]);
				if (it2 == kw2.end())
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

void TruncatedTspline_3D::run(string fn)
{
	SetProblem(fn);
	
	//RefineTest_2();
	RefineTopology();
	//RefineGeometry();
	RefineGeometry_v1();

	//AllBezierPatch();
	//CollectActives();

	string pfn("../io/hex_out/cube9"),testnum("1");
	//string pfn("../io/hex_out/cube_dense"), testnum("5");

	//VisualizeControlMesh(pfn+"_CM_"+testnum);
	//VisualizeTMesh(pfn + "_edge_" + testnum);
	//VisualizeFaceMesh(pfn+"_face_"+testnum);s
	VisualizeSolidVTK(pfn + "_geom_" + testnum);

	//CreateUniformCube("../io/hex_input/cube5");
	//VisualizeSolidVTK("../io/hex_out/cube_dense_1");
}

void TruncatedTspline_3D::runh(string fn)
{
	SetProblem(fn);

	vector<array<int, 2>> rfid;
	vector<array<int, 2>> gst;

	//Identify_Pseudo(rfid);//for cube9
	//Refine(rfid, gst);
	//Identify_Pseudo_1(rfid);
	//Refine(rfid, gst);

	//Identify_Test_2(rfid, gst);//for cube_dense
	//Refine(rfid, gst);
	//Identify_Test_3(rfid, gst);
	//Refine(rfid, gst);

	//Identify_Test_4(rfid, gst);//for cube_dense, boundary
	//Refine(rfid, gst);

	string pfn("../io/hex_out_1/cubeline0"), testnum("1"), levnum("0");
	//string pfn("../io/hex_out_1/cube_dense"), testnum("1"), levnum("0");
	//string pfn("../io/hex_out_1/cube_coarse"), testnum("1"), levnum("0");
	int lev(0);

	//PillowCube(fn+"_pillow");

	//OutputBasis(0,21,"../io/hex_out_1/cube3_basis_21");

	//OutputCM(lev, pfn + "_CM_" + testnum+"_"+levnum);
	OutputFace(lev, pfn + "_FC_" + testnum+"_"+levnum);
	OutputEdge(lev, pfn + "_ED_" + testnum+"_"+levnum);
	//OutputGeom(lev, pfn + "_geom_" + testnum+"_"+levnum);

	//VisualizeControlMesh(pfn+"_CM_"+testnum);
	//VisualizeTMesh(pfn + "_edge_" + testnum);
	//VisualizeFaceMesh(pfn+"_face_"+testnum);
	//VisualizeSolidVTK(pfn + "_geom_" + testnum);

	//CreateUniformCube("../io/hex_input/cube5");
}

//void TruncatedTspline_3D::VisualizeSurface(string fn)
//{
//	vector<array<double,3>> spt;
//	vector<array<double,3>> sval;
//	vector<array<int,4>> sele;
//	vector<array<double,3>> lpt;//visulize parameter lines
//	vector<array<int,2>> led;//line connectivity
//	int ns(5),ecount(0),loc0,loc1,loc2;
//	//vector<double> su(ns),sv(ns);
//	//for(int i=0; i<ns; i++)
//	//{
//	//	su[i]=i*1./(ns-1);
//	//	sv[i]=i*1./(ns-1);
//	//}
//
//	for(uint e=0;e<tmesh.size();e++)
//	{
//		if(tmesh[e].act==1 && (tmesh[e].type==0 || tmesh[e].type==1 || tmesh[e].type==4))
//		{
//			int loc(0);
//			vector<double> su(ns),sv(ns);
//			for(int i=0; i<ns; i++)
//			{
//				su[i]=i*tmedge[tmesh[e].edge[0]].len/(ns-1);
//				sv[i]=i*tmedge[tmesh[e].edge[3]].len/(ns-1);
//			}
//
//			for(int a=0;a<ns;a++)
//			{
//				for(int b=0;b<ns;b++)
//				{
//					array<double,3> pt;
//					array<double,3> nm;
//					SurfacePointMap(e,su[b],sv[a],pt,nm);
//					spt.push_back(pt);
//					sval.push_back(nm);
//					if(a==0||a==ns-1||b==0||b==ns-1)
//					{
//						lpt.push_back(pt);
//					}
//				}
//			}
//
//			for(int a=0;a<ns-1;a++)
//			{
//				for(int b=0;b<ns-1;b++)
//				{
//					array<int,4> el;
//					el[0]=ecount*ns*ns+a*ns+b;
//					el[1]=ecount*ns*ns+a*ns+b+1;
//					el[2]=ecount*ns*ns+(a+1)*ns+b+1;
//					el[3]=ecount*ns*ns+(a+1)*ns+b;
//					sele.push_back(el);
//				}
//			}
//			for(int a=0;a<ns-1;a++)
//			{
//				array<int,2> lc;
//				lc[0]=ecount*4*(ns-1)+a;
//				lc[1]=ecount*4*(ns-1)+a+1;
//				led.push_back(lc);
//				lc[0]=ecount*4*(ns-1)+3*ns-4+a;
//				lc[1]=ecount*4*(ns-1)+3*ns-4+a+1;
//				led.push_back(lc);
//			}
//			for(int a=0;a<ns-2;a++)
//			{
//				array<int,2> lc;
//				lc[0]=ecount*4*(ns-1)+ns+2*a;
//				lc[1]=ecount*4*(ns-1)+ns+2*a+2;
//				led.push_back(lc);
//				lc[0]=ecount*4*(ns-1)+ns+2*a-1;
//				lc[1]=ecount*4*(ns-1)+ns+2*a+1;
//				led.push_back(lc);
//			}
//			array<int,2> lc1;
//			lc1[0]=ecount*4*(ns-1);
//			lc1[1]=ecount*4*(ns-1)+ns;
//			led.push_back(lc1);
//			lc1[0]=ecount*4*(ns-1)+3*ns-5;
//			lc1[1]=ecount*4*(ns-1)+4*ns-5;
//			led.push_back(lc1);
//			ecount++;
//		}
//	}
//
//	string fname=fn+".vtk";
//	ofstream fout;
//	fout.open(fname.c_str());
//	if(fout.is_open())
//	{
//		fout<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
//		fout<<"POINTS "<<spt.size()<<" float\n";
//		for(uint i=0;i<spt.size();i++)
//		{
//			fout<<spt[i][0]<<" "<<spt[i][1]<<" "<<spt[i][2]<<"\n";
//		}
//		fout<<"\nCELLS "<<sele.size()<<" "<<5*sele.size()<<'\n';
//		for(uint i=0;i<sele.size();i++)
//		{
//			fout<<"4 "<<sele[i][0]<<" "<<sele[i][1]<<" "<<sele[i][2]<<" "<<sele[i][3]<<'\n';
//		}
//		fout<<"\nCELL_TYPES "<<sele.size()<<'\n';
//		for(uint i=0;i<sele.size();i++)
//		{
//			fout<<"9\n";
//		}
//		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
//		//for(uint i=0;i<sval.size();i++)
//		//{
//		//	fout<<sval[i]<<"\n";
//		//}
//		fout<<"\nPOINT_DATA "<<sval.size()<<"\nNORMALS Normal FLOAT\n";
//		for(uint i=0;i<sval.size();i++)
//		{
//			fout<<sval[i][0]<<" "<<sval[i][1]<<" "<<sval[i][2]<<"\n";
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

void TruncatedTspline_3D::SetProblem(string fn)
{
	//read hex vtk
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
			cp[i].act = 1;
			fin>>cp[i].coor[0]>>cp[i].coor[1]>>cp[i].coor[2];
		}
		getline(fin,stmp);
		fin>>stmp>>neles>>itmp;
		tmesh.resize(neles);
		for(int i=0;i<neles;i++)
		{
			tmesh[i].act = 1;
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3] >> 
				tmesh[i].cnct[4] >> tmesh[i].cnct[5] >> tmesh[i].cnct[6] >> tmesh[i].cnct[7];
		}
		fin.close();
	}
	else
	{
		cerr<<"Cannot open "<<fname<<"!\n";
	}

	//SetDomain();//normalize coordinates in [0,1]^3, for the cube poisson

	InitialConnect();
	//InitialRotate();
	//FindEdgeTopoDirec();
	//FindKnotInterval();
	//SetLocalCoorSystem();
	//FindIEN_PatchKV();
	//Update_IEN();

	hmesh.push_back(tmesh);
	hcp.push_back(cp);
	hface.push_back(tmface);
	hedge.push_back(tmedge);

	ReportXP();

	//construct boundary here
	SetSharpFeature();

	//string fn1("../io/pipeline/neuron1");
	//OutputEdge(0, fn1);
	////OutputCM(0, fn1);
	//cout << "done setting sharp feature\n";
	//getchar();

	AllBezierLev(0);
	//for (uint i = 0; i < hmesh[0].size(); i++)
	//{
	//	if (hmesh[0][i].type == 1)
	//	{
	//		ConstructBezierBasis_Boundary(0,i);
	//	}
	//}

	BezierPoints_ini();//represent geometry using initial control points and blending functions

	SetSupport(0);
}

void TruncatedTspline_3D::SetSharpFeature()
{
	//sharp edge and sharp corner, indicated by control point
	//sharp edge
	//double tol(.2);
	//double tol(.3);
	double tol(.57);//rod, base
	//double tol(.8);//hook
	int fc_dir[6][3] = { { 0, 3, 2 }, { 0, 1, 5 }, { 1, 2, 6 }, { 3, 7, 6 }, { 0, 4, 7 }, { 4, 5, 6 } };
	uint i, j, k;
	for (i = 0; i < hedge[0].size(); i++)
	{
		if (hedge[0][i].type == 1)
		{
			vector<array<double, 3>> vec;
			for (j = 0; j < hedge[0][i].hex.size(); j++)
			{
				if (hmesh[0][hedge[0][i].hex[j]].type == 1)
				{
					int hxid(hedge[0][i].hex[j]);
					for (k = 0; k < 6; k++)
					{
						if (hface[0][hmesh[0][hxid].face[k]].type == 1)
						{
							vector<int>::iterator it = find(hedge[0][i].face.begin(), hedge[0][i].face.end(), hmesh[0][hxid].face[k]);
							if (it != hedge[0][i].face.end())
							{
								array<double, 3> tmp1, tmp2;
								for (int dof = 0; dof < 3; dof++)
								{
									tmp1[dof] = hcp[0][hmesh[0][hxid].cnct[fc_dir[k][1]]].coor[dof] - hcp[0][hmesh[0][hxid].cnct[fc_dir[k][0]]].coor[dof];
									tmp2[dof] = hcp[0][hmesh[0][hxid].cnct[fc_dir[k][2]]].coor[dof] - hcp[0][hmesh[0][hxid].cnct[fc_dir[k][1]]].coor[dof];
								}
								array<double, 3> tmp3 = { tmp1[1] * tmp2[2] - tmp1[2] * tmp2[1], -tmp1[0] * tmp2[2] + tmp1[2] * tmp2[0], tmp1[0] * tmp2[1] - tmp1[1] * tmp2[0] };
								double dis = sqrt(tmp3[0] * tmp3[0] + tmp3[1] * tmp3[1] + tmp3[2] * tmp3[2]);
								tmp3[0] /= dis; tmp3[1] /= dis; tmp3[2] /= dis;
								vec.push_back(tmp3);
							}
						}
					}
				}
			}
			if (vec.size() == 2)
			{
				double ang(vec[0][0] * vec[1][0] + vec[0][1] * vec[1][1] + vec[0][2] * vec[1][2]);
				if (ang < tol)
				{
					hedge[0][i].sharp = 1;
					hcp[0][hedge[0][i].pt[0]].sharp = 1;
					hcp[0][hedge[0][i].pt[1]].sharp = 1;
				}
			}
			else
			{
				cerr << "Something wrong in determining sharp edge!\n";
				getchar();
			}
		}
	}

	for (i = 0; i < hcp[0].size(); i++)
	{
		if (hcp[0][i].sharp == 1)
		{
			int nshp(0);
			for (j = 0; j < hcp[0][i].edge.size(); j++)
			{
				if (hedge[0][hcp[0][i].edge[j]].sharp == 1)
				{
					nshp++;
				}
			}
			if (nshp >= 3) hcp[0][i].sharp = 2;//sharp corner
			//additional
			if (nshp == 1) hcp[0][i].sharp = 0;
		}
	}

	//additional
	for (i = 0; i < hedge[0].size(); i++)
	{
		if (hedge[0][i].sharp == 1)
		{
			if (hcp[0][hedge[0][i].pt[0]].sharp == 0 || hcp[0][hedge[0][i].pt[1]].sharp == 0)
			{
				hedge[0][i].sharp = 0;
			}
		}
	}
	//for (i = 0; i < hcp[0].size(); i++)
	//{
	//	if ((hcp[0][i].bcxp == 1 && hcp[0][i].sharp == 1) || hcp[0][i].bcxp==2)
	//	{
	//		hcp[0][i].sharp = 2;
	//	}
	//}

	//int nshp_ed(9), nshp_cn(8);
	//int shp_ed[] = { 2792, 2808, 2824, 2840, 2856, 2872, 2888, 2904, 2920};
	//int shp_cn[] = {2600,2872,2436,2276,2596,2856,2760,2468};
	//for (int i = 0; i < nshp_ed; i++)
	//{
	//	if (hcp[0][shp_ed[i]].type == 1)
	//	{
	//		hcp[0][shp_ed[i]].sharp = 1;
	//	}
	//	else
	//	{
	//		cerr << "Sharp edge not on the boundary!\n";
	//		getchar();
	//	}
	//}
	//for (int i = 0; i < nshp_cn; i++)
	//{
	//	if (hcp[0][shp_cn[i]].type == 1)
	//	{
	//		hcp[0][shp_cn[i]].sharp = 2;
	//	}
	//	else
	//	{
	//		cerr << "Sharp corner not on the boundary!\n";
	//		getchar();
	//	}
	//}
	//for (uint i = 0; i < hedge[0].size(); i++)
	//{
	//	if (hedge[0][i].type == 1 && hcp[0][hedge[0][i].pt[0]].sharp != 0 && hcp[0][hedge[0][i].pt[1]].sharp != 0)
	//	{
	//		hedge[0][i].sharp = 1;
	//	}
	//}
}

void TruncatedTspline_3D::SetDomain()
{
	double tmp(1.e6);
	double x_range[3][2] = { { tmp, -tmp }, { tmp, -tmp }, { tmp, -tmp } };
	uint i, j;
	for (i = 0; i <cp.size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (cp[i].coor[j] < x_range[j][0]) x_range[j][0] = cp[i].coor[j];
			if (cp[i].coor[j] > x_range[j][1]) x_range[j][1] = cp[i].coor[j];
		}
	}
	double xh[3] = { x_range[0][1] - x_range[0][0], x_range[1][1] - x_range[1][0], x_range[2][1] - x_range[2][0] };
	for (i = 0; i <cp.size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			cp[i].coor[j] = (cp[i].coor[j] - x_range[j][0]) / xh[j];
		}
	}
}

void TruncatedTspline_3D::BezierPoints_ini()
{
	for (uint j = 0; j < hmesh[0].size(); j++)
	{
		for (uint k = 0; k < hmesh[0][j].IEN.size(); k++)
		{
			for (int k1 = 0; k1 < 64; k1++)
			{
				hmesh[0][j].bzpt[k1][0] += hmesh[0][j].bemat[k][k1] * hcp[0][hmesh[0][j].IEN[k]].coor[0];
				hmesh[0][j].bzpt[k1][1] += hmesh[0][j].bemat[k][k1] * hcp[0][hmesh[0][j].IEN[k]].coor[1];
				hmesh[0][j].bzpt[k1][2] += hmesh[0][j].bemat[k][k1] * hcp[0][hmesh[0][j].IEN[k]].coor[2];
			}
		}
	}
}

void TruncatedTspline_3D::BezierPoints_Refine(int lev, int pos, const vector<MatrixXd>& bsmat)
{
	for (uint i = 0; i < hmesh[lev][pos].chd.size(); i++)
	{
		int cid[2] = { lev + 1, hmesh[lev][pos].chd[i] };
		for (int j = 0; j < 64; j++)//chd
		{
			hmesh[cid[0]][cid[1]].bzpt[j][0] = 0.;
			hmesh[cid[0]][cid[1]].bzpt[j][1] = 0.;
			hmesh[cid[0]][cid[1]].bzpt[j][2] = 0.;
			for (int k = 0; k < 64; k++)//father
			{
				hmesh[cid[0]][cid[1]].bzpt[j][0] += bsmat[i](j, k)*hmesh[lev][pos].bzpt[k][0];
				hmesh[cid[0]][cid[1]].bzpt[j][1] += bsmat[i](j, k)*hmesh[lev][pos].bzpt[k][1];
				hmesh[cid[0]][cid[1]].bzpt[j][2] += bsmat[i](j, k)*hmesh[lev][pos].bzpt[k][2];
			}
		}
	}
}

void TruncatedTspline_3D::BezierSubdivMatrix(vector<MatrixXd>& bsmat)
{
	bsmat.resize(8);
	for (uint i = 0; i < bsmat.size(); i++)
	{
		bsmat[i] = MatrixXd::Zero(64, 64);
	}
	double bs1d[8][4] = { { 1., 0., 0., 0. }, { .5, .5, 0., 0. }, { .25, .5, .25, 0. }, { .125, .375, .375, .125 }, 
	{ .125, .375, .375, .125 }, { 0., .25, .5, .25 }, { 0., 0., .5, .5 }, { 0., 0., 0., 1. } };
	int i0, j0, k0, i1, j1, k1, loc0(0);
	//vector<int> loc1(8, 0);
	double val;
	for (k0 = 0; k0 < 4; k0++)
	{
		for (j0 = 0; j0 < 4; j0++)
		{
			for (i0 = 0; i0 < 4; i0++)
			{
				vector<int> loc1(8, 0);
				for (k1 = 0; k1 < 8; k1++)
				{
					for (j1 = 0; j1 < 8; j1++)
					{
						for (i1 = 0; i1 < 8; i1++)
						{
							val = bs1d[i1][i0] * bs1d[j1][j0] * bs1d[k1][k0];
							if (i1 < 4 && j1 < 4 && k1 < 4)
							{
								bsmat[0](loc1[0], loc0) = val;
								loc1[0]++;
							}
							else if (i1 >= 4 && j1 < 4 && k1 < 4)
							{
								bsmat[1](loc1[1], loc0) = val;
								loc1[1]++;
							}
							else if (i1 >= 4 && j1 >= 4 && k1 < 4)
							{
								bsmat[2](loc1[2], loc0) = val;
								loc1[2]++;
							}
							else if (i1 < 4 && j1 >= 4 && k1 < 4)
							{
								bsmat[3](loc1[3], loc0) = val;
								loc1[3]++;
							}
							else if (i1 < 4 && j1 < 4 && k1 >= 4)
							{
								bsmat[4](loc1[4], loc0) = val;
								loc1[4]++;
							}
							else if (i1 >= 4 && j1 < 4 && k1 >= 4)
							{
								bsmat[5](loc1[5], loc0) = val;
								loc1[5]++;
							}
							else if (i1 >= 4 && j1 >= 4 && k1 >= 4)
							{
								bsmat[6](loc1[6], loc0) = val;
								loc1[6]++;
							}
							else if (i1 < 4 && j1 >= 4 && k1 >= 4)
							{
								bsmat[7](loc1[7], loc0) = val;
								loc1[7]++;
							}
						}
					}
				}
				loc0++;
			}
		}
	}
}

void TruncatedTspline_3D::RefineTest_1()
{
	//the input is a cube with 9 elements in each direction
	//int i(4), j(4), k(4);
	//int rfid = 81 * k + 9 * j + i;
	int rfid(-1);
	for (uint i = 0; i < tmesh.size(); i++)
	{
		for (int j = 0; j < 8; j++)
		{
			if (tmesh[i].cnct[j] == 124)
			{
				rfid = i; break;
			}
			if (rfid != -1) break;
		}
	}
	ElementSubdivide_2(rfid,2);
	//VisualizeFaceMesh("../io/hex_out/cube9_face");
	//VisualizeTMesh("../io/hex_out/cube3_edge");
	//cout << "Done\n";
	//getchar();

	//CollectActives();
	//UpdateConnect();
	//FindEdgeTopoDirec();

	////int ed(183);
	////cout << "pid: " << tmedge[eaid[ed]].pt[0] << "type: " << tmedge[eaid[ed]].pn[0][0] << " next id: " << tmedge[eaid[ed]].pn[0][1] << "\n";
	////cout << "pid: " << tmedge[eaid[ed]].pt[1] << "type: " << tmedge[eaid[ed]].pn[1][0] << " next id: " << tmedge[eaid[ed]].pn[1][1] << "\n";
	////getchar();

	//FindKnotInterval();

	////for (uint i = 0; i < cp.size(); i++)
	////{
	////	if (i>63)
	////	{
	////		cout << "pid: " << i << "\n";
	////		cout << "uv edges: " << tmedge[cp[i].uved[0]].id_act << " " << tmedge[cp[i].uved[1]].id_act << " " << tmedge[cp[i].uved[2]].id_act << "\n";
	////		cout << "kitU: ";
	////		for (int j = 0; j < 4; j++) cout << cp[i].kitvU[j] << " ";
	////		cout << "\n";
	////		cout << "kitV: ";
	////		for (int j = 0; j < 4; j++) cout << cp[i].kitvV[j] << " ";
	////		cout << "\n";
	////		cout << "kitW: ";
	////		for (int j = 0; j < 4; j++) cout << cp[i].kitvW[j] << " ";
	////		cout << "\n";
	////		getchar();
	////	}
	////}

	//SetLocalCoorSystem();
	//FindIEN_PatchKV();
	//Update_IEN();

	//Truncation();

	//int itmp[8] = {0,3,12,15,48,51,60,63};
	//for (int i = 0; i < 8; i++)
	//{
	//	cp[itmp[i]].trun = 0;
	//}
	//for (uint i = 0; i < cp.size(); i++)
	//{
	//	if (cp[i].trun == 1)
	//	{
	//		cout << "trunid: " << i << "\n";
	//		cout << "tbf: ";
	//		for (uint j = 0; j < cp[i].tbf.size(); j++) cout << cp[i].tbf[j] << " ";
	//		cout << "\n";
	//		cout << "tc: ";
	//		for (uint j = 0; j < cp[i].tc.size(); j++) cout << cp[i].tc[j] << " ";
	//		cout << "\n";
	//		getchar();
	//	}
	//}

	//for (uint i = 0; i < tmesh.size(); i++)
	//{
	//	if ( tmesh[i].act == 1 && tmesh[i].type == 0)
	//	{
	//		cout << "eid: " << i << "\n";
	//		cout << "IEN size: " << tmesh[i].IEN.size() << "\n";
	//		for (uint j = 0; j < tmesh[i].IEN.size(); j++)
	//		{
	//			cout << "pid: "<<tmesh[i].IEN[j] << "\n";
	//			//cout << "rhx: " << cp[tmesh[i].IEN[j]].rhx << "\n";
	//			//cout << "uved: " << tmedge[cp[tmesh[i].IEN[j]].uved[0]].id_act << " " << tmedge[cp[tmesh[i].IEN[j]].uved[1]].id_act << " " << tmedge[cp[tmesh[i].IEN[j]].uved[2]].id_act << "\n";
	//			//cout << "uved[0] pre: " << tmedge[cp[tmesh[i].IEN[j]].uved[0]].pn[0][0] << " " << tmedge[cp[tmesh[i].IEN[j]].uved[0]].pn[0][1] << "\n";
	//			//cout << "uved[0] next: " << tmedge[cp[tmesh[i].IEN[j]].uved[0]].pn[1][0] << " " << tmedge[cp[tmesh[i].IEN[j]].uved[0]].pn[1][1] << "\n";
	//			//cout << "uved[1] pre: " << tmedge[cp[tmesh[i].IEN[j]].uved[1]].pn[0][0] << " " << tmedge[cp[tmesh[i].IEN[j]].uved[1]].pn[0][1] << "\n";
	//			//cout << "uved[1] next: " << tmedge[cp[tmesh[i].IEN[j]].uved[1]].pn[1][0] << " " << tmedge[cp[tmesh[i].IEN[j]].uved[1]].pn[1][1] << "\n";
	//			//cout << "uved[2] pre: " << tmedge[cp[tmesh[i].IEN[j]].uved[2]].pn[0][0] << " " << tmedge[cp[tmesh[i].IEN[j]].uved[2]].pn[0][1] << "\n";
	//			//cout << "uved[2] next: " << tmedge[cp[tmesh[i].IEN[j]].uved[2]].pn[1][0] << " " << tmedge[cp[tmesh[i].IEN[j]].uved[2]].pn[1][1] << "\n";
	//			//cout << "kitU: ";
	//			//for (int k = 0; k < 4; k++)
	//			//{
	//			//	cout << cp[tmesh[i].IEN[j]].kitvU[k] << " ";
	//			//}
	//			//cout << "\n";
	//			//cout << "kitV: ";
	//			//for (int k = 0; k < 4; k++)
	//			//{
	//			//	cout << cp[tmesh[i].IEN[j]].kitvV[k] << " ";
	//			//}
	//			//cout << "\n";
	//			//cout << "kitW: ";
	//			//for (int k = 0; k < 4; k++)
	//			//{
	//			//	cout << cp[tmesh[i].IEN[j]].kitvW[k] << " ";
	//			//}
	//			//cout << "\n";

	//			cout << "patch ku: ";
	//			for (int k = 0; k < 5; k++)
	//			{
	//				cout << tmesh[i].patch_ku[j][k] << " ";
	//			}
	//			cout << "\n";
	//			cout << "patch kv: ";
	//			for (int k = 0; k < 5; k++)
	//			{
	//				cout << tmesh[i].patch_kv[j][k] << " ";
	//			}
	//			cout << "\n";
	//			cout << "patch kw: ";
	//			for (int k = 0; k < 5; k++)
	//			{
	//				cout << tmesh[i].patch_kw[j][k] << " ";
	//			}
	//			cout << "\n";
	//			getchar();
	//		}
	//		cout << "end!\n";
	//		getchar();
	//	}
	//}
}

void TruncatedTspline_3D::RefineTest_2()
{
	//the input is a cube with 9 elements in each direction
	//int i(4), j(4), k(4);
	//int rfid = 81 * k + 9 * j + i;
	int rfid(-1);
	for (uint i = 0; i < tmesh.size(); i++)
	{
		for (int j = 0; j < 8; j++)
		{
			if (tmesh[i].cnct[j] == 124)
			{
				rfid = i; break;
			}
			if (rfid != -1) break;
		}
	}
	ElementSubdivide_2(rfid, 2);
	//VisualizeFaceMesh("../io/hex_out/cube9_face");
	//VisualizeTMesh("../io/hex_out/cube3_edge");
	//cout << "Done\n";
	//getchar();

	//CollectActives();
	//UpdateConnect();
	//FindEdgeTopoDirec();

	////int ed(183);
	////cout << "pid: " << tmedge[eaid[ed]].pt[0] << "type: " << tmedge[eaid[ed]].pn[0][0] << " next id: " << tmedge[eaid[ed]].pn[0][1] << "\n";
	////cout << "pid: " << tmedge[eaid[ed]].pt[1] << "type: " << tmedge[eaid[ed]].pn[1][0] << " next id: " << tmedge[eaid[ed]].pn[1][1] << "\n";
	////getchar();

	//FindKnotInterval();

	////for (uint i = 0; i < cp.size(); i++)
	////{
	////	if (i>63)
	////	{
	////		cout << "pid: " << i << "\n";
	////		cout << "uv edges: " << tmedge[cp[i].uved[0]].id_act << " " << tmedge[cp[i].uved[1]].id_act << " " << tmedge[cp[i].uved[2]].id_act << "\n";
	////		cout << "kitU: ";
	////		for (int j = 0; j < 4; j++) cout << cp[i].kitvU[j] << " ";
	////		cout << "\n";
	////		cout << "kitV: ";
	////		for (int j = 0; j < 4; j++) cout << cp[i].kitvV[j] << " ";
	////		cout << "\n";
	////		cout << "kitW: ";
	////		for (int j = 0; j < 4; j++) cout << cp[i].kitvW[j] << " ";
	////		cout << "\n";
	////		getchar();
	////	}
	////}

	//SetLocalCoorSystem();
	//FindIEN_PatchKV();
	//Update_IEN();

	//Truncation();
}

void TruncatedTspline_3D::PatchRefine_Regular(int lev, int eid)
{
	if (hmesh.size() < lev + 2)
	{
		vector<Vertex3D> ptmp;
		vector<Element3D> etmp;
		vector<Face3D> fctmp;
		vector<Edge3D> edtmp;
		hcp.push_back(ptmp);
		hmesh.push_back(etmp);
		hface.push_back(fctmp);
		hedge.push_back(edtmp);
	}
	double range[3][2] = { { hmesh[lev][eid].dm[0][0], hmesh[lev][eid].dm[0][1] }, { hmesh[lev][eid].dm[1][0], hmesh[lev][eid].dm[1][1] },
	{ hmesh[lev][eid].dm[2][0], hmesh[lev][eid].dm[2][1] } };
	uint n, i, j, k, loc(0);
	vector<array<double, 15>> pkuvw;
	vector<array<double, 3>> pnew;
	vector<vector<int>> chdid(hmesh[lev][eid].IEN.size());
	vector<vector<double>> coef(hmesh[lev][eid].IEN.size());
	//subdivision
	for (n = 0; n < hmesh[lev][eid].IEN.size(); n++)
	{
		int pid(hmesh[lev][eid].IEN[n]);
		vector<double> ku, kv, kw;
		BisectKnotInterval(hmesh[lev][eid].patch_ku[n], ku);
		BisectKnotInterval(hmesh[lev][eid].patch_kv[n], kv);
		BisectKnotInterval(hmesh[lev][eid].patch_kw[n], kw);
		vector<vector<double>> Tu, Tv, Tw;
		vector<double> ku0(hmesh[lev][eid].patch_ku[n].begin(), hmesh[lev][eid].patch_ku[n].end());
		vector<double> kv0(hmesh[lev][eid].patch_kv[n].begin(), hmesh[lev][eid].patch_kv[n].end());
		vector<double> kw0(hmesh[lev][eid].patch_kw[n].begin(), hmesh[lev][eid].patch_kw[n].end());
		TMatrix(ku0, ku, 3, Tu);
		TMatrix(kv0, kv, 3, Tv);
		TMatrix(kw0, kw, 3, Tw);
		for (k = 0; k < kw.size() - 4; k++)
		{
			for (j = 0; j < kv.size() - 4; j++)
			{
				for (i = 0; i < ku.size() - 4; i++)
				{
					if (ku[i + 4] > range[0][0] && ku[i] < range[0][1] && kv[j + 4] > range[1][0] && kv[j] < range[1][1] && kw[k + 4] > range[2][0] && kw[k] < range[2][1])
					{
						array<double, 15> ktmp = { ku[i], ku[i + 1], ku[i + 2], ku[i + 3], ku[i + 4],
							kv[j], kv[j + 1], kv[j + 2], kv[j + 3], kv[j + 4], kw[k], kw[k + 1], kw[k + 2], kw[k + 3], kw[k + 4] };
						vector<array<double, 15>>::iterator it1 = find(pkuvw.begin(), pkuvw.end(), ktmp);
						int pos = it1 - pkuvw.begin();
						double rfc (Tu[i][0] * Tv[j][0] * Tw[k][0]);
						chdid[n].push_back(pos);
						coef[n].push_back(rfc);
						if (it1 == pkuvw.end())
						{
							pkuvw.push_back(ktmp);
							array<double, 3> ptmp = { rfc*hcp[lev][pid].coor[0], rfc*hcp[lev][pid].coor[1], rfc*hcp[lev][pid].coor[2] };
							pnew.push_back(ptmp);
						}
						else
						{
							pnew[pos][0] += rfc*hcp[lev][pid].coor[0];
							pnew[pos][1] += rfc*hcp[lev][pid].coor[1];
							pnew[pos][2] += rfc*hcp[lev][pid].coor[2];
						}
					}
				}
			}
		}
	}
	double tol(1.e-6), dis(0.);
	vector<int> pnid(pnew.size());
	uint np0(hcp[lev + 1].size());
	for (i = 0; i < pnew.size(); i++)
	{
		int flag(0);
		for (j = 0; j < np0; j++)
		{
			dis = sqrt((pnew[i][0] - hcp[lev + 1][j].coor[0])*(pnew[i][0] - hcp[lev + 1][j].coor[0]) + (pnew[i][1] - hcp[lev + 1][j].coor[1])*(pnew[i][1] - hcp[lev + 1][j].coor[1])+
				(pnew[i][2] - hcp[lev + 1][j].coor[2])*(pnew[i][2] - hcp[lev + 1][j].coor[2]));
			if (dis < tol)
			{
				flag = 1; pnid[i] = j; break;
			}
		}
		if (flag == 0)
		{
			pnid[i] = hcp[lev + 1].size();
			Vertex3D cptmp;
			cptmp.coor[0] = pnew[i][0]; cptmp.coor[1] = pnew[i][1]; cptmp.coor[2] = pnew[i][2];
			cptmp.lev = lev + 1;
			hcp[lev + 1].push_back(cptmp);
		}
	}
	//children basis functions
	for (i = 0; i < chdid.size(); i++)
	{
		int pid(hmesh[lev][eid].IEN[i]);
		for (j = 0; j < chdid[i].size(); j++)
		{
			chdid[i][j] = pnid[chdid[i][j]];
			vector<int>::iterator it1 = find(hcp[lev][pid].chd.begin(), hcp[lev][pid].chd.end(), chdid[i][j]);
			if (it1 == hcp[lev][pid].chd.end())
			{
				hcp[lev][pid].chd.push_back(chdid[i][j]);
				hcp[lev][pid].coef.push_back(coef[i][j]);
			}
		}
	}
	//new elements
	double r1[3][3] = { { range[0][0], (range[0][0] + range[0][1]) / 2., range[0][1] }, { range[1][0], (range[1][0] + range[1][1]) / 2., range[1][1] },
	{ range[2][0], (range[2][0] + range[2][1]) / 2., range[2][1] }};
	double enew[8][3][2] = { { { r1[0][0], r1[0][1] }, { r1[1][0], r1[1][1] }, { r1[2][0], r1[2][1] } }, { { r1[0][1], r1[0][2] }, { r1[1][0], r1[1][1] }, { r1[2][0], r1[2][1] } },
	{ { r1[0][1], r1[0][2] }, { r1[1][1], r1[1][2] }, { r1[2][0], r1[2][1] } }, { { r1[0][0], r1[0][1] }, { r1[1][1], r1[1][2] }, { r1[2][0], r1[2][1] } }, 
	{ { r1[0][0], r1[0][1] }, { r1[1][0], r1[1][1] }, { r1[2][1], r1[2][2] } }, { { r1[0][1], r1[0][2] }, { r1[1][0], r1[1][1] }, { r1[2][1], r1[2][2] } },
	{ { r1[0][1], r1[0][2] }, { r1[1][1], r1[1][2] }, { r1[2][1], r1[2][2] } }, { { r1[0][0], r1[0][1] }, { r1[1][1], r1[1][2] }, { r1[2][1], r1[2][2] } }	};
	for (i = 0; i < 8; i++)
	{
		Element3D etmp;
		etmp.lev = lev + 1;
		etmp.prt = eid;
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 2; k++)
			{
				etmp.dm[j][k] = enew[i][j][k];
			}
		}
		for (j = 0; j < pkuvw.size(); j++)
		{
			if (pkuvw[j][0] <= enew[i][0][1] && pkuvw[j][4] >= enew[i][0][0] && pkuvw[j][5] <= enew[i][1][1] && pkuvw[j][9] >= enew[i][1][0] &&
				pkuvw[j][10] <= enew[i][2][1] && pkuvw[j][14] >= enew[i][2][0])
			{
				etmp.IEN.push_back(pnid[j]);
				array<double, 5> utmp = { pkuvw[j][0], pkuvw[j][1], pkuvw[j][2], pkuvw[j][3], pkuvw[j][4] };
				array<double, 5> vtmp = { pkuvw[j][5], pkuvw[j][6], pkuvw[j][7], pkuvw[j][8], pkuvw[j][9] };
				array<double, 5> wtmp = { pkuvw[j][10], pkuvw[j][11], pkuvw[j][12], pkuvw[j][13], pkuvw[j][14] };
				etmp.patch_ku.push_back(utmp);
				etmp.patch_kv.push_back(vtmp);
				etmp.patch_kw.push_back(wtmp);
				if (utmp[2] == enew[i][0][0] && utmp[3] == enew[i][0][1] && vtmp[2] == enew[i][1][0] && vtmp[3] == enew[i][1][1] &&
					wtmp[2] == enew[i][2][0] && wtmp[3] == enew[i][2][1])
				{
					etmp.cnct[0] = pnid[j];
				}
				else if (utmp[1] == enew[i][0][0] && utmp[2] == enew[i][0][1] && vtmp[2] == enew[i][1][0] && vtmp[3] == enew[i][1][1] &&
					wtmp[2] == enew[i][2][0] && wtmp[3] == enew[i][2][1])
				{
					etmp.cnct[1] = pnid[j];
				}
				else if (utmp[1] == enew[i][0][0] && utmp[2] == enew[i][0][1] && vtmp[1] == enew[i][1][0] && vtmp[2] == enew[i][1][1] &&
					wtmp[2] == enew[i][2][0] && wtmp[3] == enew[i][2][1])
				{
					etmp.cnct[2] = pnid[j];
				}
				else if (utmp[2] == enew[i][0][0] && utmp[3] == enew[i][0][1] && vtmp[1] == enew[i][1][0] && vtmp[2] == enew[i][1][1] &&
					wtmp[2] == enew[i][2][0] && wtmp[3] == enew[i][2][1])
				{
					etmp.cnct[3] = pnid[j];
				}
				else if (utmp[2] == enew[i][0][0] && utmp[3] == enew[i][0][1] && vtmp[2] == enew[i][1][0] && vtmp[3] == enew[i][1][1] &&
					wtmp[1] == enew[i][2][0] && wtmp[2] == enew[i][2][1])
				{
					etmp.cnct[4] = pnid[j];
				}
				else if (utmp[1] == enew[i][0][0] && utmp[2] == enew[i][0][1] && vtmp[2] == enew[i][1][0] && vtmp[3] == enew[i][1][1] &&
					wtmp[1] == enew[i][2][0] && wtmp[2] == enew[i][2][1])
				{
					etmp.cnct[5] = pnid[j];
				}
				else if (utmp[1] == enew[i][0][0] && utmp[2] == enew[i][0][1] && vtmp[1] == enew[i][1][0] && vtmp[2] == enew[i][1][1] &&
					wtmp[1] == enew[i][2][0] && wtmp[2] == enew[i][2][1])
				{
					etmp.cnct[6] = pnid[j];
				}
				else if (utmp[2] == enew[i][0][0] && utmp[3] == enew[i][0][1] && vtmp[1] == enew[i][1][0] && vtmp[2] == enew[i][1][1] &&
					wtmp[1] == enew[i][2][0] && wtmp[2] == enew[i][2][1])
				{
					etmp.cnct[7] = pnid[j];
				}
			}
		}
		hmesh[lev][eid].chd.push_back(hmesh[lev+1].size());
		hmesh[lev + 1].push_back(etmp);
	}

	hmesh[lev][eid].act = 0;
}

void TruncatedTspline_3D::BisectKnotInterval(const array<double, 5>& kv_in, vector<double>& kv_out)
{
	kv_out.clear();
	for (int i = 0; i < 4; i++)
	{
		kv_out.push_back(kv_in[i]);
		if (kv_in[i] < kv_in[i + 1])
		{
			kv_out.push_back((kv_in[i] + kv_in[i+1]) / 2.);
		}
	}
	kv_out.push_back(kv_in[4]);
}

void TruncatedTspline_3D::PatchRefine_Irregular(int lev, int eid)
{
	if (hmesh.size() < lev + 2)
	{
		vector<Vertex3D> ptmp;
		vector<Element3D> etmp;
		vector<Face3D> fctmp;
		vector<Edge3D> edtmp;
		hcp.push_back(ptmp);
		hmesh.push_back(etmp);
		hface.push_back(fctmp);
		hedge.push_back(edtmp);
	}

	CatmullClark(lev, eid);

	//double range[3][2] = { { hmesh[lev][eid].dm[0][0], hmesh[lev][eid].dm[0][1] }, { hmesh[lev][eid].dm[1][0], hmesh[lev][eid].dm[1][1] },
	//{ hmesh[lev][eid].dm[2][0], hmesh[lev][eid].dm[2][1] } };
	//vector<double> ku0(8), kv0(8), kw0(8), ku1(9), kv1(9), kw1(9);
	//uint n, i, j, k;
	//for (i = 0; i < 4; i++)
	//{
	//	ku0[i] = range[0][0]; ku0[i + 4] = range[0][1];
	//	kv0[i] = range[1][0]; kv0[i + 4] = range[1][1];
	//	kw0[i] = range[2][0]; kw0[i + 4] = range[2][1];
	//	ku1[i] = range[0][0]; ku1[i + 5] = range[0][1];
	//	kv1[i] = range[1][0]; kv1[i + 5] = range[1][1];
	//	kw1[i] = range[2][0]; kw1[i + 5] = range[2][1];
	//}
	//ku1[4] = (range[0][0] + range[0][1]) / 2.;
	//kv1[4] = (range[1][0] + range[1][1]) / 2.;
	//kw1[4] = (range[2][0] + range[2][1]) / 2.;
	//vector<vector<double>> Tu, Tv, Tw;
	//TMatrix(ku0, ku1, 3, Tu);
	//TMatrix(kv0, kv1, 3, Tv);
	//TMatrix(kw0, kw1, 3, Tw);
	////calculate new control points using Catmull-Clark
	//vector<array<double, 3>> pnew;
	//vector<vector<int>> chdid(hmesh[lev][eid].IEN.size());
	//vector<vector<double>> coef(hmesh[lev][eid].IEN.size());


	//double r1[3][3] = { { range[0][0], (range[0][0] + range[0][1]) / 2., range[0][1] }, { range[1][0], (range[1][0] + range[1][1]) / 2., range[1][1] },
	//{ range[2][0], (range[2][0] + range[2][1]) / 2., range[2][1] } };
	//double e_dm[8][3][2] = { { { r1[0][0], r1[0][1] }, { r1[1][0], r1[1][1] }, { r1[2][0], r1[2][1] } }, { { r1[0][1], r1[0][2] }, { r1[1][0], r1[1][1] }, { r1[2][0], r1[2][1] } },
	//{ { r1[0][1], r1[0][2] }, { r1[1][1], r1[1][2] }, { r1[2][0], r1[2][1] } }, { { r1[0][0], r1[0][1] }, { r1[1][1], r1[1][2] }, { r1[2][0], r1[2][1] } },
	//{ { r1[0][0], r1[0][1] }, { r1[1][0], r1[1][1] }, { r1[2][1], r1[2][2] } }, { { r1[0][1], r1[0][2] }, { r1[1][0], r1[1][1] }, { r1[2][1], r1[2][2] } },
	//{ { r1[0][1], r1[0][2] }, { r1[1][1], r1[1][2] }, { r1[2][1], r1[2][2] } }, { { r1[0][0], r1[0][1] }, { r1[1][1], r1[1][2] }, { r1[2][1], r1[2][2] } } };
}

void TruncatedTspline_3D::CatmullClark(int lev, int eid)
{
	uint i, j, k;
	vector<int> bdid, fcid, edid;
	bdid.push_back(eid);
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].hex.size(); j++)
		{
			vector<int>::iterator it = find(bdid.begin(), bdid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].hex[j]);
			if (it == bdid.end())
			{
				bdid.push_back(hcp[lev][hmesh[lev][eid].cnct[i]].hex[j]);
			}
		}
		for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].face.size(); j++)
		{
			vector<int>::iterator it = find(fcid.begin(), fcid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].face[j]);
			if (it == fcid.end())
			{
				fcid.push_back(hcp[lev][hmesh[lev][eid].cnct[i]].face[j]);
			}
		}
		for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].edge.size(); j++)
		{
			vector<int>::iterator it = find(edid.begin(), edid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].edge[j]);
			if (it == edid.end())
			{
				edid.push_back(hcp[lev][hmesh[lev][eid].cnct[i]].edge[j]);
			}
		}
	}
	//calculate coordinates of new points
	vector<array<double, 3>> pn_bd(bdid.size()), pn_fc(fcid.size()), pn_ed(edid.size()), pn_vt(8);
	vector<array<double, 3>> pn_fc0(fcid.size()), pn_ed0(edid.size());
	for (i = 0; i < bdid.size(); i++)
	{
		pn_bd[i][0] = 0.; pn_bd[i][1] = 0.; pn_bd[i][2] = 0.;
		for (j = 0; j < 8; j++)
		{
			pn_bd[i][0] += .125*hcp[lev][hmesh[lev][bdid[i]].cnct[j]].coor[0];
			pn_bd[i][1] += .125*hcp[lev][hmesh[lev][bdid[i]].cnct[j]].coor[1];
			pn_bd[i][2] += .125*hcp[lev][hmesh[lev][bdid[i]].cnct[j]].coor[2];
		}
	}
	for (i = 0; i < fcid.size(); i++)
	{
		pn_fc0[i][0] = 0.; pn_fc0[i][1] = 0.; pn_fc0[i][2] = 0.;
		for (j = 0; j < 4; j++)
		{
			pn_fc0[i][0] += .25*hcp[lev][hface[lev][fcid[i]].cnct[j]].coor[0];
			pn_fc0[i][1] += .25*hcp[lev][hface[lev][fcid[i]].cnct[j]].coor[1];
			pn_fc0[i][2] += .25*hcp[lev][hface[lev][fcid[i]].cnct[j]].coor[2];
		}
		vector<int>::iterator it1 = find(bdid.begin(), bdid.end(), hface[lev][fcid[i]].hex[0]);
		vector<int>::iterator it2 = find(bdid.begin(), bdid.end(), hface[lev][fcid[i]].hex[1]);
		int pos1(it1 - bdid.begin()), pos2(it2 - bdid.begin());
		pn_fc[i][0] = (pn_bd[pos1][0] + pn_bd[pos2][0] + 2.*pn_fc0[i][0]) / 4.;
		pn_fc[i][1] = (pn_bd[pos1][1] + pn_bd[pos2][1] + 2.*pn_fc0[i][1]) / 4.;
		pn_fc[i][2] = (pn_bd[pos1][2] + pn_bd[pos2][2] + 2.*pn_fc0[i][2]) / 4.;
	}
	for (i = 0; i < edid.size(); i++)
	{
		pn_ed0[i][0] = (hcp[lev][hedge[lev][edid[i]].pt[0]].coor[0] + hcp[lev][hedge[lev][edid[i]].pt[1]].coor[0]) / 2.;
		pn_ed0[i][1] = (hcp[lev][hedge[lev][edid[i]].pt[0]].coor[1] + hcp[lev][hedge[lev][edid[i]].pt[1]].coor[1]) / 2.;
		pn_ed0[i][2] = (hcp[lev][hedge[lev][edid[i]].pt[0]].coor[2] + hcp[lev][hedge[lev][edid[i]].pt[1]].coor[2]) / 2.;
		array<double, 3> f_avg = { 0., 0., 0. }, b_avg = {0.,0.,0.};
		for (j = 0; j < hedge[lev][edid[i]].hex.size(); j++)
		{
			vector<int>::iterator it = find(bdid.begin(), bdid.end(), hedge[lev][edid[i]].hex[j]);
			int pos(it - bdid.begin());
			b_avg[0] += pn_bd[pos][0]; b_avg[1] += pn_bd[pos][1]; b_avg[2] += pn_bd[pos][2];
		}
		b_avg[0] /= hedge[lev][edid[i]].hex.size(); b_avg[1] /= hedge[lev][edid[i]].hex.size(); b_avg[2] /= hedge[lev][edid[i]].hex.size();
		for (j = 0; j < hedge[lev][edid[i]].face.size(); j++)
		{
			vector<int>::iterator it = find(fcid.begin(), fcid.end(), hedge[lev][edid[i]].face[j]);
			int pos(it - fcid.begin());
			f_avg[0] += pn_fc0[pos][0]; f_avg[1] += pn_fc0[pos][1]; f_avg[2] += pn_fc0[pos][2];
		}
		f_avg[0] /= hedge[lev][edid[i]].face.size(); f_avg[1] /= hedge[lev][edid[i]].face.size(); f_avg[2] /= hedge[lev][edid[i]].face.size();
		uint n = hedge[lev][edid[i]].face.size();
		pn_ed[i][0] = (b_avg[0] + 2.*f_avg[0] + (n - 3)*pn_ed0[i][0]) / n;
		pn_ed[i][1] = (b_avg[1] + 2.*f_avg[1] + (n - 3)*pn_ed0[i][1]) / n;
		pn_ed[i][2] = (b_avg[2] + 2.*f_avg[2] + (n - 3)*pn_ed0[i][2]) / n;
	}
	for (i = 0; i < 8; i++)
	{
		pn_vt[i][0] = hcp[lev][hmesh[lev][eid].cnct[i]].coor[0];
		pn_vt[i][1] = hcp[lev][hmesh[lev][eid].cnct[i]].coor[1];
		pn_vt[i][2] = hcp[lev][hmesh[lev][eid].cnct[i]].coor[2];
		array<double, 3> b_avg = { 0., 0., 0. }, f_avg = { 0., 0., 0. }, e_avg = { 0., 0., 0. };
		uint n(hcp[lev][hmesh[lev][eid].cnct[i]].hex.size());
		for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].hex.size(); j++)
		{
			vector<int>::iterator it = find(bdid.begin(), bdid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].hex[j]);
			int pos(it - bdid.begin());
			b_avg[0] += pn_bd[pos][0]; b_avg[1] += pn_bd[pos][1]; b_avg[2] += pn_bd[pos][2];
		}
		b_avg[0] /= n; b_avg[1] /= n; b_avg[2] /= n;
		n = hcp[lev][hmesh[lev][eid].cnct[i]].face.size();
		for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].face.size(); j++)
		{
			vector<int>::iterator it = find(fcid.begin(), fcid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].face[j]);
			int pos(it - fcid.begin());
			f_avg[0] += pn_fc0[pos][0]; f_avg[1] += pn_fc0[pos][1]; f_avg[2] += pn_fc0[pos][2];
		}
		f_avg[0] /= n; f_avg[1] /= n; f_avg[2] /= n;
		n = hcp[lev][hmesh[lev][eid].cnct[i]].edge.size();
		for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].edge.size(); j++)
		{
			vector<int>::iterator it = find(edid.begin(), edid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].edge[j]);
			int pos(it - edid.begin());
			e_avg[0] += pn_ed0[pos][0]; e_avg[1] += pn_ed0[pos][1]; e_avg[2] += pn_ed0[pos][2];
		}
		e_avg[0] /= n; e_avg[1] /= n; e_avg[2] /= n;
		pn_vt[i][0] = (b_avg[0] + 3.*f_avg[0] + 3.*e_avg[0] + pn_vt[i][0]) / 8.;
		pn_vt[i][1] = (b_avg[1] + 3.*f_avg[1] + 3.*e_avg[1] + pn_vt[i][1]) / 8.;
		pn_vt[i][2] = (b_avg[2] + 3.*f_avg[2] + 3.*e_avg[2] + pn_vt[i][2]) / 8.;
	}
	//find global index of new points at Level+1
	vector<int> pnid_bd(pn_bd.size()), pnid_fc(pn_fc.size()), pnid_ed(pn_ed.size()), pnid_vt(pn_vt.size());
	double tol(1.e-6), dis(0.);
	uint np0(hcp[lev + 1].size());
	for (i = 0; i < pn_bd.size(); i++)
	{
		int flag(0);
		for (j = 0; j < np0; j++)
		{
			dis = sqrt((pn_bd[i][0] - hcp[lev + 1][j].coor[0])*(pn_bd[i][0] - hcp[lev + 1][j].coor[0]) + (pn_bd[i][1] - hcp[lev + 1][j].coor[1])*(pn_bd[i][1] - hcp[lev + 1][j].coor[1]) +
				(pn_bd[i][2] - hcp[lev + 1][j].coor[2])*(pn_bd[i][2] - hcp[lev + 1][j].coor[2]));
			if (dis < tol)
			{
				flag = 1; pnid_bd[i] = j; break;
			}
		}
		if (flag == 0)
		{
			pnid_bd[i] = hcp[lev + 1].size();
			Vertex3D cptmp;
			cptmp.coor[0] = pn_bd[i][0]; cptmp.coor[1] = pn_bd[i][1]; cptmp.coor[2] = pn_bd[i][2];
			cptmp.lev = lev + 1;
			hcp[lev + 1].push_back(cptmp);
		}
	}
	for (i = 0; i < pn_fc.size(); i++)
	{
		int flag(0);
		for (j = 0; j < np0; j++)
		{
			dis = sqrt((pn_fc[i][0] - hcp[lev + 1][j].coor[0])*(pn_fc[i][0] - hcp[lev + 1][j].coor[0]) + (pn_fc[i][1] - hcp[lev + 1][j].coor[1])*(pn_fc[i][1] - hcp[lev + 1][j].coor[1]) +
				(pn_fc[i][2] - hcp[lev + 1][j].coor[2])*(pn_fc[i][2] - hcp[lev + 1][j].coor[2]));
			if (dis < tol)
			{
				flag = 1; pnid_fc[i] = j; break;
			}
		}
		if (flag == 0)
		{
			pnid_fc[i] = hcp[lev + 1].size();
			Vertex3D cptmp;
			cptmp.coor[0] = pn_fc[i][0]; cptmp.coor[1] = pn_fc[i][1]; cptmp.coor[2] = pn_fc[i][2];
			cptmp.lev = lev + 1;
			hcp[lev + 1].push_back(cptmp);
		}
	}
	for (i = 0; i < pn_ed.size(); i++)
	{
		int flag(0);
		for (j = 0; j < np0; j++)
		{
			dis = sqrt((pn_ed[i][0] - hcp[lev + 1][j].coor[0])*(pn_ed[i][0] - hcp[lev + 1][j].coor[0]) + (pn_ed[i][1] - hcp[lev + 1][j].coor[1])*(pn_ed[i][1] - hcp[lev + 1][j].coor[1]) +
				(pn_ed[i][2] - hcp[lev + 1][j].coor[2])*(pn_ed[i][2] - hcp[lev + 1][j].coor[2]));
			if (dis < tol)
			{
				flag = 1; pnid_ed[i] = j; break;
			}
		}
		if (flag == 0)
		{
			pnid_ed[i] = hcp[lev + 1].size();
			Vertex3D cptmp;
			cptmp.coor[0] = pn_ed[i][0]; cptmp.coor[1] = pn_ed[i][1]; cptmp.coor[2] = pn_ed[i][2];
			cptmp.lev = lev + 1;
			hcp[lev + 1].push_back(cptmp);
		}
	}
	for (i = 0; i < pn_vt.size(); i++)
	{
		int flag(0);
		for (j = 0; j < np0; j++)
		{
			dis = sqrt((pn_vt[i][0] - hcp[lev + 1][j].coor[0])*(pn_vt[i][0] - hcp[lev + 1][j].coor[0]) + (pn_vt[i][1] - hcp[lev + 1][j].coor[1])*(pn_vt[i][1] - hcp[lev + 1][j].coor[1]) +
				(pn_vt[i][2] - hcp[lev + 1][j].coor[2])*(pn_vt[i][2] - hcp[lev + 1][j].coor[2]));
			if (dis < tol)
			{
				flag = 1; pnid_vt[i] = j; break;
			}
		}
		if (flag == 0)
		{
			pnid_vt[i] = hcp[lev + 1].size();
			Vertex3D cptmp;
			cptmp.coor[0] = pn_vt[i][0]; cptmp.coor[1] = pn_vt[i][1]; cptmp.coor[2] = pn_vt[i][2];
			cptmp.lev = lev + 1;
			hcp[lev + 1].push_back(cptmp);
		}
	}
	//find children id and corresponding coefficients
	//vector<vector<int>> chdid(hmesh[lev][eid].IEN.size());
	//vector<vector<double>> coef(hmesh[lev][eid].IEN.size());
	vector<array<int, 8>> bdmp(bdid.size());
	vector<array<double, 8>> bdc(bdid.size());
	vector<vector<int>> fcmp(fcid.size());
	vector<vector<double>> fcc(fcid.size());
	vector<vector<int>> edmp(edid.size());
	vector<vector<double>> edc(edid.size());
	vector<vector<int>> vtmp(8);
	vector<vector<double>> vtc(8);
	//body points
	for (i = 0; i < bdid.size(); i++)
	{
		for (j = 0; j < 8; j++)
		{
			bdmp[i][j]=hmesh[lev][bdid[i]].cnct[j];
			bdc[i][j] = .125;
		}
	}
	//face points
	for (i = 0; i < fcid.size(); i++)
	{
		for (j = 0; j < 4; j++)
		{
			fcmp[i].push_back(hface[lev][fcid[i]].cnct[j]);
			fcc[i].push_back(2.*.25/4.);
		}
		for (j = 0; j < hface[lev][fcid[i]].hex.size(); j++)
		{
			int hxid(hface[lev][fcid[i]].hex[j]);
			for (k = 0; k < 8; k++)
			{
				vector<int>::iterator it = find(fcmp[i].begin(),fcmp[i].end(),hmesh[lev][hxid].cnct[k]);
				double ctmp(.125/4.);
				if (it == fcmp[i].end())
				{
					fcmp[i].push_back(hmesh[lev][hxid].cnct[k]);
					fcc[i].push_back(ctmp);
				}
				else
				{
					int pos(it - fcmp[i].begin());
					fcc[i][pos] += ctmp;
				}
			}
		}
	}
	//cout
	//for (i = 0; i < fcc.size(); i++)
	//{
	//	double sum(0.);
	//	for (j = 0; j < fcc[i].size(); j++)
	//	{
	//		sum += fcc[i][j];
	//	}
	//	if (sum != 1.)
	//	{
	//		cout << "face point sum not 1\n";
	//		getchar();
	//	}
	//}
	//edge points
	for (i = 0; i < edid.size(); i++)
	{
		uint n(hedge[lev][edid[i]].hex.size());
		for (j = 0; j < 2; j++)
		{
			edmp[i].push_back(hedge[lev][edid[i]].pt[j]);
			edc[i].push_back(.5*double(n-3)/n);
		}
		for (j = 0; j < hedge[lev][edid[i]].hex.size(); j++)
		{
			int hxid(hedge[lev][edid[i]].hex[j]);
			for (k = 0; k < 8; k++)
			{
				vector<int>::iterator it = find(edmp[i].begin(), edmp[i].end(), hmesh[lev][hxid].cnct[k]);
				double ctmp(.125 / double(n*n));
				if (it == edmp[i].end())
				{
					edmp[i].push_back(hmesh[lev][hxid].cnct[k]);
					edc[i].push_back(ctmp);
				}
				else
				{
					int pos(it - edmp[i].begin());
					edc[i][pos] += ctmp;
				}
			}
		}
		for (j = 0; j < hedge[lev][edid[i]].face.size(); j++)
		{
			int fcid(hedge[lev][edid[i]].face[j]);
			for (k = 0; k < 4; k++)
			{
				vector<int>::iterator it = find(edmp[i].begin(), edmp[i].end(), hface[lev][fcid].cnct[k]);
				double ctmp(2.*.25 / double(n*n));
				if (it == edmp[i].end())
				{
					edmp[i].push_back(hface[lev][fcid].cnct[k]);
					edc[i].push_back(ctmp);
				}
				else
				{
					int pos(it - edmp[i].begin());
					edc[i][pos] += ctmp;
				}
			}
		}
	}
	//cout
	//for (i = 0; i < edc.size(); i++)
	//{
	//	double sum(0.), tol(1.e-10);
	//	for (j = 0; j < edc[i].size(); j++)
	//	{
	//		sum += edc[i][j];
	//	}
	//	if (fabs(sum-1.)>tol)
	//	{
	//		cout << "edge point " << sum <<"\n";
	//		//cout << "edge valence hex: " << edvl_hx[i] << "\n";
	//		////cout << "edge valence face: " << edvl_fc[i] << "\n";
	//		//cout << "# contribution: " << edmp[i].size() << "\n";
	//		//for (j = 0; j < edmp[i].size(); j++)
	//		//{
	//		//	cout << setw(6)<< edmp[i][j] << " ";
	//		//}
	//		//cout << "\n";
	//		//for (j = 0; j < edc[i].size(); j++)
	//		//{
	//		//	cout << setw(6) << edc[i][j] << " ";
	//		//}
	//		//cout << "\n";
	//		getchar();
	//	}
	//}
	//vertex points
	for (i = 0; i < 8; i++)
	{
		int ptid(hmesh[lev][eid].cnct[i]);
		uint nhx(hcp[lev][ptid].hex.size()), nfc(hcp[lev][ptid].face.size()), ned(hcp[lev][ptid].edge.size());
		vtmp[i].push_back(ptid);
		vtc[i].push_back(.125);
		for (j = 0; j < nhx; j++)
		{
			int hxid(hcp[lev][ptid].hex[j]);
			for (k = 0; k < 8; k++)
			{
				vector<int>::iterator it = find(vtmp[i].begin(),vtmp[i].end(),hmesh[lev][hxid].cnct[k]);
				double ctmp(.125/(8.*nhx));
				if (it == vtmp[i].end())
				{
					vtmp[i].push_back(hmesh[lev][hxid].cnct[k]);
					vtc[i].push_back(ctmp);
				}
				else
				{
					int pos(it - vtmp[i].begin());
					vtc[i][pos] += ctmp;
				}
			}
		}
		for (j = 0; j < nfc; j++)
		{
			int fcid(hcp[lev][ptid].face[j]);
			for (k = 0; k < 4; k++)
			{
				vector<int>::iterator it = find(vtmp[i].begin(), vtmp[i].end(), hface[lev][fcid].cnct[k]);
				double ctmp(3.*.25 / (8.*nfc));
				if (it == vtmp[i].end())
				{
					vtmp[i].push_back(hface[lev][fcid].cnct[k]);
					vtc[i].push_back(ctmp);
				}
				else
				{
					int pos(it - vtmp[i].begin());
					vtc[i][pos] += ctmp;
				}
			}
		}
		for (j = 0; j < ned; j++)
		{
			int edid(hcp[lev][ptid].edge[j]);
			for (k = 0; k < 2; k++)
			{
				vector<int>::iterator it = find(vtmp[i].begin(), vtmp[i].end(), hedge[lev][edid].pt[k]);
				double ctmp(3.*.5 / (8.*ned));
				if (it == vtmp[i].end())
				{
					vtmp[i].push_back(hedge[lev][edid].pt[k]);
					vtc[i].push_back(ctmp);
				}
				else
				{
					int pos(it - vtmp[i].begin());
					vtc[i][pos] += ctmp;
				}
			}
		}
	}
	//cout
	//for (i = 0; i < vtc.size(); i++)
	//{
	//	double sum(0.);
	//	for (j = 0; j < vtc[i].size(); j++)
	//	{
	//		sum += vtc[i][j];
	//	}
	//	if (fabs(sum - 1.)>tol)
	//	{
	//		cout << "vertex point " << sum << "\n";
	//		getchar();
	//	}
	//}
	//children and coef
	for (i = 0; i < bdmp.size(); i++)
	{
		for (j = 0; j < bdmp[i].size(); j++)
		{
			vector<int>::iterator it = find(hcp[lev][bdmp[i][j]].chd.begin(), hcp[lev][bdmp[i][j]].chd.end(), pnid_bd[i]);
			if (it == hcp[lev][bdmp[i][j]].chd.end())
			{
				hcp[lev][bdmp[i][j]].chd.push_back(pnid_bd[i]);
				hcp[lev][bdmp[i][j]].coef.push_back(bdc[i][j]);
			}
		}
	}
	for (i = 0; i < fcmp.size(); i++)
	{
		for (j = 0; j < fcmp[i].size(); j++)
		{
			vector<int>::iterator it = find(hcp[lev][fcmp[i][j]].chd.begin(), hcp[lev][fcmp[i][j]].chd.end(), pnid_fc[i]);
			if (it == hcp[lev][fcmp[i][j]].chd.end())
			{
				hcp[lev][fcmp[i][j]].chd.push_back(pnid_fc[i]);
				hcp[lev][fcmp[i][j]].coef.push_back(fcc[i][j]);
			}
		}
	}
	for (i = 0; i < edmp.size(); i++)
	{
		for (j = 0; j < edmp[i].size(); j++)
		{
			vector<int>::iterator it = find(hcp[lev][edmp[i][j]].chd.begin(), hcp[lev][edmp[i][j]].chd.end(), pnid_ed[i]);
			if (it == hcp[lev][edmp[i][j]].chd.end())
			{
				hcp[lev][edmp[i][j]].chd.push_back(pnid_ed[i]);
				hcp[lev][edmp[i][j]].coef.push_back(edc[i][j]);
			}
		}
	}
	for (i = 0; i < vtmp.size(); i++)
	{
		for (j = 0; j < vtmp[i].size(); j++)
		{
			vector<int>::iterator it = find(hcp[lev][vtmp[i][j]].chd.begin(), hcp[lev][vtmp[i][j]].chd.end(), pnid_vt[i]);
			if (it == hcp[lev][vtmp[i][j]].chd.end())
			{
				hcp[lev][vtmp[i][j]].chd.push_back(pnid_vt[i]);
				hcp[lev][vtmp[i][j]].coef.push_back(vtc[i][j]);
			}
		}
	}
	//local index of new points
	vector<int> ed_loc(12), fc_loc(6), vt_loc(pnid_vt);
	int bd_loc(pnid_bd[0]);
	for (i = 0; i < 12; i++)
	{
		vector<int>::iterator it = find(edid.begin(), edid.end(), hmesh[lev][eid].edge[i]);
		int pos(it - edid.begin());
		ed_loc[i] = pnid_ed[pos];
	}
	for (i = 0; i < 6; i++)
	{
		vector<int>::iterator it = find(fcid.begin(), fcid.end(), hmesh[lev][eid].face[i]);
		int pos(it - fcid.begin());
		fc_loc[i] = pnid_fc[pos];
	}
	//construct new hex
	int enew[8][8] = { { vt_loc[0], ed_loc[0], fc_loc[0], ed_loc[3], ed_loc[4], fc_loc[1], bd_loc, fc_loc[4] },
	{ ed_loc[0], vt_loc[1], ed_loc[1], fc_loc[0], fc_loc[1], ed_loc[5], fc_loc[2], bd_loc },
	{ fc_loc[0], ed_loc[1], vt_loc[2], ed_loc[2], bd_loc, fc_loc[2], ed_loc[6], fc_loc[3] },
	{ ed_loc[3], fc_loc[0], ed_loc[2], vt_loc[3], fc_loc[4], bd_loc, fc_loc[3], ed_loc[7] },
	{ ed_loc[4], fc_loc[1], bd_loc, fc_loc[4], vt_loc[4], ed_loc[8], fc_loc[5], ed_loc[11] },
	{ fc_loc[1], ed_loc[5], fc_loc[2], bd_loc, ed_loc[8], vt_loc[5], ed_loc[9], fc_loc[5] },
	{ bd_loc, fc_loc[2], ed_loc[6], fc_loc[3], fc_loc[5], ed_loc[9], vt_loc[6], ed_loc[10] },
	{ fc_loc[4], bd_loc, fc_loc[3], ed_loc[7], ed_loc[11], fc_loc[5], ed_loc[10], vt_loc[7] } };

	for (i = 0; i < 8; i++)
	{
		Element3D etmp;
		etmp.type = 2;//tmp, determined by the topology at level (lev+1)
		etmp.lev = lev + 1;
		etmp.prt = eid;
		for (j = 0; j < 8; j++)
		{
			etmp.cnct[j] = enew[i][j];
		}
		hmesh[lev][eid].chd.push_back(hmesh[lev+1].size());
		hmesh[lev + 1].push_back(etmp);
	}

	hmesh[lev][eid].act = 0;
}

void TruncatedTspline_3D::PatchRefine_Boundary(int lev, int eid)
{
	if (hmesh.size() < lev + 2)
	{
		vector<Vertex3D> ptmp;
		vector<Element3D> etmp;
		vector<Face3D> fctmp;
		vector<Edge3D> edtmp;
		hcp.push_back(ptmp);
		hmesh.push_back(etmp);
		hface.push_back(fctmp);
		hedge.push_back(edtmp);
	}

	uint i, j, k;
	vector<int> bdid, fcid, edid;
	bdid.push_back(eid);
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].hex.size(); j++)
		{
			vector<int>::iterator it = find(bdid.begin(), bdid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].hex[j]);
			if (it == bdid.end())
			{
				bdid.push_back(hcp[lev][hmesh[lev][eid].cnct[i]].hex[j]);
			}
		}
		for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].face.size(); j++)
		{
			vector<int>::iterator it = find(fcid.begin(), fcid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].face[j]);
			if (it == fcid.end())
			{
				fcid.push_back(hcp[lev][hmesh[lev][eid].cnct[i]].face[j]);
			}
		}
		for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].edge.size(); j++)
		{
			vector<int>::iterator it = find(edid.begin(), edid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].edge[j]);
			if (it == edid.end())
			{
				edid.push_back(hcp[lev][hmesh[lev][eid].cnct[i]].edge[j]);
			}
		}
	}
	//calculate coordinates of new points
	vector<array<double, 3>> pn_bd(bdid.size()), pn_fc(fcid.size()), pn_ed(edid.size()), pn_vt(8);
	vector<array<double, 3>> pn_fc0(fcid.size()), pn_ed0(edid.size());
	vector<int> fc_type(fcid.size(),0), ed_type(edid.size(),0), vt_type(8,0);
	vector<int> ed_shp(edid.size(),0), vt_shp(8,0);
	for (i = 0; i < bdid.size(); i++)
	{
		pn_bd[i][0] = 0.; pn_bd[i][1] = 0.; pn_bd[i][2] = 0.;
		for (j = 0; j < 8; j++)
		{
			pn_bd[i][0] += .125*hcp[lev][hmesh[lev][bdid[i]].cnct[j]].coor[0];
			pn_bd[i][1] += .125*hcp[lev][hmesh[lev][bdid[i]].cnct[j]].coor[1];
			pn_bd[i][2] += .125*hcp[lev][hmesh[lev][bdid[i]].cnct[j]].coor[2];
		}
	}
	for (i = 0; i < fcid.size(); i++)
	{
		pn_fc0[i][0] = 0.; pn_fc0[i][1] = 0.; pn_fc0[i][2] = 0.;
		for (j = 0; j < 4; j++)
		{
			pn_fc0[i][0] += .25*hcp[lev][hface[lev][fcid[i]].cnct[j]].coor[0];
			pn_fc0[i][1] += .25*hcp[lev][hface[lev][fcid[i]].cnct[j]].coor[1];
			pn_fc0[i][2] += .25*hcp[lev][hface[lev][fcid[i]].cnct[j]].coor[2];
		}
		if (hface[lev][fcid[i]].type==0)
		{
			vector<int>::iterator it1 = find(bdid.begin(), bdid.end(), hface[lev][fcid[i]].hex[0]);
			vector<int>::iterator it2 = find(bdid.begin(), bdid.end(), hface[lev][fcid[i]].hex[1]);
			int pos1(it1 - bdid.begin()), pos2(it2 - bdid.begin());
			pn_fc[i][0] = (pn_bd[pos1][0] + pn_bd[pos2][0] + 2.*pn_fc0[i][0]) / 4.;
			pn_fc[i][1] = (pn_bd[pos1][1] + pn_bd[pos2][1] + 2.*pn_fc0[i][1]) / 4.;
			pn_fc[i][2] = (pn_bd[pos1][2] + pn_bd[pos2][2] + 2.*pn_fc0[i][2]) / 4.;
		}
		else if (hface[lev][fcid[i]].type == 1)
		{
			pn_fc[i][0] = pn_fc0[i][0];
			pn_fc[i][1] = pn_fc0[i][1];
			pn_fc[i][2] = pn_fc0[i][2];
			fc_type[i] = 1;
		}
	}
	for (i = 0; i < edid.size(); i++)
	{
		pn_ed0[i][0] = (hcp[lev][hedge[lev][edid[i]].pt[0]].coor[0] + hcp[lev][hedge[lev][edid[i]].pt[1]].coor[0]) / 2.;
		pn_ed0[i][1] = (hcp[lev][hedge[lev][edid[i]].pt[0]].coor[1] + hcp[lev][hedge[lev][edid[i]].pt[1]].coor[1]) / 2.;
		pn_ed0[i][2] = (hcp[lev][hedge[lev][edid[i]].pt[0]].coor[2] + hcp[lev][hedge[lev][edid[i]].pt[1]].coor[2]) / 2.;
		if (hedge[lev][edid[i]].sharp == 1)//boundary sharp edge
		{
			pn_ed[i][0] = pn_ed0[i][0]; pn_ed[i][1] = pn_ed0[i][1]; pn_ed[i][2] = pn_ed0[i][2];
			ed_type[i] = 1;
			ed_shp[i] = 1;
		}
		else if (hedge[lev][edid[i]].type==1)//boundary non-sharp edge
		{
			array<double, 3> f_avg = { 0., 0., 0. };
			int nfc(0);
			for (j = 0; j < hedge[lev][edid[i]].face.size(); j++)
			{
				if (hface[lev][hedge[lev][edid[i]].face[j]].type == 1)
				{
					vector<int>::iterator it = find(fcid.begin(), fcid.end(), hedge[lev][edid[i]].face[j]);
					int pos(it - fcid.begin());
					f_avg[0] += pn_fc0[pos][0]; f_avg[1] += pn_fc0[pos][1]; f_avg[2] += pn_fc0[pos][2];
					nfc++;
				}
			}
			if (nfc != 2)
			{
				cerr << "# of boundary faces connected to the edge is wrong!\n";
				getchar();
			}
			f_avg[0] /= nfc; f_avg[1] /= nfc; f_avg[2] /= nfc;
			pn_ed[i][0] = (f_avg[0] + pn_ed0[i][0]) / nfc;
			pn_ed[i][1] = (f_avg[1] + pn_ed0[i][1]) / nfc;
			pn_ed[i][2] = (f_avg[2] + pn_ed0[i][2]) / nfc;
			ed_type[i] = 1;
			ed_shp[i] = 0;
		}
		else//interior
		{
			array<double, 3> f_avg = { 0., 0., 0. }, b_avg = { 0., 0., 0. };
			for (j = 0; j < hedge[lev][edid[i]].hex.size(); j++)
			{
				vector<int>::iterator it = find(bdid.begin(), bdid.end(), hedge[lev][edid[i]].hex[j]);
				int pos(it - bdid.begin());
				b_avg[0] += pn_bd[pos][0]; b_avg[1] += pn_bd[pos][1]; b_avg[2] += pn_bd[pos][2];
			}
			b_avg[0] /= hedge[lev][edid[i]].hex.size(); b_avg[1] /= hedge[lev][edid[i]].hex.size(); b_avg[2] /= hedge[lev][edid[i]].hex.size();
			for (j = 0; j < hedge[lev][edid[i]].face.size(); j++)
			{
				vector<int>::iterator it = find(fcid.begin(), fcid.end(), hedge[lev][edid[i]].face[j]);
				int pos(it - fcid.begin());
				f_avg[0] += pn_fc0[pos][0]; f_avg[1] += pn_fc0[pos][1]; f_avg[2] += pn_fc0[pos][2];
			}
			f_avg[0] /= hedge[lev][edid[i]].face.size(); f_avg[1] /= hedge[lev][edid[i]].face.size(); f_avg[2] /= hedge[lev][edid[i]].face.size();
			uint n = hedge[lev][edid[i]].face.size();
			pn_ed[i][0] = (b_avg[0] + 2.*f_avg[0] + (n - 3)*pn_ed0[i][0]) / n;
			pn_ed[i][1] = (b_avg[1] + 2.*f_avg[1] + (n - 3)*pn_ed0[i][1]) / n;
			pn_ed[i][2] = (b_avg[2] + 2.*f_avg[2] + (n - 3)*pn_ed0[i][2]) / n;
		}
	}
	for (i = 0; i < 8; i++)
	{
		pn_vt[i][0] = hcp[lev][hmesh[lev][eid].cnct[i]].coor[0];
		pn_vt[i][1] = hcp[lev][hmesh[lev][eid].cnct[i]].coor[1];
		pn_vt[i][2] = hcp[lev][hmesh[lev][eid].cnct[i]].coor[2];
		if (hcp[lev][hmesh[lev][eid].cnct[i]].sharp == 2)//boundary sharp corner
		{
			vt_type[i] = 1;
			vt_shp[i] = 2;
		}
		else if (hcp[lev][hmesh[lev][eid].cnct[i]].sharp == 1)//boundary sharp edge point
		{
			array<double, 3> e_avg = { 0., 0., 0. };
			int ned(0);
			for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].edge.size(); j++)
			{
				if (hedge[lev][hcp[lev][hmesh[lev][eid].cnct[i]].edge[j]].sharp == 1)
				{
					ned++;
					vector<int>::iterator it = find(edid.begin(), edid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].edge[j]);
					int pos(it - edid.begin());
					e_avg[0] += pn_ed0[pos][0]; e_avg[1] += pn_ed0[pos][1]; e_avg[2] += pn_ed0[pos][2];
				}
			}
			if (ned != 2)
			{
				cerr << "# of sharp edges connected to the point is wrong!\n";
				cout << ned << "\n";
				getchar();
			}
			e_avg[0] /= ned; e_avg[1] /= ned; e_avg[2] /= ned;
			pn_vt[i][0] = (e_avg[0] + pn_vt[i][0]) / ned;
			pn_vt[i][1] = (e_avg[1] + pn_vt[i][1]) / ned;
			pn_vt[i][2] = (e_avg[2] + pn_vt[i][2]) / ned;
			vt_type[i] = 1;
			vt_shp[i] = 1;
		}
		else if (hcp[lev][hmesh[lev][eid].cnct[i]].sharp == 0 && hcp[lev][hmesh[lev][eid].cnct[i]].type == 1)//boundary non-sharp
		{
			array<double, 3> f_avg = { 0., 0., 0. }, e_avg = { 0., 0., 0. };
			int n(0);
			for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].face.size(); j++)
			{
				if (hface[lev][hcp[lev][hmesh[lev][eid].cnct[i]].face[j]].type == 1)
				{
					n++;
					vector<int>::iterator it = find(fcid.begin(), fcid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].face[j]);
					int pos(it - fcid.begin());
					f_avg[0] += pn_fc0[pos][0]; f_avg[1] += pn_fc0[pos][1]; f_avg[2] += pn_fc0[pos][2];
				}
			}
			f_avg[0] /= n; f_avg[1] /= n; f_avg[2] /= n;
			n = 0;
			for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].edge.size(); j++)
			{
				if (hedge[lev][hcp[lev][hmesh[lev][eid].cnct[i]].edge[j]].type == 1)
				{
					n++;
					vector<int>::iterator it = find(edid.begin(), edid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].edge[j]);
					int pos(it - edid.begin());
					e_avg[0] += pn_ed0[pos][0]; e_avg[1] += pn_ed0[pos][1]; e_avg[2] += pn_ed0[pos][2];
				}
			}
			e_avg[0] /= n; e_avg[1] /= n; e_avg[2] /= n;
			pn_vt[i][0] = (f_avg[0] + 2.*e_avg[0] + (n-3)*pn_vt[i][0]) / n;
			pn_vt[i][1] = (f_avg[1] + 2.*e_avg[1] + (n - 3)*pn_vt[i][1]) / n;
			pn_vt[i][2] = (f_avg[2] + 2.*e_avg[2] + (n - 3)*pn_vt[i][2]) / n;
			vt_type[i] = 1;
			vt_shp[i] = 0;
		}
		else if (hcp[lev][hmesh[lev][eid].cnct[i]].type != 1)//interior
		{
			array<double, 3> b_avg = { 0., 0., 0. }, f_avg = { 0., 0., 0. }, e_avg = { 0., 0., 0. };
			uint n(hcp[lev][hmesh[lev][eid].cnct[i]].hex.size());
			for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].hex.size(); j++)
			{
				vector<int>::iterator it = find(bdid.begin(), bdid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].hex[j]);
				int pos(it - bdid.begin());
				b_avg[0] += pn_bd[pos][0]; b_avg[1] += pn_bd[pos][1]; b_avg[2] += pn_bd[pos][2];
			}
			b_avg[0] /= n; b_avg[1] /= n; b_avg[2] /= n;
			n = hcp[lev][hmesh[lev][eid].cnct[i]].face.size();
			for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].face.size(); j++)
			{
				vector<int>::iterator it = find(fcid.begin(), fcid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].face[j]);
				int pos(it - fcid.begin());
				f_avg[0] += pn_fc0[pos][0]; f_avg[1] += pn_fc0[pos][1]; f_avg[2] += pn_fc0[pos][2];
			}
			f_avg[0] /= n; f_avg[1] /= n; f_avg[2] /= n;
			n = hcp[lev][hmesh[lev][eid].cnct[i]].edge.size();
			for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].edge.size(); j++)
			{
				vector<int>::iterator it = find(edid.begin(), edid.end(), hcp[lev][hmesh[lev][eid].cnct[i]].edge[j]);
				int pos(it - edid.begin());
				e_avg[0] += pn_ed0[pos][0]; e_avg[1] += pn_ed0[pos][1]; e_avg[2] += pn_ed0[pos][2];
			}
			e_avg[0] /= n; e_avg[1] /= n; e_avg[2] /= n;
			pn_vt[i][0] = (b_avg[0] + 3.*f_avg[0] + 3.*e_avg[0] + pn_vt[i][0]) / 8.;
			pn_vt[i][1] = (b_avg[1] + 3.*f_avg[1] + 3.*e_avg[1] + pn_vt[i][1]) / 8.;
			pn_vt[i][2] = (b_avg[2] + 3.*f_avg[2] + 3.*e_avg[2] + pn_vt[i][2]) / 8.;
		}
	}
	//find global index of new points at Level+1
	vector<int> pnid_bd(pn_bd.size()), pnid_fc(pn_fc.size()), pnid_ed(pn_ed.size()), pnid_vt(pn_vt.size());
	double tol(1.e-6), dis(0.);
	uint np0(hcp[lev + 1].size());
	for (i = 0; i < pn_bd.size(); i++)
	{
		int flag(0);
		for (j = 0; j < np0; j++)
		{
			dis = sqrt((pn_bd[i][0] - hcp[lev + 1][j].coor[0])*(pn_bd[i][0] - hcp[lev + 1][j].coor[0]) + (pn_bd[i][1] - hcp[lev + 1][j].coor[1])*(pn_bd[i][1] - hcp[lev + 1][j].coor[1]) +
				(pn_bd[i][2] - hcp[lev + 1][j].coor[2])*(pn_bd[i][2] - hcp[lev + 1][j].coor[2]));
			if (dis < tol)
			{
				flag = 1; pnid_bd[i] = j; break;
			}
		}
		if (flag == 0)
		{
			pnid_bd[i] = hcp[lev + 1].size();
			Vertex3D cptmp;
			cptmp.coor[0] = pn_bd[i][0]; cptmp.coor[1] = pn_bd[i][1]; cptmp.coor[2] = pn_bd[i][2];
			cptmp.lev = lev + 1;
			hcp[lev + 1].push_back(cptmp);
		}
	}
	for (i = 0; i < pn_fc.size(); i++)
	{
		int flag(0);
		for (j = 0; j < np0; j++)
		{
			dis = sqrt((pn_fc[i][0] - hcp[lev + 1][j].coor[0])*(pn_fc[i][0] - hcp[lev + 1][j].coor[0]) + (pn_fc[i][1] - hcp[lev + 1][j].coor[1])*(pn_fc[i][1] - hcp[lev + 1][j].coor[1]) +
				(pn_fc[i][2] - hcp[lev + 1][j].coor[2])*(pn_fc[i][2] - hcp[lev + 1][j].coor[2]));
			if (dis < tol)
			{
				flag = 1; pnid_fc[i] = j; break;
			}
		}
		if (flag == 0)
		{
			pnid_fc[i] = hcp[lev + 1].size();
			Vertex3D cptmp;
			cptmp.coor[0] = pn_fc[i][0]; cptmp.coor[1] = pn_fc[i][1]; cptmp.coor[2] = pn_fc[i][2];
			cptmp.lev = lev + 1;
			cptmp.type = fc_type[i];
			hcp[lev + 1].push_back(cptmp);
		}
	}
	for (i = 0; i < pn_ed.size(); i++)
	{
		int flag(0);
		for (j = 0; j < np0; j++)
		{
			dis = sqrt((pn_ed[i][0] - hcp[lev + 1][j].coor[0])*(pn_ed[i][0] - hcp[lev + 1][j].coor[0]) + (pn_ed[i][1] - hcp[lev + 1][j].coor[1])*(pn_ed[i][1] - hcp[lev + 1][j].coor[1]) +
				(pn_ed[i][2] - hcp[lev + 1][j].coor[2])*(pn_ed[i][2] - hcp[lev + 1][j].coor[2]));
			if (dis < tol)
			{
				flag = 1; pnid_ed[i] = j; break;
			}
		}
		if (flag == 0)
		{
			pnid_ed[i] = hcp[lev + 1].size();
			Vertex3D cptmp;
			cptmp.coor[0] = pn_ed[i][0]; cptmp.coor[1] = pn_ed[i][1]; cptmp.coor[2] = pn_ed[i][2];
			cptmp.lev = lev + 1;
			cptmp.type = ed_type[i];
			cptmp.sharp = ed_shp[i];
			hcp[lev + 1].push_back(cptmp);
		}
	}
	for (i = 0; i < pn_vt.size(); i++)
	{
		int flag(0);
		for (j = 0; j < np0; j++)
		{
			dis = sqrt((pn_vt[i][0] - hcp[lev + 1][j].coor[0])*(pn_vt[i][0] - hcp[lev + 1][j].coor[0]) + (pn_vt[i][1] - hcp[lev + 1][j].coor[1])*(pn_vt[i][1] - hcp[lev + 1][j].coor[1]) +
				(pn_vt[i][2] - hcp[lev + 1][j].coor[2])*(pn_vt[i][2] - hcp[lev + 1][j].coor[2]));
			if (dis < tol)
			{
				flag = 1; pnid_vt[i] = j; break;
			}
		}
		if (flag == 0)
		{
			pnid_vt[i] = hcp[lev + 1].size();
			Vertex3D cptmp;
			cptmp.coor[0] = pn_vt[i][0]; cptmp.coor[1] = pn_vt[i][1]; cptmp.coor[2] = pn_vt[i][2];
			cptmp.lev = lev + 1;
			cptmp.type = vt_type[i];
			cptmp.sharp = vt_shp[i];
			hcp[lev + 1].push_back(cptmp);
		}
	}
	//find children id and corresponding coefficients
	//vector<vector<int>> chdid(hmesh[lev][eid].IEN.size());
	//vector<vector<double>> coef(hmesh[lev][eid].IEN.size());
	vector<array<int, 8>> bdmp(bdid.size());
	vector<array<double, 8>> bdc(bdid.size());
	vector<vector<int>> fcmp(fcid.size());
	vector<vector<double>> fcc(fcid.size());
	vector<vector<int>> edmp(edid.size());
	vector<vector<double>> edc(edid.size());
	vector<vector<int>> vtmp(8);
	vector<vector<double>> vtc(8);
	//body points
	for (i = 0; i < bdid.size(); i++)
	{
		for (j = 0; j < 8; j++)
		{
			bdmp[i][j] = hmesh[lev][bdid[i]].cnct[j];
			bdc[i][j] = .125;
		}
	}
	//face points
	for (i = 0; i < fcid.size(); i++)
	{
		if (hface[lev][fcid[i]].type == 0)
		{
			for (j = 0; j < 4; j++)
			{
				fcmp[i].push_back(hface[lev][fcid[i]].cnct[j]);
				fcc[i].push_back(2.*.25 / 4.);
			}
			for (j = 0; j < hface[lev][fcid[i]].hex.size(); j++)
			{
				int hxid(hface[lev][fcid[i]].hex[j]);
				for (k = 0; k < 8; k++)
				{
					vector<int>::iterator it = find(fcmp[i].begin(), fcmp[i].end(), hmesh[lev][hxid].cnct[k]);
					double ctmp(.125 / 4.);
					if (it == fcmp[i].end())
					{
						fcmp[i].push_back(hmesh[lev][hxid].cnct[k]);
						fcc[i].push_back(ctmp);
					}
					else
					{
						int pos(it - fcmp[i].begin());
						fcc[i][pos] += ctmp;
					}
				}
			}
		}
		else if (hface[lev][fcid[i]].type == 1)
		{
			for (j = 0; j < 4; j++)
			{
				fcmp[i].push_back(hface[lev][fcid[i]].cnct[j]);
				fcc[i].push_back(.25);
			}
		}
	}
	//cout
	/*for (i = 0; i < fcc.size(); i++)
	{
		double sum(0.);
		for (j = 0; j < fcc[i].size(); j++)
		{
			sum += fcc[i][j];
		}
		if (sum != 1.)
		{
			cout << "face point sum not 1\n";
			getchar();
		}
	}*/
	//edge points
	for (i = 0; i < edid.size(); i++)
	{
		if (hedge[lev][edid[i]].sharp == 1)
		{
			for (j = 0; j < 2; j++)
			{
				edmp[i].push_back(hedge[lev][edid[i]].pt[j]);
				edc[i].push_back(.5);
			}
		}
		else if (hedge[lev][edid[i]].type == 1)
		{
			for (j = 0; j < 2; j++)
			{
				edmp[i].push_back(hedge[lev][edid[i]].pt[j]);
				edc[i].push_back(.375);
			}
			for (j = 0; j < hedge[lev][edid[i]].face.size(); j++)
			{
				if (hface[lev][hedge[lev][edid[i]].face[j]].type == 1)
				{
					int fcid(hedge[lev][edid[i]].face[j]);
					for (k = 0; k < 4; k++)
					{
						vector<int>::iterator it = find(edmp[i].begin(), edmp[i].end(), hface[lev][fcid].cnct[k]);
						if (it == edmp[i].end())
						{
							edmp[i].push_back(hface[lev][fcid].cnct[k]);
							edc[i].push_back(.0625);
						}
					}
				}
			}
		}
		else
		{
			uint n(hedge[lev][edid[i]].hex.size());
			for (j = 0; j < 2; j++)
			{
				edmp[i].push_back(hedge[lev][edid[i]].pt[j]);
				edc[i].push_back(.5*double(n - 3) / n);
			}
			for (j = 0; j < hedge[lev][edid[i]].hex.size(); j++)
			{
				int hxid(hedge[lev][edid[i]].hex[j]);
				for (k = 0; k < 8; k++)
				{
					vector<int>::iterator it = find(edmp[i].begin(), edmp[i].end(), hmesh[lev][hxid].cnct[k]);
					double ctmp(.125 / double(n*n));
					if (it == edmp[i].end())
					{
						edmp[i].push_back(hmesh[lev][hxid].cnct[k]);
						edc[i].push_back(ctmp);
					}
					else
					{
						int pos(it - edmp[i].begin());
						edc[i][pos] += ctmp;
					}
				}
			}
			for (j = 0; j < hedge[lev][edid[i]].face.size(); j++)
			{
				int fcid(hedge[lev][edid[i]].face[j]);
				for (k = 0; k < 4; k++)
				{
					vector<int>::iterator it = find(edmp[i].begin(), edmp[i].end(), hface[lev][fcid].cnct[k]);
					double ctmp(2.*.25 / double(n*n));
					if (it == edmp[i].end())
					{
						edmp[i].push_back(hface[lev][fcid].cnct[k]);
						edc[i].push_back(ctmp);
					}
					else
					{
						int pos(it - edmp[i].begin());
						edc[i][pos] += ctmp;
					}
				}
			}
		}
	}
	//cout
	//for (i = 0; i < edc.size(); i++)
	//{
	//	double sum(0.), tol(1.e-10);
	//	for (j = 0; j < edc[i].size(); j++)
	//	{
	//		sum += edc[i][j];
	//	}
	//	if (fabs(sum-1.)>tol)
	//	{
	//		cout << "edge point " << sum <<"\n";
	//		//cout << "edge valence hex: " << edvl_hx[i] << "\n";
	//		////cout << "edge valence face: " << edvl_fc[i] << "\n";
	//		//cout << "# contribution: " << edmp[i].size() << "\n";
	//		//for (j = 0; j < edmp[i].size(); j++)
	//		//{
	//		//	cout << setw(6)<< edmp[i][j] << " ";
	//		//}
	//		//cout << "\n";
	//		//for (j = 0; j < edc[i].size(); j++)
	//		//{
	//		//	cout << setw(6) << edc[i][j] << " ";
	//		//}
	//		//cout << "\n";
	//		getchar();
	//	}
	//}
	//vertex points
	for (i = 0; i < 8; i++)
	{
		if (hcp[lev][hmesh[lev][eid].cnct[i]].sharp==2)
		{
			vtmp[i].push_back(hmesh[lev][eid].cnct[i]);
			vtc[i].push_back(1.);
		}
		else if (hcp[lev][hmesh[lev][eid].cnct[i]].sharp == 1)
		{
			int ptid(hmesh[lev][eid].cnct[i]);
			vtmp[i].push_back(ptid);
			vtc[i].push_back(.75);
			for (j = 0; j < hcp[lev][ptid].edge.size(); j++)
			{
				if (hedge[lev][hcp[lev][ptid].edge[j]].sharp == 1)
				{
					int ptmp = hedge[lev][hcp[lev][ptid].edge[j]].pt[0];
					if (ptmp == ptid) ptmp = hedge[lev][hcp[lev][ptid].edge[j]].pt[1];
					vtmp[i].push_back(ptmp);
					vtc[i].push_back(.125);
				}
			}
		}
		else if (hcp[lev][hmesh[lev][eid].cnct[i]].type == 1)
		{
			int ptid(hmesh[lev][eid].cnct[i]), nv(0);
			for (j = 0; j < hcp[lev][ptid].face.size(); j++)
			{
				if (hface[lev][hcp[lev][ptid].face[j]].type == 1) nv++;
			}
			vtmp[i].push_back(ptid);
			vtc[i].push_back(1.-7./(4.*nv));
			double ctmp[2] = { 3. / (2.*nv*nv), 1. / (4.*nv*nv) };
			for (j = 0; j < hcp[lev][ptid].edge.size(); j++)
			{
				if (hedge[lev][hcp[lev][ptid].edge[j]].type == 1)
				{
					int edid(hcp[lev][ptid].edge[j]);
					int ptmp(hedge[lev][edid].pt[0]);
					if (ptmp == ptid) ptmp = hedge[lev][edid].pt[1];
					vtmp[i].push_back(ptmp);
					vtc[i].push_back(ctmp[0]);
				}
			}
			for (j = 0; j < hcp[lev][ptid].face.size(); j++)
			{
				if (hface[lev][hcp[lev][ptid].face[j]].type == 1)
				{
					int fcid(hcp[lev][ptid].face[j]);
					for (k = 0; k < 4; k++)
					{
						vector<int>::iterator it = find(vtmp[i].begin(), vtmp[i].end(), hface[lev][fcid].cnct[k]);
						if (it == vtmp[i].end())
						{
							vtmp[i].push_back(hface[lev][fcid].cnct[k]);
							vtc[i].push_back(ctmp[1]);
						}
					}
				}
			}
		}
		else
		{
			int ptid(hmesh[lev][eid].cnct[i]);
			uint nhx(hcp[lev][ptid].hex.size()), nfc(hcp[lev][ptid].face.size()), ned(hcp[lev][ptid].edge.size());
			vtmp[i].push_back(ptid);
			vtc[i].push_back(.125);
			for (j = 0; j < nhx; j++)
			{
				int hxid(hcp[lev][ptid].hex[j]);
				for (k = 0; k < 8; k++)
				{
					vector<int>::iterator it = find(vtmp[i].begin(), vtmp[i].end(), hmesh[lev][hxid].cnct[k]);
					double ctmp(.125 / (8.*nhx));
					if (it == vtmp[i].end())
					{
						vtmp[i].push_back(hmesh[lev][hxid].cnct[k]);
						vtc[i].push_back(ctmp);
					}
					else
					{
						int pos(it - vtmp[i].begin());
						vtc[i][pos] += ctmp;
					}
				}
			}
			for (j = 0; j < nfc; j++)
			{
				int fcid(hcp[lev][ptid].face[j]);
				for (k = 0; k < 4; k++)
				{
					vector<int>::iterator it = find(vtmp[i].begin(), vtmp[i].end(), hface[lev][fcid].cnct[k]);
					double ctmp(3.*.25 / (8.*nfc));
					if (it == vtmp[i].end())
					{
						vtmp[i].push_back(hface[lev][fcid].cnct[k]);
						vtc[i].push_back(ctmp);
					}
					else
					{
						int pos(it - vtmp[i].begin());
						vtc[i][pos] += ctmp;
					}
				}
			}
			for (j = 0; j < ned; j++)
			{
				int edid(hcp[lev][ptid].edge[j]);
				for (k = 0; k < 2; k++)
				{
					vector<int>::iterator it = find(vtmp[i].begin(), vtmp[i].end(), hedge[lev][edid].pt[k]);
					double ctmp(3.*.5 / (8.*ned));
					if (it == vtmp[i].end())
					{
						vtmp[i].push_back(hedge[lev][edid].pt[k]);
						vtc[i].push_back(ctmp);
					}
					else
					{
						int pos(it - vtmp[i].begin());
						vtc[i][pos] += ctmp;
					}
				}
			}
		}
	}
	//cout
	/*for (i = 0; i < vtc.size(); i++)
	{
		double sum(0.);
		for (j = 0; j < vtc[i].size(); j++)
		{
			sum += vtc[i][j];
		}
		if (fabs(sum - 1.)>tol)
		{
			cout << "vertex point " << sum << "\n";
			cout << hmesh[lev][eid].cnct[i]<< " vertex type " << hcp[lev][hmesh[lev][eid].cnct[i]].type << "\n";
			for (j = 0; j < vtmp[i].size(); j++)
			{
				cout << vtmp[i][j] << " ";
			}
			cout << "\n";
			for (j = 0; j < vtc[i].size(); j++)
			{
				cout<< vtc[i][j]<<" ";
			}
			cout << "\n";
			getchar();
		}
	}*/
	//children and coef
	for (i = 0; i < bdmp.size(); i++)
	{
		for (j = 0; j < bdmp[i].size(); j++)
		{
			vector<int>::iterator it = find(hcp[lev][bdmp[i][j]].chd.begin(), hcp[lev][bdmp[i][j]].chd.end(), pnid_bd[i]);
			if (it == hcp[lev][bdmp[i][j]].chd.end())
			{
				hcp[lev][bdmp[i][j]].chd.push_back(pnid_bd[i]);
				hcp[lev][bdmp[i][j]].coef.push_back(bdc[i][j]);
			}
		}
	}
	for (i = 0; i < fcmp.size(); i++)
	{
		for (j = 0; j < fcmp[i].size(); j++)
		{
			vector<int>::iterator it = find(hcp[lev][fcmp[i][j]].chd.begin(), hcp[lev][fcmp[i][j]].chd.end(), pnid_fc[i]);
			if (it == hcp[lev][fcmp[i][j]].chd.end())
			{
				hcp[lev][fcmp[i][j]].chd.push_back(pnid_fc[i]);
				hcp[lev][fcmp[i][j]].coef.push_back(fcc[i][j]);
			}
		}
	}
	for (i = 0; i < edmp.size(); i++)
	{
		for (j = 0; j < edmp[i].size(); j++)
		{
			vector<int>::iterator it = find(hcp[lev][edmp[i][j]].chd.begin(), hcp[lev][edmp[i][j]].chd.end(), pnid_ed[i]);
			if (it == hcp[lev][edmp[i][j]].chd.end())
			{
				hcp[lev][edmp[i][j]].chd.push_back(pnid_ed[i]);
				hcp[lev][edmp[i][j]].coef.push_back(edc[i][j]);
			}
		}
	}
	for (i = 0; i < vtmp.size(); i++)
	{
		for (j = 0; j < vtmp[i].size(); j++)
		{
			vector<int>::iterator it = find(hcp[lev][vtmp[i][j]].chd.begin(), hcp[lev][vtmp[i][j]].chd.end(), pnid_vt[i]);
			if (it == hcp[lev][vtmp[i][j]].chd.end())
			{
				hcp[lev][vtmp[i][j]].chd.push_back(pnid_vt[i]);
				hcp[lev][vtmp[i][j]].coef.push_back(vtc[i][j]);
			}
		}
	}
	//local index of new points
	vector<int> ed_loc(12), fc_loc(6), vt_loc(pnid_vt);
	int bd_loc(pnid_bd[0]);
	for (i = 0; i < 12; i++)
	{
		vector<int>::iterator it = find(edid.begin(), edid.end(), hmesh[lev][eid].edge[i]);
		int pos(it - edid.begin());
		ed_loc[i] = pnid_ed[pos];
	}
	for (i = 0; i < 6; i++)
	{
		vector<int>::iterator it = find(fcid.begin(), fcid.end(), hmesh[lev][eid].face[i]);
		int pos(it - fcid.begin());
		fc_loc[i] = pnid_fc[pos];
	}
	//construct new hex
	int enew[8][8] = { { vt_loc[0], ed_loc[0], fc_loc[0], ed_loc[3], ed_loc[4], fc_loc[1], bd_loc, fc_loc[4] },
	{ ed_loc[0], vt_loc[1], ed_loc[1], fc_loc[0], fc_loc[1], ed_loc[5], fc_loc[2], bd_loc },
	{ fc_loc[0], ed_loc[1], vt_loc[2], ed_loc[2], bd_loc, fc_loc[2], ed_loc[6], fc_loc[3] },
	{ ed_loc[3], fc_loc[0], ed_loc[2], vt_loc[3], fc_loc[4], bd_loc, fc_loc[3], ed_loc[7] },
	{ ed_loc[4], fc_loc[1], bd_loc, fc_loc[4], vt_loc[4], ed_loc[8], fc_loc[5], ed_loc[11] },
	{ fc_loc[1], ed_loc[5], fc_loc[2], bd_loc, ed_loc[8], vt_loc[5], ed_loc[9], fc_loc[5] },
	{ bd_loc, fc_loc[2], ed_loc[6], fc_loc[3], fc_loc[5], ed_loc[9], vt_loc[6], ed_loc[10] },
	{ fc_loc[4], bd_loc, fc_loc[3], ed_loc[7], ed_loc[11], fc_loc[5], ed_loc[10], vt_loc[7] } };

	for (i = 0; i < 8; i++)
	{
		Element3D etmp;
		etmp.type = 0;//tmp, determined by the topology at level (lev+1)
		etmp.lev = lev + 1;
		etmp.prt = eid;
		for (j = 0; j < 8; j++)
		{
			etmp.cnct[j] = enew[i][j];
		}
		hmesh[lev][eid].chd.push_back(hmesh[lev + 1].size());
		hmesh[lev + 1].push_back(etmp);
	}

	hmesh[lev][eid].act = 0;
}

void TruncatedTspline_3D::ConstructBezierBasis(int lev, int eid)
{
	//find IENtmp first
	hmesh[lev][eid].IEN.clear();
	uint i, j, k, hxid;
	vector<int> loc(hcp[lev].size(), -1);
	vector<int> hx1r(1, eid);
	for (i = 0; i < 8; i++)
	{
		loc[hmesh[lev][eid].cnct[i]] = hmesh[lev][eid].IEN.size();
		hmesh[lev][eid].IEN.push_back(hmesh[lev][eid].cnct[i]);
	}
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].hex.size(); j++)
		{
			hxid = hcp[lev][hmesh[lev][eid].cnct[i]].hex[j];
			vector<int>::iterator it1 = find(hx1r.begin(), hx1r.end(), hxid);
			if (it1 == hx1r.end())
			{
				hx1r.push_back(hxid);
				for (k = 0; k < 8; k++)
				{
					vector<int>::iterator it = find(hmesh[lev][eid].IEN.begin(), hmesh[lev][eid].IEN.end(), hmesh[lev][hxid].cnct[k]);
					if (it == hmesh[lev][eid].IEN.end())
					{
						loc[hmesh[lev][hxid].cnct[k]] = hmesh[lev][eid].IEN.size();
						hmesh[lev][eid].IEN.push_back(hmesh[lev][hxid].cnct[k]);
					}
				}
			}
		}
	}
	for (i = 0; i < hmesh[lev][eid].bemat.size(); i++)
	{
		hmesh[lev][eid].bemat[i].clear();
	}
	hmesh[lev][eid].bemat.clear();
	hmesh[lev][eid].bemat.resize(hmesh[lev][eid].IEN.size(), vector<double>(64, 0.));
	//8 body points, not consider boundary yet
	double w[2] = { 2. / 3., 1. / 3. };
	double a[8] = { w[0] * w[0] * w[0], w[1] * w[0] * w[0], w[1] * w[1] * w[0], w[1] * w[0] * w[0],
		w[1] * w[0] * w[0], w[1] * w[1] * w[0], w[1] * w[1] * w[1], w[1] * w[1] * w[0] };
	int bpi[8][8] = { { 0, 1, 2, 3, 4, 5, 6, 7 }, { 1, 2, 3, 0, 5, 6, 7, 4 }, { 2, 3, 0, 1, 6, 7, 4, 5 }, { 3, 0, 1, 2, 7, 4, 5, 6 },
	{ 4, 5, 6, 7, 0, 1, 2, 3 }, { 5, 6, 7, 4, 1, 2, 3, 0 }, { 6, 7, 4, 5, 2, 3, 0, 1 }, { 7, 4, 5, 6, 3, 0, 1, 2 } };
	vector<array<array<double, 8>, 8>> bpm(hx1r.size());
	vector<array<array<int, 8>, 8>> bpmap(hx1r.size());
	for (i = 0; i < hx1r.size(); i++)//which element
	{
		for (j = 0; j < 8; j++)//which body point, bezier
		{
			for (k = 0; k < 8; k++)//which local corner point, b-splines
			{
				bpm[i][j][k] = a[k];
				bpmap[i][j][k] = loc[hmesh[lev][hx1r[i]].cnct[bpi[j][k]]];
			}
		}
	}
	int layer[4] = { 0, 16, 32, 48 };
	int bpbz[8] = { 5 + layer[1], 6 + layer[1], 10 + layer[1], 9 + layer[1], 5 + layer[2], 6 + layer[2], 10 + layer[2], 9 + layer[2] };
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 8; j++)
		{
			hmesh[lev][eid].bemat[bpmap[0][i][j]][bpbz[i]] = bpm[0][i][j];
		}
	}
	//2*12 edge points
	int edi[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 } };
	int edbz[12][2] = { { 1, 2 }, { 7, 11 }, { 14, 13 }, { 8, 4 }, { 0 + layer[1], 0 + layer[2] }, { 3 + layer[1], 3 + layer[2] },
	{ 15 + layer[1], 15 + layer[2] }, { 12 + layer[1], 12 + layer[2] }, { 1 + layer[3], 2 + layer[3] }, { 7 + layer[3], 11 + layer[3] }, { 14 + layer[3], 13 + layer[3] }, { 8 + layer[3], 4 + layer[3] } };
	int pos1, pos2;
	for (i = 0; i < 12; i++)
	{
		uint nhex = hedge[lev][hmesh[lev][eid].edge[i]].hex.size();
		for (j = 0; j<hedge[lev][hmesh[lev][eid].edge[i]].hex.size(); j++)
		{
			hxid = hedge[lev][hmesh[lev][eid].edge[i]].hex[j];
			vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
			pos1 = it - hx1r.begin();
			int* it1 = find(hmesh[lev][hxid].cnct, hmesh[lev][hxid].cnct + 8, hmesh[lev][eid].cnct[edi[i][0]]);
			pos2 = it1 - hmesh[lev][hxid].cnct;
			for (k = 0; k < 8; k++)
			{
				hmesh[lev][eid].bemat[bpmap[pos1][pos2][k]][edbz[i][0]] += bpm[pos1][pos2][k] / nhex;
			}
			int* it2 = find(hmesh[lev][hxid].cnct, hmesh[lev][hxid].cnct + 8, hmesh[lev][eid].cnct[edi[i][1]]);
			pos2 = it2 - hmesh[lev][hxid].cnct;
			for (k = 0; k < 8; k++)
			{
				hmesh[lev][eid].bemat[bpmap[pos1][pos2][k]][edbz[i][1]] += bpm[pos1][pos2][k] / nhex;
			}
		}
	}
	//8 corner points
	int cnbz[8] = { 0, 3, 15, 12, 0 + layer[3], 3 + layer[3], 15 + layer[3], 12 + layer[3] };
	for (i = 0; i < 8; i++)
	{
		uint nhex = hcp[lev][hmesh[lev][eid].cnct[i]].hex.size();
		for (j = 0; j<nhex; j++)
		{
			hxid = hcp[lev][hmesh[lev][eid].cnct[i]].hex[j];
			vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
			pos1 = it - hx1r.begin();
			int* it1 = find(hmesh[lev][hxid].cnct, hmesh[lev][hxid].cnct + 8, hmesh[lev][eid].cnct[i]);
			pos2 = it1 - hmesh[lev][hxid].cnct;
			for (k = 0; k < 8; k++)
			{
				hmesh[lev][eid].bemat[bpmap[pos1][pos2][k]][cnbz[i]] += bpm[pos1][pos2][k] / nhex;
			}
		}
	}
	//4*6 face points
	int fci[6][4] = { { 0, 1, 3, 2 }, { 0, 1, 4, 5 }, { 1, 2, 5, 6 }, { 3, 2, 7, 6 }, { 0, 3, 4, 7 }, { 4, 5, 7, 6 } };
	int fcbz[6][4] = { { 5, 6, 9, 10 }, { 1 + layer[1], 2 + layer[1], 1 + layer[2], 2 + layer[2] }, { 7 + layer[1], 11 + layer[1], 7 + layer[2], 11 + layer[2] },
	{ 13 + layer[1], 14 + layer[1], 13 + layer[2], 14 + layer[2] }, { 4 + layer[1], 8 + layer[1], 4 + layer[2], 8 + layer[2] }, { 5 + layer[3], 6 + layer[3], 9 + layer[3], 10 + layer[3] } };
	for (i = 0; i < 6; i++)
	{
		uint nhex = hface[lev][hmesh[lev][eid].face[i]].hex.size();
		for (j = 0; j < nhex; j++)
		{
			hxid = hface[lev][hmesh[lev][eid].face[i]].hex[j];
			vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
			pos1 = it - hx1r.begin();
			for (int j1 = 0; j1 < 4; j1++)
			{
				int* it1 = find(hmesh[lev][hxid].cnct, hmesh[lev][hxid].cnct + 8, hmesh[lev][eid].cnct[fci[i][j1]]);
				pos2 = it1 - hmesh[lev][hxid].cnct;
				for (k = 0; k < 8; k++)
				{
					hmesh[lev][eid].bemat[bpmap[pos1][pos2][k]][fcbz[i][j1]] += bpm[pos1][pos2][k] / nhex;
				}
			}
		}
	}
}

void TruncatedTspline_3D::ConstructBezierBasis_1(int lev, int eid)//old idea to use open knot vector by pillowing layers, troublesome
{
	//find IEN first, element eid plus one-ring neighborhood
	hmesh[lev][eid].IEN.clear();
	uint i, j, k, hxid;
	vector<int> loc(hcp[lev].size(), -1);
	vector<int> hx1r(1, eid);
	for (i = 0; i < 8; i++)
	{
		loc[hmesh[lev][eid].cnct[i]] = hmesh[lev][eid].IEN.size();
		hmesh[lev][eid].IEN.push_back(hmesh[lev][eid].cnct[i]);
	}
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].hex.size(); j++)
		{
			hxid = hcp[lev][hmesh[lev][eid].cnct[i]].hex[j];
			vector<int>::iterator it1 = find(hx1r.begin(), hx1r.end(), hxid);
			if (it1 == hx1r.end())
			{
				hx1r.push_back(hxid);
				for (k = 0; k < 8; k++)
				{
					vector<int>::iterator it = find(hmesh[lev][eid].IEN.begin(), hmesh[lev][eid].IEN.end(), hmesh[lev][hxid].cnct[k]);
					if (it == hmesh[lev][eid].IEN.end())
					{
						loc[hmesh[lev][hxid].cnct[k]] = hmesh[lev][eid].IEN.size();
						hmesh[lev][eid].IEN.push_back(hmesh[lev][hxid].cnct[k]);
					}
				}
			}
		}
	}
	for (i = 0; i < hmesh[lev][eid].bemat.size(); i++)
	{
		hmesh[lev][eid].bemat[i].clear();
	}
	hmesh[lev][eid].bemat.clear();
	hmesh[lev][eid].bemat.resize(hmesh[lev][eid].IEN.size(), vector<double>(64, 0.));
	//8 body points
	//double w[2] = { 2. / 3., 1. / 3. };
	//double a[8] = { w[0] * w[0] * w[0], w[1] * w[0] * w[0], w[1] * w[1] * w[0], w[1] * w[0] * w[0],
	//	w[1] * w[0] * w[0], w[1] * w[1] * w[0], w[1] * w[1] * w[1], w[1] * w[1] * w[0] };
	vector<array<int, 6>> hx1r_type(hx1r.size());
	vector<array<double, 9>> lens(hx1r.size());
	int ed_loc[3] = { 0, 3, 4 };
	int fc_loc[3][2] = { { 4, 2 }, { 1, 3 }, { 0, 5 } };
	for (i = 0; i < hx1r.size(); i++)
	{
		for (j = 0; j < 9; j++) lens[i][j] = 0.;
		if (hmesh[lev][hx1r[i]].type != 1)
		{
			for (j = 0; j < 3; j++)
			{
				lens[i][3 * j] = hedge[lev][hmesh[lev][hx1r[i]].edge[ed_loc[j]]].len;
				for (k = 0; k < 2; k++)
				{
					int fnb = hface[lev][hmesh[lev][hx1r[i]].face[fc_loc[j][k]]].hex[0];
					if (fnb == hx1r[i]) fnb = hface[lev][hmesh[lev][hx1r[i]].face[fc_loc[j][k]]].hex[1];
					int* it = find(hmesh[lev][fnb].face, hmesh[lev][fnb].face + 6, hmesh[lev][hx1r[i]].face[fc_loc[j][k]]);
					int pos(it - hmesh[lev][fnb].face);
					lens[i][3 * j+k+1] = hedge[lev][hmesh[lev][fnb].edge[fc_ppd_ed[pos]]].len;
				}
			}
		}
		//for (j = 0; j < 3; j++)
		//{
		//	if (hedge[lev][hmesh[lev][hx1r[i]].edge[loc_tmp[j]]].len == 0.)
		//	{
		//		hx1r_type[i][2*j] = 1;
		//	}
		//}
		//if (hx1r_type[i][0] == 0 && hx1r_type[i][2] == 0 && hx1r_type[i][4] == 0)
		//{
		//	//check three directions for face neighbors
		//	//int fc_loc[3][2] = { { 2, 4 }, { 1, 3 }, {0,5} };
		//	for (j = 0; j < 3; j++)
		//	{
		//		int fnb[2] = { hface[lev][hmesh[lev][hx1r[i]].face[fc_loc[j][0]]].hex[0], hface[lev][hmesh[lev][hx1r[i]].face[fc_loc[j][1]]].hex[0]};
		//		if (hface[lev][hmesh[lev][hx1r[i]].face[fc_loc[j][0]]].hex.size() != 2 || hface[lev][hmesh[lev][hx1r[i]].face[fc_loc[j][1]]].hex.size() != 2)
		//		{
		//			cerr << "Edge interval config wrong!\n"; 
		//			getchar();
		//		}
		//		if (fnb[0] == hx1r[i]) fnb[0] = hface[lev][hmesh[lev][hx1r[i]].face[fc_loc[j][0]]].hex[1];
		//		if (fnb[1] == hx1r[i]) fnb[1] = hface[lev][hmesh[lev][hx1r[i]].face[fc_loc[j][1]]].hex[1];
		//		int* it0 = find(hmesh[lev][fnb[0]].face, hmesh[lev][fnb[0]].face + 6, hmesh[lev][hx1r[i]].face[fc_loc[j][0]]);
		//		int* it1 = find(hmesh[lev][fnb[1]].face, hmesh[lev][fnb[1]].face + 6, hmesh[lev][hx1r[i]].face[fc_loc[j][1]]);
		//		int pos[2] = { it0 - hmesh[lev][fnb[0]].face, it1 - hmesh[lev][fnb[1]].face };
		//		if (it0 == hmesh[lev][fnb[0]].face + 6 || it1 == hmesh[lev][fnb[1]].face + 6)
		//		{
		//			cerr << "Cannot find the face in a face neighbor!\n";
		//			getchar();
		//		}
		//		if (hedge[lev][hmesh[lev][fnb[0]].edge[fc_ppd_ed[pos[0]]]].len == 0.)
		//		{
		//			hx1r_type[i][2*j] == 2;
		//		}
		//		if (hedge[lev][hmesh[lev][fnb[1]].edge[fc_ppd_ed[pos[1]]]].len == 0.)
		//		{
		//			hx1r_type[i][2 * j+1] == 2;
		//		}
		//	}
		//}
	}

	int bpi[8][8] = { { 0, 1, 2, 3, 4, 5, 6, 7 }, { 1, 2, 3, 0, 5, 6, 7, 4 }, { 2, 3, 0, 1, 6, 7, 4, 5 }, { 3, 0, 1, 2, 7, 4, 5, 6 },
	{ 4, 5, 6, 7, 0, 1, 2, 3 }, { 5, 6, 7, 4, 1, 2, 3, 0 }, { 6, 7, 4, 5, 2, 3, 0, 1 }, { 7, 4, 5, 6, 3, 0, 1, 2 } };
	vector<array<array<double, 8>, 8>> bpm(hx1r.size());
	vector<array<array<int, 8>, 8>> bpmap(hx1r.size());
	for (i = 0; i < hx1r.size(); i++)//which element
	{
		if (hmesh[lev][hx1r[i]].type != 1)
		{
			double w[3][4];
			for (j = 0; j < 3; j++)
			{
				w[j][0] = (lens[i][3 * j] + lens[i][3 * j + 2]) / (lens[i][3 * j] + lens[i][3 * j + 1] + lens[i][3 * j + 2]);
				w[j][1] = 1. - w[j][0];
				w[j][2] = (lens[i][3 * j] + lens[i][3 * j + 1]) / (lens[i][3 * j] + lens[i][3 * j + 1] + lens[i][3 * j + 2]);
				w[j][3] = 1. - w[j][2];
			}
			double au[8][8] = { { w[0][0], w[0][1], w[0][1], w[0][0], w[0][0], w[0][1], w[0][1], w[0][0] }, 
			{ w[0][3], w[0][2], w[0][2], w[0][3], w[0][3], w[0][2], w[0][2], w[0][3] }, 
			{ w[0][3], w[0][2], w[0][2], w[0][3], w[0][3], w[0][2], w[0][2], w[0][3] },
			{ w[0][0], w[0][1], w[0][1], w[0][0], w[0][0], w[0][1], w[0][1], w[0][0] },
			{ w[0][0], w[0][1], w[0][1], w[0][0], w[0][0], w[0][1], w[0][1], w[0][0] },
			{ w[0][3], w[0][2], w[0][2], w[0][3], w[0][3], w[0][2], w[0][2], w[0][3] },
			{ w[0][3], w[0][2], w[0][2], w[0][3], w[0][3], w[0][2], w[0][2], w[0][3] },
			{ w[0][0], w[0][1], w[0][1], w[0][0], w[0][0], w[0][1], w[0][1], w[0][0] }};
			double av[8][8] = { { w[1][0], w[1][0], w[1][1], w[1][1], w[1][0], w[1][0], w[1][1], w[1][1] },
			{ w[1][0], w[1][0], w[1][1], w[1][1], w[1][0], w[1][0], w[1][1], w[1][1] },
			{ w[1][3], w[1][3], w[1][2], w[1][2], w[1][3], w[1][3], w[1][2], w[1][2] },
			{ w[1][3], w[1][3], w[1][2], w[1][2], w[1][3], w[1][3], w[1][2], w[1][2] },
			{ w[1][0], w[1][0], w[1][1], w[1][1], w[1][0], w[1][0], w[1][1], w[1][1] },
			{ w[1][0], w[1][0], w[1][1], w[1][1], w[1][0], w[1][0], w[1][1], w[1][1] },
			{ w[1][3], w[1][3], w[1][2], w[1][2], w[1][3], w[1][3], w[1][2], w[1][2] },
			{ w[1][3], w[1][3], w[1][2], w[1][2], w[1][3], w[1][3], w[1][2], w[1][2] }};
			double aw[8][8] = { { w[2][0], w[2][0], w[2][0], w[2][0], w[2][1], w[2][1], w[2][1], w[2][1] },
			{ w[2][0], w[2][0], w[2][0], w[2][0], w[2][1], w[2][1], w[2][1], w[2][1] },
			{ w[2][0], w[2][0], w[2][0], w[2][0], w[2][1], w[2][1], w[2][1], w[2][1] },
			{ w[2][0], w[2][0], w[2][0], w[2][0], w[2][1], w[2][1], w[2][1], w[2][1] },
			{ w[2][3], w[2][3], w[2][3], w[2][3], w[2][2], w[2][2], w[2][2], w[2][2] },
			{ w[2][3], w[2][3], w[2][3], w[2][3], w[2][2], w[2][2], w[2][2], w[2][2] },
			{ w[2][3], w[2][3], w[2][3], w[2][3], w[2][2], w[2][2], w[2][2], w[2][2] },
			{ w[2][3], w[2][3], w[2][3], w[2][3], w[2][2], w[2][2], w[2][2], w[2][2] }};
			for (j = 0; j < 8; j++)//which body point, bezier
			{
				for (k = 0; k < 8; k++)//which local corner point, b-splines
				{
					bpm[i][j][k] = au[j][k]*av[j][k]*aw[j][k];
					bpmap[i][j][k] = loc[hmesh[lev][hx1r[i]].cnct[bpi[j][k]]];
				}
			}
		}
		else
		{
			//zero-parametric-area

		}

		////double w[3][2] = { { 2. / 3., 1. / 3. }, { 2. / 3., 1. / 3. }, { 2. / 3., 1. / 3. } };
		//double a[8] = { w[0][0] * w[1][0] * w[2][0], w[0][1] * w[1][0] * w[2][0], w[0][1] * w[1][1] * w[2][0], w[0][1] * w[1][0] * w[2][0],
		//	w[0][1] * w[1][0] * w[2][0], w[0][1] * w[1][1] * w[2][0], w[0][1] * w[1][1] * w[2][1], w[0][1] * w[1][1] * w[2][0] };
		//for (j = 0; j < 8; j++)//which body point, bezier
		//{
		//	for (k = 0; k < 8; k++)//which local corner point, b-splines
		//	{
		//		//bpm[i][j][k] = a[k];
		//		bpmap[i][j][k] = loc[hmesh[lev][hx1r[i]].cnct[bpi[j][k]]];
		//	}
		//}
	}
	int layer[4] = { 0, 16, 32, 48 };
	int bpbz[8] = { 5 + layer[1], 6 + layer[1], 10 + layer[1], 9 + layer[1], 5 + layer[2], 6 + layer[2], 10 + layer[2], 9 + layer[2] };
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 8; j++)
		{
			hmesh[lev][eid].bemat[bpmap[0][i][j]][bpbz[i]] = bpm[0][i][j];
		}
	}
	//2*12 edge points
	int edi[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 } };
	int edbz[12][2] = { { 1, 2 }, { 7, 11 }, { 14, 13 }, { 8, 4 }, { 0 + layer[1], 0 + layer[2] }, { 3 + layer[1], 3 + layer[2] },
	{ 15 + layer[1], 15 + layer[2] }, { 12 + layer[1], 12 + layer[2] }, { 1 + layer[3], 2 + layer[3] }, { 7 + layer[3], 11 + layer[3] }, { 14 + layer[3], 13 + layer[3] }, { 8 + layer[3], 4 + layer[3] } };
	int pos1, pos2;
	for (i = 0; i < 12; i++)
	{
		uint nhex = hedge[lev][hmesh[lev][eid].edge[i]].hex.size();
		for (j = 0; j<hedge[lev][hmesh[lev][eid].edge[i]].hex.size(); j++)
		{
			hxid = hedge[lev][hmesh[lev][eid].edge[i]].hex[j];
			vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
			pos1 = it - hx1r.begin();
			int* it1 = find(hmesh[lev][hxid].cnct, hmesh[lev][hxid].cnct + 8, hmesh[lev][eid].cnct[edi[i][0]]);
			pos2 = it1 - hmesh[lev][hxid].cnct;
			for (k = 0; k < 8; k++)
			{
				hmesh[lev][eid].bemat[bpmap[pos1][pos2][k]][edbz[i][0]] += bpm[pos1][pos2][k] / nhex;
			}
			int* it2 = find(hmesh[lev][hxid].cnct, hmesh[lev][hxid].cnct + 8, hmesh[lev][eid].cnct[edi[i][1]]);
			pos2 = it2 - hmesh[lev][hxid].cnct;
			for (k = 0; k < 8; k++)
			{
				hmesh[lev][eid].bemat[bpmap[pos1][pos2][k]][edbz[i][1]] += bpm[pos1][pos2][k] / nhex;
			}
		}
	}
	//8 corner points
	int cnbz[8] = { 0, 3, 15, 12, 0 + layer[3], 3 + layer[3], 15 + layer[3], 12 + layer[3] };
	for (i = 0; i < 8; i++)
	{
		uint nhex = hcp[lev][hmesh[lev][eid].cnct[i]].hex.size();
		for (j = 0; j<nhex; j++)
		{
			hxid = hcp[lev][hmesh[lev][eid].cnct[i]].hex[j];
			vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
			pos1 = it - hx1r.begin();
			int* it1 = find(hmesh[lev][hxid].cnct, hmesh[lev][hxid].cnct + 8, hmesh[lev][eid].cnct[i]);
			pos2 = it1 - hmesh[lev][hxid].cnct;
			for (k = 0; k < 8; k++)
			{
				hmesh[lev][eid].bemat[bpmap[pos1][pos2][k]][cnbz[i]] += bpm[pos1][pos2][k] / nhex;
			}
		}
	}
	//4*6 face points
	int fci[6][4] = { { 0, 1, 3, 2 }, { 0, 1, 4, 5 }, { 1, 2, 5, 6 }, { 3, 2, 7, 6 }, { 0, 3, 4, 7 }, { 4, 5, 7, 6 } };
	int fcbz[6][4] = { { 5, 6, 9, 10 }, { 1 + layer[1], 2 + layer[1], 1 + layer[2], 2 + layer[2] }, { 7 + layer[1], 11 + layer[1], 7 + layer[2], 11 + layer[2] },
	{ 13 + layer[1], 14 + layer[1], 13 + layer[2], 14 + layer[2] }, { 4 + layer[1], 8 + layer[1], 4 + layer[2], 8 + layer[2] }, { 5 + layer[3], 6 + layer[3], 9 + layer[3], 10 + layer[3] } };
	for (i = 0; i < 6; i++)
	{
		uint nhex = hface[lev][hmesh[lev][eid].face[i]].hex.size();
		for (j = 0; j < nhex; j++)
		{
			hxid = hface[lev][hmesh[lev][eid].face[i]].hex[j];
			vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
			pos1 = it - hx1r.begin();
			for (int j1 = 0; j1 < 4; j1++)
			{
				int* it1 = find(hmesh[lev][hxid].cnct, hmesh[lev][hxid].cnct + 8, hmesh[lev][eid].cnct[fci[i][j1]]);
				pos2 = it1 - hmesh[lev][hxid].cnct;
				for (k = 0; k < 8; k++)
				{
					hmesh[lev][eid].bemat[bpmap[pos1][pos2][k]][fcbz[i][j1]] += bpm[pos1][pos2][k] / nhex;
				}
			}
		}
	}
}

void TruncatedTspline_3D::ConstructBezierBasis_Boundary(int lev, int eid)
{
	//find IEN first
	hmesh[lev][eid].IEN.clear();
	uint i, j, k, hxid;
	vector<int> loc(hcp[lev].size(), -1);
	vector<int> hx1r(1, eid);
	for (i = 0; i < 8; i++)
	{
		loc[hmesh[lev][eid].cnct[i]] = hmesh[lev][eid].IEN.size();
		hmesh[lev][eid].IEN.push_back(hmesh[lev][eid].cnct[i]);
	}
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].hex.size(); j++)
		{
			hxid = hcp[lev][hmesh[lev][eid].cnct[i]].hex[j];
			vector<int>::iterator it1 = find(hx1r.begin(), hx1r.end(), hxid);
			if (it1 == hx1r.end())
			{
				hx1r.push_back(hxid);
				for (k = 0; k < 8; k++)
				{
					vector<int>::iterator it = find(hmesh[lev][eid].IEN.begin(), hmesh[lev][eid].IEN.end(), hmesh[lev][hxid].cnct[k]);
					if (it == hmesh[lev][eid].IEN.end())
					{
						loc[hmesh[lev][hxid].cnct[k]] = hmesh[lev][eid].IEN.size();
						hmesh[lev][eid].IEN.push_back(hmesh[lev][hxid].cnct[k]);
					}
				}
			}
		}
	}
	for (i = 0; i < hmesh[lev][eid].bemat.size(); i++)
	{
		hmesh[lev][eid].bemat[i].clear();
	}
	hmesh[lev][eid].bemat.clear();
	hmesh[lev][eid].bemat.resize(hmesh[lev][eid].IEN.size(), vector<double>(64, 0.));
	//8 body points
	double w[2] = { 2. / 3., 1. / 3. };
	double a[8] = { w[0] * w[0] * w[0], w[1] * w[0] * w[0], w[1] * w[1] * w[0], w[1] * w[0] * w[0],
		w[1] * w[0] * w[0], w[1] * w[1] * w[0], w[1] * w[1] * w[1], w[1] * w[1] * w[0] };
	int bpi[8][8] = { { 0, 1, 2, 3, 4, 5, 6, 7 }, { 1, 2, 3, 0, 5, 6, 7, 4 }, { 2, 3, 0, 1, 6, 7, 4, 5 }, { 3, 0, 1, 2, 7, 4, 5, 6 },
	{ 4, 5, 6, 7, 0, 1, 2, 3 }, { 5, 6, 7, 4, 1, 2, 3, 0 }, { 6, 7, 4, 5, 2, 3, 0, 1 }, { 7, 4, 5, 6, 3, 0, 1, 2 } };
	vector<array<array<double, 8>, 8>> bpm(hx1r.size());
	vector<array<array<int, 8>, 8>> bpmap(hx1r.size());
	for (i = 0; i < hx1r.size(); i++)//which element
	{
		for (j = 0; j < 8; j++)//which body point, bezier
		{
			for (k = 0; k < 8; k++)//which local corner point, b-splines
			{
				bpm[i][j][k] = a[k];
				bpmap[i][j][k] = loc[hmesh[lev][hx1r[i]].cnct[bpi[j][k]]];
			}
		}
	}
	int layer[4] = { 0, 16, 32, 48 };
	int bpbz[8] = { 5 + layer[1], 6 + layer[1], 10 + layer[1], 9 + layer[1], 5 + layer[2], 6 + layer[2], 10 + layer[2], 9 + layer[2] };
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 8; j++)
		{
			hmesh[lev][eid].bemat[bpmap[0][i][j]][bpbz[i]] = bpm[0][i][j];
		}
	}
	int pos1, pos2;
	//4*6 face points
	//double af[4][4] = { { w[0] * w[0], w[1] * w[0], w[1] * w[0], w[1] * w[1] }, { w[1] * w[0], w[0] * w[0], w[1] * w[1], w[1] * w[0] },
	//{ w[1] * w[0], w[1] * w[1], w[0] * w[0], w[1] * w[0] }, { w[1] * w[1], w[1] * w[0], w[1] * w[0], w[0] * w[0] } };
	double af[4] = { w[0] * w[0], w[1] * w[0], w[1] * w[1], w[1] * w[0] };
	int f_cnct[4][4] = { { 0, 1, 2, 3 }, { 1, 2, 3, 0 }, { 2, 3, 0, 1 }, {3,0,1,2} };
	int fci[6][4] = { { 0, 1, 3, 2 }, { 0, 1, 4, 5 }, { 1, 2, 5, 6 }, { 3, 2, 7, 6 }, { 0, 3, 4, 7 }, { 4, 5, 7, 6 } };
	int fcbz[6][4] = { { 5, 6, 9, 10 }, { 1 + layer[1], 2 + layer[1], 1 + layer[2], 2 + layer[2] }, { 7 + layer[1], 11 + layer[1], 7 + layer[2], 11 + layer[2] },
	{ 13 + layer[1], 14 + layer[1], 13 + layer[2], 14 + layer[2] }, { 4 + layer[1], 8 + layer[1], 4 + layer[2], 8 + layer[2] }, { 5 + layer[3], 6 + layer[3], 9 + layer[3], 10 + layer[3] } };
	//find all boundary faces
	vector<int> fc_b;
	vector<array<array<double, 4>, 4>> fpm;
	vector<array<array<int, 4>, 4>> fpmap;//reference is face, note that reference is not solid
	for (i = 0; i < hx1r.size(); i++)
	{
		for (j = 0; j < 6; j++)
		{
			if (hface[lev][hmesh[lev][hx1r[i]].face[j]].hex.size() == 1)
			{
				int fcid(hmesh[lev][hx1r[i]].face[j]);
				fc_b.push_back(fcid);
				array<array<double, 4>, 4> fpm_tmp;
				array<array<int, 4>, 4> fpmap_tmp;
				for (int j1 = 0; j1 < 4; j1++)//Bezier
				{
					for (k = 0; k < 4; k++)//B-splines
					{
						fpm_tmp[j1][k] = af[k];
						fpmap_tmp[j1][k] = loc[hface[lev][fcid].cnct[f_cnct[j1][k]]];
					}
				}
				fpm.push_back(fpm_tmp);
				fpmap.push_back(fpmap_tmp);
			}
		}
	}
	//determine coefs
	for (i = 0; i < 6; i++)
	{
		uint nhex = hface[lev][hmesh[lev][eid].face[i]].hex.size();
		if (nhex == 2)
		{
			for (j = 0; j < nhex; j++)
			{
				hxid = hface[lev][hmesh[lev][eid].face[i]].hex[j];
				vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
				pos1 = it - hx1r.begin();
				for (int j1 = 0; j1 < 4; j1++)
				{
					int* it1 = find(hmesh[lev][hxid].cnct, hmesh[lev][hxid].cnct + 8, hmesh[lev][eid].cnct[fci[i][j1]]);
					pos2 = it1 - hmesh[lev][hxid].cnct;
					for (k = 0; k < 8; k++)
					{
						hmesh[lev][eid].bemat[bpmap[pos1][pos2][k]][fcbz[i][j1]] += bpm[pos1][pos2][k] / nhex;
					}
				}
			}
		}
		else if (nhex == 1)//new, for boundary face
		{
			int fcid = hmesh[lev][eid].face[i];
			vector<int>::iterator it = find(fc_b.begin(), fc_b.end(), fcid);
			pos1 = it - fc_b.begin();
			for (int j1 = 0; j1 < 4; j1++)
			{
				int* it1 = find(hface[lev][fcid].cnct, hface[lev][fcid].cnct + 4, hmesh[lev][eid].cnct[fci[i][j1]]);
				pos2 = it1 - hface[lev][fcid].cnct;
				for (k = 0; k < 4; k++)
				{
					hmesh[lev][eid].bemat[fpmap[pos1][pos2][k]][fcbz[i][j1]] = fpm[pos1][pos2][k];
				}
			}
			/*int* it1 = find(hface[lev][fcid].cnct, hface[lev][fcid].cnct + 4, hmesh[lev][eid].cnct[fci[i][0]]);
			int start(it1 - hface[lev][fcid].cnct);
			int fcbz1[4] = { fcbz[i][0], fcbz[i][1], fcbz[i][3], fcbz[i][2] };
			if (start == 1)
			{
				fcbz1[0] = fcbz[i][2]; fcbz1[1] = fcbz[i][0]; fcbz1[2] = fcbz[i][1]; fcbz1[3] = fcbz[i][3];
			}
			else if (start == 2)
			{
				fcbz1[0] = fcbz[i][3]; fcbz1[1] = fcbz[i][2]; fcbz1[2] = fcbz[i][0]; fcbz1[3] = fcbz[i][1];
			}
			else if (start == 3)
			{
				fcbz1[0] = fcbz[i][1]; fcbz1[1] = fcbz[i][3]; fcbz1[2] = fcbz[i][2]; fcbz1[3] = fcbz[i][0];
			}
			if (hface[lev][fcid].cnct[(start + 1) % 4] != hmesh[lev][eid].cnct[fci[i][1]])
			{
				if (start == 0)
				{
					fcbz1[0] = fcbz[i][0]; fcbz1[1] = fcbz[i][2]; fcbz1[2] = fcbz[i][3]; fcbz1[3] = fcbz[i][1];
				}
				else if (start == 1)
				{
					fcbz1[0] = fcbz[i][1]; fcbz1[1] = fcbz[i][0]; fcbz1[2] = fcbz[i][2]; fcbz1[3] = fcbz[i][3];
				}
				else if (start == 2)
				{
					fcbz1[0] = fcbz[i][3]; fcbz1[1] = fcbz[i][1]; fcbz1[2] = fcbz[i][0]; fcbz1[3] = fcbz[i][2];
				}
				else if (start == 3)
				{
					fcbz1[0] = fcbz[i][2]; fcbz1[1] = fcbz[i][3]; fcbz1[2] = fcbz[i][1]; fcbz1[3] = fcbz[i][0];
				}
			}
			for (int j1 = 0; j1 < 4; j1++)
			{
				for (k = 0; k < 4; k++)
				{
					hmesh[lev][eid].bemat[fpmap[pos1][j1][k]][fcbz1[j1]] = fpm[pos1][j1][k];
				}
			}*/
		}
	}
	//2*12 edge points
	int edi[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 } };
	int edbz[12][2] = { { 1, 2 }, { 7, 11 }, { 14, 13 }, { 8, 4 }, { 0 + layer[1], 0 + layer[2] }, { 3 + layer[1], 3 + layer[2] },
	{ 15 + layer[1], 15 + layer[2] }, { 12 + layer[1], 12 + layer[2] }, { 1 + layer[3], 2 + layer[3] }, { 7 + layer[3], 11 + layer[3] }, { 14 + layer[3], 13 + layer[3] }, { 8 + layer[3], 4 + layer[3] } };
	for (i = 0; i < 12; i++)
	{
		if (hedge[lev][hmesh[lev][eid].edge[i]].type != 1)//non-boundary
		{
			uint nhex = hedge[lev][hmesh[lev][eid].edge[i]].hex.size();
			for (j = 0; j<hedge[lev][hmesh[lev][eid].edge[i]].hex.size(); j++)
			{
				hxid = hedge[lev][hmesh[lev][eid].edge[i]].hex[j];
				vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
				pos1 = it - hx1r.begin();
				int* it1 = find(hmesh[lev][hxid].cnct, hmesh[lev][hxid].cnct + 8, hmesh[lev][eid].cnct[edi[i][0]]);
				pos2 = it1 - hmesh[lev][hxid].cnct;
				for (k = 0; k < 8; k++)
				{
					hmesh[lev][eid].bemat[bpmap[pos1][pos2][k]][edbz[i][0]] += bpm[pos1][pos2][k] / nhex;
				}
				int* it2 = find(hmesh[lev][hxid].cnct, hmesh[lev][hxid].cnct + 8, hmesh[lev][eid].cnct[edi[i][1]]);
				pos2 = it2 - hmesh[lev][hxid].cnct;
				for (k = 0; k < 8; k++)
				{
					hmesh[lev][eid].bemat[bpmap[pos1][pos2][k]][edbz[i][1]] += bpm[pos1][pos2][k] / nhex;
				}
			}
		}
		else if (hedge[lev][hmesh[lev][eid].edge[i]].type == 1 && hedge[lev][hmesh[lev][eid].edge[i]].sharp == 0)//boundary, non-sharp
		{
			int nfc_b(0);
			for (j = 0; j < hedge[lev][hmesh[lev][eid].edge[i]].face.size(); j++)
			{
				if (hface[lev][hedge[lev][hmesh[lev][eid].edge[i]].face[j]].type == 1) nfc_b++;
			}
			for (j = 0; j<hedge[lev][hmesh[lev][eid].edge[i]].face.size(); j++)
			{
				int fcid(hedge[lev][hmesh[lev][eid].edge[i]].face[j]);
				if (hface[lev][fcid].type == 1)//boundary face
				{
					vector<int>::iterator it = find(fc_b.begin(), fc_b.end(), fcid);
					pos1 = it - fc_b.begin();
					for (int j1 = 0; j1 < 2; j1++)
					{
						int* it1 = find(hface[lev][fcid].cnct, hface[lev][fcid].cnct + 4, hmesh[lev][eid].cnct[edi[i][j1]]);
						pos2 = it1 - hface[lev][fcid].cnct;
						for (k = 0; k < 4; k++)
						{
							hmesh[lev][eid].bemat[fpmap[pos1][pos2][k]][edbz[i][j1]] += fpm[pos1][pos2][k]/nfc_b;
						}
					}
				}
			}
		}
		else if (hedge[lev][hmesh[lev][eid].edge[i]].type == 1 && hedge[lev][hmesh[lev][eid].edge[i]].sharp == 1)//boundary, sharp edge
		{
			int bid[2] = { loc[hmesh[lev][eid].cnct[edi[i][0]]], loc[hmesh[lev][eid].cnct[edi[i][1]]] };
			hmesh[lev][eid].bemat[bid[0]][edbz[i][0]] = w[0];
			hmesh[lev][eid].bemat[bid[0]][edbz[i][1]] = w[1];
			hmesh[lev][eid].bemat[bid[1]][edbz[i][0]] = w[1];
			hmesh[lev][eid].bemat[bid[1]][edbz[i][1]] = w[0];
		}
	}
	//8 corner points
	int cnbz[8] = { 0, 3, 15, 12, 0 + layer[3], 3 + layer[3], 15 + layer[3], 12 + layer[3] };
	for (i = 0; i < 8; i++)
	{
		if (hcp[lev][hmesh[lev][eid].cnct[i]].type != 1)//non-boundary vertex points
		{
			uint nhex = hcp[lev][hmesh[lev][eid].cnct[i]].hex.size();
			for (j = 0; j<nhex; j++)
			{
				hxid = hcp[lev][hmesh[lev][eid].cnct[i]].hex[j];
				vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
				pos1 = it - hx1r.begin();
				int* it1 = find(hmesh[lev][hxid].cnct, hmesh[lev][hxid].cnct + 8, hmesh[lev][eid].cnct[i]);
				pos2 = it1 - hmesh[lev][hxid].cnct;
				for (k = 0; k < 8; k++)
				{
					hmesh[lev][eid].bemat[bpmap[pos1][pos2][k]][cnbz[i]] += bpm[pos1][pos2][k] / nhex;
				}
			}
		}
		else if (hcp[lev][hmesh[lev][eid].cnct[i]].type == 1 && hcp[lev][hmesh[lev][eid].cnct[i]].sharp==0)//boundary vertex points, non-sharp
		{
			int nfc_b(0);
			for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].face.size(); j++)
			{
				if (hface[lev][hcp[lev][hmesh[lev][eid].cnct[i]].face[j]].type == 1) nfc_b++;
			}
			for (j = 0; j<hcp[lev][hmesh[lev][eid].cnct[i]].face.size(); j++)
			{
				int fcid(hcp[lev][hmesh[lev][eid].cnct[i]].face[j]);
				if (hface[lev][fcid].type == 1)//boundary face
				{
					vector<int>::iterator it = find(fc_b.begin(), fc_b.end(), fcid);
					pos1 = it - fc_b.begin();
					int* it1 = find(hface[lev][fcid].cnct, hface[lev][fcid].cnct + 4, hmesh[lev][eid].cnct[i]);
					pos2 = it1 - hface[lev][fcid].cnct;
					for (k = 0; k < 4; k++)
					{
						hmesh[lev][eid].bemat[fpmap[pos1][pos2][k]][cnbz[i]] += fpm[pos1][pos2][k] / nfc_b;
					}
				}
			}
		}
		else if (hcp[lev][hmesh[lev][eid].cnct[i]].type == 1 && hcp[lev][hmesh[lev][eid].cnct[i]].sharp == 1)//boundary vertex points, sharp edge
		{
			//hmesh[lev][eid].bemat[loc[hmesh[lev][eid].cnct[i]]][cnbz[i]] = 1.;
			vector<int> shp_ed;
			for (j = 0; j < hcp[lev][hmesh[lev][eid].cnct[i]].edge.size(); j++)
			{
				if (hedge[lev][hcp[lev][hmesh[lev][eid].cnct[i]].edge[j]].sharp == 1)
				{
					shp_ed.push_back(hcp[lev][hmesh[lev][eid].cnct[i]].edge[j]);
				}
			}
			if (shp_ed.size() == 2)
			{
				for (j = 0; j < shp_ed.size(); j++)
				{
					if (hedge[lev][shp_ed[j]].pt[0] == hmesh[lev][eid].cnct[i])
					{
						hmesh[lev][eid].bemat[loc[hedge[lev][shp_ed[j]].pt[0]]][cnbz[i]] = 2. / 3.;
						hmesh[lev][eid].bemat[loc[hedge[lev][shp_ed[j]].pt[1]]][cnbz[i]] = 1. / 6.;
					}
					else
					{
						hmesh[lev][eid].bemat[loc[hedge[lev][shp_ed[j]].pt[1]]][cnbz[i]] = 2. / 3.;
						hmesh[lev][eid].bemat[loc[hedge[lev][shp_ed[j]].pt[0]]][cnbz[i]] = 1. / 6.;
					}
				}
			}
			else
			{
				cerr << "# of sharp edges of the point is wrong!\n";
				cout << shp_ed.size() << "\n";
				getchar();
			}
		}
		else if (hcp[lev][hmesh[lev][eid].cnct[i]].type == 1 && hcp[lev][hmesh[lev][eid].cnct[i]].sharp == 2)//boundary vertex points, sharp corner
		{
			hmesh[lev][eid].bemat[loc[hmesh[lev][eid].cnct[i]]][cnbz[i]] = 1.;
		}
	}
}

//void TruncatedTspline_3D::ConstructConnect(int lev)//no need for all of them
//{
//	uint i, j;
//	//construct edges
//	for (i = 0; i<hmesh[lev].size(); i++)
//	{
//		for (j = 0; j<4; j++)
//		{
//			Edge3D edtmp;
//			edtmp.pt[0] = hmesh[lev][i].cnct[j];
//			edtmp.pt[1] = hmesh[lev][i].cnct[(j + 1) % 4];
//			vector<Edge3D>::iterator it = find(hedge[lev].begin(), hedge[lev].end(), edtmp);
//			int edid(it - hedge[lev].begin());
//			if (it == hedge[lev].end())
//			{
//				hedge[lev].push_back(edtmp);
//			}
//			hmesh[lev][i].edge[j] = edid;
//		}
//		for (j = 0; j<4; j++)
//		{
//			Edge3D edtmp;
//			edtmp.pt[0] = hmesh[lev][i].cnct[j];
//			edtmp.pt[1] = hmesh[lev][i].cnct[j + 4];
//			vector<Edge3D>::iterator it = find(hedge[lev].begin(), hedge[lev].end(), edtmp);
//			int edid(it - hedge[lev].begin());
//			if (it == hedge[lev].end())
//			{
//				hedge[lev].push_back(edtmp);
//			}
//			hmesh[lev][i].edge[j + 4] = edid;
//		}
//		for (j = 0; j<4; j++)
//		{
//			Edge3D edtmp;
//			edtmp.pt[0] = hmesh[lev][i].cnct[j + 4];
//			edtmp.pt[1] = hmesh[lev][i].cnct[(j + 1) % 4 + 4];
//			vector<Edge3D>::iterator it = find(hedge[lev].begin(), hedge[lev].end(), edtmp);
//			int edid(it - hedge[lev].begin());
//			if (it == hedge[lev].end())
//			{
//				hedge[lev].push_back(edtmp);
//			}
//			hmesh[lev][i].edge[j + 8] = edid;
//		}
//	}
//	//construct faces
//	for (i = 0; i<hmesh[lev].size(); i++)
//	{
//		//one bottom face
//		Face3D fc1;
//		for (j = 0; j<4; j++)
//		{
//			fc1.cnct[j] = hmesh[lev][i].cnct[j];
//			fc1.edge[j] = hmesh[lev][i].edge[j];
//		}
//		vector<Face3D>::iterator it1 = find(hface[lev].begin(), hface[lev].end(), fc1);
//		int fc1id(it1 - hface[lev].begin());
//		if (it1 == hface[lev].end())
//		{
//			hface[lev].push_back(fc1);
//		}
//		hmesh[lev][i].face[0] = fc1id;
//		//4 side faces
//		for (j = 0; j<4; j++)
//		{
//			Face3D fc;
//			for (int k = 0; k<4; k++)
//			{
//				fc.cnct[k] = hmesh[lev][i].cnct[fc_cnct[j][k]];
//				fc.edge[k] = hmesh[lev][i].edge[ed_cnct[j][k]];
//			}
//			vector<Face3D>::iterator it = find(hface[lev].begin(), hface[lev].end(), fc);
//			int fcid(it - hface[lev].begin());
//			if (it == hface[lev].end())
//			{
//				hface[lev].push_back(fc);
//			}
//			hmesh[lev][i].face[j + 1] = fcid;
//		}
//		//one top face
//		Face3D fc2;
//		for (j = 0; j<4; j++)
//		{
//			fc2.cnct[j] = hmesh[lev][i].cnct[j + 4];
//			fc2.edge[j] = hmesh[lev][i].edge[j + 8];
//		}
//		vector<Face3D>::iterator it2 = find(hface[lev].begin(), hface[lev].end(), fc2);
//		int fc2id(it2 - hface[lev].begin());
//		if (it2 == hface[lev].end())
//		{
//			hface[lev].push_back(fc2);
//		}
//		hmesh[lev][i].face[5] = fc2id;
//	}
//	//vertex-to-hex, edge-to-hex, face-to-hex
//	for (i = 0; i<hmesh[lev].size(); i++)
//	{
//		for (j = 0; j<8; j++)
//		{
//			hcp[lev][hmesh[lev][i].cnct[j]].hex.push_back(i);
//		}
//		for (j = 0; j<12; j++)
//		{
//			hedge[lev][hmesh[lev][i].edge[j]].hex.push_back(i);
//		}
//		for (j = 0; j<6; j++)
//		{
//			hface[lev][hmesh[lev][i].face[j]].hex.push_back(i);
//		}
//	}
//	//vertex-to-face, edge-to-face
//	for (i = 0; i<hface[lev].size(); i++)
//	{
//		for (j = 0; j<4; j++)
//		{
//			hcp[lev][hface[lev][i].cnct[j]].face.push_back(i);
//			hedge[lev][hface[lev][i].edge[j]].face.push_back(i);
//		}
//	}
//	//vertex-to-edge
//	for (i = 0; i<hedge[lev].size(); i++)
//	{
//		for (j = 0; j<2; j++)
//		{
//			hcp[lev][hedge[lev][i].pt[j]].edge.push_back(i);
//		}
//	}
//	//find BC face, edge, vertex
//	int ed0[6][4] = { { 4, 5, 6, 7 }, { 1, 3, 9, 11 }, { 0, 2, 8, 10 }, { 1, 3, 9, 11 }, { 0, 2, 8, 10 }, { 4, 5, 6, 7 } };//order could be wrong, but doesn't matter
//	for (i = 0; i<hface[lev].size(); i++)
//	{
//		if (hface[lev][i].hex.size() == 1)
//		{
//			hface[lev][i].type = 1;
//			hmesh[lev][hface[lev][i].hex[0]].type = 1;
//			for (j = 0; j<4; j++)
//			{
//				hcp[lev][hface[lev][i].cnct[j]].type = 1;
//				hedge[lev][hface[lev][i].edge[j]].type = 1;
//			}
//		}
//	}
//	////find extraordinary edges and vertices
//	//for (i = 0; i<hedge[lev].size(); i++)
//	//{
//	//	if (hedge[lev][i].type != 1 && hedge[lev][i].hex.size() != 4)
//	//	{
//	//		hedge[lev][i].type = 2;
//	//		if (hcp[lev][hedge[lev][i].pt[0]].type != 1)
//	//			hcp[lev][hedge[lev][i].pt[0]].type = 3;
//	//		if (hcp[lev][hedge[lev][i].pt[1]].type != 1)
//	//			hcp[lev][hedge[lev][i].pt[1]].type = 3;
//	//	}
//	//}
//	////boundary
//	//for (i = 0; i<hcp[lev].size(); i++)
//	//{
//	//	if (hcp[lev][i].type == 1)
//	//	{
//	//		int val(0);
//	//		for (j = 0; j<hcp[lev][i].face.size(); j++)
//	//		{
//	//			if (hface[lev][hcp[lev][i].face[j]].type == 1) val++;
//	//		}
//	//		if (val == 3 || val>4) hcp[lev][i].type = 13;
//	//	}
//	//}
//	////find irregular elements
//	//for (i = 0; i<hmesh[lev].size(); i++)
//	//{
//	//	if (hmesh[lev][i].type != 1)
//	//	{
//	//		for (j = 0; j<12; j++)
//	//		{
//	//			if (hedge[lev][hmesh[lev][i].edge[j]].type == 2)
//	//			{
//	//				hmesh[lev][i].type = 2;
//	//				break;
//	//			}
//	//		}
//	//	}
//	//}
//}

void TruncatedTspline_3D::ConstructConnect(int lev)//no need for all of them
{
	uint i, j;
	//clear
	for (i = 0; i < hcp[lev].size(); i++)
	{
		//hcp[lev][i].type = 0;
		hcp[lev][i].edge.clear();
		hcp[lev][i].face.clear();
		hcp[lev][i].hex.clear();
	}
	for (i = 0; i < hedge[lev].size(); i++)
	{
		//hedge[lev][i].type = 0;
		hedge[lev][i].face.clear();
		hedge[lev][i].hex.clear();
	}
	for (i = 0; i < hface[lev].size(); i++)
	{
		//hface[lev][i].type = 0;
		hface[lev][i].hex.clear();
	}
	//for (i = 0; i < hmesh[lev].size(); i++)
	//{
	//	//hmesh[lev][i].bc_lev=0;
	//	//hmesh[lev][i].type = 0;
	//}
	//vertex-to-hex, edge-to-hex, face-to-hex
	for (i = 0; i<hmesh[lev].size(); i++)
	{
		for (j = 0; j<8; j++)
		{
			hcp[lev][hmesh[lev][i].cnct[j]].hex.push_back(i);
		}
		for (j = 0; j<12; j++)
		{
			hedge[lev][hmesh[lev][i].edge[j]].hex.push_back(i);
		}
		for (j = 0; j<6; j++)
		{
			hface[lev][hmesh[lev][i].face[j]].hex.push_back(i);
		}
	}
	//vertex-to-face, edge-to-face
	for (i = 0; i<hface[lev].size(); i++)
	{
		for (j = 0; j<4; j++)
		{
			hcp[lev][hface[lev][i].cnct[j]].face.push_back(i);
			hedge[lev][hface[lev][i].edge[j]].face.push_back(i);
		}
	}
	//vertex-to-edge
	for (i = 0; i<hedge[lev].size(); i++)
	{
		for (j = 0; j<2; j++)
		{
			hcp[lev][hedge[lev][i].pt[j]].edge.push_back(i);
		}
	}
	//check boundary type
	for (i = 0; i < hface[lev].size(); i++)
	{
		int nbc(0);
		for (j = 0; j < 4; j++)
		{
			if (hcp[lev][hface[lev][i].cnct[j]].type == 1) nbc++;
		}
		if (nbc == 4)
		{
			hface[lev][i].type = 1;
			hmesh[lev][hface[lev][i].hex[0]].type = 1;
			//for (j = 0; j < 4; j++)
			//	hedge[lev][hface[lev][i].edge[j]].type = 1;
		}
	}
	for (i = 0; i < hedge[lev].size(); i++)
	{
		int nbc(0), nshp(0);
		for (j = 0; j < 2; j++)
		{
			if (hcp[lev][hedge[lev][i].pt[j]].type == 1) nbc++;
			if (hcp[lev][hedge[lev][i].pt[j]].sharp != 0) nshp++;
		}
		if (nbc == 2)
		{
			hedge[lev][i].type = 1;
			if (nshp == 2) hedge[lev][i].sharp = 1;
		}
	}

	//new, additional boundary elements
	for (i = 0; i < hmesh[lev].size(); i++)
	{
		if (hmesh[lev][i].type != 1)
		{
			for (j = 0; j < 8; j++)
			{
				if (hcp[lev][hmesh[lev][i].cnct[j]].type == 1)
				{
					hmesh[lev][i].type = 1; break;
				}
			}
		}
	}
	for (i = 0; i<hcp[lev].size(); i++)
	{
		if (hcp[lev][i].type == 1)
		{
			int count(0);
			for (j = 0; j < hcp[lev][i].edge.size(); j++)
			{
				if (hedge[lev][hcp[lev][i].edge[j]].type == 2) count++;
			}
			if (count == 1) hcp[lev][i].bcxp = 1;
			else if (count>1) hcp[lev][i].bcxp = 2;
		}
	}
	

	////find BC face, edge, vertex, not necessary
	//for (i = 0; i<hface[lev].size(); i++)
	//{
	//	if (hface[lev][i].hex.size() == 1)
	//	{
	//		hface[lev][i].type = 1;
	//		hmesh[lev][hface[lev][i].hex[0]].type = 1;
	//		for (j = 0; j<4; j++)
	//		{
	//			hcp[lev][hface[lev][i].cnct[j]].type = 1;
	//			hedge[lev][hface[lev][i].edge[j]].type = 1;
	//		}
	//	}
	//}
	////find extraordinary edges and vertices
	//for (i = 0; i<hedge[lev].size(); i++)
	//{
	//	if (hedge[lev][i].type != 1 && hedge[lev][i].hex.size() != 4)
	//	{
	//		hedge[lev][i].type = 2;
	//		if (hcp[lev][hedge[lev][i].pt[0]].type != 1)
	//			hcp[lev][hedge[lev][i].pt[0]].type = 3;
	//		if (hcp[lev][hedge[lev][i].pt[1]].type != 1)
	//			hcp[lev][hedge[lev][i].pt[1]].type = 3;
	//	}
	//}
	////boundary
	//for (i = 0; i<hcp[lev].size(); i++)
	//{
	//	if (hcp[lev][i].type == 1)
	//	{
	//		int val(0);
	//		for (j = 0; j<hcp[lev][i].face.size(); j++)
	//		{
	//			if (hface[lev][hcp[lev][i].face[j]].type == 1) val++;
	//		}
	//		if (val == 3 || val>4) hcp[lev][i].type = 13;
	//	}
	//}
	////find irregular elements
	//for (i = 0; i<hmesh[lev].size(); i++)
	//{
	//	if (hmesh[lev][i].type != 1)
	//	{
	//		for (j = 0; j<12; j++)
	//		{
	//			if (hedge[lev][hmesh[lev][i].edge[j]].type == 2)
	//			{
	//				hmesh[lev][i].type = 2;
	//				break;
	//			}
	//		}
	//	}
	//}
}

void TruncatedTspline_3D::ConstructFaceEdge(int lev, int eid)//no need for all of them
{
	if (hmesh[lev][eid].edge[0] != -1 || hmesh[lev][eid].face[0] != -1) return;//means already constructed
	uint i(eid), j;
	//construct edges
	for (j = 0; j<4; j++)
	{
		Edge3D edtmp;
		edtmp.pt[0] = hmesh[lev][i].cnct[j];
		edtmp.pt[1] = hmesh[lev][i].cnct[(j + 1) % 4];
		vector<Edge3D>::iterator it = find(hedge[lev].begin(), hedge[lev].end(), edtmp);
		int edid(it - hedge[lev].begin());
		if (it == hedge[lev].end())
		{
			hedge[lev].push_back(edtmp);
		}
		hmesh[lev][i].edge[j] = edid;
	}
	for (j = 0; j<4; j++)
	{
		Edge3D edtmp;
		edtmp.pt[0] = hmesh[lev][i].cnct[j];
		edtmp.pt[1] = hmesh[lev][i].cnct[j + 4];
		vector<Edge3D>::iterator it = find(hedge[lev].begin(), hedge[lev].end(), edtmp);
		int edid(it - hedge[lev].begin());
		if (it == hedge[lev].end())
		{
			hedge[lev].push_back(edtmp);
		}
		hmesh[lev][i].edge[j + 4] = edid;
	}
	for (j = 0; j<4; j++)
	{
		Edge3D edtmp;
		edtmp.pt[0] = hmesh[lev][i].cnct[j + 4];
		edtmp.pt[1] = hmesh[lev][i].cnct[(j + 1) % 4 + 4];
		vector<Edge3D>::iterator it = find(hedge[lev].begin(), hedge[lev].end(), edtmp);
		int edid(it - hedge[lev].begin());
		if (it == hedge[lev].end())
		{
			hedge[lev].push_back(edtmp);
		}
		hmesh[lev][i].edge[j + 8] = edid;
	}
	//construct faces
	//one bottom face
	Face3D fc1;
	for (j = 0; j<4; j++)
	{
		fc1.cnct[j] = hmesh[lev][i].cnct[j];
		fc1.edge[j] = hmesh[lev][i].edge[j];
	}
	vector<Face3D>::iterator it1 = find(hface[lev].begin(), hface[lev].end(), fc1);
	int fc1id(it1 - hface[lev].begin());
	if (it1 == hface[lev].end())
	{
		hface[lev].push_back(fc1);
	}
	hmesh[lev][i].face[0] = fc1id;
	//4 side faces
	for (j = 0; j<4; j++)
	{
		Face3D fc;
		for (int k = 0; k<4; k++)
		{
			fc.cnct[k] = hmesh[lev][i].cnct[fc_cnct[j][k]];
			fc.edge[k] = hmesh[lev][i].edge[ed_cnct[j][k]];
		}
		vector<Face3D>::iterator it = find(hface[lev].begin(), hface[lev].end(), fc);
		int fcid(it - hface[lev].begin());
		if (it == hface[lev].end())
		{
			hface[lev].push_back(fc);
		}
		hmesh[lev][i].face[j + 1] = fcid;
	}
	//one top face
	Face3D fc2;
	for (j = 0; j<4; j++)
	{
		fc2.cnct[j] = hmesh[lev][i].cnct[j + 4];
		fc2.edge[j] = hmesh[lev][i].edge[j + 8];
	}
	vector<Face3D>::iterator it2 = find(hface[lev].begin(), hface[lev].end(), fc2);
	int fc2id(it2 - hface[lev].begin());
	if (it2 == hface[lev].end())
	{
		hface[lev].push_back(fc2);
	}
	hmesh[lev][i].face[5] = fc2id;

	////vertex-to-hex, edge-to-hex, face-to-hex
	//for (i = 0; i<hmesh[lev].size(); i++)
	//{
	//	for (j = 0; j<8; j++)
	//	{
	//		hcp[lev][hmesh[lev][i].cnct[j]].hex.push_back(i);
	//	}
	//	for (j = 0; j<12; j++)
	//	{
	//		hedge[lev][hmesh[lev][i].edge[j]].hex.push_back(i);
	//	}
	//	for (j = 0; j<6; j++)
	//	{
	//		hface[lev][hmesh[lev][i].face[j]].hex.push_back(i);
	//	}
	//}
	////vertex-to-face, edge-to-face
	//for (i = 0; i<hface[lev].size(); i++)
	//{
	//	for (j = 0; j<4; j++)
	//	{
	//		hcp[lev][hface[lev][i].cnct[j]].face.push_back(i);
	//		hedge[lev][hface[lev][i].edge[j]].face.push_back(i);
	//	}
	//}
	////vertex-to-edge
	//for (i = 0; i<hedge[lev].size(); i++)
	//{
	//	for (j = 0; j<2; j++)
	//	{
	//		hcp[lev][hedge[lev][i].pt[j]].edge.push_back(i);
	//	}
	//}
	////find BC face, edge, vertex
	//int ed0[6][4] = { { 4, 5, 6, 7 }, { 1, 3, 9, 11 }, { 0, 2, 8, 10 }, { 1, 3, 9, 11 }, { 0, 2, 8, 10 }, { 4, 5, 6, 7 } };//order could be wrong, but doesn't matter
	//for (i = 0; i<hface[lev].size(); i++)
	//{
	//	if (hface[lev][i].hex.size() == 1)
	//	{
	//		hface[lev][i].type = 1;
	//		hmesh[lev][hface[lev][i].hex[0]].type = 1;
	//		for (j = 0; j<4; j++)
	//		{
	//			hcp[lev][hface[lev][i].cnct[j]].type = 1;
	//			hedge[lev][hface[lev][i].edge[j]].type = 1;
	//		}
	//	}
	//}
	////find extraordinary edges and vertices
	//for (i = 0; i<hedge[lev].size(); i++)
	//{
	//	if (hedge[lev][i].type != 1 && hedge[lev][i].hex.size() != 4)
	//	{
	//		hedge[lev][i].type = 2;
	//		if (hcp[lev][hedge[lev][i].pt[0]].type != 1)
	//			hcp[lev][hedge[lev][i].pt[0]].type = 3;
	//		if (hcp[lev][hedge[lev][i].pt[1]].type != 1)
	//			hcp[lev][hedge[lev][i].pt[1]].type = 3;
	//	}
	//}
	////boundary
	//for (i = 0; i<hcp[lev].size(); i++)
	//{
	//	if (hcp[lev][i].type == 1)
	//	{
	//		int val(0);
	//		for (j = 0; j<hcp[lev][i].face.size(); j++)
	//		{
	//			if (hface[lev][hcp[lev][i].face[j]].type == 1) val++;
	//		}
	//		if (val == 3 || val>4) hcp[lev][i].type = 13;
	//	}
	//}
	////find irregular elements
	//for (i = 0; i<hmesh[lev].size(); i++)
	//{
	//	if (hmesh[lev][i].type != 1)
	//	{
	//		for (j = 0; j<12; j++)
	//		{
	//			if (hedge[lev][hmesh[lev][i].edge[j]].type == 2)
	//			{
	//				hmesh[lev][i].type = 2;
	//				break;
	//			}
	//		}
	//	}
	//}
}

//void TruncatedTspline_3D::Refine_Ghost(const vector<array<int, 2>>& rfid)
//{
//	//has to be irregular elements
//	for (uint i = 0; i < rfid.size(); i++)
//	{
//		if (hmesh[rfid[i][0]][rfid[i][1]].type == 0)
//		{
//			PatchRefine_Regular(rfid[i][0], rfid[i][1]);
//		}
//		else if (hmesh[rfid[i][0]][rfid[i][1]].type == 1)
//		{
//			PatchRefine_Boundary(rfid[i][0], rfid[i][1]);
//		}
//		else if (hmesh[rfid[i][0]][rfid[i][1]].type == 2)
//		{
//			PatchRefine_Irregular(rfid[i][0], rfid[i][1]);
//		}
//		for (uint j = 0; j < hmesh[rfid[i][0]][rfid[i][1]].chd.size(); j++)
//		{
//			hmesh[rfid[i][0] + 1][hmesh[rfid[i][0]][rfid[i][1]].chd[j]].act = 0;
//			hmesh[rfid[i][0] + 1][hmesh[rfid[i][0]][rfid[i][1]].chd[j]].ghost = 1;
//		}
//		hmesh[rfid[i][0]][rfid[i][1]].act = 1;
//	}
//}

void TruncatedTspline_3D::Refine_Ghost(const vector<array<int, 2>>& rfid)
{
	//has to be irregular elements
	for (uint i = 0; i < rfid.size(); i++)
	{
		if (hmesh[rfid[i][0]][rfid[i][1]].type == 0 || hmesh[rfid[i][0]][rfid[i][1]].type == 2)
		{
			if (hmesh[rfid[i][0]][rfid[i][1]].chd.size() == 0)
			{
				PatchRefine_Irregular(rfid[i][0], rfid[i][1]);
			}
		}
		else if (hmesh[rfid[i][0]][rfid[i][1]].type == 1)
		{
			if (hmesh[rfid[i][0]][rfid[i][1]].chd.size() == 0)
			{
				PatchRefine_Boundary(rfid[i][0], rfid[i][1]);
			}
		}
		//else if (hmesh[rfid[i][0]][rfid[i][1]].type == 2)
		//{
		//	PatchRefine_Irregular(rfid[i][0], rfid[i][1]);
		//}
		for (uint j = 0; j < hmesh[rfid[i][0]][rfid[i][1]].chd.size(); j++)
		{
			hmesh[rfid[i][0] + 1][hmesh[rfid[i][0]][rfid[i][1]].chd[j]].act = 0;
			hmesh[rfid[i][0] + 1][hmesh[rfid[i][0]][rfid[i][1]].chd[j]].ghost = 1;
		}
		hmesh[rfid[i][0]][rfid[i][1]].act = 1;
	}
}

void TruncatedTspline_3D::Selection(int lev)
{
	uint i, j, k, m;
	for (i = 0; i < hcp[lev].size(); i++)
	{
		hcp[lev][i].act = 0;
		hcp[lev][i].supp.clear();
	}
	for (i = 0; i < hmesh[lev].size(); i++)
	{
		if (hmesh[lev][i].ghost == 0)
		{
			for (j = 0; j < hmesh[lev][i].IEN.size(); j++)
			{
				hcp[lev][hmesh[lev][i].IEN[j]].supp.push_back(i);
			}
		}
	}
	//for (uint i = 0; i < hcp[lev].size(); i++)
	//{
	//	if (hcp[lev][i].supp.size() == 64)//tmp
	//	{
	//		int flag(0);
	//		for (uint j = 0; j < hcp[lev][i].supp.size(); j++)
	//		{
	//			if (hmesh[lev][hcp[lev][i].supp[j]].act == 1)
	//			{
	//				flag = 1; break;
	//			}
	//		}
	//		if (flag == 1)
	//		{
	//			hcp[lev][i].act = 1;
	//		}
	//	}
	//}

	for (i = 0; i < hcp[lev].size(); i++)
	{
		int flag(0);
		for (j = 0; j < hcp[lev][i].hex.size(); j++)
		{
			int hxid(hcp[lev][i].hex[j]);
			if (hmesh[lev][hxid].ghost == 1)
			{
				flag = 1; break;
			}
			for (k = 0; k < 8; k++)
			{
				int pid(hmesh[lev][hxid].cnct[k]);
				for (m = 0; m < hcp[lev][pid].hex.size(); m++)
				{
					if (hmesh[lev][hcp[lev][pid].hex[m]].ghost == 1)
					{
						flag = 1; break;
					}
				}
				if (flag == 1) break;
			}
			if (flag == 1) break;
		}
		if (flag == 0)//possibly active
		{
			int act_hx(0);
			for (j = 0; j < hcp[lev][i].supp.size(); j++)
			{
				if (hmesh[lev][hcp[lev][i].supp[j]].act == 1)
				{
					act_hx = 1; break;
				}
			}
			if (act_hx == 1)
			{
				hcp[lev][i].act = 1;
			}
		}
	}
}

void TruncatedTspline_3D::SetSupport(int lev)
{
	uint i, j;
	for (i = 0; i < hcp[lev].size(); i++)
	{
		hcp[lev][i].supp.clear();
	}
	for (i = 0; i < hmesh[lev].size(); i++)
	{
		if (hmesh[lev][i].ghost == 0)
		{
			for (j = 0; j < hmesh[lev][i].IEN.size(); j++)
			{
				hcp[lev][hmesh[lev][i].IEN[j]].supp.push_back(i);
			}
		}
	}
}

void TruncatedTspline_3D::Select()
{
	uint lev, i, j, k, m;
	for (lev = 0; lev < hcp.size(); lev++)
	{
		for (i = 0; i < hcp[lev].size(); i++)
		{
			hcp[lev][i].act = 0;
		}
		////set supp for new level
		//for (i = 0; i < hmesh[lev].size(); i++)
		//{
		//	if (hmesh[lev][i].ghost == 0)//because IEN only available for non-ghost elements
		//	{
		//		for (j = 0; j < hmesh[lev][i].IEN.size(); j++)
		//		{
		//			hcp[lev][hmesh[lev][i].IEN[j]].supp.push_back(i);
		//		}
		//	}
		//}

		for (i = 0; i < hcp[lev].size(); i++)
		{
			int flag(0);//two-ring elements no ghost
			for (j = 0; j < hcp[lev][i].hex.size(); j++)
			{
				int hxid(hcp[lev][i].hex[j]);
				if (hmesh[lev][hxid].ghost == 1)
				{
					flag = 1; break;
				}
				for (k = 0; k < 8; k++)
				{
					int pid(hmesh[lev][hxid].cnct[k]);
					for (m = 0; m < hcp[lev][pid].hex.size(); m++)
					{
						if (hmesh[lev][hcp[lev][pid].hex[m]].ghost == 1)
						{
							flag = 1; break;
						}
					}
					if (flag == 1) break;
				}
				if (flag == 1) break;
			}
			if (flag == 0)//possibly active
			{
				int act_hx(0);
				for (j = 0; j < hcp[lev][i].supp.size(); j++)
				{
					if (hmesh[lev][hcp[lev][i].supp[j]].act == 1)
					{
						act_hx = 1; break;
					}
				}
				if (act_hx == 1)
				{
					hcp[lev][i].act = 1;
				}
			}
		}
	}
}

//void TruncatedTspline_3D::Truncate(int lev)//lev is a high level or a refined level
//{
//	if (lev <= 0) return;
//	uint i, j, k;
//	//point-wise truncation
//	for (i = 0; i < hcp[lev - 1].size(); i++)
//	{
//		if (hcp[lev - 1][i].act == 1)
//		{
//			for (j = 0; j < hcp[lev - 1][i].chd.size(); j++)
//			{
//				if (hcp[lev][hcp[lev - 1][i].chd[j]].act == 1)
//				{
//					hcp[lev - 1][i].coef[j] = 0.;
//				}
//			}
//		}
//	}
//
//	//element-wise truncation
//	for (i = 0; i < hmesh[lev].size(); i++)
//	{
//		hmesh[lev][i].trun = 0;
//		if (hmesh[lev][i].act == 1 && hmesh[lev][i].prt!=-1)//must be active high-level elements
//		{
//			int flag(0);
//			for (j = 0; j < hmesh[lev][i].IEN.size(); j++)
//			{
//				if (hcp[lev][hmesh[lev][i].IEN[j]].act == 0)
//				{
//					flag++;
//				}
//			}
//			if (flag == hmesh[lev][i].IEN.size())//all the IEN are passive, something is wrong
//			{
//				cerr << "All the basis functions on " << lev << " " << i << " are passive!\n";
//				return;
//			}
//			if (flag != 0)
//			{
//				hmesh[lev][i].trun = 1;
//				for (j = 0; j < hmesh[lev][i].IEN.size(); j++)
//				{
//					if (hcp[lev][hmesh[lev][i].IEN[j]].act == 1)
//					{
//						array<int, 2> btmp = { lev, hmesh[lev][i].IEN[j] };
//						vector<double> ctmp(hmesh[lev][i].IEN.size(),0.);
//						ctmp[j] = 1.;
//						hmesh[lev][i].IEN_act.push_back(btmp);
//						hmesh[lev][i].tmat.push_back(ctmp);
//					}
//				}
//
//				int prt(hmesh[lev][i].prt);
//				int lv(lev-1);
//				MatrixXd m0 = MatrixXd::Identity(hmesh[lev][i].IEN.size(), hmesh[lev][i].IEN.size());
//				while (prt != -1 && lv>=0)
//				{
//					MatrixXd m1 = MatrixXd::Zero(hmesh[lv][prt].IEN.size(), hmesh[lev][i].IEN.size());
//					for (j = 0; j < hmesh[lv][prt].IEN.size(); j++)
//					{
//						if (hcp[lv][hmesh[lv][prt].IEN[j]].act == 1)
//						{
//							array<int, 2> btmp = { lv, hmesh[lv][prt].IEN[j] };//tmp
//							vector<double> ctmp(hmesh[lev][i].IEN.size(),0.);//tmp
//							for (k = 0; k < hcp[lv][hmesh[lv][prt].IEN[j]].chd.size(); k++)
//							{
//								int cid(hcp[lv][hmesh[lv][prt].IEN[j]].chd[k]);
//								vector<int>::iterator it = find(hmesh[lev][i].IEN.begin(), hmesh[lev][i].IEN.end(), cid);
//								if (it != hmesh[lev][i].IEN.end())
//								{
//									int pos(it - hmesh[lev][i].IEN.begin());
//									ctmp[pos] = hcp[lv][hmesh[lv][prt].IEN[j]].coef[k];
//								}
//							}
//							hmesh[lev][i].IEN_act.push_back(btmp);
//							hmesh[lev][i].tmat.push_back(ctmp);
//						}
//					}
//					prt = hmesh[lv][prt].prt;
//					lv--;
//				}
//			}
//			//else//no need for truncation
//			//{
//			//	for (j = 0; j < hmesh[lev][i].IEN.size(); j++)
//			//	{
//			//		if (hcp[lev][hmesh[lev][i].IEN[j]].act == 1)
//			//		{
//			//			array<int, 2> btmp = { lev, hmesh[lev][i].IEN[j] };
//			//			hmesh[lev][i].IEN_act.push_back(btmp);
//			//		}
//			//	}
//			//}
//		}
//	}
//}

void TruncatedTspline_3D::Truncate(int lev)//lev is a high level or a refined level
{
	if (lev <= 0) return;
	uint i, j, k;
	//point-wise truncation
	for (i = 0; i < hcp[lev - 1].size(); i++)
	{
		//if (hcp[lev - 1][i].act == 1)
		{
			for (j = 0; j < hcp[lev - 1][i].chd.size(); j++)
			{
				if (hcp[lev][hcp[lev - 1][i].chd[j]].act == 1)
				{
					hcp[lev - 1][i].coef[j] = 0.;
				}
			}
		}
	}

	//element-wise truncation
	for (i = 0; i < hmesh[lev].size(); i++)
	{
		if (hmesh[lev][i].trun == 1)
		{
			hmesh[lev][i].trun = 0;
			hmesh[lev][i].IEN_act.clear();
			for (j = 0; j < hmesh[lev][i].tmat.size(); j++)
			{
				hmesh[lev][i].tmat[j].clear();
			}
			hmesh[lev][i].tmat.clear();
		}
		if (hmesh[lev][i].act == 1 && hmesh[lev][i].prt != -1)//must be active high-level elements
		{
			//if (hmesh[lev][i].trun == 1)//initialize
			//{
			//	hmesh[lev][i].trun = 0;
			//	hmesh[lev][i].IEN_act.clear();
			//	for (j = 0; j < hmesh[lev][i].tmat.size(); j++)
			//	{
			//		hmesh[lev][i].tmat[j].clear();
			//	}
			//	hmesh[lev][i].tmat.clear();
			//}
			int flag(0);
			for (j = 0; j < hmesh[lev][i].IEN.size(); j++)
			{
				if (hcp[lev][hmesh[lev][i].IEN[j]].act == 0)
				{
					flag++;
				}
			}
			if (flag == hmesh[lev][i].IEN.size())
			{
				cerr << "All the basis functions on " << lev << " " << i << " are passive!\n";
				return;
			}
			if (flag != 0)
			{
				hmesh[lev][i].trun = 1;
				for (j = 0; j < hmesh[lev][i].IEN.size(); j++)
				{
					if (hcp[lev][hmesh[lev][i].IEN[j]].act == 1)
					{
						array<int, 2> btmp = { lev, hmesh[lev][i].IEN[j] };
						vector<double> ctmp(hmesh[lev][i].IEN.size(), 0.);
						ctmp[j] = 1.;
						hmesh[lev][i].IEN_act.push_back(btmp);
						hmesh[lev][i].tmat.push_back(ctmp);
					}
				}

				int prt(hmesh[lev][i].prt);
				int lv(lev - 1);
				int eid0(i), lv0(lev);
				MatrixXd m0 = MatrixXd::Identity(hmesh[lv0][eid0].IEN.size(), hmesh[lev][i].IEN.size());
				while (prt != -1 && lv >= 0)
				{
					MatrixXd m1 = MatrixXd::Zero(hmesh[lv][prt].IEN.size(), hmesh[lv0][eid0].IEN.size());
					for (j = 0; j < hmesh[lv][prt].IEN.size(); j++)
					{
						for (k = 0; k < hcp[lv][hmesh[lv][prt].IEN[j]].chd.size(); k++)
						{
							int cid(hcp[lv][hmesh[lv][prt].IEN[j]].chd[k]);
							vector<int>::iterator it = find(hmesh[lv0][eid0].IEN.begin(), hmesh[lv0][eid0].IEN.end(), cid);
							if (it != hmesh[lv0][eid0].IEN.end())
							{
								int pos(it - hmesh[lv0][eid0].IEN.begin());
								m1(j, pos) = hcp[lv][hmesh[lv][prt].IEN[j]].coef[k];
							}
						}
					}
					//cout
					//if (lev == 2 && lv == 0)
					//{
					//	vector<double> sum0(hmesh[lev][i].IEN.size(), 0.);
					//	cout << "sum of m0 before:\n";
					//	for (j = 0; j < hmesh[lev][i].IEN.size(); j++)
					//	{
					//		for (k = 0; k < hmesh[lv0][eid0].IEN.size(); k++)
					//		{
					//			sum0[j] += m0(k, j);
					//		}
					//		cout << sum0[j] << "\n";
					//	}
					//	//vector<double> sum1(hmesh[lv0][eid0].IEN.size(), 0.);
					//	//cout << "sum of m1:\n";
					//	//for (j = 0; j < hmesh[lv0][eid0].IEN.size(); j++)
					//	//{
					//	//	for (k = 0; k < hmesh[lv][prt].IEN.size(); k++)
					//	//	{
					//	//		sum1[j] += m1(k, j);
					//	//	}
					//	//	cout << sum1[j] << "\n";
					//	//}
					//	getchar();
					//}
					MatrixXd mtmp = m1*m0;
					m0.resize(hmesh[lv][prt].IEN.size(), hmesh[lev][i].IEN.size());
					m0 = mtmp;
					//cout
					//if (lev == 2 && lv == 0)
					//{
					//	vector<double> sum0(hmesh[lev][i].IEN.size(), 0.);
					//	cout << "sum of m0 after:\n";
					//	for (j = 0; j < hmesh[lev][i].IEN.size(); j++)
					//	{
					//		for (k = 0; k < hmesh[lv][prt].IEN.size(); k++)
					//		{
					//			sum0[j] += m0(k, j);
					//		}
					//		cout << sum0[j] << "\n";
					//	}
					//	getchar();
					//}
					for (j = 0; j < hmesh[lv][prt].IEN.size(); j++)
					{
						if (hcp[lv][hmesh[lv][prt].IEN[j]].act == 1)
						{
							array<int, 2> btmp = { lv, hmesh[lv][prt].IEN[j] };
							vector<double> ctmp(hmesh[lev][i].IEN.size());
							for (k = 0; k < hmesh[lev][i].IEN.size(); k++)
							{
								ctmp[k] = m0(j,k);
							}
							hmesh[lev][i].IEN_act.push_back(btmp);
							hmesh[lev][i].tmat.push_back(ctmp);
						}
					}
					eid0 = prt;
					lv0 = lv;
					prt = hmesh[lv][prt].prt;
					lv--;
				}
				//cout
				//if (lev == 2)
				//{
				//	vector<double> sum(hmesh[lev][i].IEN.size(), 0.);
				//	cout << "sum of trun ele: " << hmesh[lev][i].type <<"\n";
				//	for (j = 0; j < hmesh[lev][i].IEN.size(); j++)
				//	{
				//		for (k = 0; k < hmesh[lev][i].IEN_act.size(); k++)
				//		{
				//			sum[j] += hmesh[lev][i].tmat[k][j];
				//		}
				//		cout << sum[j] << "\n";
				//	}
				//	getchar();
				//}
			}
		}
	}

	//for (i = 0; i < hmesh[lev].size(); i++)
	//{
	//	if (hmesh[lev][i].trun == 1)
	//	{
	//		for (j = 0; j < hmesh[lev][i].IEN_act.size(); j++)
	//		{
	//			int idtmp[2] = { hmesh[lev][i].IEN_act[j][0], hmesh[lev][i].IEN_act[j][1] };
	//			if (hcp[idtmp[0]][idtmp[1]].act == 0)
	//			{
	//				cout << "Inactive basis function!\n";
	//				cout << idtmp[0] << " " << idtmp[1] << "\n";
	//				getchar();
	//			}
	//		}
	//	}
	//}
}

void TruncatedTspline_3D::Basis_Regular(int lev, int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt)
{
	Nt.clear();
	dNdt.clear();
	if (hmesh[lev][eid].trun == 0)
	{
		Nt.resize(hmesh[lev][eid].IEN.size());
		dNdt.resize(hmesh[lev][eid].IEN.size());
		vector<double> ku(5, 0.), kv(5, 0.), kw(5, 0.), uval, vval, wval;
		BSplineBasis bu, bv, bw;
		for (uint i = 0; i<hmesh[lev][eid].IEN.size(); i++)
		{
			ku.assign(hmesh[lev][eid].patch_ku[i].begin(), hmesh[lev][eid].patch_ku[i].end());
			kv.assign(hmesh[lev][eid].patch_kv[i].begin(), hmesh[lev][eid].patch_kv[i].end());
			kw.assign(hmesh[lev][eid].patch_kw[i].begin(), hmesh[lev][eid].patch_kw[i].end());
			bu.Set(3, ku);
			bv.Set(3, kv);
			bw.Set(3, kw);
			bu.BasisFunction(0, u[0], 1, uval);
			bv.BasisFunction(0, u[1], 1, vval);
			bw.BasisFunction(0, u[2], 1, wval);
			Nt[i] = uval[0] * vval[0] * wval[0];
			dNdt[i][0] = uval[1] * vval[0] * wval[0];
			dNdt[i][1] = uval[0] * vval[1] * wval[0];
			dNdt[i][2] = uval[0] * vval[0] * wval[1];
		}
	}
	else
	{
		Nt.resize(hmesh[lev][eid].IEN_act.size());
		dNdt.resize(hmesh[lev][eid].IEN_act.size());
		vector<double> Nt0(hmesh[lev][eid].IEN.size());
		vector<array<double, 3>> dNdt0(hmesh[lev][eid].IEN.size());
		vector<double> ku(5, 0.), kv(5, 0.), kw(5, 0.), uval, vval, wval;
		BSplineBasis bu, bv, bw;
		for (uint i = 0; i<hmesh[lev][eid].IEN.size(); i++)
		{
			ku.assign(hmesh[lev][eid].patch_ku[i].begin(), hmesh[lev][eid].patch_ku[i].end());
			kv.assign(hmesh[lev][eid].patch_kv[i].begin(), hmesh[lev][eid].patch_kv[i].end());
			kw.assign(hmesh[lev][eid].patch_kw[i].begin(), hmesh[lev][eid].patch_kw[i].end());
			bu.Set(3, ku);
			bv.Set(3, kv);
			bw.Set(3, kw);
			bu.BasisFunction(0, u[0], 1, uval);
			bv.BasisFunction(0, u[1], 1, vval);
			bw.BasisFunction(0, u[2], 1, wval);
			Nt0[i] = uval[0] * vval[0] * wval[0];
			dNdt0[i][0] = uval[1] * vval[0] * wval[0];
			dNdt0[i][1] = uval[0] * vval[1] * wval[0];
			dNdt0[i][2] = uval[0] * vval[0] * wval[1];
		}
		for (uint i = 0; i < hmesh[lev][eid].IEN_act.size(); i++)
		{
			Nt[i] = 0.; dNdt[i][0] = 0.; dNdt[i][1] = 0.; dNdt[i][2] = 0.;
			for (uint j = 0; j < hmesh[lev][eid].IEN.size(); j++)
			{
				if (hmesh[lev][eid].tmat[i][j] != 0.)
				{
					Nt[i] += hmesh[lev][eid].tmat[i][j] * Nt0[j];
					dNdt[i][0] += hmesh[lev][eid].tmat[i][j] * dNdt0[j][0];
					dNdt[i][1] += hmesh[lev][eid].tmat[i][j] * dNdt0[j][1];
					dNdt[i][2] += hmesh[lev][eid].tmat[i][j] * dNdt0[j][2];
				}
			}
		}
	}
}

void TruncatedTspline_3D::Basis_Irregular(int lev, int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt)
{
	Nt.clear();
	dNdt.clear();
	if (hmesh[lev][eid].trun == 0)
	{
		Nt.resize(hmesh[lev][eid].IEN.size());
		dNdt.resize(hmesh[lev][eid].IEN.size());
		BezierElement3D bzel;
		vector<double> Bt;
		vector<array<double, 3>> dBdt;
		bzel.Basis(u[0], u[1], u[2], Bt, dBdt);
		for (uint i = 0; i < hmesh[lev][eid].IEN.size(); i++)
		{
			Nt[i] = 0.;
			dNdt[i][0] = 0.; dNdt[i][1] = 0.; dNdt[i][2] = 0.;
			for (uint j = 0; j < 64; j++)
			{
				Nt[i] += hmesh[lev][eid].bemat[i][j] * Bt[j];
				dNdt[i][0] += hmesh[lev][eid].bemat[i][j] * dBdt[j][0];
				dNdt[i][1] += hmesh[lev][eid].bemat[i][j] * dBdt[j][1];
				dNdt[i][2] += hmesh[lev][eid].bemat[i][j] * dBdt[j][2];
			}
		}
	}
	else
	{
		Nt.resize(hmesh[lev][eid].IEN_act.size());
		dNdt.resize(hmesh[lev][eid].IEN_act.size());
		vector<double> Nt0(hmesh[lev][eid].IEN.size());
		vector<array<double, 3>> dNdt0(hmesh[lev][eid].IEN.size());
		BezierElement3D bzel;
		vector<double> Bt;
		vector<array<double, 3>> dBdt;
		bzel.Basis(u[0], u[1], u[2], Bt, dBdt);
		for (uint i = 0; i < hmesh[lev][eid].IEN.size(); i++)
		{
			Nt0[i] = 0.;
			dNdt0[i][0] = 0.; dNdt0[i][1] = 0.; dNdt0[i][2] = 0.;
			for (uint j = 0; j < 64; j++)
			{
				Nt0[i] += hmesh[lev][eid].bemat[i][j] * Bt[j];
				dNdt0[i][0] += hmesh[lev][eid].bemat[i][j] * dBdt[j][0];
				dNdt0[i][1] += hmesh[lev][eid].bemat[i][j] * dBdt[j][1];
				dNdt0[i][2] += hmesh[lev][eid].bemat[i][j] * dBdt[j][2];
			}
		}
		for (uint i = 0; i < hmesh[lev][eid].IEN_act.size(); i++)
		{
			Nt[i] = 0.; dNdt[i][0] = 0.; dNdt[i][1] = 0.; dNdt[i][2] = 0.;
			for (uint j = 0; j < hmesh[lev][eid].IEN.size(); j++)
			{
				if (hmesh[lev][eid].tmat[i][j] != 0.)
				{
					Nt[i] += hmesh[lev][eid].tmat[i][j] * Nt0[j];
					dNdt[i][0] += hmesh[lev][eid].tmat[i][j] * dNdt0[j][0];
					dNdt[i][1] += hmesh[lev][eid].tmat[i][j] * dNdt0[j][1];
					dNdt[i][2] += hmesh[lev][eid].tmat[i][j] * dNdt0[j][2];
				}
			}
		}
	}
}

void TruncatedTspline_3D::Identify_Pseudo(vector<array<int, 2>>& rfid)
{
	//int nh(9);
	//for (int i = 1; i < nh - 1; i++)
	//{
	//	array<int, 2> tmp = { 0, i*nh*nh + i*nh + i };
	//	rfid.push_back(tmp);
	//}

	rfid.clear();
	int nh(9+1), lev(0);
	vector<int> p_rf;
	for (int i = 2; i < nh - 2; i++)
	{
		p_rf.push_back(i*nh*nh + i*nh + i);
	}
	for (uint i = 0; i < p_rf.size(); i++)
	{
		for (uint j = 0; j < hcp[lev][p_rf[i]].hex.size(); j++)
		{
			array<int, 2> tmp = { lev, hcp[lev][p_rf[i]].hex[j] };
			vector<array<int, 2>>::iterator it = find(rfid.begin(),rfid.end(),tmp);
			if (it == rfid.end()) rfid.push_back(tmp);
		}
	}
	//cout <<"#rfid: "<< rfid.size() << "\n";
	//for (uint i = 0; i < rfid.size(); i++)
	//{
	//	cout << rfid[i][1] << " ";
	//}
	//cout << "\n";
	//getchar();
}

void TruncatedTspline_3D::Identify_Pseudo_1(vector<array<int, 2>>& rfid)
{
	rfid.clear();
	int lev(hmesh.size()-1);
	vector<int> p_rf;
	for (uint i = 0; i < hcp[lev].size(); i++)
	{
		if (hcp[lev][i].act == 1) p_rf.push_back(i);
	}
	for (uint i = 0; i < p_rf.size(); i++)
	{
		for (uint j = 0; j < hcp[lev][p_rf[i]].hex.size(); j++)
		{
			array<int, 2> tmp = { lev, hcp[lev][p_rf[i]].hex[j] };
			vector<array<int, 2>>::iterator it = find(rfid.begin(), rfid.end(), tmp);
			if (it == rfid.end()) rfid.push_back(tmp);
		}
	}
}

void TruncatedTspline_3D::Identify_Test_1(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst)
{
	rfid.clear();
	gst.clear();
	int lev(0);
	int pid(-1);
	for (uint i = 0; i < hmesh[lev].size(); i++)
	{
		if (hmesh[lev][i].type == 2)
		{
			for (int j = 0; j < 8; j++)
			{
				if (hcp[lev][hmesh[lev][i].cnct[j]].type == 0)
				{
					pid = hmesh[lev][i].cnct[j]; break;
				}
			}
		}
		if (pid != -1) break;
	}
	for (uint i = 0; i < hcp[lev][pid].hex.size(); i++)
	{
		array<int, 2> tmp = { lev, hcp[lev][pid].hex[i] };
		rfid.push_back(tmp);
	}
	for (uint i = 0; i < rfid.size(); i++)
	{
		for (int j = 0; j < 8; j++)
		{
			int pt(hmesh[lev][rfid[i][1]].cnct[j]);
			for (uint k = 0; k < hcp[lev][pt].hex.size(); k++)
			{
				array<int, 2> tmp = { lev, hcp[lev][pt].hex[k] };
				vector<array<int, 2>>::iterator it1 = find(gst.begin(),gst.end(),tmp);
				vector<array<int, 2>>::iterator it2 = find(rfid.begin(), rfid.end(), tmp);
				if (it1 == gst.end() && it2 == rfid.end())
				{
					gst.push_back(tmp);
				}
			}
		}
	}
}

void TruncatedTspline_3D::Identify_Test_2(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst)
{
	rfid.clear();
	gst.clear();
	int lev(0);
	int pid(-1);
	for (uint i = 0; i < hmesh[lev].size(); i++)
	{
		if (hmesh[lev][i].type == 2)
		{
			for (int j = 0; j < 8; j++)
			{
				if (hcp[lev][hmesh[lev][i].cnct[j]].type == 3)
				{
					pid = hmesh[lev][i].cnct[j]; break;
				}
			}
		}
		if (pid != -1) break;
	}
	for (uint i = 0; i < hcp[lev][pid].hex.size(); i++)
	{
		array<int, 2> tmp = { lev, hcp[lev][pid].hex[i] };
		rfid.push_back(tmp);
	}
	for (uint i = 0; i < rfid.size(); i++)
	{
		for (int j = 0; j < 8; j++)
		{
			int pt(hmesh[lev][rfid[i][1]].cnct[j]);
			for (uint k = 0; k < hcp[lev][pt].hex.size(); k++)
			{
				array<int, 2> tmp = { lev, hcp[lev][pt].hex[k] };
				vector<array<int, 2>>::iterator it1 = find(gst.begin(), gst.end(), tmp);
				vector<array<int, 2>>::iterator it2 = find(rfid.begin(), rfid.end(), tmp);
				if (it1 == gst.end() && it2 == rfid.end())
				{
					gst.push_back(tmp);
				}
			}
		}
	}
}

void TruncatedTspline_3D::Identify_Test_3(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst)
{
	rfid.clear();
	gst.clear();
	int lev(hcp.size()-1);
	int pid(-1);
	for (uint i = 0; i < hcp[lev].size(); i++)
	{
		if (hcp[lev][i].act == 1)
		{
			pid = i; break;
		}
	}
	if (pid == -1)
	{
		cerr << "Cannot find active points at level " << lev << "!\n";
		getchar();
	}
	for (uint i = 0; i < hcp[lev][pid].hex.size(); i++)
	{
		array<int, 2> tmp = { lev, hcp[lev][pid].hex[i] };
		rfid.push_back(tmp);
	}
	for (uint i = 0; i < rfid.size(); i++)
	{
		for (int j = 0; j < 8; j++)
		{
			int pt(hmesh[lev][rfid[i][1]].cnct[j]);
			for (uint k = 0; k < hcp[lev][pt].hex.size(); k++)
			{
				array<int, 2> tmp = { lev, hcp[lev][pt].hex[k] };
				vector<array<int, 2>>::iterator it1 = find(gst.begin(), gst.end(), tmp);
				vector<array<int, 2>>::iterator it2 = find(rfid.begin(), rfid.end(), tmp);
				if (it1 == gst.end() && it2 == rfid.end())
				{
					gst.push_back(tmp);
				}
			}
		}
	}
}

void TruncatedTspline_3D::Identify_Test_4(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst)
{
	rfid.clear();
	gst.clear();
	int lev(0);
	//int pid(3460);
	int pid(2792);
	//for (uint i = 0; i < hcp[lev].size(); i++)
	//{
	//	if (hcp[lev][i].type == 1)
	//	{
	//		pid=
	//		for (int j = 0; j < 8; j++)
	//		{
	//			if (hcp[lev][hmesh[lev][i].cnct[j]].type == 3)
	//			{
	//				pid = hmesh[lev][i].cnct[j]; break;
	//			}
	//		}
	//	}
	//	if (pid != -1) break;
	//}
	for (uint i = 0; i < hcp[lev][pid].hex.size(); i++)
	{
		array<int, 2> tmp = { lev, hcp[lev][pid].hex[i] };
		rfid.push_back(tmp);
	}
	for (uint i = 0; i < rfid.size(); i++)
	{
		for (int j = 0; j < 8; j++)
		{
			int pt(hmesh[lev][rfid[i][1]].cnct[j]);
			for (uint k = 0; k < hcp[lev][pt].hex.size(); k++)
			{
				array<int, 2> tmp = { lev, hcp[lev][pt].hex[k] };
				vector<array<int, 2>>::iterator it1 = find(gst.begin(), gst.end(), tmp);
				vector<array<int, 2>>::iterator it2 = find(rfid.begin(), rfid.end(), tmp);
				if (it1 == gst.end() && it2 == rfid.end())
				{
					gst.push_back(tmp);
				}
			}
		}
	}
}

//void TruncatedTspline_3D::Refine(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst)
//{
//	for (uint i = 0; i < rfid.size(); i++)
//	{
//		if (hmesh[rfid[i][0]][rfid[i][1]].type == 0)
//		{
//			PatchRefine_Regular(rfid[i][0], rfid[i][1]);
//		}
//		else if (hmesh[rfid[i][0]][rfid[i][1]].type == 1)
//		{
//			PatchRefine_Boundary(rfid[i][0], rfid[i][1]);
//		}
//		else if (hmesh[rfid[i][0]][rfid[i][1]].type == 2)
//		{
//			PatchRefine_Irregular(rfid[i][0], rfid[i][1]);
//		}
//	}
//	Refine_Ghost(gst);
//
//	int lev(rfid[0][0] + 1);//tmp
//	for (uint i = 0; i < hmesh[lev].size(); i++)
//	{
//		ConstructFaceEdge(lev, i);
//	}
//	ConstructConnect(lev);
//	for (uint i = 0; i < hmesh[lev].size(); i++)
//	{
//		if (hmesh[lev][i].ghost == 0)
//		{
//			if (hmesh[lev][i].type == 2)
//			{
//				ConstructBezierBasis(lev, i);
//			}
//			else if (hmesh[lev][i].type == 1)
//			{
//				ConstructBezierBasis_Boundary(lev, i);
//			}
//		}
//	}
//	Selection(lev);
//	Truncate(lev);
//}

void TruncatedTspline_3D::Refine(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst)
{
	vector<MatrixXd> bsmat;
	BezierSubdivMatrix(bsmat);

	cout << "# of to-be-refined elements: " << rfid.size() << "\n";
	for (uint i = 0; i < rfid.size(); i++)
	{
		if (i!=0 && i % 200 == 0) cout << i << " ";
		//cout << "Refine ID: " << i << "/" << rfid.size()<<"\n";
		if (hmesh[rfid[i][0]][rfid[i][1]].type == 0 || hmesh[rfid[i][0]][rfid[i][1]].type == 2)
		{
			if (hmesh[rfid[i][0]][rfid[i][1]].chd.size() == 0)
			{
				PatchRefine_Irregular(rfid[i][0], rfid[i][1]);
			}
			else//previously refined to ghost children elements
			{
				hmesh[rfid[i][0]][rfid[i][1]].act = 0;
				for (uint j = 0; j < hmesh[rfid[i][0]][rfid[i][1]].chd.size(); j++)
				{
					hmesh[rfid[i][0] + 1][hmesh[rfid[i][0]][rfid[i][1]].chd[j]].act = 1;
					hmesh[rfid[i][0] + 1][hmesh[rfid[i][0]][rfid[i][1]].chd[j]].ghost = 0;
				}
			}
		}
		else if (hmesh[rfid[i][0]][rfid[i][1]].type == 1)
		{
			if (hmesh[rfid[i][0]][rfid[i][1]].chd.size() == 0)
			{
				PatchRefine_Boundary(rfid[i][0], rfid[i][1]);
			}
			else
			{
				hmesh[rfid[i][0]][rfid[i][1]].act = 0;
				for (uint j = 0; j < hmesh[rfid[i][0]][rfid[i][1]].chd.size(); j++)
				{
					hmesh[rfid[i][0] + 1][hmesh[rfid[i][0]][rfid[i][1]].chd[j]].act = 1;
					hmesh[rfid[i][0] + 1][hmesh[rfid[i][0]][rfid[i][1]].chd[j]].ghost = 0;
				}
			}
		}

		BezierPoints_Refine(rfid[i][0],rfid[i][1],bsmat);

		//else if (hmesh[rfid[i][0]][rfid[i][1]].type == 2)
		//{
		//	PatchRefine_Irregular(rfid[i][0], rfid[i][1]);
		//}
	}
	cout << "# of to-be-refined ghost: " << gst.size() << "\n";
	Refine_Ghost(gst);

	//int lev(rfid[0][0] + 1);//tmp
	vector<int> lev;
	for (uint i = 0; i < rfid.size(); i++)
	{
		for (uint j = 0; j < hmesh[rfid[i][0]][rfid[i][1]].chd.size(); j++)
		{
			int eid[2] = { rfid[i][0] + 1, hmesh[rfid[i][0]][rfid[i][1]].chd[j] };
			ConstructFaceEdge(eid[0], eid[1]);
		}
		vector<int>::iterator it = find(lev.begin(), lev.end(), rfid[i][0] + 1);
		if (it == lev.end()) lev.push_back(rfid[i][0] + 1);
	}
	for (uint i = 0; i < gst.size(); i++)
	{
		for (uint j = 0; j < hmesh[gst[i][0]][gst[i][1]].chd.size(); j++)
		{
			int eid[2] = { gst[i][0] + 1, hmesh[gst[i][0]][gst[i][1]].chd[j] };
			ConstructFaceEdge(eid[0], eid[1]);
		}
	}
	//for (uint i = 0; i < hmesh[lev].size(); i++)
	//{
	//	ConstructFaceEdge(lev, i);
	//}
	for (uint i = 0; i < lev.size(); i++)
	{
		ConstructConnect(lev[i]);
	}
	//construct basis functions for new high-level elements
	for (uint i = 0; i < rfid.size(); i++)
	{
		for (uint j = 0; j < hmesh[rfid[i][0]][rfid[i][1]].chd.size(); j++)
		{
			int eid[2] = { rfid[i][0] + 1, hmesh[rfid[i][0]][rfid[i][1]].chd[j] };
			if (hmesh[eid[0]][eid[1]].ghost == 0 && hmesh[eid[0]][eid[1]].IEN.size()==0)
			{
				if (hmesh[eid[0]][eid[1]].type == 0 || hmesh[eid[0]][eid[1]].type == 2)
				{
					ConstructBezierBasis(eid[0], eid[1]);
				}
				else if (hmesh[eid[0]][eid[1]].type == 1)
				{
					ConstructBezierBasis_Boundary(eid[0], eid[1]);
				}
			}
		}
	}
	//for (uint i = 0; i < hmesh[lev].size(); i++)
	//{
	//	if (hmesh[lev][i].ghost == 0)
	//	{
	//		if (hmesh[lev][i].type == 0 || hmesh[lev][i].type == 2)
	//		{
	//			ConstructBezierBasis(lev, i);
	//		}
	//		else if (hmesh[lev][i].type == 1)
	//		{
	//			ConstructBezierBasis_Boundary(lev, i);
	//		}
	//	}
	//}

	for (uint i = 0; i < lev.size(); i++)
	{
		SetSupport(lev[i]);
	}
	Select();

	for (uint i = 1; i < hmesh.size(); i++)
	{
		Truncate(i);
	}

	//Selection(lev);
	//SetSupport(lev);
	//Select();
	//Truncate(lev);
}

void TruncatedTspline_3D::OutputCM(int lev, string fn)
{
	string fname(fn + "_CM.vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << hcp[lev].size() << " float\n";
		for (uint i = 0; i<hcp[lev].size(); i++)
		{
			fout << hcp[lev][i].coor[0] << " " << hcp[lev][i].coor[1] << " " << hcp[lev][i].coor[2] << "\n";
		}
		fout << "\nCELLS " << hmesh[lev].size() << " " << 9 * hmesh[lev].size() << '\n';
		for (uint i = 0; i<hmesh[lev].size(); i++)
		{
			fout << "8 ";
			for (int j = 0; j<8; j++)
			{
				fout << hmesh[lev][i].cnct[j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << hmesh[lev].size() << '\n';
		for (uint i = 0; i<hmesh[lev].size(); i++)
		{
			fout << "12\n";
		}
		
		//fout << "POINT_DATA " << hcp[lev].size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<hcp[lev].size(); i++)
		//{
		//	fout << hcp[lev][i].act << "\n";
		//	//if (hcp[lev][i].act == 1 && hcp[lev][i].type == 1)
		//	//{
		//	//	fout << "1\n";
		//	//}
		//	//else
		//	//{
		//	//	fout << "0\n";
		//	//}
		//}

		fout<<"\nCELL_DATA "<<hmesh[lev].size()<<"\nSCALARS eact float 1\nLOOKUP_TABLE default\n";
		for(uint i=0;i<hmesh[lev].size();i++)
		{
			//fout<<hmesh[lev][i].act<<"\n";
			fout << hmesh[lev][i].type << "\n";
		}

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline_3D::OutputFace(int lev, string fn)
{
	string fname(fn + "_face.vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << hcp[lev].size() << " float\n";
		for (uint i = 0; i<hcp[lev].size(); i++)
		{
			fout << hcp[lev][i].coor[0] << " " << hcp[lev][i].coor[1] << " " << hcp[lev][i].coor[2] << "\n";
		}
		fout << "\nCELLS " << hface[lev].size() << " " << 5 * hface[lev].size() << '\n';
		for (uint i = 0; i<hface[lev].size(); i++)
		{
			fout << "4 ";
			for (int j = 0; j<4; j++)
			{
				fout << hface[lev][i].cnct[j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << hface[lev].size() << '\n';
		for (uint i = 0; i<hface[lev].size(); i++)
		{
			fout << "9\n";
		}

		//fout << "POINT_DATA " << hcp[lev].size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<hcp[lev].size(); i++)
		//{
		//	fout << hcp[lev][i].act << "\n";
		//}

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline_3D::OutputEdge(int lev, string fn)
{
	string fname(fn + "_edge.vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << hcp[lev].size() << " float\n";
		for (uint i = 0; i<hcp[lev].size(); i++)
		{
			fout << hcp[lev][i].coor[0] << " " << hcp[lev][i].coor[1] << " " << hcp[lev][i].coor[2] << "\n";
		}
		fout << "\nCELLS " << hedge[lev].size() << " " << 3 * hedge[lev].size() << '\n';
		for (uint i = 0; i<hedge[lev].size(); i++)
		{
			fout << "2 ";
			for (int j = 0; j<2; j++)
			{
				fout << hedge[lev][i].pt[j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << hedge[lev].size() << '\n';
		for (uint i = 0; i<hedge[lev].size(); i++)
		{
			fout << "3\n";
		}

		//fout << "POINT_DATA " << hcp[lev].size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<hcp[lev].size(); i++)
		//{
		//	fout << hcp[lev][i].sharp << "\n";
		//	//fout << hcp[lev][i].bcxp << "\n";
		//}
		fout << "\nCELL_DATA " << hedge[lev].size() << "\nSCALARS eact float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<hedge[lev].size(); i++)
		{
			fout << hedge[lev][i].sharp << "\n";
			//fout << hedge[lev][i].type << "\n";
		}

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline_3D::OutputGeom(int lev, string fn)
{
	vector<array<double, 3>> spt;
	vector<array<double, 3>> sval;
	vector<double> ssum;
	vector<array<int, 8>> sele;
	vector<array<double, 3>> lpt;//visulize parameter lines
	vector<array<int, 2>> led;//line connectivity
	int ns(2), ne_ref(0), loc0, loc1, loc2;

	//for (uint lev = 0; lev < hmesh.size(); lev++)
	{
		for (uint eid = 0; eid < hmesh[lev].size(); eid++)
		{
			if (hmesh[lev][eid].act == 1 /*&& (hmesh[lev][eid].type == 0 || hmesh[lev][eid].type == 2)*/)
			{
				vector<double> su(ns), sv(ns), sw(ns);
				double ul[3] = { hmesh[lev][eid].dm[0][1] - hmesh[lev][eid].dm[0][0], hmesh[lev][eid].dm[1][1] - hmesh[lev][eid].dm[1][0], 
					hmesh[lev][eid].dm[2][1] - hmesh[lev][eid].dm[2][0]};
				for (int i = 0; i<ns; i++)
				{
					su[i] = i*ul[0] / (ns - 1) + hmesh[lev][eid].dm[0][0];
					sv[i] = i*ul[1] / (ns - 1) + hmesh[lev][eid].dm[1][0];
					sw[i] = i*ul[2] / (ns - 1) + hmesh[lev][eid].dm[2][0];
				}
				for (int a = 0; a<ns; a++)
				{
					for (int b = 0; b<ns; b++)
					{
						for (int c = 0; c < ns; c++)
						{
							array<double, 3> pt;
							array<double,3> uval = { su[c], sv[b], sw[a] };
							GeomMap(lev, eid, uval, pt);
							double sumtmp = BasisSum(lev, eid, uval);
							//GeomMap_Lev(lev, eid, uval, pt);
							//double sumtmp = BasisSum_Lev(lev, eid, uval);
							spt.push_back(pt);
							ssum.push_back(sumtmp);
							//if(a==0||a==ns-1||b==0||b==ns-1)
							//{
							//	lpt.push_back(pt);
							//}
						}
					}
				}

				for (int a = 0; a<ns - 1; a++)
				{
					for (int b = 0; b<ns - 1; b++)
					{
						for (int c = 0; c < ns - 1; c++)
						{
							array<int, 8> el;
							el[0] = ne_ref + a*ns*ns + b*ns + c;
							el[1] = ne_ref + a*ns*ns + b*ns + c + 1;
							el[2] = ne_ref + a*ns*ns + (b + 1)*ns + c + 1;
							el[3] = ne_ref + a*ns*ns + (b + 1)*ns + c;
							el[4] = ne_ref + (a + 1)*ns*ns + b*ns + c;
							el[5] = ne_ref + (a + 1)*ns*ns + b*ns + c + 1;
							el[6] = ne_ref + (a + 1)*ns*ns + (b + 1)*ns + c + 1;
							el[7] = ne_ref + (a + 1)*ns*ns + (b + 1)*ns + c;
							sele.push_back(el);
						}
					}
				}
				ne_ref += ns*ns*ns;
			}
		}
	}

	string fname = fn + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (uint i = 0; i<sele.size(); i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i<sele.size(); i++)
		{
			fout << "12\n";
		}
		fout << "\nPOINT_DATA " << ssum.size() << "\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<ssum.size(); i++)
		{
			fout << ssum[i] << "\n";
		}
		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nNORMALS Normal FLOAT\n";
		//for(uint i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i][0]<<" "<<sval[i][1]<<" "<<sval[i][2]<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline_3D::GeomMap(int lev, int eid, const array<double, 3>& u, array<double, 3>& pt)
{
	vector<double> Nt;
	vector<array<double, 3>> dNdt;
	//if (hmesh[lev][eid].type == 0)
	//{
	//	Basis_Regular(lev, eid, u, Nt, dNdt);
	//}
	//else if (hmesh[lev][eid].type == 2 || hmesh[lev][eid].type == 1)
	{
		Basis_Irregular(lev, eid, u, Nt, dNdt);
	}
	pt[0] = 0.; pt[1] = 0.; pt[2] = 0.;
	if (hmesh[lev][eid].trun == 1)
	{
		for (uint i = 0; i < hmesh[lev][eid].IEN_act.size(); i++)
		{
			int lid(hmesh[lev][eid].IEN_act[i][0]);
			int pid(hmesh[lev][eid].IEN_act[i][1]);
			pt[0] += Nt[i] * hcp[lid][pid].coor[0];
			pt[1] += Nt[i] * hcp[lid][pid].coor[1];
			pt[2] += Nt[i] * hcp[lid][pid].coor[2];
		}
	}
	else
	{
		for (uint i = 0; i < hmesh[lev][eid].IEN.size(); i++)
		{
			int pid(hmesh[lev][eid].IEN[i]);
			pt[0] += Nt[i] * hcp[lev][pid].coor[0];
			pt[1] += Nt[i] * hcp[lev][pid].coor[1];
			pt[2] += Nt[i] * hcp[lev][pid].coor[2];
		}
	}
}

void TruncatedTspline_3D::GeomMap_Bezier(int lev, int eid, const array<double, 3>& u, array<double, 3>& pt)
{
	BezierElement3D bzel;
	vector<double> Bt;
	vector<array<double, 3>> dBdt;
	bzel.Basis(u[0], u[1], u[2], Bt, dBdt);
	pt[0] = 0.; pt[1] = 0.; pt[2] = 0.;
	for (int i = 0; i < 64; i++)
	{
		pt[0] += Bt[i] * hmesh[lev][eid].bzpt[i][0];
		pt[1] += Bt[i] * hmesh[lev][eid].bzpt[i][1];
		pt[2] += Bt[i] * hmesh[lev][eid].bzpt[i][2];
	}
}

void TruncatedTspline_3D::OutputGeom_All(string fn)
{
	vector<array<double, 3>> spt;
	vector<array<double, 3>> sval;
	vector<double> ssum;
	vector<array<int, 8>> sele;
	vector<array<double, 3>> lpt;//visulize parameter lines
	vector<array<int, 2>> led;//line connectivity
	int ns(2), ne_ref(0), loc0, loc1, loc2;

	for (uint lev = 0; lev < hmesh.size(); lev++)
	{
		for (uint eid = 0; eid < hmesh[lev].size(); eid++)
		{
			if (hmesh[lev][eid].act == 1 /*&& (hmesh[lev][eid].type == 1 || hmesh[lev][eid].type == 2)*/)
			{
				vector<double> su(ns), sv(ns), sw(ns);
				double ul[3] = { hmesh[lev][eid].dm[0][1] - hmesh[lev][eid].dm[0][0], hmesh[lev][eid].dm[1][1] - hmesh[lev][eid].dm[1][0],
					hmesh[lev][eid].dm[2][1] - hmesh[lev][eid].dm[2][0] };
				for (int i = 0; i<ns; i++)
				{
					su[i] = i*ul[0] / (ns - 1) + hmesh[lev][eid].dm[0][0];//maybe not useful
					sv[i] = i*ul[1] / (ns - 1) + hmesh[lev][eid].dm[1][0];
					sw[i] = i*ul[2] / (ns - 1) + hmesh[lev][eid].dm[2][0];
				}
				for (int a = 0; a<ns; a++)
				{
					for (int b = 0; b<ns; b++)
					{
						for (int c = 0; c < ns; c++)
						{
							array<double, 3> pt;
							array<double, 3> uval = { su[c], sv[b], sw[a] };
							//GeomMap(lev, eid, uval, pt);
							GeomMap_Bezier(lev, eid, uval, pt);
							double sumtmp = BasisSum(lev, eid, uval);
							//GeomMap_Lev(lev, eid, uval, pt);
							//double sumtmp = BasisSum_Lev(lev, eid, uval);
							spt.push_back(pt);
							ssum.push_back(sumtmp);
							//if(a==0||a==ns-1||b==0||b==ns-1)
							//{
							//	lpt.push_back(pt);
							//}
						}
					}
				}

				for (int a = 0; a<ns - 1; a++)
				{
					for (int b = 0; b<ns - 1; b++)
					{
						for (int c = 0; c < ns - 1; c++)
						{
							array<int, 8> el;
							el[0] = ne_ref + a*ns*ns + b*ns + c;
							el[1] = ne_ref + a*ns*ns + b*ns + c + 1;
							el[2] = ne_ref + a*ns*ns + (b + 1)*ns + c + 1;
							el[3] = ne_ref + a*ns*ns + (b + 1)*ns + c;
							el[4] = ne_ref + (a + 1)*ns*ns + b*ns + c;
							el[5] = ne_ref + (a + 1)*ns*ns + b*ns + c + 1;
							el[6] = ne_ref + (a + 1)*ns*ns + (b + 1)*ns + c + 1;
							el[7] = ne_ref + (a + 1)*ns*ns + (b + 1)*ns + c;
							sele.push_back(el);
						}
					}
				}
				ne_ref += ns*ns*ns;
			}
		}
	}

	string fname = fn + "_geom.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (uint i = 0; i<sele.size(); i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i<sele.size(); i++)
		{
			fout << "12\n";
		}
		fout << "\nPOINT_DATA " << ssum.size() << "\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<ssum.size(); i++)
		{
			fout << ssum[i] << "\n";
		}
		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nNORMALS Normal FLOAT\n";
		//for(uint i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i][0]<<" "<<sval[i][1]<<" "<<sval[i][2]<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

double TruncatedTspline_3D::BasisSum(int lev, int eid, const array<double, 3>& u)
{
	vector<double> Nt;
	vector<array<double, 3>> dNdt;
	//if (hmesh[lev][eid].type == 0)
	//{
	//	Basis_Regular(lev, eid, u, Nt, dNdt);
	//}
	//else if (hmesh[lev][eid].type == 2 || hmesh[lev][eid].type == 1)
	{
		Basis_Irregular(lev, eid, u, Nt, dNdt);
	}
	double sum(0.);
	for (uint i = 0; i < Nt.size(); i++) sum += Nt[i];
	return sum;
}

void TruncatedTspline_3D::Basis_Regular_Lev(int lev, int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt)
{
	Nt.clear();
	dNdt.clear();
	Nt.resize(hmesh[lev][eid].IEN.size());
	dNdt.resize(hmesh[lev][eid].IEN.size());
	vector<double> ku(5, 0.), kv(5, 0.), kw(5, 0.), uval, vval, wval;
	BSplineBasis bu, bv, bw;
	for (uint i = 0; i < hmesh[lev][eid].IEN.size(); i++)
	{
		ku.assign(hmesh[lev][eid].patch_ku[i].begin(), hmesh[lev][eid].patch_ku[i].end());
		kv.assign(hmesh[lev][eid].patch_kv[i].begin(), hmesh[lev][eid].patch_kv[i].end());
		kw.assign(hmesh[lev][eid].patch_kw[i].begin(), hmesh[lev][eid].patch_kw[i].end());
		bu.Set(3, ku);
		bv.Set(3, kv);
		bw.Set(3, kw);
		bu.BasisFunction(0, u[0], 1, uval);
		bv.BasisFunction(0, u[1], 1, vval);
		bw.BasisFunction(0, u[2], 1, wval);
		Nt[i] = uval[0] * vval[0] * wval[0];
		dNdt[i][0] = uval[1] * vval[0] * wval[0];
		dNdt[i][1] = uval[0] * vval[1] * wval[0];
		dNdt[i][2] = uval[0] * vval[0] * wval[1];
	}
}

void TruncatedTspline_3D::Basis_Irregular_Lev(int lev, int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt)
{
	Nt.clear();
	dNdt.clear();
	Nt.resize(hmesh[lev][eid].IEN.size());
	dNdt.resize(hmesh[lev][eid].IEN.size());
	BezierElement3D bzel;
	vector<double> Bt;
	vector<array<double, 3>> dBdt;
	bzel.Basis(u[0], u[1], u[2], Bt, dBdt);
	for (uint i = 0; i < hmesh[lev][eid].IEN.size(); i++)
	{
		Nt[i] = 0.;
		dNdt[i][0] = 0.; dNdt[i][1] = 0.; dNdt[i][2] = 0.;
		for (uint j = 0; j < 64; j++)
		{
			Nt[i] += hmesh[lev][eid].bemat[i][j] * Bt[j];
			dNdt[i][0] += hmesh[lev][eid].bemat[i][j] * dBdt[j][0];
			dNdt[i][1] += hmesh[lev][eid].bemat[i][j] * dBdt[j][1];
			dNdt[i][2] += hmesh[lev][eid].bemat[i][j] * dBdt[j][2];
		}
	}
}

void TruncatedTspline_3D::GeomMap_Lev(int lev, int eid, const array<double, 3>& u, array<double, 3>& pt)
{
	vector<double> Nt;
	vector<array<double, 3>> dNdt;
	if (hmesh[lev][eid].type == 0)
	{
		Basis_Regular_Lev(lev, eid, u, Nt, dNdt);
	}
	else if (hmesh[lev][eid].type == 1)
	{
		Basis_Irregular_Lev(lev, eid, u, Nt, dNdt);
	}
	else if (hmesh[lev][eid].type == 2)
	{
		Basis_Irregular_Lev(lev, eid, u, Nt, dNdt);
	}
	pt[0] = 0.; pt[1] = 0.; pt[2] = 0.;
	for (uint i = 0; i < hmesh[lev][eid].IEN.size(); i++)
	{
		int pid(hmesh[lev][eid].IEN[i]);
		pt[0] += Nt[i] * hcp[lev][pid].coor[0];
		pt[1] += Nt[i] * hcp[lev][pid].coor[1];
		pt[2] += Nt[i] * hcp[lev][pid].coor[2];
	}
}

double TruncatedTspline_3D::BasisSum_Lev(int lev, int eid, const array<double, 3>& u)
{
	vector<double> Nt;
	vector<array<double, 3>> dNdt;
	if (hmesh[lev][eid].type == 0)
	{
		Basis_Regular_Lev(lev, eid, u, Nt, dNdt);
	}
	else if (hmesh[lev][eid].type == 1)
	{
		Basis_Irregular_Lev(lev, eid, u, Nt, dNdt);
	}
	else if (hmesh[lev][eid].type == 2)
	{
		Basis_Irregular_Lev(lev, eid, u, Nt, dNdt);
	}
	double sum(0.);
	for (uint i = 0; i < Nt.size(); i++) sum += Nt[i];
	return sum;
}

void TruncatedTspline_3D::AllBezierLev(int lev)
{
	cout << "# elements: " << hmesh[lev].size() << "\n";
#pragma omp parallel for
	for (int i = 0; i < hmesh[lev].size(); i++)
	{
		if (i != 0 && i % 500 == 0)
		{
			cout << i << " ";
		}
		if (hmesh[lev][i].type != 1)
		{
			ConstructBezierBasis(lev, i);
		}
		else
		{
			ConstructBezierBasis_Boundary(lev, i);
		}
	}
	//for (int i = 0; i < hmesh[lev].size(); i++)
	//{
	//	if (hmesh[lev][i].act==1 && hmesh[lev][i].type != 1)
	//	{
	//		ConstructBezierBasis(lev,i);
	//	}
	//}
}

void TruncatedTspline_3D::AnalysisInterface_Elastic(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh)//test, cube
{
	bzmesh.clear();
	IDBC.clear();
	gh.clear();
	double tmp(1.e6);
	double x_range[3][2] = { { tmp, -tmp }, { tmp, -tmp }, { tmp, -tmp } };
	uint i, j, k, k1;
	for (i = 0; i < hcp[0].size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (hcp[0][i].coor[j] < x_range[j][0]) x_range[j][0] = hcp[0][i].coor[j];
			if (hcp[0][i].coor[j] > x_range[j][1]) x_range[j][1] = hcp[0][i].coor[j];
		}
	}
	double xh[3] = { x_range[0][1] - x_range[0][0], x_range[1][1] - x_range[1][0], x_range[2][1] - x_range[2][0]};
	int loc(0);
	vector<vector<int>> aloc(hcp.size());
	for (i = 0; i < hcp.size(); i++)
	{
		aloc[i].resize(hcp[i].size(), -1);
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				aloc[i][j] = loc++;
			}
		}
	}
	IDBC.resize(3*loc);
	gh.resize(3 * loc,0.);
	loc = 0;
	int count(0);
	for (i = 0; i < hcp.size(); i++)
	{
		//aloc[i].resize(hcp[i].size(), -1);
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				if (hcp[i][j].coor[0] == x_range[0][0])
				{
					IDBC[3 * loc] = -1;
				}
				else if (hcp[i][j].coor[0] == x_range[0][1])
				{
					IDBC[3 * loc] = -1; gh[3 * loc] = .1*xh[0];
				}
				else
				{
					IDBC[3 * loc] = count++;
				}
				if (hcp[i][j].coor[1] == x_range[1][0])
				{
					IDBC[3 * loc+1] = -1;
				}
				else
				{
					IDBC[3 * loc+1] = count++;
				}
				if (hcp[i][j].coor[2] == x_range[2][0])
				{
					IDBC[3 * loc + 2] = -1;
				}
				else
				{
					IDBC[3 * loc + 2] = count++;
				}
				loc++;
			}
		}
	}
	for (i = 0; i < hmesh.size(); i++)
	{
		for (j = 0; j < hmesh[i].size(); j++)
		{
			if (hmesh[i][j].act == 1 /*&& hmesh[i][j].type!=1*/)
			{
				BezierElement3D bztmp;
				if (hmesh[i][j].trun == 0)
				{
					bztmp.IEN.resize(hmesh[i][j].IEN.size());
					bztmp.cmat.resize(hmesh[i][j].IEN.size(), vector<double>(64));
					for (k = 0; k < hmesh[i][j].IEN.size(); k++)
					{
						bztmp.IEN[k] = aloc[i][hmesh[i][j].IEN[k]];
						for (k1 = 0; k1 < 64; k1++)
						{
							bztmp.cmat[k][k1] = hmesh[i][j].bemat[k][k1];
							bztmp.pts[k1][0] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[0];
							bztmp.pts[k1][1] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[1];
							bztmp.pts[k1][2] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[2];
						}
					}
				}
				else
				{
				}
				bzmesh.push_back(bztmp);
			}
		}
	}
}

void TruncatedTspline_3D::AnalysisInterface_Poisson(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh)
{
	bzmesh.clear();
	IDBC.clear();
	gh.clear();
	double tmp(1.e6);
	double x_range[3][2] = { { tmp, -tmp }, { tmp, -tmp }, { tmp, -tmp } };
	uint i, j, k, k1, k2;
	for (i = 0; i < hcp[0].size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (hcp[0][i].coor[j] < x_range[j][0]) x_range[j][0] = hcp[0][i].coor[j];
			if (hcp[0][i].coor[j] > x_range[j][1]) x_range[j][1] = hcp[0][i].coor[j];
		}
	}
	double xh[3] = { x_range[0][1] - x_range[0][0], x_range[1][1] - x_range[1][0], x_range[2][1] - x_range[2][0] };
	int loc(0);
	vector<vector<int>> aloc(hcp.size());
	for (i = 0; i < hcp.size(); i++)
	{
		aloc[i].resize(hcp[i].size(), -1);
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				aloc[i][j] = loc++;
			}
		}
	}
	IDBC.resize(loc);
	gh.resize(loc, 0.);
	loc = 0;
	int count(0);
	for (i = 0; i < hcp.size(); i++)
	{
		//aloc[i].resize(hcp[i].size(), -1);
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				if (hcp[i][j].coor[0] == x_range[0][0] || hcp[i][j].coor[0] == x_range[0][1] ||
					hcp[i][j].coor[1] == x_range[1][0] || hcp[i][j].coor[1] == x_range[1][1] || 
					hcp[i][j].coor[2] == x_range[2][0] || hcp[i][j].coor[2] == x_range[2][1])
				{
					IDBC[loc] = -1;
					//gh[loc] = SpecifyDirichBC_2(hcp[i][j].coor);
					//gh[loc] = SpecifyDirichBC_3(hcp[i][j].coor);
					gh[loc] = SpecifyDirichBC_4(hcp[i][j].coor);
				}
				else
				{
					IDBC[loc] = count++;
				}
				loc++;
			}
		}
	}
	for (i = 0; i < hmesh.size(); i++)
	{
		for (j = 0; j < hmesh[i].size(); j++)
		{
			if (hmesh[i][j].act == 1 /*&& hmesh[i][j].type!=1*/)
			{
				BezierElement3D bztmp;
				bztmp.prt[0] = i; bztmp.prt[1] = j;
				bztmp.trun = hmesh[i][j].trun;

				for (k1 = 0; k1 < 64; k1++)
				{
					bztmp.pts[k1][0] = hmesh[i][j].bzpt[k1][0];
					bztmp.pts[k1][1] = hmesh[i][j].bzpt[k1][1];
					bztmp.pts[k1][2] = hmesh[i][j].bzpt[k1][2];
				}

				if (hmesh[i][j].trun == 0)
				{
					bztmp.IEN.resize(hmesh[i][j].IEN.size());
					bztmp.cmat.resize(hmesh[i][j].IEN.size(), vector<double>(64));
					for (k = 0; k < hmesh[i][j].IEN.size(); k++)
					{
						bztmp.IEN[k] = aloc[i][hmesh[i][j].IEN[k]];
						for (k1 = 0; k1 < 64; k1++)
						{
							bztmp.cmat[k][k1] = hmesh[i][j].bemat[k][k1];
							//bztmp.pts[k1][0] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[0];
							//bztmp.pts[k1][1] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[1];
							//bztmp.pts[k1][2] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[2];
						}
					}
				}
				else
				{
					bztmp.IEN.resize(hmesh[i][j].IEN_act.size());
					bztmp.cmat.resize(hmesh[i][j].IEN_act.size(), vector<double>(64));
					for (k = 0; k < hmesh[i][j].IEN_act.size(); k++)
					{
						int lev(hmesh[i][j].IEN_act[k][0]);
						int pid(hmesh[i][j].IEN_act[k][1]);
						if (aloc[lev][pid]==-1)
						{
							cout << "wrong aloc!\n";
							getchar();
						}
						bztmp.IEN[k] = aloc[lev][pid];
						for (k1 = 0; k1 < 64; k1++)
						{
							bztmp.cmat[k][k1] = 0.;
							for (k2 = 0; k2 < hmesh[i][j].IEN.size(); k2++)
							{
								bztmp.cmat[k][k1] += hmesh[i][j].tmat[k][k2] * hmesh[i][j].bemat[k2][k1];
							}
						}
					}
					//for (k = 0; k < hmesh[i][j].IEN.size(); k++)
					//{
					//	for (k1 = 0; k1 < 64; k1++)
					//	{
					//		bztmp.pts[k1][0] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[0];
					//		bztmp.pts[k1][1] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[1];
					//		bztmp.pts[k1][2] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[2];
					//	}
					//}
				}
				bzmesh.push_back(bztmp);
			}
		}
	}
}

void TruncatedTspline_3D::AnalysisInterface_Poisson_1(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh)
{
	bzmesh.clear();
	IDBC.clear();
	gh.clear();
	//int i, j, k, k1, k2;
	int loc(0);
	vector<vector<int>> aloc(hcp.size());
	for (uint i = 0; i < hcp.size(); i++)
	{
		aloc[i].resize(hcp[i].size(), -1);
		for (uint j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				aloc[i][j] = loc++;
			}
		}
	}
	FittingBC(IDBC, gh);
	cout << "Bezier extracting...\n";
	//cout << "# Bezier: " << hmesh[0].size() << "\n";
	vector<array<int, 2>> eact;
	for (uint i = 0; i < hmesh.size(); i++)
	{
		for (uint j = 0; j < hmesh[i].size(); j++)
		{
			if (hmesh[i][j].act == 1)
			{
				array<int, 2> tmp = { i, j };
				eact.push_back(tmp);
			}
		}
	}
	bzmesh.resize(eact.size());
//#pragma omp parallel for
	//for (int i = 0; i < hmesh.size(); i++)
	{
		//for (int j = 0; j < hmesh[i].size(); j++)
#pragma omp parallel for
		for (int eid = 0; eid < eact.size(); eid++)
		{
			int i(eact[eid][0]), j(eact[eid][1]);
			//if (hmesh[i][j].act == 1)
			{
				if (eid != 0 && eid % 500 == 0)
				{
					cout << eid << " ";
				}
				double tmp;
				//BezierElement3D bztmp;
				bzmesh[eid].prt[0] = i; bzmesh[eid].prt[1] = j;
				bzmesh[eid].trun = hmesh[i][j].trun;
				if (hmesh[i][j].type == 1) bzmesh[eid].type = 1;
				if (hmesh[i][j].trun == 0)
				{
					bzmesh[eid].IEN.resize(hmesh[i][j].IEN.size());
					bzmesh[eid].cmat.resize(hmesh[i][j].IEN.size(), vector<double>(64));
					for (uint k = 0; k < hmesh[i][j].IEN.size(); k++)
					{
						bzmesh[eid].IEN[k] = aloc[i][hmesh[i][j].IEN[k]];
						for (int k1 = 0; k1 < 64; k1++)
						{
							tmp = hmesh[i][j].bemat[k][k1];
							bzmesh[eid].cmat[k][k1] = tmp;
							bzmesh[eid].pts[k1][0] += tmp * hcp[i][hmesh[i][j].IEN[k]].coor[0];
							bzmesh[eid].pts[k1][1] += tmp * hcp[i][hmesh[i][j].IEN[k]].coor[1];
							bzmesh[eid].pts[k1][2] += tmp * hcp[i][hmesh[i][j].IEN[k]].coor[2];
						}
					}
				}
				else
				{
					bzmesh[eid].IEN.resize(hmesh[i][j].IEN_act.size());
					bzmesh[eid].cmat.resize(hmesh[i][j].IEN_act.size(), vector<double>(64));
					for (int k = 0; k < hmesh[i][j].IEN_act.size(); k++)
					{
						int lev(hmesh[i][j].IEN_act[k][0]);
						int pid(hmesh[i][j].IEN_act[k][1]);
						if (aloc[lev][pid] == -1)
						{
							cout << "wrong aloc!\n";
							getchar();
						}
						bzmesh[eid].IEN[k] = aloc[lev][pid];
						for (int k1 = 0; k1 < 64; k1++)
						{
							bzmesh[eid].cmat[k][k1] = 0.;
							for (int k2 = 0; k2 < hmesh[i][j].IEN.size(); k2++)
							{
								bzmesh[eid].cmat[k][k1] += hmesh[i][j].tmat[k][k2] * hmesh[i][j].bemat[k2][k1];
							}
						}
					}
					for (int k = 0; k < hmesh[i][j].IEN.size(); k++)
					{
						for (int k1 = 0; k1 < 64; k1++)
						{
							tmp = hmesh[i][j].bemat[k][k1];
							bzmesh[eid].pts[k1][0] += tmp * hcp[i][hmesh[i][j].IEN[k]].coor[0];
							bzmesh[eid].pts[k1][1] += tmp * hcp[i][hmesh[i][j].IEN[k]].coor[1];
							bzmesh[eid].pts[k1][2] += tmp * hcp[i][hmesh[i][j].IEN[k]].coor[2];
						}
					}
				}
				//bzmesh.push_back(bztmp);
			}
		}
	}
}

//void TruncatedTspline_3D::AnalysisInterface_Poisson_1(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh)//parallelized for global refinement
//{
//	bzmesh.clear();
//	IDBC.clear();
//	gh.clear();
//	//int i, j, k, k1, k2;
//	int loc(0);
//	vector<vector<int>> aloc(hcp.size());
//	for (uint i = 0; i < hcp.size(); i++)
//	{
//		aloc[i].resize(hcp[i].size(), -1);
//		for (uint j = 0; j < hcp[i].size(); j++)
//		{
//			if (hcp[i][j].act == 1)
//			{
//				aloc[i][j] = loc++;
//			}
//		}
//	}
//	FittingBC(IDBC, gh);
//	cout << "Bezier extracting...\n";
//	cout << "# Bezier: " << hmesh[0].size() << "\n";
//	bzmesh.resize(hmesh[0].size());
//	int i(0);
//	//#pragma omp parallel for
//	//for (int i = 0; i < hmesh.size(); i++)
//	{
//#pragma omp parallel for
//		for (int j = 0; j < hmesh[i].size(); j++)
//		{
//			//if (hmesh[i][j].act == 1)
//			{
//				if (j != 0 && j % 200 == 0)
//				{
//					cout << j << " ";
//				}
//				double tmp;
//				bzmesh[j].prt[0] = i; bzmesh[j].prt[1] = j;
//				bzmesh[j].trun = hmesh[i][j].trun;
//				if (hmesh[i][j].type == 1) bzmesh[j].type = 1;
//				if (hmesh[i][j].trun == 0)
//				{
//					bzmesh[j].IEN.resize(hmesh[i][j].IEN.size());
//					bzmesh[j].cmat.resize(hmesh[i][j].IEN.size(), vector<double>(64));
//					for (uint k = 0; k < hmesh[i][j].IEN.size(); k++)
//					{
//						bzmesh[j].IEN[k] = aloc[i][hmesh[i][j].IEN[k]];
//						for (int k1 = 0; k1 < 64; k1++)
//						{
//							tmp = hmesh[i][j].bemat[k][k1];
//							bzmesh[j].cmat[k][k1] = tmp;
//							bzmesh[j].pts[k1][0] += tmp * hcp[i][hmesh[i][j].IEN[k]].coor[0];
//							bzmesh[j].pts[k1][1] += tmp * hcp[i][hmesh[i][j].IEN[k]].coor[1];
//							bzmesh[j].pts[k1][2] += tmp * hcp[i][hmesh[i][j].IEN[k]].coor[2];
//						}
//					}
//				}
//			}
//		}
//	}
//}

void TruncatedTspline_3D::AnalysisInterface_Laplace(const vector<array<int, 2>>& pbc, const vector<double>& pdisp, vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh)
{
	bzmesh.clear();
	IDBC.clear();
	gh.clear();
	uint i, j, k, k1, k2;
	int loc(0);
	vector<vector<int>> aloc(hcp.size());
	for (i = 0; i < hcp.size(); i++)
	{
		aloc[i].resize(hcp[i].size(), -1);
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				aloc[i][j] = loc++;
			}
		}
	}
	IDBC.resize(loc);
	gh.resize(loc, 0.);
	loc = 0;
	int count(0);
	for (i = 0; i < hcp.size(); i++)
	{
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				array<int, 2> ptmp = { i, j };
				vector<array<int, 2>>::const_iterator it = find(pbc.begin(), pbc.end(), ptmp);
				if (it != pbc.end())
				{
					IDBC[loc] = -1;
					gh[loc] = pdisp[it - pbc.begin()];
					//cout << gh[loc] << "\n";
					//getchar();
				}
				else
				{
					IDBC[loc] = count++;
				}
				loc++;
			}
		}
	}
	for (i = 0; i < hmesh.size(); i++)
	{
		for (j = 0; j < hmesh[i].size(); j++)
		{
			if (hmesh[i][j].act == 1 /*&& hmesh[i][j].type!=1*/)
			{
				BezierElement3D bztmp;
				bztmp.prt[0] = i; bztmp.prt[1] = j;
				bztmp.trun = hmesh[i][j].trun;
				if (hmesh[i][j].type == 1) bztmp.type = 1;
				if (hmesh[i][j].trun == 0)
				{
					bztmp.IEN.resize(hmesh[i][j].IEN.size());
					bztmp.cmat.resize(hmesh[i][j].IEN.size(), vector<double>(64));
					for (k = 0; k < hmesh[i][j].IEN.size(); k++)
					{
						bztmp.IEN[k] = aloc[i][hmesh[i][j].IEN[k]];
						for (k1 = 0; k1 < 64; k1++)
						{
							bztmp.cmat[k][k1] = hmesh[i][j].bemat[k][k1];
							bztmp.pts[k1][0] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[0];
							bztmp.pts[k1][1] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[1];
							bztmp.pts[k1][2] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[2];
						}
					}
				}
				else
				{
					bztmp.IEN.resize(hmesh[i][j].IEN_act.size());
					bztmp.cmat.resize(hmesh[i][j].IEN_act.size(), vector<double>(64));
					for (k = 0; k < hmesh[i][j].IEN_act.size(); k++)
					{
						int lev(hmesh[i][j].IEN_act[k][0]);
						int pid(hmesh[i][j].IEN_act[k][1]);
						if (aloc[lev][pid] == -1)
						{
							cout << "wrong aloc!\n";
							getchar();
						}
						bztmp.IEN[k] = aloc[lev][pid];
						for (k1 = 0; k1 < 64; k1++)
						{
							bztmp.cmat[k][k1] = 0.;
							for (k2 = 0; k2 < hmesh[i][j].IEN.size(); k2++)
							{
								bztmp.cmat[k][k1] += hmesh[i][j].tmat[k][k2] * hmesh[i][j].bemat[k2][k1];
							}
						}
					}
					for (k = 0; k < hmesh[i][j].IEN.size(); k++)
					{
						for (k1 = 0; k1 < 64; k1++)
						{
							bztmp.pts[k1][0] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[0];
							bztmp.pts[k1][1] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[1];
							bztmp.pts[k1][2] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[2];
						}
					}
				}
				bzmesh.push_back(bztmp);
			}
		}
	}
}

void TruncatedTspline_3D::AnalysisInterface_LeastSquare(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh)
{
	bzmesh.clear();
	IDBC.clear();
	gh.clear();
	double tmp(1.e6);
	double x_range[3][2] = { { tmp, -tmp }, { tmp, -tmp }, { tmp, -tmp } };
	uint i, j, k, k1, k2;
	for (i = 0; i < hcp[0].size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (hcp[0][i].coor[j] < x_range[j][0]) x_range[j][0] = hcp[0][i].coor[j];
			if (hcp[0][i].coor[j] > x_range[j][1]) x_range[j][1] = hcp[0][i].coor[j];
		}
	}
	double xh[3] = { x_range[0][1] - x_range[0][0], x_range[1][1] - x_range[1][0], x_range[2][1] - x_range[2][0] };
	int loc(0);
	vector<vector<int>> aloc(hcp.size());
	for (i = 0; i < hcp.size(); i++)
	{
		aloc[i].resize(hcp[i].size(), -1);
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				aloc[i][j] = loc++;
			}
		}
	}
	IDBC.resize(loc);
	gh.resize(loc, 0.);
	loc = 0;
	int count(0);
	for (i = 0; i < hcp.size(); i++)
	{
		//aloc[i].resize(hcp[i].size(), -1);
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				//if (hcp[i][j].coor[0] == x_range[0][0] || hcp[i][j].coor[0] == x_range[0][1] ||
				//	hcp[i][j].coor[1] == x_range[1][0] || hcp[i][j].coor[1] == x_range[1][1] ||
				//	hcp[i][j].coor[2] == x_range[2][0] || hcp[i][j].coor[2] == x_range[2][1])
				//{
				//	IDBC[loc] = -1;
				//	//gh[loc] = SpecifyDirichBC_2(hcp[i][j].coor);
				//	//gh[loc] = SpecifyDirichBC_3(hcp[i][j].coor);
				//	gh[loc] = SpecifyDirichBC_4(hcp[i][j].coor);
				//}
				//else
				//{
					IDBC[loc] = count++;
				//}
				loc++;
			}
		}
	}
	for (i = 0; i < hmesh.size(); i++)
	{
		for (j = 0; j < hmesh[i].size(); j++)
		{
			if (hmesh[i][j].act == 1 /*&& hmesh[i][j].type!=1*/)
			{
				BezierElement3D bztmp;
				bztmp.prt[0] = i; bztmp.prt[1] = j;
				bztmp.trun = hmesh[i][j].trun;

				for (k1 = 0; k1 < 64; k1++)
				{
					bztmp.pts[k1][0] = hmesh[i][j].bzpt[k1][0];
					bztmp.pts[k1][1] = hmesh[i][j].bzpt[k1][1];
					bztmp.pts[k1][2] = hmesh[i][j].bzpt[k1][2];
				}

				if (hmesh[i][j].trun == 0)
				{
					bztmp.IEN.resize(hmesh[i][j].IEN.size());
					bztmp.cmat.resize(hmesh[i][j].IEN.size(), vector<double>(64));
					for (k = 0; k < hmesh[i][j].IEN.size(); k++)
					{
						bztmp.IEN[k] = aloc[i][hmesh[i][j].IEN[k]];
						for (k1 = 0; k1 < 64; k1++)
						{
							bztmp.cmat[k][k1] = hmesh[i][j].bemat[k][k1];
							//bztmp.pts[k1][0] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[0];
							//bztmp.pts[k1][1] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[1];
							//bztmp.pts[k1][2] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[2];
						}
					}
				}
				else
				{
					bztmp.IEN.resize(hmesh[i][j].IEN_act.size());
					bztmp.cmat.resize(hmesh[i][j].IEN_act.size(), vector<double>(64));
					for (k = 0; k < hmesh[i][j].IEN_act.size(); k++)
					{
						int lev(hmesh[i][j].IEN_act[k][0]);
						int pid(hmesh[i][j].IEN_act[k][1]);
						if (aloc[lev][pid] == -1)
						{
							cout << "wrong aloc!\n";
							getchar();
						}
						bztmp.IEN[k] = aloc[lev][pid];
						for (k1 = 0; k1 < 64; k1++)
						{
							bztmp.cmat[k][k1] = 0.;
							for (k2 = 0; k2 < hmesh[i][j].IEN.size(); k2++)
							{
								bztmp.cmat[k][k1] += hmesh[i][j].tmat[k][k2] * hmesh[i][j].bemat[k2][k1];
							}
						}
					}
					//for (k = 0; k < hmesh[i][j].IEN.size(); k++)
					//{
					//	for (k1 = 0; k1 < 64; k1++)
					//	{
					//		bztmp.pts[k1][0] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[0];
					//		bztmp.pts[k1][1] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[1];
					//		bztmp.pts[k1][2] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[2];
					//	}
					//}
				}
				bzmesh.push_back(bztmp);
			}
		}
	}
}

double TruncatedTspline_3D::SpecifyDirichBC(double x[3])
{
	//double tmp1 = SpecifyDirichBC_5(x);
	//double tmp1 = SpecifyDirichBC_6(x);
	double tmp1 = SpecifyDirichBC_7(x);
	//double tmp1 = SpecifyDirichBC_8(x);
	return tmp1;
}

double TruncatedTspline_3D::SpecifyDirichBC_2(double x[3])
{
	double tmp1 = 1. / exp((20.*x[0] - 7.)*(20.*x[0] - 7.) + (20.*x[1] - 7.)*(20.*x[1] - 7.) + (20.*x[2] - 7.)*(20.*x[2] - 7.));
	double tmp2 = 1. / exp((20.*x[0] - 13.)*(20.*x[0] - 13.) + (20.*x[1] - 13.)*(20.*x[1] - 13.) + (20.*x[2] - 13.)*(20.*x[2] - 13.));
	return 2.*(tmp1 + tmp2) / 3.;
}

double TruncatedTspline_3D::SpecifyDirichBC_3(double x[3])
{
	return 2./(3.*exp((20.*x[0] - 10.)*(20.*x[0] - 10.) + (20.*x[1] - 10.)*(20.*x[1] - 10.) + (20.*x[2] - 10.)*(20.*x[2] - 10.)));
}

double TruncatedTspline_3D::SpecifyDirichBC_4(double x[3])
{
	return tanh(1. - 50.*(-0.408248*x[0] - 0.408248*x[1] + 0.816497*x[2]));
}

double TruncatedTspline_3D::SpecifyDirichBC_5(double x[3])
{
	return (x[0] * x[0] + x[1] * x[1] + x[2] * x[2])*(x[0] + x[1] + x[2]) + x[0] * x[1] * x[2];
}

double TruncatedTspline_3D::SpecifyDirichBC_6(double x[3])
{
	return x[0] * (1. - x[0])*x[1] * (1. - x[1])*x[2] * (1. - x[2]);
}

double TruncatedTspline_3D::SpecifyDirichBC_7(double x[3])
{
	double x1[3] = { (x[0] - dmrg[0][0]), (x[1] - dmrg[1][0]), (x[2] - dmrg[2][0]) };
	return tanh(acoef*(nmpl[0] * x1[0] + nmpl[1] * x1[1] + nmpl[2] * x1[2]));
}

double TruncatedTspline_3D::SpecifyDirichBC_8(double x[3])
{
	double x1[3] = { (x[0] - dmrg[0][0]) / dmlen[0], (x[1] - dmrg[1][0]) / dmlen[1], (x[2] - dmrg[2][0]) / dmlen[2] };
	return 1. / exp(acoef*((x1[0] - nmpl[0])*(x1[0] - nmpl[0]) + (x1[1] - nmpl[1])*(x1[1] - nmpl[1]) + (x1[2] - nmpl[2])*(x1[2] - nmpl[2])));
}

void TruncatedTspline_3D::Identify_Poisson(const vector<array<double, 2>>& eh, const vector<double>& err, vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst)
{
	uint i, j, k;
	vector<vector<double>> err_bf(hcp.size());
	for (i = 0; i < hcp.size(); i++)
	{
		err_bf[i].resize(hcp[i].size(),0.);
	}
	for (i = 0; i < eh.size(); i++)
	{
		if (hmesh[eh[i][0]][eh[i][1]].trun == 0)
		{
			for (j = 0; j < hmesh[eh[i][0]][eh[i][1]].IEN.size(); j++)
			{
				err_bf[eh[i][0]][hmesh[eh[i][0]][eh[i][1]].IEN[j]] += err[i];
			}
		}
		else
		{
			for (j = 0; j < hmesh[eh[i][0]][eh[i][1]].IEN_act.size(); j++)
			{
				int tmp[2] = { hmesh[eh[i][0]][eh[i][1]].IEN_act[j][0], hmesh[eh[i][0]][eh[i][1]].IEN_act[j][1] };
				err_bf[tmp[0]][tmp[1]] += err[i];
			}
		}
	}
	double err_max(0.);
	for (i = 0; i < hcp.size(); i++)
	{
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				if (err_bf[i][j]>err_max) err_max = err_bf[i][j];
			}
		}
	}
	double eta(0.7);
	double tol(eta*err_max);
	vector<array<int, 2>> rf_pid;
	for (i = 0; i < hcp.size(); i++)
	{
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				if (err_bf[i][j]>tol)
				{
					array<int, 2> rftmp = { i, j };
					rf_pid.push_back(rftmp);
				}
			}
		}
	}
	rfid.clear();
	gst.clear();
	for (i = 0; i < rf_pid.size(); i++)
	{
		//check support, refine 1-ring as long as 1-ring elements are at the same level
		for (j = 0; j < hcp[rf_pid[i][0]][rf_pid[i][1]].hex.size(); j++)
		{
			if (hmesh[rf_pid[i][0]][hcp[rf_pid[i][0]][rf_pid[i][1]].hex[j]].act == 1)
			{
				array<int, 2> tmp = { rf_pid[i][0], hcp[rf_pid[i][0]][rf_pid[i][1]].hex[j] };
				vector<array<int, 2>>::iterator it = find(rfid.begin(), rfid.end(), tmp);
				if (it == rfid.end())
					rfid.push_back(tmp);
			}
		}
	}
	//cout << "Max error: "<< err_max << "\n";
	//cout << "# of pid: " << rf_pid.size() << "\n";
	//cout << "# of rfid: " << rfid.size() << "\n";
	//getchar();
	for (i = 0; i < rfid.size(); i++)
	{
		int lev(rfid[i][0]);
		for (j = 0; j < 8; j++)
		{
			int pt(hmesh[lev][rfid[i][1]].cnct[j]);
			for (k = 0; k < hcp[lev][pt].hex.size(); k++)
			{
				if (hmesh[lev][hcp[lev][pt].hex[k]].chd.size() == 0)
				{
					array<int, 2> tmp = { lev, hcp[lev][pt].hex[k] };
					vector<array<int, 2>>::iterator it1 = find(gst.begin(), gst.end(), tmp);
					vector<array<int, 2>>::iterator it2 = find(rfid.begin(), rfid.end(), tmp);
					if (it1 == gst.end() && it2 == rfid.end())
					{
						gst.push_back(tmp);
					}
				}
			}
		}
	}
}

void TruncatedTspline_3D::Identify_Poisson_1(const vector<array<double, 2>>& eh, const vector<double>& err, vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst)
{
	uint i, j, k;
	vector<vector<double>> err_bf(hcp.size());
	for (i = 0; i < hcp.size(); i++)
	{
		err_bf[i].resize(hcp[i].size(), 0.);
	}
	for (i = 0; i < eh.size(); i++)
	{
		if (hmesh[eh[i][0]][eh[i][1]].trun == 0)
		{
			for (j = 0; j < hmesh[eh[i][0]][eh[i][1]].IEN.size(); j++)
			{
				err_bf[eh[i][0]][hmesh[eh[i][0]][eh[i][1]].IEN[j]] += err[i];
			}
		}
		else
		{
			for (j = 0; j < hmesh[eh[i][0]][eh[i][1]].IEN_act.size(); j++)
			{
				int tmp[2] = { hmesh[eh[i][0]][eh[i][1]].IEN_act[j][0], hmesh[eh[i][0]][eh[i][1]].IEN_act[j][1] };
				err_bf[tmp[0]][tmp[1]] += err[i];
			}
		}
	}
	double err_max(0.);
	for (i = 0; i < hcp.size(); i++)
	{
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				if (err_bf[i][j]>err_max) err_max = err_bf[i][j];
			}
		}
	}
	//double eta(0.7);//hook
	double eta(0.3);//base, head
	//double eta(0.);//
	double tol(eta*err_max);
	vector<array<int, 2>> rf_pid;
	for (i = 0; i < hcp.size(); i++)
	{
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				if (err_bf[i][j]>tol)
				{
					array<int, 2> rftmp = { i, j };
					rf_pid.push_back(rftmp);
				}
			}
		}
	}
	rfid.clear();
	gst.clear();
	for (i = 0; i < rf_pid.size(); i++)
	{
		//check support, refine 1-ring as long as 1-ring elements are at the same level
		int ring1(0);
		for (j = 0; j < hcp[rf_pid[i][0]][rf_pid[i][1]].hex.size(); j++)
		{
			if (hmesh[rf_pid[i][0]][hcp[rf_pid[i][0]][rf_pid[i][1]].hex[j]].act == 1)
			{
				ring1 = 1;
				array<int, 2> tmp = { rf_pid[i][0], hcp[rf_pid[i][0]][rf_pid[i][1]].hex[j] };
				vector<array<int, 2>>::iterator it = find(rfid.begin(), rfid.end(), tmp);
				if (it == rfid.end())
					rfid.push_back(tmp);
			}
		}
		if (ring1 == 0)//new
		{
			for (j = 0; j < hcp[rf_pid[i][0]][rf_pid[i][1]].supp.size(); j++)
			{
				if (hmesh[rf_pid[i][0]][hcp[rf_pid[i][0]][rf_pid[i][1]].supp[j]].act == 1)
				{
					array<int, 2> tmp = { rf_pid[i][0], hcp[rf_pid[i][0]][rf_pid[i][1]].supp[j] };
					vector<array<int, 2>>::iterator it = find(rfid.begin(), rfid.end(), tmp);
					if (it == rfid.end())
						rfid.push_back(tmp);
				}
			}
		}
	}
	//cout << "Max error: "<< err_max << "\n";
	//cout << "# of pid: " << rf_pid.size() << "\n";
	cout << "# of rfid: " << rfid.size() << "\n";
	//getchar();
	for (i = 0; i < rfid.size(); i++)
	{
		int lev(rfid[i][0]);
		for (j = 0; j < 8; j++)
		{
			int pt(hmesh[lev][rfid[i][1]].cnct[j]);
			for (k = 0; k < hcp[lev][pt].hex.size(); k++)
			{
				if (hmesh[lev][hcp[lev][pt].hex[k]].chd.size() == 0)
				{
					array<int, 2> tmp = { lev, hcp[lev][pt].hex[k] };
					vector<array<int, 2>>::iterator it1 = find(gst.begin(), gst.end(), tmp);
					vector<array<int, 2>>::iterator it2 = find(rfid.begin(), rfid.end(), tmp);
					if (it1 == gst.end() && it2 == rfid.end())
					{
						gst.push_back(tmp);
					}
				}
			}
		}
	}
}

void TruncatedTspline_3D::Identify_Laplace(const vector<array<double, 2>>& eh, const vector<double>& err, vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst)
{
	uint i, j, k;
	vector<vector<double>> err_bf(hcp.size());
	for (i = 0; i < hcp.size(); i++)
	{
		err_bf[i].resize(hcp[i].size(), 0.);
	}
	for (i = 0; i < eh.size(); i++)
	{
		if (hmesh[eh[i][0]][eh[i][1]].trun == 0)
		{
			for (j = 0; j < hmesh[eh[i][0]][eh[i][1]].IEN.size(); j++)
			{
				err_bf[eh[i][0]][hmesh[eh[i][0]][eh[i][1]].IEN[j]] += err[i];
			}
		}
		else
		{
			for (j = 0; j < hmesh[eh[i][0]][eh[i][1]].IEN_act.size(); j++)
			{
				int tmp[2] = { hmesh[eh[i][0]][eh[i][1]].IEN_act[j][0], hmesh[eh[i][0]][eh[i][1]].IEN_act[j][1] };
				err_bf[tmp[0]][tmp[1]] += err[i];
			}
		}
	}
	double err_max(0.);
	for (i = 0; i < hcp.size(); i++)
	{
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				if (err_bf[i][j]>err_max) err_max = err_bf[i][j];
			}
		}
	}
	double eta(0.5);
	double tol(eta*err_max);
	vector<array<int, 2>> rf_pid;
	for (i = 0; i < hcp.size(); i++)
	{
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				if (err_bf[i][j]>tol)
				{
					array<int, 2> rftmp = { i, j };
					rf_pid.push_back(rftmp);
				}
			}
		}
	}
	rfid.clear();
	gst.clear();
	for (i = 0; i < rf_pid.size(); i++)
	{
		//check support, refine 1-ring as long as 1-ring elements are at the same level
		int ring1(0);
		for (j = 0; j < hcp[rf_pid[i][0]][rf_pid[i][1]].hex.size(); j++)
		{
			if (hmesh[rf_pid[i][0]][hcp[rf_pid[i][0]][rf_pid[i][1]].hex[j]].act == 1)
			{
				ring1 = 1;
				array<int, 2> tmp = { rf_pid[i][0], hcp[rf_pid[i][0]][rf_pid[i][1]].hex[j] };
				vector<array<int, 2>>::iterator it = find(rfid.begin(), rfid.end(), tmp);
				if (it == rfid.end())
					rfid.push_back(tmp);
			}
		}
		if (ring1 == 0)//new
		{
			for (j = 0; j < hcp[rf_pid[i][0]][rf_pid[i][1]].supp.size(); j++)
			{
				if (hmesh[rf_pid[i][0]][hcp[rf_pid[i][0]][rf_pid[i][1]].supp[j]].act == 1)
				{
					array<int, 2> tmp = { rf_pid[i][0], hcp[rf_pid[i][0]][rf_pid[i][1]].supp[j] };
					vector<array<int, 2>>::iterator it = find(rfid.begin(), rfid.end(), tmp);
					if (it == rfid.end())
						rfid.push_back(tmp);
				}
			}
		}
	}
	//cout << "Max error: "<< err_max << "\n";
	//cout << "# of pid: " << rf_pid.size() << "\n";
	cout << "# of rfid: " << rfid.size() << "\n";
	//getchar();
	for (i = 0; i < rfid.size(); i++)
	{
		int lev(rfid[i][0]);
		for (j = 0; j < 8; j++)
		{
			int pt(hmesh[lev][rfid[i][1]].cnct[j]);
			for (k = 0; k < hcp[lev][pt].hex.size(); k++)
			{
				if (hmesh[lev][hcp[lev][pt].hex[k]].chd.size() == 0)
				{
					array<int, 2> tmp = { lev, hcp[lev][pt].hex[k] };
					vector<array<int, 2>>::iterator it1 = find(gst.begin(), gst.end(), tmp);
					vector<array<int, 2>>::iterator it2 = find(rfid.begin(), rfid.end(), tmp);
					if (it1 == gst.end() && it2 == rfid.end())
					{
						gst.push_back(tmp);
					}
				}
			}
		}
	}
}

void TruncatedTspline_3D::Identify_LeastSquare(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst)
{
	//identify a sphere
	uint i, j, k;
	double cent[3] = { .5, .5, .5 };
	//double r(0.15);
	double r(0.23);

	vector<array<int, 2>> rf_pid;
	for (i = 0; i < hcp.size(); i++)
	{
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				int flag(0);
				for (k = 0; k<hcp[i][j].hex.size(); k++)
				{
					int hxid(hcp[i][j].hex[k]);
					if (hmesh[i][hxid].act == 1)
					{
						int count(0);
						for (int k0 = 0; k0 < 8; k0++)
						{
							int pid = hmesh[i][hxid].cnct[k0];
							double dist = sqrt((hcp[i][pid].coor[0] - cent[0])*(hcp[i][pid].coor[0] - cent[0])
								+ (hcp[i][pid].coor[1] - cent[1])*(hcp[i][pid].coor[1] - cent[1]) + (hcp[i][pid].coor[2] - cent[2])*(hcp[i][pid].coor[2] - cent[2]));
							//cout << dist << "\n"; getchar();
							if (dist < r) count++;
						}
						if (count>0 && count < 8) flag = 1;
					}
				}
				if (flag==1)
				{
					array<int, 2> rftmp = { i, j };
					rf_pid.push_back(rftmp);
				}
			}
		}
	}

	rfid.clear();
	gst.clear();
	for (i = 0; i < rf_pid.size(); i++)
	{
		//check support, refine 1-ring as long as 1-ring elements are at the same level
		int ring1(0);
		for (j = 0; j < hcp[rf_pid[i][0]][rf_pid[i][1]].hex.size(); j++)
		{
			if (hmesh[rf_pid[i][0]][hcp[rf_pid[i][0]][rf_pid[i][1]].hex[j]].act == 1)
			{
				ring1 = 1;
				array<int, 2> tmp = { rf_pid[i][0], hcp[rf_pid[i][0]][rf_pid[i][1]].hex[j] };
				vector<array<int, 2>>::iterator it = find(rfid.begin(), rfid.end(), tmp);
				if (it == rfid.end())
					rfid.push_back(tmp);
			}
		}
		if (ring1 == 0)//new
		{
			for (j = 0; j < hcp[rf_pid[i][0]][rf_pid[i][1]].supp.size(); j++)
			{
				if (hmesh[rf_pid[i][0]][hcp[rf_pid[i][0]][rf_pid[i][1]].supp[j]].act == 1)
				{
					array<int, 2> tmp = { rf_pid[i][0], hcp[rf_pid[i][0]][rf_pid[i][1]].supp[j] };
					vector<array<int, 2>>::iterator it = find(rfid.begin(), rfid.end(), tmp);
					if (it == rfid.end())
						rfid.push_back(tmp);
				}
			}
		}
	}
	//cout << "Max error: "<< err_max << "\n";
	//cout << "# of pid: " << rf_pid.size() << "\n";
	cout << "# of rfid: " << rfid.size() << "\n";
	//getchar();
	for (i = 0; i < rfid.size(); i++)
	{
		int lev(rfid[i][0]);
		for (j = 0; j < 8; j++)
		{
			int pt(hmesh[lev][rfid[i][1]].cnct[j]);
			for (k = 0; k < hcp[lev][pt].hex.size(); k++)
			{
				if (hmesh[lev][hcp[lev][pt].hex[k]].chd.size() == 0)
				{
					array<int, 2> tmp = { lev, hcp[lev][pt].hex[k] };
					vector<array<int, 2>>::iterator it1 = find(gst.begin(), gst.end(), tmp);
					vector<array<int, 2>>::iterator it2 = find(rfid.begin(), rfid.end(), tmp);
					if (it1 == gst.end() && it2 == rfid.end())
					{
						gst.push_back(tmp);
					}
				}
			}
		}
	}
}

void TruncatedTspline_3D::Identify_LeastSquare_Line(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst)
{
	//identify a sphere
	uint i, j, k;

	vector<array<int, 2>> rf_pid;
	for (i = 0; i < hcp.size(); i++)
	{
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				int flag(0);
				for (k = 0; k<hcp[i][j].hex.size(); k++)
				{
					int hxid(hcp[i][j].hex[k]);
					if (hmesh[i][hxid].act == 1)
					{
						int count(0);
						double hxdm[3][2] = { { 1.e6, -1.e6 }, { 1.e6, -1.e6 }, { 1.e6, -1.e6 } };
						for (int k0 = 0; k0 < 8; k0++)
						{
							int pid = hmesh[i][hxid].cnct[k0];
							for (int k1 = 0; k1<3; k1++)
							{
								if (hcp[i][pid].coor[k1]<hxdm[k1][0]) hxdm[k1][0] = hcp[i][pid].coor[k1];
								if (hcp[i][pid].coor[k1]>hxdm[k1][1]) hxdm[k1][1] = hcp[i][pid].coor[k1];
							}
						}
						double a[4];
						a[0] = hxdm[0][0]>hxdm[1][0] ? hxdm[0][0] : hxdm[1][0];
						a[1] = hxdm[0][1]<hxdm[1][1] ? hxdm[0][1] : hxdm[1][1];
						if (a[0] < a[1])
						{
							a[2] = a[0]>hxdm[2][0] ? a[0] : hxdm[2][0];
							a[3] = a[1]<hxdm[2][1] ? a[1] : hxdm[2][1];
							if (a[2] < a[3]) flag = 1;
						}
					}
				}
				if (flag == 1)
				{
					array<int, 2> rftmp = { i, j };
					rf_pid.push_back(rftmp);
				}
			}
		}
	}

	rfid.clear();
	gst.clear();
	for (i = 0; i < rf_pid.size(); i++)
	{
		//check support, refine 1-ring as long as 1-ring elements are at the same level
		int ring1(0);
		for (j = 0; j < hcp[rf_pid[i][0]][rf_pid[i][1]].hex.size(); j++)
		{
			if (hmesh[rf_pid[i][0]][hcp[rf_pid[i][0]][rf_pid[i][1]].hex[j]].act == 1)
			{
				ring1 = 1;
				array<int, 2> tmp = { rf_pid[i][0], hcp[rf_pid[i][0]][rf_pid[i][1]].hex[j] };
				vector<array<int, 2>>::iterator it = find(rfid.begin(), rfid.end(), tmp);
				if (it == rfid.end())
					rfid.push_back(tmp);
			}
		}
		if (ring1 == 0)//new
		{
			for (j = 0; j < hcp[rf_pid[i][0]][rf_pid[i][1]].supp.size(); j++)
			{
				if (hmesh[rf_pid[i][0]][hcp[rf_pid[i][0]][rf_pid[i][1]].supp[j]].act == 1)
				{
					array<int, 2> tmp = { rf_pid[i][0], hcp[rf_pid[i][0]][rf_pid[i][1]].supp[j] };
					vector<array<int, 2>>::iterator it = find(rfid.begin(), rfid.end(), tmp);
					if (it == rfid.end())
						rfid.push_back(tmp);
				}
			}
		}
	}
	//cout << "Max error: "<< err_max << "\n";
	//cout << "# of pid: " << rf_pid.size() << "\n";
	cout << "# of rfid: " << rfid.size() << "\n";
	//getchar();
	for (i = 0; i < rfid.size(); i++)
	{
		int lev(rfid[i][0]);
		for (j = 0; j < 8; j++)
		{
			int pt(hmesh[lev][rfid[i][1]].cnct[j]);
			for (k = 0; k < hcp[lev][pt].hex.size(); k++)
			{
				if (hmesh[lev][hcp[lev][pt].hex[k]].chd.size() == 0)
				{
					array<int, 2> tmp = { lev, hcp[lev][pt].hex[k] };
					vector<array<int, 2>>::iterator it1 = find(gst.begin(), gst.end(), tmp);
					vector<array<int, 2>>::iterator it2 = find(rfid.begin(), rfid.end(), tmp);
					if (it1 == gst.end() && it2 == rfid.end())
					{
						gst.push_back(tmp);
					}
				}
			}
		}
	}
}

void TruncatedTspline_3D::OutputBasis(int lev, int pid, string fn)
{
	vector<array<double, 3>> spt;
	vector<array<double, 3>> sval;
	vector<double> ssum;
	vector<array<int, 8>> sele;
	vector<array<double, 3>> lpt;//visulize parameter lines
	vector<array<int, 2>> led;//line connectivity
	int ns(5), ne_ref(0), loc0, loc1, loc2;

	for (uint eid = 0; eid < hmesh[lev].size(); eid++)
	{
		if (hmesh[lev][eid].act == 1 /*&& (hmesh[lev][eid].type == 0 || hmesh[lev][eid].type == 2)*/)
		{
			vector<double> su(ns), sv(ns), sw(ns);
			double ul[3] = { hmesh[lev][eid].dm[0][1] - hmesh[lev][eid].dm[0][0], hmesh[lev][eid].dm[1][1] - hmesh[lev][eid].dm[1][0],
				hmesh[lev][eid].dm[2][1] - hmesh[lev][eid].dm[2][0] };
			for (int i = 0; i<ns; i++)
			{
				su[i] = i*ul[0] / (ns - 1) + hmesh[lev][eid].dm[0][0];
				sv[i] = i*ul[1] / (ns - 1) + hmesh[lev][eid].dm[1][0];
				sw[i] = i*ul[2] / (ns - 1) + hmesh[lev][eid].dm[2][0];
			}
			int pos(-1);
			vector<int>::iterator it = find(hmesh[lev][eid].IEN.begin(), hmesh[lev][eid].IEN.end(), pid);
			if (it != hmesh[lev][eid].IEN.end()) pos = it - hmesh[lev][eid].IEN.begin();
			for (int a = 0; a<ns; a++)
			{
				for (int b = 0; b<ns; b++)
				{
					for (int c = 0; c < ns; c++)
					{
						array<double, 3> pt;
						array<double, 3> uval = { su[c], sv[b], sw[a] };
						GeomMap(lev, eid, uval, pt);
						spt.push_back(pt);
						if (pos != -1)
						{
							vector<double> Nt;
							vector<array<double, 3>> dNdt;
							Basis_Irregular(lev, eid, uval, Nt, dNdt);
							ssum.push_back(Nt[pos]);
						}
						else
						{
							ssum.push_back(0.);
						}
					}
				}
			}

			for (int a = 0; a<ns - 1; a++)
			{
				for (int b = 0; b<ns - 1; b++)
				{
					for (int c = 0; c < ns - 1; c++)
					{
						array<int, 8> el;
						el[0] = ne_ref + a*ns*ns + b*ns + c;
						el[1] = ne_ref + a*ns*ns + b*ns + c + 1;
						el[2] = ne_ref + a*ns*ns + (b + 1)*ns + c + 1;
						el[3] = ne_ref + a*ns*ns + (b + 1)*ns + c;
						el[4] = ne_ref + (a + 1)*ns*ns + b*ns + c;
						el[5] = ne_ref + (a + 1)*ns*ns + b*ns + c + 1;
						el[6] = ne_ref + (a + 1)*ns*ns + (b + 1)*ns + c + 1;
						el[7] = ne_ref + (a + 1)*ns*ns + (b + 1)*ns + c;
						sele.push_back(el);
					}
				}
			}
			ne_ref += ns*ns*ns;
		}
	}

	string fname = fn + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (uint i = 0; i<sele.size(); i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i<sele.size(); i++)
		{
			fout << "12\n";
		}
		fout << "\nPOINT_DATA " << ssum.size() << "\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<ssum.size(); i++)
		{
			fout << ssum[i] << "\n";
		}
		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nNORMALS Normal FLOAT\n";
		//for(uint i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i][0]<<" "<<sval[i][1]<<" "<<sval[i][2]<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	//int lev(0);
	//for (uint i = 0; i < hmesh[lev].size(); i++)
	//{
	//	if (hmesh[lev][i].act == 1)
	//	{
	//		cout << "eid: " << i << "\n";
	//		for (uint j = 0; j < hmesh[lev][i].IEN.size(); j++)
	//		{
	//			cout << "pid: " << hmesh[lev][i].IEN[j] << "\n";
	//			for (uint k = 0; k < 64; k++)
	//			{
	//				cout << hmesh[lev][i].bemat[j][k] << " ";
	//			}
	//			cout << "\n";
	//			getchar();
	//		}
	//	}
	//}
}

void TruncatedTspline_3D::PillowCube(string fn)
{
	double tmp(1.e6);
	double x_range[3][2] = { { tmp, -tmp }, { tmp, -tmp }, { tmp, -tmp } };
	uint i, j, k;
	for (i = 0; i < cp.size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (cp[i].coor[j] < x_range[j][0]) x_range[j][0] = cp[i].coor[j];
			if (cp[i].coor[j] > x_range[j][1]) x_range[j][1] = cp[i].coor[j];
		}
	}
	double xh[3] = { x_range[0][1] - x_range[0][0], x_range[1][1] - x_range[1][0], x_range[2][1] - x_range[2][0] };
	int nplw(2);
	int nelx(2*nplw+1);
	int eid_del(nplw*nelx*nelx+nplw*nelx+nplw);
	double plw_h[3] = {xh[0]/nplw,xh[1]/nplw,xh[2]/nplw};
	vector<array<double, 3>> psm(nelx+1);
	for (i = 0; i <= nplw; i++)
	{
		for (j = 0; j < 3; j++)
		{
			psm[i][j] = x_range[j][0] - double(nplw - i)*plw_h[j];
			psm[i+nplw+1][j] = x_range[j][1] + double(i)*plw_h[j];
		}
	}

	int nsmp(nelx);
	int npt = (nsmp + 1)*(nsmp + 1)*(nsmp + 1);
	int nel = nsmp*nsmp*nsmp;
	vector<array<double, 3>> spt(npt);
	vector<array<int, 8>> sele0(nel);
	int loc(0);
	for (i = 0; i < nsmp + 1; i++)
	{
		for (j = 0; j < nsmp + 1; j++)
		{
			for (k = 0; k < nsmp + 1; k++)
			{
				spt[loc][0] = psm[k][0];
				spt[loc][1] = psm[j][1];
				spt[loc][2] = psm[i][2];
				loc++;
			}
		}
	}
	loc = 0;
	for (i = 0; i < nsmp; i++)
	{
		for (j = 0; j < nsmp; j++)
		{
			for (k = 0; k < nsmp; k++)
			{
				sele0[loc][0] = i*(nsmp + 1)*(nsmp + 1) + j*(nsmp + 1) + k;
				sele0[loc][1] = i*(nsmp + 1)*(nsmp + 1) + j*(nsmp + 1) + k + 1;
				sele0[loc][2] = i*(nsmp + 1)*(nsmp + 1) + (j + 1)*(nsmp + 1) + k + 1;
				sele0[loc][3] = i*(nsmp + 1)*(nsmp + 1) + (j + 1)*(nsmp + 1) + k;
				sele0[loc][4] = (i + 1)*(nsmp + 1)*(nsmp + 1) + j*(nsmp + 1) + k;
				sele0[loc][5] = (i + 1)*(nsmp + 1)*(nsmp + 1) + j*(nsmp + 1) + k + 1;
				sele0[loc][6] = (i + 1)*(nsmp + 1)*(nsmp + 1) + (j + 1)*(nsmp + 1) + k + 1;
				sele0[loc][7] = (i + 1)*(nsmp + 1)*(nsmp + 1) + (j + 1)*(nsmp + 1) + k;
				loc++;
			}
		}
	}
	vector<array<int, 8>> sele(nel-1+tmesh.size());
	loc = 0;
	for (i = 0; i < sele0.size(); i++)
	{
		if (i != eid_del)
		{
			for (j = 0; j < 8; j++) sele[loc][j] = sele0[i][j];
			loc++;
		}
	}
	vector<int> id_new(cp.size());
	for (i = 0; i < cp.size(); i++)
	{
		array<double, 3> tmp = { cp[i].coor[0], cp[i].coor[1], cp[i].coor[2] };
		vector<array<double, 3>>::iterator it = find(spt.begin(),spt.end(),tmp);
		id_new[i] = it - spt.begin();
		if (it == spt.end()) spt.push_back(tmp);
	}
	for (i = 0; i < tmesh.size(); i++)
	{
		for (j = 0; j < 8; j++) sele[loc][j] = id_new[tmesh[i].cnct[j]];
		loc++;
	}

	string fname = fn + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nCube hex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (uint i = 0; i<sele.size(); i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i<sele.size(); i++)
		{
			fout << "12\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline_3D::VisualizeBezier(const vector<BezierElement3D>& bzmesh, string fn)
{
	vector<array<int, 8>> cnct(bzmesh.size()*27);
	int loc(0);
	for (uint eid = 0; eid < bzmesh.size(); eid++)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					cnct[loc][0] = 64 * eid + 16 * i + 4 * j + k;
					cnct[loc][1] = 64 * eid + 16 * i + 4 * j + k+1;
					cnct[loc][2] = 64 * eid + 16 * i + 4 * (j+1) + k + 1;
					cnct[loc][3] = 64 * eid + 16 * i + 4 * (j + 1) + k;
					cnct[loc][4] = 64 * eid + 16 * (i+1) + 4 * j + k;
					cnct[loc][5] = 64 * eid + 16 * (i+1) + 4 * j + k + 1;
					cnct[loc][6] = 64 * eid + 16 * (i+1) + 4 * (j + 1) + k + 1;
					cnct[loc][7] = 64 * eid + 16 * (i+1) + 4 * (j + 1) + k;
					loc++;
				}
			}
		}
	}

	string fname = fn + "_bezier.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << 64*bzmesh.size() << " float\n";
		for (uint i = 0; i<bzmesh.size(); i++)
		{
			for (int j = 0; j < 64; j++)
				fout << bzmesh[i].pts[j][0] << " " << bzmesh[i].pts[j][1] << " " << bzmesh[i].pts[j][2] << "\n";
		}
		fout << "\nCELLS " << cnct.size() << " " << 9 * cnct.size() << '\n';
		for (uint i = 0; i<cnct.size(); i++)
		{
			fout << "8 " << cnct[i][0] << " " << cnct[i][1] << " " << cnct[i][2] << " " << cnct[i][3] << " " << cnct[i][4] << " " << cnct[i][5] << " " << cnct[i][6] << " " << cnct[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << cnct.size() << '\n';
		for (uint i = 0; i<cnct.size(); i++)
		{
			fout << "12\n";
		}
		fout << "\nPOINT_DATA " << 64 * bzmesh.size() << "\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<bzmesh.size(); i++)
		{
			for (int j = 0; j < 64; j++)
				fout << "1\n";
		}
		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nNORMALS Normal FLOAT\n";
		//for(uint i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i][0]<<" "<<sval[i][1]<<" "<<sval[i][2]<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

//void TruncatedTspline_3D::SetInitialBC(vector<array<int, 2>>& e_BC, vector<double>& e_disp)
//{
//	//list of element
//	//rod CAD model
//	//int nbc1(10), nbc2(8);
//	//int bc1[] = {1249,1250,1251,1252,1253,1257,1258,1259,1260,1261};
//	//int bc2[] = { 775, 781, 783, 795, 801, 803, 853, 865 };
//	//cube4
//	//int nbc1(1), nbc2(1);
//	//int bc1[] = { 20,24,36,40 };
//	//int bc2[] = { 23,27,39,43 };
//	//cube9
//	int nbc1(4), nbc2(4);
//	int bc1[] = { 679,680,688,689 };
//	int bc2[] = { 40,41,49,50 };
//
//	e_BC.clear();
//	e_disp.clear();
//	for (int i = 0; i < nbc1; i++)
//	{
//		array<int, 2> tmp = { 0, bc1[i] };
//		e_BC.push_back(tmp);
//		e_disp.push_back(100.);
//	}
//
//	for (int i = 0; i < nbc2; i++)
//	{
//		array<int, 2> tmp = { 0, bc2[i] };
//		e_BC.push_back(tmp);
//		e_disp.push_back(0.);
//	}
//}
//
//void TruncatedTspline_3D::SetBC(vector<array<int, 2>>& ebc, vector<double>& edisp, vector<array<int, 2>>& pbc, vector<double>& pdisp)
//{
//	vector<array<int, 2>> ebc_new;
//	vector<double> edisp_new;
//	pbc.clear();
//	pdisp.clear();
//	for (uint i = 0; i < ebc.size(); i++)
//	{
//		ebc_new.push_back(ebc[i]);
//		edisp_new.push_back(edisp[i]);
//		//if (hmesh[ebc[i][0]][ebc[i][1]].act == 1 && hmesh[ebc[i][0]][ebc[i][1]].type == 1)
//		//{
//		//	ebc_new.push_back(ebc[i]);
//		//	edisp_new.push_back(edisp[i]);
//			for (int j = 0; j < 8; j++)
//			{
//				array<int, 2> pid1 = { ebc[i][0], hmesh[ebc[i][0]][ebc[i][1]].cnct[j] };
//				if (hcp[pid1[0]][pid1[1]].act==1 && hcp[pid1[0]][pid1[1]].type == 1)
//				{
//					for (uint k = 0; k < hcp[pid1[0]][pid1[1]].face.size(); k++)
//					{
//						array<int, 2> fcid = { ebc[i][0], hcp[pid1[0]][pid1[1]].face[k] };
//						if (hface[fcid[0]][fcid[1]].type == 1)
//						{
//							for (int k1 = 0; k1 < 4; k1++)
//							{
//								array<int, 2> pid2 = { ebc[i][0], hface[fcid[0]][fcid[1]].cnct[k1] };
//								vector<array<int, 2>>::iterator it = find(pbc.begin(), pbc.end(), pid2);
//								if (it == pbc.end())
//								{
//									pbc.push_back(pid2);
//									pdisp.push_back(edisp[i]);
//								}
//							}
//						}
//					}
//				}
//			}
//
//			/*if (hmesh[ebc[i][0]][ebc[i][1]].trun == 0)
//			{
//				for (uint j = 0; j < hmesh[ebc[i][0]][ebc[i][1]].IEN.size(); j++)
//				{
//					array<int, 2> tmp = { ebc[i][0], hmesh[ebc[i][0]][ebc[i][1]].IEN[j] };
//					if (hcp[tmp[0]][tmp[1]].type == 1)
//					{
//						vector<array<int, 2>>::iterator it = find(pbc.begin(), pbc.end(), tmp);
//						if (it == pbc.end())
//						{
//							pbc.push_back(tmp);
//							pdisp.push_back(edisp[i]);
//						}
//					}
//				}
//			}
//			else
//			{
//				for (uint j = 0; j < hmesh[ebc[i][0]][ebc[i][1]].IEN_act.size(); j++)
//				{
//					array<int, 2> tmp = { hmesh[ebc[i][0]][ebc[i][1]].IEN_act[j][0], hmesh[ebc[i][0]][ebc[i][1]].IEN_act[j][1] };
//					if (hcp[tmp[0]][tmp[1]].type == 1)
//					{
//						vector<array<int, 2>>::iterator it = find(pbc.begin(), pbc.end(), tmp);
//						if (it == pbc.end())
//						{
//							pbc.push_back(tmp);
//							pdisp.push_back(edisp[i]);
//						}
//					}
//				}
//			}*/
//		//}
//		//else if (hmesh[ebc[i][0]][ebc[i][1]].act == 0)
//		{
//			for (uint k2 = 0; k2 < hmesh[ebc[i][0]][ebc[i][1]].chd.size(); k2++)
//			{
//				array<int,2> eid = { ebc[i][0]+1, hmesh[ebc[i][0]][ebc[i][1]].chd[k2] };
//				if (hmesh[eid[0]][eid[1]].act==1 && hmesh[eid[0]][eid[1]].type == 1)
//				{
//					ebc_new.push_back(eid);
//					edisp_new.push_back(edisp[i]);
//					for (int j = 0; j < 8; j++)
//					{
//						array<int, 2> pid1 = { eid[0], hmesh[eid[0]][eid[1]].cnct[j] };
//						if (hcp[pid1[0]][pid1[1]].act == 1 && hcp[pid1[0]][pid1[1]].type == 1)
//						{
//							for (uint k = 0; k < hcp[pid1[0]][pid1[1]].face.size(); k++)
//							{
//								array<int, 2> fcid = { eid[0], hcp[pid1[0]][pid1[1]].face[k] };
//								if (hface[fcid[0]][fcid[1]].type == 1)
//								{
//									for (int k1 = 0; k1 < 4; k1++)
//									{
//										array<int, 2> pid2 = { eid[0], hface[fcid[0]][fcid[1]].cnct[k1] };
//										vector<array<int, 2>>::iterator it = find(pbc.begin(), pbc.end(), pid2);
//										if (it == pbc.end())
//										{
//											pbc.push_back(pid2);
//											pdisp.push_back(edisp[i]);
//										}
//									}
//								}
//							}
//						}
//					}
//
//					/*if (hmesh[eid[0]][eid[1]].trun == 0)
//					{
//						for (uint j = 0; j < hmesh[eid[0]][eid[1]].IEN.size(); j++)
//						{
//							array<int, 2> tmp = { ebc[i][0], hmesh[eid[0]][eid[1]].IEN[j] };
//							if (hcp[tmp[0]][tmp[1]].type == 1)
//							{
//								vector<array<int, 2>>::iterator it = find(pbc.begin(), pbc.end(), tmp);
//								if (it == pbc.end())
//								{
//									pbc.push_back(tmp);
//									pdisp.push_back(edisp[i]);
//								}
//							}
//						}
//					}
//					else
//					{
//						for (uint j = 0; j < hmesh[eid[0]][eid[1]].IEN_act.size(); j++)
//						{
//							array<int, 2> tmp = { hmesh[eid[0]][eid[1]].IEN_act[j][0], hmesh[eid[0]][eid[1]].IEN_act[j][1] };
//							if (hcp[tmp[0]][tmp[1]].type == 1)
//							{
//								vector<array<int, 2>>::iterator it = find(pbc.begin(), pbc.end(), tmp);
//								if (it == pbc.end())
//								{
//									pbc.push_back(tmp);
//									pdisp.push_back(edisp[i]);
//								}
//							}
//						}
//					}*/
//				}
//			}
//		}
//	}
//	ebc.clear();
//	edisp.clear();
//	ebc = ebc_new;
//	edisp = edisp_new;
//
//	//cout << "pbc:\n";
//	//for (uint i = 0; i < pbc.size(); i++)
//	//{
//	//	cout << pbc[i][1] << " " << pdisp[i] << "\n";
//	//}
//	//getchar();
//}

//void TruncatedTspline_3D::SetInitialBC(vector<array<int, 2>>& pbc, vector<double>& pdisp)
//{
//	//list of element
//	//rod CAD model
//	//int nbc1(10), nbc2(8);
//	//int bc1[] = {1249,1250,1251,1252,1253,1257,1258,1259,1260,1261};
//	//int bc2[] = { 775, 781, 783, 795, 801, 803, 853, 865 };
//	//cube4
//	//int nbc1(1), nbc2(1);
//	//int bc1[] = { 20,24,36,40 };
//	//int bc2[] = { 23,27,39,43 };
//	//cube9
//	int nbc1(4), nbc2(4);
//	int bc1[] = { 679, 680, 688, 689 };
//	int bc2[] = { 40, 41, 49, 50 };
//
//	pbc.clear();
//	pdisp.clear();
//	for (int i = 0; i < nbc1; i++)
//	{
//		for (int j = 0; j < 8; j++)
//		{
//			if (hcp[0][hmesh[0][bc1[i]].cnct[j]].type == 1)
//			{
//				int pid(hmesh[0][bc1[i]].cnct[j]);
//				for (uint k = 0; k < hcp[0][pid].face.size(); k++)
//				{
//					if (hface[0][hcp[0][pid].face[k]].type == 1)
//					{
//						for (int k1 = 0; k1 < 4; k1++)
//						{
//							array<int, 2> tmp = { 0, hface[0][hcp[0][pid].face[k]].cnct[k1] };
//							vector<array<int, 2>>::iterator it = find(pbc.begin(), pbc.end(), tmp);
//							if (it == pbc.end())
//							{
//								pbc.push_back(tmp);
//								pdisp.push_back(100.);
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//	for (int i = 0; i < nbc2; i++)
//	{
//		for (int j = 0; j < 8; j++)
//		{
//			if (hcp[0][hmesh[0][bc2[i]].cnct[j]].type == 1)
//			{
//				int pid(hmesh[0][bc2[i]].cnct[j]);
//				for (uint k = 0; k < hcp[0][pid].face.size(); k++)
//				{
//					if (hface[0][hcp[0][pid].face[k]].type == 1)
//					{
//						for (int k1 = 0; k1 < 4; k1++)
//						{
//							array<int, 2> tmp = { 0, hface[0][hcp[0][pid].face[k]].cnct[k1] };
//							vector<array<int, 2>>::iterator it = find(pbc.begin(), pbc.end(), tmp);
//							if (it == pbc.end())
//							{
//								pbc.push_back(tmp);
//								pdisp.push_back(0.);
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//}
//
//void TruncatedTspline_3D::SetBC(vector<array<int, 2>>& pbc, vector<double>& pdisp)
//{
//	vector<array<int, 2>> pbc_new;
//	vector<double> pdisp_new;
//	for (uint i = 0; i < pbc.size(); i++)
//	{
//		if (hcp[pbc[i][0]][pbc[i][1]].act == 1)
//		{
//			pbc_new.push_back(pbc[i]);
//			pdisp_new.push_back(pdisp[i]);
//			for (uint j = 0; j < hcp[pbc[i][0]][pbc[i][1]].chd.size(); j++)
//			{
//				array<int, 2> pid = { pbc[i][0] + 1, hcp[pbc[i][0]][pbc[i][1]].chd[j] };
//				if (hcp[pid[0]][pid[1]].act == 1 && hcp[pid[0]][pid[1]].type==1)
//				{
//					vector<array<int, 2>>::iterator it = find(pbc_new.begin(), pbc_new.end(), pid);
//					if (it == pbc_new.end())
//					{
//						pbc_new.push_back(pid);
//						pdisp_new.push_back(pdisp[i]);
//					}
//				}
//			}
//		}
//		else
//		{
//			for (uint j = 0; j < hcp[pbc[i][0]][pbc[i][1]].chd.size(); j++)
//			{
//				array<int, 2> pid = { pbc[i][0] + 1, hcp[pbc[i][0]][pbc[i][1]].chd[j] };
//				if (hcp[pid[0]][pid[1]].act == 1 && hcp[pid[0]][pid[1]].type == 1)
//				{
//					vector<array<int, 2>>::iterator it = find(pbc_new.begin(), pbc_new.end(), pid);
//					if (it == pbc_new.end())
//					{
//						pbc_new.push_back(pid);
//						pdisp_new.push_back(pdisp[i]);
//					}
//				}
//			}
//		}
//	}
//	pbc.clear();
//	pdisp.clear();
//	pbc = pbc_new;
//	pdisp = pdisp_new;
//}

void TruncatedTspline_3D::SetInitialBC(vector<array<int, 2>>& ebc, vector<double>& edisp, vector<array<int, 2>>& pbc, vector<double>& pdisp)
{
	//list of element
	//cube4
	//int nbc1(1), nbc2(1);
	//int bc1[] = { 20,24,36,40 };
	//int bc2[] = { 23,27,39,43 };
	//cube9
	//int nbc1(4), nbc2(4);
	//int bc1[] = { 679, 680, 688, 689 };
	//int bc2[] = { 40, 41, 49, 50 };
	//rod CAD model
	//int nbc1(10), nbc2(8);
	//int bc1[] = {1249,1250,1251,1252,1253,1257,1258,1259,1260,1261};
	//int bc2[] = { 775, 781, 783, 795, 801, 803, 853, 865 };
	//base CAD model
	//int nbc1(15), nbc2(16);
	//int bc1[] = { 8954,8972,8986,9086,9101,9118,9221,9236,9251,9345,9355,9370,9467,9480,9497 };
	//int bc2[] = { 7867,8095,8128,8212,8261,8301,8417,8426,8462,8495,8496,10181,10283,10364,10372,10439 };
	//cross hole CAD model
	//int nbc1(16), nbc2(16);
	//int bc1[] = { 2394,2395,2450,2451,2397,2396,2453,2452,2534,2535,2478,2479,2537,2536,2481,2480 };
	//int bc2[] = { 943,942,363,362,88,941,6,361,608,607,859,858,40,606,76,857 };
	//head model
	//int nbc1(9), nbc2(40);
	//int bc1[] = { 5579,5580,6077,5581,5582,6079,6281,6282,6553 };
	//int bc2[] = { 2658,2659,2684,2685,2686,4002,4003,4004,2655,2656,2681,2682,2683,3999,4000,4001,2283,2284,2309,2310,2311,3708,3709,3710,
	//	2280,2281,2306,2307,2308,3705,3706,3707,2277,2278,2303,2304,2305,3702,3703,3704 };
	//statue model
	int nbc1(41), nbc2(27);
	int bc1[] = { 11851,11852,11853,12181,12182,11906,11908,11910,12270,12272,11907,11909,11911,12271,12273,11897,11898,11899,12247,12248,
	11894,11887,11900,11901,11902,11903,11904,12249,12250,12252,12253,11888,11891,11892,12237,12238,11885,11889,11890,12235,12236};
	int bc2[] = { 9093,9094,9564,9565,9570,9572,9571,9096,9098,9567,9569,9578,9095,9097,9566,9568,9577,9103,9104,9582,9583,9588,9105,9107,9585,9587,9590 };
	//hook model
	//int nbc1(18), nbc2(29);
	//int bc1[] = { 5075,5077,5079,5074,5076,5078,5069,5070,5071,5059,5065,5067,5058,5064,5066,5054,5061,5063 };
	//int bc2[] = { 3612,3613,3610,3611,3766,3767,3764,3608,3609,3606,3607,3762,3763,3760,3603,3604,3601,3602,3718,3719,3715,3605,3600,3597,3598,3711,3712,3714,3717 };
	//gear model
	//int nbc1(33), nbc2(24);
	//int bc1[] = { 14589,14572,14571,14570,14569,14540,14539,14538,14534,14533,14532,14585,14568,14567,14566,14565,14537,14536,14535,14531,14530,14529,
	//14581,14564,14563,14562,14561,14528,14527,14526,14519,14518,14517};
	//int bc2[] = { 14344,14346,14345,11042,11043,11044,11048,11049,14347,14349,14348,11045,11046,11047,11050,11051,14357,14359,14358,11052,11053,11054,11058,11059 };

	ebc.clear();
	edisp.clear();
	pbc.clear();
	pdisp.clear();
	for (int i = 0; i < nbc1; i++)
	{
		array<int, 2> tmp1 = { 0, bc1[i] };
		ebc.push_back(tmp1);
		edisp.push_back(100.);
		for (int j = 0; j < 8; j++)
		{
			if (hcp[0][hmesh[0][bc1[i]].cnct[j]].type == 1)
			{
				int pid(hmesh[0][bc1[i]].cnct[j]);
				for (uint k = 0; k < hcp[0][pid].face.size(); k++)
				{
					if (hface[0][hcp[0][pid].face[k]].type == 1)
					{
						for (int k1 = 0; k1 < 4; k1++)
						{
							array<int, 2> tmp = { 0, hface[0][hcp[0][pid].face[k]].cnct[k1] };
							vector<array<int, 2>>::iterator it = find(pbc.begin(), pbc.end(), tmp);
							if (it == pbc.end())
							{
								pbc.push_back(tmp);
								pdisp.push_back(100.);
							}
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < nbc2; i++)
	{
		array<int, 2> tmp1 = { 0, bc2[i] };
		ebc.push_back(tmp1);
		edisp.push_back(0.);
		for (int j = 0; j < 8; j++)
		{
			if (hcp[0][hmesh[0][bc2[i]].cnct[j]].type == 1)
			{
				int pid(hmesh[0][bc2[i]].cnct[j]);
				for (uint k = 0; k < hcp[0][pid].face.size(); k++)
				{
					if (hface[0][hcp[0][pid].face[k]].type == 1)
					{
						for (int k1 = 0; k1 < 4; k1++)
						{
							array<int, 2> tmp = { 0, hface[0][hcp[0][pid].face[k]].cnct[k1] };
							vector<array<int, 2>>::iterator it = find(pbc.begin(), pbc.end(), tmp);
							if (it == pbc.end())
							{
								pbc.push_back(tmp);
								pdisp.push_back(0.);
							}
						}
					}
				}
			}
		}
	}
}

void TruncatedTspline_3D::SetBC(vector<array<int, 2>>& ebc, vector<double>& edisp, vector<array<int, 2>>& pbc, vector<double>& pdisp)
{
	vector<array<int, 2>> ebc_new, pbc_new;
	vector<double> edisp_new, pdisp_new;
	for (uint i = 0; i < pbc.size(); i++)
	{
		if (hcp[pbc[i][0]][pbc[i][1]].act == 1)
		{
			pbc_new.push_back(pbc[i]);
			pdisp_new.push_back(pdisp[i]);
		}
	}
	for (uint i = 0; i < ebc.size(); i++)
	{
		if (hmesh[ebc[i][0]][ebc[i][1]].act == 1)
		{
			ebc_new.push_back(ebc[i]);
			edisp_new.push_back(edisp[i]);
		}
		else if (hmesh[ebc[i][0]][ebc[i][1]].act == 0)
		{
			for (uint k2 = 0; k2 < hmesh[ebc[i][0]][ebc[i][1]].chd.size(); k2++)
			{
				array<int, 2> eid = { ebc[i][0] + 1, hmesh[ebc[i][0]][ebc[i][1]].chd[k2] };
				if (hmesh[eid[0]][eid[1]].act == 1 && hmesh[eid[0]][eid[1]].type == 1)
				{
					ebc_new.push_back(eid);
					edisp_new.push_back(edisp[i]);
					for (int j = 0; j < 8; j++)
					{
						array<int, 2> pid1 = { eid[0], hmesh[eid[0]][eid[1]].cnct[j] };
						if (hcp[pid1[0]][pid1[1]].act == 1 && hcp[pid1[0]][pid1[1]].type == 1)
						{
							for (uint k = 0; k < hcp[pid1[0]][pid1[1]].face.size(); k++)
							{
								array<int, 2> fcid = { eid[0], hcp[pid1[0]][pid1[1]].face[k] };
								if (hface[fcid[0]][fcid[1]].type == 1)
								{
									for (int k1 = 0; k1 < 4; k1++)
									{
										array<int, 2> pid2 = { eid[0], hface[fcid[0]][fcid[1]].cnct[k1] };
										if (hcp[pid2[0]][pid2[1]].act == 1 && hcp[pid2[0]][pid2[1]].type == 1)
										{
											vector<array<int, 2>>::iterator it = find(pbc_new.begin(), pbc_new.end(), pid2);
											if (it == pbc_new.end())
											{
												pbc_new.push_back(pid2);
												pdisp_new.push_back(edisp[i]);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	ebc.clear();
	edisp.clear();
	pbc.clear();
	pdisp.clear();
	ebc = ebc_new;
	edisp = edisp_new;
	pbc = pbc_new;
	pdisp = pdisp_new;

	//cout << "pbc:\n";
	//for (uint i = 0; i < pbc.size(); i++)
	//{
	//	cout << pbc[i][1] << " " << pdisp[i] << "\n";
	//}
	//getchar();
}

void TruncatedTspline_3D::FittingBC(vector<int>& IDBC1, vector<double>& gh1)
{
	IDBC1.clear();
	gh1.clear();

	uint i, j, k, k1, k2;
	int loc(0);
	//boundary points
	vector<vector<int>> flag(hcp.size());
	for (i = 0; i < hcp.size(); i++)
	{
		flag[i].resize(hcp[i].size(), 0);
	}
	for (i = 0; i < hmesh.size(); i++)
	{
		for (j = 0; j < hmesh[i].size(); j++)
		{
			if (hmesh[i][j].act == 1 && hmesh[i][j].type == 1)
			{
				if (hmesh[i][j].trun == 0)
				{
					for (k = 0; k < hmesh[i][j].IEN.size(); k++)
					{
						flag[i][hmesh[i][j].IEN[k]] = 1;
					}
				}
				else
				{
					for (k = 0; k < hmesh[i][j].IEN_act.size(); k++)
					{
						flag[hmesh[i][j].IEN_act[k][0]][hmesh[i][j].IEN_act[k][1]] = 1;
					}
				}
			}
		}
	}

	vector<vector<int>> aloc(hcp.size());
	vector<array<double, 3>> cpts;
	loc = 0;
	for (i = 0; i < hcp.size(); i++)
	{
		aloc[i].resize(hcp[i].size(), -1);
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (flag[i][j] == 1)
			{
				aloc[i][j] = loc++;
				array<double, 3> tmp = { hcp[i][j].coor[0], hcp[i][j].coor[1], hcp[i][j].coor[2] };
				cpts.push_back(tmp);
			}
		}
	}
	vector<int> IDBC(loc, -1);
	vector<double> gh(loc, 0.);
	loc = 0;
	int count(0);
	for (i = 0; i < hcp.size(); i++)
	{
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (flag[i][j] == 1)
			{
				if (hcp[i][j].type == 1)
				{
					IDBC[loc] = count++;
				}
				loc++;
			}
		}
	}
	vector<BezierElement3D> bzmesh;
	for (i = 0; i < hmesh.size(); i++)
	{
		for (j = 0; j < hmesh[i].size(); j++)
		{
			if (hmesh[i][j].act == 1 && hmesh[i][j].type==1)
			{
				BezierElement3D bztmp;
				bztmp.prt[0] = i; bztmp.prt[1] = j;
				bztmp.trun = hmesh[i][j].trun;
				if (hmesh[i][j].type == 1) bztmp.type = 1;
				for (k = 0; k < 6; k++)
				{
					if (hface[i][hmesh[i][j].face[k]].type == 1)
					{
						bztmp.bfc.push_back(k);
					}
				}
				if (hmesh[i][j].trun == 0)
				{
					bztmp.IEN.resize(hmesh[i][j].IEN.size());
					bztmp.cmat.resize(hmesh[i][j].IEN.size(), vector<double>(64));
					for (k = 0; k < hmesh[i][j].IEN.size(); k++)
					{
						bztmp.IEN[k] = aloc[i][hmesh[i][j].IEN[k]];
						for (k1 = 0; k1 < 64; k1++)
						{
							bztmp.cmat[k][k1] = hmesh[i][j].bemat[k][k1];
							bztmp.pts[k1][0] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[0];
							bztmp.pts[k1][1] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[1];
							bztmp.pts[k1][2] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[2];
						}
					}
				}
				else
				{
					bztmp.IEN.resize(hmesh[i][j].IEN_act.size());
					bztmp.cmat.resize(hmesh[i][j].IEN_act.size(), vector<double>(64));
					for (k = 0; k < hmesh[i][j].IEN_act.size(); k++)
					{
						int lev(hmesh[i][j].IEN_act[k][0]);
						int pid(hmesh[i][j].IEN_act[k][1]);
						if (aloc[lev][pid] == -1)
						{
							cout << "wrong aloc!\n";
							getchar();
						}
						bztmp.IEN[k] = aloc[lev][pid];
						for (k1 = 0; k1 < 64; k1++)
						{
							bztmp.cmat[k][k1] = 0.;
							for (k2 = 0; k2 < hmesh[i][j].IEN.size(); k2++)
							{
								bztmp.cmat[k][k1] += hmesh[i][j].tmat[k][k2] * hmesh[i][j].bemat[k2][k1];
							}
						}
					}
					for (k = 0; k < hmesh[i][j].IEN.size(); k++)
					{
						for (k1 = 0; k1 < 64; k1++)
						{
							bztmp.pts[k1][0] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[0];
							bztmp.pts[k1][1] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[1];
							bztmp.pts[k1][2] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[2];
						}
					}
				}
				bzmesh.push_back(bztmp);
			}
		}
	}

	LeastSquare ls;
	vector<double> sol;
	ls.SetProblem(IDBC, gh);
	ls.GetEqParameter(dmrg, nmpl, acoef);
	//ls.VisualizeBoundarySurface(bzmesh, cpts, "../io/complex2/rod2");
	//cout << "done output boundary surface!\n";
	//getchar();
	ls.Run_Fitting(bzmesh, "", sol);
	//ls.VisualizeBoundarySurface(bzmesh, cpts, "../io/complex2/rod2");
	//cout << "done output boundary surface!\n";
	//getchar();

	IDBC1.clear();
	gh1.clear();
	int loc1(0);
	loc = 0;
	//int count = 0;
	count = 0;
	for (i = 0; i < hcp.size(); i++)
	{
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				if (flag[i][j] == 1)
				{
					if (hcp[i][j].type == 1)
					{
						IDBC1.push_back(-1);
						gh1.push_back(sol[loc]);
						//double tmp = SpecifyDirichBC(hcp[i][j].coor);
						//gh1.push_back(tmp);
					}
					else
					{
						IDBC1.push_back(count);
						gh1.push_back(0.);
						count++;
					}
					loc++;
				}
				else
				{
					IDBC1.push_back(count);
					gh1.push_back(0.);
					count++;
				}
			}
		}
	}
}

void TruncatedTspline_3D::InputCheck(string fn)
{
	OutputGeom_All(fn);
	cout << "Output Geom done!\n";
	getchar();
}

//void TruncatedTspline_3D::MeshRepair(string fn_in, string fn_out)
//{
//	string fname(fn_in + ".vtk"), stmp;
//	int npts, neles, itmp;
//	ifstream fin;
//	fin.open(fname);
//	if (fin.is_open())
//	{
//		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
//		fin >> stmp >> npts >> stmp;
//		cp.resize(npts);
//		for (int i = 0; i<npts; i++)
//		{
//			cp[i].act = 1;
//			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
//		}
//		getline(fin, stmp);
//		fin >> stmp >> neles >> itmp;
//		tmesh.resize(neles);
//		for (int i = 0; i<neles; i++)
//		{
//			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3] >>
//				tmesh[i].cnct[4] >> tmesh[i].cnct[5] >> tmesh[i].cnct[6] >> tmesh[i].cnct[7];
//		}
//		fin.close();
//	}
//	else
//	{
//		cerr << "Cannot open " << fname << "!\n";
//	}
//	InitialConnect();
//
//	double tol(1.e-6);
//	//vector<int> move(cp.size(),0);
//	vector<array<int,2>> move;
//	for (uint i = 0; i < cp.size(); i++)
//	{
//		if (cp[i].type != 1)
//		{
//			for (uint j = 0; j < cp.size(); j++)
//			{
//				double dis = sqrt((cp[i].coor[0] - cp[j].coor[0])*(cp[i].coor[0] - cp[j].coor[0]) + (cp[i].coor[1] - cp[j].coor[1])*(cp[i].coor[1] - cp[j].coor[1]) +
//					(cp[i].coor[2] - cp[j].coor[2])*(cp[i].coor[2] - cp[j].coor[2]));
//				if (i!=j && dis < tol)
//				{
//					//move[i] = 1; break;
//					array<int, 2> tmp = { i, j };
//					move.push_back(tmp);
//					break;
//				}
//			}
//		}
//	}
//	int fc_dir[6][3] = { { 0, 3, 2 }, { 0, 1, 5 }, { 1, 2, 6 }, { 3, 7, 6 }, { 0, 4, 7 }, { 4, 5, 6 } };
//	for (uint i = 0; i < move.size(); i++)
//	{
//		vector<array<double, 3>> vec;
//		for (uint j = 0; j < cp[move[i][1]].hex.size(); j++)
//		{
//			if (tmesh[cp[move[i][1]].hex[j]].type == 1)
//			{
//				int hxid(cp[move[i][1]].hex[j]);
//				for (int k = 0; k < 6; k++)
//				{
//					if (tmface[tmesh[hxid].face[k]].type == 1)
//					{
//						array<double, 3> tmp1, tmp2;
//						for (int dof = 0; dof < 3; dof++)
//						{
//							tmp1[dof] = cp[tmesh[hxid].cnct[fc_dir[k][1]]].coor[dof] - cp[tmesh[hxid].cnct[fc_dir[k][0]]].coor[dof];
//							tmp2[dof] = cp[tmesh[hxid].cnct[fc_dir[k][2]]].coor[dof] - cp[tmesh[hxid].cnct[fc_dir[k][1]]].coor[dof];
//						}
//						array<double, 3> tmp3 = { tmp1[1] * tmp2[2] - tmp1[2] * tmp2[1], -tmp1[0] * tmp2[2] + tmp1[2] * tmp2[0], tmp1[0] * tmp2[1] - tmp1[1] * tmp2[0] };
//						double dis = sqrt(tmp3[0] * tmp3[0] + tmp3[1] * tmp3[1] + tmp3[2] * tmp3[2]);
//						tmp3[0] /= dis; tmp3[1] /= dis; tmp3[2] /= dis;
//						vec.push_back(tmp3);
//					}
//				}
//			}
//		}
//		array<double, 3> nm = {0.,0.,0.};
//		for (uint j = 0; j < vec.size(); j++)
//		{
//			nm[0] += vec[j][0]; nm[1] += vec[j][1]; nm[2] += vec[j][2];
//		}
//		nm[0] /= vec.size(); nm[1] /= vec.size(); nm[2] /= vec.size();
//		nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
//
//		double lmin(1.e6);
//		for (uint j = 0; j < cp[move[i][0]].edge.size(); j++)
//		{
//			int edid = cp[move[i][0]].edge[j];
//			double dis = sqrt((cp[tmedge[edid].pt[1]].coor[0] - cp[tmedge[edid].pt[0]].coor[0])*(cp[tmedge[edid].pt[1]].coor[0] - cp[tmedge[edid].pt[0]].coor[0]) +
//				(cp[tmedge[edid].pt[1]].coor[1] - cp[tmedge[edid].pt[0]].coor[1])*(cp[tmedge[edid].pt[1]].coor[1] - cp[tmedge[edid].pt[0]].coor[1]) +
//				(cp[tmedge[edid].pt[1]].coor[2] - cp[tmedge[edid].pt[0]].coor[2])*(cp[tmedge[edid].pt[1]].coor[2] - cp[tmedge[edid].pt[0]].coor[2]));
//			if (dis>tol && dis < lmin) lmin = dis;
//		}
//		nm[0] *= 2.*lmin; nm[1] *= 2.*lmin; nm[2] *= 2.*lmin;
//
//		cp[move[i][0]].coor[0] += nm[0]; cp[move[i][0]].coor[1] += nm[1]; cp[move[i][0]].coor[2] += nm[2];
//	}
//
//	//string fn2(fn_out + "_tmp.vtk");
//	//ofstream fout;
//	//fout.open(fn2.c_str());
//	//if (fout.is_open())
//	//{
//	//	fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
//	//	fout << "POINTS " << cp.size() << " float\n";
//	//	for (uint i = 0; i<cp.size(); i++)
//	//	{
//	//		fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
//	//	}
//	//	fout << "\nCELLS " << tmedge.size() << " " << 3 * tmedge.size() << '\n';
//	//	for (uint i = 0; i<tmedge.size(); i++)
//	//	{
//	//		fout << "2 ";
//	//		for (int j = 0; j<2; j++)
//	//		{
//	//			fout << tmedge[i].pt[j] << ' ';
//	//		}
//	//		fout << '\n';
//	//	}
//	//	fout << "\nCELL_TYPES " << tmedge.size() << '\n';
//	//	for (uint i = 0; i<tmedge.size(); i++)
//	//	{
//	//		fout << "3\n";
//	//	}
//	////	fout << "POINT_DATA " << cp.size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
//	////	for (uint i = 0; i<cp.size(); i++)
//	////	{
//	////		fout << move[i] << "\n";
//	////	}
//	////	//fout << "\nCELL_DATA " << hedge[lev].size() << "\nSCALARS eact float 1\nLOOKUP_TABLE default\n";
//	////	//for (uint i = 0; i<hedge[lev].size(); i++)
//	////	//{
//	////	//	fout << hedge[lev][i].sharp << "\n";
//	////	//}
//	//	fout.close();
//	//}
//	//else
//	//{
//	//	cout << "Cannot open " << fn2 << "!\n";
//	//}
//
//	string fn2(fn_out + ".vtk");
//	ofstream fout;
//	fout.open(fn2.c_str());
//	if (fout.is_open())
//	{
//		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
//		fout << "POINTS " << cp.size() << " float\n";
//		for (uint i = 0; i<cp.size(); i++)
//		{
//			fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
//		}
//		fout << "\nCELLS " << tmesh.size() << " " << 9 * tmesh.size() << '\n';
//		for (uint i = 0; i<tmesh.size(); i++)
//		{
//			fout << "8 ";
//			for (int j = 0; j<8; j++)
//			{
//				fout << tmesh[i].cnct[j] << ' ';
//			}
//			fout << '\n';
//		}
//		fout << "\nCELL_TYPES " << tmesh.size() << '\n';
//		for (uint i = 0; i<tmesh.size(); i++)
//		{
//			fout << "12\n";
//		}
//		//fout << "POINT_DATA " << hcp[lev].size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
//		//for (uint i = 0; i<hcp[lev].size(); i++)
//		//{
//		//	fout << hcp[lev][i].act << "\n";
//		//}
//		//fout << "\nCELL_DATA " << hmesh[lev].size() << "\nSCALARS eact float 1\nLOOKUP_TABLE default\n";
//		//for (uint i = 0; i<hmesh[lev].size(); i++)
//		//{
//		//	//fout<<hmesh[lev][i].act<<"\n";
//		//	fout << hmesh[lev][i].type << "\n";
//		//}
//		fout.close();
//	}
//	else
//	{
//		cout << "Cannot open " << fn2 << "!\n";
//	}
//
//	cout << "Done mesh repair!\n";
//	getchar();
//}

void TruncatedTspline_3D::MeshRepair(string fn_in, string fn_out)
{
	string fname(fn_in + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i<npts; i++)
		{
			cp[i].act = 1;
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i<neles; i++)
		{
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3] >>
				tmesh[i].cnct[4] >> tmesh[i].cnct[5] >> tmesh[i].cnct[6] >> tmesh[i].cnct[7];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}
	InitialConnect();

	int npr(8);
	int pr[] = {6180,6181,6182,6183,6169,6168,6171,6170};//cross hole model, point ids for bad elements
	vector<int> ebad;
	for (int i = 0; i < npr; i++)
	{
		for (uint j = 0; j < cp[pr[i]].hex.size(); j++)
		{
			int hxid(cp[pr[i]].hex[j]);
			if (tmesh[hxid].type == 1)
			{
				vector<int>::iterator it = find(ebad.begin(), ebad.end(), hxid);
				if (it == ebad.end())
				{
					ebad.push_back(hxid);
				}
			}
		}
	}
	double tol(1.e-6);
	vector<array<int, 2>> move;
	int pair[6][8] = { { 0, 1, 2, 3, 4, 5, 6, 7 }, { 0, 1, 5, 4, 3, 2, 6, 7 }, { 1, 2, 6, 5, 0, 3, 7, 4 }, { 3, 2, 6, 7, 0, 1, 5, 4 },
	{ 0, 3, 7, 4, 1, 2, 6, 5 }, { 4, 5, 6, 7, 0, 1, 2, 3 } };
	for (uint i = 0; i < ebad.size(); i++)
	{
		int fc_loc(-1);
		for (int j = 0; j < 6; j++)
		{
			if (tmface[tmesh[ebad[i]].face[j]].type == 1)
			{
				fc_loc = j; break;
			}
		}
		if (fc_loc != -1)
		{
			for (int j = 0; j < 4; j++)
			{
				array<int,2> ptmp = { tmesh[ebad[i]].cnct[pair[fc_loc][j + 4]], tmesh[ebad[i]].cnct[pair[fc_loc][j]] };
				vector<array<int, 2>>::iterator it = find(move.begin(), move.end(), ptmp);
				if (it == move.end())
				{
					move.push_back(ptmp);
				}
			}
		}
	}

	//for (uint i = 0; i < cp.size(); i++)
	//{
	//	if (cp[i].type != 1)
	//	{
	//		for (uint j = 0; j < cp.size(); j++)
	//		{
	//			double dis = sqrt((cp[i].coor[0] - cp[j].coor[0])*(cp[i].coor[0] - cp[j].coor[0]) + (cp[i].coor[1] - cp[j].coor[1])*(cp[i].coor[1] - cp[j].coor[1]) +
	//				(cp[i].coor[2] - cp[j].coor[2])*(cp[i].coor[2] - cp[j].coor[2]));
	//			if (i != j && dis < tol)
	//			{
	//				//move[i] = 1; break;
	//				array<int, 2> tmp = { i, j };
	//				move.push_back(tmp);
	//				break;
	//			}
	//		}
	//	}
	//}
	int fc_dir[6][3] = { { 0, 3, 2 }, { 0, 1, 5 }, { 1, 2, 6 }, { 3, 7, 6 }, { 0, 4, 7 }, { 4, 5, 6 } };
	for (uint i = 0; i < move.size(); i++)
	{
		vector<array<double, 3>> vec;
		for (uint j = 0; j < cp[move[i][1]].hex.size(); j++)
		{
			if (tmesh[cp[move[i][1]].hex[j]].type == 1)
			{
				int hxid(cp[move[i][1]].hex[j]);
				for (int k = 0; k < 6; k++)
				{
					if (tmface[tmesh[hxid].face[k]].type == 1)
					{
						array<double, 3> tmp1, tmp2;
						for (int dof = 0; dof < 3; dof++)
						{
							tmp1[dof] = cp[tmesh[hxid].cnct[fc_dir[k][1]]].coor[dof] - cp[tmesh[hxid].cnct[fc_dir[k][0]]].coor[dof];
							tmp2[dof] = cp[tmesh[hxid].cnct[fc_dir[k][2]]].coor[dof] - cp[tmesh[hxid].cnct[fc_dir[k][1]]].coor[dof];
						}
						array<double, 3> tmp3 = { tmp1[1] * tmp2[2] - tmp1[2] * tmp2[1], -tmp1[0] * tmp2[2] + tmp1[2] * tmp2[0], tmp1[0] * tmp2[1] - tmp1[1] * tmp2[0] };
						double dis = sqrt(tmp3[0] * tmp3[0] + tmp3[1] * tmp3[1] + tmp3[2] * tmp3[2]);
						tmp3[0] /= dis; tmp3[1] /= dis; tmp3[2] /= dis;
						vec.push_back(tmp3);
					}
				}
			}
		}
		array<double, 3> nm = { 0., 0., 0. };
		for (uint j = 0; j < vec.size(); j++)
		{
			nm[0] += vec[j][0]; nm[1] += vec[j][1]; nm[2] += vec[j][2];
		}
		nm[0] /= vec.size(); nm[1] /= vec.size(); nm[2] /= vec.size();
		nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];

		double lmin(1.e6);
		for (uint j = 0; j < cp[move[i][0]].edge.size(); j++)
		{
			int edid = cp[move[i][0]].edge[j];
			double dis = sqrt((cp[tmedge[edid].pt[1]].coor[0] - cp[tmedge[edid].pt[0]].coor[0])*(cp[tmedge[edid].pt[1]].coor[0] - cp[tmedge[edid].pt[0]].coor[0]) +
				(cp[tmedge[edid].pt[1]].coor[1] - cp[tmedge[edid].pt[0]].coor[1])*(cp[tmedge[edid].pt[1]].coor[1] - cp[tmedge[edid].pt[0]].coor[1]) +
				(cp[tmedge[edid].pt[1]].coor[2] - cp[tmedge[edid].pt[0]].coor[2])*(cp[tmedge[edid].pt[1]].coor[2] - cp[tmedge[edid].pt[0]].coor[2]));
			if (dis>tol && dis < lmin) lmin = dis;
		}
		double c_amp(2.5);
		nm[0] *= c_amp*lmin; nm[1] *= c_amp*lmin; nm[2] *= c_amp*lmin;

		cp[move[i][0]].coor[0] = cp[move[i][1]].coor[0] + nm[0];
		cp[move[i][0]].coor[1] = cp[move[i][1]].coor[1] + nm[1];
		cp[move[i][0]].coor[2] = cp[move[i][1]].coor[2] + nm[2];
	}

	//string fn2(fn_out + "_tmp.vtk");
	//ofstream fout;
	//fout.open(fn2.c_str());
	//if (fout.is_open())
	//{
	//	fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	//	fout << "POINTS " << cp.size() << " float\n";
	//	for (uint i = 0; i<cp.size(); i++)
	//	{
	//		fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
	//	}
	//	fout << "\nCELLS " << tmedge.size() << " " << 3 * tmedge.size() << '\n';
	//	for (uint i = 0; i<tmedge.size(); i++)
	//	{
	//		fout << "2 ";
	//		for (int j = 0; j<2; j++)
	//		{
	//			fout << tmedge[i].pt[j] << ' ';
	//		}
	//		fout << '\n';
	//	}
	//	fout << "\nCELL_TYPES " << tmedge.size() << '\n';
	//	for (uint i = 0; i<tmedge.size(); i++)
	//	{
	//		fout << "3\n";
	//	}
	////	fout << "POINT_DATA " << cp.size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
	////	for (uint i = 0; i<cp.size(); i++)
	////	{
	////		fout << move[i] << "\n";
	////	}
	////	//fout << "\nCELL_DATA " << hedge[lev].size() << "\nSCALARS eact float 1\nLOOKUP_TABLE default\n";
	////	//for (uint i = 0; i<hedge[lev].size(); i++)
	////	//{
	////	//	fout << hedge[lev][i].sharp << "\n";
	////	//}
	//	fout.close();
	//}
	//else
	//{
	//	cout << "Cannot open " << fn2 << "!\n";
	//}

	string fn2(fn_out + ".vtk");
	ofstream fout;
	fout.open(fn2.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << cp.size() << " float\n";
		for (uint i = 0; i<cp.size(); i++)
		{
			fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << tmesh.size() << " " << 9 * tmesh.size() << '\n';
		for (uint i = 0; i<tmesh.size(); i++)
		{
			fout << "8 ";
			for (int j = 0; j<8; j++)
			{
				fout << tmesh[i].cnct[j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << tmesh.size() << '\n';
		for (uint i = 0; i<tmesh.size(); i++)
		{
			fout << "12\n";
		}
		//fout << "POINT_DATA " << hcp[lev].size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<hcp[lev].size(); i++)
		//{
		//	fout << hcp[lev][i].act << "\n";
		//}
		//fout << "\nCELL_DATA " << hmesh[lev].size() << "\nSCALARS eact float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<hmesh[lev].size(); i++)
		//{
		//	//fout<<hmesh[lev][i].act<<"\n";
		//	fout << hmesh[lev][i].type << "\n";
		//}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fn2 << "!\n";
	}

	cout << "Done mesh repair!\n";
	getchar();
}

void TruncatedTspline_3D::ReportXP()
{
	vector<int> pflag(cp.size(),0);
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].type == 3)//interior extraordinary points
		{
			pflag[i] = 1;
			count++;
		}
		else if (cp[i].type == 1)
		{
			int nfc(0);
			for (uint j = 0; j < cp[i].face.size(); j++)
			{
				if (tmface[cp[i].face[j]].type == 1) nfc++;
			}
			if (nfc != 4)
			{
				pflag[i] = 1;
				count++;
			}
		}
	}
	cout << "# of XP: " << count << "\n";
	count = 0;
	for (uint i = 0; i < tmesh.size(); i++)
	{
		//if (tmesh[i].type == 2) count++;
		int flag(0);
		for (int j = 0; j < 8; j++)
		{
			if (pflag[tmesh[i].cnct[j]] == 1)
			{
				flag = 1; break;
			}
		}
		if (flag == 1) count++;
	}
	cout << "# of XE: " << count << "\n";
	//getchar();
}

void TruncatedTspline_3D::OutputRefineID(string fn, const vector<array<int, 2>>& rfid, const vector<array<int, 2>>& gst)
{
	string fname = fn + "_rfid.txt";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << rfid.size() << "\n";
		for (uint i = 0; i < rfid.size(); i++)
		{
			fout << rfid[i][0] << " " << rfid[i][1] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	string fn1 = fn + "_gst.txt";
	fout.open(fn1.c_str());
	if (fout.is_open())
	{
		fout << gst.size() << "\n";
		for (uint i = 0; i < gst.size(); i++)
		{
			fout << gst[i][0] << " " << gst[i][1] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fn1 << "!\n";
	}
}

void TruncatedTspline_3D::InputRefineID(string fn, vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst)
{
	int tmp;
	rfid.clear();
	gst.clear();
	string fname = fn + "_rfid.txt";
	ifstream fin;
	fin.open(fname.c_str());
	if (fin.is_open())
	{
		fin >> tmp;
		rfid.resize(tmp);
		for (uint i = 0; i < rfid.size(); i++)
		{
			fin >> rfid[i][0] >> rfid[i][1];
		}
		fin.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	string fn1 = fn + "_gst.txt";
	fin.open(fn1.c_str());
	if (fin.is_open())
	{
		fin >> tmp;
		gst.resize(tmp);
		for (uint i = 0; i < gst.size(); i++)
		{
			fin >> gst[i][0] >> gst[i][1];
		}
		fin.close();
	}
	else
	{
		cout << "Cannot open " << fn1 << "!\n";
	}
}

void TruncatedTspline_3D::GetRemoveRegion(double xy[3][2])
{
	double tmp(1.e6);
	double x_range[3][2] = { { tmp, -tmp }, { tmp, -tmp }, { tmp, -tmp } };
	uint i, j;
	for (i = 0; i <cp.size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (cp[i].coor[j] < x_range[j][0]) x_range[j][0] = cp[i].coor[j];
			if (cp[i].coor[j] > x_range[j][1]) x_range[j][1] = cp[i].coor[j];
		}
	}
	double orig[3] = { (x_range[0][0] + x_range[0][1]) / 2., (x_range[1][0] + x_range[1][1]) / 2., (x_range[2][0] + x_range[2][1]) / 2. };
	double xh[3] = { x_range[0][1] - x_range[0][0], x_range[1][1] - x_range[1][0], x_range[2][1] - x_range[2][0] };
	int dirlg(0);
	double hmax(0.);
	for (i = 0; i < 3; i++)
	{
		if (xh[i]>hmax)
		{
			hmax = xh[i];
			dirlg = i;
		}
	}
	//rod
	//if (dirlg == 0)
	//{
	//	cout << "x is the longest!\n";
	//	xy[0][0] = x_range[0][0]; xy[0][1] = x_range[0][1];
	//	xy[1][0] = orig[1]; xy[1][1] = x_range[1][1];
	//	xy[2][0] = orig[2]; xy[2][1] = x_range[2][1];
	//	//cout << x_range[0][0] << " " << x_range[0][1] << " " << xy[0][0] << " " << xy[0][1] << "\n";
	//	//cout << x_range[1][0] << " " << x_range[1][1] << " " << xy[1][0] << " " << xy[1][1] << "\n";
	//	//cout << x_range[2][0] << " " << x_range[2][1] << " " << xy[2][0] << " " << xy[2][1] << "\n";
	//	//getchar();
	//}
	//else if (dirlg == 1)
	//{
	//	cout << "y is the longest!\n";
	//	xy[0][0] = orig[0]; xy[0][1] = x_range[0][1];
	//	xy[1][0] = x_range[1][0]; xy[1][1] = x_range[1][1];
	//	xy[2][0] = orig[2]; xy[2][1] = x_range[2][1];
	//}
	//else
	//{
	//	cout << "z is the longest!\n";
	//	xy[0][0] = orig[0]; xy[0][1] = x_range[0][1];
	//	xy[1][0] = orig[1]; xy[1][1] = x_range[1][1];
	//	xy[2][0] = x_range[2][0]; xy[2][1] = x_range[2][1];
	//}

	////hook
	//xy[0][0] = x_range[0][0]; xy[0][1] = (x_range[0][0]+orig[0])/2.;
	//xy[1][0] = x_range[1][0]; xy[1][1] = x_range[1][1];
	//xy[2][0] = x_range[2][0]; xy[2][1] = x_range[2][1];

	//base
	//xy[0][0] = x_range[0][0]; xy[0][1] = orig[0];
	//xy[1][0] = x_range[1][0]; xy[1][1] = orig[1];
	//xy[2][0] = x_range[2][0]; xy[2][1] = x_range[2][1];

	//gear
	//xy[0][0] = x_range[0][0]; xy[0][1] = 0.9*orig[0] + 0.1*x_range[0][0];
	//xy[1][0] = x_range[1][0]; xy[1][1] = 0.8*orig[1] + 0.2*x_range[1][1];
	//xy[2][0] = 0.8*orig[2] + 0.2*x_range[2][1]; xy[2][1] = x_range[2][1];

	//statue
	//xy[0][0] = 0.9*orig[0]+0.1*x_range[0][1]; xy[0][1] = x_range[0][1];
	//xy[1][0] = x_range[1][0]; xy[1][1] = 0.7*orig[1] + 0.3*x_range[1][1];
	//xy[2][0] = 0.55*orig[2] + 0.45*x_range[2][1]; xy[2][1] = x_range[2][1];

	xy[0][0] = orig[0]; xy[0][1] = x_range[0][1];
	xy[1][0] = x_range[1][0]; xy[1][1] = orig[1];
	xy[2][0] = orig[2]; xy[2][1] = x_range[2][1];
}

void TruncatedTspline_3D::OutputRemoveCM(string fn, double xy[3][2])
{
	int lev(0);
	string fname(fn + "_remv_CM.vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << hcp[lev].size() << " float\n";
		for (uint i = 0; i<hcp[lev].size(); i++)
		{
			fout << hcp[lev][i].coor[0] << " " << hcp[lev][i].coor[1] << " " << hcp[lev][i].coor[2] << "\n";
		}
		vector<int> outid;
		for (uint i = 0; i < hmesh[lev].size(); i++)
		{
			array<double, 3> ctp = { 0., 0., 0. };
			for (int j = 0; j < 8; j++)
			{
				ctp[0] += hcp[lev][hmesh[lev][i].cnct[j]].coor[0];
				ctp[1] += hcp[lev][hmesh[lev][i].cnct[j]].coor[1];
				ctp[2] += hcp[lev][hmesh[lev][i].cnct[j]].coor[2];
			}
			ctp[0] /= 8.; ctp[1] /= 8.; ctp[2] /= 8.;
			if (ctp[0]>=xy[0][0] && ctp[0]<=xy[0][1] && ctp[1]>=xy[1][0] && ctp[1]<=xy[1][1] && ctp[2]>=xy[2][0] && ctp[2] <= xy[2][1])
			{
				continue;
			}
			else
			{
				outid.push_back(i);
			}
		}

		fout << "\nCELLS " << outid.size() << " " << 9 * outid.size() << '\n';
		for (uint i = 0; i<outid.size(); i++)
		{
			fout << "8 ";
			for (int j = 0; j<8; j++)
			{
				fout << hmesh[lev][outid[i]].cnct[j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << outid.size() << '\n';
		for (uint i = 0; i<outid.size(); i++)
		{
			fout << "12\n";
		}

		//fout << "POINT_DATA " << hcp[lev].size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<hcp[lev].size(); i++)
		//{
		//	fout << hcp[lev][i].act << "\n";
		//}
		//fout << "\nCELL_DATA " << hmesh[lev].size() << "\nSCALARS eact float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<hmesh[lev].size(); i++)
		//{
		//	//fout<<hmesh[lev][i].act<<"\n";
		//	fout << hmesh[lev][i].type << "\n";
		//}

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline_3D::CreateRemoveView(string fn_in, string fn_out)
{
	vector<array<double, 3>> pts;
	vector<double> val;
	vector<array<int, 8>> cube;
	vector<array<double, 3>> lpt;
	vector<array<int, 2>> line;

	string fname(fn_in + "_disp.vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		pts.resize(npts);
		val.resize(npts);
		for (int i = 0; i<npts; i++)
		{
			fin >> pts[i][0] >> pts[i][1] >> pts[i][2];
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		cube.resize(neles);
		for (int i = 0; i<neles; i++)
		{
			fin >> itmp >> cube[i][0] >> cube[i][1] >> cube[i][2] >> cube[i][3] >>
				cube[i][4] >> cube[i][5] >> cube[i][6] >> cube[i][7];
		}
		for (int i = 0; i < 3; i++)
		{
			getline(fin, stmp);//skip lines
			//cout << stmp << "\n";
			//getchar();
		}
		for (int i = 0; i<neles; i++)//skip cell types
		{
			fin >> itmp;
		}
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		for (int i = 0; i<npts; i++)
		{
			fin >> val[i];
			//cout << val[i] << "\n";
			//getchar();
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}
	string fn1(fn_in + "-lines.vtk");
	fin.open(fn1);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		lpt.resize(npts);
		for (int i = 0; i<npts; i++)
		{
			fin >> lpt[i][0] >> lpt[i][1] >> lpt[i][2];
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		line.resize(neles);
		for (int i = 0; i<neles; i++)
		{
			fin >> itmp >> line[i][0] >> line[i][1];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fn1 << "!\n";
	}

	double tmp(1.e6);
	double x_range[3][2] = { { tmp, -tmp }, { tmp, -tmp }, { tmp, -tmp } };
	uint i, j;
	for (i = 0; i <lpt.size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (lpt[i][j] < x_range[j][0]) x_range[j][0] = lpt[i][j];
			if (lpt[i][j] > x_range[j][1]) x_range[j][1] = lpt[i][j];
		}
	}
	double orig[3] = { (x_range[0][0] + x_range[0][1]) / 2., (x_range[1][0] + x_range[1][1]) / 2., (x_range[2][0] + x_range[2][1]) / 2. };
	double xh[3] = { x_range[0][1] - x_range[0][0], x_range[1][1] - x_range[1][0], x_range[2][1] - x_range[2][0] };
	int dirlg(0);
	double hmax(0.), dis;
	for (i = 0; i < 3; i++)
	{
		if (xh[i]>hmax)
		{
			hmax = xh[i];
			dirlg = i;
		}
	}
	double hmin = 1.e6;
	for (i = 0; i < line.size(); i++)
	{
		dis = (lpt[line[i][0]][0] - lpt[line[i][1]][0])*(lpt[line[i][0]][0] - lpt[line[i][1]][0]) +
			(lpt[line[i][0]][1] - lpt[line[i][1]][0])*(lpt[line[i][0]][1] - lpt[line[i][1]][1]) +
			(lpt[line[i][0]][2] - lpt[line[i][1]][0])*(lpt[line[i][0]][2] - lpt[line[i][1]][2]);
		if (dis < hmin) hmin = dis;
	}
	hmin = sqrt(hmin);
	double rmv[3][2], rmv1[3][2];
	//rod
	rmv[0][0] = x_range[0][0]; rmv[0][1] = x_range[0][1];
	rmv[1][0] = orig[1]; rmv[1][1] = x_range[1][1];
	rmv[2][0] = orig[2]; rmv[2][1] = x_range[2][1];

	rmv1[0][0] = x_range[0][0]; rmv1[0][1] = x_range[0][1];
	rmv1[1][0] = orig[1]+2.*hmin; rmv1[1][1] = x_range[1][1];
	rmv1[2][0] = orig[2] + 2.5*hmin; rmv1[2][1] = x_range[2][1];
	//hook
	//rmv[0][0] = x_range[0][0]; rmv[0][1] = (x_range[0][0] + orig[0]) / 2.;
	//rmv[1][0] = x_range[1][0]; rmv[1][1] = x_range[1][1];
	//rmv[2][0] = x_range[2][0]; rmv[2][1] = x_range[2][1];

	vector<int> ptsflag(pts.size(), 0);
	vector<int> ptsloc(pts.size(), -1);
	vector<double> val1;
	vector<array<int, 8>> cube1;
	vector<int> lptflag(lpt.size(), 0);
	vector<int> lptloc(lpt.size(), -1);
	vector<array<int, 2>> line1;
	for (i = 0; i < cube.size(); i++)
	{
		double ctp[3] = { 0., 0., 0. };
		for (j = 0; j < 8; j++)
		{
			ctp[0] += pts[cube[i][j]][0];
			ctp[1] += pts[cube[i][j]][1];
			ctp[2] += pts[cube[i][j]][2];
		}
		ctp[0] /= 8.; ctp[1] /= 8.; ctp[2] /= 8.;
		if (ctp[0] > rmv[0][0] && ctp[0] < rmv[0][1] && ctp[1] > rmv[1][0] && ctp[1] < rmv[1][1] && ctp[2] > rmv[2][0] && ctp[2] < rmv[2][1])
		{
			continue;
		}
		else
		{
			for (j = 0; j < 8; j++)
			{
				ptsflag[cube[i][j]] = 1;
			}
			cube1.push_back(cube[i]);
		}
	}
	int loc(0);
	for (i = 0; i < ptsflag.size(); i++)
	{
		if (ptsflag[i] == 1)
		{
			ptsloc[i] = loc++;
			val1.push_back(val[i]);
		}
	}
	for (i = 0; i < cube1.size(); i++)
	{
		for (j = 0; j < 8; j++) cube1[i][j] = ptsloc[cube1[i][j]];
	}

	for (i = 0; i < line.size(); i++)
	{
		double ctp[3] = { (lpt[line[i][0]][0] + lpt[line[i][1]][0]) / 2., (lpt[line[i][0]][1] + lpt[line[i][1]][1]) / 2., (lpt[line[i][0]][2] + lpt[line[i][1]][2]) / 2.};
		if (ctp[0] > rmv[0][0] && ctp[0] < rmv[0][1] && ctp[1] > rmv[1][0] && ctp[1] < rmv[1][1] && ctp[2] > rmv[2][0] && ctp[2] < rmv[2][1])
		{
			continue;
		}
		else
		{
			lptflag[line[i][0]] = 1; lptflag[line[i][1]] = 1;
			line1.push_back(line[i]);
		}
	}
	int loc1(0);
	for (i = 0; i < lptflag.size(); i++)
	{
		if (lptflag[i] == 1)
		{
			lptloc[i] = loc1++;
		}
	}
	for (i = 0; i < line1.size(); i++)
	{
		for (j = 0; j < 2; j++) line1[i][j] = lptloc[line1[i][j]];
	}
	
	string fname1 = fn_out + "_rmv_disp.vtk";
	ofstream fout;
	fout.open(fname1.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << val1.size() << " float\n";
		for (i = 0; i<pts.size(); i++)
		{
			if (ptsflag[i] == 1)
			{
				fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
			}
		}
		fout << "\nCELLS " << cube1.size() << " " << 9 * cube1.size() << '\n';
		for (i = 0; i<cube1.size(); i++)
		{
			fout << "8 " << cube1[i][0] << " " << cube1[i][1] << " " << cube1[i][2] << " " << cube1[i][3]
				<< " " << cube1[i][4] << " " << cube1[i][5] << " " << cube1[i][6] << " " << cube1[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << cube1.size() << '\n';
		for (i = 0; i<cube1.size(); i++)
		{
			fout << "12\n";
		}
		//fout<<"\nPOINT_DATA "<<sdisp.size()<<"\nVECTORS disp float\n";
		//for(i=0;i<sdisp.size();i++)
		//{
		//	fout << sdisp[i][0] << " " << sdisp[i][1] << " " << sdisp[i][2] << "\n";
		//}
		//fout << "\nPOINT_DATA " << sse.size() << "\nVECTORS strain float\n";
		//for (i = 0; i<sse.size(); i++)
		//{
		//	fout << sse[i][0] << " " << sse[i][1] << " " << sse[i][2] << "\n";
		//}
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

		fout << "\nPOINT_DATA " << val1.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<val1.size(); i++)
		{
			fout << val1[i] << "\n";
		}

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
		cout << "Cannot open " << fname1 << "!\n";
	}

	string fname2(fn_out + "_rmv-lines.vtk");
	ofstream fout1;
	fout1.open(fname2.c_str());
	if (fout1.is_open())
	{
		fout1 << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout1 << "POINTS " << loc1 << " float\n";
		for (i = 0; i<lpt.size(); i++)
		{
			if (lptflag[i] == 1)
			{
				fout1 << lpt[i][0] << " " << lpt[i][1] << " " << lpt[i][2] << "\n";
			}
		}
		fout1 << "\nCELLS " << line1.size() << " " << 3 * line1.size() << '\n';
		for (i = 0; i<line1.size(); i++)
		{
			fout1 << "2 " << line1[i][0] << " " << line1[i][1] << '\n';
		}
		fout1 << "\nCELL_TYPES " << line1.size() << '\n';
		for (i = 0; i<line1.size(); i++)
		{
			fout1 << "3\n";
		}
		fout1.close();
	}
	else
	{
		cout << "Cannot open " << fname2 << "!\n";
	}
}

void TruncatedTspline_3D::GlobalRefine(int niter)
{
	for (int it = 0; it < niter; it++)
	{
		vector<array<int, 2>> rfid, gst;
		for (uint i = 0; i < hmesh.size(); i++)
		{
			for (uint j = 0; j < hmesh[i].size(); j++)
			{
				if (hmesh[i][j].act == 1)
				{
					array<int, 2> tmp = { i, j };
					rfid.push_back(tmp);
				}
			}
		}
		Refine(rfid, gst);
		//stringstream ss;
		//ss << it + 1;
		//OutputCM(hmesh.size()-1, "../io/hex_input/cube_coarse_"+ss.str());
	}
	cout << "Global refinement done!\n";
	//getchar();
}

double TruncatedTspline_3D::MaxElementSize(const vector<BezierElement3D>& bzmesh)
{
	int ids[4][2] = { { 0, 63 }, { 3, 60 }, { 15, 48 }, { 12, 51 } };
	double max(0.), tmp;
	for (uint i = 0; i < bzmesh.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			tmp = (bzmesh[i].pts[ids[j][0]][0] - bzmesh[i].pts[ids[j][1]][0])*(bzmesh[i].pts[ids[j][0]][0] - bzmesh[i].pts[ids[j][1]][0])
				+ (bzmesh[i].pts[ids[j][0]][1] - bzmesh[i].pts[ids[j][1]][1])*(bzmesh[i].pts[ids[j][0]][1] - bzmesh[i].pts[ids[j][1]][1])
				+ (bzmesh[i].pts[ids[j][0]][2] - bzmesh[i].pts[ids[j][1]][2])*(bzmesh[i].pts[ids[j][0]][2] - bzmesh[i].pts[ids[j][1]][2]);
			max = tmp>max ? tmp : max;
		}
	}
	return sqrt(max);
}

void TruncatedTspline_3D::InitializeMesh(string fn)
{
	//read hex vtk
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i<npts; i++)
		{
			cp[i].act = 1;
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i<neles; i++)
		{
			tmesh[i].act = 1;
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3] >>
				tmesh[i].cnct[4] >> tmesh[i].cnct[5] >> tmesh[i].cnct[6] >> tmesh[i].cnct[7];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	//SetDomain();//normalize coordinates in [0,1]^3, for the cube poisson

	InitialConnect();
	SetSharpFeature_1();

	//string fn1("../io/complex4/hook1");
	//OutputEdge(fn1);
	//OutputCM(fn1);
	//cout << "done setting sharp feature\n";
	//getchar();
}

void TruncatedTspline_3D::SetSharpFeature_1()
{
	//sharp edge and sharp corner, indicated by control point
	//sharp edge
	//double tol(.3);
	double tol(.57);//rod, base
	//double tol(.8);//hook
	int fc_dir[6][3] = { { 0, 3, 2 }, { 0, 1, 5 }, { 1, 2, 6 }, { 3, 7, 6 }, { 0, 4, 7 }, { 4, 5, 6 } };
	uint i, j, k;
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].type == 1)
		{
			vector<array<double, 3>> vec;
			for (j = 0; j < tmedge[i].hex.size(); j++)
			{
				if (tmesh[tmedge[i].hex[j]].type == 1)
				{
					int hxid(tmedge[i].hex[j]);
					for (k = 0; k < 6; k++)
					{
						if (tmface[tmesh[hxid].face[k]].type == 1)
						{
							vector<int>::iterator it = find(tmedge[i].face.begin(), tmedge[i].face.end(), tmesh[hxid].face[k]);
							if (it != tmedge[i].face.end())
							{
								array<double, 3> tmp1, tmp2;
								for (int dof = 0; dof < 3; dof++)
								{
									tmp1[dof] = cp[tmesh[hxid].cnct[fc_dir[k][1]]].coor[dof] - cp[tmesh[hxid].cnct[fc_dir[k][0]]].coor[dof];
									tmp2[dof] = cp[tmesh[hxid].cnct[fc_dir[k][2]]].coor[dof] - cp[tmesh[hxid].cnct[fc_dir[k][1]]].coor[dof];
								}
								array<double, 3> tmp3 = { tmp1[1] * tmp2[2] - tmp1[2] * tmp2[1], -tmp1[0] * tmp2[2] + tmp1[2] * tmp2[0], tmp1[0] * tmp2[1] - tmp1[1] * tmp2[0] };
								double dis = sqrt(tmp3[0] * tmp3[0] + tmp3[1] * tmp3[1] + tmp3[2] * tmp3[2]);
								tmp3[0] /= dis; tmp3[1] /= dis; tmp3[2] /= dis;
								vec.push_back(tmp3);
							}
						}
					}
				}
			}
			if (vec.size() == 2)
			{
				double ang(vec[0][0] * vec[1][0] + vec[0][1] * vec[1][1] + vec[0][2] * vec[1][2]);
				if (ang < tol)
				{
					tmedge[i].sharp = 1;
					cp[tmedge[i].pt[0]].sharp = 1;
					cp[tmedge[i].pt[1]].sharp = 1;
				}
			}
			else
			{
				cerr << "Something wrong in determining sharp edge!\n";
				getchar();
			}
		}
	}

	for (i = 0; i < cp.size(); i++)
	{
		if (cp[i].sharp == 1)
		{
			int nshp(0);
			for (j = 0; j < cp[i].edge.size(); j++)
			{
				if (tmedge[cp[i].edge[j]].sharp == 1)
				{
					nshp++;
				}
			}
			if (nshp >= 3) cp[i].sharp = 2;//sharp corner
			if (nshp == 1) cp[i].sharp = 0;
		}
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].sharp == 1)
		{
			if (cp[tmedge[i].pt[0]].sharp == 0 || cp[tmedge[i].pt[1]].sharp == 0)
			{
				tmedge[i].sharp = 0;
			}
		}
	}
}

void TruncatedTspline_3D::Global_Subdivide()
{
	//body points
	vector<array<double, 3>> pts_bd(tmesh.size());
	for (int i = 0; i < tmesh.size(); i++)
	{
		pts_bd[i][0] = 0.; pts_bd[i][1] = 0.; pts_bd[i][2] = 0.;
		for (int j = 0; j < 8; j++)
		{
			pts_bd[i][0] += cp[tmesh[i].cnct[j]].coor[0];
			pts_bd[i][1] += cp[tmesh[i].cnct[j]].coor[1];
			pts_bd[i][2] += cp[tmesh[i].cnct[j]].coor[2];
		}
		pts_bd[i][0] /= 8.; pts_bd[i][1] /= 8.; pts_bd[i][2] /= 8.;
	}
	//face points
	vector<array<double, 3>> pts_fc0(tmface.size());//simply average
	vector<array<double, 3>> pts_fc(tmface.size());//desired
	for (int i = 0; i < tmface.size(); i++)
	{
		pts_fc0[i][0] = 0.; pts_fc0[i][1] = 0.; pts_fc0[i][2] = 0.;
		for (int j = 0; j < 4; j++)
		{
			pts_fc0[i][0] += cp[tmface[i].cnct[j]].coor[0];
			pts_fc0[i][1] += cp[tmface[i].cnct[j]].coor[1];
			pts_fc0[i][2] += cp[tmface[i].cnct[j]].coor[2];
		}
		pts_fc0[i][0] /= 4.; pts_fc0[i][1] /= 4.; pts_fc0[i][2] /= 4.;
		if (tmface[i].hex.size() == 2)//interior face points
		{
			int hxid[2] = { tmface[i].hex[0], tmface[i].hex[1] };
			pts_fc[i][0] = (pts_bd[hxid[0]][0] + pts_bd[hxid[1]][0] + 2.*pts_fc0[i][0]) / 4.;
			pts_fc[i][1] = (pts_bd[hxid[0]][1] + pts_bd[hxid[1]][1] + 2.*pts_fc0[i][1]) / 4.;
			pts_fc[i][2] = (pts_bd[hxid[0]][2] + pts_bd[hxid[1]][2] + 2.*pts_fc0[i][2]) / 4.;
		}
		else
		{
			pts_fc[i][0] = pts_fc0[i][0]; pts_fc[i][1] = pts_fc0[i][1]; pts_fc[i][2] = pts_fc0[i][2];
		}
	}
	//edge points
	vector<array<double, 3>> pts_ed0(tmedge.size());//simply average
	vector<array<double, 3>> pts_ed(tmedge.size());//desired
	for (int i = 0; i < tmedge.size(); i++)
	{
		pts_ed0[i][0] = (cp[tmedge[i].pt[0]].coor[0] + cp[tmedge[i].pt[1]].coor[0]) / 2.;
		pts_ed0[i][1] = (cp[tmedge[i].pt[0]].coor[1] + cp[tmedge[i].pt[1]].coor[1]) / 2.;
		pts_ed0[i][2] = (cp[tmedge[i].pt[0]].coor[2] + cp[tmedge[i].pt[1]].coor[2]) / 2.;
		if (tmedge[i].type != 1)//non-boundary, i.e. interior
		{
			double cavg[3] = { 0., 0., 0. };
			for (uint j = 0; j < tmedge[i].hex.size(); j++)
			{
				cavg[0] += pts_bd[tmedge[i].hex[j]][0];
				cavg[1] += pts_bd[tmedge[i].hex[j]][1];
				cavg[2] += pts_bd[tmedge[i].hex[j]][2];
			}
			cavg[0] /= double(tmedge[i].hex.size());
			cavg[1] /= double(tmedge[i].hex.size());
			cavg[2] /= double(tmedge[i].hex.size());
			double aavg[3] = { 0., 0., 0. };
			for (uint j = 0; j < tmedge[i].face.size(); j++)
			{
				aavg[0] += pts_fc0[tmedge[i].face[j]][0];
				aavg[1] += pts_fc0[tmedge[i].face[j]][1];
				aavg[2] += pts_fc0[tmedge[i].face[j]][2];
			}
			aavg[0] /= double(tmedge[i].face.size());
			aavg[1] /= double(tmedge[i].face.size());
			aavg[2] /= double(tmedge[i].face.size());
			pts_ed[i][0] = (cavg[0] + 2.*aavg[0] + double(tmedge[i].face.size() - 3)*pts_ed0[i][0]) / double(tmedge[i].face.size());
			pts_ed[i][1] = (cavg[1] + 2.*aavg[1] + double(tmedge[i].face.size() - 3)*pts_ed0[i][1]) / double(tmedge[i].face.size());
			pts_ed[i][2] = (cavg[2] + 2.*aavg[2] + double(tmedge[i].face.size() - 3)*pts_ed0[i][2]) / double(tmedge[i].face.size());
		}
		else
		{
			if (tmedge[i].sharp == 0)//boundary non-sharp edge
			{
				double aavg[3] = { 0., 0., 0. };
				for (uint j = 0; j < tmedge[i].face.size(); j++)
				{
					if (tmface[tmedge[i].face[j]].type == 1)
					{
						aavg[0] += pts_fc0[tmedge[i].face[j]][0];
						aavg[1] += pts_fc0[tmedge[i].face[j]][1];
						aavg[2] += pts_fc0[tmedge[i].face[j]][2];
					}
				}
				pts_ed[i][0] = (aavg[0] + 2.*pts_ed0[i][0]) / 4.;
				pts_ed[i][1] = (aavg[1] + 2.*pts_ed0[i][1]) / 4.;
				pts_ed[i][2] = (aavg[2] + 2.*pts_ed0[i][2]) / 4.;
			}
			else//sharp edge
			{
				pts_ed[i][0] = pts_ed0[i][0];
				pts_ed[i][1] = pts_ed0[i][1];
				pts_ed[i][2] = pts_ed0[i][2];
			}
		}
	}
	//vertex points
#pragma omp parallel for
	for (int i = 0; i < cp.size(); i++)
	{
		if (cp[i].type != 1)//interior
		{
			double cavg[3] = { 0., 0., 0. };
			for (uint j = 0; j < cp[i].hex.size(); j++)
			{
				cavg[0] += pts_bd[cp[i].hex[j]][0];
				cavg[1] += pts_bd[cp[i].hex[j]][1];
				cavg[2] += pts_bd[cp[i].hex[j]][2];
			}
			cavg[0] /= double(cp[i].hex.size());
			cavg[1] /= double(cp[i].hex.size());
			cavg[2] /= double(cp[i].hex.size());
			double aavg[3] = { 0., 0., 0. };
			for (uint j = 0; j < cp[i].face.size(); j++)
			{
				aavg[0] += pts_fc0[cp[i].face[j]][0];
				aavg[1] += pts_fc0[cp[i].face[j]][1];
				aavg[2] += pts_fc0[cp[i].face[j]][2];
			}
			aavg[0] /= double(cp[i].face.size());
			aavg[1] /= double(cp[i].face.size());
			aavg[2] /= double(cp[i].face.size());
			double mavg[3] = { 0., 0., 0. };
			for (uint j = 0; j < cp[i].edge.size(); j++)
			{
				mavg[0] += pts_ed0[cp[i].edge[j]][0];
				mavg[1] += pts_ed0[cp[i].edge[j]][1];
				mavg[2] += pts_ed0[cp[i].edge[j]][2];
			}
			mavg[0] /= double(cp[i].edge.size());
			mavg[1] /= double(cp[i].edge.size());
			mavg[2] /= double(cp[i].edge.size());
			cp[i].coor[0] = (cavg[0] + 3.*aavg[0] + 3.*mavg[0] + cp[i].coor[0]) / 8.;
			cp[i].coor[1] = (cavg[1] + 3.*aavg[1] + 3.*mavg[1] + cp[i].coor[1]) / 8.;
			cp[i].coor[2] = (cavg[2] + 3.*aavg[2] + 3.*mavg[2] + cp[i].coor[2]) / 8.;
		}
		else
		{
			if (cp[i].sharp == 0)
			{
				double aavg[3] = { 0., 0., 0. };
				int nfc(0);
				for (uint j = 0; j < cp[i].face.size(); j++)
				{
					if (tmface[cp[i].face[j]].type == 1)
					{
						aavg[0] += pts_fc0[cp[i].face[j]][0];
						aavg[1] += pts_fc0[cp[i].face[j]][1];
						aavg[2] += pts_fc0[cp[i].face[j]][2];
						nfc++;
					}
				}
				aavg[0] /= double(nfc); aavg[1] /= double(nfc); aavg[2] /= double(nfc);
				double mavg[3] = { 0., 0., 0. };
				int ned(0);
				for (uint j = 0; j < cp[i].edge.size(); j++)
				{
					if (tmedge[cp[i].edge[j]].type == 1)
					{
						mavg[0] += pts_ed0[cp[i].edge[j]][0];
						mavg[1] += pts_ed0[cp[i].edge[j]][1];
						mavg[2] += pts_ed0[cp[i].edge[j]][2];
						ned++;
					}
				}
				mavg[0] /= double(ned); mavg[1] /= double(ned); mavg[2] /= double(ned);
				cp[i].coor[0] = (aavg[0] + 2.*mavg[0] + double(nfc - 3)*cp[i].coor[0]) / double(nfc);
				cp[i].coor[1] = (aavg[1] + 2.*mavg[1] + double(nfc - 3)*cp[i].coor[1]) / double(nfc);
				cp[i].coor[2] = (aavg[2] + 2.*mavg[2] + double(nfc - 3)*cp[i].coor[2]) / double(nfc);
			}
			else if (cp[i].sharp == 1)//associated with sharp edge
			{
				double mavg[3] = { 0., 0., 0. };
				for (uint j = 0; j < cp[i].edge.size(); j++)
				{
					if (tmedge[cp[i].edge[j]].sharp == 1)
					{
						mavg[0] += pts_ed0[cp[i].edge[j]][0];
						mavg[1] += pts_ed0[cp[i].edge[j]][1];
						mavg[2] += pts_ed0[cp[i].edge[j]][2];
					}
				}
				cp[i].coor[0] = (mavg[0] + 2.*cp[i].coor[0]) / 4.;
				cp[i].coor[1] = (mavg[1] + 2.*cp[i].coor[1]) / 4.;
				cp[i].coor[2] = (mavg[2] + 2.*cp[i].coor[2]) / 4.;
			}
			else//sharp corner
			{
				cp[i].coor[0] = cp[i].coor[0]; cp[i].coor[1] = cp[i].coor[1]; cp[i].coor[2] = cp[i].coor[2];
			}
		}
	}
#pragma omp barrier
	//update cp
	vector<int> pid_vt(cp.size());
	vector<int> pid_ed(tmedge.size());
	vector<int> pid_fc(tmface.size());
	vector<int> pid_bd(tmesh.size());
	int count(0);
	for (int i = 0; i < cp.size(); i++)
	{
		pid_vt[i] = count++;
	}
	for (int i = 0; i < tmedge.size(); i++)
	{
		pid_ed[i] = count++;
		Vertex3D ptmp;
		ptmp.act = 1;
		if (tmedge[i].type == 1) ptmp.type = 1;//default 0
		if (tmedge[i].sharp == 1) ptmp.sharp = 1;//default 0
		ptmp.coor[0] = pts_ed[i][0]; ptmp.coor[1] = pts_ed[i][1]; ptmp.coor[2] = pts_ed[i][2];
		cp.push_back(ptmp);
	}
	for (int i = 0; i < tmface.size(); i++)
	{
		pid_fc[i] = count++;
		Vertex3D ptmp;
		ptmp.act = 1;
		if (tmface[i].type == 1) ptmp.type = 1;//default 0
		ptmp.coor[0] = pts_fc[i][0]; ptmp.coor[1] = pts_fc[i][1]; ptmp.coor[2] = pts_fc[i][2];
		cp.push_back(ptmp);
	}
	for (int i = 0; i < tmesh.size(); i++)
	{
		pid_bd[i] = count++;
		Vertex3D ptmp;
		ptmp.act = 1;
		ptmp.coor[0] = pts_bd[i][0]; ptmp.coor[1] = pts_bd[i][1]; ptmp.coor[2] = pts_bd[i][2];
		cp.push_back(ptmp);
	}
	//update edge
	vector<Edge3D> ednew(2 * tmedge.size() + 4 * tmface.size() + 6 * tmesh.size());
	int edloc(0);
	for (int i = 0; i < tmedge.size(); i++)
	{
		ednew[edloc].pt[0] = tmedge[i].pt[0];
		ednew[edloc].pt[1] = pid_ed[i];
		ednew[edloc].type = tmedge[i].type;
		ednew[edloc].sharp = tmedge[i].sharp;
		edloc++;
		ednew[edloc].pt[0] = pid_ed[i];
		ednew[edloc].pt[1] = tmedge[i].pt[1];
		ednew[edloc].type = tmedge[i].type;
		ednew[edloc].sharp = tmedge[i].sharp;
		edloc++;
	}
	//update face
	vector<Face3D> fcnew(4 * tmface.size() + 12 * tmesh.size());
	int fcloc(0);
	for (int i = 0; i < tmface.size(); i++)
	{
		int fcnct[4][4] = { { tmface[i].cnct[0], pid_ed[tmface[i].edge[0]], pid_fc[i], pid_ed[tmface[i].edge[3]] },
		{ pid_ed[tmface[i].edge[0]], tmface[i].cnct[1], pid_ed[tmface[i].edge[1]], pid_fc[i] },
		{ pid_fc[i], pid_ed[tmface[i].edge[1]], tmface[i].cnct[2], pid_ed[tmface[i].edge[2]] },
		{ pid_ed[tmface[i].edge[3]], pid_fc[i], pid_ed[tmface[i].edge[2]], tmface[i].cnct[3] } };
		int edid[12];
		for (int j = 0; j < 4; j++)
		{
			edid[2 * j] = 2 * tmface[i].edge[j]; edid[2 * j + 1] = 2 * tmface[i].edge[j] + 1;
			if (tmedge[tmface[i].edge[j]].pt[0] != tmface[i].cnct[j])
			{
				edid[2 * j] = 2 * tmface[i].edge[j] + 1; edid[2 * j + 1] = 2 * tmface[i].edge[j];
			}
		}
		//construct 4 new edges
		for (int j = 0; j < 4; j++)
		{
			ednew[edloc].pt[0] = pid_fc[i];
			ednew[edloc].pt[1] = pid_ed[tmface[i].edge[j]];
			if (tmface[i].type == 1) ednew[edloc].type = 1;
			edid[8 + j] = edloc;
			edloc++;
		}
		int fedge[4][4] = { { edid[0], edid[8], edid[11], edid[7] }, { edid[1], edid[2], edid[9], edid[8] },
		{ edid[9], edid[3], edid[4], edid[10] }, { edid[11], edid[10], edid[5], edid[6] } };
		for (int j = 0; j < 4; j++)
		{
			fcnew[fcloc].type = tmface[i].type;
			for (int k = 0; k < 4; k++)
			{
				fcnew[fcloc].cnct[k] = fcnct[j][k];
				fcnew[fcloc].edge[k] = fedge[j][k];
			}
			fcloc++;
		}
	}
	//update hex
	vector<Element3D> hxnew(8 * tmesh.size());
	int hxloc(0);
	for (int i = 0; i < tmesh.size(); i++)
	{
		//construct new edges
		for (int j = 0; j < 6; j++)
		{
			ednew[edloc].pt[0] = pid_bd[i];
			ednew[edloc].pt[1] = pid_fc[tmesh[i].face[j]];
			edloc++;
		}
		//construct new faces
		int fcnct[12][4] = { { pid_ed[tmesh[i].edge[0]], pid_fc[tmesh[i].face[0]], pid_bd[i], pid_fc[tmesh[i].face[1]] },
		{ pid_fc[tmesh[i].face[0]], pid_ed[tmesh[i].edge[2]], pid_fc[tmesh[i].face[3]], pid_bd[i] },
		{ pid_bd[i], pid_fc[tmesh[i].face[3]], pid_ed[tmesh[i].edge[10]], pid_fc[tmesh[i].face[5]] },
		{ pid_fc[tmesh[i].face[1]], pid_bd[i], pid_fc[tmesh[i].face[5]], pid_ed[tmesh[i].edge[8]] },
		{ pid_ed[tmesh[i].edge[3]], pid_fc[tmesh[i].face[0]], pid_bd[i], pid_fc[tmesh[i].face[4]] },
		{ pid_fc[tmesh[i].face[0]], pid_ed[tmesh[i].edge[1]], pid_fc[tmesh[i].face[2]], pid_bd[i] },
		{ pid_bd[i], pid_fc[tmesh[i].face[2]], pid_ed[tmesh[i].edge[9]], pid_fc[tmesh[i].face[5]] },
		{ pid_fc[tmesh[i].face[4]], pid_bd[i], pid_fc[tmesh[i].face[5]], pid_ed[tmesh[i].edge[11]] },
		{ pid_ed[tmesh[i].edge[4]], pid_fc[tmesh[i].face[1]], pid_bd[i], pid_fc[tmesh[i].face[4]] },
		{ pid_fc[tmesh[i].face[1]], pid_ed[tmesh[i].edge[5]], pid_fc[tmesh[i].face[2]], pid_bd[i] },
		{ pid_bd[i], pid_fc[tmesh[i].face[2]], pid_ed[tmesh[i].edge[6]], pid_fc[tmesh[i].face[3]] },
		{ pid_fc[tmesh[i].face[4]], pid_bd[i], pid_fc[tmesh[i].face[3]], pid_ed[tmesh[i].edge[7]] } };
		//collect all edges and faces in this element
		int edids[12 * 2 + 6 * 4 + 6], fcids[6 * 4 + 12];
		for (int j = 0; j < 12; j++)
		{
			edids[2 * j] = 2 * tmesh[i].edge[j];
			edids[2 * j + 1] = 2 * tmesh[i].edge[j] + 1;
			fcids[24 + j] = 4 * tmface.size() + 12 * i + j;
		}
		for (int j = 0; j < 6; j++)
		{
			edids[24 + 4 * j] = 2 * tmedge.size() + 4 * tmesh[i].face[j];
			edids[24 + 4 * j + 1] = 2 * tmedge.size() + 4 * tmesh[i].face[j] + 1;
			edids[24 + 4 * j + 2] = 2 * tmedge.size() + 4 * tmesh[i].face[j] + 2;
			edids[24 + 4 * j + 3] = 2 * tmedge.size() + 4 * tmesh[i].face[j] + 3;
			edids[48 + j] = 2 * tmedge.size() + 4 * tmface.size() + 6 * i + j;
			fcids[4 * j] = 4 * tmesh[i].face[j];
			fcids[4 * j + 1] = 4 * tmesh[i].face[j] + 1;
			fcids[4 * j + 2] = 4 * tmesh[i].face[j] + 2;
			fcids[4 * j + 3] = 4 * tmesh[i].face[j] + 3;
		}
		for (int j = 0; j < 12; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				fcnew[fcloc].cnct[k] = fcnct[j][k];
				for (int k0 = 24; k0 < 54; k0++)
				{
					if ((fcnct[j][k] == ednew[edids[k0]].pt[0] && fcnct[j][(k + 1) % 4] == ednew[edids[k0]].pt[1]) || 
						(fcnct[j][k] == ednew[edids[k0]].pt[1] && fcnct[j][(k + 1) % 4] == ednew[edids[k0]].pt[0]))
					{
						fcnew[fcloc].edge[k] = edids[k0]; break;
					}
				}
			}
			fcloc++;
		}
		//construct new hex
		int ecnct[8][8] = { { tmesh[i].cnct[0], pid_ed[tmesh[i].edge[0]], pid_fc[tmesh[i].face[0]], pid_ed[tmesh[i].edge[3]], pid_ed[tmesh[i].edge[4]], pid_fc[tmesh[i].face[1]], pid_bd[i], pid_fc[tmesh[i].face[4]] },
		{ pid_ed[tmesh[i].edge[0]], tmesh[i].cnct[1], pid_ed[tmesh[i].edge[1]], pid_fc[tmesh[i].face[0]], pid_fc[tmesh[i].face[1]], pid_ed[tmesh[i].edge[5]], pid_fc[tmesh[i].face[2]], pid_bd[i] },
		{ pid_fc[tmesh[i].face[0]], pid_ed[tmesh[i].edge[1]], tmesh[i].cnct[2], pid_ed[tmesh[i].edge[2]], pid_bd[i], pid_fc[tmesh[i].face[2]], pid_ed[tmesh[i].edge[6]], pid_fc[tmesh[i].face[3]] },
		{ pid_ed[tmesh[i].edge[3]], pid_fc[tmesh[i].face[0]], pid_ed[tmesh[i].edge[2]], tmesh[i].cnct[3], pid_fc[tmesh[i].face[4]], pid_bd[i], pid_fc[tmesh[i].face[3]], pid_ed[tmesh[i].edge[7]] },
		{ pid_ed[tmesh[i].edge[4]], pid_fc[tmesh[i].face[1]], pid_bd[i], pid_fc[tmesh[i].face[4]], tmesh[i].cnct[4], pid_ed[tmesh[i].edge[8]], pid_fc[tmesh[i].face[5]], pid_ed[tmesh[i].edge[11]] },
		{ pid_fc[tmesh[i].face[1]], pid_ed[tmesh[i].edge[5]], pid_fc[tmesh[i].face[2]], pid_bd[i], pid_ed[tmesh[i].edge[8]], tmesh[i].cnct[5], pid_ed[tmesh[i].edge[9]], pid_fc[tmesh[i].face[5]] },
		{ pid_bd[i], pid_fc[tmesh[i].face[2]], pid_ed[tmesh[i].edge[6]], pid_fc[tmesh[i].face[3]], pid_fc[tmesh[i].face[5]], pid_ed[tmesh[i].edge[9]], tmesh[i].cnct[6], pid_ed[tmesh[i].edge[10]] },
		{ pid_fc[tmesh[i].face[4]], pid_bd[i], pid_fc[tmesh[i].face[3]], pid_ed[tmesh[i].edge[7]], pid_ed[tmesh[i].edge[11]], pid_fc[tmesh[i].face[5]], pid_ed[tmesh[i].edge[10]], tmesh[i].cnct[7] } };
		for (int j = 0; j < 8; j++)
		{
			for (int k = 0; k < 8; k++)
			{
				hxnew[hxloc].cnct[k] = ecnct[j][k];
			}
			//find face connectivity
			for (int k = 0; k < 6; k++)
			{
				array<int, 4> tmp1 = { ecnct[j][solid_fc[k][0]], ecnct[j][solid_fc[k][1]], ecnct[j][solid_fc[k][2]], ecnct[j][solid_fc[k][3]] };
				sort(tmp1.begin(), tmp1.end());
				for (int k0 = 0; k0 < 36; k0++)
				{
					array<int, 4> tmp2 = { fcnew[fcids[k0]].cnct[0], fcnew[fcids[k0]].cnct[1], fcnew[fcids[k0]].cnct[2], fcnew[fcids[k0]].cnct[3] };
					sort(tmp2.begin(), tmp2.end());
					if (tmp1 == tmp2)
					{
						hxnew[hxloc].face[k] = fcids[k0]; break;
					}
				}
			}
			//find edge connectivity
			for (int k = 0; k < 12; k++)
			{
				array<int, 2> tmp1 = { ecnct[j][solid_ed[k][0]], ecnct[j][solid_ed[k][1]] };
				for (int k0 = 0; k0 < 54; k0++)
				{
					if ((tmp1[0] == ednew[edids[k0]].pt[0] && tmp1[1] == ednew[edids[k0]].pt[1]) ||
						(tmp1[0] == ednew[edids[k0]].pt[1] && tmp1[1] == ednew[edids[k0]].pt[0]))
					{
						hxnew[hxloc].edge[k] = edids[k0]; break;
					}
				}
			}
			hxloc++;
		}
	}
	//update connectivity
	tmedge.resize(ednew.size());
	tmface.resize(fcnew.size());
	tmesh.resize(hxnew.size());
	for (int i = 0; i < cp.size(); i++)
	{
		cp[i].edge.clear();
		cp[i].face.clear();
		cp[i].hex.clear();
	}
	for (int i = 0; i < ednew.size(); i++)
	{
		tmedge[i].act = 1;
		tmedge[i].type = ednew[i].type;
		tmedge[i].sharp = ednew[i].sharp;
		tmedge[i].pt[0] = ednew[i].pt[0]; tmedge[i].pt[1] = ednew[i].pt[1];
		tmedge[i].face.clear();
		tmedge[i].hex.clear();
	}
	for (int i = 0; i < fcnew.size(); i++)
	{
		tmface[i].act = 1;
		tmface[i].type = fcnew[i].type;
		for (int j = 0; j < 4; j++)
		{
			tmface[i].cnct[j] = fcnew[i].cnct[j];
			tmface[i].edge[j] = fcnew[i].edge[j];
		}
		tmface[i].hex.clear();
	}
	for (int i = 0; i < hxnew.size(); i++)
	{
		tmesh[i].act = 1;
		tmesh[i].type = hxnew[i].type;
		for (int j = 0; j < 8; j++)
		{
			tmesh[i].cnct[j] = hxnew[i].cnct[j];
		}
		for (int j = 0; j < 12; j++)
		{
			tmesh[i].edge[j] = hxnew[i].edge[j];
		}
		for (int j = 0; j < 6; j++)
		{
			tmesh[i].face[j] = hxnew[i].face[j];
		}
	}
	//vertex-to-hex, edge-to-hex, face-to-hex
	for (int i = 0; i<tmesh.size(); i++)
	{
		for (int j = 0; j<8; j++)
		{
			cp[tmesh[i].cnct[j]].hex.push_back(i);
		}
		for (int j = 0; j<12; j++)
		{
			tmedge[tmesh[i].edge[j]].hex.push_back(i);
		}
		for (int j = 0; j<6; j++)
		{
			tmface[tmesh[i].face[j]].hex.push_back(i);
		}
	}
	//vertex-to-face, edge-to-face
	for (int i = 0; i<tmface.size(); i++)
	{
		for (int j = 0; j<4; j++)
		{
			cp[tmface[i].cnct[j]].face.push_back(i);
			tmedge[tmface[i].edge[j]].face.push_back(i);
		}
	}
	//vertex-to-edge
	for (int i = 0; i<tmedge.size(); i++)
	{
		for (int j = 0; j<2; j++)
		{
			cp[tmedge[i].pt[j]].edge.push_back(i);
		}
	}
	//find extraordinary edges and vertices
	for (int i = 0; i<tmedge.size(); i++)
	{
		if (tmedge[i].type != 1 && tmedge[i].hex.size() != 4)
		{
			tmedge[i].type = 2;
			if (cp[tmedge[i].pt[0]].type != 1)
				cp[tmedge[i].pt[0]].type = 3;
			if (cp[tmedge[i].pt[1]].type != 1)
				cp[tmedge[i].pt[1]].type = 3;
		}
	}
	//
	//find boundary and irregular elements
	for (int i = 0; i<tmesh.size(); i++)
	{
		for (int j = 0; j < 8; j++)
		{
			if (cp[tmesh[i].cnct[j]].type == 1)
			{
				tmesh[i].type = 1; break;
			}
		}
		if (tmesh[i].type != 1)
		{
			for (int j = 0; j<12; j++)
			{
				if (tmedge[tmesh[i].edge[j]].type == 2)
				{
					tmesh[i].type = 2;
					break;
				}
			}
		}
	}

	//boundry extraordinary points
	//for (int i = 0; i<cp.size(); i++)
	//{
	//	if (cp[i].type == 1)
	//	{
	//		int count(0);
	//		for (j = 0; j < cp[i].edge.size(); j++)
	//		{
	//			if (tmedge[cp[i].edge[j]].type == 2) count++;
	//		}
	//		if (count == 1) cp[i].bcxp = 1;
	//		else if (count>1) cp[i].bcxp = 2;
	//	}
	//}
}

void TruncatedTspline_3D::Global_Subdivide_Simple()
{
	//body points
	vector<array<double, 3>> pts_bd(tmesh.size());
	for (int i = 0; i < tmesh.size(); i++)
	{
		pts_bd[i][0] = 0.; pts_bd[i][1] = 0.; pts_bd[i][2] = 0.;
		for (int j = 0; j < 8; j++)
		{
			pts_bd[i][0] += cp[tmesh[i].cnct[j]].coor[0];
			pts_bd[i][1] += cp[tmesh[i].cnct[j]].coor[1];
			pts_bd[i][2] += cp[tmesh[i].cnct[j]].coor[2];
		}
		pts_bd[i][0] /= 8.; pts_bd[i][1] /= 8.; pts_bd[i][2] /= 8.;
	}
	//face points
	vector<array<double, 3>> pts_fc0(tmface.size());//simply average
	vector<array<double, 3>> pts_fc(tmface.size());//desired
	for (int i = 0; i < tmface.size(); i++)
	{
		pts_fc0[i][0] = 0.; pts_fc0[i][1] = 0.; pts_fc0[i][2] = 0.;
		for (int j = 0; j < 4; j++)
		{
			pts_fc0[i][0] += cp[tmface[i].cnct[j]].coor[0];
			pts_fc0[i][1] += cp[tmface[i].cnct[j]].coor[1];
			pts_fc0[i][2] += cp[tmface[i].cnct[j]].coor[2];
		}
		pts_fc0[i][0] /= 4.; pts_fc0[i][1] /= 4.; pts_fc0[i][2] /= 4.;
		pts_fc[i][0] = pts_fc0[i][0]; pts_fc[i][1] = pts_fc0[i][1]; pts_fc[i][2] = pts_fc0[i][2];
		//if (tmface[i].hex.size() == 2)//interior face points
		//{
		//	int hxid[2] = { tmface[i].hex[0], tmface[i].hex[1] };
		//	pts_fc[i][0] = (pts_bd[hxid[0]][0] + pts_bd[hxid[1]][0] + 2.*pts_fc0[i][0]) / 4.;
		//	pts_fc[i][1] = (pts_bd[hxid[0]][1] + pts_bd[hxid[1]][1] + 2.*pts_fc0[i][1]) / 4.;
		//	pts_fc[i][2] = (pts_bd[hxid[0]][2] + pts_bd[hxid[1]][2] + 2.*pts_fc0[i][2]) / 4.;
		//	
		//}
		//else
		//{
		//	pts_fc[i][0] = pts_fc0[i][0]; pts_fc[i][1] = pts_fc0[i][1]; pts_fc[i][2] = pts_fc0[i][2];
		//}
	}
	//edge points
	vector<array<double, 3>> pts_ed0(tmedge.size());//simply average
	vector<array<double, 3>> pts_ed(tmedge.size());//desired
	for (int i = 0; i < tmedge.size(); i++)
	{
		pts_ed0[i][0] = (cp[tmedge[i].pt[0]].coor[0] + cp[tmedge[i].pt[1]].coor[0]) / 2.;
		pts_ed0[i][1] = (cp[tmedge[i].pt[0]].coor[1] + cp[tmedge[i].pt[1]].coor[1]) / 2.;
		pts_ed0[i][2] = (cp[tmedge[i].pt[0]].coor[2] + cp[tmedge[i].pt[1]].coor[2]) / 2.;
		pts_ed[i][0] = pts_ed0[i][0];
		pts_ed[i][1] = pts_ed0[i][1];
		pts_ed[i][2] = pts_ed0[i][2];
		//if (tmedge[i].type != 1)//non-boundary, i.e. interior
		//{
		//	double cavg[3] = { 0., 0., 0. };
		//	for (uint j = 0; j < tmedge[i].hex.size(); j++)
		//	{
		//		cavg[0] += pts_bd[tmedge[i].hex[j]][0];
		//		cavg[1] += pts_bd[tmedge[i].hex[j]][1];
		//		cavg[2] += pts_bd[tmedge[i].hex[j]][2];
		//	}
		//	cavg[0] /= double(tmedge[i].hex.size());
		//	cavg[1] /= double(tmedge[i].hex.size());
		//	cavg[2] /= double(tmedge[i].hex.size());
		//	double aavg[3] = { 0., 0., 0. };
		//	for (uint j = 0; j < tmedge[i].face.size(); j++)
		//	{
		//		aavg[0] += pts_fc0[tmedge[i].face[j]][0];
		//		aavg[1] += pts_fc0[tmedge[i].face[j]][1];
		//		aavg[2] += pts_fc0[tmedge[i].face[j]][2];
		//	}
		//	aavg[0] /= double(tmedge[i].face.size());
		//	aavg[1] /= double(tmedge[i].face.size());
		//	aavg[2] /= double(tmedge[i].face.size());
		//	pts_ed[i][0] = (cavg[0] + 2.*aavg[0] + double(tmedge[i].face.size() - 3)*pts_ed0[i][0]) / double(tmedge[i].face.size());
		//	pts_ed[i][1] = (cavg[1] + 2.*aavg[1] + double(tmedge[i].face.size() - 3)*pts_ed0[i][1]) / double(tmedge[i].face.size());
		//	pts_ed[i][2] = (cavg[2] + 2.*aavg[2] + double(tmedge[i].face.size() - 3)*pts_ed0[i][2]) / double(tmedge[i].face.size());
		//}
		//else
		//{
		//	if (tmedge[i].sharp == 0)//boundary non-sharp edge
		//	{
		//		double aavg[3] = { 0., 0., 0. };
		//		for (uint j = 0; j < tmedge[i].face.size(); j++)
		//		{
		//			if (tmface[tmedge[i].face[j]].type == 1)
		//			{
		//				aavg[0] += pts_fc0[tmedge[i].face[j]][0];
		//				aavg[1] += pts_fc0[tmedge[i].face[j]][1];
		//				aavg[2] += pts_fc0[tmedge[i].face[j]][2];
		//			}
		//		}
		//		pts_ed[i][0] = (aavg[0] + 2.*pts_ed0[i][0]) / 4.;
		//		pts_ed[i][1] = (aavg[1] + 2.*pts_ed0[i][1]) / 4.;
		//		pts_ed[i][2] = (aavg[2] + 2.*pts_ed0[i][2]) / 4.;
		//	}
		//	else//sharp edge
		//	{
		//		pts_ed[i][0] = pts_ed0[i][0];
		//		pts_ed[i][1] = pts_ed0[i][1];
		//		pts_ed[i][2] = pts_ed0[i][2];
		//	}
		//}
	}
	//vertex points
	//for (int i = 0; i < cp.size(); i++)
	//{
	//	if (cp[i].type != 1)//interior
	//	{
	//		double cavg[3] = { 0., 0., 0. };
	//		for (uint j = 0; j < cp[i].hex.size(); j++)
	//		{
	//			cavg[0] += pts_bd[cp[i].hex[j]][0];
	//			cavg[1] += pts_bd[cp[i].hex[j]][1];
	//			cavg[2] += pts_bd[cp[i].hex[j]][2];
	//		}
	//		cavg[0] /= double(cp[i].hex.size());
	//		cavg[1] /= double(cp[i].hex.size());
	//		cavg[2] /= double(cp[i].hex.size());
	//		double aavg[3] = { 0., 0., 0. };
	//		for (uint j = 0; j < cp[i].face.size(); j++)
	//		{
	//			aavg[0] += pts_fc0[cp[i].face[j]][0];
	//			aavg[1] += pts_fc0[cp[i].face[j]][1];
	//			aavg[2] += pts_fc0[cp[i].face[j]][2];
	//		}
	//		aavg[0] /= double(cp[i].face.size());
	//		aavg[1] /= double(cp[i].face.size());
	//		aavg[2] /= double(cp[i].face.size());
	//		double mavg[3] = { 0., 0., 0. };
	//		for (uint j = 0; j < cp[i].edge.size(); j++)
	//		{
	//			mavg[0] += pts_ed0[cp[i].edge[j]][0];
	//			mavg[1] += pts_ed0[cp[i].edge[j]][1];
	//			mavg[2] += pts_ed0[cp[i].edge[j]][2];
	//		}
	//		mavg[0] /= double(cp[i].edge.size());
	//		mavg[1] /= double(cp[i].edge.size());
	//		mavg[2] /= double(cp[i].edge.size());
	//		cp[i].coor[0] = (cavg[0] + 3.*aavg[0] + 3.*mavg[0] + cp[i].coor[0]) / 8.;
	//		cp[i].coor[1] = (cavg[1] + 3.*aavg[1] + 3.*mavg[1] + cp[i].coor[1]) / 8.;
	//		cp[i].coor[2] = (cavg[2] + 3.*aavg[2] + 3.*mavg[2] + cp[i].coor[2]) / 8.;
	//	}
	//	else
	//	{
	//		if (cp[i].sharp == 0)
	//		{
	//			double aavg[3] = { 0., 0., 0. };
	//			int nfc(0);
	//			for (uint j = 0; j < cp[i].face.size(); j++)
	//			{
	//				if (tmface[cp[i].face[j]].type == 1)
	//				{
	//					aavg[0] += pts_fc0[cp[i].face[j]][0];
	//					aavg[1] += pts_fc0[cp[i].face[j]][1];
	//					aavg[2] += pts_fc0[cp[i].face[j]][2];
	//					nfc++;
	//				}
	//			}
	//			aavg[0] /= double(nfc); aavg[1] /= double(nfc); aavg[2] /= double(nfc);
	//			double mavg[3] = { 0., 0., 0. };
	//			int ned(0);
	//			for (uint j = 0; j < cp[i].edge.size(); j++)
	//			{
	//				if (tmedge[cp[i].edge[j]].type == 1)
	//				{
	//					mavg[0] += pts_ed0[cp[i].edge[j]][0];
	//					mavg[1] += pts_ed0[cp[i].edge[j]][1];
	//					mavg[2] += pts_ed0[cp[i].edge[j]][2];
	//					ned++;
	//				}
	//			}
	//			mavg[0] /= double(ned); mavg[1] /= double(ned); mavg[2] /= double(ned);
	//			cp[i].coor[0] = (aavg[0] + 2.*mavg[0] + double(nfc - 3)*cp[i].coor[0]) / double(nfc);
	//			cp[i].coor[1] = (aavg[1] + 2.*mavg[1] + double(nfc - 3)*cp[i].coor[1]) / double(nfc);
	//			cp[i].coor[2] = (aavg[2] + 2.*mavg[2] + double(nfc - 3)*cp[i].coor[2]) / double(nfc);
	//		}
	//		else if (cp[i].sharp == 1)//associated with sharp edge
	//		{
	//			double mavg[3] = { 0., 0., 0. };
	//			for (uint j = 0; j < cp[i].edge.size(); j++)
	//			{
	//				if (tmedge[cp[i].edge[j]].sharp == 1)
	//				{
	//					mavg[0] += pts_ed0[cp[i].edge[j]][0];
	//					mavg[1] += pts_ed0[cp[i].edge[j]][1];
	//					mavg[2] += pts_ed0[cp[i].edge[j]][2];
	//				}
	//			}
	//			cp[i].coor[0] = (mavg[0] + 2.*cp[i].coor[0]) / 4.;
	//			cp[i].coor[1] = (mavg[1] + 2.*cp[i].coor[1]) / 4.;
	//			cp[i].coor[2] = (mavg[2] + 2.*cp[i].coor[2]) / 4.;
	//		}
	//		else//sharp corner
	//		{
	//			cp[i].coor[0] = cp[i].coor[0]; cp[i].coor[1] = cp[i].coor[1]; cp[i].coor[2] = cp[i].coor[2];
	//		}
	//	}
	//}
	//update cp
	vector<int> pid_vt(cp.size());
	vector<int> pid_ed(tmedge.size());
	vector<int> pid_fc(tmface.size());
	vector<int> pid_bd(tmesh.size());
	int count(0);
	for (int i = 0; i < cp.size(); i++)
	{
		pid_vt[i] = count++;
	}
	for (int i = 0; i < tmedge.size(); i++)
	{
		pid_ed[i] = count++;
		Vertex3D ptmp;
		ptmp.act = 1;
		if (tmedge[i].type == 1) ptmp.type = 1;//default 0
		if (tmedge[i].sharp == 1) ptmp.sharp = 1;//default 0
		ptmp.coor[0] = pts_ed[i][0]; ptmp.coor[1] = pts_ed[i][1]; ptmp.coor[2] = pts_ed[i][2];
		cp.push_back(ptmp);
	}
	for (int i = 0; i < tmface.size(); i++)
	{
		pid_fc[i] = count++;
		Vertex3D ptmp;
		ptmp.act = 1;
		if (tmface[i].type == 1) ptmp.type = 1;//default 0
		ptmp.coor[0] = pts_fc[i][0]; ptmp.coor[1] = pts_fc[i][1]; ptmp.coor[2] = pts_fc[i][2];
		cp.push_back(ptmp);
	}
	for (int i = 0; i < tmesh.size(); i++)
	{
		pid_bd[i] = count++;
		Vertex3D ptmp;
		ptmp.act = 1;
		ptmp.coor[0] = pts_bd[i][0]; ptmp.coor[1] = pts_bd[i][1]; ptmp.coor[2] = pts_bd[i][2];
		cp.push_back(ptmp);
	}
	//update edge
	vector<Edge3D> ednew(2 * tmedge.size() + 4 * tmface.size() + 6 * tmesh.size());
	int edloc(0);
	for (int i = 0; i < tmedge.size(); i++)
	{
		ednew[edloc].pt[0] = tmedge[i].pt[0];
		ednew[edloc].pt[1] = pid_ed[i];
		ednew[edloc].type = tmedge[i].type;
		ednew[edloc].sharp = tmedge[i].sharp;
		edloc++;
		ednew[edloc].pt[0] = pid_ed[i];
		ednew[edloc].pt[1] = tmedge[i].pt[1];
		ednew[edloc].type = tmedge[i].type;
		ednew[edloc].sharp = tmedge[i].sharp;
		edloc++;
	}
	//update face
	vector<Face3D> fcnew(4 * tmface.size() + 12 * tmesh.size());
	int fcloc(0);
	for (int i = 0; i < tmface.size(); i++)
	{
		int fcnct[4][4] = { { tmface[i].cnct[0], pid_ed[tmface[i].edge[0]], pid_fc[i], pid_ed[tmface[i].edge[3]] },
		{ pid_ed[tmface[i].edge[0]], tmface[i].cnct[1], pid_ed[tmface[i].edge[1]], pid_fc[i] },
		{ pid_fc[i], pid_ed[tmface[i].edge[1]], tmface[i].cnct[2], pid_ed[tmface[i].edge[2]] },
		{ pid_ed[tmface[i].edge[3]], pid_fc[i], pid_ed[tmface[i].edge[2]], tmface[i].cnct[3] } };
		int edid[12];
		for (int j = 0; j < 4; j++)
		{
			edid[2 * j] = 2 * tmface[i].edge[j]; edid[2 * j + 1] = 2 * tmface[i].edge[j] + 1;
			if (tmedge[tmface[i].edge[j]].pt[0] != tmface[i].cnct[j])
			{
				edid[2 * j] = 2 * tmface[i].edge[j] + 1; edid[2 * j + 1] = 2 * tmface[i].edge[j];
			}
		}
		//construct 4 new edges
		for (int j = 0; j < 4; j++)
		{
			ednew[edloc].pt[0] = pid_fc[i];
			ednew[edloc].pt[1] = pid_ed[tmface[i].edge[j]];
			if (tmface[i].type == 1) ednew[edloc].type = 1;
			edid[8 + j] = edloc;
			edloc++;
		}
		int fedge[4][4] = { { edid[0], edid[8], edid[11], edid[7] }, { edid[1], edid[2], edid[9], edid[8] },
		{ edid[9], edid[3], edid[4], edid[10] }, { edid[11], edid[10], edid[5], edid[6] } };
		for (int j = 0; j < 4; j++)
		{
			fcnew[fcloc].type = tmface[i].type;
			for (int k = 0; k < 4; k++)
			{
				fcnew[fcloc].cnct[k] = fcnct[j][k];
				fcnew[fcloc].edge[k] = fedge[j][k];
			}
			fcloc++;
		}
	}
	//update hex
	vector<Element3D> hxnew(8 * tmesh.size());
	int hxloc(0);
	for (int i = 0; i < tmesh.size(); i++)
	{
		//construct new edges
		for (int j = 0; j < 6; j++)
		{
			ednew[edloc].pt[0] = pid_bd[i];
			ednew[edloc].pt[1] = pid_fc[tmesh[i].face[j]];
			edloc++;
		}
		//construct new faces
		int fcnct[12][4] = { { pid_ed[tmesh[i].edge[0]], pid_fc[tmesh[i].face[0]], pid_bd[i], pid_fc[tmesh[i].face[1]] },
		{ pid_fc[tmesh[i].face[0]], pid_ed[tmesh[i].edge[2]], pid_fc[tmesh[i].face[3]], pid_bd[i] },
		{ pid_bd[i], pid_fc[tmesh[i].face[3]], pid_ed[tmesh[i].edge[10]], pid_fc[tmesh[i].face[5]] },
		{ pid_fc[tmesh[i].face[1]], pid_bd[i], pid_fc[tmesh[i].face[5]], pid_ed[tmesh[i].edge[8]] },
		{ pid_ed[tmesh[i].edge[3]], pid_fc[tmesh[i].face[0]], pid_bd[i], pid_fc[tmesh[i].face[4]] },
		{ pid_fc[tmesh[i].face[0]], pid_ed[tmesh[i].edge[1]], pid_fc[tmesh[i].face[2]], pid_bd[i] },
		{ pid_bd[i], pid_fc[tmesh[i].face[2]], pid_ed[tmesh[i].edge[9]], pid_fc[tmesh[i].face[5]] },
		{ pid_fc[tmesh[i].face[4]], pid_bd[i], pid_fc[tmesh[i].face[5]], pid_ed[tmesh[i].edge[11]] },
		{ pid_ed[tmesh[i].edge[4]], pid_fc[tmesh[i].face[1]], pid_bd[i], pid_fc[tmesh[i].face[4]] },
		{ pid_fc[tmesh[i].face[1]], pid_ed[tmesh[i].edge[5]], pid_fc[tmesh[i].face[2]], pid_bd[i] },
		{ pid_bd[i], pid_fc[tmesh[i].face[2]], pid_ed[tmesh[i].edge[6]], pid_fc[tmesh[i].face[3]] },
		{ pid_fc[tmesh[i].face[4]], pid_bd[i], pid_fc[tmesh[i].face[3]], pid_ed[tmesh[i].edge[7]] } };
		//collect all edges and faces in this element
		int edids[12 * 2 + 6 * 4 + 6], fcids[6 * 4 + 12];
		for (int j = 0; j < 12; j++)
		{
			edids[2 * j] = 2 * tmesh[i].edge[j];
			edids[2 * j + 1] = 2 * tmesh[i].edge[j] + 1;
			fcids[24 + j] = 4 * tmface.size() + 12 * i + j;
		}
		for (int j = 0; j < 6; j++)
		{
			edids[24 + 4 * j] = 2 * tmedge.size() + 4 * tmesh[i].face[j];
			edids[24 + 4 * j + 1] = 2 * tmedge.size() + 4 * tmesh[i].face[j] + 1;
			edids[24 + 4 * j + 2] = 2 * tmedge.size() + 4 * tmesh[i].face[j] + 2;
			edids[24 + 4 * j + 3] = 2 * tmedge.size() + 4 * tmesh[i].face[j] + 3;
			edids[48 + j] = 2 * tmedge.size() + 4 * tmface.size() + 6 * i + j;
			fcids[4 * j] = 4 * tmesh[i].face[j];
			fcids[4 * j + 1] = 4 * tmesh[i].face[j] + 1;
			fcids[4 * j + 2] = 4 * tmesh[i].face[j] + 2;
			fcids[4 * j + 3] = 4 * tmesh[i].face[j] + 3;
		}
		for (int j = 0; j < 12; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				fcnew[fcloc].cnct[k] = fcnct[j][k];
				for (int k0 = 24; k0 < 54; k0++)
				{
					if ((fcnct[j][k] == ednew[edids[k0]].pt[0] && fcnct[j][(k + 1) % 4] == ednew[edids[k0]].pt[1]) ||
						(fcnct[j][k] == ednew[edids[k0]].pt[1] && fcnct[j][(k + 1) % 4] == ednew[edids[k0]].pt[0]))
					{
						fcnew[fcloc].edge[k] = edids[k0]; break;
					}
				}
			}
			fcloc++;
		}
		//construct new hex
		int ecnct[8][8] = { { tmesh[i].cnct[0], pid_ed[tmesh[i].edge[0]], pid_fc[tmesh[i].face[0]], pid_ed[tmesh[i].edge[3]], pid_ed[tmesh[i].edge[4]], pid_fc[tmesh[i].face[1]], pid_bd[i], pid_fc[tmesh[i].face[4]] },
		{ pid_ed[tmesh[i].edge[0]], tmesh[i].cnct[1], pid_ed[tmesh[i].edge[1]], pid_fc[tmesh[i].face[0]], pid_fc[tmesh[i].face[1]], pid_ed[tmesh[i].edge[5]], pid_fc[tmesh[i].face[2]], pid_bd[i] },
		{ pid_fc[tmesh[i].face[0]], pid_ed[tmesh[i].edge[1]], tmesh[i].cnct[2], pid_ed[tmesh[i].edge[2]], pid_bd[i], pid_fc[tmesh[i].face[2]], pid_ed[tmesh[i].edge[6]], pid_fc[tmesh[i].face[3]] },
		{ pid_ed[tmesh[i].edge[3]], pid_fc[tmesh[i].face[0]], pid_ed[tmesh[i].edge[2]], tmesh[i].cnct[3], pid_fc[tmesh[i].face[4]], pid_bd[i], pid_fc[tmesh[i].face[3]], pid_ed[tmesh[i].edge[7]] },
		{ pid_ed[tmesh[i].edge[4]], pid_fc[tmesh[i].face[1]], pid_bd[i], pid_fc[tmesh[i].face[4]], tmesh[i].cnct[4], pid_ed[tmesh[i].edge[8]], pid_fc[tmesh[i].face[5]], pid_ed[tmesh[i].edge[11]] },
		{ pid_fc[tmesh[i].face[1]], pid_ed[tmesh[i].edge[5]], pid_fc[tmesh[i].face[2]], pid_bd[i], pid_ed[tmesh[i].edge[8]], tmesh[i].cnct[5], pid_ed[tmesh[i].edge[9]], pid_fc[tmesh[i].face[5]] },
		{ pid_bd[i], pid_fc[tmesh[i].face[2]], pid_ed[tmesh[i].edge[6]], pid_fc[tmesh[i].face[3]], pid_fc[tmesh[i].face[5]], pid_ed[tmesh[i].edge[9]], tmesh[i].cnct[6], pid_ed[tmesh[i].edge[10]] },
		{ pid_fc[tmesh[i].face[4]], pid_bd[i], pid_fc[tmesh[i].face[3]], pid_ed[tmesh[i].edge[7]], pid_ed[tmesh[i].edge[11]], pid_fc[tmesh[i].face[5]], pid_ed[tmesh[i].edge[10]], tmesh[i].cnct[7] } };
		for (int j = 0; j < 8; j++)
		{
			for (int k = 0; k < 8; k++)
			{
				hxnew[hxloc].cnct[k] = ecnct[j][k];
			}
			//find face connectivity
			for (int k = 0; k < 6; k++)
			{
				array<int, 4> tmp1 = { ecnct[j][solid_fc[k][0]], ecnct[j][solid_fc[k][1]], ecnct[j][solid_fc[k][2]], ecnct[j][solid_fc[k][3]] };
				sort(tmp1.begin(), tmp1.end());
				for (int k0 = 0; k0 < 36; k0++)
				{
					array<int, 4> tmp2 = { fcnew[fcids[k0]].cnct[0], fcnew[fcids[k0]].cnct[1], fcnew[fcids[k0]].cnct[2], fcnew[fcids[k0]].cnct[3] };
					sort(tmp2.begin(), tmp2.end());
					if (tmp1 == tmp2)
					{
						hxnew[hxloc].face[k] = fcids[k0]; break;
					}
				}
			}
			//find edge connectivity
			for (int k = 0; k < 12; k++)
			{
				array<int, 2> tmp1 = { ecnct[j][solid_ed[k][0]], ecnct[j][solid_ed[k][1]] };
				for (int k0 = 0; k0 < 54; k0++)
				{
					if ((tmp1[0] == ednew[edids[k0]].pt[0] && tmp1[1] == ednew[edids[k0]].pt[1]) ||
						(tmp1[0] == ednew[edids[k0]].pt[1] && tmp1[1] == ednew[edids[k0]].pt[0]))
					{
						hxnew[hxloc].edge[k] = edids[k0]; break;
					}
				}
			}
			hxloc++;
		}
	}
	//update connectivity
	tmedge.resize(ednew.size());
	tmface.resize(fcnew.size());
	tmesh.resize(hxnew.size());
	for (int i = 0; i < cp.size(); i++)
	{
		cp[i].edge.clear();
		cp[i].face.clear();
		cp[i].hex.clear();
	}
	for (int i = 0; i < ednew.size(); i++)
	{
		tmedge[i].act = 1;
		tmedge[i].type = ednew[i].type;
		tmedge[i].sharp = ednew[i].sharp;
		tmedge[i].pt[0] = ednew[i].pt[0]; tmedge[i].pt[1] = ednew[i].pt[1];
		tmedge[i].face.clear();
		tmedge[i].hex.clear();
	}
	for (int i = 0; i < fcnew.size(); i++)
	{
		tmface[i].act = 1;
		tmface[i].type = fcnew[i].type;
		for (int j = 0; j < 4; j++)
		{
			tmface[i].cnct[j] = fcnew[i].cnct[j];
			tmface[i].edge[j] = fcnew[i].edge[j];
		}
		tmface[i].hex.clear();
	}
	for (int i = 0; i < hxnew.size(); i++)
	{
		tmesh[i].act = 1;
		tmesh[i].type = hxnew[i].type;
		for (int j = 0; j < 8; j++)
		{
			tmesh[i].cnct[j] = hxnew[i].cnct[j];
		}
		for (int j = 0; j < 12; j++)
		{
			tmesh[i].edge[j] = hxnew[i].edge[j];
		}
		for (int j = 0; j < 6; j++)
		{
			tmesh[i].face[j] = hxnew[i].face[j];
		}
	}
	//vertex-to-hex, edge-to-hex, face-to-hex
	for (int i = 0; i<tmesh.size(); i++)
	{
		for (int j = 0; j<8; j++)
		{
			cp[tmesh[i].cnct[j]].hex.push_back(i);
		}
		for (int j = 0; j<12; j++)
		{
			tmedge[tmesh[i].edge[j]].hex.push_back(i);
		}
		for (int j = 0; j<6; j++)
		{
			tmface[tmesh[i].face[j]].hex.push_back(i);
		}
	}
	//vertex-to-face, edge-to-face
	for (int i = 0; i<tmface.size(); i++)
	{
		for (int j = 0; j<4; j++)
		{
			cp[tmface[i].cnct[j]].face.push_back(i);
			tmedge[tmface[i].edge[j]].face.push_back(i);
		}
	}
	//vertex-to-edge
	for (int i = 0; i<tmedge.size(); i++)
	{
		for (int j = 0; j<2; j++)
		{
			cp[tmedge[i].pt[j]].edge.push_back(i);
		}
	}
	//find extraordinary edges and vertices
	for (int i = 0; i<tmedge.size(); i++)
	{
		if (tmedge[i].type != 1 && tmedge[i].hex.size() != 4)
		{
			tmedge[i].type = 2;
			if (cp[tmedge[i].pt[0]].type != 1)
				cp[tmedge[i].pt[0]].type = 3;
			if (cp[tmedge[i].pt[1]].type != 1)
				cp[tmedge[i].pt[1]].type = 3;
		}
	}
	//
	//find boundary and irregular elements
	for (int i = 0; i<tmesh.size(); i++)
	{
		for (int j = 0; j < 8; j++)
		{
			if (cp[tmesh[i].cnct[j]].type == 1)
			{
				tmesh[i].type = 1; break;
			}
		}
		if (tmesh[i].type != 1)
		{
			for (int j = 0; j<12; j++)
			{
				if (tmedge[tmesh[i].edge[j]].type == 2)
				{
					tmesh[i].type = 2;
					break;
				}
			}
		}
	}
}

void TruncatedTspline_3D::BuildBasisFunction()
{
	hmesh.push_back(tmesh);
	hcp.push_back(cp);
	hface.push_back(tmface);
	hedge.push_back(tmedge);

	vector<Vertex3D>().swap(cp);
	vector<Edge3D>().swap(tmedge);
	vector<Face3D>().swap(tmface);
	vector<Element3D>().swap(tmesh);

	//ReportXP();

	AllBezierLev(0);
	//for (uint i = 0; i < hmesh[0].size(); i++)
	//{
	//	if (hmesh[0][i].type == 1)
	//	{
	//		ConstructBezierBasis_Boundary(0, i);
	//	}
	//}
}

void TruncatedTspline_3D::OutputCM(string fn)
{
	string fname(fn + "_CM.vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << cp.size() << " float\n";
		for (uint i = 0; i<cp.size(); i++)
		{
			fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << tmesh.size() << " " << 9 * tmesh.size() << '\n';
		for (uint i = 0; i<tmesh.size(); i++)
		{
			fout << "8 ";
			for (int j = 0; j<8; j++)
			{
				fout << tmesh[i].cnct[j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << tmesh.size() << '\n';
		for (uint i = 0; i<tmesh.size(); i++)
		{
			fout << "12\n";
		}

		fout << "POINT_DATA " << cp.size() << "\nVECTORS type float\n";
		for (uint i = 0; i<cp.size(); i++)
		{
			fout << cp[i].type << " " << cp[i].sharp << " 0\n";
		}

		//fout << "POINT_DATA " << cp.size() << "\nSCALARS type float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<cp.size(); i++)
		//{
		//	fout << cp[i].type << "\n";
		//}
		//fout << "\nCELL_DATA " << hmesh[lev].size() << "\nSCALARS eact float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<hmesh[lev].size(); i++)
		//{
		//	//fout<<hmesh[lev][i].act<<"\n";
		//	fout << hmesh[lev][i].type << "\n";
		//}

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline_3D::OutputFace(string fn)
{
	string fname(fn + "_face.vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << cp.size() << " float\n";
		for (uint i = 0; i<cp.size(); i++)
		{
			fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << tmface.size() << " " << 5 * tmface.size() << '\n';
		for (uint i = 0; i<tmface.size(); i++)
		{
			fout << "4 ";
			for (int j = 0; j<4; j++)
			{
				fout << tmface[i].cnct[j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << tmface.size() << '\n';
		for (uint i = 0; i<tmface.size(); i++)
		{
			fout << "9\n";
		}

		//fout << "POINT_DATA " << hcp[lev].size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<hcp[lev].size(); i++)
		//{
		//	fout << hcp[lev][i].act << "\n";
		//}

		fout << "\nCELL_DATA " << tmface.size() << "\nSCALARS type float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<tmface.size(); i++)
		{
			fout << tmface[i].type << "\n";
		}

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline_3D::OutputEdge(string fn)
{
	string fname(fn + "_edge.vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << cp.size() << " float\n";
		for (uint i = 0; i<cp.size(); i++)
		{
			fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << tmedge.size() << " " << 3 * tmedge.size() << '\n';
		for (uint i = 0; i<tmedge.size(); i++)
		{
			fout << "2 ";
			for (int j = 0; j<2; j++)
			{
				fout << tmedge[i].pt[j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << tmedge.size() << '\n';
		for (uint i = 0; i<tmedge.size(); i++)
		{
			fout << "3\n";
		}

		fout << "POINT_DATA " << cp.size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<cp.size(); i++)
		{
			fout << cp[i].sharp << "\n";
			//fout << hcp[lev][i].bcxp << "\n";
		}
		fout << "\nCELL_DATA " << tmedge.size() << "\nSCALARS type float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<tmedge.size(); i++)
		{
			fout << tmedge[i].sharp << "\n";
		}
		//fout << "\nCELL_DATA " << tmedge.size() << "\nVECTORS type float\n";
		//for (uint i = 0; i<tmedge.size(); i++)
		//{
		//	fout << tmedge[i].type << " " << tmedge[i].sharp << " 0\n";
		//}

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline_3D::run_Global(int nref, string fn)
{
	for (int i = 0; i < nref; i++)
	{
		cout << "Refine step: " << i << "\n";
		Global_Subdivide();
		//Global_Subdivide_Simple();

		stringstream ss;
		ss << (i + 1);
		//OutputCM(fn + ss.str());
		//OutputFace(fn + ss.str());
		//OutputEdge(fn + ss.str());
	}
	cout << "Building basis function...\n";
	BuildBasisFunction();
}

void TruncatedTspline_3D::SetDomainRange(double xy[3][2], double nm[3], double& a)
{
	double tmp(1.e6);
	double x_range[3][2] = { { tmp, -tmp }, { tmp, -tmp }, { tmp, -tmp } };
	uint i, j;
	for (i = 0; i <hcp[0].size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (hcp[0][i].coor[j] < x_range[j][0]) x_range[j][0] = hcp[0][i].coor[j];
			if (hcp[0][i].coor[j] > x_range[j][1]) x_range[j][1] = hcp[0][i].coor[j];
		}
	}
	double xh[3] = { x_range[0][1] - x_range[0][0], x_range[1][1] - x_range[1][0], x_range[2][1] - x_range[2][0] };

	double nmtmp[3] = { -xh[1] * xh[2], -xh[0] * xh[2], 2.*xh[0] * xh[1] };
	tmp = sqrt(nmtmp[0] * nmtmp[0] + nmtmp[1] * nmtmp[1] + nmtmp[2] * nmtmp[2]);
	nm[0] = nmtmp[0] / tmp; nm[1] = nmtmp[1] / tmp; nm[2] = nmtmp[2] / tmp;
	//nm[0] = 1.; nm[1] = 0.; nm[2] = 0.;//rod

	//nm[0] = .5; nm[1] = .5; nm[2] = .5;//cube

	xy[0][0] = x_range[0][0]; xy[1][0] = x_range[1][0]; xy[2][0] = x_range[2][0];
	xy[0][1] = x_range[0][1]; xy[1][1] = x_range[1][1]; xy[2][1] = x_range[2][1];

	dmrg[0][0] = xy[0][0]; dmrg[0][1] = xy[0][1];
	dmrg[1][0] = xy[1][0]; dmrg[1][1] = xy[1][1];
	dmrg[2][0] = xy[2][0]; dmrg[2][1] = xy[2][1];
	nmpl[0] = nm[0]; nmpl[1] = nm[1]; nmpl[2] = nm[2];
	acoef = a / sqrt(xh[0] * xh[0] + xh[1] * xh[1] + xh[2] * xh[2]);//regulation
	a = acoef;
	dmlen[0] = dmrg[0][1] - dmrg[0][0];
	dmlen[1] = dmrg[1][1] - dmrg[1][0];
	dmlen[2] = dmrg[2][1] - dmrg[2][0];
}

void TruncatedTspline_3D::PipelineDataProcess(string fn_in, string fn_out)
{
	//read hex vtk
	vector<double> disp;
	string fname(fn_in + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i<npts; i++)
		{
			//cp[i].act = 1;
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i<neles; i++)
		{
			//tmesh[i].act = 1;
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3] >>
				tmesh[i].cnct[4] >> tmesh[i].cnct[5] >> tmesh[i].cnct[6] >> tmesh[i].cnct[7];
		}
		for (int i = 0; i<neles+5; i++) getline(fin, stmp);
		disp.resize(npts);
		for (int i = 0; i < npts; i++) fin >> disp[i];
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	//remove redundant points
	vector<int> flag(cp.size(), 0);
	for (uint i = 0; i < tmesh.size(); i++)
	{
		for (int j = 0; j < 8; j++) flag[tmesh[i].cnct[j]] = 1;
	}
	vector<int> pid_new(cp.size(), -1);
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (flag[i] == 1)
		{
			pid_new[i] = count++;
		}
	}
	for (uint i = 0; i < tmesh.size(); i++)
	{
		for (int j = 0; j < 8; j++) tmesh[i].cnct[j] = pid_new[tmesh[i].cnct[j]];
	}

	//output
	string fn1(fn_out + ".vtk");
	ofstream fout;
	fout.open(fn1.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << count << " float\n";
		for (uint i = 0; i<cp.size(); i++)
		{
			if (flag[i] == 1)
			{
				fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
			}
		}
		fout << "\nCELLS " << tmesh.size() << " " << 9 * tmesh.size() << '\n';
		for (uint i = 0; i<tmesh.size(); i++)
		{
			fout << "8 ";
			for (int j = 0; j<8; j++)
			{
				fout << tmesh[i].cnct[j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << tmesh.size() << '\n';
		for (uint i = 0; i<tmesh.size(); i++)
		{
			fout << "12\n";
		}
		fout << "POINT_DATA " << count << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<cp.size(); i++)
		{
			if (flag[i] == 1)
			{
				fout << disp[i] << "\n";
			}
		}
		//fout << "\nCELL_DATA " << hmesh[lev].size() << "\nSCALARS eact float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<hmesh[lev].size(); i++)
		//{
		//	//fout<<hmesh[lev][i].act<<"\n";
		//	fout << hmesh[lev][i].type << "\n";
		//}

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fn1 << "!\n";
	}
}

void TruncatedTspline_3D::PipelineBezierExtract(vector<BezierElement3D>& bzmesh)
{
	bzmesh.clear();
	uint i, j, k, k1, k2;
	int loc(0);
	vector<vector<int>> aloc(hcp.size());
	for (i = 0; i < hcp.size(); i++)
	{
		aloc[i].resize(hcp[i].size(), -1);
		for (j = 0; j < hcp[i].size(); j++)
		{
			if (hcp[i][j].act == 1)
			{
				aloc[i][j] = loc++;
			}
		}
	}
	for (i = 0; i < hmesh.size(); i++)
	{
		for (j = 0; j < hmesh[i].size(); j++)
		{
			if (hmesh[i][j].act == 1 /*&& hmesh[i][j].type!=1*/)
			{
				BezierElement3D bztmp;
				bztmp.prt[0] = i; bztmp.prt[1] = j;
				bztmp.trun = hmesh[i][j].trun;
				if (hmesh[i][j].type == 1) bztmp.type = 1;
				if (hmesh[i][j].trun == 0)
				{
					bztmp.IEN.resize(hmesh[i][j].IEN.size());
					bztmp.cmat.resize(hmesh[i][j].IEN.size(), vector<double>(64));
					for (k = 0; k < hmesh[i][j].IEN.size(); k++)
					{
						bztmp.IEN[k] = aloc[i][hmesh[i][j].IEN[k]];
						for (k1 = 0; k1 < 64; k1++)
						{
							bztmp.cmat[k][k1] = hmesh[i][j].bemat[k][k1];
							bztmp.pts[k1][0] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[0];
							bztmp.pts[k1][1] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[1];
							bztmp.pts[k1][2] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[2];
						}
					}
				}
				else
				{
					bztmp.IEN.resize(hmesh[i][j].IEN_act.size());
					bztmp.cmat.resize(hmesh[i][j].IEN_act.size(), vector<double>(64));
					for (k = 0; k < hmesh[i][j].IEN_act.size(); k++)
					{
						int lev(hmesh[i][j].IEN_act[k][0]);
						int pid(hmesh[i][j].IEN_act[k][1]);
						if (aloc[lev][pid] == -1)
						{
							cout << "wrong aloc!\n";
							getchar();
						}
						bztmp.IEN[k] = aloc[lev][pid];
						for (k1 = 0; k1 < 64; k1++)
						{
							bztmp.cmat[k][k1] = 0.;
							for (k2 = 0; k2 < hmesh[i][j].IEN.size(); k2++)
							{
								bztmp.cmat[k][k1] += hmesh[i][j].tmat[k][k2] * hmesh[i][j].bemat[k2][k1];
							}
						}
					}
					for (k = 0; k < hmesh[i][j].IEN.size(); k++)
					{
						for (k1 = 0; k1 < 64; k1++)
						{
							bztmp.pts[k1][0] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[0];
							bztmp.pts[k1][1] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[1];
							bztmp.pts[k1][2] += hmesh[i][j].bemat[k][k1] * hcp[i][hmesh[i][j].IEN[k]].coor[2];
						}
					}
				}
				bzmesh.push_back(bztmp);
			}
		}
	}
}

void TruncatedTspline_3D::CreateBsplines(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh)
{
	int nsp(3);
	nsp *= 2;
	nsp *= 2;
	int nex[3] = { nsp, nsp, nsp };//3 by 3 by 3 cube
	int npx[3] = { nex[0] + 3, nex[1] + 3, nex[2] + 3 };
	int ncp(npx[0] * npx[1] * npx[2]), nel(nex[0] * nex[1] * nex[2]);
	double dmx[3][2] = { { 0., 1. }, { 0., 1. }, { 0., 1. } };
	double lenx[3] = { (dmx[0][1] - dmx[0][0])/double(npx[0]-1), (dmx[1][1] - dmx[1][0])/double(npx[1]-1), (dmx[2][1] - dmx[2][0])/double(npx[2]-1) };
	vector<vector<double>> ku(3);
	cp.resize(ncp);
	IDBC.resize(ncp, -1);
	gh.resize(ncp, 0.);
	bzmesh.resize(nel);
	uint i, j, k, loc(0), i0, j0, k0;
	for (i = 0; i < 4; i++)
	{
		ku[0].push_back(0.);
		ku[1].push_back(0.);
		ku[2].push_back(0.);
	}
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < nex[i]; j++) ku[i].push_back(double(j + 1));
	}
	for (i = 0; i < 3; i++)
	{
		ku[0].push_back(double(nex[0]));
		ku[1].push_back(double(nex[1]));
		ku[2].push_back(double(nex[2]));
	}
	loc = 0;
	int count(0);
	for (k = 0; k < npx[2]; k++)
	{
		for (j = 0; j < npx[1]; j++)
		{
			for (i = 0; i < npx[0]; i++)
			{
				//cp[loc].act = 1;
				cp[loc].coor[0] = dmx[0][0] + double(i) * lenx[0];
				cp[loc].coor[1] = dmx[1][0] + double(j) * lenx[1];
				cp[loc].coor[2] = dmx[2][0] + double(k) * lenx[2];
				if (i == 0 || i == npx[0] - 1 || j == 0 || j == npx[1] - 1 || k == 0 || k == npx[2] - 1)
				{
					IDBC[loc] = -1;
					gh[loc] = SpecifyDirichBC(cp[loc].coor);
				}
				else
				{
					IDBC[loc] = count++;
				}
				loc++;
			}
		}
	}
	loc = 0;
	for (k = 0; k < nex[2]; k++)
	{
		for (j = 0; j < nex[1]; j++)
		{
			for (i = 0; i < nex[0]; i++)
			{
				bzmesh[loc].IEN.resize(64);
				int loc0(0);
				for (k0 = 0; k0 < 4; k0++)
				{
					for (j0 = 0; j0 < 4; j0++)
					{
						for (i0 = 0; i0 < 4; i0++)
						{
							bzmesh[loc].IEN[loc0] = (k + k0)*npx[0] * npx[1] + (j + j0)*npx[0] + i + i0;
							loc0++;
						}
					}
				}
				bzmesh[loc].cmat.resize(64, vector<double>(64, 0.));
				vector<double> ku0(ku[0].begin() + i, ku[0].begin() + i + 8);
				vector<double> kv0(ku[1].begin() + j, ku[1].begin() + j + 8);
				vector<double> kw0(ku[2].begin() + k, ku[2].begin() + k + 8);
				vector<double> ku1, kv1, kw1;
				array<double, 2> ktu = { double(i), double(i + 1) };
				array<double, 2> ktv = { double(j), double(j + 1) };
				array<double, 2> ktw = { double(k), double(k + 1) };
				vector<vector<double>> Tu, Tv, Tw;
				int iloc[3];
				BezierInsertKnots(ku0, ktu, ku1);
				BezierInsertKnots(kv0, ktv, kv1);
				BezierInsertKnots(kw0, ktw, kw1);
				TMatrix(ku0, ku1, 3, Tu);
				TMatrix(kv0, kv1, 3, Tv);
				TMatrix(kw0, kw1, 3, Tw);
				for (i0 = 0; i0 < ku1.size() - 1; i0++)
				{
					if (ku1[i0] == double(i) && ku1[i0 + 1] == double(i + 1))
					{
						iloc[0] = i0 - 3; break;
					}
				}
				for (i0 = 0; i0 < kv1.size() - 1; i0++)
				{
					if (kv1[i0] == double(j) && kv1[i0 + 1] == double(j + 1))
					{
						iloc[1] = i0 - 3; break;
					}
				}
				for (i0 = 0; i0 < kw1.size() - 1; i0++)
				{
					if (kw1[i0] == double(k) && kw1[i0 + 1] == double(k + 1))
					{
						iloc[2] = i0 - 3; break;
					}
				}
				loc0 = 0;//Bspline
				for (k0 = 0; k0 < 4; k0++)
				{
					for (j0 = 0; j0 < 4; j0++)
					{
						for (i0 = 0; i0 < 4; i0++)
						{
							int loc1(0);//Bezier
							for (int k1 = 0; k1 < 4; k1++)
							{
								for (int j1 = 0; j1 < 4; j1++)
								{
									for (int i1 = 0; i1 < 4; i1++)
									{
										bzmesh[loc].cmat[loc0][loc1] = Tu[iloc[0] + i1][i0] * Tv[iloc[1] + j1][j0] * Tw[iloc[2] + k1][k0];
										bzmesh[loc].pts[loc1][0] += bzmesh[loc].cmat[loc0][loc1] * cp[bzmesh[loc].IEN[loc0]].coor[0];
										bzmesh[loc].pts[loc1][1] += bzmesh[loc].cmat[loc0][loc1] * cp[bzmesh[loc].IEN[loc0]].coor[1];
										bzmesh[loc].pts[loc1][2] += bzmesh[loc].cmat[loc0][loc1] * cp[bzmesh[loc].IEN[loc0]].coor[2];
										loc1++;
									}
								}
							}
							loc0++;
						}
					}
				}
				//cout << loc << "\n";
				//for (int pid = 0; pid < bzmesh[loc].pts.size(); pid++)
				//{
				//	cout << bzmesh[loc].pts[pid][0] << " " << bzmesh[loc].pts[pid][1] << " " << bzmesh[loc].pts[pid][2] << "\n";
				//	getchar();
				//}
				//cout << "done " << loc << "\n";
				loc++;
			}
		}
	}
}

void TruncatedTspline_3D::CreateBsplines(int nex[3])
{
	int npx[3] = { nex[0] + 3, nex[1] + 3, nex[2] + 3 };
	int ncp(npx[0] * npx[1] * npx[2]);
	int nel(nex[0] * nex[1] * nex[2]);
	double dmx[3][2] = { { 0., 1. }, { 0., 1. }, { 0., 1. } };
	double lenx[3] = { (dmx[0][1] - dmx[0][0]) / double(npx[0] - 1), (dmx[1][1] - dmx[1][0]) / double(npx[1] - 1), (dmx[2][1] - dmx[2][0]) / double(npx[2] - 1) };
	kvec.resize(3);
	cp.resize(ncp);
	tmesh.resize(nel);
	uint i, j, k, loc(0);
	for (i = 0; i < 4; i++)
	{
		kvec[0].push_back(0.);
		kvec[1].push_back(0.);
		kvec[2].push_back(0.);
	}
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < nex[i]; j++) kvec[i].push_back(double(j + 1));
	}
	for (i = 0; i < 3; i++)
	{
		kvec[0].push_back(double(nex[0]));
		kvec[1].push_back(double(nex[1]));
		kvec[2].push_back(double(nex[2]));
	}
	loc = 0;
	int count(0);
	for (k = 0; k < npx[2]; k++)
	{
		for (j = 0; j < npx[1]; j++)
		{
			for (i = 0; i < npx[0]; i++)
			{
				cp[loc].act = 1;
				if (i == 0 || i == npx[0] - 1 || j == 0 || j == npx[1] - 1 || k == 0 || k == npx[2] - 1)
				{
					cp[loc].type = 1;
				}
				cp[loc].coor[0] = dmx[0][0] + double(i) * lenx[0];
				cp[loc].coor[1] = dmx[1][0] + double(j) * lenx[1];
				cp[loc].coor[2] = dmx[2][0] + double(k) * lenx[2];
				loc++;
			}
		}
	}
}

void TruncatedTspline_3D::Bsplines_Refine()
{
	int nex0[3] = { kvec[0].size() - 7, kvec[1].size() - 7, kvec[2].size() - 7 };//# elements after refinement
	int npx0[3] = { nex0[0] + 3, nex0[1] + 3, nex0[2] + 3 };
	int nex[3] = { 2*(kvec[0].size() - 7), 2*(kvec[1].size() - 7), 2*(kvec[2].size() - 7) };//# elements after refinement
	int npx[3] = { nex[0] + 3, nex[1] + 3, nex[2] + 3 };
	int ncp0(npx0[0] * npx0[1] * npx0[2]), ncp(npx[0] * npx[1] * npx[2]);
	int nel0(nex0[0] * nex0[1] * nex0[2]), nel(nex[0] * nex[1] * nex[2]);
	uint i, j, k, i0, j0, k0, loc(0);
	vector<vector<double>> ku(3);//knot vectors after refinement
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < kvec[i].size() - 1; j++)
		{
			ku[i].push_back(kvec[i][j]);
			if (kvec[i][j] < kvec[i][j + 1])
			{
				ku[i].push_back((kvec[i][j] + kvec[i][j + 1]) / 2.);
			}
		}
		ku[i].push_back(kvec[i].back());
	}
	vector<vector<double>> Tu, Tv, Tw;
	TMatrix(kvec[0], ku[0], 3, Tu);
	TMatrix(kvec[1], ku[1], 3, Tv);
	TMatrix(kvec[2], ku[2], 3, Tw);
	vector<array<double, 3>> pts(ncp);//new points

	loc = 0;//new id
	for (k = 0; k < npx[2]; k++)
	{
		for (j = 0; j < npx[1]; j++)
		{
			for (i = 0; i < npx[0]; i++)
			{
				int loc0(0);//old id
				for (k0 = 0; k0 < npx0[2]; k0++)
				{
					for (j0 = 0; j0 < npx0[1]; j0++)
					{
						for (i0 = 0; i0 < npx0[0]; i0++)
						{
							if (Tu[i][i0] != 0. && Tv[j][j0] != 0. && Tw[k][k0]!=0.)
							{
								double ctmp = Tu[i][i0] * Tv[j][j0] * Tw[k][k0];
								pts[loc][0] += ctmp*cp[loc0].coor[0];
								pts[loc][1] += ctmp*cp[loc0].coor[1];
								pts[loc][2] += ctmp*cp[loc0].coor[2];
							}
							loc0++;
						}
					}
				}
				loc++;
			}
		}
	}

	for (i = 0; i < 3; i++)
	{
		kvec[i].clear();
	}
	for (i = 0; i < 4; i++)
	{
		kvec[0].push_back(0.);
		kvec[1].push_back(0.);
		kvec[2].push_back(0.);
	}
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < nex[i]; j++) kvec[i].push_back(double(j + 1));
	}
	for (i = 0; i < 3; i++)
	{
		kvec[0].push_back(double(nex[0]));
		kvec[1].push_back(double(nex[1]));
		kvec[2].push_back(double(nex[2]));
	}
	cp.clear();
	cp.resize(ncp);
	loc = 0;
	for (k = 0; k < npx[2]; k++)
	{
		for (j = 0; j < npx[1]; j++)
		{
			for (i = 0; i < npx[0]; i++)
			{
				cp[loc].act = 1;
				cp[loc].type = 0;
				if (i == 0 || i == npx[0] - 1 || j == 0 || j == npx[1] - 1 || k == 0 || k == npx[2] - 1)
				{
					cp[loc].type = 1;
				}
				cp[loc].coor[0] = pts[loc][0];
				cp[loc].coor[1] = pts[loc][1];
				cp[loc].coor[2] = pts[loc][2];
				loc++;
			}
		}
	}
	//for (i = 0; i < pts.size(); i++)
	//{
	//	cp[i].coor[0] = pts[i][0];
	//	cp[i].coor[1] = pts[i][1];
	//	cp[i].coor[2] = pts[i][2];
	//}
}

void TruncatedTspline_3D::Bspline_BezierExtract(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh)
{
	int nex[3] = { kvec[0].size() - 7, kvec[1].size() - 7, kvec[2].size() - 7 };//# elements after refinement
	int npx[3] = { nex[0] + 3, nex[1] + 3, nex[2] + 3 };
	int ncp(npx[0] * npx[1] * npx[2]);
	int nel(nex[0] * nex[1] * nex[2]);
	vector<vector<double>> ku(3);
	IDBC.resize(ncp, -1);
	gh.resize(ncp, 0.);
	bzmesh.resize(nel);
	uint i, j, k, loc(0), i0, j0, k0;
	for (i = 0; i < 3; i++)
	{
		ku[i].resize(kvec[i].size());
		for (j = 0; j < kvec[i].size(); j++)
		{
			ku[i][j] = kvec[i][j];
		}
	}
	loc = 0;
	int count(0);
	for (k = 0; k < npx[2]; k++)
	{
		for (j = 0; j < npx[1]; j++)
		{
			for (i = 0; i < npx[0]; i++)
			{
				if (i == 0 || i == npx[0] - 1 || j == 0 || j == npx[1] - 1 || k == 0 || k == npx[2] - 1)
				{
					IDBC[loc] = -1;
					gh[loc] = SpecifyDirichBC(cp[loc].coor);
				}
				else
				{
					IDBC[loc] = count++;
				}
				loc++;
			}
		}
	}
	loc = 0;
	for (k = 0; k < nex[2]; k++)
	{
		for (j = 0; j < nex[1]; j++)
		{
			for (i = 0; i < nex[0]; i++)
			{
				bzmesh[loc].IEN.resize(64);
				int loc0(0);
				for (k0 = 0; k0 < 4; k0++)
				{
					for (j0 = 0; j0 < 4; j0++)
					{
						for (i0 = 0; i0 < 4; i0++)
						{
							bzmesh[loc].IEN[loc0] = (k + k0)*npx[0] * npx[1] + (j + j0)*npx[0] + i + i0;
							loc0++;
						}
					}
				}
				bzmesh[loc].cmat.resize(64, vector<double>(64, 0.));
				vector<double> ku0(ku[0].begin() + i, ku[0].begin() + i + 8);
				vector<double> kv0(ku[1].begin() + j, ku[1].begin() + j + 8);
				vector<double> kw0(ku[2].begin() + k, ku[2].begin() + k + 8);
				vector<double> ku1, kv1, kw1;
				array<double, 2> ktu = { double(i), double(i + 1) };
				array<double, 2> ktv = { double(j), double(j + 1) };
				array<double, 2> ktw = { double(k), double(k + 1) };
				vector<vector<double>> Tu, Tv, Tw;
				int iloc[3];
				BezierInsertKnots(ku0, ktu, ku1);
				BezierInsertKnots(kv0, ktv, kv1);
				BezierInsertKnots(kw0, ktw, kw1);
				TMatrix(ku0, ku1, 3, Tu);
				TMatrix(kv0, kv1, 3, Tv);
				TMatrix(kw0, kw1, 3, Tw);
				for (i0 = 0; i0 < ku1.size() - 1; i0++)
				{
					if (ku1[i0] == double(i) && ku1[i0 + 1] == double(i + 1))
					{
						iloc[0] = i0 - 3; break;
					}
				}
				for (i0 = 0; i0 < kv1.size() - 1; i0++)
				{
					if (kv1[i0] == double(j) && kv1[i0 + 1] == double(j + 1))
					{
						iloc[1] = i0 - 3; break;
					}
				}
				for (i0 = 0; i0 < kw1.size() - 1; i0++)
				{
					if (kw1[i0] == double(k) && kw1[i0 + 1] == double(k + 1))
					{
						iloc[2] = i0 - 3; break;
					}
				}
				loc0 = 0;//Bspline
				for (k0 = 0; k0 < 4; k0++)
				{
					for (j0 = 0; j0 < 4; j0++)
					{
						for (i0 = 0; i0 < 4; i0++)
						{
							int loc1(0);//Bezier
							for (int k1 = 0; k1 < 4; k1++)
							{
								for (int j1 = 0; j1 < 4; j1++)
								{
									for (int i1 = 0; i1 < 4; i1++)
									{
										bzmesh[loc].cmat[loc0][loc1] = Tu[iloc[0] + i1][i0] * Tv[iloc[1] + j1][j0] * Tw[iloc[2] + k1][k0];
										bzmesh[loc].pts[loc1][0] += bzmesh[loc].cmat[loc0][loc1] * cp[bzmesh[loc].IEN[loc0]].coor[0];
										bzmesh[loc].pts[loc1][1] += bzmesh[loc].cmat[loc0][loc1] * cp[bzmesh[loc].IEN[loc0]].coor[1];
										bzmesh[loc].pts[loc1][2] += bzmesh[loc].cmat[loc0][loc1] * cp[bzmesh[loc].IEN[loc0]].coor[2];
										loc1++;
									}
								}
							}
							loc0++;
						}
					}
				}
				loc++;
			}
		}
	}
}

void TruncatedTspline_3D::Bspline_BezierExtract_fit(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh)
{
	int nex[3] = { kvec[0].size() - 7, kvec[1].size() - 7, kvec[2].size() - 7 };//# elements after refinement
	int npx[3] = { nex[0] + 3, nex[1] + 3, nex[2] + 3 };
	int ncp(npx[0] * npx[1] * npx[2]);
	int nel(nex[0] * nex[1] * nex[2]);
	vector<vector<double>> ku(3);
	//IDBC.resize(ncp, -1);
	//gh.resize(ncp, 0.);
	bzmesh.resize(nel);
	uint i, j, k, loc(0), i0, j0, k0;
	for (i = 0; i < 3; i++)
	{
		ku[i].resize(kvec[i].size());
		for (j = 0; j < kvec[i].size(); j++)
		{
			ku[i][j] = kvec[i][j];
		}
	}
	/*loc = 0;
	int count(0);
	for (k = 0; k < npx[2]; k++)
	{
		for (j = 0; j < npx[1]; j++)
		{
			for (i = 0; i < npx[0]; i++)
			{
				if (i == 0 || i == npx[0] - 1 || j == 0 || j == npx[1] - 1 || k == 0 || k == npx[2] - 1)
				{
					IDBC[loc] = -1;
					gh[loc] = SpecifyDirichBC(cp[loc].coor);
				}
				else
				{
					IDBC[loc] = count++;
				}
				loc++;
			}
		}
	}*/
	loc = 0;
	for (k = 0; k < nex[2]; k++)
	{
		for (j = 0; j < nex[1]; j++)
		{
			for (i = 0; i < nex[0]; i++)
			{
				if (i == 0 || i == nex[0] - 1 || j == 0 || j == nex[1] - 1 || k == 0 || k == nex[2] - 1)
				{
					bzmesh[loc].type = 1;
				}
				bzmesh[loc].IEN.resize(64);
				int loc0(0);
				for (k0 = 0; k0 < 4; k0++)
				{
					for (j0 = 0; j0 < 4; j0++)
					{
						for (i0 = 0; i0 < 4; i0++)
						{
							bzmesh[loc].IEN[loc0] = (k + k0)*npx[0] * npx[1] + (j + j0)*npx[0] + i + i0;
							loc0++;
						}
					}
				}
				bzmesh[loc].cmat.resize(64, vector<double>(64, 0.));
				vector<double> ku0(ku[0].begin() + i, ku[0].begin() + i + 8);
				vector<double> kv0(ku[1].begin() + j, ku[1].begin() + j + 8);
				vector<double> kw0(ku[2].begin() + k, ku[2].begin() + k + 8);
				vector<double> ku1, kv1, kw1;
				array<double, 2> ktu = { double(i), double(i + 1) };
				array<double, 2> ktv = { double(j), double(j + 1) };
				array<double, 2> ktw = { double(k), double(k + 1) };
				vector<vector<double>> Tu, Tv, Tw;
				int iloc[3];
				BezierInsertKnots(ku0, ktu, ku1);
				BezierInsertKnots(kv0, ktv, kv1);
				BezierInsertKnots(kw0, ktw, kw1);
				TMatrix(ku0, ku1, 3, Tu);
				TMatrix(kv0, kv1, 3, Tv);
				TMatrix(kw0, kw1, 3, Tw);
				for (i0 = 0; i0 < ku1.size() - 1; i0++)
				{
					if (ku1[i0] == double(i) && ku1[i0 + 1] == double(i + 1))
					{
						iloc[0] = i0 - 3; break;
					}
				}
				for (i0 = 0; i0 < kv1.size() - 1; i0++)
				{
					if (kv1[i0] == double(j) && kv1[i0 + 1] == double(j + 1))
					{
						iloc[1] = i0 - 3; break;
					}
				}
				for (i0 = 0; i0 < kw1.size() - 1; i0++)
				{
					if (kw1[i0] == double(k) && kw1[i0 + 1] == double(k + 1))
					{
						iloc[2] = i0 - 3; break;
					}
				}
				loc0 = 0;//Bspline
				for (k0 = 0; k0 < 4; k0++)
				{
					for (j0 = 0; j0 < 4; j0++)
					{
						for (i0 = 0; i0 < 4; i0++)
						{
							int loc1(0);//Bezier
							for (int k1 = 0; k1 < 4; k1++)
							{
								for (int j1 = 0; j1 < 4; j1++)
								{
									for (int i1 = 0; i1 < 4; i1++)
									{
										bzmesh[loc].cmat[loc0][loc1] = Tu[iloc[0] + i1][i0] * Tv[iloc[1] + j1][j0] * Tw[iloc[2] + k1][k0];
										bzmesh[loc].pts[loc1][0] += bzmesh[loc].cmat[loc0][loc1] * cp[bzmesh[loc].IEN[loc0]].coor[0];
										bzmesh[loc].pts[loc1][1] += bzmesh[loc].cmat[loc0][loc1] * cp[bzmesh[loc].IEN[loc0]].coor[1];
										bzmesh[loc].pts[loc1][2] += bzmesh[loc].cmat[loc0][loc1] * cp[bzmesh[loc].IEN[loc0]].coor[2];
										loc1++;
									}
								}
							}
							loc0++;
						}
					}
				}
				loc++;
			}
		}
	}
}

void TruncatedTspline_3D::FittingBC_Bsplines(const vector<BezierElement3D>& bzall, vector<int>& IDBC1, vector<double>& gh1)
{
	uint i, j, k, k1, k2;
	int nex[3] = { kvec[0].size() - 7, kvec[1].size() - 7, kvec[2].size() - 7 };
	//boundary points
	vector<int> flag(cp.size(), 0);
	for (i = 0; i < bzall.size(); i++)
	{
		if (bzall[i].type == 1)
		{
			for (j = 0; j < bzall[i].IEN.size(); j++)
			{
				flag[bzall[i].IEN[j]] = 1;
			}
		}
	}
	vector<int> aloc(cp.size(), -1);
	int loc = 0;
	for (i = 0; i < aloc.size(); i++)
	{
		if (flag[i] == 1)
		{
			aloc[i] = loc++;
		}
	}
	vector<int> IDBC(loc, -1);
	vector<double> gh(loc, 0.);
	loc = 0;
	int count(0);
	for (i = 0; i < cp.size(); i++)
	{
		if (flag[i] == 1)
		{
			if (cp[i].type == 1)
			{
				IDBC[loc] = count++;
			}
			loc++;
		}
	}
	vector<BezierElement3D> bzmesh;
	loc = 0;
	for (k = 0; k < nex[2]; k++)
	{
		for (j = 0; j < nex[1]; j++)
		{
			for (i = 0; i < nex[0]; i++)
			{
				if (bzall[loc].type == 1)
				{
					BezierElement3D bztmp;
					bztmp.type = 1;
					bool fbc[6] = { k == 0, j == 0, i == nex[0] - 1, j == nex[1] - 1, i == 0, k == nex[2] - 1 };
					for (k1 = 0; k1 < 6; k1++)
					{
						if (fbc[k1])
						{
							bztmp.bfc.push_back(k1);
						}
					}
					bztmp.IEN.resize(bzall[loc].IEN.size());
					bztmp.cmat.resize(bzall[loc].IEN.size(), vector<double>(64));
					for (k1 = 0; k1 < bzall[loc].IEN.size(); k1++)
					{
						bztmp.IEN[k1] = aloc[bzall[loc].IEN[k1]];
						for (k2 = 0; k2 < 64; k2++)
						{
							bztmp.cmat[k1][k2] = bzall[loc].cmat[k1][k2];
						}
					}
					for (k2 = 0; k2 < 64; k2++)
					{
						bztmp.pts[k2][0] = bzall[loc].pts[k2][0];
						bztmp.pts[k2][1] = bzall[loc].pts[k2][1];
						bztmp.pts[k2][2] = bzall[loc].pts[k2][2];
					}
					bzmesh.push_back(bztmp);
				}
				loc++;
			}
		}
	}

	LeastSquare ls;
	vector<double> sol;
	ls.SetProblem(IDBC, gh);
	ls.GetEqParameter(dmrg, nmpl, acoef);
	//ls.VisualizeBoundarySurface(bzmesh, cpts, "../io/complex2/rod2");
	//cout << "done output boundary surface!\n";
	//getchar();
	ls.Run_Fitting(bzmesh, "", sol);
	//ls.VisualizeBoundarySurface(bzmesh, cpts, "../io/complex2/rod2");
	//cout << "done output boundary surface!\n";
	//getchar();

	IDBC1.clear();
	gh1.clear();
	IDBC1.resize(cp.size(), -1);
	gh1.resize(cp.size(), 0.);
	loc = 0;
	count = 0;
	for (i = 0; i < cp.size(); i++)
	{
		if (flag[i] == 1)
		{
			if (cp[i].type == 1)
			{
				IDBC1[i] = -1;
				gh1[i] = sol[loc];
			}
			else
			{
				IDBC1[i] = count++;
			}
			loc++;
		}
		else
		{
			IDBC1[i] = count++;
		}
	}
}

void TruncatedTspline_3D::SetDomainRange_Bsplines(double xy[3][2], double nm[3], double& a)
{
	double tmp(1.e6);
	double x_range[3][2] = { { tmp, -tmp }, { tmp, -tmp }, { tmp, -tmp } };
	uint i, j;
	for (i = 0; i <cp.size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (cp[i].coor[j] < x_range[j][0]) x_range[j][0] = cp[i].coor[j];
			if (cp[i].coor[j] > x_range[j][1]) x_range[j][1] = cp[i].coor[j];
		}
	}
	double xh[3] = { x_range[0][1] - x_range[0][0], x_range[1][1] - x_range[1][0], x_range[2][1] - x_range[2][0] };

	double nmtmp[3] = { -xh[1] * xh[2], -xh[0] * xh[2], 2.*xh[0] * xh[1] };
	tmp = sqrt(nmtmp[0] * nmtmp[0] + nmtmp[1] * nmtmp[1] + nmtmp[2] * nmtmp[2]);
	nm[0] = nmtmp[0] / tmp; nm[1] = nmtmp[1] / tmp; nm[2] = nmtmp[2] / tmp;

	//nm[0] = .5; nm[1] = .5; nm[2] = .5;//cube

	xy[0][0] = x_range[0][0]; xy[1][0] = x_range[1][0]; xy[2][0] = x_range[2][0];
	xy[0][1] = x_range[0][1]; xy[1][1] = x_range[1][1]; xy[2][1] = x_range[2][1];

	dmrg[0][0] = xy[0][0]; dmrg[0][1] = xy[0][1];
	dmrg[1][0] = xy[1][0]; dmrg[1][1] = xy[1][1];
	dmrg[2][0] = xy[2][0]; dmrg[2][1] = xy[2][1];
	nmpl[0] = nm[0]; nmpl[1] = nm[1]; nmpl[2] = nm[2];
	acoef = a / sqrt(xh[0] * xh[0] + xh[1] * xh[1] + xh[2] * xh[2]);//regulation
	a = acoef;
	dmlen[0] = dmrg[0][1] - dmrg[0][0];
	dmlen[1] = dmrg[1][1] - dmrg[1][0];
	dmlen[2] = dmrg[2][1] - dmrg[2][0];
}