#include "Truncated_Tspline.h"
#include "KnotInsertion.h"
#include "BSplineBasis.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>
#include <iomanip>
//#include "Matlab_Solver_wap.h"
#include "SingularEval.h"
#include "LeastSquare.h"

//#define PI 3.141592654

typedef unsigned int uint;

const int SMOOTH_TYPE(0); //0 for non-negative, 1 for refinable
double SMOOTH_BETA(0.4);

const double lamda(0.26);

TruncatedTspline::TruncatedTspline()
{
	cp.clear();
	tmesh.clear();
}

void TruncatedTspline::FindPatchKnotVector(int eid, vector<double> &ku, vector<double> &kv)
{
	ku.clear();
	kv.clear();
}

void TruncatedTspline::Refine(const vector<int> &rf, int id)
{
	unsigned int i, j, k, l, m;
	unsigned int nold(cp.size()), neold(tmesh.size());
	const int p(3);
	vector<array<double, 10>> kuv;
	for (i = 0; i < rf.size(); i++)
	{
		vector<int> pid = tmesh[rf[i]].IEN;
		vector<vector<double>> knewU, knewV;
		double ktU = (cp[tmesh[rf[i]].cnct[0]].knotU[2] + cp[tmesh[rf[i]].cnct[1]].knotU[2]) / 2.;
		double ktV = (cp[tmesh[rf[i]].cnct[0]].knotV[2] + cp[tmesh[rf[i]].cnct[3]].knotV[2]) / 2.;
		double ktsU[45] = {cp[tmesh[rf[i]].cnct[0]].knotU[1], cp[tmesh[rf[i]].cnct[0]].knotU[2], ktU, cp[tmesh[rf[i]].cnct[0]].knotU[3], cp[tmesh[rf[i]].cnct[0]].knotU[4],
						   cp[tmesh[rf[i]].cnct[0]].knotU[0], cp[tmesh[rf[i]].cnct[0]].knotU[1], cp[tmesh[rf[i]].cnct[0]].knotU[2], ktU, cp[tmesh[rf[i]].cnct[0]].knotU[3],
						   cp[tmesh[rf[i]].cnct[0]].knotU[1], cp[tmesh[rf[i]].cnct[0]].knotU[2], ktU, cp[tmesh[rf[i]].cnct[0]].knotU[3], cp[tmesh[rf[i]].cnct[0]].knotU[4],
						   cp[tmesh[rf[i]].cnct[1]].knotU[1], ktU, cp[tmesh[rf[i]].cnct[1]].knotU[2], cp[tmesh[rf[i]].cnct[1]].knotU[3], cp[tmesh[rf[i]].cnct[1]].knotU[4],
						   cp[tmesh[rf[i]].cnct[0]].knotU[1], cp[tmesh[rf[i]].cnct[0]].knotU[2], ktU, cp[tmesh[rf[i]].cnct[0]].knotU[3], cp[tmesh[rf[i]].cnct[0]].knotU[4],

						   cp[tmesh[rf[i]].cnct[0]].knotU[0], cp[tmesh[rf[i]].cnct[0]].knotU[1], cp[tmesh[rf[i]].cnct[0]].knotU[2], ktU, cp[tmesh[rf[i]].cnct[0]].knotU[3],
						   cp[tmesh[rf[i]].cnct[1]].knotU[1], ktU, cp[tmesh[rf[i]].cnct[1]].knotU[2], cp[tmesh[rf[i]].cnct[1]].knotU[3], cp[tmesh[rf[i]].cnct[1]].knotU[4],
						   cp[tmesh[rf[i]].cnct[2]].knotU[1], ktU, cp[tmesh[rf[i]].cnct[2]].knotU[2], cp[tmesh[rf[i]].cnct[2]].knotU[3], cp[tmesh[rf[i]].cnct[2]].knotU[4],
						   cp[tmesh[rf[i]].cnct[3]].knotU[0], cp[tmesh[rf[i]].cnct[3]].knotU[1], cp[tmesh[rf[i]].cnct[3]].knotU[2], ktU, cp[tmesh[rf[i]].cnct[3]].knotU[3]};

		double ktsV[45] = {cp[tmesh[rf[i]].cnct[0]].knotV[0], cp[tmesh[rf[i]].cnct[0]].knotV[1], cp[tmesh[rf[i]].cnct[0]].knotV[2], ktV, cp[tmesh[rf[i]].cnct[0]].knotV[3],
						   cp[tmesh[rf[i]].cnct[0]].knotV[1], cp[tmesh[rf[i]].cnct[0]].knotV[2], ktV, cp[tmesh[rf[i]].cnct[0]].knotV[3], cp[tmesh[rf[i]].cnct[0]].knotV[4],
						   cp[tmesh[rf[i]].cnct[0]].knotV[1], cp[tmesh[rf[i]].cnct[0]].knotV[2], ktV, cp[tmesh[rf[i]].cnct[0]].knotV[3], cp[tmesh[rf[i]].cnct[0]].knotV[4],
						   cp[tmesh[rf[i]].cnct[0]].knotV[1], cp[tmesh[rf[i]].cnct[0]].knotV[2], ktV, cp[tmesh[rf[i]].cnct[0]].knotV[3], cp[tmesh[rf[i]].cnct[0]].knotV[4],
						   cp[tmesh[rf[i]].cnct[2]].knotV[1], ktV, cp[tmesh[rf[i]].cnct[2]].knotV[2], cp[tmesh[rf[i]].cnct[2]].knotV[3], cp[tmesh[rf[i]].cnct[2]].knotV[4],

						   cp[tmesh[rf[i]].cnct[0]].knotV[0], cp[tmesh[rf[i]].cnct[0]].knotV[1], cp[tmesh[rf[i]].cnct[0]].knotV[2], ktV, cp[tmesh[rf[i]].cnct[0]].knotV[3],
						   cp[tmesh[rf[i]].cnct[1]].knotV[1], ktV, cp[tmesh[rf[i]].cnct[1]].knotV[2], cp[tmesh[rf[i]].cnct[1]].knotV[3], cp[tmesh[rf[i]].cnct[1]].knotV[4],
						   cp[tmesh[rf[i]].cnct[2]].knotV[1], ktV, cp[tmesh[rf[i]].cnct[2]].knotV[2], cp[tmesh[rf[i]].cnct[2]].knotV[3], cp[tmesh[rf[i]].cnct[2]].knotV[4],
						   cp[tmesh[rf[i]].cnct[3]].knotV[0], cp[tmesh[rf[i]].cnct[3]].knotV[1], cp[tmesh[rf[i]].cnct[3]].knotU[2], ktV, cp[tmesh[rf[i]].cnct[3]].knotV[3]};

		for (j = 0; j < 45; j += 5)
		{
			vector<double> tmpU(5), tmpV(5);
			tmpU.assign(ktsU + j, ktsU + j + 5);
			tmpV.assign(ktsV + j, ktsV + j + 5);
			knewU.push_back(tmpU);
			knewV.push_back(tmpV);

			Vertex ptmp;
			ptmp.knotU[0] = ktsU[j];
			ptmp.knotU[1] = ktsU[j + 1];
			ptmp.knotU[2] = ktsU[j + 2];
			ptmp.knotU[3] = ktsU[j + 3];
			ptmp.knotU[4] = ktsU[j + 4];
			ptmp.knotV[0] = ktsV[j];
			ptmp.knotV[1] = ktsV[j + 1];
			ptmp.knotV[2] = ktsV[j + 2];
			ptmp.knotV[3] = ktsV[j + 3];
			ptmp.knotV[4] = ktsV[j + 4];
			cp.push_back(ptmp);
		}
		vector<double> ku(5), kv(5);
		for (j = 0; j < 9; j++)
		{
			for (k = 0; k < pid.size(); k++)
			{
				if (knewU[j][0] >= cp[pid[k]].knotU[0] && knewU[j][4] <= cp[pid[k]].knotU[4] && knewV[j][0] >= cp[pid[k]].knotV[0] && knewV[j][4] <= cp[pid[k]].knotV[4])
				{
					//3 cases
					vector<double> ku1, kv1;
					vector<vector<double>> Tu, Tv;
					ku.assign(cp[pid[k]].knotU, cp[pid[k]].knotU + 5);
					kv.assign(cp[pid[k]].knotV, cp[pid[k]].knotV + 5);
					if (knewU[j] != ku && knewV[j] == kv) //case 1
					{
						Tu.clear();
						Tv.clear();
						vector<double> insU(1, ktU);
						InsertKnots(ku, insU, ku1);
						TMatrix(ku, ku1, p, Tu);
						kv1 = kv;
						Tv.assign(1, vector<double>(1, 1.));
					}
					else if (knewV[j] != kv && knewU[j] == ku) //case 2
					{
						Tu.clear();
						Tv.clear();
						vector<double> insV(1, ktV);
						InsertKnots(kv, insV, kv1);
						TMatrix(kv, kv1, p, Tv);
						ku1 = ku;
						Tu.assign(1, vector<double>(1, 1.));
					}
					else //case 3
					{
						Tu.clear();
						Tv.clear();
						vector<double> insU(1, ktU);
						vector<double> insV(1, ktV);
						InsertKnots(ku, insU, ku1);
						InsertKnots(kv, insV, kv1);
						TMatrix(ku, ku1, p, Tu);
						TMatrix(kv, kv1, p, Tv);
					}
					for (l = 0; l < kv1.size() - 4; l++)
					{
						for (m = 0; m < ku1.size() - 4; m++)
						{
							vector<double> tmpU(ku1.begin() + m, ku1.begin() + m + 5), tmpV(kv1.begin() + l, kv1.begin() + l + 5);
							if (tmpU == knewU[j] && tmpV == knewV[j])
							{
								double coef = Tu[m][0] * Tv[l][0];
								cp[pid[k]].tbf.push_back(nold + j);
								cp[pid[k]].tc.push_back(coef);
								cp[nold + j].coor[0] += coef * cp[pid[k]].coor[0];
								cp[nold + j].coor[1] += coef * cp[pid[k]].coor[1];
								cp[nold + j].coor[2] += coef * cp[pid[k]].coor[2];
							}
						}
					}
				}
			}
		}

		//update corner control points
		for (int ii = 0; ii < 4; ii++)
		{
			for (int jj = 0; jj < 3; jj++)
			{
				cp[tmesh[rf[i]].cnct[ii]].coor[jj] = cp[nold + 5 + ii].coor[jj];
			}
			for (int jj = 0; jj < 5; jj++)
			{
				cp[tmesh[rf[i]].cnct[ii]].knotU[jj] = cp[nold + 5 + ii].knotU[jj];
				cp[tmesh[rf[i]].cnct[ii]].knotV[jj] = cp[nold + 5 + ii].knotV[jj];
			}
		}

		//add new elements
		double uanc[3] = {cp[tmesh[rf[i]].cnct[0]].knotU[2],
						  (cp[tmesh[rf[i]].cnct[0]].knotU[2] + cp[tmesh[rf[i]].cnct[1]].knotU[2]) / 2.,
						  cp[tmesh[rf[i]].cnct[1]].knotU[2]};
		double vanc[3] = {cp[tmesh[rf[i]].cnct[0]].knotV[2],
						  (cp[tmesh[rf[i]].cnct[0]].knotV[2] + cp[tmesh[rf[i]].cnct[1]].knotV[2]) / 2.,
						  cp[tmesh[rf[i]].cnct[1]].knotV[2]};
		double eanc[4][4] = {{uanc[0], vanc[0], uanc[1], vanc[1]}, {uanc[1], vanc[0], uanc[2], vanc[1]}, {uanc[1], vanc[1], uanc[2], vanc[2]}, {uanc[0], vanc[1], uanc[1], vanc[2]}};
		vector<Element> etmp(4);
		int cnctmp[4][4] = {{14, 36, 38, 37}, {36, 15, 39, 38}, {38, 39, 21, 40}, {37, 38, 40, 20}};
		for (j = 0; j < 4; j++)
		{
			etmp[j].act = 1;
			etmp[j].cnct[0] = cnctmp[j][0];
			etmp[j].cnct[1] = cnctmp[j][1];
			etmp[j].cnct[2] = cnctmp[j][2];
			etmp[j].cnct[3] = cnctmp[j][3];
			etmp[j].IEN = tmesh[rf[i]].IEN;
		}
		for (j = nold; j < cp.size(); j++)
		{
			for (k = 0; k < 4; k++)
			{
				etmp[k].IEN.push_back(j);
				//if(cp[j].knotU[2]==eanc[k][0] && cp[j].knotV[2]==eanc[k][1] && cp[j].knotU[3]==eanc[k][2] && cp[j].knotV[3]==eanc[k][3])
				//	etmp[k].cnct[0]=j;
				//if(cp[j].knotU[1]==eanc[k][0] && cp[j].knotV[2]==eanc[k][1] && cp[j].knotU[2]==eanc[k][2] && cp[j].knotV[3]==eanc[k][3])
				//	etmp[k].cnct[1]=j;
				//if(cp[j].knotU[1]==eanc[k][0] && cp[j].knotV[1]==eanc[k][1] && cp[j].knotU[2]==eanc[k][2] && cp[j].knotV[2]==eanc[k][3])
				//	etmp[k].cnct[2]=j;
				//if(cp[j].knotU[2]==eanc[k][0] && cp[j].knotV[1]==eanc[k][1] && cp[j].knotU[3]==eanc[k][2] && cp[j].knotV[2]==eanc[k][3])
				//	etmp[k].cnct[3]=j;
			}
		}
		for (j = 0; j < 4; j++)
			tmesh.push_back(etmp[j]);

		tmesh[rf[i]].act = 0;
	}
	//update elements
	for (i = 0; i < neold; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (j = nold; j < cp.size(); j++)
			{
				if (cp[j].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[j].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[j].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[j].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(j);
				}
			}
		}
	}
}

void TruncatedTspline::ElementSubdiv(int eid)
{
	int p(3); //polynomial order
	vector<int> pid = tmesh[eid].IEN;
	vector<vector<double>> knewU, knewV;
	double ktU = (cp[tmesh[eid].cnct[0]].knotU[2] + cp[tmesh[eid].cnct[1]].knotU[2]) / 2.;
	double ktV = (cp[tmesh[eid].cnct[0]].knotV[2] + cp[tmesh[eid].cnct[3]].knotV[2]) / 2.;
	double ktsU[5][5] = {{cp[tmesh[eid].cnct[0]].knotU[1], cp[tmesh[eid].cnct[0]].knotU[2], ktU, cp[tmesh[eid].cnct[0]].knotU[3], cp[tmesh[eid].cnct[0]].knotU[4]},
						 {cp[tmesh[eid].cnct[0]].knotU[0], cp[tmesh[eid].cnct[0]].knotU[1], cp[tmesh[eid].cnct[0]].knotU[2], ktU, cp[tmesh[eid].cnct[0]].knotU[3]},
						 {cp[tmesh[eid].cnct[0]].knotU[1], cp[tmesh[eid].cnct[0]].knotU[2], ktU, cp[tmesh[eid].cnct[0]].knotU[3], cp[tmesh[eid].cnct[0]].knotU[4]},
						 {cp[tmesh[eid].cnct[1]].knotU[1], ktU, cp[tmesh[eid].cnct[1]].knotU[2], cp[tmesh[eid].cnct[1]].knotU[3], cp[tmesh[eid].cnct[1]].knotU[4]},
						 {cp[tmesh[eid].cnct[0]].knotU[1], cp[tmesh[eid].cnct[0]].knotU[2], ktU, cp[tmesh[eid].cnct[0]].knotU[3], cp[tmesh[eid].cnct[0]].knotU[4]}};
	double ktsV[5][5] = {{cp[tmesh[eid].cnct[0]].knotV[0], cp[tmesh[eid].cnct[0]].knotV[1], cp[tmesh[eid].cnct[0]].knotV[2], ktV, cp[tmesh[eid].cnct[0]].knotV[3]},
						 {cp[tmesh[eid].cnct[0]].knotV[1], cp[tmesh[eid].cnct[0]].knotV[2], ktV, cp[tmesh[eid].cnct[0]].knotV[3], cp[tmesh[eid].cnct[0]].knotV[4]},
						 {cp[tmesh[eid].cnct[0]].knotV[1], cp[tmesh[eid].cnct[0]].knotV[2], ktV, cp[tmesh[eid].cnct[0]].knotV[3], cp[tmesh[eid].cnct[0]].knotV[4]},
						 {cp[tmesh[eid].cnct[0]].knotV[1], cp[tmesh[eid].cnct[0]].knotV[2], ktV, cp[tmesh[eid].cnct[0]].knotV[3], cp[tmesh[eid].cnct[0]].knotV[4]},
						 {cp[tmesh[eid].cnct[2]].knotV[1], ktV, cp[tmesh[eid].cnct[2]].knotV[2], cp[tmesh[eid].cnct[2]].knotV[3], cp[tmesh[eid].cnct[2]].knotV[4]}};

	vector<Vertex> pnew(5);
	for (int i = 0; i < 5; i++)
	{
		pnew[i].knotU[0] = ktsU[i][0];
		pnew[i].knotU[1] = ktsU[i][1];
		pnew[i].knotU[2] = ktsU[i][2];
		pnew[i].knotU[3] = ktsU[i][3];
		pnew[i].knotU[4] = ktsU[i][4];
		pnew[i].knotV[0] = ktsV[i][0];
		pnew[i].knotV[1] = ktsV[i][1];
		pnew[i].knotV[2] = ktsV[i][2];
		pnew[i].knotV[3] = ktsV[i][3];
		pnew[i].knotV[4] = ktsV[i][4];
	}
	for (int i = 0; i < 5; i++)
	{
		for (uint j = 0; j < pid.size(); j++)
		{
			if (pnew[i].knotU[0] >= cp[pid[j]].knotU[0] && pnew[i].knotU[4] <= cp[pid[j]].knotU[4] && pnew[i].knotV[0] >= cp[pid[j]].knotV[0] && pnew[i].knotV[4] <= cp[pid[j]].knotV[4])
			{
				//3 cases
				vector<double> ku1, kv1;
				vector<vector<double>> Tu, Tv;
				vector<double> ku(cp[pid[j]].knotU, cp[pid[j]].knotU + 5);
				vector<double> kv(cp[pid[j]].knotV, cp[pid[j]].knotV + 5);
				if (pnew[i].knotU[2] != ku[2] && pnew[i].knotV[2] == kv[2]) //case 1
				{
					vector<double> insU(1, ktU);
					InsertKnots(ku, insU, ku1);
					TMatrix(ku, ku1, p, Tu);
					kv1 = kv;
					Tv.assign(1, vector<double>(1, 1.));
				}
				else if (pnew[i].knotV[2] != kv[2] && pnew[i].knotU[2] == ku[2]) //case 2
				{
					vector<double> insV(1, ktV);
					InsertKnots(kv, insV, kv1);
					TMatrix(kv, kv1, p, Tv);
					ku1 = ku;
					Tu.assign(1, vector<double>(1, 1.));
				}
				else //case 3
				{
					vector<double> insU(1, ktU);
					vector<double> insV(1, ktV);
					InsertKnots(ku, insU, ku1);
					InsertKnots(kv, insV, kv1);
					TMatrix(ku, ku1, p, Tu);
					TMatrix(kv, kv1, p, Tv);
				}
				//double coef;
				for (uint iv = 0; iv < kv1.size() - 4; iv++)
				{
					for (uint iu = 0; iu < ku.size() - 4; iu++)
					{
					}
				}
				/*for(l=0;l<kv1.size()-4;l++)
				{
					for(m=0;m<ku1.size()-4;m++)
					{
						vector<double> tmpU(ku1.begin()+m,ku1.begin()+m+5),tmpV(kv1.begin()+l,kv1.begin()+l+5);
						if(tmpU==knewU[i] && tmpV==knewV[i])
						{
							double coef=Tu[m][0]*Tv[l][0];
							cp[pid[k]].tbf.push_back(nold+i);
							cp[pid[k]].tc.push_back(coef);
							cp[nold+i].coor[0]+=coef*cp[pid[k]].coor[0];
							cp[nold+i].coor[1]+=coef*cp[pid[k]].coor[1];
							cp[nold+i].coor[2]+=coef*cp[pid[k]].coor[2];
						}
					}
				}*/
			}
		}
	}

	//update all local knot vectors
}

void TruncatedTspline::Parmt2Phys(int eid, double u, double v, Vertex &pt)
{
	unsigned int i, j;
	int loc, loc1;
	double Ni;
	pt.coor[0] = 0.;
	pt.coor[1] = 0.;
	pt.coor[2] = 0.;
	vector<double> ku(5, 0.), kv(5, 0.), uval, vval;
	BSplineBasis bu, bv;
	for (i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		loc = tmesh[eid].IEN[i];
		ku.assign(cp[loc].knotU, cp[loc].knotU + 5);
		kv.assign(cp[loc].knotV, cp[loc].knotV + 5);
		bu.Set(3, ku);
		bv.Set(3, kv);
		bu.BasisFunction(0, u, 0, uval);
		bv.BasisFunction(0, v, 0, vval);
		Ni = uval[0] * vval[0];
		for (j = 0; j < cp[loc].tbf.size(); j++)
		{
			//cout<<"pid: "<<loc<<' '<<cp[loc].tbf.size()<<' '<<cp[loc].tc.size()<<'\n';
			//getchar();
			loc1 = cp[loc].tbf[j];
			ku.assign(cp[loc1].knotU, cp[loc1].knotU + 5);
			kv.assign(cp[loc1].knotV, cp[loc1].knotV + 5);
			bu.Set(3, ku);
			bv.Set(3, kv);
			bu.BasisFunction(0, u, 0, uval);
			bv.BasisFunction(0, v, 0, vval);
			Ni -= cp[loc].tc[j] * uval[0] * vval[0];
		}
		pt.coor[0] += Ni * cp[loc].coor[0];
		pt.coor[1] += Ni * cp[loc].coor[1];
		pt.coor[2] += Ni * cp[loc].coor[2];
	}
}

void TruncatedTspline::Parmt2Phys_v0(int eid, double u, double v, Vertex &pt)
{
	pt.coor[0] = 0.;
	pt.coor[1] = 0.;
	pt.coor[2] = 0.;
	vector<double> Ni, Nt;
	Ni.resize(tmesh[eid].IEN.size());
	Nt.resize(tmesh[eid].IEN.size());
	vector<double> ku(5, 0.), kv(5, 0.), uval, vval;
	BSplineBasis bu, bv;
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		int loc = tmesh[eid].IEN[i];
		ku.assign(cp[loc].knotU, cp[loc].knotU + 5);
		kv.assign(cp[loc].knotV, cp[loc].knotV + 5);
		bu.Set(3, ku);
		bv.Set(3, kv);
		bu.BasisFunction(0, u, 0, uval);
		bv.BasisFunction(0, v, 0, vval);
		Ni[i] = uval[0] * vval[0];
		Nt[i] = Ni[i];
		//for(j=0;j<cp[loc].tbf.size();j++)
		//{
		//	//cout<<"pid: "<<loc<<' '<<cp[loc].tbf.size()<<' '<<cp[loc].tc.size()<<'\n';
		//	//getchar();
		//	loc1=cp[loc].tbf[j];
		//	ku.assign(cp[loc1].knotU,cp[loc1].knotU+5);
		//	kv.assign(cp[loc1].knotV,cp[loc1].knotV+5);
		//	bu.Set(3,ku);
		//	bv.Set(3,kv);
		//	bu.BasisFunction(0,u,0,uval);
		//	bv.BasisFunction(0,v,0,vval);
		//	Ni-=cp[loc].tc[j]*uval[0]*vval[0];
		//}
		//pt.coor[0]+=Ni*cp[loc].coor[0];
		//pt.coor[1]+=Ni*cp[loc].coor[1];
		//pt.coor[2]+=Ni*cp[loc].coor[2];
	}
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		int loc = tmesh[eid].IEN[i];
		for (uint j = 0; j < cp[loc].tbf.size(); j++)
		{
			vector<int>::iterator it = find(tmesh[eid].IEN.begin(), tmesh[eid].IEN.end(), cp[loc].tbf[j]);
			if (it != tmesh[eid].IEN.end())
			{
				int pos(it - tmesh[eid].IEN.begin());
				Nt[i] -= Ni[pos] * cp[loc].tc[j];
			}
		}
		pt.coor[0] += Nt[i] * cp[loc].coor[0];
		pt.coor[1] += Nt[i] * cp[loc].coor[1];
		pt.coor[2] += Nt[i] * cp[loc].coor[2];
	}
}

double TruncatedTspline::PartitionOfUnity(int eid, double u, double v)
{
	unsigned int i, j;
	int loc, loc1;
	double Ni, sum(0.);
	vector<double> ku(5, 0.), kv(5, 0.), uval, vval;
	BSplineBasis bu, bv;
	//vector<int> bfn;
	//vector<double> bfsum;
	//double sum_old(0.);
	for (i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		loc = tmesh[eid].IEN[i];
		ku.assign(cp[loc].knotU, cp[loc].knotU + 5);
		kv.assign(cp[loc].knotV, cp[loc].knotV + 5);
		bu.Set(3, ku);
		bv.Set(3, kv);
		bu.BasisFunction(0, u, 0, uval);
		bv.BasisFunction(0, v, 0, vval);
		Ni = uval[0] * vval[0];
		//if (loc < npt_old) sum_old += Ni;
		for (j = 0; j < cp[loc].tbf.size(); j++)
		{
			loc1 = cp[loc].tbf[j];
			ku.assign(cp[loc1].knotU, cp[loc1].knotU + 5);
			kv.assign(cp[loc1].knotV, cp[loc1].knotV + 5);
			bu.Set(3, ku);
			bv.Set(3, kv);
			bu.BasisFunction(0, u, 0, uval);
			bv.BasisFunction(0, v, 0, vval);
			//if (eid >= nel_old)
			//{
			//	cout << eid<<" "<<loc<<" "<<loc1<<" "<<uval[0] * vval[0] << "\n";
			//	getchar();
			//}
			Ni -= cp[loc].tc[j] * uval[0] * vval[0];
			//vector<int>::iterator it = find(bfn.begin(),bfn.end(),loc1);
			//if (it == bfn.end())
			//{
			//	bfn.push_back(loc1);
			//	bfsum.push_back(cp[loc].tc[j]);
			//}
			//else
			//{
			//	int pos(it - bfn.begin());
			//	bfsum[pos] += cp[loc].tc[j];
			//}
		}
		sum += Ni;
	}
	//if (eid >= nel_old)
	//{
	//	cout << "sum of old: " << sum_old << "\n";
	//	cout << "sum of trun coef:\n";
	//	for (uint i = 0; i < bfn.size(); i++)
	//	{
	//		cout << bfn[i] << " " << bfsum[i] << "\n";
	//	}
	//	getchar();
	//}
	return sum;
}

void TruncatedTspline::VisualizeVTK(string fn)
{
	vector<Vertex> spt;
	vector<double> sval;
	vector<Element> sele;
	vector<Vertex> lpt;		   //visulize parameter lines
	vector<array<int, 2>> led; //line connectivity
	int ns(2), ecount(0), loc0, loc1, loc2;
	vector<double> su(ns), sv(ns);

	for (uint e = 0; e < tmesh.size(); e++)
	{
		if (tmesh[e].act == 1 && tmesh[e].type != 2 && tmesh[e].type != 3)
		{
			for (int i = 0; i < ns; i++)
			{
				loc0 = tmesh[e].cnct[0];
				loc1 = tmesh[e].cnct[1];
				loc2 = tmesh[e].cnct[3];
				su[i] = cp[loc0].knotU[2] + i * (cp[loc1].knotU[2] - cp[loc0].knotU[2]) / (ns - 1);
				sv[i] = cp[loc0].knotV[2] + i * (cp[loc2].knotV[2] - cp[loc0].knotV[2]) / (ns - 1);
			}

			int loc(0);
			for (int a = 0; a < ns; a++)
			{
				for (int b = 0; b < ns; b++)
				{
					Vertex pt;
					Parmt2Phys(e, su[b], sv[a], pt);
					spt.push_back(pt);
					if (a == 0 || a == ns - 1 || b == 0 || b == ns - 1)
					{
						lpt.push_back(pt);
					}
					double sum = PartitionOfUnity(e, su[b], sv[a]);
					sval.push_back(sum);
				}
			}
			for (int a = 0; a < ns - 1; a++)
			{
				for (int b = 0; b < ns - 1; b++)
				{
					Element el;
					el.cnct[0] = ecount * ns * ns + a * ns + b;
					el.cnct[1] = ecount * ns * ns + a * ns + b + 1;
					el.cnct[2] = ecount * ns * ns + (a + 1) * ns + b + 1;
					el.cnct[3] = ecount * ns * ns + (a + 1) * ns + b;
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
	}

	string fname = fn + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i < spt.size(); i++)
		{
			fout << spt[i].coor[0] << " " << spt[i].coor[1] << " " << spt[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i].cnct[0] << " " << sele[i].cnct[1] << " " << sele[i].cnct[2] << " " << sele[i].cnct[3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
		}
		fout << "\nPOINT_DATA " << sval.size() << "\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < sval.size(); i++)
		{
			fout << sval[i] << "\n";
		}
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
			fout1 << lpt[i].coor[0] << " " << lpt[i].coor[1] << " " << lpt[i].coor[2] << "\n";
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

void TruncatedTspline::PseudoProblem()
{
	double pts[36][3] = {{-0.1, -0.1, 0.}, {0., -0.1, 0.}, {1., -0.1, 0.}, {2., -0.1, 0.}, {3., -0.1, 0.}, {3.1, -0.1, 0.}, {-0.1, 0., 0.}, {0., 0., 0.}, {1., 0., 0.}, {2., 0., 0.}, {3., 0., 0.}, {3.1, 0., 0.}, {-0.1, 1., 0.}, {0., 1., 0.}, {1., 1., 0.}, {2., 1., 0.}, {3., 1., 0.}, {3.1, 1., 0.}, {-0.1, 2., 0.}, {0., 2., 0.}, {1., 2., 0.}, {2., 2., 0.}, {3., 2., 0.}, {3.1, 2., 0.}, {-0.1, 3., 0.}, {0., 3., 0.}, {1., 3., 0.}, {2., 3., 0.}, {3., 3., 0.}, {3.1, 3., 0.}, {-0.1, 3.1, 0.}, {0., 3.1, 0.}, {1., 3.1, 0.}, {2., 3.1, 0.}, {3., 3.1, 0.}, {3.1, 3.1, 0.}};
	int ele[9][4] = {{7, 8, 14, 13}, {8, 9, 15, 14}, {9, 10, 16, 15}, {13, 14, 20, 19}, {14, 15, 21, 20}, {15, 16, 22, 21}, {19, 20, 26, 25}, {20, 21, 27, 26}, {21, 22, 28, 27}};
	int ien[9][16] = {{0, 1, 2, 3, 6, 7, 8, 9, 12, 13, 14, 15, 18, 19, 20, 21}, {1, 2, 3, 4, 7, 8, 9, 10, 13, 14, 15, 16, 19, 20, 21, 22}, {2, 3, 4, 5, 8, 9, 10, 11, 14, 15, 16, 17, 20, 21, 22, 23}, {6, 7, 8, 9, 12, 13, 14, 15, 18, 19, 20, 21, 24, 25, 26, 27}, {7, 8, 9, 10, 13, 14, 15, 16, 19, 20, 21, 22, 25, 26, 27, 28}, {8, 9, 10, 11, 14, 15, 16, 17, 20, 21, 22, 23, 26, 27, 28, 29}, {12, 13, 14, 15, 18, 19, 20, 21, 24, 25, 26, 27, 30, 31, 32, 33}, {13, 14, 15, 16, 19, 20, 21, 22, 25, 26, 27, 28, 31, 32, 33, 34}, {14, 15, 16, 17, 20, 21, 22, 23, 26, 27, 28, 29, 32, 33, 34, 35}};
	double kU[10] = {0., 0., 0., 0., 1., 2., 3., 3., 3., 3.};
	double kV[10] = {0., 0., 0., 0., 1., 2., 3., 3., 3., 3.};

	cp.resize(36);
	tmesh.resize(9);
	int iu, iv;
	for (int i = 0; i < 36; i++)
	{
		cp[i].coor[0] = pts[i][0];
		cp[i].coor[1] = pts[i][1];
		cp[i].coor[2] = pts[i][2];
		iu = i % 6;
		iv = i / 6;
		for (int j = 0; j < 5; j++)
		{
			cp[i].knotU[j] = kU[iu + j];
			cp[i].knotV[j] = kU[iv + j];
		}
	}
	for (int i = 0; i < 9; i++)
	{
		tmesh[i].act = 1;
		;
		for (int j = 0; j < 4; j++)
			tmesh[i].cnct[j] = ele[i][j];
		tmesh[i].IEN.resize(16);
		for (int j = 0; j < 16; j++)
			tmesh[i].IEN[j] = ien[i][j];
	}
}

void TruncatedTspline::VisualizeControlMesh(string fn, int outflag)
{
	string fname(fn + "_CM_all.vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << cp.size() << " float\n";
		for (uint i = 0; i < cp.size(); i++)
		{
			fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
		}

		if (outflag == 1)
		{
			int nel_act(0);
			for (uint i = 0; i < tmesh.size(); i++)
			{
				if (tmesh[i].act == 1 /*&& tmesh[i].type != 2 && tmesh[i].type != 3*/)
					nel_act++;
			}
			fout << "\nCELLS " << nel_act << " " << 5 * nel_act << '\n';
			for (uint i = 0; i < tmesh.size(); i++)
			{
				if (tmesh[i].act == 1 /*&& tmesh[i].type != 2 && tmesh[i].type != 3*/)
				{
					fout << "4 " << tmesh[i].cnct[0] << ' ' << tmesh[i].cnct[1] << ' ' << tmesh[i].cnct[2] << ' ' << tmesh[i].cnct[3] << '\n';
				}
			}
			fout << "\nCELL_TYPES " << nel_act << '\n';
			for (uint i = 0; i < nel_act; i++)
			{
				fout << "9\n";
			}
		}
		else if (outflag == 0)
		{
			fout << "\nCELLS " << tmesh.size() << " " << 5 * tmesh.size() << '\n';
			for (uint i = 0; i < tmesh.size(); i++)
			{
				fout << "4 " << tmesh[i].cnct[0] << ' ' << tmesh[i].cnct[1] << ' ' << tmesh[i].cnct[2] << ' ' << tmesh[i].cnct[3] << '\n';
			}
			fout << "\nCELL_TYPES " << tmesh.size() << '\n';
			for (uint i = 0; i < tmesh.size(); i++)
			{
				fout << "9\n";
			}
			fout << "\nCELL_DATA " << tmesh.size() << '\n';
			fout << "SCALARS cell_scalars int 1\n";
			fout << "LOOKUP_TABLE default\n";
			for (uint i = 0; i < tmesh.size(); i++)
			{
				fout << tmesh[i].act << "\n";
			}
		}

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
		cout << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline::Test()
{
	//uint i,j;
	//for(i=0;i<cp.size();i++)
	//{
	//	cout<<"pid "<<i<<": ";
	//	for(j=0;j<cp[i].tbf.size();j++)
	//	{
	//		cout<<cp[i].tbf[j]<<" ";
	//	}
	//	cout<<'\n';
	//}
	//cout<<'\n';

	//for(i=0;i<tmesh.size();i++)
	//{
	//	cout<<"eid: "<<i<<" ";
	//	for(j=0;j<tmesh[i].IEN.size();j++)
	//	{
	//		cout<<tmesh[i].IEN[j]<<" ";
	//	}
	//	cout<<'\n';
	//}
	//cout<<'\n';
	//for(i=tmesh.size()-5;i<tmesh.size();i++)
	//{
	//	cout<<"eid "<<i<<": ";
	//	for(j=0;j<4;j++)
	//	{
	//		cout<<tmesh[i].cnct[j]<<" ";
	//	}
	//	cout<<'\n';
	//}
}

void TruncatedTspline::CreateInputMesh(string fn, int nu, int nv)
{
	double du(1.), dv(1.);
	vector<array<double, 3>> pts((nu + 1) * (nv + 1));
	vector<array<int, 4>> eles(nu * nv);
	int count(0);
	for (int i = 0; i < nv + 1; i++)
	{
		for (int j = 0; j < nu + 1; j++)
		{
			pts[count][0] = du * j;
			pts[count][1] = dv * i;
			pts[count][2] = du * j * (du * nu - du * j) / (du * nu);
			count++;
		}
	}
	count = 0;
	for (int i = 0; i < nv; i++)
	{
		for (int j = 0; j < nu; j++)
		{
			eles[count][0] = i * (nu + 1) + j;
			eles[count][1] = i * (nu + 1) + j + 1;
			eles[count][2] = (i + 1) * (nu + 1) + j + 1;
			eles[count][3] = (i + 1) * (nu + 1) + j;
			count++;
		}
	}

	string fname = fn + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		for (uint i = 0; i < pts.size(); i++)
		{
			fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
		}
		fout << "\nCELLS " << eles.size() << " " << 5 * eles.size() << '\n';
		for (uint i = 0; i < eles.size(); i++)
		{
			fout << "4 " << eles[i][0] << " " << eles[i][1] << " " << eles[i][2] << " " << eles[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << eles.size() << '\n';
		for (uint i = 0; i < eles.size(); i++)
		{
			fout << "9\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	//knot vectors
	int npu(nu + 1), npv(nv + 1);
	vector<double> kU(npu + 4);
	vector<double> kV(npv + 4);
	for (int i = 0; i < 3; i++)
	{
		kU[i] = 0.;
		kV[i] = 0.;
	}
	for (int i = 3; i < npu + 1; i++)
	{
		kU[i] = double(i - 3);
	}
	for (int i = npu + 1; i < npu + 4; i++)
	{
		kU[i] = double(npu - 3);
	}
	for (int i = 3; i < npv + 1; i++)
	{
		kV[i] = double(i - 3);
	}
	for (int i = npv + 1; i < npv + 4; i++)
	{
		kV[i] = double(npv - 3);
	}
	//IEN
	vector<array<int, 16>> ien;
	for (int j = 1; j < nv - 1; j++)
	{
		for (int i = 1; i < nu - 1; i++)
		{
			array<int, 16> tmp = {(j - 1) * npu + i - 1, (j - 1) * npu + i, (j - 1) * npu + i + 1, (j - 1) * npu + i + 2, j * npu + i - 1, j * npu + i, j * npu + i + 1, j * npu + i + 2,
								  (j + 1) * npu + i - 1, (j + 1) * npu + i, (j + 1) * npu + i + 1, (j + 1) * npu + i + 2, (j + 2) * npu + i - 1, (j + 2) * npu + i, (j + 2) * npu + i + 1, (j + 2) * npu + i + 2};
			ien.push_back(tmp);
		}
	}

	string fname1 = fn + "_appendix.txt";
	ofstream fout1;
	fout1.open(fname1.c_str());
	if (fout1.is_open())
	{
		fout1 << "LocalKnotVectors: " << (nu + 1) * (nv + 1) << '\n';
		count = 0;
		for (int j = 0; j < npv; j++)
		{
			for (int i = 0; i < npu; i++)
			{
				fout1 << count << ' ';
				fout1 << kU[i] << ' ' << kU[i + 1] << ' ' << kU[i + 2] << ' ' << kU[i + 3] << ' ' << kU[i + 4] << ' ';
				fout1 << kV[j] << ' ' << kV[j + 1] << ' ' << kV[j + 2] << ' ' << kV[j + 3] << ' ' << kV[j + 4] << '\n';
				count++;
			}
		}
		fout1 << "IEN: " << (nu - 2) * (nv - 2) << '\n';
		count = 0;
		for (int j = 1; j < nv - 1; j++)
		{
			for (int i = 1; i < nu - 1; i++)
			{
				fout1 << j * nu + i << ' ';
				for (int k = 0; k < 16; k++)
					fout1 << ien[count][k] << ' ';
				fout1 << '\n';
				count++;
			}
		}
		fout1.close();
	}
	else
	{
		cout << "Cannot open " << fname1 << "!\n";
	}
}

void TruncatedTspline::SetPseudoProblem(string fn)
{
	//read quad vtk
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++)
			getline(fin, stmp); //skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
			//cp[i].coor[0]/=10.; cp[i].coor[1]/=10.; cp[i].coor[2]=0.;
			//cp[i].coor[2]=0.;
			cp[i].act = 1;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	//read local knot vectors and IEN
	string fname1(fn + "_appendix.txt");
	int nlkv, nien;
	ifstream fin1;
	fin1.open(fname1);
	if (fin1.is_open())
	{
		fin1 >> stmp >> nlkv;
		for (int i = 0; i < nlkv; i++)
		{
			fin1 >> itmp;
			fin1 >> cp[itmp].knotU[0] >> cp[itmp].knotU[1] >> cp[itmp].knotU[2] >> cp[itmp].knotU[3] >> cp[itmp].knotU[4];
			fin1 >> cp[itmp].knotV[0] >> cp[itmp].knotV[1] >> cp[itmp].knotV[2] >> cp[itmp].knotV[3] >> cp[itmp].knotV[4];
			for (int j = 0; j < 4; j++)
			{
				cp[itmp].kitvU[j] = cp[itmp].knotU[j + 1] - cp[itmp].knotU[j];
				cp[itmp].kitvV[j] = cp[itmp].knotV[j + 1] - cp[itmp].knotV[j];
			}
			cp[itmp].pm[0] = cp[itmp].knotU[2];
			cp[itmp].pm[1] = cp[itmp].knotV[2];
		}
		fin1 >> stmp >> nien;
		for (int i = 0; i < nien; i++)
		{
			fin1 >> itmp;
			tmesh[itmp].act = 1;
			tmesh[itmp].IEN.resize(16);
			for (int j = 0; j < 16; j++)
				fin1 >> tmesh[itmp].IEN[j];
		}
		fin1.close();
	}
	else
	{
		cerr << "Cannot open " << fname1 << "!\n";
	}

	int nx(11), ny(11);
	double kntu[11] = {0., 0., 1., 2., 3., 4., 5., 6., 7., 8., 8.};
	double kntv[11] = {0., 0., 1., 2., 3., 4., 5., 6., 7., 8., 8.};
	int loc(0);
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			cp[loc].index[0] = i;
			cp[loc].index[1] = j;
			cp[loc].pm[0] = kntu[i];
			cp[loc].pm[1] = kntv[j];
			loc++;
		}
	}
	uanc.clear();
	vanc.clear();
	for (int i = 0; i < nx; i++)
	{
		pair<int, double> utmp(i, kntu[i]);
		uanc.push_back(utmp);
	}
	for (int i = 0; i < ny; i++)
	{
		pair<int, double> vtmp(i, kntv[i]);
		vanc.push_back(vtmp);
	}
}

void TruncatedTspline::ElementSubdivTest(int eid) //hierarchical level basis functions
{
	//int lev=tmesh[eid].lev;
	//if(lev+2>phid.size())
	//{
	//	vector<int> tmp;
	//	phid.push_back(tmp);
	//}
	int p(3), npold(cp.size()), neold(tmesh.size());
	uint i, j, k;
	vector<int> pid = tmesh[eid].IEN;
	int cid[4] = {tmesh[eid].cnct[0], tmesh[eid].cnct[1], tmesh[eid].cnct[2], tmesh[eid].cnct[3]};

	//define 5 new points
	//double ktsU[5][5]={{(cp[cid[0]].knotU[1]+cp[cid[0]].knotU[2])/2.,cp[cid[0]].knotU[2],(cp[cid[0]].knotU[2]+cp[cid[0]].knotU[3])/2.,cp[cid[0]].knotU[3],(cp[cid[0]].knotU[3]+cp[cid[0]].knotU[4])/2.},
	//	{cp[cid[0]].knotU[1],(cp[cid[0]].knotU[1]+cp[cid[0]].knotU[2])/2.,cp[cid[0]].knotU[2],(cp[cid[0]].knotU[2]+cp[cid[0]].knotU[3])/2.,cp[cid[0]].knotU[3]},
	//	{(cp[cid[0]].knotU[1]+cp[cid[0]].knotU[2])/2.,cp[cid[0]].knotU[2],(cp[cid[0]].knotU[2]+cp[cid[0]].knotU[3])/2.,cp[cid[0]].knotU[3],(cp[cid[0]].knotU[3]+cp[cid[0]].knotU[4])/2.},
	//	{cp[cid[1]].knotU[1],(cp[cid[1]].knotU[1]+cp[cid[1]].knotU[2])/2.,cp[cid[1]].knotU[2],(cp[cid[1]].knotU[2]+cp[cid[1]].knotU[3])/2.,cp[cid[1]].knotU[3]},
	//	{(cp[cid[0]].knotU[1]+cp[cid[0]].knotU[2])/2.,cp[cid[0]].knotU[2],(cp[cid[0]].knotU[2]+cp[cid[0]].knotU[3])/2.,cp[cid[0]].knotU[3],(cp[cid[0]].knotU[3]+cp[cid[0]].knotU[4])/2.}};
	//double ktsV[5][5]={{cp[cid[0]].knotV[1],(cp[cid[0]].knotV[1]+cp[cid[0]].knotV[2])/2.,cp[cid[0]].knotV[2],(cp[cid[0]].knotV[2]+cp[cid[0]].knotV[3])/2.,cp[cid[0]].knotV[3]},
	//	{(cp[cid[0]].knotV[1]+cp[cid[0]].knotV[2])/2.,cp[cid[0]].knotV[2],(cp[cid[0]].knotV[2]+cp[cid[0]].knotV[3])/2.,cp[cid[0]].knotV[3],(cp[cid[0]].knotV[3]+cp[cid[0]].knotV[4])/2.},
	//	{(cp[cid[0]].knotV[1]+cp[cid[0]].knotV[2])/2.,cp[cid[0]].knotV[2],(cp[cid[0]].knotV[2]+cp[cid[0]].knotV[3])/2.,cp[cid[0]].knotV[3],(cp[cid[0]].knotV[3]+cp[cid[0]].knotV[4])/2.},
	//	{(cp[cid[0]].knotV[1]+cp[cid[0]].knotV[2])/2.,cp[cid[0]].knotV[2],(cp[cid[0]].knotV[2]+cp[cid[0]].knotV[3])/2.,cp[cid[0]].knotV[3],(cp[cid[0]].knotV[3]+cp[cid[0]].knotV[4])/2.},
	//	{cp[cid[3]].knotV[1],(cp[cid[3]].knotV[1]+cp[cid[3]].knotV[2])/2.,cp[cid[3]].knotV[2],(cp[cid[3]].knotV[2]+cp[cid[3]].knotV[3])/2.,cp[cid[3]].knotV[3]}};
	//vector<Vertex> pnew(5);
	//for(i=0; i<5; i++)
	//{
	//	for(j=0; j<5; j++)
	//	{
	//		pnew[i].knotU[j]=ktsU[i][j];
	//		pnew[i].knotV[j]=ktsV[i][j];
	//	}
	//}
	double intvU((cp[cid[1]].knotU[2] - cp[cid[0]].knotU[2]) / 2.);
	double intvV((cp[cid[3]].knotV[2] - cp[cid[0]].knotV[2]) / 2.);
	double ktsU[7] = {cp[cid[0]].knotU[2] - 2. * intvU, cp[cid[0]].knotU[2] - intvU, cp[cid[0]].knotU[2], cp[cid[0]].knotU[2] + intvU, cp[cid[1]].knotU[2], cp[cid[1]].knotU[2] + intvU, cp[cid[1]].knotU[2] + 2. * intvU};
	double ktsV[7] = {cp[cid[0]].knotV[2] - 2. * intvV, cp[cid[0]].knotV[2] - intvV, cp[cid[0]].knotV[2], cp[cid[0]].knotV[2] + intvV, cp[cid[3]].knotV[2], cp[cid[3]].knotV[2] + intvV, cp[cid[3]].knotV[2] + 2. * intvV};
	if (cp[cid[0]].knotU[1] == cp[cid[0]].knotU[2])
	{
		ktsU[1] = cp[cid[0]].knotU[1];
		if (cp[cid[0]].knotU[0] != cp[cid[0]].knotU[1])
		{
			ktsU[0] = cp[cid[0]].knotU[1] - intvU;
		}
		else
		{
			ktsU[0] = cp[cid[0]].knotU[1];
		}
	}
	if (cp[cid[1]].knotU[2] == cp[cid[1]].knotU[3])
	{
		ktsU[5] = cp[cid[1]].knotU[3];
		if (cp[cid[1]].knotU[3] != cp[cid[1]].knotU[4])
		{
			ktsU[6] = cp[cid[1]].knotU[3] + intvU;
		}
		else
		{
			ktsU[6] = cp[cid[1]].knotU[3];
		}
	}
	if (cp[cid[0]].knotV[1] == cp[cid[0]].knotV[2])
	{
		ktsV[1] = cp[cid[0]].knotV[1];
		if (cp[cid[0]].knotV[0] != cp[cid[0]].knotV[1])
		{
			ktsV[0] = cp[cid[0]].knotV[1] - intvV;
		}
		else
		{
			ktsV[0] = cp[cid[0]].knotV[1];
		}
	}
	if (cp[cid[3]].knotV[2] == cp[cid[3]].knotV[3])
	{
		ktsV[5] = cp[cid[3]].knotV[3];
		if (cp[cid[3]].knotV[3] != cp[cid[3]].knotV[4])
		{
			ktsV[6] = cp[cid[3]].knotV[3] + intvV;
		}
		else
		{
			ktsV[6] = cp[cid[3]].knotV[3];
		}
	}
	int setU[3][3] = {{-1, 0, -1}, {1, 2, 3}, {-1, 4, -1}};
	int setV[3][3] = {{-1, 1, -1}, {0, 2, 4}, {-1, 3, -1}};
	vector<Vertex> pnew(5);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (setU[i][j] != -1)
			{
				for (k = 0; k < 5; k++)
					pnew[setU[i][j]].knotU[k] = ktsU[j + k];
				for (k = 0; k < 4; k++)
					pnew[setU[i][j]].kitvU[k] = ktsU[j + k + 1] - ktsU[j + k];
			}
			if (setV[i][j] != -1)
			{
				for (k = 0; k < 5; k++)
					pnew[setV[i][j]].knotV[k] = ktsV[j + k];
				for (k = 0; k < 4; k++)
					pnew[setV[i][j]].kitvV[k] = ktsV[j + k + 1] - ktsV[j + k];
			}
		}
	}
	//for(i=0;i<5;i++)
	//{
	//	//cout<<pnew[i].knotU[0]<<' '<<pnew[i].knotU[1]<<' '<<pnew[i].knotU[2]<<' '<<pnew[i].knotU[3]<<' '<<pnew[i].knotU[4]<<'\n';
	//	//cout<<is_sorted(pnew[i].knotU,pnew[i].knotU+5)<<'\n';
	//	cout<<pnew[i].knotV[0]<<' '<<pnew[i].knotV[1]<<' '<<pnew[i].knotV[2]<<' '<<pnew[i].knotV[3]<<' '<<pnew[i].knotV[4]<<'\n';
	//	cout<<is_sorted(pnew[i].knotV,pnew[i].knotV+5)<<'\n';
	//}
	//getchar();

	vector<int> flag(5, 1);
	vector<int> pgn(5, -1);
	int count(0);
	for (i = 0; i < 5; i++)
	{
		for (j = 0; j < pid.size(); j++)
		{
			if (equal(pnew[i].knotU, pnew[i].knotU + 5, cp[pid[j]].knotU) && equal(pnew[i].knotV, pnew[i].knotV + 5, cp[pid[j]].knotV))
			{
				flag[i] = 0;
				pgn[i] = pid[j];
			}
		}
		if (pgn[i] == -1)
		{
			pgn[i] = npold + count;
			count++;
		}
	}
	vector<vector<double>> coef(5, vector<double>(pid.size(), 0.));
	for (i = 0; i < 5; i++)
	{
		if (flag[i] == 1)
		{
			for (j = 0; j < pid.size(); j++)
			{
				if (pnew[i].knotU[0] >= cp[pid[j]].knotU[0] && pnew[i].knotU[4] <= cp[pid[j]].knotU[4] && pnew[i].knotV[0] >= cp[pid[j]].knotV[0] && pnew[i].knotV[4] <= cp[pid[j]].knotV[4])
				{
					vector<double> ku1(10), kv1(10);
					vector<vector<double>> Tu, Tv;
					vector<double> ku(cp[pid[j]].knotU, cp[pid[j]].knotU + 5);
					vector<double> kv(cp[pid[j]].knotV, cp[pid[j]].knotV + 5);
					vector<double>::iterator it1, it2;
					it1 = set_union(cp[pid[j]].knotU, cp[pid[j]].knotU + 5, pnew[i].knotU, pnew[i].knotU + 5, ku1.begin());
					it2 = set_union(cp[pid[j]].knotV, cp[pid[j]].knotV + 5, pnew[i].knotV, pnew[i].knotV + 5, kv1.begin());
					ku1.resize(it1 - ku1.begin());
					kv1.resize(it2 - kv1.begin());
					TMatrix(ku, ku1, p, Tu);
					TMatrix(kv, kv1, p, Tv);
					it1 = search(ku1.begin(), ku1.end(), pnew[i].knotU, pnew[i].knotU + 5);
					it2 = search(kv1.begin(), kv1.end(), pnew[i].knotV, pnew[i].knotV + 5);
					int loc1 = it1 - ku1.begin();
					int loc2 = it2 - kv1.begin();
					coef[i][j] = Tu[loc1][0] * Tv[loc2][0];
					//double tmp=Tu[loc1][0]*Tv[loc2][0];
					//double coef=tmp;
					//for(k=0;k<cp[pid[j]].tbf.size();k++)
					//{
					//	for(uint a=0;a<pid.size();a++)
					//	{
					//		vector<int>::iterator it=find(cp[pid[j]].tbf.begin(),cp[pid[j]].tbf.end(),pid[a]);
					//		if(it!=cp[pid[j]].tbf.end())
					//		{
					//			coef-=tmp*cp[pid[j]].tc[it-cp[pid[j]].tbf.begin()];
					//		}
					//	}
					//}
					//pnew[i].coor[0]+=coef*cp[pid[j]].coor[0];
					//pnew[i].coor[1]+=coef*cp[pid[j]].coor[1];
					//pnew[i].coor[2]+=coef*cp[pid[j]].coor[2];
					//cp[pid[j]].tc.push_back(coef);
					//cp[pid[j]].tbf.push_back(pgn[i]);
				}
			}
		}
	}
	for (i = 0; i < 5; i++)
	{
		if (flag[i] == 1)
		{
			for (j = 0; j < pid.size(); j++)
			{
				if (pnew[i].knotU[0] >= cp[pid[j]].knotU[0] && pnew[i].knotU[4] <= cp[pid[j]].knotU[4] && pnew[i].knotV[0] >= cp[pid[j]].knotV[0] && pnew[i].knotV[4] <= cp[pid[j]].knotV[4])
				{
					double tmp = coef[i][j];
					for (k = 0; k < cp[pid[j]].tbf.size(); k++)
					{
						vector<int>::iterator it = find(pid.begin(), pid.end(), cp[pid[j]].tbf[k]);
						if (it != pid.end())
							tmp -= coef[i][it - pid.begin()] * cp[pid[j]].tc[k];
						//for(uint a=0;a<pid.size();a++)
						//{
						//	vector<int>::iterator it=find(cp[pid[j]].tbf.begin(),cp[pid[j]].tbf.end(),pid[a]);
						//	if(it!=cp[pid[j]].tbf.end())
						//	{
						//
						//	}
						//}
					}
					pnew[i].coor[0] += tmp * cp[pid[j]].coor[0];
					pnew[i].coor[1] += tmp * cp[pid[j]].coor[1];
					pnew[i].coor[2] += tmp * cp[pid[j]].coor[2];
					cp[pid[j]].tc.push_back(tmp);
				}
			}
		}
	}
	for (i = 0; i < 5; i++)
	{
		if (pgn[i] >= npold)
			cp.push_back(pnew[i]);
		if (flag[i] == 1)
			for (j = 0; j < pid.size(); j++)
				if (pnew[i].knotU[0] >= cp[pid[j]].knotU[0] && pnew[i].knotU[4] <= cp[pid[j]].knotU[4] && pnew[i].knotV[0] >= cp[pid[j]].knotV[0] && pnew[i].knotV[4] <= cp[pid[j]].knotV[4])
				{
					/*cp[pid[j]].tc.push_back(coef[i][j]);*/
					cp[pid[j]].tbf.push_back(pgn[i]);
				}
	}

	//new elements
	vector<Element> enew(4);
	int elen[4][4] = {{cid[0], pgn[0], pgn[2], pgn[1]}, {pgn[0], cid[1], pgn[3], pgn[2]}, {pgn[1], pgn[2], pgn[4], cid[3]}, {pgn[2], pgn[3], cid[2], pgn[4]}};
	for (i = 0; i < 4; i++)
	{
		enew[i].act = 1;
		enew[i].IEN = pid;
		for (j = 0; j < 4; j++)
		{
			enew[i].cnct[j] = elen[i][j];
		}
		for (j = 0; j < 5; j++)
		{
			if (flag[j] == 1)
				enew[i].IEN.push_back(pgn[j]);
		}
		tmesh.push_back(enew[i]);
	}
	tmesh[eid].act = 0;

	//update other elements
	for (i = 0; i < neold; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (j = 0; j < 5; j++)
			{
				if (flag[j] == 1 && pnew[j].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && pnew[j].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					pnew[j].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && pnew[j].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(pgn[j]);
				}
			}
		}
	}
}

void TruncatedTspline::CollectActives()
{
	paid.clear();
	paid.resize(cp.size());
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			paid[i] = count;
			count++;
		}
		else
		{
			paid[i] = -1;
		}
	}
	npta = count;

	eaid.clear();
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].type != 2 && tmesh[i].type != 3)
		{
			eaid.push_back(i);
		}
	}
}

void TruncatedTspline::BezierExtract(vector<BezierElement> &bzmesh)
{
	bzmesh.clear();
	//for(uint i=0; i<cp.size(); i++)
	//{
	//	cp[i].coor[0]*=cp[i].w;
	//	cp[i].coor[1]*=cp[i].w;
	//	cp[i].coor[2]*=cp[i].w;
	//}
	//for(uint i=69; i<cp.size(); i++)
	//{
	//	cp[i].w=1.;
	//}
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		double plen[2] = {cp[tmesh[eid].cnct[1]].knotU[2] - cp[tmesh[eid].cnct[0]].knotU[2], cp[tmesh[eid].cnct[3]].knotV[2] - cp[tmesh[eid].cnct[0]].knotV[2]};
		if (tmesh[eid].act == 1 && plen[0] != 0. && plen[1] != 0.)
		{
			//cout<<"eid: "<<eid<<'\n';
			BezierElementExtract(eid, bzmesh);
		}
	}
	//for(uint i=0; i<cp.size(); i++)
	//{
	//	cp[i].coor[0]/=cp[i].w;
	//	cp[i].coor[1]/=cp[i].w;
	//	cp[i].coor[2]/=cp[i].w;
	//}
}

void TruncatedTspline::BezierElementExtract(int eid, vector<BezierElement> &bzmesh)
{
	vector<array<double, 4>> be;
	BezierFinder_1(eid, be);
	uint nbzold = bzmesh.size();
	for (uint i = 0; i < be.size(); i++)
	{
		BezierUnit(tmesh[eid].IEN, be[i], bzmesh);
		//BezierUnit_glb(tmesh[eid].IEN,be[i],bzmesh);
	}
	for (uint i = nbzold; i < bzmesh.size(); i++)
	{
		bzmesh[i].prt = eid;
	}
	//int cid[4]={tmesh[eid].cnct[0],tmesh[eid].cnct[1],tmesh[eid].cnct[2],tmesh[eid].cnct[3]};
	//double intv[2]={cp[cid[1]].knotU[2]-cp[cid[0]].knotU[2],cp[cid[3]].knotV[2]-cp[cid[0]].knotV[2]};
	//vector<int> pid=tmesh[eid].IEN;
	//vector<vector<double>> coef(pid.size(),vector<double>(16));
	//vector<int> bid;
	//for(uint i=0;i<pid.size();i++)
	//{
	//	for(int j=0;j<4;j++)
	//	{
	//		if((cp[pid[i]].kitvU[j]!=0. && cp[pid[i]].kitvU[j]<intv[0]) || (cp[pid[i]].kitvV[j]!=0. && cp[pid[i]].kitvV[j]<intv[1]))
	//		{
	//			bid.push_back(i);
	//			break;
	//		}
	//	}
	//}
	//if(bid.size()==0)//no need to subdivide the element
	//{
	//	//BezierElement bzel;
	//	//bzel.cmat.resize(pid.size(),vector<double>(16,0.));
	//	//double bzku[6]={cp[cid[0]].knotU[2],cp[cid[0]].knotU[2],cp[cid[0]].knotU[2],cp[cid[1]].knotU[2],cp[cid[1]].knotU[2],cp[cid[1]].knotU[2]};
	//	//double bzkv[6]={cp[cid[0]].knotV[2],cp[cid[0]].knotV[2],cp[cid[0]].knotV[2],cp[cid[3]].knotV[2],cp[cid[3]].knotV[2],cp[cid[3]].knotV[2]};
	//	array<double,2> ktsU={cp[cid[0]].knotU[2],cp[cid[1]].knotU[2]};
	//	array<double,2> ktsV={cp[cid[0]].knotV[2],cp[cid[3]].knotV[2]};
	//	BezierUnit(pid,ktsU,ktsV);
	//	/*for(uint i=0;i<pid.size();i++)
	//	{
	//		if(cp[pid[i]].knotU[0]<ktsU[1] && cp[pid[i]].knotU[4]>ktsU[0] && cp[pid[i]].knotV[0]<ktsV[1] && cp[pid[i]].knotV[4]>ktsV[0])
	//		{
	//			vector<double> ku(cp[pid[i]].knotU,cp[pid[i]].knotU+5);
	//			vector<double> kv(cp[pid[i]].knotV,cp[pid[i]].knotV+5);
	//			vector<double> ku1,kv1;
	//			vector<vector<double>> Tu,Tv;
	//			BezierInsertKnots(ku,ktsU,ku1);
	//			BezierInsertKnots(kv,ktsV,kv1);
	//			TMatrix(ku,ku1,3,Tu);
	//			TMatrix(kv,kv1,3,Tv);
	//			vector<double>::iterator it1=search(ku1.begin(),ku1.end(),bzku,bzku+6);
	//			vector<double>::iterator it2=search(kv1.begin(),kv1.end(),bzkv,bzkv+6);
	//			int loc1(it1-ku1.begin()-1),loc2(it2-kv1.begin()-1),count(0);
	//			for(int j=0;j<4;j++)
	//			{
	//				for(int k=0;k<4;k++)
	//				{
	//					coef[i][count]=Tu[loc1+k][0]*Tv[loc2+j][0];
	//					count++;
	//				}
	//			}
	//		}
	//	}
	//	for(uint i=0;i<pid.size();i++)
	//	{
	//		if(cp[pid[i]].knotU[0]<ktsU[1] && cp[pid[i]].knotU[4]>ktsU[0] && cp[pid[i]].knotV[0]<ktsV[1] && cp[pid[i]].knotV[4]>ktsV[0])
	//		{
	//			for(int j=0;j<16;j++)
	//			{
	//				double tmp(coef[i][j]);
	//				for(uint k=0;k<cp[pid[i]].tbf.size();k++)
	//				{
	//					vector<int>::iterator it=find(pid.begin(),pid.end(),cp[pid[i]].tbf[k]);
	//					if(it!=pid.end())
	//					{
	//						int loc=it-pid.begin();
	//						tmp-=cp[pid[i]].tc[k]*coef[loc][j];
	//					}
	//				}
	//				bzel.pts[j][0]+=tmp*cp[pid[i]].coor[0];
	//				bzel.pts[j][1]+=tmp*cp[pid[i]].coor[1];
	//				bzel.pts[j][2]+=tmp*cp[pid[i]].coor[2];
	//				bzel.cmat[i][j]=tmp;
	//			}
	//		}
	//	}
	//	bzmesh.push_back(bzel);*/
	//}
	//else
	//{

	//}
}

void TruncatedTspline::BezierUnit(vector<int> &pid, array<double, 4> &kts, vector<BezierElement> &bzmesh)
{
	BezierElement bzel;
	bzel.cmat.resize(pid.size());
	for (uint i = 0; i < pid.size(); i++)
	{
		for (int j = 0; j < 16; j++)
		{
			bzel.cmat[i][j] = 0.;
		}
	}
	vector<vector<double>> coef(pid.size(), vector<double>(16));
	array<double, 2> ktsU = {kts[0], kts[1]};
	array<double, 2> ktsV = {kts[2], kts[3]};
	double bzku[6] = {ktsU[0], ktsU[0], ktsU[0], ktsU[1], ktsU[1], ktsU[1]};
	double bzkv[6] = {ktsV[0], ktsV[0], ktsV[0], ktsV[1], ktsV[1], ktsV[1]};
	for (uint i = 0; i < pid.size(); i++)
	{
		if (cp[pid[i]].knotU[0] < ktsU[1] && cp[pid[i]].knotU[4] > ktsU[0] && cp[pid[i]].knotV[0] < ktsV[1] && cp[pid[i]].knotV[4] > ktsV[0])
		{
			vector<double> ku(cp[pid[i]].knotU, cp[pid[i]].knotU + 5);
			vector<double> kv(cp[pid[i]].knotV, cp[pid[i]].knotV + 5);
			vector<double> ku1, kv1;
			vector<vector<double>> Tu, Tv;
			BezierInsertKnots(ku, ktsU, ku1);
			BezierInsertKnots(kv, ktsV, kv1);
			TMatrix(ku, ku1, 3, Tu);
			TMatrix(kv, kv1, 3, Tv);
			vector<double>::iterator it1 = search(ku1.begin(), ku1.end(), bzku, bzku + 6);
			vector<double>::iterator it2 = search(kv1.begin(), kv1.end(), bzkv, bzkv + 6);
			int loc1(it1 - ku1.begin() - 1), loc2(it2 - kv1.begin() - 1), count(0);
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					coef[i][count] = Tu[loc1 + k][0] * Tv[loc2 + j][0];
					count++;
				}
			}
		}
	}
	for (uint i = 0; i < pid.size(); i++)
	{
		if (cp[pid[i]].knotU[0] < ktsU[1] && cp[pid[i]].knotU[4] > ktsU[0] && cp[pid[i]].knotV[0] < ktsV[1] && cp[pid[i]].knotV[4] > ktsV[0])
		{
			for (int j = 0; j < 16; j++)
			{
				double tmp(coef[i][j]);
				for (uint k = 0; k < cp[pid[i]].tbf.size(); k++)
				{
					vector<int>::iterator it = find(pid.begin(), pid.end(), cp[pid[i]].tbf[k]);
					if (it != pid.end())
					{
						int loc = it - pid.begin();
						tmp -= cp[pid[i]].tc[k] * coef[loc][j];
					}
				}
				bzel.pts[j][0] += tmp * cp[pid[i]].coor[0];
				bzel.pts[j][1] += tmp * cp[pid[i]].coor[1];
				bzel.pts[j][2] += tmp * cp[pid[i]].coor[2];
				bzel.cmat[i][j] = tmp;
			}
		}
	}
	for (uint i = 0; i < pid.size(); i++)
	{
		if (paid[pid[i]] != -1)
			bzel.IEN.push_back(paid[pid[i]]); //problematic
	}
	bzmesh.push_back(bzel);
}

void TruncatedTspline::BezierUnit_glb(vector<int> &pid, array<double, 4> &kts, vector<BezierElement> &bzmesh)
{
	BezierElement bzel;
	bzel.cmat.resize(pid.size());
	for (uint i = 0; i < pid.size(); i++)
	{
		for (int j = 0; j < 16; j++)
		{
			bzel.cmat[i][j] = 0.;
		}
	}
	//vector<vector<double>> coef(pid.size(),vector<double>(16));
	array<double, 2> ktsU = {kts[0], kts[1]};
	array<double, 2> ktsV = {kts[2], kts[3]};
	double bzku[6] = {ktsU[0], ktsU[0], ktsU[0], ktsU[1], ktsU[1], ktsU[1]};
	double bzkv[6] = {ktsV[0], ktsV[0], ktsV[0], ktsV[1], ktsV[1], ktsV[1]};
	for (uint i = 0; i < pid.size(); i++)
	{
		if (cp[pid[i]].knotU[0] < ktsU[1] && cp[pid[i]].knotU[4] > ktsU[0] && cp[pid[i]].knotV[0] < ktsV[1] && cp[pid[i]].knotV[4] > ktsV[0])
		{
			vector<double> ku(cp[pid[i]].knotU, cp[pid[i]].knotU + 5);
			vector<double> kv(cp[pid[i]].knotV, cp[pid[i]].knotV + 5);
			vector<double> ku1, kv1;
			vector<vector<double>> Tu, Tv;
			BezierInsertKnots(ku, ktsU, ku1);
			BezierInsertKnots(kv, ktsV, kv1);
			TMatrix(ku, ku1, 3, Tu);
			TMatrix(kv, kv1, 3, Tv);
			vector<double>::iterator it1 = search(ku1.begin(), ku1.end(), bzku, bzku + 6);
			vector<double>::iterator it2 = search(kv1.begin(), kv1.end(), bzkv, bzkv + 6);
			int loc1(it1 - ku1.begin() - 1), loc2(it2 - kv1.begin() - 1), count(0);
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					bzel.cmat[i][count] = Tu[loc1 + k][0] * Tv[loc2 + j][0];
					bzel.pts[count][0] += bzel.cmat[i][count] * cp[pid[i]].coor[0] * cp[pid[i]].w;
					bzel.pts[count][1] += bzel.cmat[i][count] * cp[pid[i]].coor[1] * cp[pid[i]].w;
					bzel.pts[count][2] += bzel.cmat[i][count] * cp[pid[i]].coor[2] * cp[pid[i]].w;
					bzel.w[count] += bzel.cmat[i][count] * cp[pid[i]].w;
					//bzel.pts[count][0]+=bzel.cmat[i][count]*cp[pid[i]].coor[0]*1.;
					//bzel.pts[count][1]+=bzel.cmat[i][count]*cp[pid[i]].coor[1]*1.;
					//bzel.pts[count][2]+=bzel.cmat[i][count]*cp[pid[i]].coor[2]*1.;
					//bzel.w[count]+=bzel.cmat[i][count]*1.;
					count++;
				}
			}
		}
	}
	for (int i = 0; i < 16; i++)
	{
		bzel.pts[i][0] /= bzel.w[i];
		bzel.pts[i][1] /= bzel.w[i];
		bzel.pts[i][2] /= bzel.w[i];
	}
	//if(bzmesh.size()<68)
	//{
	//	double cir_x[2]={0.,0.};
	//	double u(.5), sum_w(0.);
	//	double Nu[4]={(1.-u)*(1.-u)*(1.-u),3.*(1.-u)*(1.-u)*u,3.*(1.-u)*u*u,u*u*u};
	//	for(int i=0; i<4; i++)
	//	{
	//		//cir_x[0]+=bzel.pts[i][0]*bzel.w[i]*Nu[i];
	//		//cir_x[1]+=bzel.pts[i][1]*bzel.w[i]*Nu[i];
	//		//sum_w+=bzel.w[i]*Nu[i];
	//		cir_x[0]+=bzel.pts[i][0]*Nu[i];
	//		cir_x[1]+=bzel.pts[i][1]*Nu[i];
	//	}
	//	//cir_x[0]/=sum_w;
	//	//cir_x[1]/=sum_w;
	//	double err=sqrt(cir_x[0]*cir_x[0]+cir_x[1]*cir_x[1])-1.;
	//	cout<<err<<"\n";
	//}

	//for(uint i=0;i<pid.size();i++)
	//{
	//	if(cp[pid[i]].knotU[0]<ktsU[1] && cp[pid[i]].knotU[4]>ktsU[0] && cp[pid[i]].knotV[0]<ktsV[1] && cp[pid[i]].knotV[4]>ktsV[0])
	//	{
	//		for(int j=0;j<16;j++)
	//		{
	//			double tmp(coef[i][j]);
	//			for(uint k=0;k<cp[pid[i]].tbf.size();k++)
	//			{
	//				vector<int>::iterator it=find(pid.begin(),pid.end(),cp[pid[i]].tbf[k]);
	//				if(it!=pid.end())
	//				{
	//					int loc=it-pid.begin();
	//					tmp-=cp[pid[i]].tc[k]*coef[loc][j];
	//				}
	//			}
	//			bzel.pts[j][0]+=tmp*cp[pid[i]].coor[0];
	//			bzel.pts[j][1]+=tmp*cp[pid[i]].coor[1];
	//			bzel.pts[j][2]+=tmp*cp[pid[i]].coor[2];
	//			bzel.cmat[i][j]=tmp;
	//		}
	//	}
	//}
	for (uint i = 0; i < pid.size(); i++)
	{
		if (paid[pid[i]] != -1)
		{
			bzel.IEN.push_back(paid[pid[i]]); //problematic
			bzel.w_nurbs.push_back(cp[pid[i]].w);
			//cout<<cp[pid[i]].w<<" ";
		}
	}
	//cout<<"\n";
	//getchar();
	bzmesh.push_back(bzel);
}

void TruncatedTspline::BezierVTK(string fn, vector<BezierElement> &bzmesh)
{
	vector<Vertex> spt;
	vector<double> sval;
	vector<Element> sele;
	vector<Vertex> lpt;		   //visulize parameter lines
	vector<array<int, 2>> led; //line connectivity
	int ns(5), ecount(0);
	vector<double> su(ns), sv(ns);
	for (int i = 0; i < ns; i++)
	{
		su[i] = double(i) / double(ns - 1);
		sv[i] = double(i) / double(ns - 1);
	}

	for (uint e = 0; e < bzmesh.size(); e++)
	{
		int loc(0);
		for (int a = 0; a < ns; a++)
		{
			for (int b = 0; b < ns; b++)
			{
				Vertex pt;
				double ptmp[3];
				bzmesh[e].Para2Phys(su[b], sv[a], ptmp);
				pt.coor[0] = ptmp[0];
				pt.coor[1] = ptmp[1];
				pt.coor[2] = ptmp[2];
				spt.push_back(pt);
				if (a == 0 || a == ns - 1 || b == 0 || b == ns - 1)
				{
					lpt.push_back(pt);
				}
				//double sum=PartitionOfUnity(e,su[b],sv[a]);
				//sval.push_back(sum);
			}
		}
		for (int a = 0; a < ns - 1; a++)
		{
			for (int b = 0; b < ns - 1; b++)
			{
				Element el;
				el.cnct[0] = ecount * ns * ns + a * ns + b;
				el.cnct[1] = ecount * ns * ns + a * ns + b + 1;
				el.cnct[2] = ecount * ns * ns + (a + 1) * ns + b + 1;
				el.cnct[3] = ecount * ns * ns + (a + 1) * ns + b;
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

	string fname = fn + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i < spt.size(); i++)
		{
			fout << spt[i].coor[0] << " " << spt[i].coor[1] << " " << spt[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i].cnct[0] << " " << sele[i].cnct[1] << " " << sele[i].cnct[2] << " " << sele[i].cnct[3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
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
			fout1 << lpt[i].coor[0] << " " << lpt[i].coor[1] << " " << lpt[i].coor[2] << "\n";
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

void TruncatedTspline::BezierFinder(int eid, vector<array<double, 4>> &be)
{
	int cid[4] = {tmesh[eid].cnct[0], tmesh[eid].cnct[1], tmesh[eid].cnct[2], tmesh[eid].cnct[3]};
	vector<int> pid = tmesh[eid].IEN;
	array<double, 4> e0 = {cp[cid[0]].knotU[2], cp[cid[1]].knotU[2], cp[cid[0]].knotV[2], cp[cid[3]].knotV[2]};
	be.clear();
	be.push_back(e0);
	int flag(1);
	while (flag)
	{
		flag = 0;
		vector<array<double, 4>> betmp;
		for (uint i = 0; i < be.size(); i++)
		{
			int eflag(0);
			for (uint j = 0; j < pid.size(); j++)
			{
				if (cp[pid[j]].knotU[0] < be[i][1] && cp[pid[j]].knotU[4] > be[i][0] && cp[pid[j]].knotV[0] < be[i][3] && cp[pid[j]].knotV[4] > be[i][2])
				{
					for (int k = 0; k < 4; k++)
					{
						if ((cp[pid[j]].kitvU[k] != 0. && cp[pid[j]].kitvU[k] < be[i][1] - be[i][0]) || (cp[pid[j]].kitvV[k] != 0. && cp[pid[j]].kitvV[k] < be[i][3] - be[i][2]))
						{
							flag = 1;
							eflag = 1;
						}
					}
				}
			}
			if (eflag)
			{
				array<double, 4> tmp1 = {be[i][0], (be[i][0] + be[i][1]) / 2., be[i][2], (be[i][2] + be[i][3]) / 2.};
				array<double, 4> tmp2 = {(be[i][0] + be[i][1]) / 2., be[i][1], be[i][2], (be[i][2] + be[i][3]) / 2.};
				array<double, 4> tmp3 = {be[i][0], (be[i][0] + be[i][1]) / 2., (be[i][2] + be[i][3]) / 2., be[i][3]};
				array<double, 4> tmp4 = {(be[i][0] + be[i][1]) / 2., be[i][1], (be[i][2] + be[i][3]) / 2., be[i][3]};
				betmp.push_back(tmp1);
				betmp.push_back(tmp2);
				betmp.push_back(tmp3);
				betmp.push_back(tmp4);
			}
			else
			{
				betmp.push_back(be[i]);
			}
		}
		be.clear();
		be = betmp;
	}
}

void TruncatedTspline::BezierFinder_1(int eid, vector<array<double, 4>> &be)
{
	int cid[4] = {tmesh[eid].cnct[0], tmesh[eid].cnct[1], tmesh[eid].cnct[2], tmesh[eid].cnct[3]};
	vector<int> pid = tmesh[eid].IEN;
	array<double, 4> e0 = {cp[cid[0]].knotU[2], cp[cid[1]].knotU[2], cp[cid[0]].knotV[2], cp[cid[3]].knotV[2]};
	be.clear();
	//be.push_back(e0);

	vector<double> ukt, vkt;
	for (uint i = 0; i < pid.size(); i++)
	{
		for (int j = 0; j < 5; j++)
		{
			if (cp[pid[i]].knotU[j] > e0[0] && cp[pid[i]].knotU[j] < e0[1])
			{
				vector<double>::iterator it = find(ukt.begin(), ukt.end(), cp[pid[i]].knotU[j]);
				if (it == ukt.end())
					ukt.push_back(cp[pid[i]].knotU[j]);
			}
			if (cp[pid[i]].knotV[j] > e0[2] && cp[pid[i]].knotV[j] < e0[3])
			{
				vector<double>::iterator it = find(vkt.begin(), vkt.end(), cp[pid[i]].knotV[j]);
				if (it == vkt.end())
					vkt.push_back(cp[pid[i]].knotV[j]);
			}
		}
	}
	if (ukt.size() != 0 && vkt.size() == 0)
	{
		sort(ukt.begin(), ukt.end());
		vector<double> ukt1;
		ukt1.push_back(e0[0]);
		for (uint i = 0; i < ukt.size(); i++)
		{
			ukt1.push_back(ukt[i]);
		}
		ukt1.push_back(e0[1]);
		for (uint i = 0; i < ukt1.size() - 1; i++)
		{
			array<double, 4> betmp = {ukt1[i], ukt1[i + 1], e0[2], e0[3]};
			be.push_back(betmp);
		}
	}
	else if (ukt.size() == 0 && vkt.size() != 0)
	{
		sort(vkt.begin(), vkt.end());
		vector<double> vkt1;
		vkt1.push_back(e0[2]);
		for (uint i = 0; i < vkt.size(); i++)
		{
			vkt1.push_back(vkt[i]);
		}
		vkt1.push_back(e0[3]);
		for (uint i = 0; i < vkt1.size() - 1; i++)
		{
			array<double, 4> betmp = {e0[0], e0[1], vkt1[i], vkt1[i + 1]};
			be.push_back(betmp);
		}
	}
	else if (ukt.size() == 0 && vkt.size() == 0)
	{
		be.push_back(e0);
	}
	else
	{
		cerr << "Not available for face intersection!\n";
		getchar();
	}

	/*int flag(1);
	while(flag)
	{
		flag=0;
		vector<array<double,4>> betmp;
		for(uint i=0;i<be.size();i++)
		{
			int eflag(0);
			for(uint j=0;j<pid.size();j++)
			{
				if(cp[pid[j]].knotU[0]<be[i][1] && cp[pid[j]].knotU[4]>be[i][0] && cp[pid[j]].knotV[0]<be[i][3] && cp[pid[j]].knotV[4]>be[i][2])
				{
					for(int k=0;k<4;k++)
					{
						if((cp[pid[j]].kitvU[k]!=0. && cp[pid[j]].kitvU[k]<be[i][1]-be[i][0]) || (cp[pid[j]].kitvV[k]!=0. && cp[pid[j]].kitvV[k]<be[i][3]-be[i][2]))
						{
							flag=1;
							eflag=1;
						}
					}
				}
			}
			if(eflag)
			{
				array<double,4> tmp1={be[i][0],(be[i][0]+be[i][1])/2.,be[i][2],(be[i][2]+be[i][3])/2.};
				array<double,4> tmp2={(be[i][0]+be[i][1])/2.,be[i][1],be[i][2],(be[i][2]+be[i][3])/2.};
				array<double,4> tmp3={be[i][0],(be[i][0]+be[i][1])/2.,(be[i][2]+be[i][3])/2.,be[i][3]};
				array<double,4> tmp4={(be[i][0]+be[i][1])/2.,be[i][1],(be[i][2]+be[i][3])/2.,be[i][3]};
				betmp.push_back(tmp1);
				betmp.push_back(tmp2);
				betmp.push_back(tmp3);
				betmp.push_back(tmp4);
			}
			else
			{
				betmp.push_back(be[i]);
			}
		}
		be.clear();
		be=betmp;
	}*/
}

bool TruncatedTspline::CheckSameKnotVector(const Vertex &p1, const Vertex &p2)
{
	for (int i = 0; i < 5; i++)
	{
		if (p1.knotU[i] != p2.knotU[i])
		{
			return false;
		}
		if (p1.knotV[i] != p2.knotV[i])
		{
			return false;
		}
	}
	return true;
}

void TruncatedTspline::RefineTest()
{
	int rid[2] = {54, 55};
	int trun1_id[6] = {59, 60, 61, 70, 71, 72};
	int trun2_id[10] = {48, 49, 50, 58, 62, 69, 73, 81, 82, 83};
	int ptt[6] = {121, 126, 122, 128, 125, 129}; //T-junctions
	int pte[1] = {124};							 //inner edge
	int ptf[2] = {123, 127};					 //face
	int ton[6][4] = {{59, 60, 48, 49}, {60, 61, 49, 50}, {59, 70, 58, 69}, {61, 72, 62, 73}, {70, 71, 81, 82}, {71, 72, 82, 83}};
	double tonc[4] = {5. / 12., 5. / 12., 1. / 12., 1. / 12.};
	int eon[1][6] = {{60, 71, 59, 70, 61, 72}};
	double eonc[6] = {3. / 8., 3. / 8., 1. / 16., 1. / 16., 1. / 16., 1. / 16.};
	int fon[2][4] = {{59, 60, 71, 70}, {60, 61, 72, 71}};
	double fonc[4] = {1. / 4., 1. / 4., 1. / 4., 1. / 4.};
	int enew[8][4] = {{59, 121, 123, 122}, {121, 60, 124, 123}, {123, 124, 71, 125}, {122, 123, 125, 70}, {60, 126, 127, 124}, {126, 61, 128, 127}, {127, 128, 72, 129}, {124, 127, 129, 71}};
	double ktsU[9] = {cp[59].knotU[0], cp[59].knotU[1], cp[59].knotU[2], (cp[59].knotU[2] + cp[59].knotU[3]) / 2., cp[59].knotU[3], (cp[59].knotU[3] + cp[59].knotU[4]) / 2., cp[59].knotU[4], cp[61].knotU[3], cp[61].knotU[4]};
	double ktsV[7] = {cp[59].knotV[0], cp[59].knotV[1], cp[59].knotV[2], (cp[59].knotV[2] + cp[59].knotV[3]) / 2., cp[59].knotV[3], cp[59].knotV[4], cp[70].knotV[4]};
	for (int i = 0; i < 9; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int setKU[3][5] = {{-1, 121, -1, 126, -1}, {122, 123, 124, 127, 128}, {-1, 125, -1, 129, -1}};
	int setKV[5][3] = {{-1, 122, -1}, {121, 123, 125}, {-1, 124, -1}, {126, 127, 129}, {-1, 128, -1}};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				if (setKU[i][j] != -1)
					cp[setKU[i][j]].knotU[k] = ktsU[j + k];
			}
		}
	}
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				if (setKV[i][j] != -1)
					cp[setKV[i][j]].knotV[k] = ktsV[j + k];
			}
		}
	}
	for (int i = 0; i < 6; i++)
	{
		cp[trun1_id[i]].type = 1;
		cp[trun1_id[i]].trun = 1;
	}
	for (int i = 0; i < 10; i++)
	{
		cp[trun2_id[i]].type = 2;
		cp[trun2_id[i]].trun = 1;
	}
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptt[i]].coor[0] += tonc[j] * cp[ton[i][j]].coor[0];
			cp[ptt[i]].coor[1] += tonc[j] * cp[ton[i][j]].coor[1];
			cp[ptt[i]].coor[2] += tonc[j] * cp[ton[i][j]].coor[2];
		}
	}
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cp[pte[i]].coor[0] += eonc[j] * cp[eon[i][j]].coor[0];
			cp[pte[i]].coor[1] += eonc[j] * cp[eon[i][j]].coor[1];
			cp[pte[i]].coor[2] += eonc[j] * cp[eon[i][j]].coor[2];
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptf[i]].coor[0] += fonc[j] * cp[fon[i][j]].coor[0];
			cp[ptf[i]].coor[1] += fonc[j] * cp[fon[i][j]].coor[1];
			cp[ptf[i]].coor[2] += fonc[j] * cp[fon[i][j]].coor[2];
		}
	}

	int t1id[8] = {48, 50, 58, 62, 69, 73, 81, 83};
	int t1cid[8] = {121, 126, 122, 128, 122, 128, 125, 129};
	int t2id[2] = {49, 82};
	int t2cid[2][2] = {{121, 126}, {125, 129}};
	int t3id[4] = {59, 61, 70, 72};
	int t3cid[4][4] = {{121, 122, 123, 124}, {126, 128, 127, 124}, {122, 125, 123, 124}, {128, 129, 127, 124}};
	double t3c[4] = {5. / 12., 5. / 12., 1. / 4., 1. / 16.};
	int t4id[2] = {60, 71};
	int t4cid[2][5] = {{121, 126, 124, 123, 127}, {125, 129, 124, 123, 127}};
	double t4c[5] = {5. / 12., 5. / 12., 3. / 8., 1. / 4., 1. / 4.};
	for (int i = 0; i < 8; i++)
	{
		cp[t1id[i]].tbf.push_back(t1cid[i]);
		cp[t1id[i]].tc.push_back(1. / 12.);
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			cp[t2id[i]].tbf.push_back(t2cid[i][j]);
			cp[t2id[i]].tc.push_back(1. / 12.);
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[t3id[i]].tbf.push_back(t3cid[i][j]);
			cp[t3id[i]].tc.push_back(t3c[j]);
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			cp[t4id[i]].tbf.push_back(t4cid[i][j]);
			cp[t4id[i]].tc.push_back(t4c[j]);
		}
	}

	//new elements
	int idtmp[9] = {121, 122, 123, 124, 125, 126, 127, 128, 129};
	for (int i = 0; i < 8; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		if (i < 4)
			etmp.IEN = tmesh[rid[0]].IEN;
		else
			etmp.IEN = tmesh[rid[1]].IEN;
		for (int j = 0; j < 9; j++)
			etmp.IEN.push_back(idtmp[j]);
		tmesh.push_back(etmp);
	}

	for (int i = 0; i < 2; i++)
		tmesh[rid[i]].act = 0;

	for (int i = 0; i < 100; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (int j = 0; j < 9; j++)
			{
				if (cp[idtmp[j]].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[idtmp[j]].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[idtmp[j]].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[idtmp[j]].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(idtmp[j]);
				}
			}
		}
	}
}

void TruncatedTspline::RefineTest1()
{
	int rid[1] = {54};
	int ptt[4] = {121, 122, 124, 125}; //T-junctions
	int ptf[1] = {123};				   //face
	int ptc[4] = {59, 60, 70, 71};	   //corner
	int ton[4][4] = {{59, 60, 48, 49}, {59, 70, 58, 69}, {60, 71, 61, 72}, {70, 71, 81, 82}};
	double tonc[4] = {5. / 12., 5. / 12., 1. / 12., 1. / 12.};
	int fon[1][4] = {{59, 60, 71, 70}};
	double fonc[4] = {1. / 4., 1. / 4., 1. / 4., 1. / 4.};
	int con[4][4] = {{59, 48, 58, 47}, {60, 49, 61, 50}, {70, 69, 81, 80}, {71, 72, 82, 83}};
	double conc[4] = {25. / 36., 5. / 36., 5. / 36., 1. / 36.};
	int enew[4][4] = {{59, 121, 123, 122}, {121, 60, 124, 123}, {123, 124, 71, 125}, {122, 123, 125, 70}};
	double ktsU[7] = {cp[59].knotU[0], cp[59].knotU[1], cp[59].knotU[2], (cp[59].knotU[2] + cp[59].knotU[3]) / 2., cp[59].knotU[3], cp[59].knotU[4], cp[60].knotU[4]};
	double ktsV[7] = {cp[59].knotV[0], cp[59].knotV[1], cp[59].knotV[2], (cp[59].knotV[2] + cp[59].knotV[3]) / 2., cp[59].knotV[3], cp[59].knotV[4], cp[70].knotV[4]};
	for (int i = 0; i < 5; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int setKU[3][3] = {{59, 121, 60}, {122, 123, 124}, {70, 125, 71}};
	int setKV[3][3] = {{59, 122, 70}, {121, 123, 125}, {60, 124, 71}};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKU[i][j]].knotU[k] = ktsU[j + k];
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKV[i][j]].knotV[k] = ktsV[j + k];
			}
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptt[i]].coor[0] += tonc[j] * cp[ton[i][j]].coor[0];
			cp[ptt[i]].coor[1] += tonc[j] * cp[ton[i][j]].coor[1];
			cp[ptt[i]].coor[2] += tonc[j] * cp[ton[i][j]].coor[2];
		}
	}
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptf[i]].coor[0] += fonc[j] * cp[fon[i][j]].coor[0];
			cp[ptf[i]].coor[1] += fonc[j] * cp[fon[i][j]].coor[1];
			cp[ptf[i]].coor[2] += fonc[j] * cp[fon[i][j]].coor[2];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 4; j++)
		{
			tmp.coor[0] += conc[j] * cp[con[i][j]].coor[0];
			tmp.coor[1] += conc[j] * cp[con[i][j]].coor[1];
			tmp.coor[2] += conc[j] * cp[con[i][j]].coor[2];
		}
		cp[ptc[i]].coor[0] = tmp.coor[0];
		cp[ptc[i]].coor[1] = tmp.coor[1];
		cp[ptc[i]].coor[2] = tmp.coor[2];
	}

	//int t1id[8]={48,49,58,61,69,72,81,82};
	//int t1cid[8][2]={{121,59},{121,60},{122,59},{124,60},{122,70},{124,71},{125,70},{125,71}};
	//double t1c[2]={1./12.,5./36.};
	int t2id[4] = {47, 50, 80, 83};
	int t2cid[4][1] = {{59}, {60}, {70}, {71}};
	double t2c[1] = {1. / 36.};
	//for(int i=0;i<8;i++)
	//{
	//	for(int j=0;j<2;j++)
	//	{
	//		cp[t1id[i]].tbf.push_back(t1cid[i][j]);
	//		cp[t1id[i]].tc.push_back(t1c[j]);
	//	}
	//}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t2id[i]].tbf.push_back(t2cid[i][j]);
			cp[t2id[i]].tc.push_back(t2c[j]);
		}
	}

	cp[48].knotV[4] = ktsV[3];
	cp[49].knotV[4] = ktsV[3];
	cp[81].knotV[0] = ktsV[3];
	cp[82].knotV[0] = ktsV[3];
	cp[58].knotU[4] = ktsU[3];
	cp[69].knotU[4] = ktsU[3];
	cp[61].knotU[0] = ktsU[3];
	cp[72].knotU[0] = ktsU[3];

	//new elements
	int idtmp[5] = {121, 122, 123, 124, 125};
	for (int i = 0; i < 4; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		etmp.IEN = tmesh[rid[0]].IEN;
		for (int j = 0; j < 5; j++)
			etmp.IEN.push_back(idtmp[j]);
		tmesh.push_back(etmp);
	}

	for (int i = 0; i < 1; i++)
		tmesh[rid[i]].act = 0;

	for (int i = 0; i < 100; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (int j = 0; j < 5; j++)
			{
				if (cp[idtmp[j]].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[idtmp[j]].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[idtmp[j]].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[idtmp[j]].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(idtmp[j]);
				}
			}
		}
	}
}

void TruncatedTspline::RefineTest2()
{
	int rid[1] = {54};
	int ptt[4] = {121, 122, 124, 125}; //T-junctions
	int ptf[1] = {123};				   //face
	//int ptc[4]={59,60,70,71};//corner
	int ton[4][4] = {{59, 60, 48, 49}, {59, 70, 58, 69}, {60, 71, 61, 72}, {70, 71, 81, 82}};
	double tonc[4] = {5. / 12., 5. / 12., 1. / 12., 1. / 12.};
	int fon[1][4] = {{59, 60, 71, 70}};
	double fonc[4] = {1. / 4., 1. / 4., 1. / 4., 1. / 4.};
	//int con[4][4]={{59,48,58,47},{60,49,61,50},{70,69,81,80},{71,72,82,83}};
	//double conc[4]={25./36.,5./36.,5./36.,1./36.};
	int enew[4][4] = {{59, 121, 123, 122}, {121, 60, 124, 123}, {123, 124, 71, 125}, {122, 123, 125, 70}};
	double ktsU[7] = {cp[59].knotU[0], cp[59].knotU[1], cp[59].knotU[2], (cp[59].knotU[2] + cp[59].knotU[3]) / 2., cp[59].knotU[3], cp[59].knotU[4], cp[60].knotU[4]};
	double ktsV[7] = {cp[59].knotV[0], cp[59].knotV[1], cp[59].knotV[2], (cp[59].knotV[2] + cp[59].knotV[3]) / 2., cp[59].knotV[3], cp[59].knotV[4], cp[70].knotV[4]};
	for (int i = 0; i < 5; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int setKU[3][3] = {{-1, 121, -1}, {122, 123, 124}, {-1, 125, -1}};
	int setKV[3][3] = {{-1, 122, -1}, {121, 123, 125}, {-1, 124, -1}};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				if (setKU[i][j] != -1)
					cp[setKU[i][j]].knotU[k] = ktsU[j + k];
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				if (setKV[i][j] != -1)
					cp[setKV[i][j]].knotV[k] = ktsV[j + k];
			}
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptt[i]].coor[0] += tonc[j] * cp[ton[i][j]].coor[0];
			cp[ptt[i]].coor[1] += tonc[j] * cp[ton[i][j]].coor[1];
			cp[ptt[i]].coor[2] += tonc[j] * cp[ton[i][j]].coor[2];
		}
	}
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptf[i]].coor[0] += fonc[j] * cp[fon[i][j]].coor[0];
			cp[ptf[i]].coor[1] += fonc[j] * cp[fon[i][j]].coor[1];
			cp[ptf[i]].coor[2] += fonc[j] * cp[fon[i][j]].coor[2];
		}
	}
	//for(int i=0;i<4;i++)
	//{
	//	Vertex tmp;
	//	for(int j=0;j<4;j++)
	//	{
	//		tmp.coor[0]+=conc[j]*cp[con[i][j]].coor[0];
	//		tmp.coor[1]+=conc[j]*cp[con[i][j]].coor[1];
	//		tmp.coor[2]+=conc[j]*cp[con[i][j]].coor[2];
	//	}
	//	cp[ptc[i]].coor[0]=tmp.coor[0]; cp[ptc[i]].coor[1]=tmp.coor[1]; cp[ptc[i]].coor[2]=tmp.coor[2];
	//}

	int t1id[8] = {48, 49, 58, 61, 69, 72, 81, 82};
	int t1cid[8][1] = {{121}, {121}, {122}, {124}, {122}, {124}, {125}, {125}};
	double t1c[1] = {1. / 12.};
	int t2id[4] = {59, 60, 70, 71};
	int t2cid[4][3] = {{121, 122, 123}, {121, 124, 123}, {122, 125, 123}, {124, 125, 123}};
	double t2c[3] = {5. / 12., 5. / 12., 1. / 4.};
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t1id[i]].tbf.push_back(t1cid[i][j]);
			cp[t1id[i]].tc.push_back(t1c[j]);
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			cp[t2id[i]].tbf.push_back(t2cid[i][j]);
			cp[t2id[i]].tc.push_back(t2c[j]);
		}
	}

	//new elements
	int idtmp[5] = {121, 122, 123, 124, 125};
	for (int i = 0; i < 4; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		etmp.IEN = tmesh[rid[0]].IEN;
		for (int j = 0; j < 5; j++)
			etmp.IEN.push_back(idtmp[j]);
		tmesh.push_back(etmp);
	}

	for (int i = 0; i < 1; i++)
		tmesh[rid[i]].act = 0;

	for (int i = 0; i < 100; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (int j = 0; j < 5; j++)
			{
				if (cp[idtmp[j]].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[idtmp[j]].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[idtmp[j]].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[idtmp[j]].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(idtmp[j]);
				}
			}
		}
	}
}

void TruncatedTspline::RefineTest3()
{
	int rid[2] = {54, 55};
	int ptt[6] = {121, 126, 122, 128, 125, 129}; //T-junctions
	int ptf[2] = {123, 127};					 //face
	int ptc[4] = {59, 61, 70, 72};				 //corner
	int pte[1] = {124};							 //edge
	int ptec[2] = {60, 71};
	int ton[6][4] = {{59, 60, 48, 49}, {60, 61, 49, 50}, {59, 70, 58, 69}, {61, 72, 62, 73}, {70, 71, 81, 82}, {71, 72, 82, 83}};
	double tonc[4] = {5. / 12., 5. / 12., 1. / 12., 1. / 12.};
	int fon[2][4] = {{59, 60, 71, 70}, {60, 61, 72, 71}};
	double fonc[4] = {1. / 4., 1. / 4., 1. / 4., 1. / 4.};
	int con[4][4] = {{59, 48, 58, 47}, {61, 50, 62, 51}, {70, 69, 81, 80}, {72, 73, 83, 84}};
	double conc[4] = {25. / 36., 5. / 36., 5. / 36., 1. / 36.};
	int eon[1][6] = {{60, 71, 59, 70, 61, 72}};
	double eonc[6] = {3. / 8., 3. / 8., 1. / 16., 1. / 16., 1. / 16., 1. / 16.};
	int econ[2][6] = {{60, 49, 59, 48, 61, 50}, {71, 82, 70, 81, 72, 83}};
	double econc[6] = {15. / 24., 3. / 24., 5. / 48., 1. / 48., 5. / 48., 1. / 48.};
	int enew[8][4] = {{59, 121, 123, 122}, {121, 60, 124, 123}, {123, 124, 71, 125}, {122, 123, 125, 70}, {60, 126, 127, 124}, {126, 61, 128, 127}, {127, 128, 72, 129}, {124, 127, 129, 71}};
	double ktsU[9] = {cp[59].knotU[0], cp[59].knotU[1], cp[59].knotU[2], (cp[59].knotU[2] + cp[59].knotU[3]) / 2., cp[59].knotU[3], (cp[59].knotU[3] + cp[59].knotU[4]) / 2., cp[59].knotU[4], cp[61].knotU[3], cp[61].knotU[4]};
	double ktsV[7] = {cp[59].knotV[0], cp[59].knotV[1], cp[59].knotV[2], (cp[59].knotV[2] + cp[59].knotV[3]) / 2., cp[59].knotV[3], cp[59].knotV[4], cp[70].knotV[4]};
	for (int i = 0; i < 9; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int setKU[3][5] = {{59, 121, 60, 126, 61}, {122, 123, 124, 127, 128}, {70, 125, 71, 129, 72}};
	int setKV[5][3] = {{59, 122, 70}, {121, 123, 125}, {60, 124, 71}, {126, 127, 129}, {61, 128, 72}};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKU[i][j]].knotU[k] = ktsU[j + k];
			}
		}
	}
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKV[i][j]].knotV[k] = ktsV[j + k];
			}
		}
	}
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptt[i]].coor[0] += tonc[j] * cp[ton[i][j]].coor[0];
			cp[ptt[i]].coor[1] += tonc[j] * cp[ton[i][j]].coor[1];
			cp[ptt[i]].coor[2] += tonc[j] * cp[ton[i][j]].coor[2];
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptf[i]].coor[0] += fonc[j] * cp[fon[i][j]].coor[0];
			cp[ptf[i]].coor[1] += fonc[j] * cp[fon[i][j]].coor[1];
			cp[ptf[i]].coor[2] += fonc[j] * cp[fon[i][j]].coor[2];
		}
	}
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cp[pte[i]].coor[0] += eonc[j] * cp[eon[i][j]].coor[0];
			cp[pte[i]].coor[1] += eonc[j] * cp[eon[i][j]].coor[1];
			cp[pte[i]].coor[2] += eonc[j] * cp[eon[i][j]].coor[2];
		}
	}
	for (int i = 0; i < 2; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 6; j++)
		{
			tmp.coor[0] += econc[j] * cp[econ[i][j]].coor[0];
			tmp.coor[1] += econc[j] * cp[econ[i][j]].coor[1];
			tmp.coor[2] += econc[j] * cp[econ[i][j]].coor[2];
		}
		cp[ptec[i]].coor[0] = tmp.coor[0];
		cp[ptec[i]].coor[1] = tmp.coor[1];
		cp[ptec[i]].coor[2] = tmp.coor[2];
	}
	for (int i = 0; i < 4; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 4; j++)
		{
			tmp.coor[0] += conc[j] * cp[con[i][j]].coor[0];
			tmp.coor[1] += conc[j] * cp[con[i][j]].coor[1];
			tmp.coor[2] += conc[j] * cp[con[i][j]].coor[2];
		}
		cp[ptc[i]].coor[0] = tmp.coor[0];
		cp[ptc[i]].coor[1] = tmp.coor[1];
		cp[ptc[i]].coor[2] = tmp.coor[2];
	}

	int t1id[4] = {48, 50, 81, 83};
	int t1cid[4][3] = {{59, 121, 60}, {61, 126, 60}, {70, 125, 71}, {72, 129, 71}};
	double t1c[3] = {5. / 36., 1. / 12., 1. / 48.};
	int t2id[4] = {47, 51, 80, 84};
	int t2cid[4][1] = {{59}, {61}, {70}, {72}};
	double t2c[1] = {1. / 36.};
	int t3id[2] = {49, 82};
	int t3cid[2][3] = {{60, 121, 126}, {71, 125, 129}};
	double t3c[3] = {3. / 24., 1. / 12., 1. / 12.};
	int t4id[4] = {58, 62, 69, 73};
	int t4cid[4][2] = {{59, 122}, {61, 128}, {70, 122}, {72, 128}};
	double t4c[2] = {5. / 36., 1. / 12.};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			cp[t1id[i]].tbf.push_back(t1cid[i][j]);
			cp[t1id[i]].tc.push_back(t1c[j]);
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t2id[i]].tbf.push_back(t2cid[i][j]);
			cp[t2id[i]].tc.push_back(t2c[j]);
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			cp[t3id[i]].tbf.push_back(t3cid[i][j]);
			cp[t3id[i]].tc.push_back(t3c[j]);
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			cp[t4id[i]].tbf.push_back(t4cid[i][j]);
			cp[t4id[i]].tc.push_back(t4c[j]);
		}
	}

	//new elements
	int idtmp[9] = {121, 122, 123, 124, 125, 126, 127, 128, 129};
	for (int i = 0; i < 8; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		if (i < 4)
			etmp.IEN = tmesh[rid[0]].IEN;
		else
			etmp.IEN = tmesh[rid[1]].IEN;
		for (int j = 0; j < 9; j++)
			etmp.IEN.push_back(idtmp[j]);
		tmesh.push_back(etmp);
	}

	for (int i = 0; i < 2; i++)
		tmesh[rid[i]].act = 0;

	for (int i = 0; i < 100; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (int j = 0; j < 9; j++)
			{
				if (cp[idtmp[j]].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[idtmp[j]].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[idtmp[j]].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[idtmp[j]].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(idtmp[j]);
				}
			}
		}
	}
}

void TruncatedTspline::RefineTest4()
{
	int rid[3] = {54, 55, 45};
	for (int i = 0; i < 13; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int ptt[8] = {121, 130, 122, 128, 125, 129, 131, 133}; //T-junctions
	int ton[8][4] = {{59, 60, 48, 49}, {49, 50, 38, 39}, {59, 70, 58, 69}, {61, 72, 62, 73}, {70, 71, 81, 82}, {71, 72, 82, 83}, {49, 60, 48, 59}, {50, 61, 51, 62}};
	double tonc[4] = {5. / 12., 5. / 12., 1. / 12., 1. / 12.};
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptt[i]].coor[0] += tonc[j] * cp[ton[i][j]].coor[0];
			cp[ptt[i]].coor[1] += tonc[j] * cp[ton[i][j]].coor[1];
			cp[ptt[i]].coor[2] += tonc[j] * cp[ton[i][j]].coor[2];
		}
	}
	int ptf[3] = {123, 127, 132}; //face
	int fon[3][4] = {{59, 60, 71, 70}, {60, 61, 72, 71}, {49, 50, 61, 60}};
	double fonc[4] = {1. / 4., 1. / 4., 1. / 4., 1. / 4.};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptf[i]].coor[0] += fonc[j] * cp[fon[i][j]].coor[0];
			cp[ptf[i]].coor[1] += fonc[j] * cp[fon[i][j]].coor[1];
			cp[ptf[i]].coor[2] += fonc[j] * cp[fon[i][j]].coor[2];
		}
	}
	int pte[2] = {124, 126}; //edge
	int eon[2][6] = {{60, 71, 59, 70, 61, 72}, {60, 61, 49, 50, 71, 72}};
	double eonc[6] = {3. / 8., 3. / 8., 1. / 16., 1. / 16., 1. / 16., 1. / 16.};
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cp[pte[i]].coor[0] += eonc[j] * cp[eon[i][j]].coor[0];
			cp[pte[i]].coor[1] += eonc[j] * cp[eon[i][j]].coor[1];
			cp[pte[i]].coor[2] += eonc[j] * cp[eon[i][j]].coor[2];
		}
	}
	int ptcc[1] = {60}; //corner corner
	int ccon[1][9] = {60, 49, 59, 61, 71, 48, 50, 70, 72};
	double cconc[9] = {9. / 16., 3. / 32., 3. / 32., 3. / 32., 3. / 32., 1. / 64., 1. / 64., 1. / 64., 1. / 64.};
	for (int i = 0; i < 1; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 9; j++)
		{
			tmp.coor[0] += cconc[j] * cp[ccon[i][j]].coor[0];
			tmp.coor[1] += cconc[j] * cp[ccon[i][j]].coor[1];
			tmp.coor[2] += cconc[j] * cp[ccon[i][j]].coor[2];
		}
		cp[ptcc[i]].coor[0] = tmp.coor[0];
		cp[ptcc[i]].coor[1] = tmp.coor[1];
		cp[ptcc[i]].coor[2] = tmp.coor[2];
	}
	int ptec[2] = {61, 71}; //edge corner
	int econ[2][6] = {{61, 62, 50, 51, 72, 73}, {71, 82, 70, 81, 72, 83}};
	double econc[6] = {15. / 24., 3. / 24., 5. / 48., 1. / 48., 5. / 48., 1. / 48.};
	for (int i = 0; i < 2; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 6; j++)
		{
			tmp.coor[0] += econc[j] * cp[econ[i][j]].coor[0];
			tmp.coor[1] += econc[j] * cp[econ[i][j]].coor[1];
			tmp.coor[2] += econc[j] * cp[econ[i][j]].coor[2];
		}
		cp[ptec[i]].coor[0] = tmp.coor[0];
		cp[ptec[i]].coor[1] = tmp.coor[1];
		cp[ptec[i]].coor[2] = tmp.coor[2];
	}
	int ptc[5] = {59, 49, 70, 72, 50}; //corner
	int con[5][4] = {{59, 48, 58, 47}, {49, 38, 48, 37}, {70, 69, 81, 80}, {72, 73, 83, 84}, {50, 39, 51, 40}};
	double conc[4] = {25. / 36., 5. / 36., 5. / 36., 1. / 36.};
	for (int i = 0; i < 5; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 4; j++)
		{
			tmp.coor[0] += conc[j] * cp[con[i][j]].coor[0];
			tmp.coor[1] += conc[j] * cp[con[i][j]].coor[1];
			tmp.coor[2] += conc[j] * cp[con[i][j]].coor[2];
		}
		cp[ptc[i]].coor[0] = tmp.coor[0];
		cp[ptc[i]].coor[1] = tmp.coor[1];
		cp[ptc[i]].coor[2] = tmp.coor[2];
	}

	//update knot vectors
	double ktsU[9] = {cp[59].knotU[0], cp[59].knotU[1], cp[59].knotU[2], (cp[59].knotU[2] + cp[59].knotU[3]) / 2., cp[59].knotU[3], (cp[59].knotU[3] + cp[59].knotU[4]) / 2., cp[59].knotU[4], cp[61].knotU[3], cp[61].knotU[4]};
	double ktsV[7] = {cp[59].knotV[0], cp[59].knotV[1], cp[59].knotV[2], (cp[59].knotV[2] + cp[59].knotV[3]) / 2., cp[59].knotV[3], cp[59].knotV[4], cp[70].knotV[4]};
	double ktsU1[7] = {cp[49].knotU[0], cp[49].knotU[1], cp[49].knotU[2], (cp[49].knotU[2] + cp[49].knotU[3]) / 2., cp[49].knotU[3], cp[49].knotU[4], cp[50].knotU[4]};
	double ktsV1[9] = {cp[49].knotV[0], cp[49].knotV[1], cp[49].knotV[2], (cp[49].knotV[2] + cp[49].knotV[3]) / 2., cp[49].knotV[3], (cp[49].knotV[3] + cp[49].knotV[4]) / 2., cp[49].knotV[4], cp[71].knotV[3], cp[71].knotV[4]};
	int setKU[3][5] = {{59, 121, 60, 126, 61}, {122, 123, 124, 127, 128}, {70, 125, 71, 129, 72}};
	int setKV[2][3] = {{59, 122, 70}, {121, 123, 125}};
	int setKU1[2][3] = {{49, 130, 50}, {131, 132, 133}};
	int setKV1[3][5] = {{49, 131, 60, 124, 71}, {130, 132, 126, 127, 129}, {50, 133, 61, 128, 72}};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKU[i][j]].knotU[k] = ktsU[j + k];
			}
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKV[i][j]].knotV[k] = ktsV[j + k];
			}
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKU1[i][j]].knotU[k] = ktsU1[j + k];
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKV1[i][j]].knotV[k] = ktsV1[j + k];
			}
		}
	}
	cp[38].knotV[4] = cp[131].knotV[2];
	cp[39].knotV[4] = cp[133].knotV[2];
	//cp[48].knotU[4]=cp[130].knotU[2]; cp[48].knotV[4]=cp[122].knotV[2];
	cp[51].knotU[0] = cp[130].knotU[2];
	cp[58].knotU[4] = cp[121].knotU[2];
	cp[62].knotU[0] = cp[126].knotU[2];
	cp[69].knotU[4] = cp[125].knotU[2];
	cp[73].knotU[0] = cp[129].knotU[2];
	cp[81].knotV[0] = cp[122].knotV[2];
	cp[82].knotV[0] = cp[124].knotV[2];
	cp[83].knotV[0] = cp[128].knotV[2];
	cp[131].knotU[0] = cp[60].knotU[0];
	cp[131].knotU[1] = cp[60].knotU[1];
	cp[121].knotV[0] = cp[60].knotV[0];
	cp[121].knotV[1] = cp[60].knotV[1];

	//int t1id[4]={48,50,81,83};
	//int t1cid[4][3]={{59,121,60},{61,126,60},{70,125,71},{72,129,71}};
	//double t1c[3]={5./36.,1./12.,1./48.};
	int t2id[5] = {47, 37, 40, 80, 84};
	int t2cid[5][1] = {{59}, {49}, {50}, {70}, {72}};
	double t2c[1] = {1. / 36.};
	//int t3id[2]={49,82};
	//int t3cid[2][3]={{60,121,126},{71,125,129}};
	//double t3c[3]={3./24.,1./12.,1./12.};
	//int t4id[4]={58,62,69,73};
	//int t4cid[4][2]={{59,122},{61,128},{70,122},{72,128}};
	//double t4c[2]={5./36.,1./12.};
	int t5id[1] = {48};
	int t5cid[1][3] = {{60, 121, 131}};
	double t5c[3] = {1. / 64., 1. / 12., 1. / 12.};
	//for(int i=0;i<4;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t1id[i]].tbf.push_back(t1cid[i][j]);
	//		cp[t1id[i]].tc.push_back(t1c[j]);
	//	}
	//}
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t2id[i]].tbf.push_back(t2cid[i][j]);
			cp[t2id[i]].tc.push_back(t2c[j]);
		}
	}
	//for(int i=0;i<2;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t3id[i]].tbf.push_back(t3cid[i][j]);
	//		cp[t3id[i]].tc.push_back(t3c[j]);
	//	}
	//}
	//for(int i=0;i<4;i++)
	//{
	//	for(int j=0;j<2;j++)
	//	{
	//		cp[t4id[i]].tbf.push_back(t4cid[i][j]);
	//		cp[t4id[i]].tc.push_back(t4c[j]);
	//	}
	//}
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			cp[t5id[i]].tbf.push_back(t5cid[i][j]);
			cp[t5id[i]].tc.push_back(t5c[j]);
		}
	}

	//new elements
	int enew[12][4] = {{59, 121, 123, 122}, {121, 60, 124, 123}, {123, 124, 71, 125}, {122, 123, 125, 70}, {60, 126, 127, 124}, {126, 61, 128, 127}, {127, 128, 72, 129}, {124, 127, 129, 71}, {49, 130, 132, 131}, {130, 50, 133, 132}, {132, 133, 61, 126}, {131, 132, 126, 60}};
	int idtmp[13] = {121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133};
	for (int i = 0; i < 12; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		if (i < 4)
			etmp.IEN = tmesh[rid[0]].IEN;
		else if (i >= 4 && i < 8)
			etmp.IEN = tmesh[rid[1]].IEN;
		else
			etmp.IEN = tmesh[rid[2]].IEN;
		for (int j = 0; j < 13; j++)
			etmp.IEN.push_back(idtmp[j]);
		tmesh.push_back(etmp);
	}

	for (int i = 0; i < 3; i++)
		tmesh[rid[i]].act = 0;

	for (int i = 0; i < 100; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (int j = 0; j < 13; j++)
			{
				if (cp[idtmp[j]].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[idtmp[j]].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[idtmp[j]].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[idtmp[j]].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(idtmp[j]);
				}
			}
		}
	}
}

void TruncatedTspline::RefineTest3_1()
{
	int rid[2] = {54, 55};
	for (int i = 0; i < 9; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int ptt[6] = {121, 126, 122, 128, 125, 129}; //T-junctions
	int ton[6][4] = {{59, 60, 48, 49}, {60, 61, 49, 50}, {59, 70, 58, 69}, {61, 72, 62, 73}, {70, 71, 81, 82}, {71, 72, 82, 83}};
	double tonc[4] = {5. / 12., 5. / 12., 1. / 12., 1. / 12.};
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptt[i]].coor[0] += tonc[j] * cp[ton[i][j]].coor[0];
			cp[ptt[i]].coor[1] += tonc[j] * cp[ton[i][j]].coor[1];
			cp[ptt[i]].coor[2] += tonc[j] * cp[ton[i][j]].coor[2];
		}
	}
	int ptf[2] = {123, 127}; //face
	int fon[2][4] = {{59, 60, 71, 70}, {60, 61, 72, 71}};
	double fonc[4] = {1. / 4., 1. / 4., 1. / 4., 1. / 4.};
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptf[i]].coor[0] += fonc[j] * cp[fon[i][j]].coor[0];
			cp[ptf[i]].coor[1] += fonc[j] * cp[fon[i][j]].coor[1];
			cp[ptf[i]].coor[2] += fonc[j] * cp[fon[i][j]].coor[2];
		}
	}
	int pte[1] = {124}; //edge
	int eon[1][6] = {{60, 71, 59, 70, 61, 72}};
	double eonc[6] = {3. / 8., 3. / 8., 1. / 16., 1. / 16., 1. / 16., 1. / 16.};
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cp[pte[i]].coor[0] += eonc[j] * cp[eon[i][j]].coor[0];
			cp[pte[i]].coor[1] += eonc[j] * cp[eon[i][j]].coor[1];
			cp[pte[i]].coor[2] += eonc[j] * cp[eon[i][j]].coor[2];
		}
	}
	//int ptcc[1]={60};//corner corner
	//int ccon[1][9]={60,49,59,61,71,48,50,70,72};
	//double cconc[9]={9./16.,3./32.,3./32.,3./32.,3./32.,1./64.,1./64.,1./64.,1./64.};
	//for(int i=0;i<1;i++)
	//{
	//	Vertex tmp;
	//	for(int j=0;j<9;j++)
	//	{
	//		tmp.coor[0]+=cconc[j]*cp[ccon[i][j]].coor[0];
	//		tmp.coor[1]+=cconc[j]*cp[ccon[i][j]].coor[1];
	//		tmp.coor[2]+=cconc[j]*cp[ccon[i][j]].coor[2];
	//	}
	//	cp[ptcc[i]].coor[0]=tmp.coor[0]; cp[ptcc[i]].coor[1]=tmp.coor[1]; cp[ptcc[i]].coor[2]=tmp.coor[2];
	//}
	int ptec[2] = {60, 71}; //edge corner
	int econ[2][6] = {{60, 49, 59, 48, 61, 50}, {71, 82, 70, 81, 72, 83}};
	double econc[6] = {15. / 24., 3. / 24., 5. / 48., 1. / 48., 5. / 48., 1. / 48.};
	for (int i = 0; i < 2; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 6; j++)
		{
			tmp.coor[0] += econc[j] * cp[econ[i][j]].coor[0];
			tmp.coor[1] += econc[j] * cp[econ[i][j]].coor[1];
			tmp.coor[2] += econc[j] * cp[econ[i][j]].coor[2];
		}
		cp[ptec[i]].coor[0] = tmp.coor[0];
		cp[ptec[i]].coor[1] = tmp.coor[1];
		cp[ptec[i]].coor[2] = tmp.coor[2];
	}
	int ptc[4] = {59, 61, 70, 72}; //corner
	int con[4][4] = {{59, 48, 58, 47}, {61, 50, 62, 51}, {70, 69, 81, 80}, {72, 73, 83, 84}};
	double conc[4] = {25. / 36., 5. / 36., 5. / 36., 1. / 36.};
	for (int i = 0; i < 4; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 4; j++)
		{
			tmp.coor[0] += conc[j] * cp[con[i][j]].coor[0];
			tmp.coor[1] += conc[j] * cp[con[i][j]].coor[1];
			tmp.coor[2] += conc[j] * cp[con[i][j]].coor[2];
		}
		cp[ptc[i]].coor[0] = tmp.coor[0];
		cp[ptc[i]].coor[1] = tmp.coor[1];
		cp[ptc[i]].coor[2] = tmp.coor[2];
	}

	//update knot vectors
	double ktsU[9] = {cp[59].knotU[0], cp[59].knotU[1], cp[59].knotU[2], (cp[59].knotU[2] + cp[59].knotU[3]) / 2., cp[59].knotU[3], (cp[59].knotU[3] + cp[59].knotU[4]) / 2., cp[59].knotU[4], cp[61].knotU[3], cp[61].knotU[4]};
	double ktsV[7] = {cp[59].knotV[0], cp[59].knotV[1], cp[59].knotV[2], (cp[59].knotV[2] + cp[59].knotV[3]) / 2., cp[59].knotV[3], cp[59].knotV[4], cp[70].knotV[4]};
	//double ktsU1[7]={cp[49].knotU[0],cp[49].knotU[1],cp[49].knotU[2],(cp[49].knotU[2]+cp[49].knotU[3])/2.,cp[49].knotU[3],cp[49].knotU[4],cp[50].knotU[4]};
	//double ktsV1[9]={cp[49].knotV[0],cp[49].knotV[1],cp[49].knotV[2],(cp[49].knotV[2]+cp[49].knotV[3])/2.,cp[49].knotV[3],(cp[49].knotV[3]+cp[49].knotV[4])/2.,cp[49].knotV[4],cp[71].knotV[3],cp[71].knotV[4]};
	int setKU[3][5] = {{59, 121, 60, 126, 61}, {122, 123, 124, 127, 128}, {70, 125, 71, 129, 72}};
	int setKV[5][3] = {{59, 122, 70}, {121, 123, 125}, {60, 124, 71}, {126, 127, 129}, {61, 128, 72}};
	//int setKU1[2][3]={{49,130,50},{131,132,133}};
	//int setKV1[3][5]={{49,131,60,124,71},{130,132,126,127,129},{50,133,61,128,72}};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKU[i][j]].knotU[k] = ktsU[j + k];
			}
		}
	}
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKV[i][j]].knotV[k] = ktsV[j + k];
			}
		}
	}
	//for(int i=0;i<2;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		for(int k=0;k<5;k++)
	//		{
	//			cp[setKU1[i][j]].knotU[k]=ktsU1[j+k];
	//		}
	//	}
	//}
	//for(int i=0;i<3;i++)
	//{
	//	for(int j=0;j<5;j++)
	//	{
	//		for(int k=0;k<5;k++)
	//		{
	//			cp[setKV1[i][j]].knotV[k]=ktsV1[j+k];
	//		}
	//	}
	//}

	cp[48].knotV[4] = cp[122].knotV[2];
	cp[49].knotV[4] = cp[124].knotV[2];
	cp[50].knotV[4] = cp[128].knotV[2];
	cp[58].knotU[4] = cp[121].knotU[2];
	cp[62].knotU[0] = cp[126].knotU[2];
	cp[69].knotU[4] = cp[125].knotU[2];
	cp[73].knotU[0] = cp[129].knotU[2];
	cp[81].knotV[0] = cp[122].knotV[2];
	cp[82].knotV[0] = cp[124].knotV[2];
	cp[83].knotV[0] = cp[128].knotV[2];

	//int t1id[4]={48,50,81,83};
	//int t1cid[4][3]={{59,121,60},{61,126,60},{70,125,71},{72,129,71}};
	//double t1c[3]={5./36.,1./12.,1./48.};
	int t2id[4] = {47, 51, 80, 84};
	int t2cid[4][1] = {{59}, {61}, {70}, {72}};
	double t2c[1] = {1. / 36.};
	//int t3id[2]={49,82};
	//int t3cid[2][3]={{60,121,126},{71,125,129}};
	//double t3c[3]={3./24.,1./12.,1./12.};
	//int t4id[4]={58,62,69,73};
	//int t4cid[4][2]={{59,122},{61,128},{70,122},{72,128}};
	//double t4c[2]={5./36.,1./12.};
	//int t5id[1]={48};
	//int t5cid[1][3]={{60,121,131}};
	//double t5c[3]={1./64.,1./12.,1./12.};
	//for(int i=0;i<4;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t1id[i]].tbf.push_back(t1cid[i][j]);
	//		cp[t1id[i]].tc.push_back(t1c[j]);
	//	}
	//}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t2id[i]].tbf.push_back(t2cid[i][j]);
			cp[t2id[i]].tc.push_back(t2c[j]);
		}
	}
	//for(int i=0;i<2;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t3id[i]].tbf.push_back(t3cid[i][j]);
	//		cp[t3id[i]].tc.push_back(t3c[j]);
	//	}
	//}
	//for(int i=0;i<4;i++)
	//{
	//	for(int j=0;j<2;j++)
	//	{
	//		cp[t4id[i]].tbf.push_back(t4cid[i][j]);
	//		cp[t4id[i]].tc.push_back(t4c[j]);
	//	}
	//}
	//for(int i=0;i<1;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t5id[i]].tbf.push_back(t5cid[i][j]);
	//		cp[t5id[i]].tc.push_back(t5c[j]);
	//	}
	//}

	//new elements
	int enew[8][4] = {{59, 121, 123, 122}, {121, 60, 124, 123}, {123, 124, 71, 125}, {122, 123, 125, 70}, {60, 126, 127, 124}, {126, 61, 128, 127}, {127, 128, 72, 129}, {124, 127, 129, 71}};
	int idtmp[9] = {121, 122, 123, 124, 125, 126, 127, 128, 129};
	for (int i = 0; i < 8; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		if (i < 4)
			etmp.IEN = tmesh[rid[0]].IEN;
		else if (i >= 4 && i < 8)
			etmp.IEN = tmesh[rid[1]].IEN;
		//else etmp.IEN=tmesh[rid[2]].IEN;
		for (int j = 0; j < 9; j++)
			etmp.IEN.push_back(idtmp[j]);
		tmesh.push_back(etmp);
	}

	for (int i = 0; i < 2; i++)
		tmesh[rid[i]].act = 0;

	for (int i = 0; i < 100; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (int j = 0; j < 9; j++)
			{
				if (cp[idtmp[j]].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[idtmp[j]].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[idtmp[j]].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[idtmp[j]].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(idtmp[j]);
				}
			}
		}
	}
}

void TruncatedTspline::RefineTest5()
{
	int rid[3] = {54, 55, 45};
	for (int i = 0; i < 13; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int ptt[6] = {130, 122, 128, 125, 129, 133}; //T-junctions
	int ton[6][4] = {{49, 50, 38, 39}, {59, 70, 58, 69}, {61, 72, 62, 73}, {70, 71, 81, 82}, {71, 72, 82, 83}, {50, 61, 51, 62}};
	double tonc[4] = {5. / 12., 5. / 12., 1. / 12., 1. / 12.};
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptt[i]].coor[0] += tonc[j] * cp[ton[i][j]].coor[0];
			cp[ptt[i]].coor[1] += tonc[j] * cp[ton[i][j]].coor[1];
			cp[ptt[i]].coor[2] += tonc[j] * cp[ton[i][j]].coor[2];
		}
	}
	int ptf[3] = {123, 127, 132}; //face
	int fon[3][4] = {{59, 60, 71, 70}, {60, 61, 72, 71}, {49, 50, 61, 60}};
	double fonc[4] = {1. / 4., 1. / 4., 1. / 4., 1. / 4.};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptf[i]].coor[0] += fonc[j] * cp[fon[i][j]].coor[0];
			cp[ptf[i]].coor[1] += fonc[j] * cp[fon[i][j]].coor[1];
			cp[ptf[i]].coor[2] += fonc[j] * cp[fon[i][j]].coor[2];
		}
	}
	int pte[4] = {124, 126, 121, 131}; //edge
	int eon[4][6] = {{60, 71, 59, 70, 61, 72}, {60, 61, 49, 50, 71, 72}, {59, 60, 48, 49, 70, 71}, {49, 60, 48, 59, 50, 61}};
	double eonc[6] = {3. / 8., 3. / 8., 1. / 16., 1. / 16., 1. / 16., 1. / 16.};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cp[pte[i]].coor[0] += eonc[j] * cp[eon[i][j]].coor[0];
			cp[pte[i]].coor[1] += eonc[j] * cp[eon[i][j]].coor[1];
			cp[pte[i]].coor[2] += eonc[j] * cp[eon[i][j]].coor[2];
		}
	}
	int ptcc[1] = {60}; //corner corner
	int ccon[1][9] = {60, 49, 59, 61, 71, 48, 50, 70, 72};
	double cconc[9] = {9. / 16., 3. / 32., 3. / 32., 3. / 32., 3. / 32., 1. / 64., 1. / 64., 1. / 64., 1. / 64.};
	for (int i = 0; i < 1; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 9; j++)
		{
			tmp.coor[0] += cconc[j] * cp[ccon[i][j]].coor[0];
			tmp.coor[1] += cconc[j] * cp[ccon[i][j]].coor[1];
			tmp.coor[2] += cconc[j] * cp[ccon[i][j]].coor[2];
		}
		cp[ptcc[i]].coor[0] = tmp.coor[0];
		cp[ptcc[i]].coor[1] = tmp.coor[1];
		cp[ptcc[i]].coor[2] = tmp.coor[2];
	}
	int ptec[2] = {61, 71}; //edge corner
	int econ[2][6] = {{61, 62, 50, 51, 72, 73}, {71, 82, 70, 81, 72, 83}};
	double econc[6] = {15. / 24., 3. / 24., 5. / 48., 1. / 48., 5. / 48., 1. / 48.};
	for (int i = 0; i < 2; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 6; j++)
		{
			tmp.coor[0] += econc[j] * cp[econ[i][j]].coor[0];
			tmp.coor[1] += econc[j] * cp[econ[i][j]].coor[1];
			tmp.coor[2] += econc[j] * cp[econ[i][j]].coor[2];
		}
		cp[ptec[i]].coor[0] = tmp.coor[0];
		cp[ptec[i]].coor[1] = tmp.coor[1];
		cp[ptec[i]].coor[2] = tmp.coor[2];
	}
	int ptc[3] = {70, 72, 50}; //corner
	int con[3][4] = {{70, 69, 81, 80}, {72, 73, 83, 84}, {50, 39, 51, 40}};
	double conc[4] = {25. / 36., 5. / 36., 5. / 36., 1. / 36.};
	for (int i = 0; i < 3; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 4; j++)
		{
			tmp.coor[0] += conc[j] * cp[con[i][j]].coor[0];
			tmp.coor[1] += conc[j] * cp[con[i][j]].coor[1];
			tmp.coor[2] += conc[j] * cp[con[i][j]].coor[2];
		}
		cp[ptc[i]].coor[0] = tmp.coor[0];
		cp[ptc[i]].coor[1] = tmp.coor[1];
		cp[ptc[i]].coor[2] = tmp.coor[2];
	}

	//update knot vectors
	double ktsU[9] = {cp[59].knotU[0], cp[59].knotU[1], cp[59].knotU[2], (cp[59].knotU[2] + cp[59].knotU[3]) / 2., cp[59].knotU[3], (cp[59].knotU[3] + cp[59].knotU[4]) / 2., cp[59].knotU[4], cp[61].knotU[3], cp[61].knotU[4]};
	double ktsV[7] = {cp[59].knotV[0], cp[59].knotV[1], cp[59].knotV[2], (cp[59].knotV[2] + cp[59].knotV[3]) / 2., cp[59].knotV[3], cp[59].knotV[4], cp[70].knotV[4]};
	double ktsU1[7] = {cp[49].knotU[0], cp[49].knotU[1], cp[49].knotU[2], (cp[49].knotU[2] + cp[49].knotU[3]) / 2., cp[49].knotU[3], cp[49].knotU[4], cp[50].knotU[4]};
	double ktsV1[9] = {cp[49].knotV[0], cp[49].knotV[1], cp[49].knotV[2], (cp[49].knotV[2] + cp[49].knotV[3]) / 2., cp[49].knotV[3], (cp[49].knotV[3] + cp[49].knotV[4]) / 2., cp[49].knotV[4], cp[71].knotV[3], cp[71].knotV[4]};
	int setKU[3][5] = {{-1, 121, 60, 126, 61}, {122, 123, 124, 127, 128}, {70, 125, 71, 129, 72}};
	int setKV[2][3] = {{-1, 122, 70}, {121, 123, 125}};
	int setKU1[2][3] = {{-1, 130, 50}, {131, 132, 133}};
	int setKV1[3][5] = {{-1, 131, 60, 124, 71}, {130, 132, 126, 127, 129}, {50, 133, 61, 128, 72}};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				if (setKU[i][j] != -1)
					cp[setKU[i][j]].knotU[k] = ktsU[j + k];
			}
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				if (setKV[i][j] != -1)
					cp[setKV[i][j]].knotV[k] = ktsV[j + k];
			}
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				if (setKU1[i][j] != -1)
					cp[setKU1[i][j]].knotU[k] = ktsU1[j + k];
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				if (setKV1[i][j] != -1)
					cp[setKV1[i][j]].knotV[k] = ktsV1[j + k];
			}
		}
	}
	//cp[38].knotV[4]=cp[131].knotV[2];
	//cp[39].knotV[4]=cp[133].knotV[2];
	//cp[48].knotU[4]=cp[130].knotU[2]; cp[48].knotV[4]=cp[122].knotV[2];
	cp[51].knotU[0] = cp[130].knotU[2];
	//cp[58].knotU[4]=cp[121].knotU[2];
	cp[62].knotU[0] = cp[126].knotU[2];
	//cp[69].knotU[4]=cp[125].knotU[2];
	cp[73].knotU[0] = cp[129].knotU[2];
	cp[81].knotV[0] = cp[122].knotV[2];
	cp[82].knotV[0] = cp[124].knotV[2];
	cp[83].knotV[0] = cp[128].knotV[2];
	//cp[131].knotU[0]=cp[60].knotU[0]; cp[131].knotU[1]=cp[60].knotU[1];
	//cp[121].knotV[0]=cp[60].knotV[0]; cp[121].knotV[1]=cp[60].knotV[1];

	//cp[59].knotV[0]=cp[48].knotV[2]; cp[59].knotV[1]=cp[131].knotV[2];
	cp[121].knotV[0] = cp[48].knotV[2];
	cp[121].knotV[1] = cp[131].knotV[2];
	//cp[49].knotU[0]=cp[48].knotU[2]; cp[49].knotU[1]=cp[121].knotU[2];
	cp[131].knotU[0] = cp[48].knotU[2];
	cp[131].knotU[1] = cp[121].knotU[2];
	cp[122].knotV[0] = cp[131].knotV[2];
	cp[123].knotV[0] = cp[131].knotV[2];
	//cp[130].knotU[0]=cp[121].knotU[2];
	cp[132].knotU[0] = cp[121].knotU[2];

	int t1id[1] = {48};
	int t1cid[1][3] = {{121, 131, 60}};
	double t1c[3] = {1. / 16., 1. / 16., 1. / 64.};
	int t2id[2] = {49, 59};
	int t2cid[2][5] = {{131, 60, 130, 132, 126}, {121, 60, 122, 123, 124}};
	double t2c[6] = {3. / 8., 3. / 32., 5. / 12., 1. / 4., 1. / 16.};
	int t3id[2] = {38, 58};
	int t3cid[2][1] = {{130}, {122}};
	double t3c[1] = {1. / 12.};
	int t4id[3] = {40, 80, 84};
	int t4cid[3][1] = {{50}, {70}, {72}};
	double t4c[1] = {1. / 36.};
	int t5id[2] = {39, 69};
	int t5cid[2][2] = {{50, 130}, {70, 122}};
	double t5c[2] = {5. / 36., 1. / 12.};
	//int t6id[2]={47,37};
	//int t6cid[2][1]={{59},{49}};
	//double t6c[1]={1./48.};
	//int t7id[2]={58,38};
	//int t7cid[2][2]={{59,122},{49,130}};
	//double t7c[2]={3./24.,1./12.};
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			cp[t1id[i]].tbf.push_back(t1cid[i][j]);
			cp[t1id[i]].tc.push_back(t1c[j]);
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			cp[t2id[i]].tbf.push_back(t2cid[i][j]);
			cp[t2id[i]].tc.push_back(t2c[j]);
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t3id[i]].tbf.push_back(t3cid[i][j]);
			cp[t3id[i]].tc.push_back(t3c[j]);
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t4id[i]].tbf.push_back(t4cid[i][j]);
			cp[t4id[i]].tc.push_back(t4c[j]);
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			cp[t5id[i]].tbf.push_back(t5cid[i][j]);
			cp[t5id[i]].tc.push_back(t5c[j]);
		}
	}
	//for(int i=0;i<2;i++)
	//{
	//	for(int j=0;j<1;j++)
	//	{
	//		cp[t6id[i]].tbf.push_back(t6cid[i][j]);
	//		cp[t6id[i]].tc.push_back(t6c[j]);
	//	}
	//}
	//for(int i=0;i<2;i++)
	//{
	//	for(int j=0;j<2;j++)
	//	{
	//		cp[t7id[i]].tbf.push_back(t7cid[i][j]);
	//		cp[t7id[i]].tc.push_back(t7c[j]);
	//	}
	//}

	//new elements
	int enew[12][4] = {{59, 121, 123, 122}, {121, 60, 124, 123}, {123, 124, 71, 125}, {122, 123, 125, 70}, {60, 126, 127, 124}, {126, 61, 128, 127}, {127, 128, 72, 129}, {124, 127, 129, 71}, {49, 130, 132, 131}, {130, 50, 133, 132}, {132, 133, 61, 126}, {131, 132, 126, 60}};
	int idtmp[13] = {121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133};
	for (int i = 0; i < 12; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		if (i < 4)
			etmp.IEN = tmesh[rid[0]].IEN;
		else if (i >= 4 && i < 8)
			etmp.IEN = tmesh[rid[1]].IEN;
		else
			etmp.IEN = tmesh[rid[2]].IEN;
		for (int j = 0; j < 13; j++)
			etmp.IEN.push_back(idtmp[j]);
		tmesh.push_back(etmp);
	}

	for (int i = 0; i < 3; i++)
		tmesh[rid[i]].act = 0;

	for (int i = 0; i < 100; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (int j = 0; j < 13; j++)
			{
				if (cp[idtmp[j]].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[idtmp[j]].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[idtmp[j]].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[idtmp[j]].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(idtmp[j]);
				}
			}
		}
	}
}

void TruncatedTspline::RefineTest6()
{
	int rid[3] = {54, 55, 45};
	for (int i = 0; i < 13; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int ptf[3] = {123, 127, 132}; //face node
	int fon[3][4] = {{59, 60, 71, 70}, {60, 61, 72, 71}, {49, 50, 61, 60}};
	double f(1. / 4.);
	double fonc[4] = {f, f, f, f};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptf[i]].coor[0] += fonc[j] * cp[fon[i][j]].coor[0];
			cp[ptf[i]].coor[1] += fonc[j] * cp[fon[i][j]].coor[1];
			cp[ptf[i]].coor[2] += fonc[j] * cp[fon[i][j]].coor[2];
		}
	}
	int pte[10] = {121, 122, 124, 125, 126, 128, 129, 130, 131, 133}; //edge node
	int eon[10][6] = {{59, 60, 48, 49, 70, 71}, {59, 70, 58, 69, 60, 71}, {60, 71, 59, 70, 61, 72}, {70, 71, 59, 60, 81, 82}, {60, 61, 49, 50, 71, 72}, {61, 72, 60, 71, 62, 73}, {71, 72, 60, 61, 82, 83}, {49, 50, 38, 39, 60, 61}, {49, 60, 48, 59, 50, 61}, {50, 61, 49, 60, 51, 62}};
	double e1(3. / 8.), e2(1. / 16.);
	double eonc[6] = {e1, e1, e2, e2, e2, e2};
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cp[pte[i]].coor[0] += eonc[j] * cp[eon[i][j]].coor[0];
			cp[pte[i]].coor[1] += eonc[j] * cp[eon[i][j]].coor[1];
			cp[pte[i]].coor[2] += eonc[j] * cp[eon[i][j]].coor[2];
		}
	}

	//update knot vectors
	double ktsU[9] = {cp[58].knotU[2], (cp[58].knotU[2] + cp[59].knotU[2]) / 2., cp[59].knotU[2], (cp[59].knotU[2] + cp[60].knotU[2]) / 2.,
					  cp[60].knotU[2], (cp[60].knotU[2] + cp[61].knotU[2]) / 2., cp[61].knotU[2], (cp[61].knotU[2] + cp[62].knotU[2]) / 2., cp[62].knotU[2]};
	double ktsV[9] = {cp[37].knotV[2], (cp[37].knotV[2] + cp[48].knotV[2]) / 2., cp[48].knotV[2], (cp[48].knotV[2] + cp[59].knotV[2]) / 2.,
					  cp[59].knotV[2], (cp[59].knotV[2] + cp[70].knotV[2]) / 2., cp[70].knotV[2], (cp[70].knotV[2] + cp[81].knotV[2]) / 2., cp[81].knotV[2]};
	int setKU[5][5] = {{-1, -1, -1, 130, -1}, {-1, -1, 131, 132, 133}, {-1, 121, -1, 126, -1}, {122, 123, 124, 127, 128}, {-1, 125, -1, 129, -1}};
	int setKV[5][5] = {{-1, -1, -1, 122, -1}, {-1, -1, 121, 123, 125}, {-1, 131, -1, 124, -1}, {130, 132, 126, 127, 129}, {-1, 133, -1, 128, -1}};
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				if (setKU[i][j] != -1)
					cp[setKU[i][j]].knotU[k] = ktsU[j + k];
			}
		}
	}
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				if (setKV[i][j] != -1)
					cp[setKV[i][j]].knotV[k] = ktsV[j + k];
			}
		}
	}

	int t1id[3] = {70, 72, 50};
	int t1cid[3][5] = {{122, 125, 123, 121, 124}, {128, 129, 127, 126, 124}, {130, 133, 132, 131, 126}};
	double t1c[5] = {e1, e1, f, e2, e2};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			cp[t1id[i]].tbf.push_back(t1cid[i][j]);
			cp[t1id[i]].tc.push_back(t1c[j]);
		}
	}
	int t2id[3] = {48, 62, 82};
	int t2cid[3][2] = {{121, 131}, {133, 128}, {125, 129}};
	double t2c[2] = {e2, e2};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			cp[t2id[i]].tbf.push_back(t2cid[i][j]);
			cp[t2id[i]].tc.push_back(t2c[j]);
		}
	}
	int t3id[2] = {61, 71};
	int t3cid[2][9] = {{133, 128, 126, 132, 127, 130, 131, 124, 129}, {125, 124, 129, 123, 127, 122, 121, 126, 128}};
	double t3c[9] = {e1, e1, e1, f, f, e2, e2, e2, e2};
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 9; j++)
		{
			cp[t3id[i]].tbf.push_back(t3cid[i][j]);
			cp[t3id[i]].tc.push_back(t3c[j]);
		}
	}
	int t4id[1] = {60};
	int t4cid[1][13] = {{123, 127, 132, 121, 124, 126, 131, 122, 125, 129, 128, 133, 130}};
	double t4c[13] = {f, f, f, e1, e1, e1, e1, e2, e2, e2, e2, e2, e2};
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 13; j++)
		{
			cp[t4id[i]].tbf.push_back(t4cid[i][j]);
			cp[t4id[i]].tc.push_back(t4c[j]);
		}
	}
	int t5id[8] = {38, 39, 51, 73, 83, 81, 69, 58};
	int t5cid[8][1] = {{130}, {130}, {133}, {128}, {129}, {125}, {122}, {122}};
	double t5c[1] = {e2};
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t5id[i]].tbf.push_back(t5cid[i][j]);
			cp[t5id[i]].tc.push_back(t5c[j]);
		}
	}
	int t6id[2] = {49, 59};
	int t6cid[2][6] = {{130, 131, 132, 133, 126, 121}, {121, 122, 123, 124, 125, 131}};
	double t6c[6] = {e1, e1, f, e2, e2, e2};
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cp[t6id[i]].tbf.push_back(t6cid[i][j]);
			cp[t6id[i]].tc.push_back(t6c[j]);
		}
	}

	//new elements
	int enew[12][4] = {{59, 121, 123, 122}, {121, 60, 124, 123}, {123, 124, 71, 125}, {122, 123, 125, 70}, {60, 126, 127, 124}, {126, 61, 128, 127}, {127, 128, 72, 129}, {124, 127, 129, 71}, {49, 130, 132, 131}, {130, 50, 133, 132}, {132, 133, 61, 126}, {131, 132, 126, 60}};
	int idtmp[13] = {121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133};
	for (int i = 0; i < 12; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		if (i < 4)
			etmp.IEN = tmesh[rid[0]].IEN;
		else if (i >= 4 && i < 8)
			etmp.IEN = tmesh[rid[1]].IEN;
		else
			etmp.IEN = tmesh[rid[2]].IEN;
		for (int j = 0; j < 13; j++)
			etmp.IEN.push_back(idtmp[j]);
		tmesh.push_back(etmp);
	}

	for (int i = 0; i < 3; i++)
		tmesh[rid[i]].act = 0;

	for (int i = 0; i < 100; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (int j = 0; j < 13; j++)
			{
				if (cp[idtmp[j]].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[idtmp[j]].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[idtmp[j]].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[idtmp[j]].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(idtmp[j]);
				}
			}
		}
	}
}

void TruncatedTspline::RefineTest7()
{
	int rid[3] = {54, 55, 45};
	for (int i = 0; i < 14; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int ptt[7] = {121, 130, 122, 128, 125, 129, 133}; //T-junctions
	int ton[7][4] = {{59, 60, 48, 49}, {49, 50, 38, 39}, {59, 70, 58, 69}, {61, 72, 62, 73}, {70, 71, 81, 82}, {71, 72, 82, 83}, {50, 61, 51, 62}};
	double tonc[4] = {5. / 12., 5. / 12., 1. / 12., 1. / 12.};
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptt[i]].coor[0] += tonc[j] * cp[ton[i][j]].coor[0];
			cp[ptt[i]].coor[1] += tonc[j] * cp[ton[i][j]].coor[1];
			cp[ptt[i]].coor[2] += tonc[j] * cp[ton[i][j]].coor[2];
		}
	}
	int ptf[3] = {123, 127, 132}; //face
	int fon[3][4] = {{59, 60, 71, 70}, {60, 61, 72, 71}, {49, 50, 61, 60}};
	double fonc[4] = {1. / 4., 1. / 4., 1. / 4., 1. / 4.};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptf[i]].coor[0] += fonc[j] * cp[fon[i][j]].coor[0];
			cp[ptf[i]].coor[1] += fonc[j] * cp[fon[i][j]].coor[1];
			cp[ptf[i]].coor[2] += fonc[j] * cp[fon[i][j]].coor[2];
		}
	}
	//int pit[1]={134};//inserted T-junction
	//int iton[1][2]={48,49};
	//double itonc[2]={.5,.5};
	//for(int i=0;i<1;i++)
	//{
	//	for(int j=0;j<2;j++)
	//	{
	//		cp[pit[i]].coor[0]+=itonc[j]*cp[iton[i][j]].coor[0];
	//		cp[pit[i]].coor[1]+=itonc[j]*cp[iton[i][j]].coor[1];
	//		cp[pit[i]].coor[2]+=itonc[j]*cp[iton[i][j]].coor[2];
	//	}
	//}
	int pte[3] = {124, 126, 131}; //edge
	int eon[3][6] = {{60, 71, 59, 70, 61, 72}, {60, 61, 49, 50, 71, 72}, {49, 60, 48, 59, 50, 61}};
	double eonc[6] = {3. / 8., 3. / 8., 1. / 16., 1. / 16., 1. / 16., 1. / 16.};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cp[pte[i]].coor[0] += eonc[j] * cp[eon[i][j]].coor[0];
			cp[pte[i]].coor[1] += eonc[j] * cp[eon[i][j]].coor[1];
			cp[pte[i]].coor[2] += eonc[j] * cp[eon[i][j]].coor[2];
		}
	}

	cp[134].coor[0] = (cp[48].coor[0] + cp[49].coor[0]) / 2.;
	cp[134].coor[1] = (cp[48].coor[1] + cp[49].coor[1]) / 2.;
	cp[134].coor[2] = (cp[48].coor[2] + cp[49].coor[2]) / 2.;

	int ptcc[1] = {60}; //corner corner
	int ccon[1][9] = {60, 49, 59, 61, 71, 48, 50, 70, 72};
	double cconc[9] = {9. / 16., 3. / 32., 3. / 32., 3. / 32., 3. / 32., 1. / 64., 1. / 64., 1. / 64., 1. / 64.};
	for (int i = 0; i < 1; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 9; j++)
		{
			tmp.coor[0] += cconc[j] * cp[ccon[i][j]].coor[0];
			tmp.coor[1] += cconc[j] * cp[ccon[i][j]].coor[1];
			tmp.coor[2] += cconc[j] * cp[ccon[i][j]].coor[2];
		}
		cp[ptcc[i]].coor[0] = tmp.coor[0];
		cp[ptcc[i]].coor[1] = tmp.coor[1];
		cp[ptcc[i]].coor[2] = tmp.coor[2];
	}
	int ptec[3] = {61, 71, 49}; //edge corner
	int econ[3][6] = {{61, 62, 50, 51, 72, 73}, {71, 82, 70, 81, 72, 83}, {49, 38, 48, 37, 50, 39}};
	double econc[6] = {15. / 24., 3. / 24., 5. / 48., 1. / 48., 5. / 48., 1. / 48.};
	for (int i = 0; i < 3; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 6; j++)
		{
			tmp.coor[0] += econc[j] * cp[econ[i][j]].coor[0];
			tmp.coor[1] += econc[j] * cp[econ[i][j]].coor[1];
			tmp.coor[2] += econc[j] * cp[econ[i][j]].coor[2];
		}
		cp[ptec[i]].coor[0] = tmp.coor[0];
		cp[ptec[i]].coor[1] = tmp.coor[1];
		cp[ptec[i]].coor[2] = tmp.coor[2];
	}
	int ptc[4] = {59, 70, 72, 50}; //corner
	int con[4][4] = {{59, 48, 58, 47}, {70, 69, 81, 80}, {72, 73, 83, 84}, {50, 39, 51, 40}};
	double conc[4] = {25. / 36., 5. / 36., 5. / 36., 1. / 36.};
	for (int i = 0; i < 4; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 4; j++)
		{
			tmp.coor[0] += conc[j] * cp[con[i][j]].coor[0];
			tmp.coor[1] += conc[j] * cp[con[i][j]].coor[1];
			tmp.coor[2] += conc[j] * cp[con[i][j]].coor[2];
		}
		cp[ptc[i]].coor[0] = tmp.coor[0];
		cp[ptc[i]].coor[1] = tmp.coor[1];
		cp[ptc[i]].coor[2] = tmp.coor[2];
	}
	double vtmp[3] = {(cp[47].coor[0] + 5. * cp[48].coor[0]) / 6., (cp[47].coor[1] + 5. * cp[48].coor[1]) / 6., (cp[47].coor[2] + 5. * cp[48].coor[2]) / 6.};
	cp[48].coor[0] = vtmp[0];
	cp[48].coor[1] = vtmp[1];
	cp[48].coor[2] = vtmp[2];

	//update knot vectors
	double ktsU[9] = {cp[59].knotU[0], cp[59].knotU[1], cp[59].knotU[2], (cp[59].knotU[2] + cp[59].knotU[3]) / 2., cp[59].knotU[3], (cp[59].knotU[3] + cp[59].knotU[4]) / 2., cp[59].knotU[4], cp[61].knotU[3], cp[61].knotU[4]};
	double ktsV[7] = {cp[59].knotV[0], cp[59].knotV[1], cp[59].knotV[2], (cp[59].knotV[2] + cp[59].knotV[3]) / 2., cp[59].knotV[3], cp[59].knotV[4], cp[70].knotV[4]};
	double ktsU1[7] = {cp[49].knotU[1], (cp[49].knotU[1] + cp[49].knotU[2]) / 2., cp[49].knotU[2], (cp[49].knotU[2] + cp[49].knotU[3]) / 2., cp[49].knotU[3], cp[49].knotU[4], cp[50].knotU[4]};
	double ktsV1[9] = {cp[49].knotV[0], cp[49].knotV[1], cp[49].knotV[2], (cp[49].knotV[2] + cp[49].knotV[3]) / 2., cp[49].knotV[3], (cp[49].knotV[3] + cp[49].knotV[4]) / 2., cp[49].knotV[4], cp[71].knotV[3], cp[71].knotV[4]};
	int setKU[3][5] = {{59, 121, 60, 126, 61}, {122, 123, 124, 127, 128}, {70, 125, 71, 129, 72}};
	int setKV[2][3] = {{59, 122, 70}, {121, 123, 125}};
	int setKU1[2][3] = {{49, 130, 50}, {131, 132, 133}};
	int setKV1[3][5] = {{49, 131, 60, 124, 71}, {130, 132, 126, 127, 129}, {50, 133, 61, 128, 72}};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKU[i][j]].knotU[k] = ktsU[j + k];
			}
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKV[i][j]].knotV[k] = ktsV[j + k];
			}
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKU1[i][j]].knotU[k] = ktsU1[j + k];
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKV1[i][j]].knotV[k] = ktsV1[j + k];
			}
		}
	}

	//cp[38].knotV[4]=cp[131].knotV[2];
	cp[39].knotV[4] = cp[133].knotV[2];
	cp[47].knotU[4] = cp[121].knotU[2];
	cp[48].knotU[3] = cp[121].knotU[2];
	cp[48].knotU[4] = cp[131].knotU[2];
	cp[48].knotV[4] = cp[122].knotV[2];
	cp[51].knotU[0] = cp[130].knotU[2];
	cp[58].knotU[4] = cp[121].knotU[2];
	cp[62].knotU[0] = cp[126].knotU[2];
	cp[69].knotU[4] = cp[125].knotU[2];
	cp[73].knotU[0] = cp[129].knotU[2];
	cp[81].knotV[0] = cp[122].knotV[2];
	cp[82].knotV[0] = cp[124].knotV[2];
	cp[83].knotV[0] = cp[128].knotV[2];

	for (int i = 0; i < 5; i++)
	{
		cp[134].knotU[i] = cp[121].knotU[i];
		cp[134].knotV[i] = cp[48].knotV[i];
	}

	//int t1id[4]={48,50,81,83};
	//int t1cid[4][3]={{59,121,60},{61,126,60},{70,125,71},{72,129,71}};
	//double t1c[3]={5./36.,1./12.,1./48.};
	int t2id[3] = {40, 80, 84};
	int t2cid[3][1] = {{50}, {70}, {72}};
	double t2c[1] = {1. / 36.};
	//int t3id[2]={49,82};
	//int t3cid[2][3]={{60,121,126},{71,125,129}};
	//double t3c[3]={3./24.,1./12.,1./12.};
	//int t4id[4]={58,62,69,73};
	//int t4cid[4][2]={{59,122},{61,128},{70,122},{72,128}};
	//double t4c[2]={5./36.,1./12.};
	int t4id[1] = {37};
	int t4cid[1][1] = {{49}};
	double t4c[1] = {1. / 48.};
	int t5id[1] = {38};
	int t5cid[1][2] = {{49, 130}};
	double t5c[2] = {3. / 24., 1. / 12.};
	//for(int i=0;i<4;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t1id[i]].tbf.push_back(t1cid[i][j]);
	//		cp[t1id[i]].tc.push_back(t1c[j]);
	//	}
	//}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t2id[i]].tbf.push_back(t2cid[i][j]);
			cp[t2id[i]].tc.push_back(t2c[j]);
		}
	}
	//for(int i=0;i<2;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t3id[i]].tbf.push_back(t3cid[i][j]);
	//		cp[t3id[i]].tc.push_back(t3c[j]);
	//	}
	//}
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t4id[i]].tbf.push_back(t4cid[i][j]);
			cp[t4id[i]].tc.push_back(t4c[j]);
		}
	}
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			cp[t5id[i]].tbf.push_back(t5cid[i][j]);
			cp[t5id[i]].tc.push_back(t5c[j]);
		}
	}

	//new elements
	int enew[14][4] = {{59, 121, 123, 122}, {121, 60, 124, 123}, {123, 124, 71, 125}, {122, 123, 125, 70}, {60, 126, 127, 124}, {126, 61, 128, 127}, {127, 128, 72, 129}, {124, 127, 129, 71}, {49, 130, 132, 131}, {130, 50, 133, 132}, {132, 133, 61, 126}, {131, 132, 126, 60}, {48, 134, 121, 59}, {134, 49, 60, 121}};
	int idtmp[14] = {121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134};
	for (int i = 0; i < 14; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		if (i < 4)
			etmp.IEN = tmesh[rid[0]].IEN;
		else if (i >= 4 && i < 8)
			etmp.IEN = tmesh[rid[1]].IEN;
		else if (i >= 8 && i < 12)
			etmp.IEN = tmesh[rid[2]].IEN;
		else
			etmp.IEN = tmesh[44].IEN;
		for (int j = 0; j < 14; j++)
			etmp.IEN.push_back(idtmp[j]);
		tmesh.push_back(etmp);
	}

	for (int i = 0; i < 3; i++)
		tmesh[rid[i]].act = 0;
	tmesh[44].act = 0;

	for (int i = 0; i < 100; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (int j = 0; j < 14; j++)
			{
				if (cp[idtmp[j]].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[idtmp[j]].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[idtmp[j]].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[idtmp[j]].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(idtmp[j]);
				}
			}
		}
	}
}

void TruncatedTspline::RefineTest8()
{
	int rid[2] = {54, 55};
	for (int i = 0; i < 10; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int ptt[6] = {121, 126, 122, 128, 125, 129}; //T-junctions
	int ton[6][4] = {{59, 60, 48, 49}, {60, 61, 49, 50}, {59, 70, 58, 69}, {61, 72, 62, 73}, {70, 71, 81, 82}, {71, 72, 82, 83}};
	double tonc[4] = {5. / 12., 5. / 12., 1. / 12., 1. / 12.};
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptt[i]].coor[0] += tonc[j] * cp[ton[i][j]].coor[0];
			cp[ptt[i]].coor[1] += tonc[j] * cp[ton[i][j]].coor[1];
			cp[ptt[i]].coor[2] += tonc[j] * cp[ton[i][j]].coor[2];
		}
	}
	int ptf[2] = {123, 127}; //face
	int fon[2][4] = {{59, 60, 71, 70}, {60, 61, 72, 71}};
	double fonc[4] = {1. / 4., 1. / 4., 1. / 4., 1. / 4.};
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptf[i]].coor[0] += fonc[j] * cp[fon[i][j]].coor[0];
			cp[ptf[i]].coor[1] += fonc[j] * cp[fon[i][j]].coor[1];
			cp[ptf[i]].coor[2] += fonc[j] * cp[fon[i][j]].coor[2];
		}
	}
	int pte[1] = {124}; //edge
	int eon[1][6] = {{60, 71, 59, 70, 61, 72}};
	double eonc[6] = {3. / 8., 3. / 8., 1. / 16., 1. / 16., 1. / 16., 1. / 16.};
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cp[pte[i]].coor[0] += eonc[j] * cp[eon[i][j]].coor[0];
			cp[pte[i]].coor[1] += eonc[j] * cp[eon[i][j]].coor[1];
			cp[pte[i]].coor[2] += eonc[j] * cp[eon[i][j]].coor[2];
		}
	}
	//int ptcc[1]={60};//corner corner
	//int ccon[1][9]={60,49,59,61,71,48,50,70,72};
	//double cconc[9]={9./16.,3./32.,3./32.,3./32.,3./32.,1./64.,1./64.,1./64.,1./64.};
	//for(int i=0;i<1;i++)
	//{
	//	Vertex tmp;
	//	for(int j=0;j<9;j++)
	//	{
	//		tmp.coor[0]+=cconc[j]*cp[ccon[i][j]].coor[0];
	//		tmp.coor[1]+=cconc[j]*cp[ccon[i][j]].coor[1];
	//		tmp.coor[2]+=cconc[j]*cp[ccon[i][j]].coor[2];
	//	}
	//	cp[ptcc[i]].coor[0]=tmp.coor[0]; cp[ptcc[i]].coor[1]=tmp.coor[1]; cp[ptcc[i]].coor[2]=tmp.coor[2];
	//}
	int ptec[2] = {60, 71}; //edge corner
	int econ[2][6] = {{60, 49, 59, 48, 61, 50}, {71, 82, 70, 81, 72, 83}};
	double econc[6] = {15. / 24., 3. / 24., 5. / 48., 1. / 48., 5. / 48., 1. / 48.};
	for (int i = 0; i < 2; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 6; j++)
		{
			tmp.coor[0] += econc[j] * cp[econ[i][j]].coor[0];
			tmp.coor[1] += econc[j] * cp[econ[i][j]].coor[1];
			tmp.coor[2] += econc[j] * cp[econ[i][j]].coor[2];
		}
		cp[ptec[i]].coor[0] = tmp.coor[0];
		cp[ptec[i]].coor[1] = tmp.coor[1];
		cp[ptec[i]].coor[2] = tmp.coor[2];
	}
	int ptc[4] = {59, 61, 70, 72}; //corner
	int con[4][4] = {{59, 48, 58, 47}, {61, 50, 62, 51}, {70, 69, 81, 80}, {72, 73, 83, 84}};
	double conc[4] = {25. / 36., 5. / 36., 5. / 36., 1. / 36.};
	for (int i = 0; i < 4; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 4; j++)
		{
			tmp.coor[0] += conc[j] * cp[con[i][j]].coor[0];
			tmp.coor[1] += conc[j] * cp[con[i][j]].coor[1];
			tmp.coor[2] += conc[j] * cp[con[i][j]].coor[2];
		}
		cp[ptc[i]].coor[0] = tmp.coor[0];
		cp[ptc[i]].coor[1] = tmp.coor[1];
		cp[ptc[i]].coor[2] = tmp.coor[2];
	}

	cp[130].coor[0] = (cp[48].coor[0] + cp[49].coor[0]) / 2.;
	cp[130].coor[1] = (cp[48].coor[1] + cp[49].coor[1]) / 2.;
	cp[130].coor[2] = (cp[48].coor[2] + cp[49].coor[2]) / 2.;
	double vtmp[3] = {(cp[47].coor[0] + 5. * cp[48].coor[0]) / 6., (cp[47].coor[1] + 5. * cp[48].coor[1]) / 6., (cp[47].coor[2] + 5. * cp[48].coor[2]) / 6.};
	cp[48].coor[0] = vtmp[0];
	cp[48].coor[1] = vtmp[1];
	cp[48].coor[2] = vtmp[2];
	double vtmp1[3] = {(cp[50].coor[0] + 5. * cp[49].coor[0]) / 6., (cp[50].coor[1] + 5. * cp[49].coor[1]) / 6., (cp[50].coor[2] + 5. * cp[49].coor[2]) / 6.};
	cp[49].coor[0] = vtmp1[0];
	cp[49].coor[1] = vtmp1[1];
	cp[49].coor[2] = vtmp1[2];

	//update knot vectors
	double ktsU[9] = {cp[59].knotU[0], cp[59].knotU[1], cp[59].knotU[2], (cp[59].knotU[2] + cp[59].knotU[3]) / 2., cp[59].knotU[3], (cp[59].knotU[3] + cp[59].knotU[4]) / 2., cp[59].knotU[4], cp[61].knotU[3], cp[61].knotU[4]};
	double ktsV[7] = {cp[59].knotV[0], cp[59].knotV[1], cp[59].knotV[2], (cp[59].knotV[2] + cp[59].knotV[3]) / 2., cp[59].knotV[3], cp[59].knotV[4], cp[70].knotV[4]};
	//double ktsU1[7]={cp[49].knotU[0],cp[49].knotU[1],cp[49].knotU[2],(cp[49].knotU[2]+cp[49].knotU[3])/2.,cp[49].knotU[3],cp[49].knotU[4],cp[50].knotU[4]};
	//double ktsV1[9]={cp[49].knotV[0],cp[49].knotV[1],cp[49].knotV[2],(cp[49].knotV[2]+cp[49].knotV[3])/2.,cp[49].knotV[3],(cp[49].knotV[3]+cp[49].knotV[4])/2.,cp[49].knotV[4],cp[71].knotV[3],cp[71].knotV[4]};
	int setKU[3][5] = {{59, 121, 60, 126, 61}, {122, 123, 124, 127, 128}, {70, 125, 71, 129, 72}};
	int setKV[5][3] = {{59, 122, 70}, {121, 123, 125}, {60, 124, 71}, {126, 127, 129}, {61, 128, 72}};
	//int setKU1[2][3]={{49,130,50},{131,132,133}};
	//int setKV1[3][5]={{49,131,60,124,71},{130,132,126,127,129},{50,133,61,128,72}};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKU[i][j]].knotU[k] = ktsU[j + k];
			}
		}
	}
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[setKV[i][j]].knotV[k] = ktsV[j + k];
			}
		}
	}
	//for(int i=0;i<2;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		for(int k=0;k<5;k++)
	//		{
	//			cp[setKU1[i][j]].knotU[k]=ktsU1[j+k];
	//		}
	//	}
	//}
	//for(int i=0;i<3;i++)
	//{
	//	for(int j=0;j<5;j++)
	//	{
	//		for(int k=0;k<5;k++)
	//		{
	//			cp[setKV1[i][j]].knotV[k]=ktsV1[j+k];
	//		}
	//	}
	//}

	cp[48].knotV[4] = cp[122].knotV[2];
	cp[49].knotV[4] = cp[124].knotV[2];
	cp[50].knotV[4] = cp[128].knotV[2];
	cp[58].knotU[4] = cp[121].knotU[2];
	cp[62].knotU[0] = cp[126].knotU[2];
	cp[69].knotU[4] = cp[125].knotU[2];
	cp[73].knotU[0] = cp[129].knotU[2];
	cp[81].knotV[0] = cp[122].knotV[2];
	cp[82].knotV[0] = cp[124].knotV[2];
	cp[83].knotV[0] = cp[128].knotV[2];

	for (int i = 0; i < 5; i++)
	{
		cp[130].knotU[i] = cp[121].knotU[i];
		cp[130].knotV[i] = cp[48].knotV[i];
	}
	cp[130].knotU[4] = cp[50].knotU[2];
	cp[47].knotU[4] = cp[130].knotU[2];
	cp[48].knotU[3] = cp[130].knotU[2];
	cp[48].knotU[4] = cp[49].knotU[2];
	cp[49].knotU[0] = cp[48].knotU[2];
	cp[49].knotU[1] = cp[130].knotU[2];
	cp[50].knotU[0] = cp[130].knotU[2];

	//int t1id[4]={48,50,81,83};
	//int t1cid[4][3]={{59,121,60},{61,126,60},{70,125,71},{72,129,71}};
	//double t1c[3]={5./36.,1./12.,1./48.};
	int t2id[3] = {51, 80, 84};
	int t2cid[3][1] = {{61}, {70}, {72}};
	double t2c[1] = {1. / 36.};
	//int t3id[2]={49,82};
	//int t3cid[2][3]={{60,121,126},{71,125,129}};
	//double t3c[3]={3./24.,1./12.,1./12.};
	//int t4id[4]={58,62,69,73};
	//int t4cid[4][2]={{59,122},{61,128},{70,122},{72,128}};
	//double t4c[2]={5./36.,1./12.};
	//int t5id[1]={48};
	//int t5cid[1][3]={{60,121,131}};
	//double t5c[3]={1./64.,1./12.,1./12.};
	//for(int i=0;i<4;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t1id[i]].tbf.push_back(t1cid[i][j]);
	//		cp[t1id[i]].tc.push_back(t1c[j]);
	//	}
	//}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t2id[i]].tbf.push_back(t2cid[i][j]);
			cp[t2id[i]].tc.push_back(t2c[j]);
		}
	}
	//for(int i=0;i<2;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t3id[i]].tbf.push_back(t3cid[i][j]);
	//		cp[t3id[i]].tc.push_back(t3c[j]);
	//	}
	//}
	//for(int i=0;i<4;i++)
	//{
	//	for(int j=0;j<2;j++)
	//	{
	//		cp[t4id[i]].tbf.push_back(t4cid[i][j]);
	//		cp[t4id[i]].tc.push_back(t4c[j]);
	//	}
	//}
	//for(int i=0;i<1;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t5id[i]].tbf.push_back(t5cid[i][j]);
	//		cp[t5id[i]].tc.push_back(t5c[j]);
	//	}
	//}

	//new elements
	int enew[10][4] = {{59, 121, 123, 122}, {121, 60, 124, 123}, {123, 124, 71, 125}, {122, 123, 125, 70}, {60, 126, 127, 124}, {126, 61, 128, 127}, {127, 128, 72, 129}, {124, 127, 129, 71}, {48, 130, 121, 59}, {130, 49, 60, 121}};
	int idtmp[10] = {121, 122, 123, 124, 125, 126, 127, 128, 129, 130};
	for (int i = 0; i < 10; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		if (i < 4)
			etmp.IEN = tmesh[rid[0]].IEN;
		else if (i >= 4 && i < 8)
			etmp.IEN = tmesh[rid[1]].IEN;
		else
			etmp.IEN = tmesh[44].IEN;
		for (int j = 0; j < 10; j++)
			etmp.IEN.push_back(idtmp[j]);
		tmesh.push_back(etmp);
	}

	for (int i = 0; i < 2; i++)
		tmesh[rid[i]].act = 0;
	tmesh[44].act = 0;

	for (int i = 0; i < 100; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (int j = 0; j < 10; j++)
			{
				if (cp[idtmp[j]].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[idtmp[j]].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[idtmp[j]].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[idtmp[j]].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(idtmp[j]);
				}
			}
		}
	}
}

void TruncatedTspline::RefineTest9() //not complete yet, stop at change knot vector
{
	int rid[3] = {54, 55, 45};
	for (int i = 0; i < 13; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int ptt[6] = {130, 122, 128, 125, 129, 133}; //T-junctions
	int ton[6][4] = {{49, 50, 38, 39}, {59, 70, 58, 69}, {61, 72, 62, 73}, {70, 71, 81, 82}, {71, 72, 82, 83}, {50, 61, 51, 62}};
	double tonc[4] = {5. / 12., 5. / 12., 1. / 12., 1. / 12.};
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptt[i]].coor[0] += tonc[j] * cp[ton[i][j]].coor[0];
			cp[ptt[i]].coor[1] += tonc[j] * cp[ton[i][j]].coor[1];
			cp[ptt[i]].coor[2] += tonc[j] * cp[ton[i][j]].coor[2];
		}
	}
	int ptf[3] = {123, 127, 132}; //face
	int fon[3][4] = {{59, 60, 71, 70}, {60, 61, 72, 71}, {49, 50, 61, 60}};
	double fonc[4] = {1. / 4., 1. / 4., 1. / 4., 1. / 4.};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptf[i]].coor[0] += fonc[j] * cp[fon[i][j]].coor[0];
			cp[ptf[i]].coor[1] += fonc[j] * cp[fon[i][j]].coor[1];
			cp[ptf[i]].coor[2] += fonc[j] * cp[fon[i][j]].coor[2];
		}
	}
	int pte[4] = {121, 124, 126, 131}; //edge
	int eon[4][6] = {{59, 60, 48, 49, 70, 71}, {60, 71, 59, 70, 61, 72}, {60, 61, 49, 50, 71, 72}, {49, 60, 48, 59, 50, 61}};
	double eonc[6] = {3. / 8., 3. / 8., 1. / 16., 1. / 16., 1. / 16., 1. / 16.};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cp[pte[i]].coor[0] += eonc[j] * cp[eon[i][j]].coor[0];
			cp[pte[i]].coor[1] += eonc[j] * cp[eon[i][j]].coor[1];
			cp[pte[i]].coor[2] += eonc[j] * cp[eon[i][j]].coor[2];
		}
	}
	//int ptcc[1]={60};//corner corner
	//int ccon[1][9]={60,49,59,61,71,48,50,70,72};
	//double cconc[9]={9./16.,3./32.,3./32.,3./32.,3./32.,1./64.,1./64.,1./64.,1./64.};
	//for(int i=0;i<1;i++)
	//{
	//	Vertex tmp;
	//	for(int j=0;j<9;j++)
	//	{
	//		tmp.coor[0]+=cconc[j]*cp[ccon[i][j]].coor[0];
	//		tmp.coor[1]+=cconc[j]*cp[ccon[i][j]].coor[1];
	//		tmp.coor[2]+=cconc[j]*cp[ccon[i][j]].coor[2];
	//	}
	//	cp[ptcc[i]].coor[0]=tmp.coor[0]; cp[ptcc[i]].coor[1]=tmp.coor[1]; cp[ptcc[i]].coor[2]=tmp.coor[2];
	//}
	int ptec[2] = {61, 71}; //edge corner
	int econ[2][6] = {{61, 62, 50, 51, 72, 73}, {71, 82, 70, 81, 72, 83}};
	double econc[6] = {15. / 24., 3. / 24., 5. / 48., 1. / 48., 5. / 48., 1. / 48.};
	for (int i = 0; i < 2; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 6; j++)
		{
			tmp.coor[0] += econc[j] * cp[econ[i][j]].coor[0];
			tmp.coor[1] += econc[j] * cp[econ[i][j]].coor[1];
			tmp.coor[2] += econc[j] * cp[econ[i][j]].coor[2];
		}
		cp[ptec[i]].coor[0] = tmp.coor[0];
		cp[ptec[i]].coor[1] = tmp.coor[1];
		cp[ptec[i]].coor[2] = tmp.coor[2];
	}
	int ptc[3] = {70, 72, 50}; //corner
	int con[3][4] = {{70, 69, 81, 80}, {72, 73, 83, 84}, {50, 39, 51, 40}};
	double conc[4] = {25. / 36., 5. / 36., 5. / 36., 1. / 36.};
	for (int i = 0; i < 3; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 4; j++)
		{
			tmp.coor[0] += conc[j] * cp[con[i][j]].coor[0];
			tmp.coor[1] += conc[j] * cp[con[i][j]].coor[1];
			tmp.coor[2] += conc[j] * cp[con[i][j]].coor[2];
		}
		cp[ptc[i]].coor[0] = tmp.coor[0];
		cp[ptc[i]].coor[1] = tmp.coor[1];
		cp[ptc[i]].coor[2] = tmp.coor[2];
	}

	//update knot vectors
	double ktsU[9] = {cp[59].knotU[0], cp[59].knotU[1], cp[59].knotU[2], (cp[59].knotU[2] + cp[59].knotU[3]) / 2., cp[59].knotU[3], (cp[59].knotU[3] + cp[59].knotU[4]) / 2., cp[59].knotU[4], cp[61].knotU[3], cp[61].knotU[4]};
	//double ktsV[7]={cp[59].knotV[1],(cp[59].knotV[1]+cp[59].knotV[2])/2.,cp[59].knotV[2],(cp[59].knotV[2]+cp[59].knotV[3])/2.,cp[59].knotV[3],cp[59].knotV[4],cp[70].knotV[4]};
	//double ktsU1[7]={cp[49].knotU[1],(cp[49].knotU[1]+cp[49].knotU[2])/2.,cp[49].knotU[2],(cp[49].knotU[2]+cp[49].knotU[3])/2.,cp[49].knotU[3],cp[49].knotU[4],cp[50].knotU[4]};
	double ktsV[9] = {cp[49].knotV[0], cp[49].knotV[1], cp[49].knotV[2], (cp[49].knotV[2] + cp[49].knotV[3]) / 2., cp[49].knotV[3], (cp[49].knotV[3] + cp[49].knotV[4]) / 2., cp[49].knotV[4], cp[71].knotV[3], cp[71].knotV[4]};
	int setKU[5][5] = {{-1, -1, -1, 130, 50}, {-1, -1, 131, 132, 133}, {-1, 121, -1, 126, 61}, {122, 123, 124, 127, 128}, {70, 125, 71, 129, 72}};
	int setKV[5][5] = {{-1, -1, -1, 122, 70}, {-1, -1, 121, 123, 125}, {-1, 131, -1, 124, 71}, {130, 132, 126, 127, 129}, {50, 133, 61, 128, 72}};
	//int setKU1[2][3]={{49,130,50},{131,132,133}};
	//int setKV1[3][5]={{49,131,60,124,71},{130,132,126,127,129},{50,133,61,128,72}};
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			if (setKU[i][j] != 0)
			{
				for (int k = 0; k < 5; k++)
				{
					cp[setKU[i][j]].knotU[k] = ktsU[j + k];
				}
			}
		}
	}
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			if (setKV[i][j] != 0)
			{
				for (int k = 0; k < 5; k++)
				{
					cp[setKV[i][j]].knotV[k] = ktsV[j + k];
				}
			}
		}
	}
	//for(int i=0;i<2;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		for(int k=0;k<5;k++)
	//		{
	//			cp[setKU1[i][j]].knotU[k]=ktsU1[j+k];
	//		}
	//	}
	//}
	//for(int i=0;i<3;i++)
	//{
	//	for(int j=0;j<5;j++)
	//	{
	//		for(int k=0;k<5;k++)
	//		{
	//			cp[setKV1[i][j]].knotV[k]=ktsV1[j+k];
	//		}
	//	}
	//}

	//cp[38].knotV[4]=cp[131].knotV[2];
	cp[39].knotV[4] = cp[133].knotV[2];
	cp[47].knotU[4] = cp[121].knotU[2];
	cp[48].knotU[3] = cp[121].knotU[2];
	cp[48].knotU[4] = cp[131].knotU[2];
	cp[48].knotV[4] = cp[122].knotV[2];
	cp[51].knotU[0] = cp[130].knotU[2];
	cp[58].knotU[4] = cp[121].knotU[2];
	cp[62].knotU[0] = cp[126].knotU[2];
	cp[69].knotU[4] = cp[125].knotU[2];
	cp[73].knotU[0] = cp[129].knotU[2];
	cp[81].knotV[0] = cp[122].knotV[2];
	cp[82].knotV[0] = cp[124].knotV[2];
	cp[83].knotV[0] = cp[128].knotV[2];

	//int t1id[4]={48,50,81,83};
	//int t1cid[4][3]={{59,121,60},{61,126,60},{70,125,71},{72,129,71}};
	//double t1c[3]={5./36.,1./12.,1./48.};
	int t2id[3] = {40, 80, 84};
	int t2cid[3][1] = {{50}, {70}, {72}};
	double t2c[1] = {1. / 36.};
	//int t3id[2]={49,82};
	//int t3cid[2][3]={{60,121,126},{71,125,129}};
	//double t3c[3]={3./24.,1./12.,1./12.};
	//int t4id[4]={58,62,69,73};
	//int t4cid[4][2]={{59,122},{61,128},{70,122},{72,128}};
	//double t4c[2]={5./36.,1./12.};
	int t4id[1] = {37};
	int t4cid[1][1] = {{49}};
	double t4c[1] = {1. / 48.};
	int t5id[1] = {38};
	int t5cid[1][2] = {{49, 130}};
	double t5c[2] = {3. / 24., 1. / 12.};
	//for(int i=0;i<4;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t1id[i]].tbf.push_back(t1cid[i][j]);
	//		cp[t1id[i]].tc.push_back(t1c[j]);
	//	}
	//}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t2id[i]].tbf.push_back(t2cid[i][j]);
			cp[t2id[i]].tc.push_back(t2c[j]);
		}
	}
	//for(int i=0;i<2;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t3id[i]].tbf.push_back(t3cid[i][j]);
	//		cp[t3id[i]].tc.push_back(t3c[j]);
	//	}
	//}
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t4id[i]].tbf.push_back(t4cid[i][j]);
			cp[t4id[i]].tc.push_back(t4c[j]);
		}
	}
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			cp[t5id[i]].tbf.push_back(t5cid[i][j]);
			cp[t5id[i]].tc.push_back(t5c[j]);
		}
	}

	//new elements
	int enew[14][4] = {{59, 121, 123, 122}, {121, 60, 124, 123}, {123, 124, 71, 125}, {122, 123, 125, 70}, {60, 126, 127, 124}, {126, 61, 128, 127}, {127, 128, 72, 129}, {124, 127, 129, 71}, {49, 130, 132, 131}, {130, 50, 133, 132}, {132, 133, 61, 126}, {131, 132, 126, 60}, {48, 134, 121, 59}, {134, 49, 60, 121}};
	int idtmp[14] = {121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134};
	for (int i = 0; i < 14; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		if (i < 4)
			etmp.IEN = tmesh[rid[0]].IEN;
		else if (i >= 4 && i < 8)
			etmp.IEN = tmesh[rid[1]].IEN;
		else if (i >= 8 && i < 12)
			etmp.IEN = tmesh[rid[2]].IEN;
		else
			etmp.IEN = tmesh[44].IEN;
		for (int j = 0; j < 14; j++)
			etmp.IEN.push_back(idtmp[j]);
		tmesh.push_back(etmp);
	}

	for (int i = 0; i < 3; i++)
		tmesh[rid[i]].act = 0;
	tmesh[44].act = 0;

	for (int i = 0; i < 100; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (int j = 0; j < 14; j++)
			{
				if (cp[idtmp[j]].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[idtmp[j]].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[idtmp[j]].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[idtmp[j]].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(idtmp[j]);
				}
			}
		}
	}
}

void TruncatedTspline::RefineTest10()
{
	int rid[3] = {54, 55, 45};
	for (int i = 0; i < 15; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int ptt[7] = {134, 130, 122, 128, 125, 129, 133}; //T-junctions
	int ton[7][4] = {{48, 49, 37, 38}, {49, 50, 38, 39}, {59, 70, 58, 69}, {61, 72, 62, 73}, {70, 71, 81, 82}, {71, 72, 82, 83}, {50, 61, 51, 62}};
	double tonc[4] = {5. / 12., 5. / 12., 1. / 12., 1. / 12.};
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptt[i]].coor[0] += tonc[j] * cp[ton[i][j]].coor[0];
			cp[ptt[i]].coor[1] += tonc[j] * cp[ton[i][j]].coor[1];
			cp[ptt[i]].coor[2] += tonc[j] * cp[ton[i][j]].coor[2];
		}
	}
	int ptf[4] = {123, 127, 132, 135}; //face
	int fon[4][4] = {{59, 60, 71, 70}, {60, 61, 72, 71}, {49, 50, 61, 60}, {48, 49, 60, 59}};
	double fonc[4] = {1. / 4., 1. / 4., 1. / 4., 1. / 4.};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptf[i]].coor[0] += fonc[j] * cp[fon[i][j]].coor[0];
			cp[ptf[i]].coor[1] += fonc[j] * cp[fon[i][j]].coor[1];
			cp[ptf[i]].coor[2] += fonc[j] * cp[fon[i][j]].coor[2];
		}
	}
	//int pit[1]={134};//inserted T-junction
	//int iton[1][2]={48,49};
	//double itonc[2]={.5,.5};
	//for(int i=0;i<1;i++)
	//{
	//	for(int j=0;j<2;j++)
	//	{
	//		cp[pit[i]].coor[0]+=itonc[j]*cp[iton[i][j]].coor[0];
	//		cp[pit[i]].coor[1]+=itonc[j]*cp[iton[i][j]].coor[1];
	//		cp[pit[i]].coor[2]+=itonc[j]*cp[iton[i][j]].coor[2];
	//	}
	//}
	int pte[4] = {124, 126, 131, 121}; //edge
	int eon[4][6] = {{60, 71, 59, 70, 61, 72}, {60, 61, 49, 50, 71, 72}, {49, 60, 48, 59, 50, 61}, {59, 60, 48, 49, 70, 71}};
	double eonc[6] = {3. / 8., 3. / 8., 1. / 16., 1. / 16., 1. / 16., 1. / 16.};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cp[pte[i]].coor[0] += eonc[j] * cp[eon[i][j]].coor[0];
			cp[pte[i]].coor[1] += eonc[j] * cp[eon[i][j]].coor[1];
			cp[pte[i]].coor[2] += eonc[j] * cp[eon[i][j]].coor[2];
		}
	}

	int ptcc[1] = {60}; //corner corner
	int ccon[1][9] = {60, 49, 59, 61, 71, 48, 50, 70, 72};
	double cconc[9] = {9. / 16., 3. / 32., 3. / 32., 3. / 32., 3. / 32., 1. / 64., 1. / 64., 1. / 64., 1. / 64.};
	for (int i = 0; i < 1; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 9; j++)
		{
			tmp.coor[0] += cconc[j] * cp[ccon[i][j]].coor[0];
			tmp.coor[1] += cconc[j] * cp[ccon[i][j]].coor[1];
			tmp.coor[2] += cconc[j] * cp[ccon[i][j]].coor[2];
		}
		cp[ptcc[i]].coor[0] = tmp.coor[0];
		cp[ptcc[i]].coor[1] = tmp.coor[1];
		cp[ptcc[i]].coor[2] = tmp.coor[2];
	}
	int ptec[3] = {61, 71, 49}; //edge corner
	int econ[3][6] = {{61, 62, 50, 51, 72, 73}, {71, 82, 70, 81, 72, 83}, {49, 38, 48, 37, 50, 39}};
	double econc[6] = {15. / 24., 3. / 24., 5. / 48., 1. / 48., 5. / 48., 1. / 48.};
	for (int i = 0; i < 3; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 6; j++)
		{
			tmp.coor[0] += econc[j] * cp[econ[i][j]].coor[0];
			tmp.coor[1] += econc[j] * cp[econ[i][j]].coor[1];
			tmp.coor[2] += econc[j] * cp[econ[i][j]].coor[2];
		}
		cp[ptec[i]].coor[0] = tmp.coor[0];
		cp[ptec[i]].coor[1] = tmp.coor[1];
		cp[ptec[i]].coor[2] = tmp.coor[2];
	}
	int ptc[4] = {59, 70, 72, 50}; //corner
	int con[4][4] = {{59, 48, 58, 47}, {70, 69, 81, 80}, {72, 73, 83, 84}, {50, 39, 51, 40}};
	double conc[4] = {25. / 36., 5. / 36., 5. / 36., 1. / 36.};
	for (int i = 0; i < 4; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 4; j++)
		{
			tmp.coor[0] += conc[j] * cp[con[i][j]].coor[0];
			tmp.coor[1] += conc[j] * cp[con[i][j]].coor[1];
			tmp.coor[2] += conc[j] * cp[con[i][j]].coor[2];
		}
		cp[ptc[i]].coor[0] = tmp.coor[0];
		cp[ptc[i]].coor[1] = tmp.coor[1];
		cp[ptc[i]].coor[2] = tmp.coor[2];
	}
	double vtmp[3] = {(cp[47].coor[0] + 5. * cp[48].coor[0]) / 6., (cp[47].coor[1] + 5. * cp[48].coor[1]) / 6., (cp[47].coor[2] + 5. * cp[48].coor[2]) / 6.};
	cp[48].coor[0] = vtmp[0];
	cp[48].coor[1] = vtmp[1];
	cp[48].coor[2] = vtmp[2];

	//update knot vectors
	double ktsU[9] = {cp[59].knotU[0], cp[59].knotU[1], cp[59].knotU[2], (cp[59].knotU[2] + cp[59].knotU[3]) / 2., cp[59].knotU[3], (cp[59].knotU[3] + cp[59].knotU[4]) / 2., cp[59].knotU[4], cp[61].knotU[3], cp[61].knotU[4]};
	double ktsV[9] = {cp[49].knotV[0], cp[49].knotV[1], cp[49].knotV[2], (cp[49].knotV[2] + cp[49].knotV[3]) / 2., cp[49].knotV[3], (cp[49].knotV[3] + cp[49].knotV[4]) / 2., cp[49].knotV[4], cp[71].knotV[3], cp[71].knotV[4]};
	int setKU[5][5] = {{48, 134, 49, 130, 50}, {-1, 135, 131, 132, 133}, {59, 121, 60, 126, 61}, {122, 123, 124, 127, 128}, {70, 125, 71, 129, 72}};
	int setKV[5][5] = {{-1, -1, 59, 122, 70}, {134, 135, 121, 123, 125}, {49, 131, 60, 124, 71}, {130, 132, 126, 127, 129}, {50, 133, 61, 128, 72}};
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			if (setKU[i][j] != -1)
			{
				for (int k = 0; k < 5; k++)
				{
					cp[setKU[i][j]].knotU[k] = ktsU[j + k];
				}
			}
		}
	}
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			if (setKV[i][j] != -1)
			{
				for (int k = 0; k < 5; k++)
				{
					cp[setKV[i][j]].knotV[k] = ktsV[j + k];
				}
			}
		}
	}

	cp[38].knotV[4] = cp[131].knotV[2];
	cp[39].knotV[4] = cp[133].knotV[2];
	cp[47].knotU[4] = cp[121].knotU[2];
	//cp[48].knotU[3]=cp[121].knotU[2]; cp[48].knotU[4]=cp[131].knotU[2];
	cp[48].knotV[4] = cp[122].knotV[2];
	cp[51].knotU[0] = cp[130].knotU[2];
	cp[58].knotU[4] = cp[121].knotU[2];
	cp[62].knotU[0] = cp[126].knotU[2];
	cp[69].knotU[4] = cp[125].knotU[2];
	cp[73].knotU[0] = cp[129].knotU[2];
	cp[81].knotV[0] = cp[122].knotV[2];
	cp[82].knotV[0] = cp[124].knotV[2];
	cp[83].knotV[0] = cp[128].knotV[2];

	cp[59].knotV[0] = cp[37].knotV[2];
	cp[59].knotV[1] = cp[48].knotV[2];
	cp[122].knotV[0] = cp[48].knotV[2];

	//int t1id[4]={48,50,81,83};
	//int t1cid[4][3]={{59,121,60},{61,126,60},{70,125,71},{72,129,71}};
	//double t1c[3]={5./36.,1./12.,1./48.};
	int t2id[3] = {40, 80, 84};
	int t2cid[3][1] = {{50}, {70}, {72}};
	double t2c[1] = {1. / 36.};
	//int t3id[2]={49,82};
	//int t3cid[2][3]={{60,121,126},{71,125,129}};
	//double t3c[3]={3./24.,1./12.,1./12.};
	//int t4id[4]={58,62,69,73};
	//int t4cid[4][2]={{59,122},{61,128},{70,122},{72,128}};
	//double t4c[2]={5./36.,1./12.};
	int t4id[1] = {37};
	int t4cid[1][2] = {{49, 134}};
	double t4c[2] = {1. / 48., 1. / 12.};
	//int t5id[1]={38};
	//int t5cid[1][2]={{49,130}};
	//double t5c[2]={3./24.,1./12.};
	//for(int i=0;i<4;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t1id[i]].tbf.push_back(t1cid[i][j]);
	//		cp[t1id[i]].tc.push_back(t1c[j]);
	//	}
	//}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t2id[i]].tbf.push_back(t2cid[i][j]);
			cp[t2id[i]].tc.push_back(t2c[j]);
		}
	}
	//for(int i=0;i<2;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t3id[i]].tbf.push_back(t3cid[i][j]);
	//		cp[t3id[i]].tc.push_back(t3c[j]);
	//	}
	//}
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			cp[t4id[i]].tbf.push_back(t4cid[i][j]);
			cp[t4id[i]].tc.push_back(t4c[j]);
		}
	}
	//for(int i=0;i<1;i++)
	//{
	//	for(int j=0;j<2;j++)
	//	{
	//		cp[t5id[i]].tbf.push_back(t5cid[i][j]);
	//		cp[t5id[i]].tc.push_back(t5c[j]);
	//	}
	//}

	//new elements
	int ne_new(15), npt_new(15);
	int enew[15][4] = {{59, 121, 123, 122}, {121, 60, 124, 123}, {123, 124, 71, 125}, {122, 123, 125, 70}, {60, 126, 127, 124}, {126, 61, 128, 127}, {127, 128, 72, 129}, {124, 127, 129, 71}, {49, 130, 132, 131}, {130, 50, 133, 132}, {132, 133, 61, 126}, {131, 132, 126, 60}, {48, 134, 121, 59}, {134, 49, 131, 135}, {135, 131, 60, 121}};
	int idtmp[15] = {121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135};
	for (int i = 0; i < ne_new; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		if (i < 4)
			etmp.IEN = tmesh[rid[0]].IEN;
		else if (i >= 4 && i < 8)
			etmp.IEN = tmesh[rid[1]].IEN;
		else if (i >= 8 && i < 12)
			etmp.IEN = tmesh[rid[2]].IEN;
		else
			etmp.IEN = tmesh[44].IEN;
		for (int j = 0; j < npt_new; j++)
			etmp.IEN.push_back(idtmp[j]);
		tmesh.push_back(etmp);
	}

	for (int i = 0; i < 3; i++)
		tmesh[rid[i]].act = 0;
	tmesh[44].act = 0;

	for (int i = 0; i < 100; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (int j = 0; j < npt_new; j++)
			{
				if (cp[idtmp[j]].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[idtmp[j]].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[idtmp[j]].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[idtmp[j]].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(idtmp[j]);
				}
			}
		}
	}
}

void TruncatedTspline::RefineTest10_1()
{
	int rid[1] = {104};
	for (int i = 0; i < 5; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int ptt[4] = {136, 137, 139, 140}; //T-junctions
	int ton[4][4] = {{60, 126, 131, 132}, {60, 124, 121, 123}, {126, 127, 61, 128}, {124, 127, 71, 129}};
	double tonc[4] = {7. / 16., 7. / 16., 1. / 16., 1. / 16.};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptt[i]].coor[0] += tonc[j] * cp[ton[i][j]].coor[0];
			cp[ptt[i]].coor[1] += tonc[j] * cp[ton[i][j]].coor[1];
			cp[ptt[i]].coor[2] += tonc[j] * cp[ton[i][j]].coor[2];
		}
	}
	int ptf[1] = {138}; //face
	int fon[1][4] = {{60, 126, 127, 124}};
	double fonc[4] = {1. / 4., 1. / 4., 1. / 4., 1. / 4.};
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[ptf[i]].coor[0] += fonc[j] * cp[fon[i][j]].coor[0];
			cp[ptf[i]].coor[1] += fonc[j] * cp[fon[i][j]].coor[1];
			cp[ptf[i]].coor[2] += fonc[j] * cp[fon[i][j]].coor[2];
		}
	}

	//int pte[4]={124,126,131,121};//edge
	//int eon[4][6]={{60,71,59,70,61,72},{60,61,49,50,71,72},{49,60,48,59,50,61},{59,60,48,49,70,71}};
	//double eonc[6]={3./8.,3./8.,1./16.,1./16.,1./16.,1./16.};
	//for(int i=0;i<4;i++)
	//{
	//	for(int j=0;j<6;j++)
	//	{
	//		cp[pte[i]].coor[0]+=eonc[j]*cp[eon[i][j]].coor[0];
	//		cp[pte[i]].coor[1]+=eonc[j]*cp[eon[i][j]].coor[1];
	//		cp[pte[i]].coor[2]+=eonc[j]*cp[eon[i][j]].coor[2];
	//	}
	//}

	//int ptcc[1]={60};//corner corner
	//int ccon[1][9]={60,49,59,61,71,48,50,70,72};
	//double cconc[9]={9./16.,3./32.,3./32.,3./32.,3./32.,1./64.,1./64.,1./64.,1./64.};
	//for(int i=0;i<1;i++)
	//{
	//	Vertex tmp;
	//	for(int j=0;j<9;j++)
	//	{
	//		tmp.coor[0]+=cconc[j]*cp[ccon[i][j]].coor[0];
	//		tmp.coor[1]+=cconc[j]*cp[ccon[i][j]].coor[1];
	//		tmp.coor[2]+=cconc[j]*cp[ccon[i][j]].coor[2];
	//	}
	//	cp[ptcc[i]].coor[0]=tmp.coor[0]; cp[ptcc[i]].coor[1]=tmp.coor[1]; cp[ptcc[i]].coor[2]=tmp.coor[2];
	//}
	//int ptec[3]={61,71,49};//edge corner
	//int econ[3][6]={{61,62,50,51,72,73},{71,82,70,81,72,83},{49,38,48,37,50,39}};
	//double econc[6]={15./24.,3./24.,5./48.,1./48.,5./48.,1./48.};
	//for(int i=0;i<3;i++)
	//{
	//	Vertex tmp;
	//	for(int j=0;j<6;j++)
	//	{
	//		tmp.coor[0]+=econc[j]*cp[econ[i][j]].coor[0];
	//		tmp.coor[1]+=econc[j]*cp[econ[i][j]].coor[1];
	//		tmp.coor[2]+=econc[j]*cp[econ[i][j]].coor[2];
	//	}
	//	cp[ptec[i]].coor[0]=tmp.coor[0]; cp[ptec[i]].coor[1]=tmp.coor[1]; cp[ptec[i]].coor[2]=tmp.coor[2];
	//}
	int ptc[4] = {60, 126, 124, 127}; //corner
	int con[4][4] = {{60, 121, 131, 135}, {126, 61, 132, 133}, {124, 123, 71, 125}, {127, 128, 129, 72}};
	double conc[4] = {49. / 64., 7. / 64., 7. / 64., 1. / 64.};
	for (int i = 0; i < 4; i++)
	{
		Vertex tmp;
		for (int j = 0; j < 4; j++)
		{
			tmp.coor[0] += conc[j] * cp[con[i][j]].coor[0];
			tmp.coor[1] += conc[j] * cp[con[i][j]].coor[1];
			tmp.coor[2] += conc[j] * cp[con[i][j]].coor[2];
		}
		cp[ptc[i]].coor[0] = tmp.coor[0];
		cp[ptc[i]].coor[1] = tmp.coor[1];
		cp[ptc[i]].coor[2] = tmp.coor[2];
	}

	//update knot vectors
	double ktsU[7] = {cp[60].knotU[0], cp[60].knotU[1], cp[60].knotU[2], (cp[60].knotU[2] + cp[60].knotU[3]) / 2., cp[60].knotU[3], cp[60].knotU[4], cp[126].knotU[4]};
	double ktsV[7] = {cp[60].knotV[0], cp[60].knotV[1], cp[60].knotV[2], (cp[60].knotV[2] + cp[60].knotV[3]) / 2., cp[60].knotV[3], cp[60].knotV[4], cp[124].knotV[4]};
	int setKU[3][3] = {{60, 136, 126}, {137, 138, 139}, {124, 140, 127}};
	int setKV[3][3] = {{60, 137, 124}, {136, 138, 140}, {126, 139, 127}};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (setKU[i][j] != -1)
			{
				for (int k = 0; k < 5; k++)
				{
					cp[setKU[i][j]].knotU[k] = ktsU[j + k];
				}
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (setKV[i][j] != -1)
			{
				for (int k = 0; k < 5; k++)
				{
					cp[setKV[i][j]].knotV[k] = ktsV[j + k];
				}
			}
		}
	}

	cp[131].knotV[4] = cp[137].knotV[2];
	cp[132].knotV[4] = cp[139].knotV[2];
	cp[121].knotU[4] = cp[136].knotU[2];
	cp[61].knotU[0] = cp[136].knotU[2];
	cp[123].knotU[4] = cp[140].knotU[2];
	cp[128].knotU[0] = cp[140].knotU[2];
	cp[71].knotV[0] = cp[137].knotV[2];
	cp[129].knotV[0] = cp[139].knotV[2];

	//int t1id[4]={48,50,81,83};
	//int t1cid[4][3]={{59,121,60},{61,126,60},{70,125,71},{72,129,71}};
	//double t1c[3]={5./36.,1./12.,1./48.};
	int t2id[4] = {60, 126, 124, 127};
	int t2cid[4][1] = {{135}, {133}, {125}, {72}};
	double t2c[1] = {1. / 64.};
	//int t3id[2]={49,82};
	//int t3cid[2][3]={{60,121,126},{71,125,129}};
	//double t3c[3]={3./24.,1./12.,1./12.};
	//int t4id[4]={58,62,69,73};
	//int t4cid[4][2]={{59,122},{61,128},{70,122},{72,128}};
	//double t4c[2]={5./36.,1./12.};
	//int t4id[1]={37};
	//int t4cid[1][2]={{49,134}};
	//double t4c[2]={1./48.,1./12.};
	//int t5id[1]={38};
	//int t5cid[1][2]={{49,130}};
	//double t5c[2]={3./24.,1./12.};
	//for(int i=0;i<4;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t1id[i]].tbf.push_back(t1cid[i][j]);
	//		cp[t1id[i]].tc.push_back(t1c[j]);
	//	}
	//}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			cp[t2id[i]].tbf.push_back(t2cid[i][j]);
			cp[t2id[i]].tc.push_back(t2c[j]);
		}
	}
	//for(int i=0;i<2;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cp[t3id[i]].tbf.push_back(t3cid[i][j]);
	//		cp[t3id[i]].tc.push_back(t3c[j]);
	//	}
	//}
	//for(int i=0;i<1;i++)
	//{
	//	for(int j=0;j<2;j++)
	//	{
	//		cp[t4id[i]].tbf.push_back(t4cid[i][j]);
	//		cp[t4id[i]].tc.push_back(t4c[j]);
	//	}
	//}
	//for(int i=0;i<1;i++)
	//{
	//	for(int j=0;j<2;j++)
	//	{
	//		cp[t5id[i]].tbf.push_back(t5cid[i][j]);
	//		cp[t5id[i]].tc.push_back(t5c[j]);
	//	}
	//}

	//new elements
	int ne_old = tmesh.size();
	int ne_new(4), npt_new(5);
	int enew[4][4] = {{60, 136, 138, 137}, {136, 126, 139, 138}, {137, 138, 140, 124}, {138, 139, 127, 140}};
	int idtmp[5] = {136, 137, 138, 139, 140};
	for (int i = 0; i < ne_new; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		//if(i<4)
		etmp.IEN = tmesh[rid[0]].IEN;
		//else if(i>=4 && i<8) etmp.IEN=tmesh[rid[1]].IEN;
		//else if(i>=8 && i<12) etmp.IEN=tmesh[rid[2]].IEN;
		//else etmp.IEN=tmesh[44].IEN;
		for (int j = 0; j < npt_new; j++)
			etmp.IEN.push_back(idtmp[j]);
		tmesh.push_back(etmp);
	}

	for (int i = 0; i < 1; i++)
		tmesh[rid[i]].act = 0;
	//tmesh[44].act=0;

	for (int i = 0; i < ne_old; i++)
	{
		if (tmesh[i].act == 1)
		{
			for (int j = 0; j < npt_new; j++)
			{
				if (cp[idtmp[j]].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[idtmp[j]].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[idtmp[j]].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[idtmp[j]].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(idtmp[j]);
				}
			}
		}
	}
}

void TruncatedTspline::Refine1()
{
	//int rid[3]={54,55,45};
	//for(int i=0; i<3; i++)
	//{
	//	ElementSubdivTest(rid[i]);
	//}
	//ElementSubdivTest(104);
	//ElementSubdivTest(107);
	const int niter(3);
	for (int i = 0; i < niter; i++)
	{
		cout << "step: " << i << '\n';
		vector<int> rid;
		IdentifyTest(rid);
		cout << "Refining...\n";
		for (uint j = 0; j < rid.size(); j++)
		{
			ElementSubdivTest(rid[j]);
		}
		//cout<<"Visualizing...\n";
		//if(i==niter-1)
		VisualizeControlMesh("test6/tsp_test_CM_" + to_string(i));
		//VisualizeVTK("test6/tsp_test_"+to_string(long long(i)));
	}
	//CollectActives();
	//cout<<"Bezier extracting...\n";
	//vector<BezierElement> bzmesh;
	//BezierExtract(bzmesh);
	//cout<<"Visualizing...\n";
	//BezierVTK("test2/tsp_test_bezier2",bzmesh);
}

void TruncatedTspline::Refine2()
{
	int rid[3] = {54, 55, 45};
	for (int i = 0; i < 3; i++)
	{
		ElementSubdivTest(rid[i]);
	}
	ElementSubdivTest(104);

	//VisualizeVTK("test6/patchtest_2_surf");

	CollectActives();
	//vector<BezierElement> bzmesh;
	//BezierExtract(bzmesh);
	//BezierVTK("test2/tsp_test_bezier1",bzmesh);
}

void TruncatedTspline::Refine3()
{
	int erid[4] = {44, 45, 54, 55};
	TopologyRefineTest();
	for (int i = 0; i < 4; i++)
	{
		ElementRefine(erid[i]);
	}
	Truncation();
	UpdateIEN();

	TopologyRefineTest1();
	ElementRefine1(104);
	Truncation();
	UpdateIEN();
}

void TruncatedTspline::Truncation()
{
	//update knot vectors
	uint i, j;
	//for(i=0; i<cp.size(); i++)
	//{
	//	if(cp[i].aff==1)
	//	{
	//		for(j=0; j<5; j++)
	//		{
	//			cp[i].knotU[j]=cp[i].kutmp[j];
	//			cp[i].knotV[j]=cp[i].kvtmp[j];
	//		}
	//		if(cp[i].update==1)
	//		{
	//			for(j=0; j<3; j++) cp[i].coor[j]=cp[i].coortmp[j];
	//			cp[i].update=0;
	//		}
	//		cp[i].aff=0;
	//	}
	//}

	//check support and truncation
	for (i = 0; i < cp.size(); i++)
	{
		for (j = i + 1; j < cp.size(); j++)
		{
			if (cp[j].knotU[0] >= cp[i].knotU[0] && cp[j].knotU[4] <= cp[i].knotU[4] && cp[j].knotV[0] >= cp[i].knotV[0] && cp[j].knotV[4] <= cp[i].knotV[4] && cp[j].min_itv[0] <= cp[i].min_itv[0] && cp[j].min_itv[1] <= cp[i].min_itv[1]) //i cover j
			{
				vector<int>::iterator it = find(cp[i].tbf.begin(), cp[i].tbf.end(), j);
				if (it == cp[i].tbf.end())
				{
					vector<double> ku1(10), kv1(10);
					vector<vector<double>> Tu, Tv;
					vector<double> ku(cp[i].knotU, cp[i].knotU + 5);
					vector<double> kv(cp[i].knotV, cp[i].knotV + 5);
					vector<double>::iterator it1, it2;
					it1 = set_union(cp[i].knotU, cp[i].knotU + 5, cp[j].knotU, cp[j].knotU + 5, ku1.begin());
					it2 = set_union(cp[i].knotV, cp[i].knotV + 5, cp[j].knotV, cp[j].knotV + 5, kv1.begin());
					ku1.resize(it1 - ku1.begin());
					kv1.resize(it2 - kv1.begin());
					TMatrix(ku, ku1, 3, Tu);
					TMatrix(kv, kv1, 3, Tv);
					it1 = search(ku1.begin(), ku1.end(), cp[j].knotU, cp[j].knotU + 5);
					it2 = search(kv1.begin(), kv1.end(), cp[j].knotV, cp[j].knotV + 5);
					int loc1 = it1 - ku1.begin();
					int loc2 = it2 - kv1.begin();
					double coef(Tu[loc1][0] * Tv[loc2][0]);
					if (coef > 0.)
					{
						//cout<<i<<' '<<j<<'\n';
						//getchar();
						cp[i].trun = 1;
						cp[i].tc.push_back(coef);
						cp[i].tbf.push_back(j);
					}
					/*if(it1==ku1.end() || it2==kv1.end())
					{
						cout<<"U:\n";
						for(uint k=0; k<5; k++)
						{
							cout<<cp[i].knotU[k]<<" ";
						}
						cout<<'\n';
						for(uint k=0; k<5; k++)
						{
							cout<<cp[j].knotU[k]<<" ";
						}
						cout<<'\n';
						for(uint k=0; k<ku1.size(); k++)
						{
							cout<<ku1[k]<<" ";
						}
						cout<<'\n';
						cout<<"V:\n";
						for(uint k=0; k<5; k++)
						{
							cout<<cp[i].knotV[k]<<" ";
						}
						cout<<'\n';
						for(uint k=0; k<5; k++)
						{
							cout<<cp[j].knotV[k]<<" ";
						}
						cout<<'\n';
						for(uint k=0; k<kv1.size(); k++)
						{
							cout<<kv1[k]<<" ";
						}
						cout<<'\n';
						getchar();
					}*/
				}
			}
			if (cp[i].knotU[0] >= cp[j].knotU[0] && cp[i].knotU[4] <= cp[j].knotU[4] && cp[i].knotV[0] >= cp[j].knotV[0] && cp[i].knotV[4] <= cp[j].knotV[4] && cp[i].min_itv[0] <= cp[j].min_itv[0] && cp[i].min_itv[1] <= cp[j].min_itv[1]) //j cover i
			{
				vector<int>::iterator it = find(cp[j].tbf.begin(), cp[j].tbf.end(), i);
				if (it == cp[j].tbf.end())
				{
					vector<double> ku1(10), kv1(10);
					vector<vector<double>> Tu, Tv;
					vector<double> ku(cp[j].knotU, cp[j].knotU + 5);
					vector<double> kv(cp[j].knotV, cp[j].knotV + 5);
					vector<double>::iterator it1, it2;
					it1 = set_union(cp[j].knotU, cp[j].knotU + 5, cp[i].knotU, cp[i].knotU + 5, ku1.begin());
					it2 = set_union(cp[j].knotV, cp[j].knotV + 5, cp[i].knotV, cp[i].knotV + 5, kv1.begin());
					ku1.resize(it1 - ku1.begin());
					kv1.resize(it2 - kv1.begin());
					TMatrix(ku, ku1, 3, Tu);
					TMatrix(kv, kv1, 3, Tv);
					it1 = search(ku1.begin(), ku1.end(), cp[i].knotU, cp[i].knotU + 5);
					it2 = search(kv1.begin(), kv1.end(), cp[i].knotV, cp[i].knotV + 5);
					int loc1 = it1 - ku1.begin();
					int loc2 = it2 - kv1.begin();
					double coef(Tu[loc1][0] * Tv[loc2][0]);
					if (coef > 0.)
					{
						//cout<<i<<' '<<j<<'\n';
						//getchar();
						cp[j].trun = 1;
						cp[j].tc.push_back(coef);
						cp[j].tbf.push_back(i);
					}
				}
			}
		}
	}
}

void TruncatedTspline::ElementRefine(int eid) //calculate control points
{
	uint i, j;
	for (i = 0; i < tmesh[eid].newpID.size(); i++) //new nodes
	{
		if (tmesh[eid].newpID[i] != -1)
		{
			int pnew(tmesh[eid].newpID[i]);
			double ptmp[3] = {0., 0., 0.};
			for (j = 0; j < tmesh[eid].IEN.size(); j++)
			{
				int pold(tmesh[eid].IEN[j]);
				if (cp[pnew].knotU[0] >= cp[pold].knotU[0] && cp[pnew].knotU[4] <= cp[pold].knotU[4] && cp[pnew].knotV[0] >= cp[pold].knotV[0] && cp[pnew].knotV[4] <= cp[pold].knotV[4] && cp[pold].trun == 0)
				{
					vector<double> ku1(10), kv1(10);
					vector<vector<double>> Tu, Tv;
					vector<double> ku(cp[pold].knotU, cp[pold].knotU + 5);
					vector<double> kv(cp[pold].knotV, cp[pold].knotV + 5);
					vector<double>::iterator it1, it2;
					it1 = set_union(cp[pold].knotU, cp[pold].knotU + 5, cp[pnew].knotU, cp[pnew].knotU + 5, ku1.begin());
					it2 = set_union(cp[pold].knotV, cp[pold].knotV + 5, cp[pnew].knotV, cp[pnew].knotV + 5, kv1.begin());
					ku1.resize(it1 - ku1.begin());
					kv1.resize(it2 - kv1.begin());
					TMatrix(ku, ku1, 3, Tu);
					TMatrix(kv, kv1, 3, Tv);
					it1 = search(ku1.begin(), ku1.end(), cp[pnew].knotU, cp[pnew].knotU + 5);
					it2 = search(kv1.begin(), kv1.end(), cp[pnew].knotV, cp[pnew].knotV + 5);
					int loc1 = it1 - ku1.begin();
					int loc2 = it2 - kv1.begin();
					double coef = Tu[loc1][0] * Tv[loc2][0];
					ptmp[0] += coef * cp[pold].coor[0];
					ptmp[1] += coef * cp[pold].coor[1];
					ptmp[2] += coef * cp[pold].coor[2];
				}
			}
			cp[pnew].coor[0] = ptmp[0];
			cp[pnew].coor[1] = ptmp[1];
			cp[pnew].coor[2] = ptmp[2];
		}
	}

	for (i = 0; i < 4; i++) //corner nodes
	{
		int pcn(tmesh[eid].cnct[i]);
		double ptmp[3] = {0., 0., 0.};
		for (j = 0; j < tmesh[eid].IEN.size(); j++)
		{
			int pold(tmesh[eid].IEN[j]);
			if (cp[pcn].kutmp[0] >= cp[pold].knotU[0] && cp[pcn].kutmp[4] <= cp[pold].knotU[4] && cp[pcn].kvtmp[0] >= cp[pold].knotV[0] && cp[pcn].kvtmp[4] <= cp[pold].knotV[4] && cp[pold].trun == 0)
			{
				vector<double> ku1(10), kv1(10);
				vector<vector<double>> Tu, Tv;
				vector<double> ku(cp[pold].knotU, cp[pold].knotU + 5);
				vector<double> kv(cp[pold].knotV, cp[pold].knotV + 5);
				vector<double>::iterator it1, it2;
				it1 = set_union(cp[pold].knotU, cp[pold].knotU + 5, cp[pcn].kutmp, cp[pcn].kutmp + 5, ku1.begin());
				it2 = set_union(cp[pold].knotV, cp[pold].knotV + 5, cp[pcn].kvtmp, cp[pcn].kvtmp + 5, kv1.begin());
				ku1.resize(it1 - ku1.begin());
				kv1.resize(it2 - kv1.begin());
				TMatrix(ku, ku1, 3, Tu);
				TMatrix(kv, kv1, 3, Tv);
				it1 = search(ku1.begin(), ku1.end(), cp[pcn].kutmp, cp[pcn].kutmp + 5);
				it2 = search(kv1.begin(), kv1.end(), cp[pcn].kvtmp, cp[pcn].kvtmp + 5);
				int loc1 = it1 - ku1.begin();
				int loc2 = it2 - kv1.begin();
				double coef = Tu[loc1][0] * Tv[loc2][0];
				ptmp[0] += coef * cp[pold].coor[0];
				ptmp[1] += coef * cp[pold].coor[1];
				ptmp[2] += coef * cp[pold].coor[2];
			}
		}
		cp[pcn].update = 1;
		cp[pcn].coortmp[0] = ptmp[0];
		cp[pcn].coortmp[1] = ptmp[1];
		cp[pcn].coortmp[2] = ptmp[2];
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::ElementRefine1(int eid) //calculate control points
{
	uint i, j;
	for (i = 0; i < tmesh[eid].newpID.size(); i++) //new nodes
	{
		if (tmesh[eid].newpID[i] != -1)
		{
			int pnew(tmesh[eid].newpID[i]);
			double ptmp[3] = {0., 0., 0.};
			for (j = 0; j < tmesh[eid].IEN.size(); j++)
			{
				int pold(tmesh[eid].IEN[j]);
				if (cp[pnew].knotU[0] >= cp[pold].knotU[0] && cp[pnew].knotU[4] <= cp[pold].knotU[4] && cp[pnew].knotV[0] >= cp[pold].knotV[0] && cp[pnew].knotV[4] <= cp[pold].knotV[4] && cp[pold].trun == 0)
				{
					vector<double> ku1(10), kv1(10);
					vector<vector<double>> Tu, Tv;
					vector<double> ku(cp[pold].knotU, cp[pold].knotU + 5);
					vector<double> kv(cp[pold].knotV, cp[pold].knotV + 5);
					vector<double>::iterator it1, it2;
					it1 = set_union(cp[pold].knotU, cp[pold].knotU + 5, cp[pnew].knotU, cp[pnew].knotU + 5, ku1.begin());
					it2 = set_union(cp[pold].knotV, cp[pold].knotV + 5, cp[pnew].knotV, cp[pnew].knotV + 5, kv1.begin());
					ku1.resize(it1 - ku1.begin());
					kv1.resize(it2 - kv1.begin());
					TMatrix(ku, ku1, 3, Tu);
					TMatrix(kv, kv1, 3, Tv);
					it1 = search(ku1.begin(), ku1.end(), cp[pnew].knotU, cp[pnew].knotU + 5);
					it2 = search(kv1.begin(), kv1.end(), cp[pnew].knotV, cp[pnew].knotV + 5);
					int loc1 = it1 - ku1.begin();
					int loc2 = it2 - kv1.begin();
					//cout<<pnew<<' '<<pold<<'\n';
					//cout<<loc1<<' '<<loc2<<'\n';
					//for(int jl=0; jl<5; jl++) cout<<cp[pnew].knotV[jl]<<' ';
					//cout<<'\n';
					//for(uint il=0; il<kv1.size(); il++)
					//{
					//	cout<<kv1[il]<<' ';
					//}
					//cout<<'\n';
					//getchar();
					double coef = Tu[loc1][0] * Tv[loc2][0];
					ptmp[0] += coef * cp[pold].coor[0];
					ptmp[1] += coef * cp[pold].coor[1];
					ptmp[2] += coef * cp[pold].coor[2];
				}
			}
			cp[pnew].coor[0] = ptmp[0];
			cp[pnew].coor[1] = ptmp[1];
			cp[pnew].coor[2] = ptmp[2];
		}
	}

	for (i = 0; i < 4; i++) //corner nodes
	{
		int pcn(tmesh[eid].cnct[i]);
		double ptmp[3] = {0., 0., 0.};
		for (j = 0; j < tmesh[eid].IEN.size(); j++)
		{
			int pold(tmesh[eid].IEN[j]);
			if (cp[pcn].kutmp[0] >= cp[pold].knotU[0] && cp[pcn].kutmp[4] <= cp[pold].knotU[4] && cp[pcn].kvtmp[0] >= cp[pold].knotV[0] && cp[pcn].kvtmp[4] <= cp[pold].knotV[4] && cp[pold].trun == 0)
			{
				vector<double> ku1(10), kv1(10);
				vector<vector<double>> Tu, Tv;
				vector<double> ku(cp[pold].knotU, cp[pold].knotU + 5);
				vector<double> kv(cp[pold].knotV, cp[pold].knotV + 5);
				vector<double>::iterator it1, it2;
				it1 = set_union(cp[pold].knotU, cp[pold].knotU + 5, cp[pcn].kutmp, cp[pcn].kutmp + 5, ku1.begin());
				it2 = set_union(cp[pold].knotV, cp[pold].knotV + 5, cp[pcn].kvtmp, cp[pcn].kvtmp + 5, kv1.begin());
				ku1.resize(it1 - ku1.begin());
				kv1.resize(it2 - kv1.begin());
				TMatrix(ku, ku1, 3, Tu);
				TMatrix(kv, kv1, 3, Tv);
				it1 = search(ku1.begin(), ku1.end(), cp[pcn].kutmp, cp[pcn].kutmp + 5);
				it2 = search(kv1.begin(), kv1.end(), cp[pcn].kvtmp, cp[pcn].kvtmp + 5);
				int loc1 = it1 - ku1.begin();
				int loc2 = it2 - kv1.begin();
				double coef = Tu[loc1][0] * Tv[loc2][0];
				ptmp[0] += coef * cp[pold].coor[0];
				ptmp[1] += coef * cp[pold].coor[1];
				ptmp[2] += coef * cp[pold].coor[2];
			}
		}
		cp[pcn].update = 1;
		cp[pcn].coortmp[0] = ptmp[0];
		cp[pcn].coortmp[1] = ptmp[1];
		cp[pcn].coortmp[2] = ptmp[2];
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::TopologyRefine(const vector<int> &rid) //work for regular mesh
{
	int npold(cp.size() - 1), pid(cp.size() - 1); //start numbering
	int neold(tmesh.size() - 1), eid(tmesh.size() - 1);
	int nedold(tmedge.size() - 1), edid(tmedge.size() - 1);
	vector<int> ridsb, ridtp;
	StrongBalanceCheck(rid, ridsb);
	StrongBalanceRefine(ridsb);
	OneTjunctionCheck(ridtp);
	OneTjunctionRefine(ridtp);
	TargetRefine(rid);
	//OneTjunctionCheck(ridtp);
	//OneTjunctionRefine(ridtp);
	UpdateTopology();
	//VisualizeTMesh("test4/CM_6_0");
	//getchar();
	//topology check stage: One T-junction
	int niter(0);
	while (1)
	{
		niter++;
		if (niter > 20)
			break;
		vector<int> ridtp1, ridsb1;
		OneTjunctionCheck(ridtp1);
		if (ridtp1.size() == 0)
			break;
		StrongBalanceCheck(ridtp1, ridsb1);
		StrongBalanceRefine(ridsb1);
		OneTjunctionRefine(ridtp1);
		UpdateTopology();
		//VisualizeTMesh("test4/CM_6_"+to_string(long long(niter)));
		//getchar();
	}
	TjuncExtentCheck();
	UpdateTopology();

	FindLocalKnotVectors();

	//		tmesh[i].ref=1;
	//		//int crnid[4]={tmesh[i].cnct[0],tmesh[i].cnct[1],tmesh[i].cnct[2],tmesh[i].cnct[3]};
	//		tmesh[i].nT=0;//if this element is identified as refined, there will be no T-junctions
	//		if(tmesh[i].square==1)
	//		{
	//			double umid((cp[tmesh[i].cnct[0]].knotU[2]+cp[tmesh[i].cnct[1]].knotU[2])/2.),vmid((cp[tmesh[i].cnct[0]].knotV[2]+cp[tmesh[i].cnct[3]].knotV[2])/2.);
	//			double ucoor[5]={umid,cp[tmesh[i].cnct[1]].knotU[2],umid,cp[tmesh[i].cnct[0]].knotU[2],umid};
	//			double vcoor[5]={cp[tmesh[i].cnct[0]].knotV[2],vmid,cp[tmesh[i].cnct[3]].knotV[2],vmid,vmid};
	//			for(j=0; j<4; j++)
	//			{
	//				if(tmesh[i].Tjunc[j]==-1)
	//				{
	//					tmesh[i].newpID[j]=++pid;
	//					Vertex vtmp;
	//					vtmp.pmcoor[0]=ucoor[j]; vtmp.pmcoor[1]=vcoor[j];
	//					cp.push_back(vtmp);
	//				}
	//				else
	//				{
	//					tmesh[i].newpID[j]=tmesh[i].Tjunc[j];
	//					tmesh[i].Tjunc[j]=-1;
	//				}
	//			}
	//			tmesh[i].newpID[4]=++pid;
	//			Vertex vtmp;
	//			vtmp.pmcoor[0]=ucoor[4]; vtmp.pmcoor[1]=vcoor[4];
	//			cp.push_back(vtmp);
	//			for(j=0; j<4; j++)//4 edge neighbors, which means a strongly balanced structure
	//			{
	//				if(tmesh[i].edge_nb[j][0]!=-1 && tmesh[tmesh[i].edge_nb[j][0]].ref==0)//edge neighbor is not yet refined
	//				{
	//					tmesh[tmesh[i].edge_nb[j][0]].Tjunc[tmesh[i].edge_nb[j][1]]=tmesh[i].newpID[j];
	//					tmesh[tmesh[i].edge_nb[j][0]].nT++;
	//				}
	//			}
	//			vector<Element> etmp(4);
	//			int enid[4][4]={{tmesh[i].cnct[0],tmesh[i].newpID[0],tmesh[i].newpID[4],tmesh[i].newpID[3]},{tmesh[i].newpID[0],tmesh[i].cnct[1],tmesh[i].newpID[1],tmesh[i].newpID[4]},
	//							{tmesh[i].newpID[3],tmesh[i].newpID[4],tmesh[i].newpID[2],tmesh[i].cnct[3]},{tmesh[i].newpID[4],tmesh[i].newpID[1],tmesh[i].cnct[2],tmesh[i].newpID[2]}};
	//			for(j=0; j<4; j++)
	//			{
	//				for(int k=0; k<4; k++)
	//				{
	//					etmp[j].cnct[k]=enid[j][k];
	//				}
	//				tmesh[i].chd[j]=++eid;
	//				etmp[j].prt=i;
	//				tmesh.push_back[etmp[j]];
	//			}
	//		}
	//		else if(tmesh[i].square==0)
	//		{
	//			int Tloc;
	//			for(j=0; j<4; j++)
	//			{
	//				if(tmesh[i].Tjunc[j]!=-1)
	//				{
	//					Tloc=j;
	//					break;
	//				}
	//			}
	//			if(tmesh[i].Tjunc[(Tloc+2)%4]!=-1)
	//			{
	//				int enid[2][4]={};
	//			}
	//			else
	//			{
	//			}
	//		}
	//	}
	//}

	//for(i=0; i<tmesh.size(); i++)
	//{
	//	if(tmesh[i].ref==0)
	//	{
	//		if(tmesh[i].nT==2)
	//		{
	//			int loc(0);
	//			for(j=0; j<4; j++)
	//			{
	//				if(tmesh[i].Tjunc[j]!=-1)
	//				{
	//					loc=j;
	//					break;
	//				}
	//			}
	//			if(tmesh[i].Tjunc[loc+1]!=0)//2 T-junctions next to each other
	//			{
	//				tmesh[i].ref=2;
	//				tmesh[i].newpID.resize(5,-1);
	//				//tmesh[i].newpID[loc]=;
	//			}
	//			else//2 T-junctions opposite to each other
	//			{
	//			}
	//		}
	//		else if(tmesh[i].ref==3)
	//		{
	//		}
	//		else if(tmesh[i].ref==4)
	//		{
	//		}
	//	}
	//}
}

void TruncatedTspline::TopologyRefine_1(const vector<int> &rid)
{
	vector<int> ridsb, ridtp;
	StrongBalanceCheck(rid, ridsb);
	StrongBalanceRefine(ridsb);
	OneTjunctionCheck(ridtp);
	OneTjunctionRefine(ridtp);
	TargetRefine(rid);
	UpdateTopology();
	//VisualizeTMesh("test6/Topo_CM_1_0");
	//topology check stage: One T-junction
	int niter(0);
	while (1)
	{
		niter++;
		if (niter > 20)
			break;
		vector<int> ridtp1, ridsb1;
		OneTjunctionCheck(ridtp1);
		if (ridtp1.size() == 0)
			break;
		StrongBalanceCheck(ridtp1, ridsb1);
		StrongBalanceRefine(ridsb1);
		OneTjunctionRefine(ridtp1);
		UpdateTopology();
		//VisualizeTMesh("test6/Topo_CM_1_"+to_string(long long(niter)));
		//getchar();
	}

	//extend edges
	niter = 0;
	while (1)
	{
		niter++;
		if (niter > 20)
			break;
		vector<int> ridtjx, ridsb2;
		TjuncExtentCheck_2(ridtjx);
		if (ridtjx.size() == 0)
			break;
		StrongBalanceCheck(ridtjx, ridsb2);
		StrongBalanceRefine(ridsb2);
		TjuncExtentRefine(ridtjx);
		UpdateTopology();
		//VisualizeTMesh("test4/CM_6_"+to_string(long long(niter)));
		//getchar();
	}

	FindLocalKnotVectors();
}

void TruncatedTspline::TopologyRefine_2(const vector<int> &rid)
{
	vector<int> ridsb, ridtp;
	StrongBalanceCheck(rid, ridsb);
	StrongBalanceRefine(ridsb);
	OneTjunctionCheck(ridtp);
	OneTjunctionRefine(ridtp);
	TargetRefine(rid);
	UpdateTopology();
	//VisualizeTMesh("test6/Topo_CM_1_0");
	//topology check stage: One T-junction
	int niter(0);
	while (1)
	{
		niter++;
		if (niter > 20)
			break;
		vector<int> ridtp1, ridsb1;
		OneTjunctionCheck(ridtp1);
		if (ridtp1.size() == 0)
			break;
		StrongBalanceCheck(ridtp1, ridsb1);
		StrongBalanceRefine(ridsb1);
		OneTjunctionRefine(ridtp1);
		UpdateTopology();
		//VisualizeTMesh("test6/Topo_CM_1_"+to_string(long long(niter)));
		//getchar();
	}

	//extend edges
	niter = 0;
	while (1)
	{
		niter++;
		if (niter > 20)
			break;
		vector<int> ridtjx, ridsb2;
		TjuncExtentCheck_2(ridtjx);
		if (ridtjx.size() == 0)
			break;
		StrongBalanceCheck(ridtjx, ridsb2);
		StrongBalanceRefine(ridsb2);
		TjuncExtentRefine(ridtjx);
		UpdateTopology();
		//VisualizeTMesh("test4/CM_6_"+to_string(long long(niter)));
		//getchar();
	}

	FindLocalKnotVectors();
}

void TruncatedTspline::ElementTopologyRefine(int eid)
{
}

void TruncatedTspline::StrongBalanceCheck(const vector<int> &rid, vector<int> &rid2)
{
	uint i, j, k;
	rid2.clear();
	for (i = 0; i < rid.size(); i++)
	{
		if (tmesh[rid[i]].act == 1)
		{
			for (j = 0; j < 4; j++)
			{
				for (k = 0; k < cp[tmesh[rid[i]].cnct[j]].face.size(); k++)
				{
					int ktmp = cp[tmesh[rid[i]].cnct[j]].face[k];
					if (tmesh[ktmp].act == 1 && ktmp != rid[i])
					{
						double lvdf = tmesh[ktmp].lv - tmesh[rid[i]].lv;
						if (lvdf != 0. && lvdf != 0.5 && lvdf != 1.)
						{
							vector<int>::iterator it = find(rid2.begin(), rid2.end(), ktmp);
							if (it == rid2.end())
								rid2.push_back(ktmp);
						}
					}
				}
			}
		}
	}
}

void TruncatedTspline::OneTjunctionCheck(vector<int> &rid2)
{
	rid2.clear();
	for (uint i = 0; i < tmesh.size(); i++)
	{
		tmesh[i].ref = 0;
		if (tmesh[i].act == 1)
		{
			int ntjc(0);
			for (int j = 0; j < 4; j++)
			{
				if (tmedge[tmesh[i].edge[j]].act == 0)
				{
					ntjc++;
				}
			}
			if (ntjc == 1 && tmesh[i].type == 2) //boundary element
			{
				tmesh[i].ref = 10;
				rid2.push_back(i);
			}
			else if (ntjc == 2)
			{
				int pos(0);
				for (int j = 0; j < 4; j++)
				{
					if (tmedge[tmesh[i].edge[j]].act == 0)
					{
						pos = j;
						break;
					}
				}
				if (tmedge[tmesh[i].edge[(pos + 1) % 4]].act == 0)
				{
					tmesh[i].ref = 20;
				}
				else if (tmedge[tmesh[i].edge[(pos + 2) % 4]].act == 0)
				{
					tmesh[i].ref = 21;
				}
				rid2.push_back(i);
			}
			else if (ntjc == 3)
			{
				tmesh[i].ref = 3;
				rid2.push_back(i);
			}
			else if (ntjc == 4)
			{
				tmesh[i].ref = 4;
				rid2.push_back(i);
			}
		}
	}
}

void TruncatedTspline::UpdateTopology()
{
	uint i, j, k;
	//first clear all active
	for (i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			cp[i].face.clear();
			cp[i].edge.clear();
		}
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].act == 1)
		{
			tmedge[i].face.clear();
		}
	}
	//loop all faces
	for (i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 /* || (tmesh[i].type==2 && tmesh[i].chd[0]==-1 && tmesh[i].chd[1]==-1) || tmesh[i].type==3*/)
		{
			for (j = 0; j < 4; j++)
			{
				cp[tmesh[i].cnct[j]].face.push_back(i);
				if (tmedge[tmesh[i].edge[j]].act == 1)
				{
					tmedge[tmesh[i].edge[j]].face.push_back(i);
				}
				else
				{
					cp[tmedge[tmesh[i].edge[j]].midpt].face.push_back(i);
					int chdid[2] = {tmedge[tmesh[i].edge[j]].chd[0], tmedge[tmesh[i].edge[j]].chd[1]};
					if (tmedge[chdid[0]].act == 1 && tmedge[chdid[1]].act == 1)
					{
						tmedge[chdid[0]].face.push_back(i);
						tmedge[chdid[1]].face.push_back(i);
					}
					else
					{
						cerr << "Configuration not recognized!\n";
					}
				}
			}
		}
	}
	//loop all edges
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].act == 1)
		{
			cp[tmedge[i].pt[0]].edge.push_back(i);
			cp[tmedge[i].pt[1]].edge.push_back(i);
		}
	}

	//update corner element level
	for (i = 0; i < cornerEID.size(); i++)
	{
		double lev(tmesh[cornerEID[i]].lv);
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < cp[tmesh[cornerEID[i]].cnct[j]].face.size(); k++)
			{
				if (cp[tmesh[cornerEID[i]].cnct[j]].face[k] != cornerEID[k])
					lev = tmesh[cp[tmesh[cornerEID[i]].cnct[j]].face[k]].lv;
			}
		}
		tmesh[cornerEID[i]].lv = lev;
	}
}

void TruncatedTspline::InitialTopology()
{
	uint i, j;
	tmedge.clear();
	for (i = 0; i < tmesh.size(); i++)
	{
		tmesh[i].act = 1;
		tmesh[i].type = 0;
		double lentmp[4];
		for (j = 0; j < 4; j++)
		{
			cp[tmesh[i].cnct[j]].face.push_back(i);
			Edge edtmp;
			edtmp.act = 1;
			edtmp.pt[0] = tmesh[i].cnct[j];
			edtmp.pt[1] = tmesh[i].cnct[(j + 1) % 4];
			if (cp[tmesh[i].cnct[j]].index[0] == cp[tmesh[i].cnct[(j + 1) % 4]].index[0])
			{
				edtmp.len = cp[tmesh[i].cnct[(j + 1) % 4]].pm[1] - cp[tmesh[i].cnct[j]].pm[1];
				if (edtmp.len < 0.)
					edtmp.len = -edtmp.len;
				lentmp[j] = edtmp.len;
			}
			else
			{
				edtmp.len = cp[tmesh[i].cnct[(j + 1) % 4]].pm[0] - cp[tmesh[i].cnct[j]].pm[0];
				if (edtmp.len < 0.)
					edtmp.len = -edtmp.len;
				lentmp[j] = edtmp.len;
			}
			vector<Edge>::iterator it = find(tmedge.begin(), tmedge.end(), edtmp);
			int edid(it - tmedge.begin());
			if (it == tmedge.end())
			{
				tmedge.push_back(edtmp);
			}
			tmesh[i].edge[j] = edid;
			tmedge[edid].face.push_back(i);
		}
		if ((lentmp[0] == 0. && lentmp[1] > 0.) || (lentmp[0] > 0. && lentmp[1] == 0.))
		{
			tmesh[i].type = 2;
		}
		else if (lentmp[0] == 0. && lentmp[1] == 0.)
		{
			tmesh[i].type = 3;
			vector<int>::iterator it = find(cornerEID.begin(), cornerEID.end(), i);
			if (it == cornerEID.end())
				cornerEID.push_back(i);
		}
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		cp[tmedge[i].pt[0]].edge.push_back(i);
		cp[tmedge[i].pt[1]].edge.push_back(i);
	}
	/*for(i=0; i<tmedge.size(); i++)
	{
		if(tmedge[i].face.size()==1)
		{
			int eid(tmedge[i].face[0]);
			tmesh[eid].act=1;
			tmesh[eid].type=2;
			int* it=find(tmesh[eid].edge,tmesh[eid].edge+3,i);
			int loc(it-tmesh[eid].edge);
			tmedge[tmesh[eid].edge[(loc+1)%4]].len=0.;
			tmedge[tmesh[eid].edge[(loc+3)%4]].len=0.;
		}
	}
	cornerEID.clear();
	for(i=0; i<cp.size(); i++)
	{
		if(cp[i].face.size()==1)
		{
			int eid(cp[i].face[0]);
			tmesh[eid].act=1;
			tmesh[eid].type=3;
			vector<int>::iterator it=find(cornerEID.begin(),cornerEID.end(),eid);
			if(it==cornerEID.end())
				cornerEID.push_back(eid);
		}
	}*/

	FindLocalKnotVectors();
	UpdateKnotVectors();
	FindIENglb();
}

void TruncatedTspline::StrongBalanceRefine(const vector<int> &ridsb)
{
	for (uint i = 0; i < ridsb.size(); i++)
	{
		if (tmesh[ridsb[i]].act == 1)
		{
			if (tmesh[ridsb[i]].type == 0)
			{
				ElementRefine_Square_4(ridsb[i]);
			}
			else if (tmesh[ridsb[i]].type == 1)
			{
				ElementRefine_Rectangular(ridsb[i]);
			}
			else if (tmesh[ridsb[i]].type == 2)
			{
				ElementRefine_Boundary(ridsb[i]);
			}
		}
	}
}

void TruncatedTspline::TargetRefine(const vector<int> &rid)
{
	for (uint i = 0; i < rid.size(); i++)
	{
		if (tmesh[rid[i]].act == 1)
		{
			if (tmesh[rid[i]].type == 0)
			{
				ElementRefine_Square_4(rid[i]);
			}
			else if (tmesh[rid[i]].type == 1)
			{
				ElementRefine_Rectangular(rid[i]);
			}
		}
	}
}

void TruncatedTspline::OneTjunctionRefine(const vector<int> &ridtp)
{
	for (uint i = 0; i < ridtp.size(); i++)
	{
		if (tmesh[ridtp[i]].act == 1)
		{
			if (tmesh[ridtp[i]].ref == 10)
			{
				ElementRefine_Boundary(ridtp[i]);
			}
			if (tmesh[ridtp[i]].ref == 20 || tmesh[ridtp[i]].ref == 3)
			{
				//ElementRefine_Square_3(ridtp[i]);
				ElementRefine_Square_2(ridtp[i]);
			}
			else if (tmesh[ridtp[i]].ref == 21)
			{
				ElementRefine_Square_2(ridtp[i]);
			}
			else if (tmesh[ridtp[i]].ref == 4)
			{
				ElementRefine_Square_4(ridtp[i]);
			}
		}
	}
}

void TruncatedTspline::TjuncExtentCheck()
{
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].type == 0)
		{
			int pos(-1);
			for (int j = 0; j < 4; j++)
			{
				if (tmedge[tmesh[i].edge[j]].act == 0)
				{
					pos = j;
					break;
				}
			}
			if (pos != -1)
			{
				int ed[2] = {(pos + 1) % 4, (pos + 3) % 4}, ref(0);
				for (int k = 0; k < 2; k++)
				{
					int ednb(tmedge[tmesh[i].edge[ed[k]]].face[0]);
					if (ednb == i)
						ednb = tmedge[tmesh[i].edge[ed[k]]].face[1];
					int *it = find(tmesh[ednb].edge, tmesh[ednb].edge + 4, ed[k]);
					int loc = it - tmesh[ednb].edge;
					loc = (loc + 2) % 4;
					if (tmedge[tmesh[ednb].edge[loc]].act == 0)
						ref = 1;
				}
				if (ref == 1)
				{
					ElementRefine_Square_2(i);
				}
			}
		}
	}
}

void TruncatedTspline::TjuncExtentCheck_1(vector<int> &ridtjx)
{
	ridtjx.clear();
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			int pos(-1);
			for (int j = 0; j < 4; j++)
			{
				if (tmedge[tmesh[i].edge[j]].act == 0)
				{
					pos = j;
					break;
				}
			}
			if (pos != -1)
			{
				int ed[2] = {(pos + 1) % 4, (pos + 3) % 4}, ref(0);
				for (int k = 0; k < 2; k++)
				{
					int ednb(tmedge[tmesh[i].edge[ed[k]]].face[0]);
					if (ednb == i)
						ednb = tmedge[tmesh[i].edge[ed[k]]].face[1];
					int *it = find(tmesh[ednb].edge, tmesh[ednb].edge + 4, tmesh[i].edge[ed[k]]);
					if (it == tmesh[ednb].edge + 4)
					{
						vector<int>::iterator it1 = find(ridtjx.begin(), ridtjx.end(), ednb);
						if (it1 == ridtjx.end())
							ridtjx.push_back(ednb);
					}
				}
			}
		}
	}
}

void TruncatedTspline::TjuncExtentCheck_2(vector<int> &ridtjx)
{
	ridtjx.clear();
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].type == 0)
		{
			int pos(-1);
			for (int j = 0; j < 4; j++)
			{
				if (tmedge[tmesh[i].edge[j]].act == 0)
				{
					pos = j;
					break;
				}
			}
			if (pos != -1)
			{
				int ed[2] = {(pos + 1) % 4, (pos + 3) % 4}, ref(0);
				for (int k = 0; k < 2; k++)
				{
					int ednb(tmedge[tmesh[i].edge[ed[k]]].face[0]);
					if (ednb == i)
						ednb = tmedge[tmesh[i].edge[ed[k]]].face[1];
					int *it = find(tmesh[ednb].edge, tmesh[ednb].edge + 4, tmesh[i].edge[ed[k]]);
					int loc = it - tmesh[ednb].edge;
					loc = (loc + 2) % 4;
					if (tmedge[tmesh[ednb].edge[loc]].act == 0)
					{
						vector<int>::iterator it = find(ridtjx.begin(), ridtjx.end(), i);
						if (it == ridtjx.end())
							ridtjx.push_back(i);
					}
				}
			}
		}
	}
}

void TruncatedTspline::TjuncExtentRefine(const vector<int> &ridtjx)
{
	for (uint i = 0; i < ridtjx.size(); i++)
	{
		if (tmesh[ridtjx[i]].act == 1)
		{
			if (tmesh[ridtjx[i]].type == 0)
			{
				ElementRefine_Square_2(ridtjx[i]);
			}
			else if (tmesh[ridtjx[i]].type == 1)
			{
				ElementRefine_Rectangular(ridtjx[i]);
			}
			else if (tmesh[ridtjx[i]].type == 2)
			{
				ElementRefine_Boundary(ridtjx[i]);
			}
		}
	}
}

void TruncatedTspline::ElementRefine_Square_4(int eid)
{
	int pid[5];
	int edid[12];
	Vertex ptmp1;
	ptmp1.pm[0] = (cp[tmesh[eid].cnct[0]].pm[0] + cp[tmesh[eid].cnct[1]].pm[0]) / 2.;
	ptmp1.pm[1] = (cp[tmesh[eid].cnct[0]].pm[1] + cp[tmesh[eid].cnct[3]].pm[1]) / 2.;
	AssignIndex_NewFacePoint(eid, ptmp1);
	cp.push_back(ptmp1);
	pid[0] = cp.size() - 1;
	for (int j = 0; j < 4; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[j], tmesh[eid].cnct[(j + 1) % 4]};
		if (tmedge[tmesh[eid].edge[j]].act == 1)
		{
			Vertex ptmp;
			ptmp.pm[0] = (cp[itmp[0]].pm[0] + cp[itmp[1]].pm[0]) / 2.;
			ptmp.pm[1] = (cp[itmp[0]].pm[1] + cp[itmp[1]].pm[1]) / 2.;
			AssignIndex_NewEdgePoint(tmesh[eid].edge[j], ptmp);
			cp.push_back(ptmp);
			pid[j + 1] = cp.size() - 1;
			tmedge[tmesh[eid].edge[j]].midpt = pid[j + 1];
			//int ednb(tmedge[tmesh[eid].edge[j]].face[0]);
			//if(ednb==eid) ednb=tmedge[tmesh[eid].edge[j]].face[1];
			//cp[pid[j+1]].face.push_back(ednb);
			Edge edtmp1, edtmp2, edtmp3;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j + 1];
			edtmp1.len = tmedge[tmesh[eid].edge[j]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[j];
			//edtmp1.face.push_back(ednb);
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j + 1];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[j]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[j];
			//edtmp2.face.push_back(ednb);
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[(j + 1) % 4]].len / 2.;
			tmedge.push_back(edtmp1);
			edid[3 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].chd[0] = edid[3 * j];
			tmedge.push_back(edtmp2);
			edid[3 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].chd[1] = edid[3 * j + 1];
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].act = 0;
		}
		else
		{
			pid[j + 1] = tmedge[tmesh[eid].edge[j]].midpt;
			Edge edtmp3;
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[(j + 1) % 4]].len / 2.;
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			int ied(tmedge[tmesh[eid].edge[j]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[3 * j] = ied;
				edid[3 * j + 1] = tmedge[tmesh[eid].edge[j]].chd[1];
			}
			else
			{
				edid[3 * j] = tmedge[tmesh[eid].edge[j]].chd[1];
				edid[3 * j + 1] = ied;
			}
		}
	}

	int e_cnct[4][4] = {{tmesh[eid].cnct[0], pid[1], pid[0], pid[4]}, {pid[1], tmesh[eid].cnct[1], pid[2], pid[0]}, {pid[0], pid[2], tmesh[eid].cnct[2], pid[3]}, {pid[4], pid[0], pid[3], tmesh[eid].cnct[3]}};
	int e_edge[4][4] = {{edid[0], edid[2], edid[11], edid[10]}, {edid[1], edid[3], edid[5], edid[2]}, {edid[5], edid[4], edid[6], edid[8]}, {edid[11], edid[8], edid[7], edid[9]}};
	int enewid[4];
	vector<Element> etmp(4);
	for (int i = 0; i < 4; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 0;
		etmp[i].lv = tmesh[eid].lv + 1.;
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[j] = e_cnct[i][j];
			etmp[i].edge[j] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
	}

	tmesh[eid].act = 0;

	//int vf[4];
	//for(int i=0; i<4; i++)
	//{
	//	vector<int>::iterator it=find(cp[tmesh[eid].cnct[i]].face.begin(),cp[tmesh[eid].cnct[i]].face.end(),eid);
	//	vf[4]=it-cp[tmesh[eid].cnct[i]].face.begin();
	//}
	//for(int i=0; i<4; i++)
	//{
	//	cp[tmesh[eid].cnct[i]].face[vf[i]]=enewid[i];
	//	//cp[pid[0]].face.push_back(enewid[i]);
	//}
}

void TruncatedTspline::ElementRefine_Square_3(int eid)
{
	int pid[4], edid[11], pos(0);
	for (int i = 0; i < 4; i++)
	{
		if (tmedge[tmesh[eid].edge[i]].act == 0)
		{
			pos = i;
			break;
		}
	}
	int cnid[4] = {pos, (pos + 1) % 4, (pos + 2) % 4, (pos + 3) % 4};
	Vertex ptmp1;
	ptmp1.pm[0] = (cp[tmesh[eid].cnct[0]].pm[0] + cp[tmesh[eid].cnct[1]].pm[0]) / 2.;
	ptmp1.pm[1] = (cp[tmesh[eid].cnct[0]].pm[1] + cp[tmesh[eid].cnct[3]].pm[1]) / 2.;
	AssignIndex_NewFacePoint(eid, ptmp1);
	cp.push_back(ptmp1);
	pid[0] = cp.size() - 1;
	for (int j = 0; j < 3; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[cnid[j]], tmesh[eid].cnct[cnid[j + 1]]};
		if (tmedge[tmesh[eid].edge[cnid[j]]].act == 1)
		{
			Vertex ptmp;
			ptmp.pm[0] = (cp[itmp[0]].pm[0] + cp[itmp[1]].pm[0]) / 2.;
			ptmp.pm[1] = (cp[itmp[0]].pm[1] + cp[itmp[1]].pm[1]) / 2.;
			AssignIndex_NewEdgePoint(tmesh[eid].edge[cnid[j]], ptmp);
			cp.push_back(ptmp);
			pid[j + 1] = cp.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].midpt = pid[j + 1];
			Edge edtmp1, edtmp2, edtmp3;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j + 1];
			edtmp1.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[cnid[j]];
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j + 1];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[cnid[j]];
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[cnid[j + 1]]].len / 2.;
			tmedge.push_back(edtmp1);
			edid[3 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[0] = edid[3 * j];
			tmedge.push_back(edtmp2);
			edid[3 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[1] = edid[3 * j + 1];
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].act = 0;
		}
		else
		{
			pid[j + 1] = tmedge[tmesh[eid].edge[cnid[j]]].midpt;
			Edge edtmp3;
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[cnid[j + 1]]].len / 2.;
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			int ied(tmedge[tmesh[eid].edge[cnid[j]]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[3 * j] = ied;
				edid[3 * j + 1] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
			}
			else
			{
				edid[3 * j] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
				edid[3 * j + 1] = ied;
			}
		}
	}
	edid[9] = tmesh[eid].edge[cnid[3]];
	Edge edtmp;
	edtmp.act = 0;
	edtmp.pt[0] = pid[1];
	edtmp.pt[1] = pid[3];
	edtmp.len = tmedge[edid[9]].len;
	edtmp.chd[0] = edid[2];
	edtmp.chd[1] = edid[8];
	edtmp.midpt = pid[0];
	tmedge.push_back(edtmp);
	edid[10] = tmedge.size() - 1;
	tmedge[edid[2]].prt = edid[10];
	tmedge[edid[8]].prt = edid[10];

	int e_cnct[3][4] = {{tmesh[eid].cnct[cnid[0]], pid[1], pid[3], tmesh[eid].cnct[cnid[3]]}, {pid[1], tmesh[eid].cnct[cnid[1]], pid[2], pid[0]}, {pid[0], pid[2], tmesh[eid].cnct[cnid[2]], pid[3]}};
	int e_edge[3][4] = {{edid[0], edid[10], edid[7], edid[9]}, {edid[1], edid[3], edid[5], edid[2]}, {edid[5], edid[4], edid[6], edid[8]}};
	int enewid[3];
	vector<Element> etmp(3);
	for (int i = 0; i < 3; i++)
	{
		etmp[i].act = 1;
		if (i == 0)
		{
			etmp[i].type = 1;
			etmp[i].lv = tmesh[eid].lv + 0.5;
		}
		else
		{
			etmp[i].type = 0;
			etmp[i].lv = tmesh[eid].lv + 1.;
		}
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[(pos + j) % 4] = e_cnct[i][j];
			etmp[i].edge[(pos + j) % 4] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
	}

	tmesh[eid].act = 0;

	//int vf[4];
	//for(int i=0; i<4; i++)
	//{
	//	vector<int>::iterator it=find(cp[tmesh[eid].cnct[i]].face.begin(),cp[tmesh[eid].cnct[i]].face.end(),eid);
	//	vf[4]=it-cp[tmesh[eid].cnct[i]].face.begin();
	//}
	//for(int i=0; i<3; i++)
	//{
	//	cp[tmesh[eid].cnct[cnid[i]]].face[vf[cnid[i]]]=enewid[i];
	//}
	//cp[tmesh[eid].cnct[cnid[3]]].face[vf[cnid[3]]]=enewid[0];
}

void TruncatedTspline::ElementRefine_Square_2(int eid)
{
	int pid[2], edid[7], pos(0);
	for (int i = 0; i < 4; i++)
	{
		if (tmedge[tmesh[eid].edge[i]].act == 0)
		{
			pos = i;
			break;
		}
	}
	int cnid[2] = {pos, (pos + 2) % 4};
	for (int j = 0; j < 2; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[cnid[j]], tmesh[eid].cnct[(cnid[j] + 1) % 4]};
		if (tmedge[tmesh[eid].edge[cnid[j]]].act == 1)
		{
			Vertex ptmp;
			ptmp.pm[0] = (cp[itmp[0]].pm[0] + cp[itmp[1]].pm[0]) / 2.;
			ptmp.pm[1] = (cp[itmp[0]].pm[1] + cp[itmp[1]].pm[1]) / 2.;
			AssignIndex_NewEdgePoint(tmesh[eid].edge[cnid[j]], ptmp);
			cp.push_back(ptmp);
			pid[j] = cp.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].midpt = pid[j];
			Edge edtmp1, edtmp2, edtmp3;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j];
			edtmp1.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[cnid[j]];
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[cnid[j]];
			tmedge.push_back(edtmp1);
			edid[2 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[0] = edid[2 * j];
			tmedge.push_back(edtmp2);
			edid[2 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[1] = edid[2 * j + 1];
			tmedge[tmesh[eid].edge[cnid[j]]].act = 0;
		}
		else
		{
			pid[j] = tmedge[tmesh[eid].edge[cnid[j]]].midpt;
			int ied(tmedge[tmesh[eid].edge[cnid[j]]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[2 * j] = ied;
				edid[2 * j + 1] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
			}
			else
			{
				edid[2 * j] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
				edid[2 * j + 1] = ied;
			}
		}
	}
	Edge edtmp;
	edtmp.act = 1;
	edtmp.pt[0] = pid[0];
	edtmp.pt[1] = pid[1];
	edtmp.len = tmedge[tmesh[eid].edge[(pos + 1) % 4]].len;
	tmedge.push_back(edtmp);
	edid[4] = tmedge.size() - 1;
	edid[5] = tmesh[eid].edge[(pos + 1) % 4];
	edid[6] = tmesh[eid].edge[(pos + 3) % 4];

	int e_cnct[2][4] = {{tmesh[eid].cnct[pos], pid[0], pid[1], tmesh[eid].cnct[(pos + 3) % 4]}, {pid[0], tmesh[eid].cnct[(pos + 1) % 4], tmesh[eid].cnct[(pos + 2) % 4], pid[1]}};
	int e_edge[2][4] = {{edid[0], edid[4], edid[3], edid[6]}, {edid[1], edid[5], edid[2], edid[4]}};
	int enewid[2];
	vector<Element> etmp(2);
	for (int i = 0; i < 2; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 1;
		etmp[i].lv = tmesh[eid].lv + 0.5;
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[(pos + j) % 4] = e_cnct[i][j];
			etmp[i].edge[(pos + j) % 4] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
	}

	tmesh[eid].act = 0;

	//int vf[4];
	//for(int i=0; i<4; i++)
	//{
	//	vector<int>::iterator it=find(cp[tmesh[eid].cnct[i]].face.begin(),cp[tmesh[eid].cnct[i]].face.end(),eid);
	//	vf[4]=it-cp[tmesh[eid].cnct[i]].face.begin();
	//}
	//for(int i=0; i<2; i++)
	//{
	//	cp[tmesh[eid].cnct[(pos+i)%4]].face[vf[(pos+i)%4]]=enewid[i];
	//}
	//cp[tmesh[eid].cnct[(pos+2)%4]].face[vf[(pos+2)%4]]=enewid[1];
	//cp[tmesh[eid].cnct[(pos+3)%4]].face[vf[(pos+3)%4]]=enewid[0];
}

void TruncatedTspline::ElementRefine_Rectangular(int eid)
{
	int pos(0); //long edge position
	if (tmedge[tmesh[eid].edge[0]].len < tmedge[tmesh[eid].edge[1]].len)
		pos = 1;
	int ie[2] = {pos, pos + 2};
	int pid[2];
	int edid[7];
	for (int i = 0; i < 2; i++)
	{
		int itmp[2] = {tmedge[tmesh[eid].edge[ie[i]]].pt[0], tmedge[tmesh[eid].edge[ie[i]]].pt[1]};
		if (tmedge[tmesh[eid].edge[ie[i]]].act == 1)
		{
			Vertex ptmp1;
			ptmp1.pm[0] = (cp[itmp[0]].pm[0] + cp[itmp[1]].pm[0]) / 2.;
			ptmp1.pm[1] = (cp[itmp[0]].pm[1] + cp[itmp[1]].pm[1]) / 2.;
			AssignIndex_NewEdgePoint(tmesh[eid].edge[ie[i]], ptmp1);
			cp.push_back(ptmp1);
			pid[i] = cp.size() - 1;
			vector<Edge> edtmp(2);
			edtmp[0].act = 1;
			edtmp[0].prt = tmesh[eid].edge[ie[i]];
			edtmp[0].len = tmedge[tmesh[eid].edge[ie[i]]].len / 2.;
			edtmp[0].pt[0] = tmesh[eid].cnct[ie[i]];
			edtmp[0].pt[1] = pid[i];
			tmedge.push_back(edtmp[0]);
			edid[2 * i] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[ie[i]]].chd[0] = edid[2 * i];
			edtmp[1].act = 1;
			edtmp[1].prt = tmesh[eid].edge[ie[i]];
			edtmp[1].len = tmedge[tmesh[eid].edge[ie[i]]].len / 2.;
			edtmp[1].pt[0] = pid[i];
			edtmp[1].pt[1] = tmesh[eid].cnct[(ie[i] + 1) % 4];
			tmedge.push_back(edtmp[1]);
			edid[2 * i + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[ie[i]]].chd[1] = edid[2 * i + 1];
			tmedge[tmesh[eid].edge[ie[i]]].act = 0;
			tmedge[tmesh[eid].edge[ie[i]]].midpt = pid[i];
		}
		else
		{
			pid[i] = tmedge[tmesh[eid].edge[ie[i]]].midpt;
			int ied(tmedge[tmesh[eid].edge[ie[i]]].chd[0]);
			if (tmedge[ied].pt[0] == tmesh[eid].cnct[ie[i]] || tmedge[ied].pt[1] == tmesh[eid].cnct[ie[i]])
			{
				edid[2 * i] = ied;
				edid[2 * i + 1] = tmedge[tmesh[eid].edge[ie[i]]].chd[1];
			}
			else
			{
				edid[2 * i] = tmedge[tmesh[eid].edge[ie[i]]].chd[1];
				edid[2 * i + 1] = ied;
			}
		}
	}
	Edge edtmp;
	edtmp.act = 1;
	edtmp.len = tmedge[tmesh[eid].edge[ie[0] + 1]].len;
	edtmp.pt[0] = pid[0];
	edtmp.pt[1] = pid[1];
	tmedge.push_back(edtmp);
	edid[4] = tmedge.size() - 1;
	edid[5] = tmesh[eid].edge[ie[0] + 1];
	edid[6] = tmesh[eid].edge[(ie[1] + 1) % 4];

	int e_cnct[2][4] = {{tmesh[eid].cnct[ie[0]], pid[0], pid[1], tmesh[eid].cnct[(ie[1] + 1) % 4]}, {pid[0], tmesh[eid].cnct[ie[0] + 1], tmesh[eid].cnct[ie[1]], pid[1]}};
	int e_edge[2][4] = {{edid[0], edid[4], edid[3], edid[6]}, {edid[1], edid[5], edid[2], edid[4]}};
	int enewid[2];
	vector<Element> etmp(2);
	for (int i = 0; i < 2; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 0;
		etmp[i].lv = tmesh[eid].lv + 0.5;
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[(pos + j) % 4] = e_cnct[i][j];
			etmp[i].edge[(pos + j) % 4] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
	}

	tmesh[eid].act = 0;

	//int vf[4];
	//for(int i=0; i<4; i++)
	//{
	//	vector<int>::iterator it=find(cp[tmesh[eid].cnct[i]].face.begin(),cp[tmesh[eid].cnct[i]].face.end(),eid);
	//	vf[4]=it-cp[tmesh[eid].cnct[i]].face.begin();
	//}
	//for(int i=0; i<2; i++)
	//{
	//	cp[tmesh[eid].cnct[(pos+i)%4]].face[vf[(pos+i)%4]]=enewid[i];
	//}
	//cp[tmesh[eid].cnct[(pos+2)%4]].face[vf[(pos+2)%4]]=enewid[1];
	//cp[tmesh[eid].cnct[(pos+3)%4]].face[vf[(pos+3)%4]]=enewid[0];
}

void TruncatedTspline::ElementRefine_Boundary(int eid)
{
	int pos(0); //long edge position
	if (tmedge[tmesh[eid].edge[0]].len < tmedge[tmesh[eid].edge[1]].len)
		pos = 1;
	int ie[2] = {pos, pos + 2};
	int pid[2];
	int edid[7];
	for (int i = 0; i < 2; i++)
	{
		int itmp[2] = {tmedge[tmesh[eid].edge[ie[i]]].pt[0], tmedge[tmesh[eid].edge[ie[i]]].pt[1]};
		if (tmedge[tmesh[eid].edge[ie[i]]].act == 1)
		{
			Vertex ptmp1;
			ptmp1.pm[0] = (cp[itmp[0]].pm[0] + cp[itmp[1]].pm[0]) / 2.;
			ptmp1.pm[1] = (cp[itmp[0]].pm[1] + cp[itmp[1]].pm[1]) / 2.;
			AssignIndex_NewEdgePoint(tmesh[eid].edge[ie[i]], ptmp1);
			cp.push_back(ptmp1);
			pid[i] = cp.size() - 1;
			vector<Edge> edtmp(2);
			edtmp[0].act = 1;
			edtmp[0].prt = tmesh[eid].edge[ie[i]];
			edtmp[0].len = tmedge[tmesh[eid].edge[ie[i]]].len / 2.;
			edtmp[0].pt[0] = tmesh[eid].cnct[ie[i]];
			edtmp[0].pt[1] = pid[i];
			tmedge.push_back(edtmp[0]);
			edid[2 * i] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[ie[i]]].chd[0] = edid[2 * i];
			edtmp[1].act = 1;
			edtmp[1].prt = tmesh[eid].edge[ie[i]];
			edtmp[1].len = tmedge[tmesh[eid].edge[ie[i]]].len / 2.;
			edtmp[1].pt[0] = pid[i];
			edtmp[1].pt[1] = tmesh[eid].cnct[(ie[i] + 1) % 4];
			tmedge.push_back(edtmp[1]);
			edid[2 * i + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[ie[i]]].chd[1] = edid[2 * i + 1];
			tmedge[tmesh[eid].edge[ie[i]]].act = 0;
			tmedge[tmesh[eid].edge[ie[i]]].midpt = pid[i];
		}
		else
		{
			pid[i] = tmedge[tmesh[eid].edge[ie[i]]].midpt;
			int ied(tmedge[tmesh[eid].edge[ie[i]]].chd[0]);
			if (tmedge[ied].pt[0] == tmesh[eid].cnct[ie[i]] || tmedge[ied].pt[1] == tmesh[eid].cnct[ie[i]])
			{
				edid[2 * i] = ied;
				edid[2 * i + 1] = tmedge[tmesh[eid].edge[ie[i]]].chd[1];
			}
			else
			{
				edid[2 * i] = tmedge[tmesh[eid].edge[ie[i]]].chd[1];
				edid[2 * i + 1] = ied;
			}
		}
	}
	Edge edtmp;
	edtmp.act = 1;
	edtmp.len = tmedge[tmesh[eid].edge[ie[0] + 1]].len;
	edtmp.pt[0] = pid[0];
	edtmp.pt[1] = pid[1];
	tmedge.push_back(edtmp);
	edid[4] = tmedge.size() - 1;
	edid[5] = tmesh[eid].edge[ie[0] + 1];
	edid[6] = tmesh[eid].edge[(ie[1] + 1) % 4];

	int e_cnct[2][4] = {{tmesh[eid].cnct[ie[0]], pid[0], pid[1], tmesh[eid].cnct[(ie[1] + 1) % 4]}, {pid[0], tmesh[eid].cnct[ie[0] + 1], tmesh[eid].cnct[ie[1]], pid[1]}};
	int e_edge[2][4] = {{edid[0], edid[4], edid[3], edid[6]}, {edid[1], edid[5], edid[2], edid[4]}};
	int enewid[2];
	vector<Element> etmp(2);
	for (int i = 0; i < 2; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 2;
		etmp[i].lv = tmesh[eid].lv + 1.;
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[(pos + j) % 4] = e_cnct[i][j];
			etmp[i].edge[(pos + j) % 4] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
	}

	tmesh[eid].act = 0;

	//int vf[4];
	//for(int i=0; i<4; i++)
	//{
	//	vector<int>::iterator it=find(cp[tmesh[eid].cnct[i]].face.begin(),cp[tmesh[eid].cnct[i]].face.end(),eid);
	//	vf[4]=it-cp[tmesh[eid].cnct[i]].face.begin();
	//}
	//for(int i=0; i<2; i++)
	//{
	//	cp[tmesh[eid].cnct[(pos+i)%4]].face[vf[(pos+i)%4]]=enewid[i];
	//}
	//cp[tmesh[eid].cnct[(pos+2)%4]].face[vf[(pos+2)%4]]=enewid[1];
	//cp[tmesh[eid].cnct[(pos+3)%4]].face[vf[(pos+3)%4]]=enewid[0];
}

void TruncatedTspline::ElementRefine_Square_2(int eid, int dir)
{
	int pid[2], edid[7], pos(dir);
	//for (int i = 0; i<4; i++)
	//{
	//	if (tmedge[tmesh[eid].edge[i]].act == 0)
	//	{
	//		pos = i;
	//		break;
	//	}
	//}
	int cnid[2] = {pos, (pos + 2) % 4};
	for (int j = 0; j < 2; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[cnid[j]], tmesh[eid].cnct[(cnid[j] + 1) % 4]};
		if (tmedge[tmesh[eid].edge[cnid[j]]].act == 1)
		{
			Vertex ptmp;
			ptmp.pm[0] = (cp[itmp[0]].pm[0] + cp[itmp[1]].pm[0]) / 2.;
			ptmp.pm[1] = (cp[itmp[0]].pm[1] + cp[itmp[1]].pm[1]) / 2.;
			AssignIndex_NewEdgePoint(tmesh[eid].edge[cnid[j]], ptmp);
			cp.push_back(ptmp);
			pid[j] = cp.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].midpt = pid[j];
			Edge edtmp1, edtmp2, edtmp3;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j];
			edtmp1.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[cnid[j]];
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[cnid[j]];
			tmedge.push_back(edtmp1);
			edid[2 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[0] = edid[2 * j];
			tmedge.push_back(edtmp2);
			edid[2 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[1] = edid[2 * j + 1];
			tmedge[tmesh[eid].edge[cnid[j]]].act = 0;
		}
		else
		{
			pid[j] = tmedge[tmesh[eid].edge[cnid[j]]].midpt;
			int ied(tmedge[tmesh[eid].edge[cnid[j]]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[2 * j] = ied;
				edid[2 * j + 1] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
			}
			else
			{
				edid[2 * j] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
				edid[2 * j + 1] = ied;
			}
		}
	}
	Edge edtmp;
	edtmp.act = 1;
	edtmp.pt[0] = pid[0];
	edtmp.pt[1] = pid[1];
	edtmp.len = tmedge[tmesh[eid].edge[(pos + 1) % 4]].len;
	tmedge.push_back(edtmp);
	edid[4] = tmedge.size() - 1;
	edid[5] = tmesh[eid].edge[(pos + 1) % 4];
	edid[6] = tmesh[eid].edge[(pos + 3) % 4];

	int e_cnct[2][4] = {{tmesh[eid].cnct[pos], pid[0], pid[1], tmesh[eid].cnct[(pos + 3) % 4]}, {pid[0], tmesh[eid].cnct[(pos + 1) % 4], tmesh[eid].cnct[(pos + 2) % 4], pid[1]}};
	int e_edge[2][4] = {{edid[0], edid[4], edid[3], edid[6]}, {edid[1], edid[5], edid[2], edid[4]}};
	int enewid[2];
	vector<Element> etmp(2);
	for (int i = 0; i < 2; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 1;
		etmp[i].lv = tmesh[eid].lv + 0.5;
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[(pos + j) % 4] = e_cnct[i][j];
			etmp[i].edge[(pos + j) % 4] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
	}

	tmesh[eid].act = 0;

	//int vf[4];
	//for(int i=0; i<4; i++)
	//{
	//	vector<int>::iterator it=find(cp[tmesh[eid].cnct[i]].face.begin(),cp[tmesh[eid].cnct[i]].face.end(),eid);
	//	vf[4]=it-cp[tmesh[eid].cnct[i]].face.begin();
	//}
	//for(int i=0; i<2; i++)
	//{
	//	cp[tmesh[eid].cnct[(pos+i)%4]].face[vf[(pos+i)%4]]=enewid[i];
	//}
	//cp[tmesh[eid].cnct[(pos+2)%4]].face[vf[(pos+2)%4]]=enewid[1];
	//cp[tmesh[eid].cnct[(pos+3)%4]].face[vf[(pos+3)%4]]=enewid[0];
}

void TruncatedTspline::VisualizeTMesh(string fn)
{
	string fname(fn + ".vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << cp.size() << " float\n";
		for (uint i = 0; i < cp.size(); i++)
		{
			fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
			//fout<<cp[i].pm[0]<<" "<<cp[i].pm[1]<<" 0\n";
			//fout<<cp[i].knotU[2]<<" "<<cp[i].knotV[2]<<" 0\n";
		}
		int ned(0);
		for (uint i = 0; i < tmedge.size(); i++)
		{
			if (tmedge[i].act == 1)
				ned++;
		}
		fout << "\nCELLS " << ned << " " << 3 * ned << '\n';
		for (uint i = 0; i < tmedge.size(); i++)
		{
			if (tmedge[i].act == 1)
				fout << "2 " << tmedge[i].pt[0] << " " << tmedge[i].pt[1] << '\n';
		}
		fout << "\nCELL_TYPES " << ned << '\n';
		for (uint i = 0; i < tmedge.size(); i++)
		{
			if (tmedge[i].act == 1)
				fout << "3\n";
		}
		fout << "\nCELL_DATA " << ned << "\nSCALARS len float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < tmedge.size(); i++)
		{
			if (tmedge[i].act == 1)
				fout << tmedge[i].len << "\n";
		}
		fout << "POINT_DATA " << cp.size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < cp.size(); i++)
		{
			fout << cp[i].trun << "\n";
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
		cout << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline::InitializeTopologyDirect()
{
	for (uint i = 0; i < tmedge.size(); i++)
	{
		tmedge[i].pn[0][0] = 3;
		tmedge[i].pn[1][0] = 3;
		tmedge[i].pn[0][1] = -1;
		tmedge[i].pn[1][1] = -1;
	}
	for (uint i = 0; i < tmesh.size(); i++)
	{
		for (uint j = 0; j < 4; j++)
		{
			tmesh[i].pn[j][0] = 3;
			tmesh[i].pn[j][1] = -1;
		}
	}
}

void TruncatedTspline::FindTopologyDirect()
{
	InitializeTopologyDirect();
	uint i, j, k, i1, j1;
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].act == 1)
		{
			for (j = 0; j < 2; j++)
			{
				if (cp[tmedge[i].pt[j]].type == 0) //regular
				{
					for (k = 0; k < cp[tmedge[i].pt[j]].edge.size(); k++)
					{
						if (cp[tmedge[i].pt[j]].edge[k] != i)
						{
							int flag(0);
							for (i1 = 0; i1 < tmedge[i].face.size(); i1++)
							{
								for (j1 = 0; j1 < 4; j1++)
								{
									int edid(tmesh[tmedge[i].face[i1]].edge[j1]);
									if (tmedge[edid].act == 1 && edid == cp[tmedge[i].pt[j]].edge[k])
									{
										flag = 1;
									}
									else if (tmedge[edid].act == 0 && (tmedge[edid].chd[0] == cp[tmedge[i].pt[j]].edge[k] || tmedge[edid].chd[1] == cp[tmedge[i].pt[j]].edge[k]))
									{
										flag = 1;
									}
								}
							}
							if (flag == 0)
							{
								tmedge[i].pn[j][0] = 0;
								tmedge[i].pn[j][1] = cp[tmedge[i].pt[j]].edge[k];
							}
						}
					}
				}
				else if (cp[tmedge[i].pt[j]].type == 1) //T-junction
				{
					int itf(cp[tmedge[i].pt[j]].face[0]);
					for (k = 0; k < cp[tmedge[i].pt[j]].face.size(); k++)
					{
						int found(0);
						for (i1 = 0; i1 < 4; i1++)
						{
							int fcid(cp[tmedge[i].pt[j]].face[k]);
							if (tmedge[tmesh[fcid].edge[i1]].act == 0 && tmedge[tmesh[fcid].edge[i1]].midpt == tmedge[i].pt[j])
							{
								itf = fcid;
								found = 1;
								break;
							}
						}
						if (found == 1)
							break;
					}
					vector<int>::iterator it = find(tmedge[i].face.begin(), tmedge[i].face.end(), itf);
					if (it != tmedge[i].face.end())
					{
						for (k = 0; k < cp[tmedge[i].pt[j]].edge.size(); k++)
						{
							int edid(cp[tmedge[i].pt[j]].edge[k]);
							vector<int>::iterator it1 = find(tmedge[edid].face.begin(), tmedge[edid].face.end(), itf);
							if (edid != i && it1 != tmedge[edid].face.end())
							{
								tmedge[i].pn[j][0] = 0;
								tmedge[i].pn[j][1] = edid;
								break;
							}
						}
					}
					else
					{
						tmedge[i].pn[j][0] = 1;
						tmedge[i].pn[j][1] = itf;
						int ied(0);
						for (k = 0; k < 4; k++)
						{
							if (tmedge[tmesh[itf].edge[k]].act == 0 && tmedge[tmesh[itf].edge[k]].midpt == tmedge[i].pt[j])
							{
								ied = k;
								break;
							}
						}
						if (ied == 0 || ied == 2)
							tmedge[i].pn[j][0] = 2;
					}
				}
				else //extraordinary
				{
				}
			}
		}
	}

	for (i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			for (j = 0; j < 4; j++)
			{
				if (tmedge[tmesh[i].edge[j]].act == 1)
				{
					if (tmedge[tmesh[i].edge[j]].face.size() == 2)
					{
						int fcid(tmedge[tmesh[i].edge[j]].face[0]);
						if (fcid == i)
							fcid = tmedge[tmesh[i].edge[j]].face[1];
						tmesh[i].pn[j][0] = 1;
						tmesh[i].pn[j][1] = fcid;
						int *it = find(tmesh[fcid].edge, tmesh[fcid].edge + 4, tmesh[i].edge[j]);
						int pos = it - tmesh[fcid].edge;
						if (pos == 0 || pos == 2)
							tmesh[i].pn[j][0] = 2;
					}
				}
				else
				{
					tmesh[i].pn[j][0] = 0;
					int pid(tmedge[tmesh[i].edge[j]].midpt);
					for (k = 0; k < cp[pid].edge.size(); k++)
					{
						if (cp[pid].edge[k] != tmedge[tmesh[i].edge[j]].chd[0] && cp[pid].edge[k] != tmedge[tmesh[i].edge[j]].chd[1])
						{
							tmesh[i].pn[j][1] = cp[pid].edge[k];
							break;
						}
					}
				}
			}
		}
	}
}

void TruncatedTspline::FindKnotInterval()
{
	FindTopologyDirect();
	uint i, j;
	for (i = 0; i < cp.size(); i++)
	{
		if (cp[i].type != 2)
		{
			int reffc(cp[i].face[0]);
			int *it = find(tmesh[reffc].cnct, tmesh[reffc].cnct + 4, i);
			if (it != tmesh[reffc].cnct + 4)
			{
				int pos = it - tmesh[reffc].cnct;
				int ued(tmesh[reffc].edge[pos]), ved(tmesh[reffc].edge[(pos + 3) % 4]);
				if (tmedge[ued].act == 0)
				{
					if (i == tmedge[tmedge[ued].chd[0]].pt[0] || i == tmedge[tmedge[ued].chd[0]].pt[1])
						ued = tmedge[ued].chd[0];
					else
						ued = tmedge[ued].chd[1];
				}
				if (tmedge[ved].act == 0)
				{
					if (i == tmedge[tmedge[ved].chd[0]].pt[0] || i == tmedge[tmedge[ved].chd[0]].pt[1])
						ved = tmedge[ved].chd[0];
					else
						ved = tmedge[ved].chd[1];
				}
				double ku[4], kv[4];
				ShootRay_Edge(ued, i, ku);
				ShootRay_Edge(ved, i, kv);
				if (pos == 0)
				{
					for (j = 0; j < 4; j++)
					{
						cp[i].kitvU[j] = ku[j];
						cp[i].kitvV[j] = kv[j];
					}
				}
				else if (pos == 1)
				{
					for (j = 0; j < 4; j++)
					{
						cp[i].kitvU[j] = kv[3 - j];
						cp[i].kitvV[j] = ku[j];
					}
				}
				else if (pos == 2)
				{
					for (j = 0; j < 4; j++)
					{
						cp[i].kitvU[j] = ku[3 - j];
						cp[i].kitvV[j] = kv[3 - j];
					}
				}
				else
				{
					for (j = 0; j < 4; j++)
					{
						cp[i].kitvU[j] = kv[j];
						cp[i].kitvV[j] = ku[3 - j];
					}
				}
			}
			else
			{
				int pos(0);
				for (j = 0; j < 4; j++)
				{
					if (tmedge[tmesh[reffc].edge[j]].act == 0 && tmedge[tmesh[reffc].edge[j]].midpt == i)
					{
						pos = j;
						break;
					}
				}
				int ued(tmedge[tmesh[reffc].edge[pos]].chd[0]);
				if (tmedge[tmedge[ued].chd[1]].pt[0] == tmesh[reffc].cnct[(pos + 1) % 4] || tmedge[tmedge[ued].chd[1]].pt[1] == tmesh[reffc].cnct[(pos + 1) % 4])
					ued = tmedge[tmesh[reffc].edge[pos]].chd[1];
				double ku[4], kv[4];
				ShootRay_Edge(ued, i, ku);
				ShootRay_Face(reffc, i, kv);
				if (pos == 0)
				{
					for (j = 0; j < 4; j++)
					{
						cp[i].kitvU[j] = ku[j];
						cp[i].kitvV[j] = kv[j];
					}
				}
				else if (pos == 1)
				{
					for (j = 0; j < 4; j++)
					{
						cp[i].kitvU[j] = kv[3 - j];
						cp[i].kitvV[j] = ku[j];
					}
				}
				else if (pos == 2)
				{
					for (j = 0; j < 4; j++)
					{
						cp[i].kitvU[j] = ku[3 - j];
						cp[i].kitvV[j] = kv[3 - j];
					}
				}
				else
				{
					for (j = 0; j < 4; j++)
					{
						cp[i].kitvU[j] = kv[j];
						cp[i].kitvV[j] = ku[3 - j];
					}
				}
			}
		}
	}
}

void TruncatedTspline::ShootRay_Edge(int edid, int pid, double kv[4])
{
	kv[2] = tmedge[edid].len;
	int start(0), end(1);
	if (tmedge[edid].pt[1] == pid)
	{
		start = 1;
		end = 0;
	}
	if (tmedge[edid].pn[end][0] == 0)
		kv[3] = tmedge[tmedge[edid].pn[end][1]].len;
	else if (tmedge[edid].pn[end][0] == 1)
	{
		int fcid(tmedge[edid].pn[end][1]);
		kv[3] = tmedge[tmesh[fcid].edge[0]].len;
	}
	else if (tmedge[edid].pn[end][0] == 2)
	{
		int fcid(tmedge[edid].pn[end][1]);
		kv[3] = tmedge[tmesh[fcid].edge[1]].len;
	}
	else
		kv[3] = 0.;

	if (tmedge[edid].pn[start][0] == 0)
	{
		kv[1] = tmedge[tmedge[edid].pn[start][1]].len;
		int edpre(tmedge[edid].pn[start][1]);
		int start1(0), end1(1);
		if (tmedge[edpre].pt[0] == pid)
		{
			start1 = 1;
			end1 = 0;
		}
		if (tmedge[edpre].pn[start1][0] == 0)
			kv[0] = tmedge[tmedge[edpre].pn[start1][1]].len;
		else if (tmedge[edpre].pn[start1][0] == 1)
		{
			int fcid(tmedge[edpre].pn[start1][1]);
			kv[0] = tmedge[tmesh[fcid].edge[0]].len;
		}
		else if (tmedge[edpre].pn[start1][0] == 2)
		{
			int fcid(tmedge[edpre].pn[start1][1]);
			kv[0] = tmedge[tmesh[fcid].edge[1]].len;
		}
		else
			kv[0] = 0.;
	}
	else if (tmedge[edid].pn[start][0] == 1)
	{
		int fcid(tmedge[edid].pn[start][1]);
		kv[1] = tmedge[tmesh[fcid].edge[0]].len;
		int loc;
		for (int j = 0; j < 4; j++)
		{
			if (tmesh[fcid].pn[j][0] == 0 && tmesh[fcid].pn[j][1] == edid)
			{
				loc = j;
				break;
			}
		}
		loc = (loc + 2) % 4;
		if (tmesh[fcid].pn[loc][0] == 0)
		{
			kv[0] = tmedge[tmesh[fcid].pn[loc][1]].len;
		}
		else if (tmesh[fcid].pn[loc][0] == 1)
		{
			kv[0] = tmedge[tmesh[tmesh[fcid].pn[loc][0]].edge[0]].len;
		}
		else if (tmesh[fcid].pn[loc][0] == 2)
		{
			kv[0] = tmedge[tmesh[tmesh[fcid].pn[loc][0]].edge[1]].len;
		}
		else
		{
			kv[0] = 0.;
		}
	}
	else if (tmedge[edid].pn[start][0] == 2)
	{
		int fcid(tmedge[edid].pn[start][1]);
		kv[1] = tmedge[tmesh[fcid].edge[1]].len;
		int loc;
		for (int j = 0; j < 4; j++)
		{
			if (tmesh[fcid].pn[j][0] == 0 && tmesh[fcid].pn[j][1] == edid)
			{
				loc = j;
				break;
			}
		}
		loc = (loc + 2) % 4;
		if (tmesh[fcid].pn[loc][0] == 0)
		{
			kv[0] = tmedge[tmesh[fcid].pn[loc][1]].len;
		}
		else if (tmesh[fcid].pn[loc][0] == 1)
		{
			kv[0] = tmedge[tmesh[tmesh[fcid].pn[loc][0]].edge[0]].len;
		}
		else if (tmesh[fcid].pn[loc][0] == 2)
		{
			kv[0] = tmedge[tmesh[tmesh[fcid].pn[loc][0]].edge[1]].len;
		}
		else
		{
			kv[0] = 0.;
		}
	}
	else
	{
		kv[1] = 0.;
		kv[0] = 0.;
	}
}

void TruncatedTspline::ShootRay_Face(int fcid, int pid, double kv[4])
{
	int dir(0), loc;
	for (int j = 0; j < 4; j++)
	{
		if (tmedge[tmesh[fcid].edge[j]].act == 0 && tmedge[tmesh[fcid].edge[j]].midpt == pid)
		{
			loc = j;
			break;
		}
	}
	if (loc == 0 || loc == 2)
		dir = 1;
	kv[2] = tmedge[tmesh[fcid].edge[dir]].len;
	int loc1 = (loc + 2) % 4;
	if (tmesh[fcid].pn[loc1][0] == 0)
	{
		kv[3] = tmedge[tmesh[fcid].pn[loc1][1]].len;
	}
	else if (tmesh[fcid].pn[loc1][0] == 1)
	{
		int fcid1(tmesh[fcid].pn[loc1][1]);
		kv[3] = tmedge[tmesh[fcid1].edge[0]].len;
	}
	else if (tmesh[fcid].pn[loc1][0] == 2)
	{
		int fcid1(tmesh[fcid].pn[loc1][1]);
		kv[3] = tmedge[tmesh[fcid1].edge[1]].len;
	}
	else
	{
		kv[3] = 0.;
	}

	if (tmesh[fcid].pn[loc][0] == 0)
	{
		int edid(tmesh[fcid].pn[loc][1]);
		kv[1] = tmedge[edid].len;
		int start(0), end(1);
		if (tmedge[edid].pt[1] == pid)
		{
			start = 1;
			end = 0;
		}
		if (tmedge[edid].pn[end][0] == 0)
		{
			kv[0] = tmedge[tmedge[edid].pn[end][1]].len;
		}
		else if (tmedge[edid].pn[end][0] == 1)
		{
			int fcid1(tmedge[edid].pn[end][1]);
			kv[0] = tmedge[tmesh[fcid1].edge[0]].len;
		}
		else if (tmedge[edid].pn[end][0] == 2)
		{
			int fcid1(tmedge[edid].pn[end][1]);
			kv[0] = tmedge[tmesh[fcid1].edge[0]].len;
		}
		else
		{
			kv[0] = 0.;
		}
	}
	else
	{
		cerr << "Wrong edge connectivity!\n";
	}
}

void TruncatedTspline::FindIEN() //some points missing
{
	uint i, j, k;
	for (i = 0; i < tmesh.size(); i++)
	{
		tmesh[i].IEN.clear();
		tmesh[i].UV.clear();
		tmesh[i].pmcoor.clear();
	}

	for (i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			double crtmp[8][2] = {{0., 0.}, {tmedge[tmesh[i].edge[0]].len, 0.}, {tmedge[tmesh[i].edge[0]].len, tmedge[tmesh[i].edge[1]].len}, {0., tmedge[tmesh[i].edge[1]].len}, {tmedge[tmesh[i].edge[0]].len / 2., 0.}, {tmedge[tmesh[i].edge[0]].len, tmedge[tmesh[i].edge[1]].len / 2.}, {tmedge[tmesh[i].edge[0]].len / 2., tmedge[tmesh[i].edge[1]].len}, {0., tmedge[tmesh[i].edge[1]].len / 2.}};
			vector<int> one_ring;
			one_ring.push_back(i);
			//loop edge neighbors
			for (j = 0; j < 4; j++)
			{
				tmesh[i].IEN.push_back(tmesh[i].cnct[j]);
				if (tmedge[tmesh[i].edge[j]].act == 1)
				{
					int ednb(tmedge[tmesh[i].edge[j]].face[0]);
					if (tmedge[tmesh[i].edge[j]].face.size() == 2 && ednb == i)
						ednb = tmedge[tmesh[i].edge[j]].face[1];
					one_ring.push_back(ednb);
					int *it = find(tmesh[ednb].edge, tmesh[ednb].edge + 4, tmesh[i].edge[j]);
					if (it != tmesh[ednb].edge + 4)
					{
						int pos = it - tmesh[ednb].edge;
						int ed0((pos + 3) % 4), ed1((pos + 1) % 4);
						if (tmedge[tmesh[ednb].edge[ed0]].act == 1)
						{
							vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), tmesh[ednb].cnct[ed0]);
							if (it1 == tmesh[i].IEN.end())
								tmesh[i].IEN.push_back(tmesh[ednb].cnct[ed0]);
						}
						else
						{
							vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), tmedge[tmesh[ednb].edge[ed0]].midpt);
							if (it1 == tmesh[i].IEN.end())
								tmesh[i].IEN.push_back(tmedge[tmesh[ednb].edge[ed0]].midpt);
						}
						if (tmedge[tmesh[ednb].edge[ed1]].act == 1)
						{
							vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), tmesh[ednb].cnct[(ed1 + 1) % 4]);
							if (it1 == tmesh[i].IEN.end())
								tmesh[i].IEN.push_back(tmesh[ednb].cnct[(ed1 + 1) % 4]);
						}
						else
						{
							vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), tmedge[tmesh[ednb].edge[ed1]].midpt);
							if (it1 == tmesh[i].IEN.end())
								tmesh[i].IEN.push_back(tmedge[tmesh[ednb].edge[ed1]].midpt);
						}
					}
					else
					{
						for (k = 0; k < 4; k++)
						{
							vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), tmesh[ednb].cnct[k]);
							if (it1 == tmesh[i].IEN.end())
								tmesh[i].IEN.push_back(tmesh[ednb].cnct[k]);
							if (tmedge[tmesh[ednb].edge[k]].act == 0)
							{
								vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), tmedge[tmesh[ednb].edge[k]].midpt);
								if (it1 == tmesh[i].IEN.end())
									tmesh[i].IEN.push_back(tmedge[tmesh[ednb].edge[k]].midpt);
							}
						}
					}
				}
				else
				{
					tmesh[i].IEN.push_back(tmedge[tmesh[i].edge[j]].midpt);
					for (k = 0; k < 2; k++)
					{
						int ied(tmedge[tmesh[i].edge[j]].chd[k]);
						int ednb(tmedge[ied].face[0]);
						if (tmedge[ied].face.size() == 2 && ednb == i)
							ednb = tmedge[ied].face[1];
						one_ring.push_back(ednb);
						int *it = find(tmesh[ednb].edge, tmesh[ednb].edge + 4, ied);
						if (it != tmesh[ednb].edge + 4)
						{
							int pos = it - tmesh[ednb].edge;
							int ed0((pos + 3) % 4), ed1((pos + 1) % 4);
							if (tmedge[tmesh[ednb].edge[ed0]].act == 1)
							{
								vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), tmesh[ednb].cnct[ed0]);
								if (it1 == tmesh[i].IEN.end())
									tmesh[i].IEN.push_back(tmesh[ednb].cnct[ed0]);
							}
							else
							{
								vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), tmedge[tmesh[ednb].edge[ed0]].midpt);
								if (it1 == tmesh[i].IEN.end())
									tmesh[i].IEN.push_back(tmedge[tmesh[ednb].edge[ed0]].midpt);
							}
							if (tmedge[tmesh[ednb].edge[ed1]].act == 1)
							{
								vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), tmesh[ednb].cnct[(ed1 + 1) % 4]);
								if (it1 == tmesh[i].IEN.end())
									tmesh[i].IEN.push_back(tmesh[ednb].cnct[(ed1 + 1) % 4]);
							}
							else
							{
								vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), tmedge[tmesh[ednb].edge[ed1]].midpt);
								if (it1 == tmesh[i].IEN.end())
									tmesh[i].IEN.push_back(tmedge[tmesh[ednb].edge[ed1]].midpt);
							}
						}
						else
						{
							for (k = 0; k < 4; k++)
							{
								vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), tmesh[ednb].cnct[k]);
								if (it1 == tmesh[i].IEN.end())
									tmesh[i].IEN.push_back(tmesh[ednb].cnct[k]);
								if (tmedge[tmesh[ednb].edge[k]].act == 0)
								{
									vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), tmedge[tmesh[ednb].edge[k]].midpt);
									if (it1 == tmesh[i].IEN.end())
										tmesh[i].IEN.push_back(tmedge[tmesh[ednb].edge[k]].midpt);
								}
							}
						}
					}
				}
			}

			//then loop corner neighbors
			for (j = 0; j < 4; j++)
			{
				int crnd(cp[tmesh[i].cnct[j]].face[0]);
				for (k = 0; k < cp[tmesh[i].cnct[j]].face.size(); k++)
				{
					vector<int>::iterator it = find(one_ring.begin(), one_ring.end(), cp[tmesh[i].cnct[j]].face[k]);
					if (it == one_ring.end())
					{
						crnd = cp[tmesh[i].cnct[j]].face[k];
						break;
					}
				}
				int *it = find(tmesh[crnd].cnct, tmesh[crnd].cnct + 4, tmesh[i].cnct[j]);
				int pos = it - tmesh[crnd].cnct;
				int ed[4] = {pos, (pos + 3) % 4, (pos + 1) % 4, (pos + 2) % 4};
				int itmp[4] = {(pos + 1) % 4, (pos + 3) % 4, (pos + 2) % 4, (pos + 2) % 4};
				for (k = 0; k < 4; k++)
				{
					if (tmedge[tmesh[crnd].edge[ed[k]]].act == 1)
					{
						vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), tmesh[crnd].cnct[itmp[k]]);
						if (it1 == tmesh[i].IEN.end())
							tmesh[i].IEN.push_back(tmesh[crnd].cnct[itmp[k]]);
					}
					else
					{
						vector<int>::iterator it1 = find(tmesh[i].IEN.begin(), tmesh[i].IEN.end(), tmedge[tmesh[crnd].edge[ed[k]]].midpt);
						if (it1 == tmesh[i].IEN.end())
							tmesh[i].IEN.push_back(tmedge[tmesh[crnd].edge[ed[k]]].midpt);
					}
				}
			}

			for (j = 0; j < tmesh[i].IEN.size(); i++)
			{
				int uvtmp = FindLocalUV(i, tmesh[i].IEN[j]);
				tmesh[i].UV.push_back(uvtmp);
			}
		}
	}
}

int TruncatedTspline::FindLocalUV(int eid, int pid) //some problem
{
	int pfc(cp[pid].face[0]);
	if (pfc == eid)
		return 0;
	vector<int>::iterator it = find(tmesh[eid].nb.begin(), tmesh[eid].nb.end(), pfc);
	if (it != tmesh[eid].nb.end())
	{
		int pos = it - tmesh[eid].nb.begin();
		return tmesh[eid].nbrot[pos];
	}
	else
	{
		int fctmp(tmesh[pfc].nb[0]), pos1(0);
		for (uint i = 0; i < tmesh[pfc].nb.size(); i++)
		{
			vector<int>::iterator it1 = find(tmesh[eid].nb.begin(), tmesh[eid].nb.end(), tmesh[pfc].nb[i]);
			if (it1 != tmesh[eid].nb.begin())
			{
				pos1 = it1 - tmesh[eid].nb.begin();
				fctmp = tmesh[pfc].nb[i];
				break;
			}
		}
		vector<int>::iterator it1 = find(tmesh[fctmp].nb.begin(), tmesh[fctmp].nb.end(), pfc);
		int pos2 = it1 - tmesh[fctmp].nb.begin();
		int rottmp = (tmesh[fctmp].nbrot[pos2] + tmesh[eid].nbrot[pos1]) % 4;
		return rottmp;
	}
}

void TruncatedTspline::RelativeElementDirect()
{
	double ijv[4][2][2] = {{{1., 0.}, {0., 1.}}, {{0., 1.}, {-1., 0.}}, {{-1., 0.}, {0., -1.}}, {{0, -1.}, {1., 0.}}};
	int rotcase[4][4] = {{2, 3, 0, 1}, {1, 2, 3, 0}, {0, 1, 2, 3}, {3, 0, 1, 2}};
	uint i, j, k;
	for (i = 0; i < tmesh.size(); i++)
	{
		tmesh[i].nb.clear();
		tmesh[i].nbrot.clear();
		tmesh[i].nbog.clear();
		//tmesh[i].rotM.clear();
	}
	for (i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			//edge neighbor first
			for (j = 0; j < 4; j++)
			{
				if (tmedge[tmesh[i].edge[j]].act == 1)
				{
					int fcid(tmedge[tmesh[i].edge[j]].face[0]);
					if (tmedge[tmesh[i].edge[j]].face.size() == 2 && fcid == i)
						fcid = tmedge[tmesh[i].edge[j]].face[1];
					vector<int>::iterator it = find(tmesh[i].nb.begin(), tmesh[i].nb.end(), fcid);
					if (it == tmesh[i].nb.end())
					{
						tmesh[i].nb.push_back(fcid);
						int *it1 = find(tmesh[fcid].edge, tmesh[fcid].edge + 4, tmesh[i].edge[j]);
						if (it1 != tmesh[fcid].edge + 4)
						{
							int pos = it1 - tmesh[fcid].edge;
							if ((j + 2) % 4 == pos)
								tmesh[i].nbrot.push_back(0);
							if ((j + 1) % 4 == pos)
								tmesh[i].nbrot.push_back(1);
							if (j == pos)
								tmesh[i].nbrot.push_back(2);
							if ((j + 3) % 4 == pos)
								tmesh[i].nbrot.push_back(3);
						}
						else
						{
							int pos(0);
							for (k = 0; k < 4; k++)
							{
								if (tmedge[tmesh[fcid].edge[k]].act == 0 && (tmedge[tmesh[fcid].edge[k]].chd[0] == tmesh[i].edge[j] || tmedge[tmesh[fcid].edge[k]].chd[1] == tmesh[i].edge[j]))
								{
									pos = k;
									break;
								}
							}
							if ((j + 2) % 4 == pos)
								tmesh[i].nbrot.push_back(0);
							if ((j + 1) % 4 == pos)
								tmesh[i].nbrot.push_back(1);
							if (j == pos)
								tmesh[i].nbrot.push_back(2);
							if ((j + 3) % 4 == pos)
								tmesh[i].nbrot.push_back(3);
						}
					}
				}
				else
				{
					for (k = 0; k < 2; k++)
					{
						int edid(tmedge[tmesh[i].edge[j]].chd[k]);
						int fcid(tmedge[edid].face[0]);
						if (tmedge[tmesh[i].edge[j]].face.size() == 2 && fcid == i)
							fcid = tmedge[tmesh[i].edge[j]].face[1];
						vector<int>::iterator it = find(tmesh[i].nb.begin(), tmesh[i].nb.end(), fcid);
						if (it == tmesh[i].nb.end())
						{
							tmesh[i].nb.push_back(fcid);
							int *it1 = find(tmesh[fcid].edge, tmesh[fcid].edge + 4, tmesh[i].edge[j]);
							if (it1 != tmesh[fcid].edge + 4)
							{
								int pos = it1 - tmesh[fcid].edge;
								if ((j + 2) % 4 == pos)
									tmesh[i].nbrot.push_back(0);
								if ((j + 1) % 4 == pos)
									tmesh[i].nbrot.push_back(1);
								if (j == pos)
									tmesh[i].nbrot.push_back(2);
								if ((j + 3) % 4 == pos)
									tmesh[i].nbrot.push_back(3);
							}
							else
							{
								int pos(0);
								for (int k1 = 0; k1 < 4; k1++)
								{
									if (tmedge[tmesh[fcid].edge[k]].act == 0 && (tmedge[tmesh[fcid].edge[k1]].chd[0] == tmesh[i].edge[j] || tmedge[tmesh[fcid].edge[k1]].chd[1] == tmesh[i].edge[j]))
									{
										pos = k1;
										break;
									}
								}
								if ((j + 2) % 4 == pos)
									tmesh[i].nbrot.push_back(0);
								if ((j + 1) % 4 == pos)
									tmesh[i].nbrot.push_back(1);
								if (j == pos)
									tmesh[i].nbrot.push_back(2);
								if ((j + 3) % 4 == pos)
									tmesh[i].nbrot.push_back(3);
							}
						}
					}
				}
			}

			//corner neighbor
			for (j = 0; j < 1; j++)
			{
				if (cp[tmesh[i].cnct[j]].type != 2)
				{
					for (k = 0; k < cp[tmesh[i].cnct[j]].face.size(); k++)
					{
						int fcid(cp[tmesh[i].cnct[j]].face[k]);
						vector<int>::iterator it = find(tmesh[i].nb.begin(), tmesh[i].nb.end(), fcid);
						if (it == tmesh[i].nb.end())
						{
							tmesh[i].nb.push_back(fcid);
							int *it1 = find(tmesh[fcid].cnct, tmesh[fcid].cnct + 4, tmesh[i].cnct[j]);
							int pos = it1 - tmesh[fcid].cnct;
							int caseid(0);
							for (int k1 = 0; k1 < 4; k1++)
							{
								if (rotcase[k1][j] == pos)
								{
									caseid = k1;
									break;
								}
							}
							tmesh[i].nbrot.push_back(caseid);
						}
					}
				}
			}

			//relative origin
			for (j = 0; j < tmesh[i].nb.size(); j++)
			{
				int fid(tmesh[i].nb[j]);
				double pmc[4][2] = {{0., 0.}, {tmedge[tmesh[fid].edge[0]].len, 0.}, {tmedge[tmesh[fid].edge[0]].len, tmedge[tmesh[fid].edge[1]].len}, {0., tmedge[tmesh[fid].edge[1]].len}};
				double pmt[4][2] = {{tmedge[tmesh[fid].edge[0]].len / 2., 0.}, {tmedge[tmesh[fid].edge[0]].len, tmedge[tmesh[fid].edge[1]].len / 2.}, {tmedge[tmesh[fid].edge[0]].len / 2., tmedge[tmesh[fid].edge[1]].len}, {0., tmedge[tmesh[fid].edge[1]].len / 2.}};
			}
		}
	}

	//relative origin
	for (i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			double pmc[4][2] = {{0., 0.}, {tmedge[tmesh[i].edge[0]].len, 0.}, {tmedge[tmesh[i].edge[0]].len, tmedge[tmesh[i].edge[1]].len}, {0., tmedge[tmesh[i].edge[1]].len}};
			//double pmt[4][2]={{tmedge[tmesh[i].edge[0]].len/2.,0.},{tmedge[tmesh[i].edge[0]].len,tmedge[tmesh[i].edge[1]].len/2.},
			//{tmedge[tmesh[i].edge[0]].len/2.,tmedge[tmesh[i].edge[1]].len},{0.,tmedge[tmesh[i].edge[1]].len/2.}};
			for (j = 0; j < 4; j++)
			{
				for (k = 0; k < cp[tmesh[i].cnct[j]].face.size(); k++)
				{
					int fid(cp[tmesh[i].cnct[j]].face[k]);
					if (fid != i)
					{
						vector<int>::iterator it = find(tmesh[i].nb.begin(), tmesh[i].nb.end(), fid);
						int pos = it - tmesh[i].nb.begin();
						double pmc1[4][2] = {{0., 0.}, {tmedge[tmesh[fid].edge[0]].len, 0.}, {tmedge[tmesh[fid].edge[0]].len, tmedge[tmesh[fid].edge[1]].len}, {0., tmedge[tmesh[fid].edge[1]].len}};
						double pmt1[4][2] = {{tmedge[tmesh[fid].edge[0]].len / 2., 0.}, {tmedge[tmesh[fid].edge[0]].len, tmedge[tmesh[fid].edge[1]].len / 2.}, {tmedge[tmesh[fid].edge[0]].len / 2., tmedge[tmesh[fid].edge[1]].len}, {0., tmedge[tmesh[fid].edge[1]].len / 2.}};
						int *it1 = find(tmesh[fid].cnct, tmesh[fid].cnct + 4, tmesh[i].cnct[j]);
						if (it1 != tmesh[fid].cnct + 4)
						{
							int loc = it1 - tmesh[fid].cnct;
						}
						else
						{
							int loc(0);
							for (int k2 = 0; k2 < 4; k2++)
							{
								if (tmedge[tmesh[fid].edge[k2]].act == 0 && tmedge[tmesh[fid].edge[k2]].midpt == tmesh[i].cnct[j])
								{
									loc = k2;
									break;
								}
							}
						}
					}
				}
			}
		}
	}
}

void TruncatedTspline::ElementKnotVectors(int eid)
{
	uint i;
	for (i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		int pid(tmesh[eid].IEN[i]);
		if (tmesh[eid].UV[i] == 0)
		{
			cp[pid].knotU[0] = tmesh[i].pmcoor[i][0] - cp[pid].kitvU[1] - cp[pid].kitvU[0];
			cp[pid].knotU[1] = tmesh[i].pmcoor[i][0] - cp[pid].kitvU[1];
			cp[pid].knotU[2] = tmesh[i].pmcoor[i][0];
			cp[pid].knotU[3] = tmesh[i].pmcoor[i][0] + cp[pid].kitvU[2];
			cp[pid].knotU[4] = tmesh[i].pmcoor[i][0] + cp[pid].kitvU[2] + cp[pid].kitvU[3];
			cp[pid].knotV[0] = tmesh[i].pmcoor[i][1] - cp[pid].kitvV[1] - cp[pid].kitvV[0];
			cp[pid].knotV[1] = tmesh[i].pmcoor[i][1] - cp[pid].kitvV[1];
			cp[pid].knotV[2] = tmesh[i].pmcoor[i][1];
			cp[pid].knotV[3] = tmesh[i].pmcoor[i][1] + cp[pid].kitvV[2];
			cp[pid].knotV[4] = tmesh[i].pmcoor[i][1] + cp[pid].kitvV[2] + cp[pid].kitvV[3];
		}
		else if (tmesh[eid].UV[i] == 1)
		{
			cp[pid].knotU[0] = tmesh[i].pmcoor[i][0] - cp[pid].kitvV[2] - cp[pid].kitvV[3];
			cp[pid].knotU[1] = tmesh[i].pmcoor[i][0] - cp[pid].kitvV[2];
			cp[pid].knotU[2] = tmesh[i].pmcoor[i][0];
			cp[pid].knotU[3] = tmesh[i].pmcoor[i][0] + cp[pid].kitvV[1];
			cp[pid].knotU[4] = tmesh[i].pmcoor[i][0] + cp[pid].kitvV[1] + cp[pid].kitvV[0];
			cp[pid].knotV[0] = tmesh[i].pmcoor[i][1] - cp[pid].kitvU[1] - cp[pid].kitvU[0];
			cp[pid].knotV[1] = tmesh[i].pmcoor[i][1] - cp[pid].kitvU[1];
			cp[pid].knotV[2] = tmesh[i].pmcoor[i][1];
			cp[pid].knotV[3] = tmesh[i].pmcoor[i][1] + cp[pid].kitvU[2];
			cp[pid].knotV[4] = tmesh[i].pmcoor[i][1] + cp[pid].kitvU[2] + cp[pid].kitvU[3];
		}
		else if (tmesh[eid].UV[i] == 2)
		{
			cp[pid].knotU[0] = tmesh[i].pmcoor[i][0] - cp[pid].kitvU[2] - cp[pid].kitvU[3];
			cp[pid].knotU[1] = tmesh[i].pmcoor[i][0] - cp[pid].kitvU[2];
			cp[pid].knotU[2] = tmesh[i].pmcoor[i][0];
			cp[pid].knotU[3] = tmesh[i].pmcoor[i][0] + cp[pid].kitvU[1];
			cp[pid].knotU[4] = tmesh[i].pmcoor[i][0] + cp[pid].kitvU[1] + cp[pid].kitvU[0];
			cp[pid].knotV[0] = tmesh[i].pmcoor[i][1] - cp[pid].kitvV[2] - cp[pid].kitvV[3];
			cp[pid].knotV[1] = tmesh[i].pmcoor[i][1] - cp[pid].kitvV[2];
			cp[pid].knotV[2] = tmesh[i].pmcoor[i][1];
			cp[pid].knotV[3] = tmesh[i].pmcoor[i][1] + cp[pid].kitvV[1];
			cp[pid].knotV[4] = tmesh[i].pmcoor[i][1] + cp[pid].kitvV[1] + cp[pid].kitvV[0];
		}
		else
		{
			cp[pid].knotU[0] = tmesh[i].pmcoor[i][0] - cp[pid].kitvV[1] - cp[pid].kitvV[0];
			cp[pid].knotU[1] = tmesh[i].pmcoor[i][0] - cp[pid].kitvV[1];
			cp[pid].knotU[2] = tmesh[i].pmcoor[i][0];
			cp[pid].knotU[3] = tmesh[i].pmcoor[i][0] + cp[pid].kitvV[2];
			cp[pid].knotU[4] = tmesh[i].pmcoor[i][0] + cp[pid].kitvV[2] + cp[pid].kitvV[3];
			cp[pid].knotV[0] = tmesh[i].pmcoor[i][1] - cp[pid].kitvU[2] - cp[pid].kitvU[3];
			cp[pid].knotV[1] = tmesh[i].pmcoor[i][1] - cp[pid].kitvU[2];
			cp[pid].knotV[2] = tmesh[i].pmcoor[i][1];
			cp[pid].knotV[3] = tmesh[i].pmcoor[i][1] + cp[pid].kitvU[1];
			cp[pid].knotV[4] = tmesh[i].pmcoor[i][1] + cp[pid].kitvU[1] + cp[pid].kitvU[0];
		}
	}
}

void TruncatedTspline::AssignIndex_NewFacePoint(int fcid, Vertex &pt)
{
	int flag(0);
	for (uint i = 0; i < uanc.size(); i++)
	{
		if (pt.pm[0] == uanc[i].second)
		{
			flag = 1;
			pt.index[0] = uanc[i].first;
			break;
		}
	}
	if (flag == 0)
	{
		pt.index[0] = uanc.size();
		pair<int, double> upair(pt.index[0], pt.pm[0]);
		uanc.push_back(upair);
	}
	flag = 0;
	for (uint i = 0; i < vanc.size(); i++)
	{
		if (pt.pm[1] == vanc[i].second)
		{
			flag = 1;
			pt.index[1] = vanc[i].first;
			break;
		}
	}
	if (flag == 0)
	{
		pt.index[1] = vanc.size();
		pair<int, double> vpair(pt.index[1], pt.pm[1]);
		vanc.push_back(vpair);
	}
}

void TruncatedTspline::AssignIndex_NewEdgePoint(int edid, Vertex &pt)
{
	if (pt.pm[0] != cp[tmedge[edid].pt[0]].pm[0] && pt.pm[0] != cp[tmedge[edid].pt[1]].pm[0])
	{
		int flag(0);
		for (uint i = 0; i < uanc.size(); i++)
		{
			if (pt.pm[0] == uanc[i].second)
			{
				flag = 1;
				pt.index[0] = uanc[i].first;
				break;
			}
		}
		if (flag == 0)
		{
			pt.index[0] = uanc.size();
			pair<int, double> upair(pt.index[0], pt.pm[0]);
			uanc.push_back(upair);
		}
		pt.index[1] = cp[tmedge[edid].pt[0]].index[1];
	}
	else
	{
		int flag = 0;
		for (uint i = 0; i < vanc.size(); i++)
		{
			if (pt.pm[1] == vanc[i].second)
			{
				flag = 1;
				pt.index[1] = vanc[i].first;
				break;
			}
		}
		if (flag == 0)
		{
			pt.index[1] = vanc.size();
			pair<int, double> vpair(pt.index[1], pt.pm[1]);
			vanc.push_back(vpair);
		}
		pt.index[0] = cp[tmedge[edid].pt[0]].index[0];
	}
}

void TruncatedTspline::SortEdge()
{
	//vector<double> uanc,vanc;
	//vector<vector<int>> uedge,vedge;
	//uanc.clear();
	//vanc.clear();
	for (uint i = 0; i < uedge.size(); i++)
	{
		uedge[i].clear();
	}
	for (uint i = 0; i < vedge.size(); i++)
	{
		vedge[i].clear();
	}
	uedge.clear();
	vedge.clear();
	//uanc.push_back(0.);//boundary
	//vanc.push_back(0.);//boundary
	//for(uint i=0; i<cp.size(); i++)
	//{
	//	pair<int,double> upair(cp[i].index[0],cp[i].pm[0]);
	//	pair<int,double> vpair(cp[i].index[1],cp[i].pm[1]);
	//	vector<pair<int,double>>::iterator it1=find(uanc.begin(),uanc.end(),upair);
	//	vector<pair<int,double>>::iterator it2=find(vanc.begin(),vanc.end(),vpair);
	//	if(it1==uanc.end()) uanc.push_back(upair);
	//	if(it2==vanc.end()) vanc.push_back(vpair);
	//}
	sort(uanc.begin(), uanc.end(), CompareAnchor);
	sort(vanc.begin(), vanc.end(), CompareAnchor);
	for (uint i = 0; i < uanc.size() - 1; i++)
	{
		if (uanc[i].second == uanc[i + 1].second && uanc[i].first > uanc[i].first)
		{
			int tmp = uanc[i].first;
			uanc[i].first = uanc[i + 1].first;
			uanc[i + 1].first = tmp;
		}
	}
	for (uint i = 0; i < vanc.size() - 1; i++)
	{
		if (vanc[i].second == vanc[i + 1].second && vanc[i].first > vanc[i].first)
		{
			int tmp = vanc[i].first;
			vanc[i].first = vanc[i + 1].first;
			vanc[i + 1].first = tmp;
		}
	}
	//for(uint i=0; i<uanc.size(); i++)
	//{
	//	cout<<uanc[i].first<<" ";
	//}
	//cout<<"\n";
	//for(uint i=0; i<uanc.size(); i++)
	//{
	//	cout<<uanc[i].second<<" ";
	//}
	//cout<<"\n";
	//for(uint i=0; i<vanc.size(); i++)
	//{
	//	cout<<vanc[i].first<<" ";
	//}
	//cout<<"\n";
	//for(uint i=0; i<vanc.size(); i++)
	//{
	//	cout<<vanc[i].second<<" ";
	//}
	//cout<<"\n";
	//getchar();

	/*if(uanc[0].first > uanc[1].first)
	{
		int tmp=uanc[0].first;
		uanc[0].first=uanc[1].first;
		uanc[1].first=tmp;
	}
	if(vanc[0].first > vanc[1].first)
	{
		int tmp=vanc[0].first;
		vanc[0].first=vanc[1].first;
		vanc[1].first=tmp;
	}
	if(uanc[uanc.size()-2].first > uanc[uanc.size()-1].first)
	{
		int tmp=uanc[uanc.size()-2].first;
		uanc[uanc.size()-2].first=uanc[uanc.size()-1].first;
		uanc[uanc.size()-1].first=tmp;
	}
	if(vanc[vanc.size()-2].first > vanc[vanc.size()-1].first)
	{
		int tmp=vanc[vanc.size()-2].first;
		vanc[vanc.size()-2].first=vanc[vanc.size()-1].first;
		vanc[vanc.size()-1].first=tmp;
	}*/
	//uanc.push_back(uanc.back());//boundary
	//vanc.push_back(vanc.back());//boundary
	uedge.resize(uanc.size());
	vedge.resize(vanc.size());

	for (uint i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].act == 1)
		{
			if (cp[tmedge[i].pt[0]].index[0] == cp[tmedge[i].pt[1]].index[0])
			{
				pair<int, double> upair(cp[tmedge[i].pt[0]].index[0], cp[tmedge[i].pt[0]].pm[0]);
				vector<pair<int, double>>::iterator it = find(uanc.begin(), uanc.end(), upair);
				int pos = it - uanc.begin();
				uedge[pos].push_back(i);
			}
			else
			{
				pair<int, double> vpair(cp[tmedge[i].pt[0]].index[1], cp[tmedge[i].pt[0]].pm[1]);
				vector<pair<int, double>>::iterator it = find(vanc.begin(), vanc.end(), vpair);
				int pos = it - vanc.begin();
				//if(it==vanc.end())
				//{
				//	cout<<tmedge[i].pt[0]<<'\n';
				//	for(uint i0=0; i0<vanc.size(); i0++)
				//	{
				//		cout<<vanc[i0].first<<' ';
				//	}
				//	cout<<'\n';
				//	for(uint i0=0; i0<vanc.size(); i0++)
				//	{
				//		cout<<vanc[i0].second<<' ';
				//	}
				//	cout<<'\n';
				//	cout<<vpair.first<<' '<<vpair.second<<'\n';
				//	getchar();
				//}
				vedge[pos].push_back(i);
			}

			//if(tmedge[i].face.size()==2)//inner edge
			//{
			//	if(cp[tmedge[i].pt[0]].pm[0]==cp[tmedge[i].pt[1]].pm[0])//uedge, same u
			//	{
			//		vector<double>::iterator it=find(uanc.begin(),uanc.end(),cp[tmedge[i].pt[0]].pm[0]);
			//		int pos=it-uanc.begin();
			//		uedge[pos].push_back(i);
			//	}
			//	else
			//	{
			//		vector<double>::iterator it=find(vanc.begin(),vanc.end(),cp[tmedge[i].pt[0]].pm[1]);
			//		int pos=it-vanc.begin();
			//		vedge[pos].push_back(i);
			//	}
			//}
			//else//boundary edge
			//{
			//	if(cp[tmedge[i].pt[0]].pm[0]==cp[tmedge[i].pt[1]].pm[0])//uedge
			//	{
			//		if(cp[tmedge[i].pt[0]].pm[0]==uanc[1])
			//		{
			//			uedge[0].push_back(i);
			//		}
			//		else
			//		{
			//			uedge[uedge.size()-1].push_back(i);
			//		}
			//	}
			//	else
			//	{
			//		if(cp[tmedge[i].pt[0]].pm[1]==vanc[1])
			//		{
			//			vedge[0].push_back(i);
			//		}
			//		else
			//		{
			//			vedge[vedge.size()-1].push_back(i);
			//		}
			//	}
			//}
		}
	}
}

void TruncatedTspline::FindLocalKnotVectors()
{
	//cout<<"before sort!\n";
	//for(uint i=0; i<vanc.size(); i++)
	//{
	//	cout<<vanc[i].first<<' '<<vanc[i].second<<'\n';
	//}
	//cout<<"sort\n";
	SortEdge();
	//cout<<"after sort!\n";
	//for(uint i=0; i<vanc.size(); i++)
	//{
	//	cout<<vanc[i].first<<' '<<vanc[i].second<<'\n';
	//}
	//cout<<'\n';
	for (uint i = 0; i < cp.size(); i++)
	{
		cp[i].kutmp[2] = cp[i].pm[0];
		cp[i].kvtmp[2] = cp[i].pm[1];
		//u direction
		pair<int, double> pu(cp[i].index[0], cp[i].pm[0]);
		vector<pair<int, double>>::iterator it = find(uanc.begin(), uanc.end(), pu);
		int pos = it - uanc.begin();
		int count(1), loc(pos);
		while (count < 3) //two intersections along +u
		{
			loc++;
			if (loc < uanc.size())
			{
				for (uint j = 0; j < uedge[loc].size(); j++)
				{
					int eid(uedge[loc][j]);
					double edpm[2] = {cp[tmedge[eid].pt[0]].pm[1], cp[tmedge[eid].pt[1]].pm[1]};
					sort(edpm, edpm + 2);
					if (cp[i].pm[1] >= edpm[0] && cp[i].pm[1] <= edpm[1])
					{
						cp[i].kutmp[2 + count] = uanc[loc].second;
						//if(i==118) cout<<cp[i].kutmp[2+count]<<'\n';
						count++;
						break;
					}
				}
			}
			else
			{
				cp[i].kutmp[2 + count] = uanc[uanc.size() - 1].second;
				//if(i==118) cout<<cp[i].kutmp[2+count]<<'\n';
				count++;
			}
		}
		count = 1;
		loc = pos;
		while (count < 3) //two intersections along -u
		{
			loc--;
			if (loc > 0)
			{
				for (uint j = 0; j < uedge[loc].size(); j++)
				{
					int eid(uedge[loc][j]);
					double edpm[2] = {cp[tmedge[eid].pt[0]].pm[1], cp[tmedge[eid].pt[1]].pm[1]};
					sort(edpm, edpm + 2);
					if (cp[i].pm[1] >= edpm[0] && cp[i].pm[1] <= edpm[1])
					{
						cp[i].kutmp[2 - count] = uanc[loc].second;
						//if(i==118) cout<<cp[i].kutmp[2-count]<<'\n';
						count++;
						break;
					}
				}
			}
			else
			{
				cp[i].kutmp[2 - count] = uanc[0].second;
				//if(i==118) cout<<cp[i].kutmp[2-count]<<'\n';
				count++;
			}
		}
		//v direction
		pair<int, double> pv(cp[i].index[1], cp[i].pm[1]);
		it = find(vanc.begin(), vanc.end(), pv);
		pos = it - vanc.begin();
		//if((cp[i].face.size()==1 || cp[i].face.size()==2) && pos==1)
		//	pos=0;
		//if((cp[i].face.size()==1 || cp[i].face.size()==2) && pos==vanc.size()-2)
		//	pos=vanc.size()-1;
		count = 1;
		loc = pos;
		while (count < 3) //two intersections along +v
		{
			loc++;
			if (loc < vanc.size())
			{
				for (uint j = 0; j < vedge[loc].size(); j++)
				{
					int eid(vedge[loc][j]);
					double edpm[2] = {cp[tmedge[eid].pt[0]].pm[0], cp[tmedge[eid].pt[1]].pm[0]};
					sort(edpm, edpm + 2);
					if (cp[i].pm[0] >= edpm[0] && cp[i].pm[0] <= edpm[1])
					{
						cp[i].kvtmp[2 + count] = vanc[loc].second;
						//if(i==118) cout<<cp[i].kvtmp[2+count]<<'\n';
						count++;
						break;
					}
				}
			}
			else
			{
				cp[i].kvtmp[2 + count] = vanc[vanc.size() - 1].second;
				//if(i==118) cout<<cp[i].kvtmp[2+count]<<'\n';
				count++;
			}
		}
		count = 1;
		loc = pos;
		while (count < 3) //two intersections along -v
		{
			loc--;
			if (loc > 0)
			{
				for (uint j = 0; j < vedge[loc].size(); j++)
				{
					int eid(vedge[loc][j]);
					double edpm[2] = {cp[tmedge[eid].pt[0]].pm[0], cp[tmedge[eid].pt[1]].pm[0]};
					sort(edpm, edpm + 2);
					if (cp[i].pm[0] >= edpm[0] && cp[i].pm[0] <= edpm[1])
					{
						cp[i].kvtmp[2 - count] = vanc[loc].second;
						//if(i==118)
						//{
						//	cout<<"knot: "<<cp[i].kvtmp[2-count]<<'\n';
						//	cout<<"loc: "<<loc<<'\n';
						//	cout<<"index: "<<vanc[loc].first<<'\n';
						//	cout<<'\n';
						//}
						count++;
						break;
					}
				}
			}
			else
			{
				cp[i].kvtmp[2 - count] = vanc[0].second;
				if (i == 118)
					cout << cp[i].kvtmp[2 - count] << '\n';
				count++;
			}
		}

		//cout<<cp[i].knotU[0]<<" "<<cp[i].knotU[1]<<" "<<cp[i].knotU[2]<<" "<<cp[i].knotU[3]<<" "<<cp[i].knotU[4]<<"\n";
		//cout<<cp[i].knotV[0]<<" "<<cp[i].knotV[1]<<" "<<cp[i].knotV[2]<<" "<<cp[i].knotV[3]<<" "<<cp[i].knotV[4]<<"\n";
		//getchar();
	}
}

void TruncatedTspline::FindIENglb()
{
	for (uint i = 0; i < tmesh.size(); i++)
	{
		tmesh[i].IEN.clear();
	}

	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			for (uint j = 0; j < cp.size(); j++)
			{
				if (cp[j].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[j].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[j].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[j].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(j);
				}
			}
		}
	}
}

void TruncatedTspline::UpdateControlPoints() //calculate control points
{
	for (uint i = 0; i < npt_old; i++)
	{
		cp[i].aff = 0;
		array<double, 10> kvnew = {cp[i].kutmp[0], cp[i].kutmp[1], cp[i].kutmp[2], cp[i].kutmp[3], cp[i].kutmp[4], cp[i].kvtmp[0], cp[i].kvtmp[1], cp[i].kvtmp[2], cp[i].kvtmp[3], cp[i].kvtmp[4]};
		array<double, 10> kvold = {cp[i].knotU[0], cp[i].knotU[1], cp[i].knotU[2], cp[i].knotU[3], cp[i].knotU[4], cp[i].knotV[0], cp[i].knotV[1], cp[i].knotV[2], cp[i].knotV[3], cp[i].knotV[4]};
		if (kvnew != kvold)
			cp[i].aff = 1;
	}
	for (uint i = 0; i < cp.size(); i++) //to be updated
	{
		cp[i].coortmp[0] = 0.;
		cp[i].coortmp[1] = 0.;
		cp[i].coortmp[2] = 0.;
		for (uint j = 0; j < npt_old; j++) //old
		{
			if (cp[i].kutmp[0] >= cp[j].knotU[0] && cp[i].kutmp[4] <= cp[j].knotU[4] && cp[i].kvtmp[0] >= cp[j].knotV[0] && cp[i].kvtmp[4] <= cp[j].knotV[4])
			{
				array<double, 10> kvnew = {cp[i].kutmp[0], cp[i].kutmp[1], cp[i].kutmp[2], cp[i].kutmp[3], cp[i].kutmp[4], cp[i].kvtmp[0], cp[i].kvtmp[1], cp[i].kvtmp[2], cp[i].kvtmp[3], cp[i].kvtmp[4]};
				array<double, 10> kvold = {cp[j].knotU[0], cp[j].knotU[1], cp[j].knotU[2], cp[j].knotU[3], cp[j].knotU[4], cp[j].knotV[0], cp[j].knotV[1], cp[j].knotV[2], cp[j].knotV[3], cp[j].knotV[4]};
				if (kvnew != kvold)
				{
					vector<double> ku1(10), kv1(10);
					vector<vector<double>> Tu, Tv;
					vector<double> ku(cp[j].knotU, cp[j].knotU + 5);
					vector<double> kv(cp[j].knotV, cp[j].knotV + 5);
					vector<double>::iterator it1, it2;
					it1 = set_union(cp[j].knotU, cp[j].knotU + 5, cp[i].kutmp, cp[i].kutmp + 5, ku1.begin());
					it2 = set_union(cp[j].knotV, cp[j].knotV + 5, cp[i].kvtmp, cp[i].kvtmp + 5, kv1.begin());
					ku1.resize(it1 - ku1.begin());
					kv1.resize(it2 - kv1.begin());
					TMatrix(ku, ku1, 3, Tu);
					TMatrix(kv, kv1, 3, Tv);
					it1 = search(ku1.begin(), ku1.end(), cp[i].kutmp, cp[i].kutmp + 5);
					it2 = search(kv1.begin(), kv1.end(), cp[i].kvtmp, cp[i].kvtmp + 5);
					int loc1 = it1 - ku1.begin();
					int loc2 = it2 - kv1.begin();
					double coef = Tu[loc1][0] * Tv[loc2][0];
					if (i < npt_old)
					{
						cp[i].coortmp[0] += coef * cp[j].coor[0];
						cp[i].coortmp[1] += coef * cp[j].coor[1];
						cp[i].coortmp[2] += coef * cp[j].coor[2];
					}
					else
					{
						cp[i].coor[0] += coef * cp[j].coor[0];
						cp[i].coor[1] += coef * cp[j].coor[1];
						cp[i].coor[2] += coef * cp[j].coor[2];
					}
				}
			}
		}
	}
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].aff == 1)
		{
			cp[i].coor[0] = cp[i].coortmp[0];
			cp[i].coor[1] = cp[i].coortmp[1];
			cp[i].coor[2] = cp[i].coortmp[2];
		}
	}
}

void TruncatedTspline::UpdateControlPoints_1() //calculate control points
{
	double tol(1.e-6);
	vector<vector<double>> cmat(cp.size(), vector<double>(npt_old, 0.));
	//for(uint i=0; i<npt_old; i++)
	//{
	//	cp[i].aff=0;
	//	array<double,10> kvnew={cp[i].kutmp[0],cp[i].kutmp[1],cp[i].kutmp[2],cp[i].kutmp[3],cp[i].kutmp[4],cp[i].kvtmp[0],cp[i].kvtmp[1],cp[i].kvtmp[2],cp[i].kvtmp[3],cp[i].kvtmp[4]};
	//	array<double,10> kvold={cp[i].knotU[0],cp[i].knotU[1],cp[i].knotU[2],cp[i].knotU[3],cp[i].knotU[4],cp[i].knotV[0],cp[i].knotV[1],cp[i].knotV[2],cp[i].knotV[3],cp[i].knotV[4]};
	//	if(kvnew!=kvold)
	//		cp[i].aff=1;
	//}
	for (uint i = 0; i < cp.size(); i++) //to be updated
	{
		cp[i].coortmp[0] = 0.;
		cp[i].coortmp[1] = 0.;
		cp[i].coortmp[2] = 0.;
		double sum(0.);
		for (uint j = 0; j < npt_old; j++) //old
		{
			if (cp[i].kutmp[0] >= cp[j].knotU[0] && cp[i].kutmp[4] <= cp[j].knotU[4] && cp[i].kvtmp[0] >= cp[j].knotV[0] && cp[i].kvtmp[4] <= cp[j].knotV[4])
			{
				array<double, 10> kvnew = {cp[i].kutmp[0], cp[i].kutmp[1], cp[i].kutmp[2], cp[i].kutmp[3], cp[i].kutmp[4], cp[i].kvtmp[0], cp[i].kvtmp[1], cp[i].kvtmp[2], cp[i].kvtmp[3], cp[i].kvtmp[4]};
				array<double, 10> kvold = {cp[j].knotU[0], cp[j].knotU[1], cp[j].knotU[2], cp[j].knotU[3], cp[j].knotU[4], cp[j].knotV[0], cp[j].knotV[1], cp[j].knotV[2], cp[j].knotV[3], cp[j].knotV[4]};
				if (kvnew != kvold)
				{
					vector<double> ku1(10), kv1(10);
					vector<vector<double>> Tu, Tv;
					vector<double> ku(cp[j].knotU, cp[j].knotU + 5);
					vector<double> kv(cp[j].knotV, cp[j].knotV + 5);
					vector<double>::iterator it1, it2;
					it1 = set_union(cp[j].knotU, cp[j].knotU + 5, cp[i].kutmp, cp[i].kutmp + 5, ku1.begin());
					it2 = set_union(cp[j].knotV, cp[j].knotV + 5, cp[i].kvtmp, cp[i].kvtmp + 5, kv1.begin());
					ku1.resize(it1 - ku1.begin());
					kv1.resize(it2 - kv1.begin());
					TMatrix(ku, ku1, 3, Tu);
					TMatrix(kv, kv1, 3, Tv);
					it1 = search(ku1.begin(), ku1.end(), cp[i].kutmp, cp[i].kutmp + 5);
					it2 = search(kv1.begin(), kv1.end(), cp[i].kvtmp, cp[i].kvtmp + 5);
					int loc1 = it1 - ku1.begin();
					int loc2 = it2 - kv1.begin();
					double coef = Tu[loc1][0] * Tv[loc2][0];
					cmat[i][j] = coef;
					sum += coef;
					//if(i<npt_old)
					//{
					//	cp[i].coortmp[0]+=coef*cp[j].coor[0]; cp[i].coortmp[1]+=coef*cp[j].coor[1]; cp[i].coortmp[2]+=coef*cp[j].coor[2];
					//}
					//else
					//{
					//	cp[i].coor[0]+=coef*cp[j].coor[0]; cp[i].coor[1]+=coef*cp[j].coor[1]; cp[i].coor[2]+=coef*cp[j].coor[2];
					//}
				}
				else
				{
					cmat[i][j] = 1.;
					sum += 1.;
				}
			}
		}
		sum -= 1.;
		if (sum < 0.)
			sum = -sum;
		if (sum < tol)
		{
			cp[i].trun = 0;
			for (uint j = 0; j < npt_old; j++)
			{
				if (cmat[i][j] != 0.)
				{
					cp[i].coortmp[0] += cmat[i][j] * cp[j].coor[0];
					cp[i].coortmp[1] += cmat[i][j] * cp[j].coor[1];
					cp[i].coortmp[2] += cmat[i][j] * cp[j].coor[2];
				}
			}
		}
		else
		{
			cp[i].trun = 1;
			if (i >= npt_old)
			{
				cout << "New points are not weighted average!\n";
				getchar();
			}
		}
	}
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].trun == 0)
		{
			cp[i].coor[0] = cp[i].coortmp[0];
			cp[i].coor[1] = cp[i].coortmp[1];
			cp[i].coor[2] = cp[i].coortmp[2];
			for (uint j = 0; j < 5; j++)
			{
				cp[i].knotU[j] = cp[i].kutmp[j];
				cp[i].knotV[j] = cp[i].kvtmp[j];
			}
		}
	}
	for (uint j = 0; j < npt_old; j++)
	{
		if (cp[j].trun == 1)
		{
			for (uint i = 0; i < cp.size(); i++)
			{
				if (cmat[i][j] != 0. && cp[i].trun == 0)
				{
					cp[j].tbf.push_back(i);
					cp[j].tc.push_back(cmat[i][j]);
				}
			}
		}
		else if (cp[j].trun == 0)
		{
			for (uint i = 0; i < cp.size(); i++)
			{
				if (i != j && cmat[i][j] != 0. && cp[i].trun == 0 && cp[i].knotU[0] >= cp[j].knotU[0] && cp[i].knotU[4] <= cp[j].knotU[4] && cp[i].knotV[0] >= cp[j].knotV[0] && cp[i].knotV[4] <= cp[j].knotV[4])
				{
					cp[j].trun = 1;
					cp[j].tbf.push_back(i);
					cp[j].tc.push_back(cmat[i][j]);
				}
			}
		}
	}
}

void TruncatedTspline::UpdateControlPoints_2() //calculate control points, this
{
	for (uint i = 0; i < npt_old; i++)
	{
		//cp[i].trun=0;
		//cp[i].tbf.clear();
		//cp[i].tc.clear();
		cp[i].aff = 0;
		//array<double,10> kvnew={cp[i].kutmp[0],cp[i].kutmp[1],cp[i].kutmp[2],cp[i].kutmp[3],cp[i].kutmp[4],cp[i].kvtmp[0],cp[i].kvtmp[1],cp[i].kvtmp[2],cp[i].kvtmp[3],cp[i].kvtmp[4]};
		//array<double,10> kvold={cp[i].knotU[0],cp[i].knotU[1],cp[i].knotU[2],cp[i].knotU[3],cp[i].knotU[4],cp[i].knotV[0],cp[i].knotV[1],cp[i].knotV[2],cp[i].knotV[3],cp[i].knotV[4]};
		if (CheckSubKnotVector(cp[i].kutmp, cp[i].kvtmp, cp[i].knotU, cp[i].knotV))
			cp[i].aff = 1;
	}
	for (uint i = 0; i < cp.size(); i++) //to be updated
	{
		cp[i].coortmp[0] = 0.;
		cp[i].coortmp[1] = 0.;
		cp[i].coortmp[2] = 0.;
		if (cp[i].aff == 1 || i >= npt_old)
		{
			for (uint j = 0; j < npt_old; j++) //old
			{
				//if(cp[i].kutmp[0]>=cp[j].knotU[0] && cp[i].kutmp[4]<=cp[j].knotU[4] && cp[i].kvtmp[0]>=cp[j].knotV[0] && cp[i].kvtmp[4]<=cp[j].knotV[4])
				if (CheckSubKnotVector(cp[i].kutmp, cp[i].kvtmp, cp[j].knotU, cp[j].knotV))
				{
					//array<double,10> kvnew={cp[i].kutmp[0],cp[i].kutmp[1],cp[i].kutmp[2],cp[i].kutmp[3],cp[i].kutmp[4],cp[i].kvtmp[0],cp[i].kvtmp[1],cp[i].kvtmp[2],cp[i].kvtmp[3],cp[i].kvtmp[4]};
					//array<double,10> kvold={cp[j].knotU[0],cp[j].knotU[1],cp[j].knotU[2],cp[j].knotU[3],cp[j].knotU[4],cp[j].knotV[0],cp[j].knotV[1],cp[j].knotV[2],cp[j].knotV[3],cp[j].knotV[4]};
					//if(kvnew!=kvold)
					{
						vector<double> ku1(10), kv1(10);
						vector<vector<double>> Tu, Tv;
						vector<double> ku(cp[j].knotU, cp[j].knotU + 5);
						vector<double> kv(cp[j].knotV, cp[j].knotV + 5);
						vector<double>::iterator it1, it2;
						it1 = set_union(cp[j].knotU, cp[j].knotU + 5, cp[i].kutmp, cp[i].kutmp + 5, ku1.begin());
						it2 = set_union(cp[j].knotV, cp[j].knotV + 5, cp[i].kvtmp, cp[i].kvtmp + 5, kv1.begin());
						ku1.resize(it1 - ku1.begin());
						kv1.resize(it2 - kv1.begin());
						TMatrix(ku, ku1, 3, Tu);
						TMatrix(kv, kv1, 3, Tv);
						it1 = search(ku1.begin(), ku1.end(), cp[i].kutmp, cp[i].kutmp + 5);
						it2 = search(kv1.begin(), kv1.end(), cp[i].kvtmp, cp[i].kvtmp + 5);
						int loc1 = it1 - ku1.begin();
						int loc2 = it2 - kv1.begin();
						/*if(it1==ku1.end() || it2==kv1.end())
						{
							cout<<"kutmp, kvtmp: "<<'\n';
							for(uint k=0; k<5; k++)
							{
								cout<<cp[i].kutmp[k]<<' ';
							}
							cout<<'\n';
							for(uint k=0; k<5; k++)
							{
								cout<<cp[i].kvtmp[k]<<' ';
							}
							cout<<'\n';
							cout<<"ku, kv: "<<'\n';
							for(uint k=0; k<5; k++)
							{
								cout<<cp[j].knotU[k]<<' ';
							}
							cout<<'\n';
							for(uint k=0; k<5; k++)
							{
								cout<<cp[j].knotV[k]<<' ';
							}
							cout<<'\n';
							cout<<"ku1, kv1: "<<'\n';
							for(uint k=0; k<ku1.size(); k++)
							{
								cout<<ku1[k]<<' ';
							}
							cout<<'\n';
							for(uint k=0; k<kv1.size(); k++)
							{
								cout<<kv1[k]<<' ';
							}
							cout<<'\n';
							getchar();
						}*/
						double coef = Tu[loc1][0] * Tv[loc2][0];
						if (i < npt_old)
						{
							cp[i].coortmp[0] += coef * cp[j].coor[0];
							cp[i].coortmp[1] += coef * cp[j].coor[1];
							cp[i].coortmp[2] += coef * cp[j].coor[2];
						}
						else
						{
							//if(i!=131 || j!=84)//tmp
							{
								cp[i].coor[0] += coef * cp[j].coor[0];
								cp[i].coor[1] += coef * cp[j].coor[1];
								cp[i].coor[2] += coef * cp[j].coor[2];
							}
						}
						cp[j].trun = 1;
						vector<int>::iterator it3 = find(cp[j].tbf.begin(), cp[j].tbf.end(), i);
						if (it3 == cp[j].tbf.end() /*&& (i!=131 || j!=84)*/)
						{
							cp[j].tbf.push_back(i);
							cp[j].tc.push_back(coef);
						}
					}
				}
			}
		}
	}
	//for(uint i=0; i<cp.size(); i++)
	//{
	//	if(cp[i].aff==1)
	//	{
	//		cp[i].coor[0]=cp[i].coortmp[0]; cp[i].coor[1]=cp[i].coortmp[1]; cp[i].coor[2]=cp[i].coortmp[2];
	//	}
	//}

	//for(uint i=0; i<cp[58].tbf.size(); i++)
	//{
	//	cout<<cp[58].tbf[i]<<" ";
	//}
	//cout<<'\n';
	//getchar();
}

void TruncatedTspline::UpdateControlPoints_3() //calculate control points, support of multiple levels
{
	for (uint i = 0; i < npt_old; i++)
	{
		cp[i].aff = 0;
		if (!equal(cp[i].kutmp, cp[i].kutmp + 5, cp[i].knotU) || !equal(cp[i].kvtmp, cp[i].kvtmp + 5, cp[i].knotV))
		{
			cp[i].aff = 1;
		}
	}
	vector<vector<double>> cmat(cp.size(), vector<double>(npt_old, 0.));
	vector<vector<double>> tmat(cp.size(), vector<double>(npt_old, 0.));
	for (uint i = 0; i < cp.size(); i++) //to be updated
	{
		//cp[i].truntmp=0;
		if (cp[i].trun == 0)
		{
			cp[i].coortmp[0] = 0.;
			cp[i].coortmp[1] = 0.;
			cp[i].coortmp[2] = 0.;
		}
		if (cp[i].aff == 1 || i >= npt_old)
		{
			for (uint j = 0; j < npt_old; j++) //old
			{
				if (CheckSubKnotVector(cp[i].kutmp, cp[i].kvtmp, cp[j].knotU, cp[j].knotV))
				{
					vector<double> ku1(10), kv1(10);
					vector<vector<double>> Tu, Tv;
					vector<double> ku(cp[j].knotU, cp[j].knotU + 5);
					vector<double> kv(cp[j].knotV, cp[j].knotV + 5);
					vector<double>::iterator it1, it2;
					it1 = set_union(cp[j].knotU, cp[j].knotU + 5, cp[i].kutmp, cp[i].kutmp + 5, ku1.begin());
					it2 = set_union(cp[j].knotV, cp[j].knotV + 5, cp[i].kvtmp, cp[i].kvtmp + 5, kv1.begin());
					ku1.resize(it1 - ku1.begin());
					kv1.resize(it2 - kv1.begin());
					TMatrix(ku, ku1, 3, Tu);
					TMatrix(kv, kv1, 3, Tv);
					it1 = search(ku1.begin(), ku1.end(), cp[i].kutmp, cp[i].kutmp + 5);
					it2 = search(kv1.begin(), kv1.end(), cp[i].kvtmp, cp[i].kvtmp + 5);
					if (it1 != ku1.end() && it2 != kv1.end())
					{
						int loc1 = it1 - ku1.begin();
						int loc2 = it2 - kv1.begin();
						double coef = Tu[loc1][0] * Tv[loc2][0];
						cmat[i][j] = coef;
						tmat[i][j] = coef;
						//if(i==49 && j==39)
						//{
						//	cout<<cmat[i][j]<<'\n';
						//	getchar();
						//}
						//if(coef!=0.)
						//{
						//	if(i<npt_old)
						//	{
						//		cp[i].coortmp[0]+=coef*cp[j].coor[0]; cp[i].coortmp[1]+=coef*cp[j].coor[1]; cp[i].coortmp[2]+=coef*cp[j].coor[2];
						//	}
						//	else
						//	{
						//		cp[i].coor[0]+=coef*cp[j].coor[0]; cp[i].coor[1]+=coef*cp[j].coor[1]; cp[i].coor[2]+=coef*cp[j].coor[2];
						//	}
						//	cp[j].trun=1;
						//	vector<int>::iterator it3=find(cp[j].tbf.begin(),cp[j].tbf.end(),i);
						//	if(it3==cp[j].tbf.end())
						//	{
						//		cp[j].tbf.push_back(i);
						//		cp[j].tc.push_back(coef);
						//	}
						//}
					}
				}
				/*int flag=CheckKnotInsertion(cp[i].kutmp,cp[i].kvtmp,cp[j].knotU,cp[j].knotV);
				double coef(0.);
				if(flag==0)
				{
					vector<double> ku1(10);
					vector<vector<double>> Tu;
					vector<double> ku(cp[j].knotU,cp[j].knotU+5);
					vector<double>::iterator it1;
					it1=set_union(cp[j].knotU,cp[j].knotU+5,cp[i].kutmp,cp[i].kutmp+5,ku1.begin());
					ku1.resize(it1-ku1.begin());
					TMatrix(ku,ku1,3,Tu);
					it1=search(ku1.begin(),ku1.end(),cp[i].kutmp,cp[i].kutmp+5);
					int loc1=it1-ku1.begin();
					coef=Tu[loc1][0];
				}
				else if(flag==1)
				{
					vector<double> kv1(10);
					vector<vector<double>> Tv;
					vector<double> kv(cp[j].knotV,cp[j].knotV+5);
					vector<double>::iterator it2;
					it2=set_union(cp[j].knotV,cp[j].knotV+5,cp[i].kvtmp,cp[i].kvtmp+5,kv1.begin());
					kv1.resize(it2-kv1.begin());
					TMatrix(kv,kv1,3,Tv);
					it2=search(kv1.begin(),kv1.end(),cp[i].kvtmp,cp[i].kvtmp+5);
					int loc2=it2-kv1.begin();
					coef=Tv[loc2][0];
				}
				else if(flag==2)
				{
					vector<double> ku1(10),kv1(10);
					vector<vector<double>> Tu,Tv;
					vector<double> ku(cp[j].knotU,cp[j].knotU+5);
					vector<double> kv(cp[j].knotV,cp[j].knotV+5);
					vector<double>::iterator it1,it2;
					it1=set_union(cp[j].knotU,cp[j].knotU+5,cp[i].kutmp,cp[i].kutmp+5,ku1.begin());
					it2=set_union(cp[j].knotV,cp[j].knotV+5,cp[i].kvtmp,cp[i].kvtmp+5,kv1.begin());
					ku1.resize(it1-ku1.begin());
					kv1.resize(it2-kv1.begin());
					TMatrix(ku,ku1,3,Tu);
					TMatrix(kv,kv1,3,Tv);
					it1=search(ku1.begin(),ku1.end(),cp[i].kutmp,cp[i].kutmp+5);
					it2=search(kv1.begin(),kv1.end(),cp[i].kvtmp,cp[i].kvtmp+5);
					int loc1=it1-ku1.begin();
					int loc2=it2-kv1.begin();
					coef=Tu[loc1][0]*Tv[loc2][0];
				}*/
				//if(flag!=-1 && coef!=0.)
				//{
				//	if(i<npt_old)
				//	{
				//		cp[i].coortmp[0]+=coef*cp[j].coor[0]; cp[i].coortmp[1]+=coef*cp[j].coor[1]; cp[i].coortmp[2]+=coef*cp[j].coor[2];
				//	}
				//	else
				//	{
				//		cp[i].coor[0]+=coef*cp[j].coor[0]; cp[i].coor[1]+=coef*cp[j].coor[1]; cp[i].coor[2]+=coef*cp[j].coor[2];
				//	}
				//	cp[j].trun=1;
				//	vector<int>::iterator it3=find(cp[j].tbf.begin(),cp[j].tbf.end(),i);
				//	if(it3==cp[j].tbf.end())
				//	{
				//		cp[j].tbf.push_back(i);
				//		cp[j].tc.push_back(coef);
				//	}
				//}
				//if(cp[i].kutmp[0]>=cp[j].knotU[0] && cp[i].kutmp[4]<=cp[j].knotU[4] &&
				//	cp[i].kvtmp[0]>=cp[j].knotV[0] && cp[i].kvtmp[4]<=cp[j].knotV[4])
				//{
				//if(!equal(cp[i].kutmp,cp[i].kutmp+5,cp[i].knotU) && equal(cp[i].kvtmp,cp[i].kvtmp+5,cp[i].knotV))
				//{
				//}
				//else if(equal(cp[i].kutmp,cp[i].kutmp+5,cp[i].knotU) && !equal(cp[i].kvtmp,cp[i].kvtmp+5,cp[i].knotV))
				//{
				//}
				//else if(!equal(cp[i].kutmp,cp[i].kutmp+5,cp[i].knotU) && !equal(cp[i].kvtmp,cp[i].kvtmp+5,cp[i].knotV))
				//{
				//}
				//}
				//if(CheckSubKnotVector(cp[i].kutmp,cp[i].kvtmp,cp[j].knotU,cp[j].knotV))
				//{
				//	{
				//		vector<double> ku1(10),kv1(10);
				//		vector<vector<double>> Tu,Tv;
				//		vector<double> ku(cp[j].knotU,cp[j].knotU+5);
				//		vector<double> kv(cp[j].knotV,cp[j].knotV+5);
				//		vector<double>::iterator it1,it2;
				//		it1=set_union(cp[j].knotU,cp[j].knotU+5,cp[i].kutmp,cp[i].kutmp+5,ku1.begin());
				//		it2=set_union(cp[j].knotV,cp[j].knotV+5,cp[i].kvtmp,cp[i].kvtmp+5,kv1.begin());
				//		ku1.resize(it1-ku1.begin());
				//		kv1.resize(it2-kv1.begin());
				//		TMatrix(ku,ku1,3,Tu);
				//		TMatrix(kv,kv1,3,Tv);
				//		it1=search(ku1.begin(),ku1.end(),cp[i].kutmp,cp[i].kutmp+5);
				//		it2=search(kv1.begin(),kv1.end(),cp[i].kvtmp,cp[i].kvtmp+5);
				//		int loc1=it1-ku1.begin();
				//		int loc2=it2-kv1.begin();
				//		/*if(it1==ku1.end() || it2==kv1.end())
				//		{
				//			cout<<"kutmp, kvtmp: "<<'\n';
				//			for(uint k=0; k<5; k++)
				//			{
				//				cout<<cp[i].kutmp[k]<<' ';
				//			}
				//			cout<<'\n';
				//			for(uint k=0; k<5; k++)
				//			{
				//				cout<<cp[i].kvtmp[k]<<' ';
				//			}
				//			cout<<'\n';
				//			cout<<"ku, kv: "<<'\n';
				//			for(uint k=0; k<5; k++)
				//			{
				//				cout<<cp[j].knotU[k]<<' ';
				//			}
				//			cout<<'\n';
				//			for(uint k=0; k<5; k++)
				//			{
				//				cout<<cp[j].knotV[k]<<' ';
				//			}
				//			cout<<'\n';
				//			cout<<"ku1, kv1: "<<'\n';
				//			for(uint k=0; k<ku1.size(); k++)
				//			{
				//				cout<<ku1[k]<<' ';
				//			}
				//			cout<<'\n';
				//			for(uint k=0; k<kv1.size(); k++)
				//			{
				//				cout<<kv1[k]<<' ';
				//			}
				//			cout<<'\n';
				//			getchar();
				//		}*/
				//		double coef=Tu[loc1][0]*Tv[loc2][0];
				//		if(i<npt_old)
				//		{
				//			cp[i].coortmp[0]+=coef*cp[j].coor[0]; cp[i].coortmp[1]+=coef*cp[j].coor[1]; cp[i].coortmp[2]+=coef*cp[j].coor[2];
				//		}
				//		else
				//		{
				//			//if(i!=131 || j!=84)//tmp
				//			{
				//			cp[i].coor[0]+=coef*cp[j].coor[0]; cp[i].coor[1]+=coef*cp[j].coor[1]; cp[i].coor[2]+=coef*cp[j].coor[2];
				//			}
				//		}
				//		cp[j].trun=1;
				//		vector<int>::iterator it3=find(cp[j].tbf.begin(),cp[j].tbf.end(),i);
				//		if(it3==cp[j].tbf.end() /*&& (i!=131 || j!=84)*/)
				//		{
				//			cp[j].tbf.push_back(i);
				//			cp[j].tc.push_back(coef);
				//		}
				//	}
				//}
			}
		}
	}

	for (uint i = 0; i < cp.size(); i++)
	{
		for (uint j = 0; j < npt_old; j++)
		{
			if (cmat[i][j] != 0.)
			{
				for (uint k = 0; k < cp[j].tbf.size(); k++)
				{
					tmat[i][j] -= cmat[i][cp[j].tbf[k]] * cp[j].tc[k];
				}
			}
		}
	}

	for (uint i = 0; i < cp.size(); i++)
	{
		for (uint j = 0; j < npt_old; j++)
		{
			if (tmat[i][j] != 0.) //control points
			{
				if (i < npt_old)
				{
					cp[i].coortmp[0] += tmat[i][j] * cp[j].coor[0];
					cp[i].coortmp[1] += tmat[i][j] * cp[j].coor[1];
					cp[i].coortmp[2] += tmat[i][j] * cp[j].coor[2];
				}
				else
				{
					cp[i].coor[0] += tmat[i][j] * cp[j].coor[0];
					cp[i].coor[1] += tmat[i][j] * cp[j].coor[1];
					cp[i].coor[2] += tmat[i][j] * cp[j].coor[2];
				}
				//cp[j].truntmp=1;
				/*cp[j].trun=1;
				vector<int>::iterator it3=find(cp[j].tbf.begin(),cp[j].tbf.end(),i);
				if(it3==cp[j].tbf.end())
				{
					cp[j].tbf.push_back(i);
					cp[j].tc.push_back(cmat[i][j]);
				}
				else
				{
					int loc=it3-cp[j].tbf.begin();
					cp[j].tc[loc]=cmat[i][j];
				}*/
			}
			if (cmat[i][j] != 0.) //truncation
			{
				cp[j].trun = 1;
				vector<int>::iterator it3 = find(cp[j].tbf.begin(), cp[j].tbf.end(), i);
				if (it3 == cp[j].tbf.end())
				{
					cp[j].tbf.push_back(i);
					cp[j].tc.push_back(cmat[i][j]);
				}
				else
				{
					int loc = it3 - cp[j].tbf.begin();
					cp[j].tc[loc] = cmat[i][j];
				}
			}
		}
	}
}

void TruncatedTspline::UpdateControlPoints_4() //calculate control points, support of multiple levels
{
	for (uint i = 0; i < npt_old; i++)
	{
		cp[i].aff = 0;
		if (!equal(cp[i].kutmp, cp[i].kutmp + 5, cp[i].knotU) || !equal(cp[i].kvtmp, cp[i].kvtmp + 5, cp[i].knotV))
		{
			cp[i].aff = 1;
		}
	}
	vector<vector<double>> cmat(cp.size(), vector<double>(npt_old, 0.));
	vector<vector<double>> wmat(cp.size(), vector<double>(npt_old, 0.));
	//vector<vector<double>> tmat(cp.size(),vector<double>(npt_old,0.));
	for (uint i = 0; i < cp.size(); i++) //to be updated
	{
		cp[i].coortmp[0] = 0.;
		cp[i].coortmp[1] = 0.;
		cp[i].coortmp[2] = 0.;
		cp[i].wtmp = 0.;
		if (cp[i].aff == 1 || i >= npt_old)
		{
			for (uint j = 0; j < npt_old; j++) //old
			{
				if (CheckSubKnotVector(cp[i].kutmp, cp[i].kvtmp, cp[j].knotU, cp[j].knotV))
				{
					vector<double> ku1(10), kv1(10);
					vector<vector<double>> Tu, Tv;
					vector<double> ku(cp[j].knotU, cp[j].knotU + 5);
					vector<double> kv(cp[j].knotV, cp[j].knotV + 5);
					vector<double>::iterator it1, it2;
					it1 = set_union(cp[j].knotU, cp[j].knotU + 5, cp[i].kutmp, cp[i].kutmp + 5, ku1.begin());
					it2 = set_union(cp[j].knotV, cp[j].knotV + 5, cp[i].kvtmp, cp[i].kvtmp + 5, kv1.begin());
					ku1.resize(it1 - ku1.begin());
					kv1.resize(it2 - kv1.begin());
					TMatrix(ku, ku1, 3, Tu);
					TMatrix(kv, kv1, 3, Tv);
					it1 = search(ku1.begin(), ku1.end(), cp[i].kutmp, cp[i].kutmp + 5);
					it2 = search(kv1.begin(), kv1.end(), cp[i].kvtmp, cp[i].kvtmp + 5);
					if (it1 != ku1.end() && it2 != kv1.end())
					{
						int loc1 = it1 - ku1.begin();
						int loc2 = it2 - kv1.begin();
						double coef = Tu[loc1][0] * Tv[loc2][0];
						cmat[i][j] = coef;
						wmat[i][j] = coef;
						//tmat[i][j]=coef;
					}
				}
			}
		}
	}

	for (uint i = 0; i < cp.size(); i++)
	{
		for (uint j = 0; j < npt_old; j++)
		{
			if (cmat[i][j] != 0.)
			{
				for (uint k = 0; k < cp[j].tbf.size(); k++)
				{
					wmat[i][j] -= cmat[i][cp[j].tbf[k]] * cp[j].tc[k];
					//tmat[i][j]+=cmat[i][cp[j].tbf[k]]*cp[j].tc[k];
				}
			}
		}
	}

	for (uint i = npt_old; i < cp.size(); i++)
	{
		cp[i].coor[0] = 0.;
		cp[i].coor[1] = 0.;
		cp[i].coor[0] = 0.;
		cp[i].w = 0.;
	}
	for (uint i = 0; i < cp.size(); i++)
	{
		for (uint j = 0; j < npt_old; j++)
		{
			if (wmat[i][j] != 0.) //control points
			{
				if (i < npt_old)
				{
					cp[i].coortmp[0] += wmat[i][j] * cp[j].coor[0] * cp[j].w;
					cp[i].coortmp[1] += wmat[i][j] * cp[j].coor[1] * cp[j].w;
					cp[i].coortmp[2] += wmat[i][j] * cp[j].coor[2] * cp[j].w;
					cp[i].wtmp += wmat[i][j] * cp[j].w;
				}
				else
				{
					cp[i].coor[0] += wmat[i][j] * cp[j].coor[0] * cp[j].w;
					cp[i].coor[1] += wmat[i][j] * cp[j].coor[1] * cp[j].w;
					cp[i].coor[2] += wmat[i][j] * cp[j].coor[2] * cp[j].w;
					cp[i].w += wmat[i][j] * cp[j].w;
				}
			}
			if (cmat[i][j] != 0.) //truncation
			{
				cp[j].trun = 1;
				vector<int>::iterator it3 = find(cp[j].tbf.begin(), cp[j].tbf.end(), i);
				if (it3 == cp[j].tbf.end())
				{
					cp[j].tbf.push_back(i);
					cp[j].tc.push_back(cmat[i][j]);
				}
				else
				{
					int loc = it3 - cp[j].tbf.begin();
					cp[j].tc[loc] = cmat[i][j];
				}
			}
		}
		//rational
		if (i < npt_old)
		{
			cp[i].coortmp[0] /= cp[i].wtmp;
			cp[i].coortmp[1] /= cp[i].wtmp;
			cp[i].coortmp[2] /= cp[i].wtmp;
		}
		else
		{
			cp[i].coor[0] /= cp[i].w;
			cp[i].coor[1] /= cp[i].w;
			cp[i].coor[2] /= cp[i].w;
		}
	}
}

void TruncatedTspline::CheckTruncation()
{
	double tol(1.e-8);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].trun == 1)
		{
			double Ni = EvaluateTrunBF(i);
			if (Ni < 0.)
				Ni = -Ni;
			if (Ni < tol)
			{
				cp[i].trun = 0;
				cp[i].tbf.clear();
				cp[i].tc.clear();
				cp[i].coor[0] = cp[i].coortmp[0];
				cp[i].coor[1] = cp[i].coortmp[1];
				cp[i].coor[2] = cp[i].coortmp[2];
			}
			if (i == 47 || i == 73 || i == 39 || 81)
				cp[i].trun = 1; //tmp for square surface
		}
	}

	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].trun == 1)
		{
			vector<int> tbf_tmp(cp[i].tbf);
			vector<double> tc_tmp(cp[i].tc);
			cp[i].tbf.clear();
			cp[i].tc.clear();
			for (uint j = 0; j < tbf_tmp.size(); j++)
			{
				if (tbf_tmp[j] != i && cp[tbf_tmp[j]].trun == 0)
				{
					cp[i].tbf.push_back(tbf_tmp[j]);
					cp[i].tc.push_back(tc_tmp[j]);
				}
			}
		}
	}
}

void TruncatedTspline::CheckTruncation_1()
{
	double tol(1.e-8);
	//for(uint i=0; i<cp.size(); i++)
	//{
	//	if(cp[i].truntmp==1 && cp[i].trun==1)
	//	{
	//		for(uint j=0; j<cp[i].tbf.size(); j++)
	//		{
	//			if(cp[cp[i].tbf[j]].aff==0)
	//			{
	//				vector<int>::iterator it=find(cp[i].tbftmp.begin(),cp[i].tbftmp.end(),cp[i].tbf[j]);
	//				if(it==cp[i].tbftmp.end())
	//				{
	//					cp[i].tbftmp.push_back(cp[i].tbf[j]);
	//					cp[i].tctmp.push_back(cp[i].tc[j]);
	//				}
	//			}
	//		}
	//		vector<int>::iterator it1=find(cp[i].tbftmp.begin(),cp[i].tbftmp.end(),i);
	//		if(it1==cp[i].tbftmp.end())
	//		{
	//			cp[i].tbftmp.push_back(i);
	//			cp[i].tctmp.push_back(1.);
	//		}
	//	}
	//}

	//cout<<"before!\n";
	//for(uint i=0; i<cp[50].tbf.size(); i++)
	//{
	//	cout<<cp[50].tbf[i]<<" ";
	//}
	//cout<<"\n";
	//for(uint i=0; i<cp[72].tc.size(); i++)
	//{
	//	cout<<cp[72].tc[i]<<" ";
	//}
	//cout<<"\n";

	for (uint i = 0; i < cp.size(); i++) //I_1^1 and I_2^1
	{
		if (cp[i].trun == 1)
		{
			double Ni = EvaluateTrunBF(i);
			if (Ni < 0.)
				Ni = -Ni;
			if (Ni < tol)
			{
				cp[i].trun = 0;
				cp[i].tbf.clear();
				cp[i].tc.clear();
				if (cp[i].aff == 1)
				{
					cp[i].coor[0] = cp[i].coortmp[0];
					cp[i].coor[1] = cp[i].coortmp[1];
					cp[i].coor[2] = cp[i].coortmp[2];
				}
			}
		}
	}

	for (uint i = 0; i < cp.size(); i++) //T_1^1
	{
		if (cp[i].trun == 1)
		{
			//double Ni=EvaluateTrunBF1(i);
			//if(Ni<0.) Ni=-Ni;
			//if(Ni<tol)
			//{
			//	cp[i].trun=0;
			//	cp[i].tbf.clear();
			//	cp[i].tc.clear();
			//	if(cp[i].aff==1)
			//	{
			//		cp[i].coor[0]=cp[i].coortmp[0]; cp[i].coor[1]=cp[i].coortmp[1]; cp[i].coor[2]=cp[i].coortmp[2];
			//	}
			//}
			//else
			//{
			vector<int> tbf_tmp(cp[i].tbf);
			vector<double> tc_tmp(cp[i].tc);
			cp[i].tbf.clear();
			cp[i].tc.clear();
			for (uint j = 0; j < tbf_tmp.size(); j++)
			{
				//if(tbf_tmp[j]!=i && tc_tmp[j]!=0.)
				if (tbf_tmp[j] != i && tc_tmp[j] != 0. && cp[tbf_tmp[j]].trun == 0)
				{
					cp[i].tbf.push_back(tbf_tmp[j]);
					cp[i].tc.push_back(tc_tmp[j]);
				}
			}
			//}
		}
	}

	//for(uint i=0; i<cp.size(); i++)
	//{
	//	if(cp[i].trun==1)
	//	{
	//		double Ni=EvaluateTrunBF1(i);
	//		if(Ni<0.) Ni=-Ni;
	//		if(Ni<tol)
	//		{
	//			cp[i].trun=0;
	//			cp[i].tbf.clear();
	//			cp[i].tc.clear();
	//			if(cp[i].aff==1)
	//			{
	//				cp[i].coor[0]=cp[i].coortmp[0]; cp[i].coor[1]=cp[i].coortmp[1]; cp[i].coor[2]=cp[i].coortmp[2];
	//			}
	//		}
	//	}
	//}

	//cout<<"after!\n";
	//for(uint i=0; i<cp[50].tbf.size(); i++)
	//{
	//	cout<<cp[50].tbf[i]<<" ";
	//}
	//cout<<"\n";
}

double TruncatedTspline::EvaluateTrunBF(int loc)
{
	int nsp(5);
	double sp[5][2] = {{cp[loc].knotU[2], cp[loc].knotV[2]}, {cp[loc].knotU[2], cp[loc].knotV[1]}, {cp[loc].knotU[3], cp[loc].knotV[2]}, {cp[loc].knotU[2], cp[loc].knotV[3]}, {cp[loc].knotU[1], cp[loc].knotV[2]}};
	double sum(0.);
	for (int i = 0; i < nsp; i++)
	{
		vector<double> ku(cp[loc].knotU, cp[loc].knotU + 5);
		vector<double> kv(cp[loc].knotV, cp[loc].knotV + 5);
		vector<double> uval, vval;
		BSplineBasis bu, bv;
		bu.Set(3, ku);
		bv.Set(3, kv);
		bu.BasisFunction(0, sp[i][0], 0, uval);
		bv.BasisFunction(0, sp[i][1], 0, vval);
		double Ni = uval[0] * vval[0];
		for (uint j = 0; j < cp[loc].tbf.size(); j++)
		{
			int loc1 = cp[loc].tbf[j];
			ku.assign(cp[loc1].kutmp, cp[loc1].kutmp + 5);
			kv.assign(cp[loc1].kvtmp, cp[loc1].kvtmp + 5);
			bu.Set(3, ku);
			bv.Set(3, kv);
			bu.BasisFunction(0, sp[i][0], 0, uval);
			bv.BasisFunction(0, sp[i][1], 0, vval);
			Ni -= cp[loc].tc[j] * uval[0] * vval[0];
		}
		sum += Ni;
	}
	return (sum / nsp);
}

double TruncatedTspline::EvaluateTrunBF1(int loc)
{
	int nsp(5);
	double sp[5][2] = {{cp[loc].knotU[2], cp[loc].knotV[2]}, {cp[loc].knotU[2], cp[loc].knotV[1]}, {cp[loc].knotU[3], cp[loc].knotV[2]}, {cp[loc].knotU[2], cp[loc].knotV[3]}, {cp[loc].knotU[1], cp[loc].knotV[2]}};
	double sum(0.);
	for (int i = 0; i < nsp; i++)
	{
		vector<double> ku(cp[loc].knotU, cp[loc].knotU + 5);
		vector<double> kv(cp[loc].knotV, cp[loc].knotV + 5);
		vector<double> uval, vval;
		BSplineBasis bu, bv;
		bu.Set(3, ku);
		bv.Set(3, kv);
		bu.BasisFunction(0, sp[i][0], 0, uval);
		bv.BasisFunction(0, sp[i][1], 0, vval);
		double Ni = uval[0] * vval[0];
		for (uint j = 0; j < cp[loc].tbf.size(); j++)
		{
			int loc1 = cp[loc].tbf[j];
			ku.assign(cp[loc1].kutmp, cp[loc1].kutmp + 5);
			kv.assign(cp[loc1].kvtmp, cp[loc1].kvtmp + 5);
			bu.Set(3, ku);
			bv.Set(3, kv);
			bu.BasisFunction(0, sp[i][0], 0, uval);
			bv.BasisFunction(0, sp[i][1], 0, vval);
			double Nitmp = uval[0] * vval[0];
			for (uint k = 0; k < cp[loc1].tbf.size(); k++)
			{
				int loc2 = cp[loc1].tbf[k];
				if (loc1 != loc2)
				{
					ku.assign(cp[loc2].kutmp, cp[loc2].kutmp + 5);
					kv.assign(cp[loc2].kvtmp, cp[loc2].kvtmp + 5);
					bu.Set(3, ku);
					bv.Set(3, kv);
					bu.BasisFunction(0, sp[i][0], 0, uval);
					bv.BasisFunction(0, sp[i][1], 0, vval);
					Nitmp -= cp[loc1].tc[k] * uval[0] * vval[0];
				}
			}
			Ni -= cp[loc].tc[j] * Nitmp;
		}
		sum += Ni;
	}
	return (sum / nsp);
}

void TruncatedTspline::UpdateKnotVectors()
{
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].trun == 0)
		{
			for (int j = 0; j < 5; j++)
			{
				cp[i].knotU[j] = cp[i].kutmp[j];
				cp[i].knotV[j] = cp[i].kvtmp[j];
			}
			double minu(10.), minv(10.);
			for (int j = 0; j < 4; j++)
			{
				cp[i].kitvU[j] = cp[i].knotU[j + 1] - cp[i].knotU[j];
				cp[i].kitvV[j] = cp[i].knotV[j + 1] - cp[i].knotV[j];
				if (minu > cp[i].kitvU[j] && cp[i].kitvU[j] != 0)
					minu = cp[i].kitvU[j];
				if (minv > cp[i].kitvV[j] && cp[i].kitvV[j] != 0)
					minv = cp[i].kitvV[j];
			}
			cp[i].min_itv[0] = minu;
			cp[i].min_itv[1] = minv;
		}
	}
}

void TruncatedTspline::Refine_Addition(const vector<int> &rid)
{
	npt_old = cp.size();
	nel_old = tmesh.size();
	vector<int> ridsb, ridtp;
	StrongBalanceCheck(rid, ridsb);
	StrongBalanceRefine(ridsb);
	OneTjunctionCheck(ridtp);
	OneTjunctionRefine(ridtp);
	UpdateTopology();
	FindLocalKnotVectors();
	GeometryRefine();
}

void TruncatedTspline::Refine_Target(const vector<int> &rid)
{
	npt_old = cp.size();
	nel_old = tmesh.size();
	TargetRefine(rid);
	UpdateTopology();
	int niter(0);
	while (1)
	{
		niter++;
		if (niter > 20)
			break;
		vector<int> ridtp1, ridsb1;
		OneTjunctionCheck(ridtp1);
		if (ridtp1.size() == 0)
			break;
		StrongBalanceCheck(ridtp1, ridsb1);
		StrongBalanceRefine(ridsb1);
		OneTjunctionRefine(ridtp1);
		UpdateTopology();
	}

	//extend edges
	niter = 0;
	while (1)
	{
		niter++;
		if (niter > 20)
			break;
		vector<int> ridtjx, ridsb2;
		TjuncExtentCheck_2(ridtjx);
		if (ridtjx.size() == 0)
			break;
		StrongBalanceCheck(ridtjx, ridsb2);
		StrongBalanceRefine(ridsb2);
		TjuncExtentRefine(ridtjx);
		UpdateTopology();
	}

	FindLocalKnotVectors();

	GeometryRefine();
}

void TruncatedTspline::GeometryRefine()
{
	UpdateControlPoints_4();
	CheckTruncation();
	UpdateKnotVectors();
	FindIENglb();

	//Truncation();

	//for(uint i=0; i<cp.size(); i++)
	//{
	//	if(cp[i].trun==1)
	//	{
	//		cout<<"pid: "<<i<<"\ntbf: ";
	//		for(uint j=0; j<cp[i].tbf.size(); j++)
	//		{
	//			cout<<cp[i].tbf[j]<<' ';
	//			//cout<<"tbf: "<<cp[i].tbf[j]<<' '<<cp[i].tc[j]<<'\n';
	//			//for(int k=0; k<5; k++)
	//			//{
	//			//	cout<<cp[cp[i].tbf[j]].knotU[k]<<' ';
	//			//}
	//			//cout<<'\n';
	//			//for(int k=0; k<5; k++)
	//			//{
	//			//	cout<<cp[cp[i].tbf[j]].knotV[k]<<' ';
	//			//}
	//			//cout<<'\n';
	//		}
	//		cout<<"\ntc: ";
	//		for(uint j=0; j<cp[i].tbf.size(); j++)
	//		{
	//			cout<<cp[i].tc[j]<<' ';
	//		}
	//		cout<<'\n';
	//		//getchar();
	//	}
	//}
	//getchar();
}

void TruncatedTspline::TopologyRefineTest()
{
	int npt_new(15);
	for (int i = 0; i < npt_new; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int eids[4] = {44, 45, 54, 55};
	int ref_ids[4] = {2, 1, 1, 1};
	int pnewids[4][5] = {{134, 131, 121, -1, 135}, {130, 133, 126, 131, 132}, {121, 124, 125, 122, 123}, {126, 128, 129, 124, 127}};
	for (int i = 0; i < 4; i++)
	{
		tmesh[eids[i]].ref = ref_ids[i];
		tmesh[eids[i]].newpID.resize(5);
		for (int j = 0; j < 5; j++)
		{
			tmesh[eids[i]].newpID[j] = pnewids[i][j];
		}
	}

	//set knot vectors
	int affid[20] = {38, 39, 47, 48, 49, 50, 51, 58, 59, 60, 61, 62, 69, 70, 71, 72, 73, 81, 82, 83};
	for (int i = 0; i < 20; i++)
	{
		cp[affid[i]].aff = 1;
		for (int j = 0; j < 5; j++)
		{
			cp[affid[i]].kutmp[j] = cp[affid[i]].knotU[j];
			cp[affid[i]].kvtmp[j] = cp[affid[i]].knotV[j];
		}
	}
	double ktsU[9] = {cp[59].knotU[0], cp[59].knotU[1], cp[59].knotU[2], (cp[59].knotU[2] + cp[59].knotU[3]) / 2., cp[59].knotU[3], (cp[59].knotU[3] + cp[59].knotU[4]) / 2., cp[59].knotU[4], cp[61].knotU[3], cp[61].knotU[4]};
	double ktsV[9] = {cp[49].knotV[0], cp[49].knotV[1], cp[49].knotV[2], (cp[49].knotV[2] + cp[49].knotV[3]) / 2., cp[49].knotV[3], (cp[49].knotV[3] + cp[49].knotV[4]) / 2., cp[49].knotV[4], cp[71].knotV[3], cp[71].knotV[4]};
	int setKU[5][5] = {{-1, 134, -1, 130, -1}, {-1, 135, 131, 132, 133}, {-1, 121, -1, 126, -1}, {122, 123, 124, 127, 128}, {-1, 125, -1, 129, -1}};
	int setKV[5][5] = {{-1, -1, -1, 122, -1}, {134, 135, 121, 123, 125}, {-1, 131, -1, 124, -1}, {130, 132, 126, 127, 129}, {-1, 133, -1, 128, -1}};
	int setKU1[5][5] = {{48, -1, 49, -1, 50}, {-1, -1, -1, -1, -1}, {59, -1, 60, -1, 61}, {-1, -1, -1, -1, -1}, {70, -1, 71, -1, 72}};
	int setKV1[5][5] = {{-1, -1, 59, -1, 70}, {-1, -1, -1, -1, -1}, {49, -1, 60, -1, 71}, {-1, -1, -1, -1, -1}, {50, -1, 61, -1, 72}};
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			if (setKU[i][j] != -1)
			{
				for (int k = 0; k < 5; k++)
				{
					cp[setKU[i][j]].knotU[k] = ktsU[j + k];
				}
			}
		}
	}
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			if (setKV[i][j] != -1)
			{
				for (int k = 0; k < 5; k++)
				{
					cp[setKV[i][j]].knotV[k] = ktsV[j + k];
				}
			}
		}
	}
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			if (setKU1[i][j] != -1)
			{
				for (int k = 0; k < 5; k++)
				{
					cp[setKU1[i][j]].kutmp[k] = ktsU[j + k];
				}
			}
		}
	}
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			if (setKV1[i][j] != -1)
			{
				for (int k = 0; k < 5; k++)
				{
					cp[setKV1[i][j]].kvtmp[k] = ktsV[j + k];
				}
			}
		}
	}

	cp[38].kvtmp[4] = cp[131].knotV[2];
	cp[39].kvtmp[4] = cp[133].knotV[2];
	cp[47].kutmp[4] = cp[121].knotU[2];
	cp[48].kvtmp[4] = cp[122].knotV[2];
	cp[51].kutmp[0] = cp[130].knotU[2];
	cp[58].kutmp[4] = cp[121].knotU[2];
	cp[62].kutmp[0] = cp[126].knotU[2];
	cp[69].kutmp[4] = cp[125].knotU[2];
	cp[73].kutmp[0] = cp[129].knotU[2];
	cp[81].kvtmp[0] = cp[122].knotV[2];
	cp[82].kvtmp[0] = cp[124].knotV[2];
	cp[83].kvtmp[0] = cp[128].knotV[2];
	cp[59].kvtmp[0] = cp[37].knotV[2];
	cp[59].kvtmp[1] = cp[48].knotV[2];
	cp[122].knotV[0] = cp[48].knotV[2];

	//new elements
	int ne_new(15);
	int enew[15][4] = {{59, 121, 123, 122}, {121, 60, 124, 123}, {123, 124, 71, 125}, {122, 123, 125, 70}, {60, 126, 127, 124}, {126, 61, 128, 127}, {127, 128, 72, 129}, {124, 127, 129, 71}, {49, 130, 132, 131}, {130, 50, 133, 132}, {132, 133, 61, 126}, {131, 132, 126, 60}, {48, 134, 121, 59}, {134, 49, 131, 135}, {135, 131, 60, 121}};
	//int idtmp[15]={121,122,123,124,125,126,127,128,129,130,131,132,133,134,135};
	for (int i = 0; i < ne_new; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		tmesh.push_back(etmp);
	}
}

void TruncatedTspline::TopologyRefineTest1()
{
	int npt_new(5);
	for (int i = 0; i < npt_new; i++)
	{
		Vertex ptmp;
		cp.push_back(ptmp);
	}
	int eids[1] = {104};
	int ref_ids[1] = {1};
	int pnewids[1][5] = {{136, 139, 140, 137, 138}};
	for (int i = 0; i < 1; i++)
	{
		tmesh[eids[i]].ref = ref_ids[i];
		tmesh[eids[i]].newpID.resize(5);
		for (int j = 0; j < 5; j++)
		{
			tmesh[eids[i]].newpID[j] = pnewids[i][j];
		}
	}

	//set knot vectors
	int affid[12] = {131, 132, 121, 61, 123, 128, 71, 129, 60, 126, 124, 127};
	for (int i = 0; i < 12; i++)
	{
		cp[affid[i]].aff = 1;
		for (int j = 0; j < 5; j++)
		{
			cp[affid[i]].kutmp[j] = cp[affid[i]].knotU[j];
			cp[affid[i]].kvtmp[j] = cp[affid[i]].knotV[j];
		}
	}
	double ktsU[7] = {cp[60].knotU[0], cp[60].knotU[1], cp[60].knotU[2], (cp[60].knotU[2] + cp[60].knotU[3]) / 2., cp[60].knotU[3], cp[60].knotU[4], cp[126].knotU[4]};
	double ktsV[7] = {cp[60].knotV[0], cp[60].knotV[1], cp[60].knotV[2], (cp[60].knotV[2] + cp[60].knotV[3]) / 2., cp[60].knotV[3], cp[60].knotV[4], cp[124].knotV[4]};
	int setKU[3][3] = {{-1, 136, -1}, {137, 138, 139}, {-1, 140, -1}};
	int setKV[3][3] = {{-1, 137, -1}, {136, 138, 140}, {-1, 139, -1}};
	int setKU1[3][3] = {{60, -1, 126}, {-1, -1, -1}, {124, -1, 127}};
	int setKV1[3][3] = {{60, -1, 124}, {-1, -1, -1}, {126, -1, 127}};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (setKU[i][j] != -1)
			{
				for (int k = 0; k < 5; k++)
				{
					cp[setKU[i][j]].knotU[k] = ktsU[j + k];
				}
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (setKV[i][j] != -1)
			{
				for (int k = 0; k < 5; k++)
				{
					cp[setKV[i][j]].knotV[k] = ktsV[j + k];
				}
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (setKU1[i][j] != -1)
			{
				for (int k = 0; k < 5; k++)
				{
					cp[setKU1[i][j]].kutmp[k] = ktsU[j + k];
				}
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (setKV1[i][j] != -1)
			{
				for (int k = 0; k < 5; k++)
				{
					cp[setKV1[i][j]].kvtmp[k] = ktsV[j + k];
				}
			}
		}
	}

	cp[131].kvtmp[4] = cp[137].knotV[2];
	cp[132].kvtmp[4] = cp[139].knotV[2];
	cp[121].kutmp[4] = cp[136].knotU[2];
	cp[61].kutmp[0] = cp[136].knotU[2];
	cp[123].kutmp[4] = cp[140].knotU[2];
	cp[128].kutmp[0] = cp[140].knotU[2];
	cp[71].kvtmp[0] = cp[137].knotV[2];
	cp[129].kvtmp[0] = cp[139].knotV[2];

	//new elements
	int ne_new(4);
	int enew[4][4] = {{60, 136, 138, 137}, {136, 126, 139, 138}, {137, 138, 140, 124}, {138, 139, 127, 140}};
	for (int i = 0; i < ne_new; i++)
	{
		Element etmp;
		etmp.act = 1;
		etmp.cnct[0] = enew[i][0];
		etmp.cnct[1] = enew[i][1];
		etmp.cnct[2] = enew[i][2];
		etmp.cnct[3] = enew[i][3];
		tmesh.push_back(etmp);
	}
}

void TruncatedTspline::TopologyRefineTest2()
{
	InitialTopology();
	//VisualizeTMesh("test4/CM_1");
	const int iter(3);
	for (int i = 0; i < iter; i++)
	{
		//cout<<i<<'\n';
		vector<int> rid;
		IdentifyTest1(rid);
		TopologyRefine(rid);
		//VisualizeTMesh("test4/CM_5_"+to_string(long long(i)));
	}
}

void TruncatedTspline::TopologyRefineTest3()
{
	InitialTopology();

	//int rid_tmp[2]={44,55};
	//int rid_tmp1[1]={104};
	//vector<int> rid1(rid_tmp,rid_tmp+2);
	//vector<int> rid2(rid_tmp1,rid_tmp1+1);
	////npt_old=cp.size();
	////TopologyRefine_1(rid1);
	//Refine_Addition(rid1);
	//Refine_Target(rid1);
	////GeometryRefine();
	////VisualizeTMesh("test6/surf_14_CM");
	////VisualizeVTK("test6/surf_8_surf");
	////npt_old=cp.size();
	////TopologyRefine_1(rid2);
	//Refine_Addition(rid2);
	//Refine_Target(rid2);
	////GeometryRefine();
	//VisualizeTMesh("test6/surf_15_CM");
	//VisualizeVTK("test6/surf_15_surf");

	//const int iter(7);
	//for(uint i=0; i<iter; i++)
	//{
	//	vector<int> rid1=cp[60].face;
	//	npt_old=cp.size();
	//	TopologyRefine_1(rid1);
	//	GeometryRefine();
	//}
	//VisualizeTMesh("test6/surf_13_CM");
	//VisualizeVTK("test6/surf_13_surf");

	//npt_old=cp.size();
	//vector<int> rid;
	//IdentifyTest2(rid);
	//TopologyRefine_1(rid);
	//GeometryRefine();
	////VisualizeTMesh("test6/surf_10_CM");
	////VisualizeVTK("test6/surf_10_surf");
	//npt_old=cp.size();
	//IdentifyTest(rid);
	//TopologyRefine_1(rid);
	//GeometryRefine();
	//VisualizeTMesh("test6/surf_11_CM");
	//VisualizeVTK("test6/surf_11_surf");

	//vector<int> rid;
	//IdentifyTest(rid);
	//Refine_Addition(rid);
	//Refine_Target(rid);
	//IdentifyTest(rid);
	//Refine_Addition(rid);
	//VisualizeTMesh("test6/surf_CM_16_1");
	//VisualizeVTK("test6/surf_surf_16_1");

	//for(uint i=0; i<cp[82].tbf.size(); i++)
	//{
	//	cout<<cp[82].tbf[i]<<" ";
	//}
	//cout<<"\n";
	//for(uint i=0; i<cp[72].tc.size(); i++)
	//{
	//	cout<<cp[72].tc[i]<<" ";
	//}
	//cout<<"\n";

	//const int iter(2);
	//for(int i=0; i<iter; i++)
	//{
	//	//npt_old=cp.size();
	//	cout<<"Iteration setp: "<<i<<'\n';
	//	vector<int> rid;
	//	IdentifyTest(rid);
	//	Refine_Addition(rid);
	//	VisualizeTMesh("test6/surf_CM_"+to_string(long long(16+i))+"_0");
	//	Refine_Target(rid);
	//	//TopologyRefine_1(rid);
	//	//GeometryRefine();
	//	//VisualizeTMesh("test6/surf_CM_"+to_string(long long(16+i)));
	//}
	//VisualizeTMesh("test6/surf_16_CM");
	//VisualizeVTK("test6/surf_14_surf");

	CollectActives();
	//cout<<"Bezier extracting...\n";
	//vector<BezierElement> bzmesh;
	//BezierExtract(bzmesh);
	//cout<<"Visualizing...\n";
	//BezierVTK("test6/patchtest_1_bezier",bzmesh);

	//npt_old=cp.size();
	////int refid[4]={35,45,54,55};
	////vector<int> rid(refid,refid+4);
	////TopologyRefine_1(rid);

	//ElementRefine_Square_4(45);
	////ElementRefine_Square_4(53);
	//ElementRefine_Square_4(54);
	//ElementRefine_Square_4(55);
	////ElementRefine_Square_4(35);
	//ElementRefine_Square_2(44);
	////ElementRefine_Square_2(34);
	////ElementRefine_Square_2(33);
	//UpdateTopology();
	//FindLocalKnotVectors();

	//GeometryRefine();

	//VisualizeTMesh("test4/CM_1");
	//const int iter(3);
	//for(int i=0; i<iter; i++)
	//{
	//	//cout<<i<<'\n';
	//	vector<int> rid;
	//	IdentifyTest1(rid);
	//	TopologyRefine(rid);
	//	//VisualizeTMesh("test4/CM_5_"+to_string(long long(i)));
	//}
}

void TruncatedTspline::TopologyRefineTest4()
{
	InitialTopology();

	npt_old = cp.size();
	ElementRefine_Square_4(45);
	ElementRefine_Square_4(54);
	ElementRefine_Square_4(55);
	ElementRefine_Square_3(44);
	UpdateTopology();
	FindLocalKnotVectors();
	GeometryRefine();
	double tmp = cp[84].tc[0];
	npt_old = cp.size();
	ElementRefine_Square_4(108);
	UpdateTopology();
	FindLocalKnotVectors();
	GeometryRefine();
	cp[84].tbf.push_back(72);
	cp[84].tc.push_back(tmp);

	//vector<int>::iterator it=find(cp[84].tbf.begin(),cp[84].tbf.end(),131);
	//int loc=it-cp[84].tbf.begin();
	//cp[84].tc[loc]=0.;

	//VisualizeTMesh("test6/patchtest_2_CM");
	//VisualizeVTK("test6/patchtest_2_surf");

	CollectActives();
}

void TruncatedTspline::UpdateIEN()
{
	uint i, j;
	for (i = 0; i < tmesh.size(); i++)
	{
		tmesh[i].IEN.clear();
		if (tmesh[i].act == 1)
		{
			for (j = 0; j < cp.size(); j++)
			{
				if (cp[j].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[j].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[j].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[j].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(j);
				}
			}
		}
	}
}

void TruncatedTspline::IdentifyTest(vector<int> &rid)
{
	rid.clear();
	double tol(1.e-6);
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			double tmp = (cp[tmesh[i].cnct[0]].knotU[2] + cp[tmesh[i].cnct[1]].knotU[2]) / 2. - (cp[tmesh[i].cnct[0]].knotV[2] + cp[tmesh[i].cnct[3]].knotV[2]) / 2.;
			if (tmp > -1. * tol && tmp < tol)
				rid.push_back(i);
		}
	}
}

void TruncatedTspline::IdentifyTest2(vector<int> &rid)
{
	rid.clear();
	vector<int> ridtmp, ridtmp1;
	double tol(1.e-6);
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].type == 0)
		{
			double tmp = (cp[tmesh[i].cnct[0]].knotU[2] + cp[tmesh[i].cnct[1]].knotU[2]) / 2. - (cp[tmesh[i].cnct[0]].knotV[2] + cp[tmesh[i].cnct[3]].knotV[2]) / 2.;
			if (tmp > -1. * tol && tmp < tol)
				ridtmp.push_back(i);
		}
	}
	for (uint i = 0; i < ridtmp.size(); i++)
	{
		ridtmp1.push_back(ridtmp[i]);
		for (int j = 0; j < 4; j++)
		{
			int fcid = tmedge[tmesh[ridtmp[i]].edge[j]].face[0];
			if (fcid == ridtmp[i])
				fcid = tmedge[tmesh[ridtmp[i]].edge[j]].face[1];
			vector<int>::iterator it = find(ridtmp1.begin(), ridtmp1.end(), fcid);
			if (it == ridtmp1.end() && tmesh[fcid].type == 0)
				ridtmp1.push_back(fcid);
		}
	}
	for (uint i = 0; i < ridtmp1.size(); i++)
	{
		rid.push_back(ridtmp1[i]);
		for (int j = 0; j < 4; j++)
		{
			int fcid = tmedge[tmesh[ridtmp1[i]].edge[j]].face[0];
			if (fcid == ridtmp1[i])
				fcid = tmedge[tmesh[ridtmp1[i]].edge[j]].face[1];
			vector<int>::iterator it = find(rid.begin(), rid.end(), fcid);
			if (it == rid.end() && tmesh[fcid].type == 0)
				rid.push_back(fcid);
		}
	}
}

void TruncatedTspline::IdentifyTest1(vector<int> &rid)
{
	rid.clear();
	double tol(1.e-3);
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			double sumx(0.), sumy(0.);
			for (int j = 0; j < 4; j++)
			{
				sumx += cp[tmesh[i].cnct[j]].pm[0];
				sumy += cp[tmesh[i].cnct[j]].pm[1];
			}
			sumx /= 4.;
			sumy /= 4.;
			double df(sumx - sumy);
			if (df < 0.)
				df = -df;
			if (df < tol && sumx > 2. && sumx < 6. && sumy > 2. && sumy < 6.)
				//if(df<tol)
				rid.push_back(i);
		}
	}
	//rid.resize(3);
	//rid[0]=44; rid[1]=45; rid[2]=55;
}

void TruncatedTspline::OutputNumbers()
{
	//for(uint i=0;i<tmesh.size();i++)
	//{
	//	cout<<"eid: "<<i<<"\n";
	//	for(uint j=0;j<tmesh[i].IEN.size();j++)
	//	{
	//		cout<<tmesh[i].IEN[j]<<" ";
	//	}
	//	cout<<'\n';
	//	getchar();
	//}
	//for(int i=121;i<129;i++)
	//{
	//	cout<<"pid: "<<i<<'\n';
	//	cout<<"KU: \n";
	//	for(int j=0;j<5;j++)
	//	{
	//		cout<<cp[i].knotU[j]<<' ';
	//	}
	//	cout<<"\nKU: \n";
	//	for(int j=0;j<5;j++)
	//	{
	//		cout<<cp[i].knotV[j]<<' ';
	//	}
	//	cout<<'\n';
	//}
	//int eid(25);
	//cout<<"#IEN: "<<tmesh[eid].IEN.size()<<'\n';
	//for(uint i=0;i<tmesh[eid].IEN.size();i++)
	//{
	//	int pid=tmesh[eid].IEN[i];
	//	cout<<pid<<" U: "<<cp[pid].knotU[0]<<' '<<cp[pid].knotU[1]<<' '<<cp[pid].knotU[2]<<' '<<cp[pid].knotU[3]<<' '<<cp[pid].knotU[4]<<'\n';
	//	cout<<pid<<" V: "<<cp[pid].knotV[0]<<' '<<cp[pid].knotV[1]<<' '<<cp[pid].knotV[2]<<' '<<cp[pid].knotV[3]<<' '<<cp[pid].knotV[4]<<'\n';
	//	if(cp[pid].tbf.size()!=0)
	//	{
	//		cout<<"child: ";
	//		for(uint j=0;j<cp[pid].tbf.size();j++)
	//		{
	//			cout<<cp[pid].tbf[j]<<' ';
	//		}
	//		cout<<'\n';
	//	}
	//}
	//cout<<'\n';
	//for(uint pid=121;pid<cp.size();pid++)
	//{
	//	cout<<"pid: "<<pid<<'\n';
	//	cout<<"a: "<<cp[pid].knotU[0]<<' '<<cp[pid].knotU[1]<<' '<<cp[pid].knotU[2]<<' '<<cp[pid].knotU[3]<<' '<<cp[pid].knotU[4]<<'\n';
	//	cout<<"b: "<<cp[pid].knotV[0]<<' '<<cp[pid].knotV[1]<<' '<<cp[pid].knotV[2]<<' '<<cp[pid].knotV[3]<<' '<<cp[pid].knotV[4]<<'\n';
	//}
	//for(uint i=0;i<cp[38].tbf.size();i++)
	//{
	//	int pid=cp[38].tbf[i];
	//	cout<<"a: "<<cp[pid].knotU[0]<<' '<<cp[pid].knotU[1]<<' '<<cp[pid].knotU[2]<<' '<<cp[pid].knotU[3]<<' '<<cp[pid].knotU[4]<<'\n';
	//	cout<<"b: "<<cp[pid].knotV[0]<<' '<<cp[pid].knotV[1]<<' '<<cp[pid].knotV[2]<<' '<<cp[pid].knotV[3]<<' '<<cp[pid].knotV[4]<<'\n';
	//	cout<<cp[38].tc[i]<<'\n';
	//}
	//for(int i=104; i<108; i++)
	//{
	//	cout<<"pid "<<i<<": ";
	//	for(int j=0;j<4;j++)
	//		cout<<tmesh[i].cnct[j]<<' ';
	//	cout<<'\n';
	//	for(uint j=0;j<tmesh[i].IEN.size();j++)
	//		cout<<tmesh[i].IEN[j]<<' ';
	//	cout<<'\n';
	//}

	for (uint j = 0; j < tmesh[11].IEN.size(); j++)
	{
		//cout<<tmesh[11].IEN[j]<<' ';
		if (cp[tmesh[11].IEN[j]].tbf.size() != 0)
		{
			cout << tmesh[11].IEN[j] << '\n';
			for (uint i = 0; i < cp[tmesh[11].IEN[j]].tbf.size(); i++)
			{
				cout << cp[tmesh[11].IEN[j]].tbf[i] << ' ';
			}
			cout << '\n';
		}
	}
	cout << '\n';
}

void TruncatedTspline::ElementSubdivide_4(int eid)
{
	//find knots
	uint i, j;
	int cid[4] = {tmesh[eid].cnct[0], tmesh[eid].cnct[1], tmesh[eid].cnct[2], tmesh[eid].cnct[3]};
	double intvU((cp[cid[1]].knotU[2] - cp[cid[0]].knotU[2]) / 2.);
	double intvV((cp[cid[3]].knotV[2] - cp[cid[0]].knotV[2]) / 2.);
	double ktsU[7] = {cp[cid[0]].knotU[2] - 2. * intvU, cp[cid[0]].knotU[2] - intvU, cp[cid[0]].knotU[2], cp[cid[0]].knotU[2] + intvU, cp[cid[1]].knotU[2], cp[cid[1]].knotU[2] + intvU, cp[cid[1]].knotU[2] + 2. * intvU};
	double ktsV[7] = {cp[cid[0]].knotV[2] - 2. * intvV, cp[cid[0]].knotV[2] - intvV, cp[cid[0]].knotV[2], cp[cid[0]].knotV[2] + intvV, cp[cid[3]].knotV[2], cp[cid[3]].knotV[2] + intvV, cp[cid[3]].knotV[2] + 2. * intvV};
	if (cp[cid[0]].knotU[1] == cp[cid[0]].knotU[2])
	{
		ktsU[1] = cp[cid[0]].knotU[1];
		if (cp[cid[0]].knotU[0] != cp[cid[0]].knotU[1])
		{
			ktsU[0] = cp[cid[0]].knotU[1] - intvU;
		}
		else
		{
			ktsU[0] = cp[cid[0]].knotU[1];
		}
	}
	if (cp[cid[1]].knotU[2] == cp[cid[1]].knotU[3])
	{
		ktsU[5] = cp[cid[1]].knotU[3];
		if (cp[cid[1]].knotU[3] != cp[cid[1]].knotU[4])
		{
			ktsU[6] = cp[cid[1]].knotU[3] + intvU;
		}
		else
		{
			ktsU[6] = cp[cid[1]].knotU[3];
		}
	}
	if (cp[cid[0]].knotV[1] == cp[cid[0]].knotV[2])
	{
		ktsV[1] = cp[cid[0]].knotV[1];
		if (cp[cid[0]].knotV[0] != cp[cid[0]].knotV[1])
		{
			ktsV[0] = cp[cid[0]].knotV[1] - intvV;
		}
		else
		{
			ktsV[0] = cp[cid[0]].knotV[1];
		}
	}
	if (cp[cid[3]].knotV[2] == cp[cid[3]].knotV[3])
	{
		ktsV[5] = cp[cid[3]].knotV[3];
		if (cp[cid[3]].knotV[3] != cp[cid[3]].knotV[4])
		{
			ktsV[6] = cp[cid[3]].knotV[3] + intvV;
		}
		else
		{
			ktsV[6] = cp[cid[3]].knotV[3];
		}
	}
	int setkv[5][2] = {{3, 3}, {3, 2}, {4, 3}, {3, 4}, {2, 3}};
	vector<Vertex> pnew(5);
	for (i = 0; i < 5; i++)
	{
		for (j = 0; j < 5; j++)
		{
			pnew[i].knotU[j] = ktsU[setkv[i][0] + j - 2];
			pnew[i].knotV[j] = ktsU[setkv[i][1] + j - 2];
		}
	}

	int pid[5];
	int edid[12];
	cp.push_back(pnew[0]);
	pid[0] = cp.size() - 1;
	for (j = 0; j < 4; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[j], tmesh[eid].cnct[(j + 1) % 4]};
		if (tmedge[tmesh[eid].edge[j]].act == 1)
		{
			cp.push_back(pnew[j + 1]);
			pid[j + 1] = cp.size() - 1;
			tmedge[tmesh[eid].edge[j]].midpt = pid[j + 1];
			Edge edtmp1, edtmp2, edtmp3;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j + 1];
			edtmp1.len = tmedge[tmesh[eid].edge[j]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[j];
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j + 1];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[j]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[j];
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[(j + 1) % 4]].len / 2.;
			tmedge.push_back(edtmp1);
			edid[3 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].chd[0] = edid[3 * j];
			tmedge.push_back(edtmp2);
			edid[3 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].chd[1] = edid[3 * j + 1];
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].act = 0;
		}
		else
		{
			pid[j + 1] = tmedge[tmesh[eid].edge[j]].midpt;
			Edge edtmp3;
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[(j + 1) % 4]].len / 2.;
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			int ied(tmedge[tmesh[eid].edge[j]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[3 * j] = ied;
				edid[3 * j + 1] = tmedge[tmesh[eid].edge[j]].chd[1];
			}
			else
			{
				edid[3 * j] = tmedge[tmesh[eid].edge[j]].chd[1];
				edid[3 * j + 1] = ied;
			}
		}
	}

	int e_cnct[4][4] = {{tmesh[eid].cnct[0], pid[1], pid[0], pid[4]}, {pid[1], tmesh[eid].cnct[1], pid[2], pid[0]}, {pid[0], pid[2], tmesh[eid].cnct[2], pid[3]}, {pid[4], pid[0], pid[3], tmesh[eid].cnct[3]}};
	int e_edge[4][4] = {{edid[0], edid[2], edid[11], edid[10]}, {edid[1], edid[3], edid[5], edid[2]}, {edid[5], edid[4], edid[6], edid[8]}, {edid[11], edid[8], edid[7], edid[9]}};
	int enewid[4];
	vector<Element> etmp(4);
	for (i = 0; i < 4; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 0;
		etmp[i].lv = tmesh[eid].lv + 1.;
		etmp[i].prt = eid;
		for (j = 0; j < 4; j++)
		{
			etmp[i].cnct[j] = e_cnct[i][j];
			etmp[i].edge[j] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::Topo_Refine_Struct(const vector<int> &rfid, const vector<int> &rftype, vector<int> &rfid_more, vector<int> &rftype_more)
{
	//initialize for refinement
	npt_old = cp.size();
	nel_old = tmesh.size();
	rfid_more.clear();
	rftype_more.clear();
	rfid_more = rfid;
	rftype_more = rftype;
	for (uint i = 0; i < cp.size(); i++)
	{
		cp[i].update = 0;
		cp[i].aff = 0;
		cp[i].coortmp[0] = 0.;
		cp[i].coortmp[1] = 0.;
		cp[i].coortmp[2] = 0.;
		for (int j = 0; j < 4; j++)
		{
			cp[i].kitvUtmp[j] = 0.;
			cp[i].kitvVtmp[j] = 0.;
		}
		cp[i].truntmp = 0;
		vector<int>().swap(cp[i].tbftmp);
		vector<double>().swap(cp[i].tctmp);
	}
	for (uint i = 0; i < tmesh.size(); i++)
	{
		//tmesh[i].aff=0;
		vector<int>().swap(tmesh[i].IENtmp);
		vector<array<double, 5>>().swap(tmesh[i].patch_kutmp);
		vector<array<double, 5>>().swap(tmesh[i].patch_kvtmp);
	}
	//refine
	for (uint i = 0; i < rfid.size(); i++)
	{
		if (rftype[i] == 0)
		{
			ElementRefine_Unstruct_4(rfid[i]);
		}
		else if (rftype[i] == 1)
		{
			ElementRefine_Unstruct_2(rfid[i], rftype[i] - 1);
		}
		else if (rftype[i] == 2)
		{
			ElementRefine_Unstruct_2(rfid[i], rftype[i] - 1);
		}
		else if (rftype[i] == 3)
		{
			ElementRefine_Unstruct_b(rfid[i]);
		}
		else if (rftype[i] == 4)
		{
			ElementRefine_Unstruct_b(rfid[i]);
		}
		else
		{
			cout << "Other types of elements are not supported to be refined!\n";
			getchar();
		}
	}

	vector<int> rid_b;
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			if (tmesh[i].type == 2)
			{
				for (int j = 0; j < 4; j++)
				{
					if (tmedge[tmesh[i].edge[j]].act == 0)
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
	for (uint i = 0; i < rid_b.size(); i++)
	{
		ElementRefine_Unstruct_b(rid_b[i]);
	}
}

void TruncatedTspline::Geom_Refine_Struct(const vector<int> &rfid, const vector<int> &rftype)
{
}

void TruncatedTspline::ElementRefine_Struct_4(int eid)
{
	int pid[5];
	int edid[12];
	Vertex ptmp1;
	ptmp1.pm[0] = (cp[tmesh[eid].cnct[0]].pm[0] + cp[tmesh[eid].cnct[1]].pm[0]) / 2.;
	ptmp1.pm[1] = (cp[tmesh[eid].cnct[0]].pm[1] + cp[tmesh[eid].cnct[3]].pm[1]) / 2.;
	AssignIndex_NewFacePoint(eid, ptmp1);
	cp.push_back(ptmp1);
	pid[0] = cp.size() - 1;
	for (int j = 0; j < 4; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[j], tmesh[eid].cnct[(j + 1) % 4]};
		if (tmedge[tmesh[eid].edge[j]].act == 1)
		{
			Vertex ptmp;
			ptmp.pm[0] = (cp[itmp[0]].pm[0] + cp[itmp[1]].pm[0]) / 2.;
			ptmp.pm[1] = (cp[itmp[0]].pm[1] + cp[itmp[1]].pm[1]) / 2.;
			AssignIndex_NewEdgePoint(tmesh[eid].edge[j], ptmp);
			cp.push_back(ptmp);
			pid[j + 1] = cp.size() - 1;
			tmedge[tmesh[eid].edge[j]].midpt = pid[j + 1];
			Edge edtmp1, edtmp2, edtmp3;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j + 1];
			edtmp1.len = tmedge[tmesh[eid].edge[j]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[j];
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j + 1];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[j]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[j];
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[(j + 1) % 4]].len / 2.;
			tmedge.push_back(edtmp1);
			edid[3 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].chd[0] = edid[3 * j];
			tmedge.push_back(edtmp2);
			edid[3 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].chd[1] = edid[3 * j + 1];
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].act = 0;
		}
		else
		{
			pid[j + 1] = tmedge[tmesh[eid].edge[j]].midpt;
			Edge edtmp3;
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[(j + 1) % 4]].len / 2.;
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			int ied(tmedge[tmesh[eid].edge[j]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[3 * j] = ied;
				edid[3 * j + 1] = tmedge[tmesh[eid].edge[j]].chd[1];
			}
			else
			{
				edid[3 * j] = tmedge[tmesh[eid].edge[j]].chd[1];
				edid[3 * j + 1] = ied;
			}
		}
	}

	int e_cnct[4][4] = {{tmesh[eid].cnct[0], pid[1], pid[0], pid[4]}, {pid[1], tmesh[eid].cnct[1], pid[2], pid[0]}, {pid[0], pid[2], tmesh[eid].cnct[2], pid[3]}, {pid[4], pid[0], pid[3], tmesh[eid].cnct[3]}};
	int e_edge[4][4] = {{edid[0], edid[2], edid[11], edid[10]}, {edid[1], edid[3], edid[5], edid[2]}, {edid[5], edid[4], edid[6], edid[8]}, {edid[11], edid[8], edid[7], edid[9]}};
	int enewid[4];
	vector<Element> etmp(4);
	for (int i = 0; i < 4; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 0;
		etmp[i].lv = tmesh[eid].lv + 1.;
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[j] = e_cnct[i][j];
			etmp[i].edge[j] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::ElementRefine_Struct_3(int eid, int dir)
{
	int pid[4], edid[11], pos(dir);
	int cnid[4] = {pos, (pos + 1) % 4, (pos + 2) % 4, (pos + 3) % 4};
	Vertex ptmp1;
	ptmp1.pm[0] = (cp[tmesh[eid].cnct[0]].pm[0] + cp[tmesh[eid].cnct[1]].pm[0]) / 2.;
	ptmp1.pm[1] = (cp[tmesh[eid].cnct[0]].pm[1] + cp[tmesh[eid].cnct[3]].pm[1]) / 2.;
	AssignIndex_NewFacePoint(eid, ptmp1);
	cp.push_back(ptmp1);
	pid[0] = cp.size() - 1;
	for (int j = 0; j < 3; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[cnid[j]], tmesh[eid].cnct[cnid[j + 1]]};
		if (tmedge[tmesh[eid].edge[cnid[j]]].act == 1)
		{
			Vertex ptmp;
			ptmp.pm[0] = (cp[itmp[0]].pm[0] + cp[itmp[1]].pm[0]) / 2.;
			ptmp.pm[1] = (cp[itmp[0]].pm[1] + cp[itmp[1]].pm[1]) / 2.;
			AssignIndex_NewEdgePoint(tmesh[eid].edge[cnid[j]], ptmp);
			cp.push_back(ptmp);
			pid[j + 1] = cp.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].midpt = pid[j + 1];
			Edge edtmp1, edtmp2, edtmp3;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j + 1];
			edtmp1.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[cnid[j]];
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j + 1];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[cnid[j]];
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[cnid[j + 1]]].len / 2.;
			tmedge.push_back(edtmp1);
			edid[3 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[0] = edid[3 * j];
			tmedge.push_back(edtmp2);
			edid[3 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[1] = edid[3 * j + 1];
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].act = 0;
		}
		else
		{
			pid[j + 1] = tmedge[tmesh[eid].edge[cnid[j]]].midpt;
			Edge edtmp3;
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[cnid[j + 1]]].len / 2.;
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			int ied(tmedge[tmesh[eid].edge[cnid[j]]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[3 * j] = ied;
				edid[3 * j + 1] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
			}
			else
			{
				edid[3 * j] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
				edid[3 * j + 1] = ied;
			}
		}
	}
	edid[9] = tmesh[eid].edge[cnid[3]];
	Edge edtmp;
	edtmp.act = 0;
	edtmp.pt[0] = pid[1];
	edtmp.pt[1] = pid[3];
	edtmp.len = tmedge[edid[9]].len;
	edtmp.chd[0] = edid[2];
	edtmp.chd[1] = edid[8];
	edtmp.midpt = pid[0];
	tmedge.push_back(edtmp);
	edid[10] = tmedge.size() - 1;
	tmedge[edid[2]].prt = edid[10];
	tmedge[edid[8]].prt = edid[10];

	int e_cnct[3][4] = {{tmesh[eid].cnct[cnid[0]], pid[1], pid[3], tmesh[eid].cnct[cnid[3]]}, {pid[1], tmesh[eid].cnct[cnid[1]], pid[2], pid[0]}, {pid[0], pid[2], tmesh[eid].cnct[cnid[2]], pid[3]}};
	int e_edge[3][4] = {{edid[0], edid[10], edid[7], edid[9]}, {edid[1], edid[3], edid[5], edid[2]}, {edid[5], edid[4], edid[6], edid[8]}};
	int enewid[3];
	vector<Element> etmp(3);
	for (int i = 0; i < 3; i++)
	{
		etmp[i].act = 1;
		if (i == 0)
		{
			etmp[i].type = 0;
			etmp[i].lv = tmesh[eid].lv + 0.5;
		}
		else
		{
			etmp[i].type = 0;
			etmp[i].lv = tmesh[eid].lv + 1.;
		}
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[(pos + j) % 4] = e_cnct[i][j];
			etmp[i].edge[(pos + j) % 4] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::ElementRefine_Struct_2(int eid, int dir)
{
	int pid[2], edid[7], pos(dir);
	int cnid[2] = {pos, (pos + 2) % 4};
	for (int j = 0; j < 2; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[cnid[j]], tmesh[eid].cnct[(cnid[j] + 1) % 4]};
		if (tmedge[tmesh[eid].edge[cnid[j]]].act == 1)
		{
			Vertex ptmp;
			ptmp.pm[0] = (cp[itmp[0]].pm[0] + cp[itmp[1]].pm[0]) / 2.;
			ptmp.pm[1] = (cp[itmp[0]].pm[1] + cp[itmp[1]].pm[1]) / 2.;
			AssignIndex_NewEdgePoint(tmesh[eid].edge[cnid[j]], ptmp);
			cp.push_back(ptmp);
			pid[j] = cp.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].midpt = pid[j];
			Edge edtmp1, edtmp2, edtmp3;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j];
			edtmp1.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[cnid[j]];
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[cnid[j]];
			tmedge.push_back(edtmp1);
			edid[2 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[0] = edid[2 * j];
			tmedge.push_back(edtmp2);
			edid[2 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[1] = edid[2 * j + 1];
			tmedge[tmesh[eid].edge[cnid[j]]].act = 0;
		}
		else
		{
			pid[j] = tmedge[tmesh[eid].edge[cnid[j]]].midpt;
			int ied(tmedge[tmesh[eid].edge[cnid[j]]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[2 * j] = ied;
				edid[2 * j + 1] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
			}
			else
			{
				edid[2 * j] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
				edid[2 * j + 1] = ied;
			}
		}
	}
	Edge edtmp;
	edtmp.act = 1;
	edtmp.pt[0] = pid[0];
	edtmp.pt[1] = pid[1];
	edtmp.len = tmedge[tmesh[eid].edge[(pos + 1) % 4]].len;
	tmedge.push_back(edtmp);
	edid[4] = tmedge.size() - 1;
	edid[5] = tmesh[eid].edge[(pos + 1) % 4];
	edid[6] = tmesh[eid].edge[(pos + 3) % 4];

	int e_cnct[2][4] = {{tmesh[eid].cnct[pos], pid[0], pid[1], tmesh[eid].cnct[(pos + 3) % 4]}, {pid[0], tmesh[eid].cnct[(pos + 1) % 4], tmesh[eid].cnct[(pos + 2) % 4], pid[1]}};
	int e_edge[2][4] = {{edid[0], edid[4], edid[3], edid[6]}, {edid[1], edid[5], edid[2], edid[4]}};
	int enewid[2];
	vector<Element> etmp(2);
	for (int i = 0; i < 2; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 0;
		etmp[i].lv = tmesh[eid].lv + 0.5;
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[(pos + j) % 4] = e_cnct[i][j];
			etmp[i].edge[(pos + j) % 4] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::ElementRefine_Struct_b(int eid)
{
	int pos(0); //long edge position
	if (tmedge[tmesh[eid].edge[0]].len < tmedge[tmesh[eid].edge[1]].len)
		pos = 1;
	int ie[2] = {pos, pos + 2};
	int pid[2];
	int edid[7];
	for (int i = 0; i < 2; i++)
	{
		int itmp[2] = {tmedge[tmesh[eid].edge[ie[i]]].pt[0], tmedge[tmesh[eid].edge[ie[i]]].pt[1]};
		if (tmedge[tmesh[eid].edge[ie[i]]].act == 1)
		{
			Vertex ptmp1;
			ptmp1.pm[0] = (cp[itmp[0]].pm[0] + cp[itmp[1]].pm[0]) / 2.;
			ptmp1.pm[1] = (cp[itmp[0]].pm[1] + cp[itmp[1]].pm[1]) / 2.;
			AssignIndex_NewEdgePoint(tmesh[eid].edge[ie[i]], ptmp1);
			cp.push_back(ptmp1);
			pid[i] = cp.size() - 1;
			vector<Edge> edtmp(2);
			edtmp[0].act = 1;
			edtmp[0].prt = tmesh[eid].edge[ie[i]];
			edtmp[0].len = tmedge[tmesh[eid].edge[ie[i]]].len / 2.;
			edtmp[0].pt[0] = tmesh[eid].cnct[ie[i]];
			edtmp[0].pt[1] = pid[i];
			tmedge.push_back(edtmp[0]);
			edid[2 * i] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[ie[i]]].chd[0] = edid[2 * i];
			edtmp[1].act = 1;
			edtmp[1].prt = tmesh[eid].edge[ie[i]];
			edtmp[1].len = tmedge[tmesh[eid].edge[ie[i]]].len / 2.;
			edtmp[1].pt[0] = pid[i];
			edtmp[1].pt[1] = tmesh[eid].cnct[(ie[i] + 1) % 4];
			tmedge.push_back(edtmp[1]);
			edid[2 * i + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[ie[i]]].chd[1] = edid[2 * i + 1];
			tmedge[tmesh[eid].edge[ie[i]]].act = 0;
			tmedge[tmesh[eid].edge[ie[i]]].midpt = pid[i];
		}
		else
		{
			pid[i] = tmedge[tmesh[eid].edge[ie[i]]].midpt;
			int ied(tmedge[tmesh[eid].edge[ie[i]]].chd[0]);
			if (tmedge[ied].pt[0] == tmesh[eid].cnct[ie[i]] || tmedge[ied].pt[1] == tmesh[eid].cnct[ie[i]])
			{
				edid[2 * i] = ied;
				edid[2 * i + 1] = tmedge[tmesh[eid].edge[ie[i]]].chd[1];
			}
			else
			{
				edid[2 * i] = tmedge[tmesh[eid].edge[ie[i]]].chd[1];
				edid[2 * i + 1] = ied;
			}
		}
	}
	Edge edtmp;
	edtmp.act = 1;
	edtmp.len = tmedge[tmesh[eid].edge[ie[0] + 1]].len;
	edtmp.pt[0] = pid[0];
	edtmp.pt[1] = pid[1];
	tmedge.push_back(edtmp);
	edid[4] = tmedge.size() - 1;
	edid[5] = tmesh[eid].edge[ie[0] + 1];
	edid[6] = tmesh[eid].edge[(ie[1] + 1) % 4];

	int e_cnct[2][4] = {{tmesh[eid].cnct[ie[0]], pid[0], pid[1], tmesh[eid].cnct[(ie[1] + 1) % 4]}, {pid[0], tmesh[eid].cnct[ie[0] + 1], tmesh[eid].cnct[ie[1]], pid[1]}};
	int e_edge[2][4] = {{edid[0], edid[4], edid[3], edid[6]}, {edid[1], edid[5], edid[2], edid[4]}};
	int enewid[2];
	vector<Element> etmp(2);
	for (int i = 0; i < 2; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 2;
		etmp[i].lv = tmesh[eid].lv + 1.;
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[(pos + j) % 4] = e_cnct[i][j];
			etmp[i].edge[(pos + j) % 4] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::UpdateControlPoints_v0()
{
	int npt_new(cp.size() - npt_old);
	vector<vector<double>> cmat(npt_new, vector<double>(npt_old, 0.));
	vector<vector<double>> tmat(npt_new, vector<double>(npt_old, 0.));
	for (uint i = npt_old; i < cp.size(); i++)
	{
		for (uint j = 0; j < npt_old; j++)
		{
			if (cp[i].knotU[0] >= cp[j].knotU[0] && cp[i].knotU[4] <= cp[j].knotU[4] && cp[i].knotV[0] >= cp[j].knotV[0] && cp[i].knotV[4] <= cp[j].knotV[4])
			{
				vector<double> ku1(10), kv1(10);
				vector<vector<double>> Tu, Tv;
				vector<double> ku(cp[j].knotU, cp[j].knotU + 5);
				vector<double> kv(cp[j].knotV, cp[j].knotV + 5);
				vector<double>::iterator it1, it2;
				it1 = set_union(cp[j].knotU, cp[j].knotU + 5, cp[i].knotU, cp[i].knotU + 5, ku1.begin());
				it2 = set_union(cp[j].knotV, cp[j].knotV + 5, cp[i].knotV, cp[i].knotV + 5, kv1.begin());
				ku1.resize(it1 - ku1.begin());
				kv1.resize(it2 - kv1.begin());
				TMatrix(ku, ku1, 3, Tu);
				TMatrix(kv, kv1, 3, Tv);
				it1 = search(ku1.begin(), ku1.end(), cp[i].knotU, cp[i].knotU + 5);
				it2 = search(kv1.begin(), kv1.end(), cp[i].knotV, cp[i].knotV + 5);
				int loc1 = it1 - ku1.begin();
				int loc2 = it2 - kv1.begin();
				cmat[i - npt_old][j] = Tu[loc1][0] * Tv[loc2][0];
				tmat[i - npt_old][j] = cmat[i - npt_old][j];
			}
		}
	}

	for (uint i = 0; i < npt_new; i++)
	{
		for (uint j = 0; j < npt_old; j++)
		{
			if (cmat[i][j] != 0.)
			{
				for (uint k = 0; k < cp[j].tbf.size(); k++)
				{
					tmat[i][j] -= cmat[i][cp[j].tbf[k]] * cp[j].tc[k];
				}
			}
		}
	}

	for (uint i = 0; i < npt_new; i++)
	{
		int loc(npt_old + i);
		for (uint j = 0; j < npt_old; j++)
		{
			if (tmat[i][j] != 0.)
			{
				cp[loc].coor[0] += tmat[i][j] * cp[j].coor[0];
				cp[loc].coor[1] += tmat[i][j] * cp[j].coor[1];
				cp[loc].coor[2] += tmat[i][j] * cp[j].coor[2];
				cp[j].tbf.push_back(loc);
				cp[j].tc.push_back(tmat[i][j]);
			}
		}
	}

	//for(uint i=npt_old; i<cp.size(); i++)
	//{
	//	for(uint j=0; j<npt_old; j++)
	//	{
	//		if(cp[i].knotU[0]>=cp[j].knotU[0] && cp[i].knotU[4]<=cp[j].knotU[4] && cp[i].knotV[0]>=cp[j].knotV[0] && cp[i].knotV[4]<=cp[j].knotV[4])
	//		{
	//			double tmp=cmat[i-npt_old][j];
	//			for(uint k=0; k<cp[j].tbf.size(); k++)
	//			{
	//				vector<int>::iterator it=find(cp.begin(),cp.begin()+npt_old,cp[j].tbf[k]);
	//				if(it!=(cp.begin()+npt_old))
	//					tmp-=cmat[i-npt_old][it-cp.begin()]*cp[j].tc[k];
	//			}
	//			pnew[i].coor[0]+=tmp*cp[pid[j]].coor[0];
	//			pnew[i].coor[1]+=tmp*cp[pid[j]].coor[1];
	//			pnew[i].coor[2]+=tmp*cp[pid[j]].coor[2];
	//			cp[pid[j]].tc.push_back(tmp);
	//		}
	//	}
	//}
	//for(i=0; i<5; i++)
	//{
	//	if(pgn[i]>=npold)
	//		cp.push_back(pnew[i]);
	//	if(flag[i]==1)
	//		for(j=0;j<pid.size();j++)
	//			if(pnew[i].knotU[0]>=cp[pid[j]].knotU[0] && pnew[i].knotU[4]<=cp[pid[j]].knotU[4] && pnew[i].knotV[0]>=cp[pid[j]].knotV[0] && pnew[i].knotV[4]<=cp[pid[j]].knotV[4])
	//			{
	//				/*cp[pid[j]].tc.push_back(coef[i][j]);*/
	//				cp[pid[j]].tbf.push_back(pgn[i]);
	//			}
	//}
}

void TruncatedTspline::UpdateIEN_v0()
{
	uint i, j;
	for (i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			for (j = npt_old; j < cp.size(); j++)
			{
				if (cp[j].knotU[0] < cp[tmesh[i].cnct[1]].knotU[2] && cp[j].knotU[4] > cp[tmesh[i].cnct[0]].knotU[2] &&
					cp[j].knotV[0] < cp[tmesh[i].cnct[3]].knotV[2] && cp[j].knotV[4] > cp[tmesh[i].cnct[0]].knotV[2])
				{
					tmesh[i].IEN.push_back(j);
				}
			}
		}
	}
}

void TruncatedTspline::FindEdge_v0()
{
	uint i, j;
	tmedge.clear();
	for (i = 0; i < tmesh.size(); i++)
	{
		for (j = 0; j < 4; j++)
		{
			Edge edtmp;
			edtmp.act = 1;
			edtmp.pt[0] = tmesh[i].cnct[j];
			edtmp.pt[1] = tmesh[i].cnct[(j + 1) % 4];
			vector<Edge>::iterator it = find(tmedge.begin(), tmedge.end(), edtmp);
			int edid(it - tmedge.begin());
			if (it == tmedge.end())
			{
				tmedge.push_back(edtmp);
			}
			tmesh[i].edge[j] = edid;
			tmedge[edid].face.push_back(i);
		}
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].face.size() == 1)
		{
			int eid(tmedge[i].face[0]);
			tmesh[eid].act = 1;
			tmesh[eid].type = 2;
			int *it = find(tmesh[eid].edge, tmesh[eid].edge + 3, i);
			int loc(it - tmesh[eid].edge);
			tmedge[tmesh[eid].edge[(loc + 1) % 4]].len = 0.;
			tmedge[tmesh[eid].edge[(loc + 3) % 4]].len = 0.;
		}
	}
	//cornerEID.clear();
	//for(i=0; i<cp.size(); i++)
	//{
	//	if(cp[i].face.size()==1)
	//	{
	//		int eid(cp[i].face[0]);
	//		tmesh[eid].act=1;
	//		tmesh[eid].type=3;
	//		vector<int>::iterator it=find(cornerEID.begin(),cornerEID.end(),eid);
	//		if(it==cornerEID.end())
	//			cornerEID.push_back(eid);
	//	}
	//}
}

void TruncatedTspline::RefineTest_v0_1()
{
	FindEdge_v0();
	const int niter(3);
	for (int i = 0; i < niter; i++)
	{
		npt_old = cp.size();
		cout << "step: " << i << '\n';
		vector<int> rid;
		IdentifyTest(rid);
		cout << "Refining...\n";
		for (uint j = 0; j < rid.size(); j++)
		{
			ElementSubdivide_4(rid[j]);
		}
		UpdateControlPoints_v0();
		FindIENglb();

		cout << "Visualizing...\n";
		VisualizeTMesh("test6/tsp_test_CM_" + to_string(i));
		//cout<<"Visualizing...\n";
		//if(i==niter-1)
		//VisualizeVTK("test6/tsp_test_"+to_string(long long(i)));
	}
	//CollectActives();
	//cout<<"Bezier extracting...\n";
	//vector<BezierElement> bzmesh;
	//BezierExtract(bzmesh);
	//cout<<"Visualizing...\n";
	//BezierVTK("test2/tsp_test_bezier2",bzmesh);
}

bool TruncatedTspline::CheckSubKnotVector(const double ku1[5], const double kv1[5], const double ku2[5], const double kv2[5])
{
	if (ku1[0] >= ku2[0] && ku1[4] <= ku2[4] && kv1[0] >= kv2[0] && kv1[4] <= kv2[4])
	{
		double min_u1(10.), min_v1(10.), min_u2(10.), min_v2(10.);
		for (int i = 0; i < 4; i++)
		{
			double tmp = ku1[i + 1] - ku1[i];
			if (tmp != 0. && tmp < min_u1)
				min_u1 = tmp;
			tmp = kv1[i + 1] - kv1[i];
			if (tmp != 0. && tmp < min_v1)
				min_v1 = tmp;
			tmp = ku2[i + 1] - ku2[i];
			if (tmp != 0. && tmp < min_u2)
				min_u2 = tmp;
			tmp = kv2[i + 1] - kv2[i];
			if (tmp != 0. && tmp < min_v2)
				min_v2 = tmp;
		}
		if (min_u1 <= min_u2 && min_v1 <= min_v2)
		{
			for (int i = 0; i < 5; i++)
			{
				const double *it1 = find(ku2, ku2 + 5, ku1[i]);
				if (it1 == ku2 + 5)
				{
					return true;
				}
				const double *it2 = find(kv2, kv2 + 5, kv1[i]);
				if (it2 == kv2 + 5)
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

bool TruncatedTspline::CheckSubKnotVector(const array<double, 5> &ku1, const array<double, 5> &kv1, const array<double, 5> &ku2, const array<double, 5> &kv2)
{
	if (ku1[0] >= ku2[0] && ku1[4] <= ku2[4] && kv1[0] >= kv2[0] && kv1[4] <= kv2[4])
	{
		double min_u1(10.), min_v1(10.), min_u2(10.), min_v2(10.);
		for (int i = 0; i < 4; i++)
		{
			double tmp = ku1[i + 1] - ku1[i];
			if (tmp != 0. && tmp < min_u1)
				min_u1 = tmp;
			tmp = kv1[i + 1] - kv1[i];
			if (tmp != 0. && tmp < min_v1)
				min_v1 = tmp;
			tmp = ku2[i + 1] - ku2[i];
			if (tmp != 0. && tmp < min_u2)
				min_u2 = tmp;
			tmp = kv2[i + 1] - kv2[i];
			if (tmp != 0. && tmp < min_v2)
				min_v2 = tmp;
		}
		if (min_u1 <= min_u2 && min_v1 <= min_v2)
		{
			for (int i = 0; i < 5; i++)
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

int TruncatedTspline::CheckKnotInsertion(const double ku1[5], const double kv1[5], const double ku2[5], const double kv2[5])
{
	if (ku1[0] >= ku2[0] && ku1[4] <= ku2[4] && kv1[0] >= kv2[0] && kv1[4] <= kv2[4])
	{
		if (!equal(ku1, ku1 + 5, ku2) && equal(kv1, kv1 + 5, kv2))
		{
			if (NewKnots(ku1, ku2))
				return 0;
			else
				return -1;
		}
		else if (equal(ku1, ku1 + 5, ku2) && !equal(kv1, kv1 + 5, kv2))
		{
			if (NewKnots(kv1, kv2))
				return 1;
			else
				return -1;
		}
		else if (!equal(ku1, ku1 + 5, ku2) && !equal(kv1, kv1 + 5, kv2))
		{
			if (NewKnots(ku1, ku2) && !NewKnots(kv1, kv2))
				return 0;
			else if (!NewKnots(ku1, ku2) && NewKnots(kv1, kv2))
				return 1;
			else if (NewKnots(ku1, ku2) && NewKnots(kv1, kv2))
				return 2;
			else
				return -1;
		}
		else
		{
			return -1;
		}
	}
	else
	{
		return -1;
	}
}

bool TruncatedTspline::NewKnots(const double kv1[5], const double kv2[5])
{
	for (int i = 0; i < 5; i++)
	{
		const double *it = find(kv2, kv2 + 5, kv1[i]);
		if (it == kv2 + 5)
		{
			return true;
		}
	}
	return false;
}

void TruncatedTspline::PlotBasisFunctions(string fn, vector<int> &pid)
{
	vector<Vertex> spt;
	vector<double> sval;
	vector<Element> sele;
	vector<Vertex> lpt;		   //visulize parameter lines
	vector<array<int, 2>> led; //line connectivity
	int ns(5), ecount(0), loc0, loc1, loc2;
	vector<double> su(ns), sv(ns);

	for (uint e = 0; e < tmesh.size(); e++)
	{
		if (tmesh[e].act == 1)
		{
			int bfid(-1);
			for (uint i = 0; i < tmesh[e].IEN.size(); i++)
			{
				vector<int>::iterator it = find(pid.begin(), pid.end(), tmesh[e].IEN[i]);
				if (it != pid.end())
				{
					bfid = tmesh[e].IEN[i];
					break;
				}
			}
			for (int i = 0; i < ns; i++)
			{
				loc0 = tmesh[e].cnct[0];
				loc1 = tmesh[e].cnct[1];
				loc2 = tmesh[e].cnct[3];
				su[i] = cp[loc0].knotU[2] + i * (cp[loc1].knotU[2] - cp[loc0].knotU[2]) / (ns - 1);
				sv[i] = cp[loc0].knotV[2] + i * (cp[loc2].knotV[2] - cp[loc0].knotV[2]) / (ns - 1);
			}

			int loc(0);
			for (int a = 0; a < ns; a++)
			{
				for (int b = 0; b < ns; b++)
				{
					Vertex pt;
					FunctionMagnitude(e, bfid, su[b], sv[a], pt);
					spt.push_back(pt);
					if (a == 0 || a == ns - 1 || b == 0 || b == ns - 1)
					{
						lpt.push_back(pt);
					}
					//double sum=PartitionOfUnity(e,su[b],sv[a]);
					//sval.push_back(sum);
				}
			}
			for (int a = 0; a < ns - 1; a++)
			{
				for (int b = 0; b < ns - 1; b++)
				{
					Element el;
					el.cnct[0] = ecount * ns * ns + a * ns + b;
					el.cnct[1] = ecount * ns * ns + a * ns + b + 1;
					el.cnct[2] = ecount * ns * ns + (a + 1) * ns + b + 1;
					el.cnct[3] = ecount * ns * ns + (a + 1) * ns + b;
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
	}

	string fname = fn + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i < spt.size(); i++)
		{
			fout << spt[i].coor[0] << " " << spt[i].coor[1] << " " << spt[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i].cnct[0] << " " << sele[i].cnct[1] << " " << sele[i].cnct[2] << " " << sele[i].cnct[3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
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
			fout1 << lpt[i].coor[0] << " " << lpt[i].coor[1] << " " << lpt[i].coor[2] << "\n";
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

void TruncatedTspline::PlotTrunBasisFunctions(string fn, vector<int> &pid)
{
	vector<Vertex> spt;
	vector<double> sval;
	vector<Element> sele;
	vector<Vertex> lpt;		   //visulize parameter lines
	vector<array<int, 2>> led; //line connectivity
	int ns(5), ecount(0), loc0, loc1, loc2;
	vector<double> su(ns), sv(ns);

	for (uint e = 0; e < tmesh.size(); e++)
	{
		if (tmesh[e].act == 1)
		{
			int bfid(-1);
			for (uint i = 0; i < tmesh[e].IEN.size(); i++)
			{
				vector<int>::iterator it = find(pid.begin(), pid.end(), tmesh[e].IEN[i]);
				if (it != pid.end())
				{
					bfid = tmesh[e].IEN[i];
					break;
				}
			}
			for (int i = 0; i < ns; i++)
			{
				loc0 = tmesh[e].cnct[0];
				loc1 = tmesh[e].cnct[1];
				loc2 = tmesh[e].cnct[3];
				su[i] = cp[loc0].knotU[2] + i * (cp[loc1].knotU[2] - cp[loc0].knotU[2]) / (ns - 1);
				sv[i] = cp[loc0].knotV[2] + i * (cp[loc2].knotV[2] - cp[loc0].knotV[2]) / (ns - 1);
			}

			int loc(0);
			for (int a = 0; a < ns; a++)
			{
				for (int b = 0; b < ns; b++)
				{
					Vertex pt;
					TrunFunctionMagnitude(e, bfid, su[b], sv[a], pt);
					spt.push_back(pt);
					if (a == 0 || a == ns - 1 || b == 0 || b == ns - 1)
					{
						lpt.push_back(pt);
					}
					//double sum=PartitionOfUnity(e,su[b],sv[a]);
					//sval.push_back(sum);
				}
			}
			for (int a = 0; a < ns - 1; a++)
			{
				for (int b = 0; b < ns - 1; b++)
				{
					Element el;
					el.cnct[0] = ecount * ns * ns + a * ns + b;
					el.cnct[1] = ecount * ns * ns + a * ns + b + 1;
					el.cnct[2] = ecount * ns * ns + (a + 1) * ns + b + 1;
					el.cnct[3] = ecount * ns * ns + (a + 1) * ns + b;
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
	}

	string fname = fn + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i < spt.size(); i++)
		{
			fout << spt[i].coor[0] << " " << spt[i].coor[1] << " " << spt[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i].cnct[0] << " " << sele[i].cnct[1] << " " << sele[i].cnct[2] << " " << sele[i].cnct[3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
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
			fout1 << lpt[i].coor[0] << " " << lpt[i].coor[1] << " " << lpt[i].coor[2] << "\n";
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

void TruncatedTspline::FunctionMagnitude(int eid, int pid, double u, double v, Vertex &pt)
{
	unsigned int i, j;
	int loc, loc1;
	double Ni;
	pt.coor[0] = u;
	pt.coor[1] = v;
	pt.coor[2] = 0.;
	vector<double> ku(5, 0.), kv(5, 0.), uval, vval;
	BSplineBasis bu, bv;
	for (i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		if (tmesh[eid].IEN[i] == pid)
		{
			loc = tmesh[eid].IEN[i];
			ku.assign(cp[loc].knotU, cp[loc].knotU + 5);
			kv.assign(cp[loc].knotV, cp[loc].knotV + 5);
			bu.Set(3, ku);
			bv.Set(3, kv);
			bu.BasisFunction(0, u, 0, uval);
			bv.BasisFunction(0, v, 0, vval);
			Ni = uval[0] * vval[0];
			//for(j=0;j<cp[loc].tbf.size();j++)
			//{
			//	loc1=cp[loc].tbf[j];
			//	ku.assign(cp[loc1].knotU,cp[loc1].knotU+5);
			//	kv.assign(cp[loc1].knotV,cp[loc1].knotV+5);
			//	bu.Set(3,ku);
			//	bv.Set(3,kv);
			//	bu.BasisFunction(0,u,0,uval);
			//	bv.BasisFunction(0,v,0,vval);
			//	Ni-=cp[loc].tc[j]*uval[0]*vval[0];
			//}
			//pt.coor[0]+=Ni*cp[loc].coor[0];
			//pt.coor[1]+=Ni*cp[loc].coor[1];
			//pt.coor[2]+=Ni*cp[loc].coor[2];
			pt.coor[2] = Ni;
		}
	}
}

void TruncatedTspline::TrunFunctionMagnitude(int eid, int pid, double u, double v, Vertex &pt)
{
	unsigned int i, j;
	int loc, loc1;
	double Ni;
	pt.coor[0] = u;
	pt.coor[1] = v;
	pt.coor[2] = 0.;
	vector<double> ku(5, 0.), kv(5, 0.), uval, vval;
	BSplineBasis bu, bv;
	for (i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		if (tmesh[eid].IEN[i] == pid)
		{
			loc = tmesh[eid].IEN[i];
			ku.assign(cp[loc].knotU, cp[loc].knotU + 5);
			kv.assign(cp[loc].knotV, cp[loc].knotV + 5);
			bu.Set(3, ku);
			bv.Set(3, kv);
			bu.BasisFunction(0, u, 0, uval);
			bv.BasisFunction(0, v, 0, vval);
			Ni = uval[0] * vval[0];
			for (j = 0; j < cp[loc].tbf.size(); j++)
			{
				loc1 = cp[loc].tbf[j];
				ku.assign(cp[loc1].knotU, cp[loc1].knotU + 5);
				kv.assign(cp[loc1].knotV, cp[loc1].knotV + 5);
				bu.Set(3, ku);
				bv.Set(3, kv);
				bu.BasisFunction(0, u, 0, uval);
				bv.BasisFunction(0, v, 0, vval);
				Ni -= cp[loc].tc[j] * uval[0] * vval[0];
			}
			//pt.coor[0]+=Ni*cp[loc].coor[0];
			//pt.coor[1]+=Ni*cp[loc].coor[1];
			//pt.coor[2]+=Ni*cp[loc].coor[2];
			pt.coor[2] = Ni;
		}
	}
}

void TruncatedTspline::InitialConnect()
{
	uint i, j;
	tmedge.clear();
	for (i = 0; i < tmesh.size(); i++)
	{
		tmesh[i].act = 1;
		tmesh[i].type = 0;
		for (j = 0; j < 4; j++)
		{
			cp[tmesh[i].cnct[j]].face.push_back(i);
		}
	}
	for (i = 0; i < cp.size(); i++)
	{
		if (cp[i].face.size() == 3 || cp[i].face.size() > 4)
		{
			cp[i].type = 2;
		}
	}
	for (i = 0; i < tmesh.size(); i++)
	{
		int flag(0), pos(0);
		for (j = 0; j < 4; j++)
		{
			if (cp[tmesh[i].cnct[j]].type == 2)
			{
				flag = 1;
				pos = j;
				break;
			}
		}
		if (flag == 1)
		{
			tmesh[i].type = 4;
			int cnctnew[4] = {tmesh[i].cnct[pos], tmesh[i].cnct[(pos + 1) % 4], tmesh[i].cnct[(pos + 2) % 4], tmesh[i].cnct[(pos + 3) % 4]};
			for (j = 0; j < 4; j++)
			{
				tmesh[i].cnct[j] = cnctnew[j];
			}
		}
	}
	for (i = 0; i < tmesh.size(); i++)
	{
		for (j = 0; j < 4; j++)
		{
			Edge edtmp;
			edtmp.act = 1;
			edtmp.pt[0] = tmesh[i].cnct[j];
			edtmp.pt[1] = tmesh[i].cnct[(j + 1) % 4];
			edtmp.len = 1.;
			vector<Edge>::iterator it = find(tmedge.begin(), tmedge.end(), edtmp);
			int edid(it - tmedge.begin());
			if (it == tmedge.end())
			{
				tmedge.push_back(edtmp);
			}
			tmesh[i].edge[j] = edid;
			tmedge[edid].face.push_back(i);
		}
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		cp[tmedge[i].pt[0]].edge.push_back(i);
		cp[tmedge[i].pt[1]].edge.push_back(i);
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].face.size() == 1)
		{
			int eid(tmedge[i].face[0]);
			if (tmesh[eid].type == 0)
				tmesh[eid].type = 2;
			else if (tmesh[eid].type == 2)
				tmesh[eid].type = 3;
			int *it = find(tmesh[eid].edge, tmesh[eid].edge + 4, i);
			int loc(it - tmesh[eid].edge);
			tmedge[tmesh[eid].edge[(loc + 1) % 4]].len = 0.;
			tmedge[tmesh[eid].edge[(loc + 3) % 4]].len = 0.;
		}
	}

	FindEdgeTopoDirec_1();
	FindKnotInterval_1();
	UpdateKnotInterval_1();
	SetLocalCoorSystem();
	FindIEN_1();
}

void TruncatedTspline::InitialConnect_1()
{
	uint i, j;
	tmedge.clear();
	for (i = 0; i < tmesh.size(); i++)
	{
		tmesh[i].act = 1;
		tmesh[i].type = 0;
		for (j = 0; j < 4; j++)
		{
			cp[tmesh[i].cnct[j]].face.push_back(i);
		}
	}
	for (i = 0; i < cp.size(); i++)
	{
		if (cp[i].face.size() == 3 || cp[i].face.size() > 4)
		{
			cp[i].type = 2;
		}
	}
	for (i = 0; i < tmesh.size(); i++)
	{
		int flag(0), pos(0);
		for (j = 0; j < 4; j++)
		{
			if (cp[tmesh[i].cnct[j]].type == 2)
			{
				flag = 1;
				pos = j;
				break;
			}
		}
		if (flag == 1)
		{
			tmesh[i].type = 4;
			int cnctnew[4] = {tmesh[i].cnct[pos], tmesh[i].cnct[(pos + 1) % 4], tmesh[i].cnct[(pos + 2) % 4], tmesh[i].cnct[(pos + 3) % 4]};
			for (j = 0; j < 4; j++)
			{
				tmesh[i].cnct[j] = cnctnew[j];
			}
		}
	}
	for (i = 0; i < tmesh.size(); i++)
	{
		for (j = 0; j < 4; j++)
		{
			Edge edtmp;
			edtmp.act = 1;
			edtmp.pt[0] = tmesh[i].cnct[j];
			edtmp.pt[1] = tmesh[i].cnct[(j + 1) % 4];
			edtmp.len = 1.;
			vector<Edge>::iterator it = find(tmedge.begin(), tmedge.end(), edtmp);
			int edid(it - tmedge.begin());
			if (it == tmedge.end())
			{
				tmedge.push_back(edtmp);
			}
			tmesh[i].edge[j] = edid;
			tmedge[edid].face.push_back(i);
		}
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		cp[tmedge[i].pt[0]].edge.push_back(i);
		cp[tmedge[i].pt[1]].edge.push_back(i);
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].face.size() == 1)
		{
			int eid(tmedge[i].face[0]);
			if (tmesh[eid].type == 0)
				tmesh[eid].type = 2;
			else if (tmesh[eid].type == 2)
				tmesh[eid].type = 3;
			int *it = find(tmesh[eid].edge, tmesh[eid].edge + 4, i);
			int loc(it - tmesh[eid].edge);
			tmedge[tmesh[eid].edge[(loc + 1) % 4]].len = 0.;
			tmedge[tmesh[eid].edge[(loc + 3) % 4]].len = 0.;
		}
	}

	//FindEdgeTopoDirec_1();
	//FindKnotInterval_1();
	//UpdateKnotInterval_1();
	//SetLocalCoorSystem();
	//FindIEN_1();
}

void TruncatedTspline::InitialConnect_2()
{
	uint i, j;
	tmedge.clear();
	for (i = 0; i < tmesh.size(); i++)
	{
		tmesh[i].act = 1;
		tmesh[i].type = 0;
		for (j = 0; j < 4; j++)
		{
			cp[tmesh[i].cnct[j]].face.push_back(i);
		}
	}
	for (i = 0; i < cp.size(); i++)
	{
		if (cp[i].face.size() == 3 || cp[i].face.size() > 4)
		{
			cp[i].type = 2;
		}
	}
	for (i = 0; i < tmesh.size(); i++)
	{
		int flag(0), pos(0);
		for (j = 0; j < 4; j++)
		{
			if (cp[tmesh[i].cnct[j]].type == 2)
			{
				flag = 1;
				pos = j;
				break;
			}
		}
		if (flag == 1)
		{
			tmesh[i].type = 4;
			int cnctnew[4] = {tmesh[i].cnct[pos], tmesh[i].cnct[(pos + 1) % 4], tmesh[i].cnct[(pos + 2) % 4], tmesh[i].cnct[(pos + 3) % 4]};
			for (j = 0; j < 4; j++)
			{
				tmesh[i].cnct[j] = cnctnew[j];
			}
		}
		for (j = 1; j < 4; j++)
		{
			if (cp[tmesh[i].cnct[j]].type == 2)
			{
				tmesh[i].type = 5;
			}
		}
	}
	for (i = 0; i < tmesh.size(); i++)
	{
		for (j = 0; j < 4; j++)
		{
			Edge edtmp;
			edtmp.act = 1;
			edtmp.pt[0] = tmesh[i].cnct[j];
			edtmp.pt[1] = tmesh[i].cnct[(j + 1) % 4];
			edtmp.len = 1.;
			vector<Edge>::iterator it = find(tmedge.begin(), tmedge.end(), edtmp);
			int edid(it - tmedge.begin());
			if (it == tmedge.end())
			{
				tmedge.push_back(edtmp);
			}
			tmesh[i].edge[j] = edid;
			tmedge[edid].face.push_back(i);
		}
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		cp[tmedge[i].pt[0]].edge.push_back(i);
		cp[tmedge[i].pt[1]].edge.push_back(i);
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].face.size() == 1)
		{
			int eid(tmedge[i].face[0]);
			if (tmesh[eid].type == 0)
				tmesh[eid].type = 2;
			else if (tmesh[eid].type == 2)
				tmesh[eid].type = 3;
			int *it = find(tmesh[eid].edge, tmesh[eid].edge + 4, i);
			int loc(it - tmesh[eid].edge);
			tmedge[tmesh[eid].edge[(loc + 1) % 4]].len = 0.;
			tmedge[tmesh[eid].edge[(loc + 3) % 4]].len = 0.;
		}
	}
}

void TruncatedTspline::UpdateConnect()
{
	uint i, j, k;
	for (i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			cp[i].face.clear();
			cp[i].edge.clear();
		}
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].act == 1)
		{
			tmedge[i].face.clear();
		}
	}
	//loop all faces
	for (i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			for (j = 0; j < 4; j++)
			{
				cp[tmesh[i].cnct[j]].face.push_back(i);
				if (tmedge[tmesh[i].edge[j]].act == 1)
				{
					tmedge[tmesh[i].edge[j]].face.push_back(i);
				}
				else
				{
					cp[tmedge[tmesh[i].edge[j]].midpt].face.push_back(i);
					int chdid[2] = {tmedge[tmesh[i].edge[j]].chd[0], tmedge[tmesh[i].edge[j]].chd[1]};
					if (tmedge[chdid[0]].act == 1 && tmedge[chdid[1]].act == 1)
					{
						tmedge[chdid[0]].face.push_back(i);
						tmedge[chdid[1]].face.push_back(i);
					}
					else
					{
						cerr << "Configuration not recognized!\n";
					}
				}
			}
		}
	}
	//loop all edges
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].act == 1)
		{
			cp[tmedge[i].pt[0]].edge.push_back(i);
			cp[tmedge[i].pt[1]].edge.push_back(i);
		}
	}

	for (i = 0; i < cp.size(); i++)
	{
		if (cp[i].face.size() == 3 && cp[i].type != 2)
		{
			cp[i].type = 1;
		}
	}

	//update corner element level
	//for(i=0; i<cornerEID.size(); i++)
	//{
	//	double lev(tmesh[cornerEID[i]].lv);
	//	for(j=0; j<4; j++)
	//	{
	//		for(k=0; k<cp[tmesh[cornerEID[i]].cnct[j]].face.size(); k++)
	//		{
	//			if(cp[tmesh[cornerEID[i]].cnct[j]].face[k] != cornerEID[k])
	//				lev=tmesh[cp[tmesh[cornerEID[i]].cnct[j]].face[k]].lv;
	//		}
	//	}
	//	tmesh[cornerEID[i]].lv=lev;
	//}

	FindEdgeTopoDirec_1();
	FindKnotInterval_1(); //find kitvtmp
	UpdateKnotInterval_1();
	SetLocalCoorSystem();
	FindIEN_1();
}

void TruncatedTspline::UpdateConnect_1()
{
	uint i, j, k;
	for (i = 0; i < cp.size(); i++)
	{
		//if(cp[i].act==1)
		{
			cp[i].face.clear();
			cp[i].edge.clear();
		}
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		//if(tmedge[i].act==1)
		{
			tmedge[i].face.clear();
		}
	}
	//loop all faces
	for (i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			for (j = 0; j < 4; j++)
			{
				cp[tmesh[i].cnct[j]].face.push_back(i);
				if (tmedge[tmesh[i].edge[j]].act == 1)
				{
					tmedge[tmesh[i].edge[j]].face.push_back(i);
				}
				else
				{
					cp[tmedge[tmesh[i].edge[j]].midpt].face.push_back(i);
					int chdid[2] = {tmedge[tmesh[i].edge[j]].chd[0], tmedge[tmesh[i].edge[j]].chd[1]};
					if (tmedge[chdid[0]].act == 1 && tmedge[chdid[1]].act == 1)
					{
						tmedge[chdid[0]].face.push_back(i);
						tmedge[chdid[1]].face.push_back(i);
					}
					else
					{
						cerr << "Configuration not recognized!\n";
					}
				}
			}
		}
	}
	//loop all edges
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].act == 1)
		{
			cp[tmedge[i].pt[0]].edge.push_back(i);
			cp[tmedge[i].pt[1]].edge.push_back(i);
		}
	}

	for (i = 0; i < cp.size(); i++)
	{
		if (cp[i].face.size() == 3 && cp[i].type != 2)
		{
			cp[i].type = 1;
		}
	}

	//FindEdgeTopoDirec_1();
	//FindKnotInterval_1();//find kitvtmp
	//UpdateKnotInterval_1();
	//SetLocalCoorSystem();
	//FindIEN_1();
}

void TruncatedTspline::FindEdgeTopoDirec()
{
	for (uint i = 0; i < tmedge.size(); i++)
	{
		tmedge[i].pn[0][0] = 0;
		tmedge[i].pn[0][1] = -1;
		tmedge[i].pn[1][0] = 0;
		tmedge[i].pn[1][1] = -1; //initialize as end
		if (tmedge[i].act == 1)
		{
			for (int j = 0; j < 2; j++)
			{
				if (cp[tmedge[i].pt[j]].type == 0) //regular or T-junctions
				{
					for (uint k = 0; k < cp[tmedge[i].pt[j]].edge.size(); k++)
					{
						if (cp[tmedge[i].pt[j]].edge[k] != i)
						{
							int flag(0);
							for (uint i1 = 0; i1 < tmedge[i].face.size(); i1++)
							{
								for (uint j1 = 0; j1 < tmesh[tmedge[i].face[i1]].edge_act.size(); j1++)
								{
									for (uint k1 = 0; k1 < tmesh[tmedge[i].face[i1]].edge_act[j1].size(); k1++)
									{
										if (tmesh[tmedge[i].face[i1]].edge_act[j1][k1] == cp[tmedge[i].pt[j]].edge[k])
										{
											flag = 1;
											break;
										}
									}
								}
							}
							if (flag == 0)
							{
								tmedge[i].pn[j][0] = 0;
								tmedge[i].pn[j][1] = cp[tmedge[i].pt[j]].edge[k];
								break;
							}
						}
					}
				}
				else if (cp[tmedge[i].pt[j]].type == 1) //T-junctions
				{
					int fid;
					for (uint k = 0; k < cp[tmedge[i].pt[j]].face.size(); k++)
					{
						int flag(0);
						for (uint i1 = 0; i1 < tmedge[i].face.size(); i1++)
						{
							if (cp[tmedge[i].pt[j]].face[k] == tmedge[i].face[i1])
							{
								flag = 1;
								break;
							}
						}
						if (flag == 0)
						{
							fid = cp[tmedge[i].pt[j]].face[k];
							break;
						}
					}
					int loc, flag(0);
					for (uint k = 0; k < tmesh[fid].edge_act.size(); k++)
					{
						vector<int>::iterator it = find(tmesh[fid].edge_act[k].begin(), tmesh[fid].edge_act[k].end(), i);
						if (it != tmesh[fid].edge_act[k].end())
						{
							flag = 1;
							loc = k;
							break;
						}
					}
					if (flag == 1)
					{
						for (uint k = 0; k < cp[tmedge[i].pt[j]].edge.size(); k++)
						{
							if (cp[tmedge[i].pt[j]].edge[k] != i)
							{
								vector<int>::iterator it = find(tmesh[fid].edge_act[loc].begin(), tmesh[fid].edge_act[loc].end(), cp[tmedge[i].pt[j]].edge[k]);
								if (it != tmesh[fid].edge_act[loc].end())
								{
									tmedge[i].pn[j][0] = 0;
									tmedge[i].pn[j][0] = cp[tmedge[i].pt[j]].edge[k];
									break;
								}
							}
						}
					}
					else if (flag == 0)
					{
						tmedge[i].pn[j][0] = 1;
						tmedge[i].pn[j][1] = fid;
					}
				}
				else if (cp[tmedge[i].pt[j]].type == 2) //extraordinary
				{
					tmedge[i].pn[j][0] = 2;
				}
			}
		}
	}
}

void TruncatedTspline::FindEdgeTopoDirec_1()
{
	for (uint i = 0; i < tmedge.size(); i++)
	{
		tmedge[i].pn[0][0] = 3;
		tmedge[i].pn[0][1] = -1;
		tmedge[i].pn[1][0] = 3;
		tmedge[i].pn[1][1] = -1; //initialize as end
		if (tmedge[i].act == 1)
		{
			for (int j = 0; j < 2; j++)
			{
				if (cp[tmedge[i].pt[j]].type == 0) //regular
				{
					for (uint k = 0; k < cp[tmedge[i].pt[j]].edge.size(); k++)
					{
						if (cp[tmedge[i].pt[j]].edge[k] != i)
						{
							int flag(0);
							for (uint i1 = 0; i1 < tmedge[i].face.size(); i1++)
							{
								for (uint j1 = 0; j1 < 4; j1++)
								{
									int ed(tmesh[tmedge[i].face[i1]].edge[j1]);
									if (tmedge[ed].act == 1 && ed == cp[tmedge[i].pt[j]].edge[k])
									{
										flag = 1;
										break;
									}
									else if (tmedge[ed].act == 0 && (tmedge[ed].chd[0] == cp[tmedge[i].pt[j]].edge[k] || tmedge[ed].chd[1] == cp[tmedge[i].pt[j]].edge[k]))
									{
										flag = 1;
										break;
									}
								}
							}
							if (flag == 0)
							{
								tmedge[i].pn[j][0] = 0;
								tmedge[i].pn[j][1] = cp[tmedge[i].pt[j]].edge[k];
								break;
							}
						}
					}
				}
				else if (cp[tmedge[i].pt[j]].type == 1) //T-junctions
				{
					int fid(-1);
					/*for(uint k=0; k<cp[tmedge[i].pt[j]].face.size(); k++)
					{
						int flag(0);
						for(uint i1=0; i1<tmedge[i].face.size(); i1++)
						{
							if(cp[tmedge[i].pt[j]].face[k]==tmedge[i].face[i1])
							{
								flag=1; break;
							}
						}
						if(flag==0)
						{
							fid=cp[tmedge[i].pt[j]].face[k]; break;
						}
					}
					int flag(0);
					for(uint k=0; k<4; k++)
					{
						if(tmedge[tmesh[fid].edge[k]].act==0 && (tmedge[tmesh[fid].edge[k]].chd[0]==i || tmedge[tmesh[fid].edge[k]].chd[1]==i))
						{
							flag=1;
							tmedge[i].pn[j][0]=0;
							if(tmedge[tmesh[fid].edge[k]].chd[0]==i)
								tmedge[i].pn[j][1]=tmedge[tmesh[fid].edge[k]].chd[1];
							else
								tmedge[i].pn[j][1]=tmedge[tmesh[fid].edge[k]].chd[0];
							break;
						}
					}
					if(flag==0)
					{
						tmedge[i].pn[j][0]=1;
						tmedge[i].pn[j][1]=fid;
					}*/
					int loc(0);
					for (uint k = 0; k < cp[tmedge[i].pt[j]].face.size(); k++)
					{
						int ftmp(cp[tmedge[i].pt[j]].face[k]);
						for (int k1 = 0; k1 < 4; k1++)
						{
							if (tmedge[tmesh[ftmp].edge[k1]].act == 0 && tmedge[tmesh[ftmp].edge[k1]].midpt == tmedge[i].pt[j])
							{
								loc = k1;
								fid = ftmp;
								break;
							}
						}
						if (fid != -1)
						{
							break;
						}
					}
					if (fid == -1)
					{
						cout << "edge id: " << i << "\n";
						cout << tmedge[i].pt[0] << " " << tmedge[i].pt[1] << "\n";
						cout << tmedge[i].prt << "\n";
						//cout<<cp[tmedge[i].pt[0]].face.size()<<" "<<cp[tmedge[i].pt[1]].face.size()<<"\n";
						cout << cp[tmedge[i].pt[0]].face[0] << " " << cp[tmedge[i].pt[0]].face[1] << " " << cp[tmedge[i].pt[0]].face[2] << "\n";
						cerr << "T-junction cannot be found in any neighboring elements!\n";
						getchar();
					}
					if (tmedge[tmesh[fid].edge[loc]].chd[0] == i)
					{
						tmedge[i].pn[j][0] = 0;
						tmedge[i].pn[j][1] = tmedge[tmesh[fid].edge[loc]].chd[1];
					}
					else if (tmedge[tmesh[fid].edge[loc]].chd[1] == i)
					{
						tmedge[i].pn[j][0] = 0;
						tmedge[i].pn[j][1] = tmedge[tmesh[fid].edge[loc]].chd[0];
					}
					else
					{
						tmedge[i].pn[j][0] = 1;
						tmedge[i].pn[j][1] = fid;
					}
				}
				else if (cp[tmedge[i].pt[j]].type == 2) //extraordinary
				{
					tmedge[i].pn[j][0] = 2;
				}
			}
		}
	}
}

void TruncatedTspline::FindKnotInterval_1()
{
	for (uint i = 0; i < cp.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[i].kitvUtmp[j] = 1.;
			cp[i].kitvVtmp[j] = 1.;
		}
		if (cp[i].type != 2)
		{
			int pos(0);
			cp[i].rfc = -1;
			for (uint j = 0; j < cp[i].face.size(); j++)
			{
				int *it = find(tmesh[cp[i].face[j]].cnct, tmesh[cp[i].face[j]].cnct + 4, i);
				if (it != tmesh[cp[i].face[j]].cnct + 4)
				{
					cp[i].rfc = cp[i].face[j];
					pos = it - tmesh[cp[i].face[j]].cnct;
					break;
				}
			}
			if (cp[i].rfc == -1)
			{
				cerr << "Cannot find correct reference face!\n";
				getchar();
			}
			cp[i].uved[0] = tmesh[cp[i].rfc].edge[pos];
			cp[i].uved[1] = tmesh[cp[i].rfc].edge[(pos + 3) % 4];
			if (tmedge[cp[i].uved[0]].act == 0)
			{
				int edtmp(tmedge[cp[i].uved[0]].chd[0]);
				if (tmedge[edtmp].pt[0] == i || tmedge[edtmp].pt[1] == i)
				{
					cp[i].uved[0] = tmedge[cp[i].uved[0]].chd[0];
				}
				else
				{
					cp[i].uved[0] = tmedge[cp[i].uved[0]].chd[1];
				}
			}
			if (tmedge[cp[i].uved[1]].act == 0)
			{
				int edtmp(tmedge[cp[i].uved[1]].chd[0]);
				if (tmedge[edtmp].pt[0] == i || tmedge[edtmp].pt[1] == i)
				{
					cp[i].uved[1] = tmedge[cp[i].uved[1]].chd[0];
				}
				else
				{
					cp[i].uved[1] = tmedge[cp[i].uved[1]].chd[1];
				}
			}
			vector<int>::iterator it1 = find(cp[i].edge.begin(), cp[i].edge.end(), cp[i].uved[0]);
			vector<int>::iterator it2 = find(cp[i].edge.begin(), cp[i].edge.end(), cp[i].uved[1]);
			if (it1 == cp[i].edge.end() || it2 == cp[i].edge.end())
			{
				cerr << "Cannot find correct uv edges!\n";
				getchar();
			}
			ShootRay(i, cp[i].uved[0], cp[i].kitvUtmp);
			ShootRay(i, cp[i].uved[1], cp[i].kitvVtmp);
		}
		else
		{
			for (int j = 0; j < 4; j++)
			{
				cp[i].kitvUtmp[j] = tmedge[cp[i].edge[0]].len;
				cp[i].kitvVtmp[j] = tmedge[cp[i].edge[0]].len;
			}
		}
	}
}

void TruncatedTspline::UpdateKnotInterval_1()
{
	for (uint i = 0; i < cp.size(); i++)
	{
		//if(cp[i].trun==0)
		{
			for (int j = 0; j < 4; j++)
			{
				cp[i].kitvU[j] = cp[i].kitvUtmp[j];
				cp[i].kitvV[j] = cp[i].kitvVtmp[j];
			}
		}
	}
}

void TruncatedTspline::ShootRay(int pid, int edid, double kv[4])
{
	int loc0(0), loc1(1);
	//positive direction
	if (pid != tmedge[edid].pt[0])
	{
		loc0 = 1;
		loc1 = 0;
	}
	kv[2] = tmedge[edid].len;
	if (tmedge[edid].pn[loc1][0] == 0) //next is edge
	{
		kv[3] = tmedge[tmedge[edid].pn[loc1][1]].len;
	}
	else if (tmedge[edid].pn[loc1][0] == 1) //next is face
	{
		int fid(tmedge[edid].pn[loc1][1]), pos(0);
		for (int i = 0; i < 4; i++)
		{
			if (tmedge[tmesh[fid].edge[i]].act == 0 && tmedge[tmesh[fid].edge[i]].midpt == tmedge[edid].pt[loc1])
			{
				pos = i;
				break;
			}
		}
		kv[3] = tmedge[tmesh[fid].edge[(pos + 1) % 4]].len;
	}
	else if (tmedge[edid].pn[loc1][0] == 2) //next is XP
	{
		kv[3] = kv[2];
	}
	else if (tmedge[edid].pn[loc1][0] == 3) //end
	{
		kv[3] = 0.;
	}
	//negative direction
	if (tmedge[edid].pn[loc0][0] == 0) //previous is edge
	{
		int ed0 = tmedge[edid].pn[loc0][1];
		kv[1] = tmedge[ed0].len;
		int a0 = 0, a1 = 1;
		if (tmedge[ed0].pt[0] != pid)
		{
			a0 = 1;
			a1 = 0;
		}
		if (tmedge[ed0].pn[a1][0] == 0)
		{
			kv[0] = tmedge[tmedge[ed0].pn[a1][1]].len;
		}
		else if (tmedge[ed0].pn[a1][0] == 1)
		{
			int pt0(tmedge[ed0].pt[a1]), fid(tmedge[ed0].pn[a1][1]), pos(0);
			for (int i = 0; i < 4; i++)
			{
				if (tmedge[tmesh[fid].edge[i]].act == 0 && tmedge[tmesh[fid].edge[i]].midpt == pt0)
				{
					pos = i;
					break;
				}
			}
			kv[0] = tmedge[tmesh[fid].edge[(pos + 1) % 4]].len;
		}
		else if (tmedge[ed0].pn[a1][0] == 2)
		{
			kv[0] = kv[1];
		}
		else if (tmedge[ed0].pn[a1][0] == 3)
		{
			kv[0] = 0.;
		}
	}
	else if (tmedge[edid].pn[loc0][0] == 1)
	{
		int fid0(tmedge[edid].pn[loc0][1]), pos(0);
		for (int i = 0; i < 4; i++)
		{
			if (tmedge[tmesh[fid0].edge[i]].act == 0 && tmedge[tmesh[fid0].edge[i]].midpt == pid)
			{
				pos = i;
				break;
			}
		}
		kv[1] = tmedge[tmesh[fid0].edge[(pos + 1) % 4]].len;
		int ed0(tmesh[fid0].edge[(pos + 2) % 4]);
		if (tmedge[ed0].act == 1)
		{
			if (tmedge[ed0].face.size() == 2)
			{
				int fid1(tmedge[ed0].face[0]);
				if (fid1 == fid0)
					fid1 = tmedge[ed0].face[1];
				int *it = find(tmesh[fid1].edge, tmesh[fid1].edge + 4, ed0);
				int pos1(it - tmesh[fid1].edge);
				kv[0] = tmedge[tmesh[fid1].edge[(pos1 + 1) % 4]].len;
			}
			else
			{
				kv[0] = 0.;
			}
		}
		else
		{
			int pt0(tmedge[ed0].midpt), ed1;
			for (uint i = 0; i < cp[pt0].edge.size(); i++)
			{
				if (cp[pt0].edge[i] != tmedge[ed0].chd[0] && cp[pt0].edge[i] != tmedge[ed0].chd[1])
				{
					ed1 = cp[pt0].edge[i];
					break;
				}
			}
			kv[0] = tmedge[ed1].len;
		}
	}
	else if (tmedge[edid].pn[loc0][0] == 2)
	{
		kv[1] = kv[2];
		kv[0] = kv[2];
	}
	else
	{
		kv[1] = 0.;
		kv[0] = 0.;
	}
}

void TruncatedTspline::FindIEN_1()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		tmesh[eid].IEN.clear();
		tmesh[eid].patch_ku.clear();
		tmesh[eid].patch_kv.clear();
		if (tmesh[eid].act == 1 && (tmesh[eid].type == 0 || tmesh[eid].type == 1)) //find two ring neighorhood
		{
			array<double, 2> urang = {0., tmedge[tmesh[eid].edge[0]].len};
			array<double, 2> vrang = {0., tmedge[tmesh[eid].edge[3]].len};
			int count(0);
			//initial
			vector<int> pr0(tmesh[eid].node), er0(1, eid), pr1, er1, pr1_pref, pr1_eref;
			vector<int> rot_ref(tmesh[eid].node.size());
			vector<array<double, 2>> uv_ref(tmesh[eid].node.size());
			for (uint i = 0; i < tmesh[eid].node.size(); i++)
			{
				rot_ref[i] = tmesh[eid].lcs[i].rot;
				uv_ref[i][0] = tmesh[eid].lcs[i].u[0];
				uv_ref[i][1] = tmesh[eid].lcs[i].u[1];
			}
			for (uint i = 0; i < tmesh[eid].node.size(); i++)
			{
				array<double, 5> kui, kvi;
				FindLocalKnotVector_1(tmesh[eid].node[i], rot_ref[i], uv_ref[i], kui, kvi);
				if (CheckSupport(urang, vrang, kui, kvi))
				{
					tmesh[eid].IEN.push_back(tmesh[eid].node[i]);
					tmesh[eid].patch_ku.push_back(kui);
					tmesh[eid].patch_kv.push_back(kvi);
				}
			}
			while (count < 2)
			{
				FindNextRing(pr0, er0, pr1, er1, pr1_pref, pr1_eref);
				vector<int> rot_tmp(pr1.size());
				vector<array<double, 2>> uv_tmp(pr1.size());
				for (uint i = 0; i < pr1.size(); i++)
				{
					array<double, 5> kui, kvi;
					FindRotateAndUVCoor(pr0[pr1_pref[i]], rot_ref[pr1_pref[i]], uv_ref[pr1_pref[i]], pr1_eref[i], pr1[i], rot_tmp[i], uv_tmp[i]);
					FindLocalKnotVector_1(pr1[i], rot_tmp[i], uv_tmp[i], kui, kvi);
					if (CheckSupport(urang, vrang, kui, kvi))
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
				pr0 = pr1;
				er0 = er1;
				rot_ref = rot_tmp;
				uv_ref = uv_tmp;
				count++;
			}
		}
		else if (tmesh[eid].act == 1 && tmesh[eid].type == 4)
		{
			tmesh[eid].IEN.push_back(tmesh[eid].cnct[0]);
			int fc_pre(tmedge[tmesh[eid].edge[3]].face[0]);
			if (fc_pre == eid)
				fc_pre = tmedge[tmesh[eid].edge[3]].face[1];
			tmesh[eid].IEN.push_back(tmesh[fc_pre].cnct[3]);
			tmesh[eid].IEN.push_back(tmesh[fc_pre].cnct[2]);
			tmesh[eid].IEN.push_back(tmesh[eid].cnct[3]);
			tmesh[eid].IEN.push_back(tmesh[eid].cnct[2]);
			int fc_next(tmedge[tmesh[eid].edge[0]].face[0]);
			if (fc_next == eid)
				fc_next = tmedge[tmesh[eid].edge[0]].face[1];
			int fc_next0 = fc_next;
			while (fc_next != fc_pre)
			{
				tmesh[eid].IEN.push_back(tmesh[fc_next].cnct[3]);
				tmesh[eid].IEN.push_back(tmesh[fc_next].cnct[2]);
				int fc_nn(tmedge[tmesh[fc_next].edge[0]].face[0]);
				if (fc_nn == fc_next)
					fc_nn = tmedge[tmesh[fc_next].edge[0]].face[1];
				fc_next = fc_nn;
			}
			for (int j = 1; j < 4; j++)
			{
				for (uint k = 0; k < cp[tmesh[eid].cnct[j]].face.size(); k++)
				{
					int fcid(cp[tmesh[eid].cnct[j]].face[k]);
					if (fcid != eid && fcid != fc_pre && fcid != fc_next0)
					{
						for (uint k1 = 0; k1 < tmesh[fcid].node.size(); k1++)
						{
							vector<int>::iterator it = find(tmesh[eid].IEN.begin(), tmesh[eid].IEN.end(), tmesh[fcid].node[k1]);
							if (it == tmesh[eid].IEN.end())
							{
								array<double, 2> uv_ref = {tmesh[eid].lcs[j].u[0], tmesh[eid].lcs[j].u[1]}, uv;
								int rot;
								array<double, 5> kui, kvi;
								FindRotateAndUVCoor(tmesh[eid].cnct[j], tmesh[eid].lcs[j].rot, uv_ref, fcid, tmesh[fcid].node[k1], rot, uv);
								FindLocalKnotVector_1(tmesh[fcid].node[k1], rot, uv, kui, kvi);
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

void TruncatedTspline::FindIEN_2()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		tmesh[eid].IENtmp.clear();
		tmesh[eid].patch_kutmp.clear();
		tmesh[eid].patch_kvtmp.clear();
		if (tmesh[eid].act == 1 && (tmesh[eid].type == 0 || tmesh[eid].type == 1)) //find two ring neighorhood
		{
			array<double, 2> urang = {0., tmedge[tmesh[eid].edge[0]].len};
			array<double, 2> vrang = {0., tmedge[tmesh[eid].edge[3]].len};
			int count(0);
			//initial
			vector<int> pr0(tmesh[eid].node), er0(1, eid), pr1, er1, pr1_pref, pr1_eref;
			vector<int> rot_ref(tmesh[eid].node.size());
			vector<array<double, 2>> uv_ref(tmesh[eid].node.size());
			for (uint i = 0; i < tmesh[eid].node.size(); i++)
			{
				rot_ref[i] = tmesh[eid].lcs[i].rot;
				uv_ref[i][0] = tmesh[eid].lcs[i].u[0];
				uv_ref[i][1] = tmesh[eid].lcs[i].u[1];
			}
			for (uint i = 0; i < tmesh[eid].node.size(); i++)
			{
				array<double, 5> kui, kvi;
				FindLocalKnotVector_1(tmesh[eid].node[i], rot_ref[i], uv_ref[i], kui, kvi);
				if (CheckSupport(urang, vrang, kui, kvi))
				{
					//if(tmesh[eid].node[i]>=npt_old) tmesh[eid].aff=1;
					tmesh[eid].IENtmp.push_back(tmesh[eid].node[i]);
					tmesh[eid].patch_kutmp.push_back(kui);
					tmesh[eid].patch_kvtmp.push_back(kvi);
				}
			}
			while (count < 2)
			{
				FindNextRing(pr0, er0, pr1, er1, pr1_pref, pr1_eref);
				vector<int> rot_tmp(pr1.size());
				vector<array<double, 2>> uv_tmp(pr1.size());
				for (uint i = 0; i < pr1.size(); i++)
				{
					array<double, 5> kui, kvi;
					FindRotateAndUVCoor(pr0[pr1_pref[i]], rot_ref[pr1_pref[i]], uv_ref[pr1_pref[i]], pr1_eref[i], pr1[i], rot_tmp[i], uv_tmp[i]);
					FindLocalKnotVector_1(pr1[i], rot_tmp[i], uv_tmp[i], kui, kvi);
					if (CheckSupport(urang, vrang, kui, kvi))
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
				pr0 = pr1;
				er0 = er1;
				rot_ref = rot_tmp;
				uv_ref = uv_tmp;
				count++;
			}
		}
		else if (tmesh[eid].act == 1 && tmesh[eid].type == 4)
		{
			tmesh[eid].IENtmp.push_back(tmesh[eid].cnct[0]);
			int fc_pre(tmedge[tmesh[eid].edge[3]].face[0]);
			if (fc_pre == eid)
				fc_pre = tmedge[tmesh[eid].edge[3]].face[1];
			tmesh[eid].IENtmp.push_back(tmesh[fc_pre].cnct[3]);
			tmesh[eid].IENtmp.push_back(tmesh[fc_pre].cnct[2]);
			tmesh[eid].IENtmp.push_back(tmesh[eid].cnct[3]);
			tmesh[eid].IENtmp.push_back(tmesh[eid].cnct[2]);
			int fc_next(tmedge[tmesh[eid].edge[0]].face[0]);
			if (fc_next == eid)
				fc_next = tmedge[tmesh[eid].edge[0]].face[1];
			int fc_next0 = fc_next;
			while (fc_next != fc_pre)
			{
				tmesh[eid].IENtmp.push_back(tmesh[fc_next].cnct[3]);
				tmesh[eid].IENtmp.push_back(tmesh[fc_next].cnct[2]);
				int fc_nn(tmedge[tmesh[fc_next].edge[0]].face[0]);
				if (fc_nn == fc_next)
					fc_nn = tmedge[tmesh[fc_next].edge[0]].face[1];
				fc_next = fc_nn;
			}
			for (int j = 1; j < 4; j++)
			{
				for (uint k = 0; k < cp[tmesh[eid].cnct[j]].face.size(); k++)
				{
					int fcid(cp[tmesh[eid].cnct[j]].face[k]);
					if (fcid != eid && fcid != fc_pre && fcid != fc_next0)
					{
						for (uint k1 = 0; k1 < tmesh[fcid].node.size(); k1++)
						{
							vector<int>::iterator it = find(tmesh[eid].IENtmp.begin(), tmesh[eid].IENtmp.end(), tmesh[fcid].node[k1]);
							if (it == tmesh[eid].IENtmp.end())
							{
								array<double, 2> uv_ref = {tmesh[eid].lcs[j].u[0], tmesh[eid].lcs[j].u[1]}, uv;
								int rot;
								array<double, 5> kui, kvi;
								FindRotateAndUVCoor(tmesh[eid].cnct[j], tmesh[eid].lcs[j].rot, uv_ref, fcid, tmesh[fcid].node[k1], rot, uv);
								FindLocalKnotVector_1(tmesh[fcid].node[k1], rot, uv, kui, kvi);
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

void TruncatedTspline::FindIEN_3()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		tmesh[eid].IENtmp.clear();
		tmesh[eid].patch_kutmp.clear();
		tmesh[eid].patch_kvtmp.clear();
		if (tmesh[eid].act == 1 && (tmesh[eid].type == 0 || tmesh[eid].type == 1)) //find two ring neighorhood
		{
			array<double, 2> urang = {0., tmedge[tmesh[eid].edge[0]].len};
			array<double, 2> vrang = {0., tmedge[tmesh[eid].edge[3]].len};
			int count(0);
			//initial
			vector<int> pr0(tmesh[eid].node), er0(1, eid), pr1, er1, pr1_pref, pr1_eref;
			vector<int> rot_ref(tmesh[eid].node.size());
			vector<array<double, 2>> uv_ref(tmesh[eid].node.size());
			for (uint i = 0; i < tmesh[eid].node.size(); i++)
			{
				rot_ref[i] = tmesh[eid].lcs[i].rot;
				uv_ref[i][0] = tmesh[eid].lcs[i].u[0];
				uv_ref[i][1] = tmesh[eid].lcs[i].u[1];
			}
			for (uint i = 0; i < tmesh[eid].node.size(); i++)
			{
				array<double, 5> kui, kvi;
				FindLocalKnotVector_1(tmesh[eid].node[i], rot_ref[i], uv_ref[i], kui, kvi);
				if (CheckSupport(urang, vrang, kui, kvi))
				{
					//if(tmesh[eid].node[i]>=npt_old) tmesh[eid].aff=1;
					tmesh[eid].IENtmp.push_back(tmesh[eid].node[i]);
					tmesh[eid].patch_kutmp.push_back(kui);
					tmesh[eid].patch_kvtmp.push_back(kvi);
				}
			}
			while (count < 2)
			{
				FindNextRing(pr0, er0, pr1, er1, pr1_pref, pr1_eref);
				vector<int> rot_tmp(pr1.size());
				vector<array<double, 2>> uv_tmp(pr1.size());
				for (uint i = 0; i < pr1.size(); i++)
				{
					array<double, 5> kui, kvi;
					FindRotateAndUVCoor(pr0[pr1_pref[i]], rot_ref[pr1_pref[i]], uv_ref[pr1_pref[i]], pr1_eref[i], pr1[i], rot_tmp[i], uv_tmp[i]);
					FindLocalKnotVector_1(pr1[i], rot_tmp[i], uv_tmp[i], kui, kvi);
					if (CheckSupport(urang, vrang, kui, kvi))
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
				pr0 = pr1;
				er0 = er1;
				rot_ref = rot_tmp;
				uv_ref = uv_tmp;
				count++;
			}
		}
		else if (tmesh[eid].act == 1 && tmesh[eid].type == 4)
		{
			tmesh[eid].IENtmp.push_back(tmesh[eid].cnct[0]);
			int fc_pre(tmedge[tmesh[eid].edge[3]].face[0]);
			if (fc_pre == eid)
				fc_pre = tmedge[tmesh[eid].edge[3]].face[1];
			int pt_pre[2] = {tmesh[fc_pre].cnct[3], tmesh[fc_pre].cnct[2]};
			if (tmesh[fc_pre].type == 5)
			{
				int *it1 = find(tmesh[fc_pre].cnct, tmesh[fc_pre].cnct + 4, tmesh[eid].cnct[0]);
				int loc1(it1 - tmesh[fc_pre].cnct);
				pt_pre[0] = tmesh[fc_pre].cnct[(loc1 + 3) % 4];
				pt_pre[1] = tmesh[fc_pre].cnct[(loc1 + 2) % 4];
			}
			tmesh[eid].IENtmp.push_back(pt_pre[0]);
			tmesh[eid].IENtmp.push_back(pt_pre[1]);
			tmesh[eid].IENtmp.push_back(tmesh[eid].cnct[3]);
			tmesh[eid].IENtmp.push_back(tmesh[eid].cnct[2]);
			int fc_next(tmedge[tmesh[eid].edge[0]].face[0]);
			if (fc_next == eid)
				fc_next = tmedge[tmesh[eid].edge[0]].face[1];
			int fc_next0 = fc_next;
			int count(0);
			while (fc_next != fc_pre)
			{
				if (fc_next >= tmesh.size())
				{
					//cout<<eid<<" "<<fc_next<<" "<<tmesh.size()<<"\n";
					//cout<<tmesh[eid].cnct[0]<<" "<<tmesh[eid].cnct[1]<<" "<<tmesh[eid].cnct[2]<<" "<<tmesh[eid].cnct[3]<<"\n";
					cout << tmesh[eid].type << "\n";
					cout << cp[tmesh[eid].cnct[0]].face.size() << "\n";
					cout << cp[tmesh[eid].cnct[1]].face.size() << "\n";
					cout << cp[tmesh[eid].cnct[2]].face.size() << "\n";
					cout << cp[tmesh[eid].cnct[3]].face.size() << "\n";
					getchar();
				}
				int pt_next[2] = {tmesh[fc_next].cnct[3], tmesh[fc_next].cnct[2]};
				int loc2(0);
				if (tmesh[fc_next].type == 5)
				{
					int *it2 = find(tmesh[fc_next].cnct, tmesh[fc_next].cnct + 4, tmesh[eid].cnct[0]);
					loc2 = it2 - tmesh[fc_next].cnct;
					pt_next[0] = tmesh[fc_next].cnct[(loc2 + 3) % 4];
					pt_next[1] = tmesh[fc_next].cnct[(loc2 + 2) % 4];
				}
				tmesh[eid].IENtmp.push_back(pt_next[0]);
				tmesh[eid].IENtmp.push_back(pt_next[1]);
				int fc_nn(tmedge[tmesh[fc_next].edge[loc2]].face[0]);
				if (fc_nn == fc_next)
					fc_nn = tmedge[tmesh[fc_next].edge[loc2]].face[1];
				fc_next = fc_nn;
				count++;
				if (count > 20)
				{
					cout << "eid:" << eid << "\n";
					cerr << "Loop more than 20 times!\n";
					getchar();
					break;
				}
			}
			for (int j = 1; j < 4; j++)
			{
				for (uint k = 0; k < cp[tmesh[eid].cnct[j]].face.size(); k++)
				{
					int fcid(cp[tmesh[eid].cnct[j]].face[k]);
					if (fcid != eid && fcid != fc_pre && fcid != fc_next0)
					{
						for (uint k1 = 0; k1 < tmesh[fcid].node.size(); k1++)
						{
							vector<int>::iterator it = find(tmesh[eid].IENtmp.begin(), tmesh[eid].IENtmp.end(), tmesh[fcid].node[k1]);
							if (it == tmesh[eid].IENtmp.end())
							{
								array<double, 2> uv_ref = {tmesh[eid].lcs[j].u[0], tmesh[eid].lcs[j].u[1]}, uv;
								int rot;
								array<double, 5> kui, kvi;
								FindRotateAndUVCoor(tmesh[eid].cnct[j], tmesh[eid].lcs[j].rot, uv_ref, fcid, tmesh[fcid].node[k1], rot, uv);
								FindLocalKnotVector_1(tmesh[fcid].node[k1], rot, uv, kui, kvi);
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

void TruncatedTspline::FindIEN_4()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		tmesh[eid].IENtmp.clear();
		tmesh[eid].patch_kutmp.clear();
		tmesh[eid].patch_kvtmp.clear();
		if (tmesh[eid].act == 1 && (tmesh[eid].type == 0 || tmesh[eid].type == 1)) //find two ring neighorhood
		{
			array<double, 2> urang = {0., tmedge[tmesh[eid].edge[0]].len};
			array<double, 2> vrang = {0., tmedge[tmesh[eid].edge[3]].len};
			int count(0);
			//initial
			vector<int> pr0(tmesh[eid].node), er0(1, eid), pr1, er1, pr1_pref, pr1_eref;
			vector<int> rot_ref(tmesh[eid].node.size());
			vector<array<double, 2>> uv_ref(tmesh[eid].node.size());
			for (uint i = 0; i < tmesh[eid].node.size(); i++)
			{
				rot_ref[i] = tmesh[eid].lcs[i].rot;
				uv_ref[i][0] = tmesh[eid].lcs[i].u[0];
				uv_ref[i][1] = tmesh[eid].lcs[i].u[1];
			}
			for (uint i = 0; i < tmesh[eid].node.size(); i++)
			{
				array<double, 5> kui, kvi;
				FindLocalKnotVector_1(tmesh[eid].node[i], rot_ref[i], uv_ref[i], kui, kvi);
				if (CheckSupport(urang, vrang, kui, kvi))
				{
					//if(tmesh[eid].node[i]>=npt_old) tmesh[eid].aff=1;
					tmesh[eid].IENtmp.push_back(tmesh[eid].node[i]);
					tmesh[eid].patch_kutmp.push_back(kui);
					tmesh[eid].patch_kvtmp.push_back(kvi);
				}
			}
			while (count < 2)
			{
				FindNextRing(pr0, er0, pr1, er1, pr1_pref, pr1_eref);
				vector<int> rot_tmp(pr1.size());
				vector<array<double, 2>> uv_tmp(pr1.size());
				for (uint i = 0; i < pr1.size(); i++)
				{
					array<double, 5> kui, kvi;
					FindRotateAndUVCoor(pr0[pr1_pref[i]], rot_ref[pr1_pref[i]], uv_ref[pr1_pref[i]], pr1_eref[i], pr1[i], rot_tmp[i], uv_tmp[i]);
					FindLocalKnotVector_1(pr1[i], rot_tmp[i], uv_tmp[i], kui, kvi);
					if (CheckSupport(urang, vrang, kui, kvi))
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
				pr0 = pr1;
				er0 = er1;
				rot_ref = rot_tmp;
				uv_ref = uv_tmp;
				count++;
			}
		}
		else if (tmesh[eid].act == 1 && tmesh[eid].type == 4)
		{
			tmesh[eid].IENtmp.push_back(tmesh[eid].cnct[0]);
			int fc_pre(tmedge[tmesh[eid].edge[3]].face[0]);
			if (fc_pre == eid)
				fc_pre = tmedge[tmesh[eid].edge[3]].face[1];
			int pt_pre[2] = {tmesh[fc_pre].cnct[3], tmesh[fc_pre].cnct[2]};
			if (tmesh[fc_pre].type == 5)
			{
				int *it1 = find(tmesh[fc_pre].cnct, tmesh[fc_pre].cnct + 4, tmesh[eid].cnct[0]);
				int loc1(it1 - tmesh[fc_pre].cnct);
				pt_pre[0] = tmesh[fc_pre].cnct[(loc1 + 3) % 4];
				pt_pre[1] = tmesh[fc_pre].cnct[(loc1 + 2) % 4];
			}
			tmesh[eid].IENtmp.push_back(pt_pre[0]);
			tmesh[eid].IENtmp.push_back(pt_pre[1]);
			tmesh[eid].IENtmp.push_back(tmesh[eid].cnct[3]);
			tmesh[eid].IENtmp.push_back(tmesh[eid].cnct[2]);
			int fc_next(tmedge[tmesh[eid].edge[0]].face[0]);
			if (fc_next == eid)
				fc_next = tmedge[tmesh[eid].edge[0]].face[1];
			int fc_next0 = fc_next;
			int count(0);
			while (fc_next != fc_pre)
			{
				if (fc_next >= tmesh.size())
				{
					//cout<<eid<<" "<<fc_next<<" "<<tmesh.size()<<"\n";
					//cout<<tmesh[eid].cnct[0]<<" "<<tmesh[eid].cnct[1]<<" "<<tmesh[eid].cnct[2]<<" "<<tmesh[eid].cnct[3]<<"\n";
					cout << tmesh[eid].type << "\n";
					cout << cp[tmesh[eid].cnct[0]].face.size() << "\n";
					cout << cp[tmesh[eid].cnct[1]].face.size() << "\n";
					cout << cp[tmesh[eid].cnct[2]].face.size() << "\n";
					cout << cp[tmesh[eid].cnct[3]].face.size() << "\n";
					getchar();
				}
				int pt_next[2] = {tmesh[fc_next].cnct[3], tmesh[fc_next].cnct[2]};
				int loc2(0);
				if (tmesh[fc_next].type == 5)
				{
					int *it2 = find(tmesh[fc_next].cnct, tmesh[fc_next].cnct + 4, tmesh[eid].cnct[0]);
					loc2 = it2 - tmesh[fc_next].cnct;
					pt_next[0] = tmesh[fc_next].cnct[(loc2 + 3) % 4];
					pt_next[1] = tmesh[fc_next].cnct[(loc2 + 2) % 4];
				}
				tmesh[eid].IENtmp.push_back(pt_next[0]);
				tmesh[eid].IENtmp.push_back(pt_next[1]);
				int fc_nn(tmedge[tmesh[fc_next].edge[loc2]].face[0]);
				if (fc_nn == fc_next)
					fc_nn = tmedge[tmesh[fc_next].edge[loc2]].face[1];
				fc_next = fc_nn;
				count++;
				//if (eid == 101)
				//{
				//	cout << tmesh[eid].cnct[0] << " " << tmesh[eid].cnct[1] << " " << tmesh[eid].cnct[2] << " " << tmesh[eid].cnct[3] << "\n";
				//	cout << tmesh[eid].edge[0] << " " << tmesh[eid].edge[1] << " " << tmesh[eid].edge[2] << " " << tmesh[eid].edge[3] << "\n";
				//	int eid1(5624);
				//	cout << tmesh[eid1].cnct[0] << " " << tmesh[eid1].cnct[1] << " " << tmesh[eid1].cnct[2] << " " << tmesh[eid1].cnct[3] << "\n";
				//	cout << tmesh[eid1].edge[0] << " " << tmesh[eid1].edge[1] << " " << tmesh[eid1].edge[2] << " " << tmesh[eid1].edge[3] << "\n";
				//	//cout << cp[134].face.size() << "\n";
				//	//cout << fc_next << "\n";
				//	getchar();
				//}
				if (count > 20)
				{
					cout << "eid: " << eid << "\n";
					cerr << "Loop more than 20 times!\n";
					getchar();
					break;
				}
			}
			for (int j = 0; j < 7; j++)
				tmesh[eid].IENtmp.push_back(-1);
			int nv(cp[tmesh[eid].cnct[0]].face.size());
			for (int j = 1; j < 4; j++) //has to be regular
			{
				for (uint k = 0; k < cp[tmesh[eid].cnct[j]].face.size(); k++)
				{
					int fcid(cp[tmesh[eid].cnct[j]].face[k]);
					if (fcid == eid)
						continue;
					int nshare(0);
					for (int k1 = 0; k1 < 4; k1++)
					{
						if (tmesh[fcid].cnct[k1] == tmesh[eid].cnct[0] || tmesh[fcid].cnct[k1] == tmesh[eid].cnct[1] ||
							tmesh[fcid].cnct[k1] == tmesh[eid].cnct[2] || tmesh[fcid].cnct[k1] == tmesh[eid].cnct[3])
						{
							nshare++;
						}
					}
					if (nshare == 1)
					{
						int *it = find(tmesh[fcid].cnct, tmesh[fcid].cnct + 4, tmesh[eid].cnct[j]);
						int loc3 = it - tmesh[fcid].cnct;
						if (j == 1)
						{
							tmesh[eid].IENtmp[2 * nv + 4] = tmesh[fcid].cnct[(loc3 + 2) % 4];
							tmesh[eid].IENtmp[2 * nv + 3] = tmesh[fcid].cnct[(loc3 + 3) % 4];
						}
						else if (j == 2)
						{
							tmesh[eid].IENtmp[2 * nv + 2] = tmesh[fcid].cnct[(loc3 + 1) % 4];
							tmesh[eid].IENtmp[2 * nv + 1] = tmesh[fcid].cnct[(loc3 + 2) % 4];
							tmesh[eid].IENtmp[2 * nv + 5] = tmesh[fcid].cnct[(loc3 + 3) % 4];
						}
						else
						{
							tmesh[eid].IENtmp[2 * nv + 6] = tmesh[fcid].cnct[(loc3 + 1) % 4];
							tmesh[eid].IENtmp[2 * nv + 7] = tmesh[fcid].cnct[(loc3 + 2) % 4];
						}
					}
				}
				/*for (uint k = 0; k<cp[tmesh[eid].cnct[j]].face.size(); k++)
				{
					int fcid(cp[tmesh[eid].cnct[j]].face[k]);
					if (fcid != eid && fcid != fc_pre && fcid != fc_next0)
					{
						for (uint k1 = 0; k1<tmesh[fcid].node.size(); k1++)
						{
							vector<int>::iterator it = find(tmesh[eid].IENtmp.begin(), tmesh[eid].IENtmp.end(), tmesh[fcid].node[k1]);
							if (it == tmesh[eid].IENtmp.end())
							{
								array<double, 2> uv_ref = { tmesh[eid].lcs[j].u[0],tmesh[eid].lcs[j].u[1] }, uv;
								int rot;
								array<double, 5> kui, kvi;
								FindRotateAndUVCoor(tmesh[eid].cnct[j], tmesh[eid].lcs[j].rot, uv_ref, fcid, tmesh[fcid].node[k1], rot, uv);
								FindLocalKnotVector_1(tmesh[fcid].node[k1], rot, uv, kui, kvi);
								tmesh[eid].IENtmp.push_back(tmesh[fcid].node[k1]);
								tmesh[eid].patch_kutmp.push_back(kui);
								tmesh[eid].patch_kvtmp.push_back(kvi);
							}
						}
					}
				}*/
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

void TruncatedTspline::FindNextRing(const vector<int> &pr0, const vector<int> &er0, vector<int> &pr1, vector<int> &er1, vector<int> &pr1_pref, vector<int> &pr1_eref)
{
	pr1.clear();
	er1.clear();
	pr1_pref.clear();
	pr1_eref.clear();
	for (uint i = 0; i < pr0.size(); i++)
	{
		if (cp[pr0[i]].type != 2)
		{
			for (uint j = 0; j < cp[pr0[i]].face.size(); j++)
			{
				int fc(cp[pr0[i]].face[j]);
				vector<int>::const_iterator it1 = find(er0.begin(), er0.end(), fc);
				if (it1 == er0.end())
				{
					er1.push_back(fc);
					for (uint k = 0; k < tmesh[fc].node.size(); k++)
					{
						vector<int>::const_iterator it2 = find(pr0.begin(), pr0.end(), tmesh[fc].node[k]);
						vector<int>::iterator it3 = find(pr1.begin(), pr1.end(), tmesh[fc].node[k]);
						if (it2 == pr0.end() && it3 == pr1.end())
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

void TruncatedTspline::FindRotateAndUVCoor(int pref, int rot_ref, const array<double, 2> &uv_ref, int eid, int pid, int &rot, array<double, 2> &uv) //uv_ref is "global"
{
	vector<int>::iterator it0 = find(tmesh[eid].node.begin(), tmesh[eid].node.end(), pref);
	vector<int>::iterator it1 = find(tmesh[eid].node.begin(), tmesh[eid].node.end(), pid);
	int loc0 = it0 - tmesh[eid].node.begin();
	int loc1 = it1 - tmesh[eid].node.begin();
	int rot1 = (tmesh[eid].lcs[loc1].rot + 4 - tmesh[eid].lcs[loc0].rot) % 4;
	rot = (rot1 + rot_ref) % 4;
	double tmp[2] = {tmesh[eid].lcs[loc1].u[0] - tmesh[eid].lcs[loc0].u[0], tmesh[eid].lcs[loc1].u[1] - tmesh[eid].lcs[loc0].u[1]}; //direction vecton in element eid
	int rot2 = (rot_ref - tmesh[eid].lcs[loc0].rot + 4) % 4;
	double tmp1[2] = {tmp[0], tmp[1]};
	if (rot2 == 1)
	{
		tmp1[0] = -tmp[1];
		tmp1[1] = tmp[0];
	}
	else if (rot2 == 2)
	{
		tmp1[0] = -tmp[0];
		tmp1[1] = -tmp[1];
	}
	else if (rot2 == 3)
	{
		tmp1[0] = tmp[1];
		tmp1[1] = -tmp[0];
	}
	uv[0] = uv_ref[0] + tmp1[0];
	uv[1] = uv_ref[1] + tmp1[1];
}

void TruncatedTspline::FindLocalKnotVector_1(int id, int rot, const array<double, 2> &uv, array<double, 5> &ku, array<double, 5> &kv)
{
	ku[2] = uv[0];
	kv[2] = uv[1];
	if (rot == 0)
	{
		for (int i = 0; i < 2; i++)
		{
			ku[i + 3] = ku[i + 2] + cp[id].kitvU[i + 2];
			kv[i + 3] = kv[i + 2] + cp[id].kitvV[i + 2];
			ku[1 - i] = ku[2 - i] - cp[id].kitvU[1 - i];
			kv[1 - i] = kv[2 - i] - cp[id].kitvV[1 - i];
		}
	}
	else if (rot == 1)
	{
		for (int i = 0; i < 2; i++)
		{
			ku[i + 3] = ku[i + 2] + cp[id].kitvV[1 - i];
			kv[i + 3] = kv[i + 2] + cp[id].kitvU[i + 2];
			ku[1 - i] = ku[2 - i] - cp[id].kitvV[i + 2];
			kv[1 - i] = kv[2 - i] - cp[id].kitvU[1 - i];
		}
	}
	else if (rot == 2)
	{
		for (int i = 0; i < 2; i++)
		{
			ku[i + 3] = ku[i + 2] + cp[id].kitvU[1 - i];
			kv[i + 3] = kv[i + 2] + cp[id].kitvV[1 - i];
			ku[1 - i] = ku[2 - i] - cp[id].kitvU[2 + i];
			kv[1 - i] = kv[2 - i] - cp[id].kitvV[2 + i];
		}
	}
	else if (rot == 3)
	{
		for (int i = 0; i < 2; i++)
		{
			ku[i + 3] = ku[i + 2] + cp[id].kitvV[i + 2];
			kv[i + 3] = kv[i + 2] + cp[id].kitvU[1 - i];
			ku[1 - i] = ku[2 - i] - cp[id].kitvV[1 - i];
			kv[1 - i] = kv[2 - i] - cp[id].kitvU[i + 2];
		}
	}
	else
	{
		cerr << "rot cannot be other numbers!\n";
		getchar();
	}
}

bool TruncatedTspline::CheckSupport(const array<double, 2> &u, const array<double, 2> &v, const array<double, 5> &ku, const array<double, 5> &kv)
{
	if (ku[0] < u[1] && ku[4] > u[0] && kv[0] < v[1] && kv[4] > v[0])
	{
		return true;
	}
	else
	{
		return false;
	}
}

void TruncatedTspline::SetLocalCoorSystem()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		tmesh[eid].node.clear();
		tmesh[eid].lcs.clear();
		if (tmesh[eid].act == 1)
		{
			double ul[2] = {tmedge[tmesh[eid].edge[0]].len, tmedge[tmesh[eid].edge[1]].len};
			double uvcoor[8][2] = {{0., 0.}, {ul[0] / 2., 0.}, {ul[0], 0.}, {ul[0], ul[1] / 2.}, {ul[0], ul[1]}, {ul[0] / 2., ul[1]}, {0., ul[1]}, {0., ul[1] / 2.}};
			for (int i = 0; i < 4; i++)
			{
				int uved[2] = {tmesh[eid].edge[i], tmesh[eid].edge[(i + 3) % 4]};
				for (int j = 0; j < 2; j++)
				{
					if (tmedge[uved[j]].act == 0)
					{
						int edtmp(tmedge[uved[j]].chd[0]);
						if (tmedge[edtmp].pt[0] != tmesh[eid].cnct[i] && tmedge[edtmp].pt[1] != tmesh[eid].cnct[i])
						{
							edtmp = tmedge[uved[j]].chd[1];
						}
						uved[j] = edtmp;
					}
				}
				tmesh[eid].node.push_back(tmesh[eid].cnct[i]);
				ELCS tmp;
				tmp.u[0] = uvcoor[2 * i][0];
				tmp.u[1] = uvcoor[2 * i][1];
				if (cp[tmesh[eid].cnct[i]].uved[0] == uved[0])
				{
					tmp.rot = i;
				}
				else if (cp[tmesh[eid].cnct[i]].uved[0] == uved[1])
				{
					tmp.rot = (i + 1) % 4;
				}
				else if (cp[tmesh[eid].cnct[i]].uved[1] == uved[0])
				{
					tmp.rot = (i + 3) % 4;
				}
				else
				{
					tmp.rot = (i + 2) % 4;
				}
				tmesh[eid].lcs.push_back(tmp);

				if (tmedge[tmesh[eid].edge[i]].act == 0)
				{
					int ued(tmedge[tmesh[eid].edge[i]].chd[1]);
					if (tmedge[ued].pt[0] == tmesh[eid].cnct[i] || tmedge[ued].pt[1] == tmesh[eid].cnct[i])
					{
						ued = tmedge[tmesh[eid].edge[i]].chd[0];
					}
					tmesh[eid].node.push_back(tmedge[tmesh[eid].edge[i]].midpt);
					ELCS tmp1;
					tmp1.u[0] = uvcoor[2 * i + 1][0];
					tmp1.u[1] = uvcoor[2 * i + 1][1];
					if (cp[tmedge[tmesh[eid].edge[i]].midpt].uved[1] == ued)
					{
						tmp1.rot = (i + 3) % 4;
					}
					else
					{
						tmp1.rot = (i + 2) % 4;
					}
					tmesh[eid].lcs.push_back(tmp1);
				}
			}
		}
	}
}

void TruncatedTspline::UpdatePatchCP_Unstruct(int eid)
{
	if (tmesh[eid].act == 0 && (tmesh[eid].type == 0 || tmesh[eid].type == 1 || tmesh[eid].type == 2))
	{
		int chdid(tmesh[eid].chd[0]);
		//vector<int> upid;
		vector<vector<double>> cmat(tmesh[chdid].IENtmp.size(), vector<double>(tmesh[eid].IEN.size(), 0.)), tmat(tmesh[chdid].IENtmp.size(), vector<double>(tmesh[eid].IEN.size(), 0.));
		for (uint i = 0; i < tmesh[chdid].IENtmp.size(); i++)
		{
			//if(tmesh[chdid].IENtmp[i] >= npt_old)
			{
				//upid.push_back(i);
				//int pos(tmesh[chdid].IENtmp[i]);
				//vector<double> cvec(tmesh[eid].IEN.size(),0.);
				for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
				{
					if (tmesh[chdid].patch_kutmp[i] != tmesh[eid].patch_ku[j] || tmesh[chdid].patch_kvtmp[i] != tmesh[eid].patch_kv[j])
					{
						if (CheckSubKnotVector(tmesh[chdid].patch_kutmp[i], tmesh[chdid].patch_kvtmp[i], tmesh[eid].patch_ku[j], tmesh[eid].patch_kv[j]))
						{
							vector<double> ku1(10), kv1(10);
							vector<vector<double>> Tu, Tv;
							vector<double> ku(tmesh[eid].patch_ku[j].begin(), tmesh[eid].patch_ku[j].end());
							vector<double> kv(tmesh[eid].patch_kv[j].begin(), tmesh[eid].patch_kv[j].end());
							vector<double>::iterator it1, it2;
							it1 = set_union(tmesh[eid].patch_ku[j].begin(), tmesh[eid].patch_ku[j].end(), tmesh[chdid].patch_kutmp[i].begin(), tmesh[chdid].patch_kutmp[i].end(), ku1.begin());
							it2 = set_union(tmesh[eid].patch_kv[j].begin(), tmesh[eid].patch_kv[j].end(), tmesh[chdid].patch_kvtmp[i].begin(), tmesh[chdid].patch_kvtmp[i].end(), kv1.begin());
							ku1.resize(it1 - ku1.begin());
							kv1.resize(it2 - kv1.begin());
							TMatrix(ku, ku1, 3, Tu);
							TMatrix(kv, kv1, 3, Tv);
							it1 = search(ku1.begin(), ku1.end(), tmesh[chdid].patch_kutmp[i].begin(), tmesh[chdid].patch_kutmp[i].end());
							it2 = search(kv1.begin(), kv1.end(), tmesh[chdid].patch_kvtmp[i].begin(), tmesh[chdid].patch_kvtmp[i].end());
							if (it1 != ku1.end() && it2 != kv1.end())
							{
								int loc1 = it1 - ku1.begin();
								int loc2 = it2 - kv1.begin();
								//double coef=Tu[loc1][0]*Tv[loc2][0];
								//cvec[j]=Tu[loc1][0]*Tv[loc2][0];
								cmat[i][j] = Tu[loc1][0] * Tv[loc2][0];
								tmat[i][j] = cmat[i][j];
							}
						}
					}
					else
					{
						cmat[i][j] = 1.;
						tmat[i][j] = 1.;
					}
				}
				//cmat.push_back(cvec);
				//tmat.push_back(cvec);
			}
			//else
			//{
			//	vector<int>::iterator it=find(tmesh[eid].IEN.begin(),tmesh[eid].IEN.end(),tmesh[chdid].IENtmp[i]);
			//	if(it!=tmesh[eid].IEN.end())
			//	{
			//		int i0(it-tmesh[eid].IEN.begin());
			//		if(tmesh[chdid].patch_kutmp[i]!=tmesh[eid].patch_ku[i0] || tmesh[chdid].patch_kvtmp[i]!=tmesh[eid].patch_kv[i0])
			//		{
			//			//upid.push_back(i);
			//			//vector<double> cvec(tmesh[eid].IEN.size(),0.);
			//			for(uint j=0; j<tmesh[eid].IEN.size(); j++)
			//			{
			//				if(CheckSubKnotVector(tmesh[chdid].patch_kutmp[i],tmesh[chdid].patch_kvtmp[i],tmesh[eid].patch_ku[j],tmesh[eid].patch_kv[j]))
			//				{
			//					vector<double> ku1(10),kv1(10);
			//					vector<vector<double>> Tu,Tv;
			//					vector<double> ku(tmesh[eid].patch_ku[j].begin(),tmesh[eid].patch_ku[j].end());
			//					vector<double> kv(tmesh[eid].patch_kv[j].begin(),tmesh[eid].patch_kv[j].end());
			//					vector<double>::iterator it1,it2;
			//					it1=set_union(tmesh[eid].patch_ku[j].begin(),tmesh[eid].patch_ku[j].end(),tmesh[chdid].patch_kutmp[i].begin(),tmesh[chdid].patch_kutmp[i].end(),ku1.begin());
			//					it2=set_union(tmesh[eid].patch_kv[j].begin(),tmesh[eid].patch_kv[j].end(),tmesh[chdid].patch_kvtmp[i].begin(),tmesh[chdid].patch_kvtmp[i].end(),kv1.begin());
			//					ku1.resize(it1-ku1.begin());
			//					kv1.resize(it2-kv1.begin());
			//					TMatrix(ku,ku1,3,Tu);
			//					TMatrix(kv,kv1,3,Tv);
			//					it1=search(ku1.begin(),ku1.end(),tmesh[chdid].patch_kutmp[i].begin(),tmesh[chdid].patch_kutmp[i].end());
			//					it2=search(kv1.begin(),kv1.end(),tmesh[chdid].patch_kvtmp[i].begin(),tmesh[chdid].patch_kvtmp[i].end());
			//					if(it1!=ku1.end() && it2!=kv1.end())
			//					{
			//						int loc1=it1-ku1.begin();
			//						int loc2=it2-kv1.begin();
			//						//double coef=Tu[loc1][0]*Tv[loc2][0];
			//						//cvec[j]=Tu[loc1][0]*Tv[loc2][0];
			//						cmat[i][j]=Tu[loc1][0]*Tv[loc2][0];
			//						tmat[i][j]=cmat[i][j];
			//					}
			//				}
			//			}
			//			//cmat.push_back(cvec);
			//			//tmat.push_back(cvec);
			//		}
			//	}
			//}
		}
		for (uint i = 0; i < tmat.size(); i++)
		{
			int ptid(tmesh[chdid].IENtmp[i]);
			cp[ptid].coortmp[0] = 0.;
			cp[ptid].coortmp[1] = 0.;
			cp[ptid].coortmp[2] = 0.;
			for (uint j = 0; j < tmat[i].size(); j++)
			{
				if (tmat[i][j] != 0.)
				{
					for (uint k = 0; k < cp[tmesh[eid].IEN[j]].tbf.size(); k++)
					{
						vector<int>::iterator it = find(tmesh[eid].IEN.begin(), tmesh[eid].IEN.end(), cp[tmesh[eid].IEN[j]].tbf[k]);
						if (it != tmesh[eid].IEN.end())
						{
							int loc(it - tmesh[eid].IEN.begin());
							tmat[i][j] -= cmat[i][loc] * cp[tmesh[eid].IEN[j]].tc[k];
						}
					}
					cp[ptid].coortmp[0] += tmat[i][j] * cp[tmesh[eid].IEN[j]].coor[0];
					cp[ptid].coortmp[1] += tmat[i][j] * cp[tmesh[eid].IEN[j]].coor[1];
					cp[ptid].coortmp[2] += tmat[i][j] * cp[tmesh[eid].IEN[j]].coor[2];
				}
			}
		}
		//check truncation
		vector<int> pid_old(tmesh[eid].IEN);
		vector<int> pid_new(tmesh[chdid].IENtmp);
		for (uint i = 0; i < pid_new.size(); i++)
		{
			if (cp[pid_new[i]].update == 0)
			{
				if (pid_new[i] >= npt_old)
				{
					cp[pid_new[i]].trun = 0;
					cp[pid_new[i]].coor[0] = cp[pid_new[i]].coortmp[0];
					cp[pid_new[i]].coor[1] = cp[pid_new[i]].coortmp[1];
					cp[pid_new[i]].coor[2] = cp[pid_new[i]].coortmp[2];
					cp[pid_new[i]].update = 1;
				}
				else
				{
					vector<int>::iterator it = find(pid_old.begin(), pid_old.end(), pid_new[i]);
					if (it != pid_old.end())
					{
						int i0(it - pid_old.begin());
						if (tmesh[chdid].patch_kutmp[i] != tmesh[eid].patch_ku[i0] || tmesh[chdid].patch_kvtmp[i] != tmesh[eid].patch_kv[i0])
						{
							if (cp[pid_new[i]].trun == 1)
							{
								cout << "not ready yet!\n";
								getchar();
							}
							else
							{
								vector<array<double, 5>> chdu, chdv;
								vector<double> coef, tc0;
								vector<int> tbf0;
								for (uint j = 0; j < pid_new.size(); j++)
								{
									if (cmat[j][i0] != 0.)
									{
										chdu.push_back(tmesh[chdid].patch_kutmp[j]);
										chdv.push_back(tmesh[chdid].patch_kvtmp[j]);
										coef.push_back(cmat[j][i0]);
										if (j != i)
										{
											tbf0.push_back(pid_new[j]);
											tc0.push_back(cmat[j][i0]);
										}
									}
								}
								if (CheckFullChildren(tmesh[eid].patch_ku[i0], tmesh[eid].patch_kv[i0], chdu, chdv, coef))
								{
									cp[pid_new[i]].trun = 0;
									cp[pid_new[i]].update = 2;
								}
								else
								{
									cp[pid_new[i]].truntmp = 1;
									for (uint j = 0; j < tbf0.size(); j++)
									{
										cp[pid_new[i]].tbftmp.push_back(tbf0[j]);
										cp[pid_new[i]].tctmp.push_back(tc0[j]);
									}
									cp[pid_new[i]].update = 3;
								}
							}
						}
					}
				}
			}
		}
	}
	else if (tmesh[eid].act == 0 && tmesh[eid].type == 4)
	{
		//don't develop truncated basis funcions here, just Catmull-Clark subdivision
	}
	else
	{
		cerr << "Unknown element to update CP!\n";
		getchar();
	}
}

void TruncatedTspline::UpdatePatchCP_Unstruct_1(int eid)
{
	vector<int> IEN_new_all;
	vector<array<double, 5>> KU_new_all;
	vector<array<double, 5>> KV_new_all;
	for (int i = 0; i < 4; i++)
	{
		if (tmesh[eid].chd[i] != -1)
		{
			int cid(tmesh[eid].chd[i]);
			for (uint j = 0; j < tmesh[cid].IENtmp.size(); j++)
			{
				vector<int>::iterator it = find(IEN_new_all.begin(), IEN_new_all.end(), tmesh[cid].IENtmp[j]);
				if (it == IEN_new_all.end())
				{
					IEN_new_all.push_back(tmesh[cid].IENtmp[j]);
					array<double, 5> kutmp, kvtmp;
					for (int k = 0; k < 5; k++)
					{
						kutmp[k] = tmesh[cid].patch_kutmp[j][k] + tmesh[eid].chd_o[i][0];
						kvtmp[k] = tmesh[cid].patch_kvtmp[j][k] + tmesh[eid].chd_o[i][1];
					}
					KU_new_all.push_back(kutmp);
					KV_new_all.push_back(kvtmp);
				}
			}
		}
	}

	//cout<<"old:\n";
	//for(uint i=0; i<tmesh[eid].IEN.size(); i++)
	//{
	//	cout<<"pid: "<<tmesh[eid].IEN[i]<<"\n";
	//	for(int j=0; j<5; j++) cout<<tmesh[eid].patch_ku[i][j]<<" ";
	//	cout<<"\n";
	//	for(int j=0; j<5; j++) cout<<tmesh[eid].patch_kv[i][j]<<" ";
	//	cout<<"\n\n";
	//}

	//cout<<"new:\n";
	//for(uint i=0; i<IEN_new_all.size(); i++)
	//{
	//	cout<<"pid: "<<IEN_new_all[i]<<"\n";
	//	for(int j=0; j<5; j++) cout<<KU_new_all[i][j]<<" ";
	//	cout<<"\n";
	//	for(int j=0; j<5; j++) cout<<KV_new_all[i][j]<<" ";
	//	cout<<"\n\n";
	//}
	//getchar();

	vector<int> IEN_old(tmesh[eid].IEN);
	vector<vector<double>> cmat(IEN_new_all.size(), vector<double>(IEN_old.size(), 0.));
	//vector<vector<double>> tmat(IEN_new_all.size(),vector<double>(IEN_old.size(),0.));
	for (uint i = 0; i < IEN_old.size(); i++)
	{
		vector<int>::iterator it = find(IEN_new_all.begin(), IEN_new_all.end(), IEN_old[i]);
		if (it != IEN_new_all.end())
		{
			int loc(it - IEN_new_all.begin());
			if (tmesh[eid].patch_ku[i] != KU_new_all[loc] || tmesh[eid].patch_kv[i] != KV_new_all[loc])
			{
				cp[IEN_old[i]].aff = 1;
			}
			else
			{
				cmat[loc][i] = 1.;
			}
		}
	}

	for (uint i = 0; i < IEN_new_all.size(); i++)
	{
		for (uint j = 0; j < IEN_old.size(); j++)
		{
			if (CheckSubKnotVector(KU_new_all[i], KV_new_all[i], tmesh[eid].patch_ku[j], tmesh[eid].patch_kv[j]))
			{
				vector<double> ku1(10), kv1(10);
				vector<vector<double>> Tu, Tv;
				vector<double> ku(tmesh[eid].patch_ku[j].begin(), tmesh[eid].patch_ku[j].end());
				vector<double> kv(tmesh[eid].patch_kv[j].begin(), tmesh[eid].patch_kv[j].end());
				vector<double>::iterator it1, it2;
				it1 = set_union(tmesh[eid].patch_ku[j].begin(), tmesh[eid].patch_ku[j].end(), KU_new_all[i].begin(), KU_new_all[i].end(), ku1.begin());
				it2 = set_union(tmesh[eid].patch_kv[j].begin(), tmesh[eid].patch_kv[j].end(), KV_new_all[i].begin(), KV_new_all[i].end(), kv1.begin());
				ku1.resize(it1 - ku1.begin());
				kv1.resize(it2 - kv1.begin());
				TMatrix(ku, ku1, 3, Tu);
				TMatrix(kv, kv1, 3, Tv);
				it1 = search(ku1.begin(), ku1.end(), KU_new_all[i].begin(), KU_new_all[i].end());
				it2 = search(kv1.begin(), kv1.end(), KV_new_all[i].begin(), KV_new_all[i].end());
				if (it1 != ku1.end() && it2 != kv1.end())
				{
					int loc1 = it1 - ku1.begin();
					int loc2 = it2 - kv1.begin();
					cmat[i][j] = Tu[loc1][0] * Tv[loc2][0];
					//tmat[i][j]=cmat[i][j];
				}
			}
		}
	}
	for (uint i = 0; i < IEN_new_all.size(); i++)
	{
		cp[IEN_new_all[i]].coortmp[0] = 0.;
		cp[IEN_new_all[i]].coortmp[1] = 0.;
		cp[IEN_new_all[i]].coortmp[2] = 0.;
		cp[IEN_new_all[i]].wtmp = 0.;
		for (uint j = 0; j < IEN_old.size(); j++)
		{
			if (cmat[i][j] != 0.)
			{
				//for(uint k=0; k<cp[tmesh[eid].IEN[j]].tbf.size(); k++)
				//{
				//	vector<int>::iterator it=find(tmesh[eid].IEN.begin(),tmesh[eid].IEN.end(),cp[tmesh[eid].IEN[j]].tbf[k]);
				//	if(it!=tmesh[eid].IEN.end())
				//	{
				//		int loc(it-tmesh[eid].IEN.begin());
				//		tmat[i][j]-=cmat[i][loc]*cp[tmesh[eid].IEN[j]].tc[k];
				//	}
				//}
				cp[IEN_new_all[i]].coortmp[0] += cmat[i][j] * cp[IEN_old[j]].coor[0] * cp[IEN_old[j]].w;
				cp[IEN_new_all[i]].coortmp[1] += cmat[i][j] * cp[IEN_old[j]].coor[1] * cp[IEN_old[j]].w;
				cp[IEN_new_all[i]].coortmp[2] += cmat[i][j] * cp[IEN_old[j]].coor[2] * cp[IEN_old[j]].w;
				cp[IEN_new_all[i]].wtmp += cmat[i][j] * cp[IEN_old[j]].w;
				//cp[IEN_old[j]].truntmp=1;
				//vector<int>::iterator it=find(cp[IEN_old[j]].tbftmp.begin(),cp[IEN_old[j]].tbftmp.end(),IEN_new_all[i]);
				//if(it==cp[IEN_old[j]].tbftmp.end())
				//{
				//	cp[IEN_old[j]].tbftmp.push_back(IEN_new_all[i]);
				//	cp[IEN_old[j]].tctmp.push_back(cmat[i][j]);
				//}
			}
		}
		cp[IEN_new_all[i]].coortmp[0] /= cp[IEN_new_all[i]].wtmp;
		cp[IEN_new_all[i]].coortmp[1] /= cp[IEN_new_all[i]].wtmp;
		cp[IEN_new_all[i]].coortmp[2] /= cp[IEN_new_all[i]].wtmp;
	}
	//check truncation
	double sp0[4][2] = {{0., 0.}, {tmedge[tmesh[eid].edge[0]].len, 0.}, {tmedge[tmesh[eid].edge[0]].len, tmedge[tmesh[eid].edge[1]].len}, {0., tmedge[tmesh[eid].edge[1]].len}};
	vector<array<double, 2>> spt(9);
	for (int i = 0; i < 4; i++)
	{
		spt[i][0] = sp0[i][0];
		spt[i][1] = sp0[i][1];
		spt[i + 4][0] = (sp0[i][0] + sp0[(i + 1) % 4][0]) / 2.;
		spt[i + 4][1] = (sp0[i][1] + sp0[(i + 1) % 4][1]) / 2.;
	}
	spt[8][0] = spt[4][0];
	spt[8][1] = spt[5][1];
	for (uint i = 0; i < IEN_old.size(); i++)
	{
		vector<double> coef;
		for (uint j = 0; j < IEN_new_all.size(); j++)
		{
			coef.push_back(cmat[j][i]);
		}
		if (!CheckFullChildren_1(spt, tmesh[eid].patch_ku[i], tmesh[eid].patch_kv[i], KU_new_all, KV_new_all, coef))
		{
			//if(IEN_old[i]==1)
			//{
			//	cout<<"old knot vector:\n";
			//	cout<<tmesh[eid].patch_ku[i][0]<<" "<<tmesh[eid].patch_ku[i][1]<<" "<<tmesh[eid].patch_ku[i][2]<<" "<<tmesh[eid].patch_ku[i][3]<<" "<<tmesh[eid].patch_ku[i][4]<<"\n";
			//	cout<<tmesh[eid].patch_kv[i][0]<<" "<<tmesh[eid].patch_kv[i][1]<<" "<<tmesh[eid].patch_kv[i][2]<<" "<<tmesh[eid].patch_kv[i][3]<<" "<<tmesh[eid].patch_kv[i][4]<<"\n";
			//	cout<<"new knot vector:\n";
			//	for(uint j=0; j<IEN_new_all.size(); j++)
			//	{
			//		if(cmat[j][i]!=0.)
			//		{
			//			cout<<"ID: "<<IEN_new_all[j]<<"\n";
			//			cout<<KU_new_all[j][0]<<" "<<KU_new_all[j][1]<<" "<<KU_new_all[j][2]<<" "<<KU_new_all[j][3]<<" "<<KU_new_all[j][4]<<"\n";
			//			cout<<KV_new_all[j][0]<<" "<<KV_new_all[j][1]<<" "<<KV_new_all[j][2]<<" "<<KV_new_all[j][3]<<" "<<KV_new_all[j][4]<<"\n";
			//			cout<<coef[j]<<"\n";
			//		}
			//	}
			//	//double tmp=CheckFullChildren_2(spt,tmesh[eid].patch_ku[i],tmesh[eid].patch_kv[i],KU_new_all,KV_new_all,coef);
			//	//cout<<tmp<<"\n";
			//	getchar();
			//}
			cp[IEN_old[i]].truntmp = 1;
			for (uint j = 0; j < IEN_new_all.size(); j++)
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

	/*vector<int> pid_old(tmesh[eid].IEN);
	vector<int> pid_new(tmesh[chid].IENtmp);
	for(uint i=0; i<pid_new.size(); i++)
	{
		if(cp[pid_new[i]].update==0)
		{
			if(pid_new[i]>=npt_old)
			{
				cp[pid_new[i]].trun=0;
				cp[pid_new[i]].coor[0]=cp[pid_new[i]].coortmp[0];
				cp[pid_new[i]].coor[1]=cp[pid_new[i]].coortmp[1];
				cp[pid_new[i]].coor[2]=cp[pid_new[i]].coortmp[2];
				cp[pid_new[i]].update=1;
			}
			else
			{
				vector<int>::iterator it=find(pid_old.begin(),pid_old.end(),pid_new[i]);
				if(it!=pid_old.end())
				{
					int i0(it-pid_old.begin());
					if(tmesh[chdid].patch_kutmp[i]!=tmesh[eid].patch_ku[i0] || tmesh[chdid].patch_kvtmp[i]!=tmesh[eid].patch_kv[i0])
					{
						if(cp[pid_new[i]].trun==1)
						{
							cout<<"not ready yet!\n";
							getchar();
						}
						else
						{
							vector<array<double,5>> chdu,chdv;
							vector<double> coef, tc0;
							vector<int> tbf0;
							for(uint j=0; j<pid_new.size(); j++)
							{
								if(cmat[j][i0]!=0.)
								{
									chdu.push_back(tmesh[chdid].patch_kutmp[j]);
									chdv.push_back(tmesh[chdid].patch_kvtmp[j]);
									coef.push_back(cmat[j][i0]);
									if(j!=i)
									{
										tbf0.push_back(pid_new[j]);
										tc0.push_back(cmat[j][i0]);
									}
								}
							}
							if(CheckFullChildren(tmesh[eid].patch_ku[i0],tmesh[eid].patch_kv[i0],chdu,chdv,coef))
							{
								cp[pid_new[i]].trun=0;
								cp[pid_new[i]].update=2;
							}
							else
							{
								cp[pid_new[i]].truntmp=1;
								for(uint j=0; j<tbf0.size(); j++)
								{
									cp[pid_new[i]].tbftmp.push_back(tbf0[j]);
									cp[pid_new[i]].tctmp.push_back(tc0[j]);
								}
								cp[pid_new[i]].update=3;
							}
						}
					}
				}
			}
		}
	}*/
}

void TruncatedTspline::UpdatePatchCP_Unstruct_2(int eid)
{
	vector<int> IEN_new_all;
	vector<array<double, 5>> KU_new_all;
	vector<array<double, 5>> KV_new_all;
	for (int i = 0; i < 4; i++)
	{
		if (tmesh[eid].chd[i] != -1)
		{
			int cid(tmesh[eid].chd[i]);
			for (uint j = 0; j < tmesh[cid].IENtmp.size(); j++)
			{
				vector<int>::iterator it = find(IEN_new_all.begin(), IEN_new_all.end(), tmesh[cid].IENtmp[j]);
				if (it == IEN_new_all.end())
				{
					IEN_new_all.push_back(tmesh[cid].IENtmp[j]);
					array<double, 5> kutmp, kvtmp;
					for (int k = 0; k < 5; k++)
					{
						kutmp[k] = tmesh[cid].patch_kutmp[j][k] + tmesh[eid].chd_o[i][0];
						kvtmp[k] = tmesh[cid].patch_kvtmp[j][k] + tmesh[eid].chd_o[i][1];
					}
					KU_new_all.push_back(kutmp);
					KV_new_all.push_back(kvtmp);
				}
			}
		}
	}

	vector<int> IEN_old(tmesh[eid].IEN);

	for (uint i = 0; i < IEN_new_all.size(); i++)
	{
		if (cp[IEN_new_all[i]].update == 0) //not updated yet
		{
			for (uint j = 0; j < IEN_old.size(); j++)
			{
				if (CheckSubKnotVector(KU_new_all[i], KV_new_all[i], tmesh[eid].patch_ku[j], tmesh[eid].patch_kv[j]))
				{
					vector<double> ku1(10), kv1(10);
					vector<vector<double>> Tu, Tv;
					vector<double> ku(tmesh[eid].patch_ku[j].begin(), tmesh[eid].patch_ku[j].end());
					vector<double> kv(tmesh[eid].patch_kv[j].begin(), tmesh[eid].patch_kv[j].end());
					vector<double>::iterator it1, it2;
					it1 = set_union(tmesh[eid].patch_ku[j].begin(), tmesh[eid].patch_ku[j].end(), KU_new_all[i].begin(), KU_new_all[i].end(), ku1.begin());
					it2 = set_union(tmesh[eid].patch_kv[j].begin(), tmesh[eid].patch_kv[j].end(), KV_new_all[i].begin(), KV_new_all[i].end(), kv1.begin());
					ku1.resize(it1 - ku1.begin());
					kv1.resize(it2 - kv1.begin());
					TMatrix(ku, ku1, 3, Tu);
					TMatrix(kv, kv1, 3, Tv);
					it1 = search(ku1.begin(), ku1.end(), KU_new_all[i].begin(), KU_new_all[i].end());
					it2 = search(kv1.begin(), kv1.end(), KV_new_all[i].begin(), KV_new_all[i].end());
					if (it1 != ku1.end() && it2 != kv1.end())
					{
						int loc1 = it1 - ku1.begin();
						int loc2 = it2 - kv1.begin();
						double tmp = Tu[loc1][0] * Tv[loc2][0];
						cp[IEN_new_all[i]].coortmp[0] += tmp * cp[IEN_old[j]].coor[0];
						cp[IEN_new_all[i]].coortmp[1] += tmp * cp[IEN_old[j]].coor[1];
						cp[IEN_new_all[i]].coortmp[2] += tmp * cp[IEN_old[j]].coor[2];
						cp[IEN_new_all[i]].wtmp += tmp * cp[IEN_old[j]].w;
					}
				}
			}
			cp[IEN_new_all[i]].update = 1;
		}
	}
}

bool TruncatedTspline::CheckFullChildren(const array<double, 5> &motu, const array<double, 5> &motv, const vector<array<double, 5>> &chdu, const vector<array<double, 5>> &chdv, const vector<double> &coef)
{
	double tol(1.e-8);
	int nsp(5);
	double sp[5][2] = {{motu[2], motv[2]}, {motu[2], motv[1]}, {motu[3], motv[2]}, {motu[2], motv[3]}, {motu[1], motv[2]}};
	double sum(0.);
	for (int i = 0; i < nsp; i++)
	{
		vector<double> ku(motu.begin(), motu.end());
		vector<double> kv(motv.begin(), motv.end());
		vector<double> uval, vval;
		BSplineBasis bu, bv;
		bu.Set(3, ku);
		bv.Set(3, kv);
		bu.BasisFunction(0, sp[i][0], 0, uval);
		bv.BasisFunction(0, sp[i][1], 0, vval);
		double Ni = uval[0] * vval[0];
		for (uint j = 0; j < chdu.size(); j++)
		{
			ku.assign(chdu[j].begin(), chdu[j].end());
			kv.assign(chdv[j].begin(), chdv[j].end());
			bu.Set(3, ku);
			bv.Set(3, kv);
			bu.BasisFunction(0, sp[i][0], 0, uval);
			bv.BasisFunction(0, sp[i][1], 0, vval);
			Ni -= coef[j] * uval[0] * vval[0];
		}
		if (Ni < 0.)
			Ni = -Ni;
		sum += Ni;
	}
	if (sum / nsp < tol)
		return true;
	else
		return false;
}

bool TruncatedTspline::CheckFullChildren_1(const vector<array<double, 2>> &spt, const array<double, 5> &motu, const array<double, 5> &motv, const vector<array<double, 5>> &chdu, const vector<array<double, 5>> &chdv, const vector<double> &coef)
{
	double tol(1.e-8);
	double sum(0.);
	for (uint i = 0; i < spt.size(); i++)
	{
		vector<double> ku(motu.begin(), motu.end());
		vector<double> kv(motv.begin(), motv.end());
		vector<double> uval, vval;
		BSplineBasis bu, bv;
		bu.Set(3, ku);
		bv.Set(3, kv);
		bu.BasisFunction(0, spt[i][0], 0, uval);
		bv.BasisFunction(0, spt[i][1], 0, vval);
		double Ni = uval[0] * vval[0];
		for (uint j = 0; j < chdu.size(); j++)
		{
			if (coef[j] != 0.)
			{
				ku.assign(chdu[j].begin(), chdu[j].end());
				kv.assign(chdv[j].begin(), chdv[j].end());
				bu.Set(3, ku);
				bv.Set(3, kv);
				bu.BasisFunction(0, spt[i][0], 0, uval);
				bv.BasisFunction(0, spt[i][1], 0, vval);
				Ni -= coef[j] * uval[0] * vval[0];
			}
		}
		if (Ni < 0.)
			Ni = -Ni;
		sum += Ni;
	}
	if (sum / spt.size() < tol)
		return true;
	else
		return false;
}

double TruncatedTspline::CheckFullChildren_2(const vector<array<double, 2>> &spt, const array<double, 5> &motu, const array<double, 5> &motv, const vector<array<double, 5>> &chdu, const vector<array<double, 5>> &chdv, const vector<double> &coef)
{
	double tol(1.e-8);
	double sum(0.);
	for (uint i = 0; i < spt.size(); i++)
	{
		vector<double> ku(motu.begin(), motu.end());
		vector<double> kv(motv.begin(), motv.end());
		vector<double> uval, vval;
		BSplineBasis bu, bv;
		bu.Set(3, ku);
		bv.Set(3, kv);
		bu.BasisFunction(0, spt[i][0], 0, uval);
		bv.BasisFunction(0, spt[i][1], 0, vval);
		double Ni = uval[0] * vval[0];
		for (uint j = 0; j < chdu.size(); j++)
		{
			if (coef[j] != 0.)
			{
				ku.assign(chdu[j].begin(), chdu[j].end());
				kv.assign(chdv[j].begin(), chdv[j].end());
				bu.Set(3, ku);
				bv.Set(3, kv);
				bu.BasisFunction(0, spt[i][0], 0, uval);
				bv.BasisFunction(0, spt[i][1], 0, vval);
				Ni -= coef[j] * uval[0] * vval[0];
			}
		}
		if (Ni < 0.)
			Ni = -Ni;
		sum += Ni;
	}
	return sum;
}

void TruncatedTspline::SetProblem2(string fn)
{
	//read quad vtk
	string fname(fn), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++)
			getline(fin, stmp); //skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
			//cp[i].coor[0]/=12.; cp[i].coor[1]/=12.;
			//cp[i].coor[0]/=6.8; cp[i].coor[1]/=6.8;
			//cp[i].coor[2]=0.;
			cp[i].act = 1;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			tmesh[i].act = 1;
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline::ElementBasis(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	if (tmesh[eid].act == 1)
	{
		if (tmesh[eid].type == 0 || tmesh[eid].type == 1)
		{
			ElementBasis_Regular(eid, u, v, Nt, dNdt);
		}
		else if (tmesh[eid].type == 4)
		{
			//ElementBasis_Irregular(eid,u,v,Nt,dNdt);
			ElementBasis_Irregular_DPatch(eid, u, v, Nt, dNdt);
			//ElementBasis_Irregular_DPatch22(eid, u, v, Nt, dNdt);
		}
	}
}

void TruncatedTspline::ElementBasis_Regular(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	Nt.clear();
	dNdt.clear();
	Nt.resize(tmesh[eid].IEN.size());
	dNdt.resize(tmesh[eid].IEN.size());
	vector<double> Nt0(tmesh[eid].IEN.size());
	vector<array<double, 2>> dNdt0(tmesh[eid].IEN.size());
	vector<double> ku(5, 0.), kv(5, 0.), uval, vval;
	BSplineBasis bu, bv;
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		ku.assign(tmesh[eid].patch_ku[i].begin(), tmesh[eid].patch_ku[i].end());
		kv.assign(tmesh[eid].patch_kv[i].begin(), tmesh[eid].patch_kv[i].end());
		bu.Set(3, ku);
		bv.Set(3, kv);
		bu.BasisFunction(0, u, 1, uval);
		bv.BasisFunction(0, v, 1, vval);
		Nt0[i] = uval[0] * vval[0];
		dNdt0[i][0] = uval[1] * vval[0];
		dNdt0[i][1] = uval[0] * vval[1];
	}
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		Nt[i] = Nt0[i];
		dNdt[i][0] = dNdt0[i][0];
		dNdt[i][1] = dNdt0[i][1];
		if (cp[tmesh[eid].IEN[i]].trun == 1)
		{
			int pid(tmesh[eid].IEN[i]);
			for (uint j = 0; j < cp[pid].tbf.size(); j++)
			{
				vector<int>::iterator it = find(tmesh[eid].IEN.begin(), tmesh[eid].IEN.end(), cp[pid].tbf[j]);
				if (it != tmesh[eid].IEN.end())
				{
					int loc(it - tmesh[eid].IEN.begin());
					Nt[i] -= cp[pid].tc[j] * Nt0[loc];
					dNdt[i][0] -= cp[pid].tc[j] * dNdt0[loc][0];
					dNdt[i][1] -= cp[pid].tc[j] * dNdt0[loc][1];
				}
			}
		}
	}
}

void TruncatedTspline::ElementBasis_Irregular(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	//need to change variabe for Bezier, (u,v) -> (u_b,v_b)
	Nt.clear();
	dNdt.clear();
	Nt.resize(tmesh[eid].IEN.size(), 0.);
	dNdt.resize(tmesh[eid].IEN.size());
	vector<double> Nt1(tmesh[eid].IEN.size(), 0.);
	vector<array<double, 2>> dNdt1(tmesh[eid].IEN.size());
	uint nv(cp[tmesh[eid].cnct[0]].face.size());
	vector<vector<double>> bmat(tmesh[eid].bemat.size(), vector<double>(tmesh[eid].bemat[0].size()));
	for (uint i = 0; i < tmesh[eid].bemat.size(); i++)
	{
		for (uint j = 0; j < tmesh[eid].bemat[i].size(); j++)
		{
			bmat[i][j] = tmesh[eid].bemat[i][j];
		}
	}
	//SetBezier4TranMatOP(nv,bmat);
	BezierElement be;
	double Nt0[25], dNdt0[25][2];
	double u_b(u / tmedge[tmesh[eid].edge[0]].len), v_b(v / tmedge[tmesh[eid].edge[3]].len);
	be.Basis4(u_b, v_b, Nt0, dNdt0);
	for (uint i = 0; i < 2 * nv + 1; i++)
	{
		dNdt1[i][0] = 0.;
		dNdt1[i][1] = 0.;
		for (int j = 0; j < 25; j++)
		{
			Nt1[i] += bmat[i][j] * Nt0[j];
			dNdt1[i][0] += bmat[i][j] * dNdt0[j][0];
			dNdt1[i][1] += bmat[i][j] * dNdt0[j][1];
		}
		dNdt1[i][0] /= tmedge[tmesh[eid].edge[0]].len;
		dNdt1[i][1] /= tmedge[tmesh[eid].edge[3]].len;
	}
	vector<double> ku(5, 0.), kv(5, 0.), uval, vval;
	BSplineBasis bu, bv;
	for (uint i = 2 * nv + 1; i < tmesh[eid].IEN.size(); i++)
	{
		ku.assign(tmesh[eid].patch_ku[i - (2 * nv + 1)].begin(), tmesh[eid].patch_ku[i - (2 * nv + 1)].end());
		kv.assign(tmesh[eid].patch_kv[i - (2 * nv + 1)].begin(), tmesh[eid].patch_kv[i - (2 * nv + 1)].end());
		bu.Set(3, ku);
		bv.Set(3, kv);
		bu.BasisFunction(0, u, 1, uval);
		bv.BasisFunction(0, v, 1, vval);
		Nt1[i] = uval[0] * vval[0];
		dNdt1[i][0] = uval[1] * vval[0];
		dNdt1[i][1] = uval[0] * vval[1];
	}
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		Nt[i] = Nt1[i];
		dNdt[i][0] = dNdt1[i][0];
		dNdt[i][1] = dNdt1[i][1];
		if (cp[tmesh[eid].IEN[i]].trun == 1)
		{
			int pid(tmesh[eid].IEN[i]);
			for (uint j = 0; j < cp[pid].tbf.size(); j++)
			{
				vector<int>::iterator it = find(tmesh[eid].IEN.begin(), tmesh[eid].IEN.end(), cp[pid].tbf[j]);
				if (it != tmesh[eid].IEN.end())
				{
					int loc(it - tmesh[eid].IEN.begin());
					Nt[i] -= cp[pid].tc[j] * Nt1[loc];
					dNdt[i][0] -= cp[pid].tc[j] * dNdt1[loc][0];
					dNdt[i][1] -= cp[pid].tc[j] * dNdt1[loc][1];
				}
			}
		}
	}
}

void TruncatedTspline::SetBezierMatIrrPatch()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1 && tmesh[eid].type == 4)
		{
			tmesh[eid].bemat.clear();
			uint nv(cp[tmesh[eid].cnct[0]].face.size());
			//vector<vector<double>> bmat;
			SetBezier4TranMatOP(nv, tmesh[eid].bemat);
			//SetBezier4TranMat(nv,tmesh[eid].bemat);
			vector<array<double, 5>> patch_ku, patch_kv;
			FindPatchKnotVector_Irr(eid, patch_ku, patch_kv);
			vector<vector<double>> coef(patch_ku.size(), vector<double>(16, 0.));
			vector<vector<double>> coef1(patch_ku.size(), vector<double>(25, 0.));
			array<double, 2> ktsU = {0., tmedge[tmesh[eid].edge[0]].len};
			array<double, 2> ktsV = {0., tmedge[tmesh[eid].edge[3]].len};
			double bzku[6] = {ktsU[0], ktsU[0], ktsU[0], ktsU[1], ktsU[1], ktsU[1]};
			double bzkv[6] = {ktsV[0], ktsV[0], ktsV[0], ktsV[1], ktsV[1], ktsV[1]};
			for (uint i = 0; i < patch_ku.size(); i++)
			{
				if (patch_ku[i][0] < ktsU[1] && patch_ku[i][4] > ktsU[0] && patch_kv[i][0] < ktsV[1] && patch_kv[i][4] > ktsV[0])
				{
					vector<double> ku(patch_ku[i].begin(), patch_ku[i].end());
					vector<double> kv(patch_kv[i].begin(), patch_kv[i].end());
					vector<double> ku1, kv1;
					vector<vector<double>> Tu, Tv;
					BezierInsertKnots(ku, ktsU, ku1);
					BezierInsertKnots(kv, ktsV, kv1);
					TMatrix(ku, ku1, 3, Tu);
					TMatrix(kv, kv1, 3, Tv);
					vector<double>::iterator it1 = search(ku1.begin(), ku1.end(), bzku, bzku + 6);
					vector<double>::iterator it2 = search(kv1.begin(), kv1.end(), bzkv, bzkv + 6);
					//if(it1!=ku1.end() && it2!=kv1.end())
					//{
					int loc1(it1 - ku1.begin() - 1), loc2(it2 - kv1.begin() - 1), count(0);
					for (int j = 0; j < 4; j++)
					{
						for (int k = 0; k < 4; k++)
						{
							coef[i][count] = Tu[loc1 + k][0] * Tv[loc2 + j][0];
							count++;
						}
					}
					//}
				}
			}
			vector<vector<double>> demat;
			DegreeElevate(demat);
			for (uint i = 0; i < patch_ku.size(); i++)
			{
				for (int j = 0; j < 16; j++)
				{
					for (int k = 0; k < 25; k++)
					{
						coef1[i][k] += coef[i][j] * demat[k][j];
					}
				}
			}
			int ploc[5] = {2, 3, 4, 5, 6};
			int bzloc[16] = {3, 4, 8, 9, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
			for (int i = 0; i < 5; i++)
			{
				if (cp[tmesh[eid].IEN[ploc[i]]].trun == 0)
				{
					for (int j = 0; j < 16; j++)
					{
						tmesh[eid].bemat[ploc[i]][bzloc[j]] = coef1[i][bzloc[j]];
					}
				}
			}
		}
	}
}

void TruncatedTspline::FindPatchKnotVector_Irr(int eid, vector<array<double, 5>> &patch_ku, vector<array<double, 5>> &patch_kv)
{
	patch_ku.clear();
	patch_kv.clear();
	patch_ku.resize(5);
	patch_kv.resize(5);
	//find knot vectors for node 2 to 6
	if (tmesh[eid].act == 1 && tmesh[eid].type == 4)
	{
		//first for node 2
		array<double, 2> uv_ref0 = {tmesh[eid].lcs[3].u[0], tmesh[eid].lcs[3].u[1]}, uv0;
		int fcid0(tmedge[tmesh[eid].edge[3]].face[0]), rot0;
		if (fcid0 == eid)
			fcid0 = tmedge[tmesh[eid].edge[3]].face[1];
		FindRotateAndUVCoor(tmesh[eid].cnct[3], tmesh[eid].lcs[3].rot, uv_ref0, fcid0, tmesh[fcid0].cnct[2], rot0, uv0);
		FindLocalKnotVector_1(tmesh[fcid0].cnct[2], rot0, uv0, patch_ku[0], patch_kv[0]);
		//node 3, 4, 5
		for (int i = 3; i > 0; i--)
		{
			array<double, 2> uv1 = {tmesh[eid].lcs[i].u[0], tmesh[eid].lcs[i].u[1]};
			FindLocalKnotVector_1(tmesh[eid].cnct[i], tmesh[eid].lcs[i].rot, uv1, patch_ku[4 - i], patch_kv[4 - i]);
		}
		//node 6
		array<double, 2> uv_ref2 = {tmesh[eid].lcs[1].u[0], tmesh[eid].lcs[1].u[1]}, uv2;
		int fcid2(tmedge[tmesh[eid].edge[0]].face[0]), rot2;
		if (fcid2 == eid)
			fcid2 = tmedge[tmesh[eid].edge[0]].face[1];
		FindRotateAndUVCoor(tmesh[eid].cnct[1], tmesh[eid].lcs[1].rot, uv_ref2, fcid2, tmesh[fcid2].cnct[2], rot2, uv2);
		FindLocalKnotVector_1(tmesh[fcid2].cnct[2], rot2, uv2, patch_ku[4], patch_kv[4]);
	}
}

void TruncatedTspline::SurfacePointMap(int eid, double u, double v, array<double, 3> &pt, array<double, 3> &norm)
{
	vector<double> Nt;
	vector<array<double, 2>> dNdt;
	ElementBasis(eid, u, v, Nt, dNdt);
	//ElementBasisDual(eid,u,v,Nt,dNdt);
	//ElementBasisCC(eid, u, v, Nt, dNdt);
	pt[0] = 0.;
	pt[1] = 0.;
	pt[2] = 0.;
	double nmtmp[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		pt[0] += cp[tmesh[eid].IEN[i]].coor[0] * Nt[i];
		pt[1] += cp[tmesh[eid].IEN[i]].coor[1] * Nt[i];
		pt[2] += cp[tmesh[eid].IEN[i]].coor[2] * Nt[i];
		nmtmp[0][0] += cp[tmesh[eid].IEN[i]].coor[0] * dNdt[i][0];
		nmtmp[0][1] += cp[tmesh[eid].IEN[i]].coor[1] * dNdt[i][0];
		nmtmp[0][2] += cp[tmesh[eid].IEN[i]].coor[2] * dNdt[i][0];
		nmtmp[1][0] += cp[tmesh[eid].IEN[i]].coor[0] * dNdt[i][1];
		nmtmp[1][1] += cp[tmesh[eid].IEN[i]].coor[1] * dNdt[i][1];
		nmtmp[1][2] += cp[tmesh[eid].IEN[i]].coor[2] * dNdt[i][1];
	}
	//norm[0]=nmtmp[0][1]*nmtmp[1][2]-nmtmp[0][2]*nmtmp[1][1];
	//norm[1]=nmtmp[0][2]*nmtmp[1][0]-nmtmp[0][0]*nmtmp[1][2];
	//norm[2]=nmtmp[0][0]*nmtmp[1][1]-nmtmp[0][1]*nmtmp[1][0];
	//double len=sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
	//norm[0]/=len; norm[1]/=len; norm[2]/=len;
}

void TruncatedTspline::ElementRefine_Unstruct_4(int eid)
{
	int pid[5];
	int edid[12];
	Vertex ptmp1;
	//ptmp1.update = 1;
	//ptmp1.coortmp[0]=(cp[tmesh[eid].cnct[0]].coor[0]+cp[tmesh[eid].cnct[1]].coor[0]+cp[tmesh[eid].cnct[2]].coor[0]+cp[tmesh[eid].cnct[3]].coor[0])/4.;
	//ptmp1.coortmp[1]=(cp[tmesh[eid].cnct[0]].coor[1]+cp[tmesh[eid].cnct[1]].coor[1]+cp[tmesh[eid].cnct[2]].coor[1]+cp[tmesh[eid].cnct[3]].coor[1])/4.;
	//ptmp1.coortmp[2]=(cp[tmesh[eid].cnct[0]].coor[2]+cp[tmesh[eid].cnct[1]].coor[2]+cp[tmesh[eid].cnct[2]].coor[2]+cp[tmesh[eid].cnct[3]].coor[2])/4.;
	cp.push_back(ptmp1);
	pid[0] = cp.size() - 1;
	for (int j = 0; j < 4; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[j], tmesh[eid].cnct[(j + 1) % 4]};
		if (tmedge[tmesh[eid].edge[j]].act == 1)
		{
			Vertex ptmp;
			//ptmp.update = 1;
			//ptmp.coortmp[0]=(cp[itmp[0]].coor[0]+cp[itmp[1]].coor[0])/2.;
			//ptmp.coortmp[1]=(cp[itmp[0]].coor[1]+cp[itmp[1]].coor[1])/2.;
			//ptmp.coortmp[2]=(cp[itmp[0]].coor[2]+cp[itmp[1]].coor[2])/2.;
			cp.push_back(ptmp);
			pid[j + 1] = cp.size() - 1;
			tmedge[tmesh[eid].edge[j]].midpt = pid[j + 1];
			Edge edtmp1, edtmp2, edtmp3;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j + 1];
			edtmp1.len = tmedge[tmesh[eid].edge[j]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[j];
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j + 1];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[j]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[j];
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[(j + 1) % 4]].len / 2.;
			tmedge.push_back(edtmp1);
			edid[3 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].chd[0] = edid[3 * j];
			tmedge.push_back(edtmp2);
			edid[3 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].chd[1] = edid[3 * j + 1];
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].act = 0;
		}
		else
		{
			pid[j + 1] = tmedge[tmesh[eid].edge[j]].midpt;
			Edge edtmp3;
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[(j + 1) % 4]].len / 2.;
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			int ied(tmedge[tmesh[eid].edge[j]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[3 * j] = ied;
				edid[3 * j + 1] = tmedge[tmesh[eid].edge[j]].chd[1];
			}
			else
			{
				edid[3 * j] = tmedge[tmesh[eid].edge[j]].chd[1];
				edid[3 * j + 1] = ied;
			}
		}
	}

	int e_cnct[4][4] = {{tmesh[eid].cnct[0], pid[1], pid[0], pid[4]}, {pid[1], tmesh[eid].cnct[1], pid[2], pid[0]}, {pid[0], pid[2], tmesh[eid].cnct[2], pid[3]}, {pid[4], pid[0], pid[3], tmesh[eid].cnct[3]}};
	int e_edge[4][4] = {{edid[0], edid[2], edid[11], edid[10]}, {edid[1], edid[3], edid[5], edid[2]}, {edid[5], edid[4], edid[6], edid[8]}, {edid[11], edid[8], edid[7], edid[9]}};
	int enewid[4];
	double chd_org[4][2] = {{0., 0.}, {tmedge[tmesh[eid].edge[0]].len / 2., 0.}, {tmedge[tmesh[eid].edge[0]].len / 2., tmedge[tmesh[eid].edge[1]].len / 2.}, {0., tmedge[tmesh[eid].edge[1]].len / 2.}};
	vector<Element> etmp(4);
	for (int i = 0; i < 4; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 0;
		etmp[i].lv = tmesh[eid].lv + 1.; //lv is used in all kinds of elements in T-mesh
		//etmp[i].lev = tmesh[eid].lev + 1;//lev is used in type-0 elements in T-mesh
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[j] = e_cnct[i][j];
			etmp[i].edge[j] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
		tmesh[eid].chd_o[i][0] = chd_org[i][0];
		tmesh[eid].chd_o[i][1] = chd_org[i][1];
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::ElementRefine_Irregular_4(int eid)
{
	int pid[5];
	int edid[12];
	Vertex ptmp1;
	ptmp1.update = 1;
	ptmp1.coortmp[0] = (cp[tmesh[eid].cnct[0]].coor[0] + cp[tmesh[eid].cnct[1]].coor[0] + cp[tmesh[eid].cnct[2]].coor[0] + cp[tmesh[eid].cnct[3]].coor[0]) / 4.;
	ptmp1.coortmp[1] = (cp[tmesh[eid].cnct[0]].coor[1] + cp[tmesh[eid].cnct[1]].coor[1] + cp[tmesh[eid].cnct[2]].coor[1] + cp[tmesh[eid].cnct[3]].coor[1]) / 4.;
	ptmp1.coortmp[2] = (cp[tmesh[eid].cnct[0]].coor[2] + cp[tmesh[eid].cnct[1]].coor[2] + cp[tmesh[eid].cnct[2]].coor[2] + cp[tmesh[eid].cnct[3]].coor[2]) / 4.;
	ptmp1.wtmp = 1.;
	cp.push_back(ptmp1);
	pid[0] = cp.size() - 1;
	cp[pid[0]].update = 1;
	for (int j = 0; j < 4; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[j], tmesh[eid].cnct[(j + 1) % 4]};
		if (tmedge[tmesh[eid].edge[j]].act == 1)
		{
			Vertex ptmp;
			int id0[4] = {tmesh[eid].cnct[j], tmesh[eid].cnct[(j + 1) % 4], tmesh[eid].cnct[(j + 2) % 4], tmesh[eid].cnct[(j + 3) % 4]};
			int fcnb(tmedge[tmesh[eid].edge[j]].face[0]);
			if (fcnb == eid)
				fcnb = tmedge[tmesh[eid].edge[j]].face[1];
			int *it1 = find(tmesh[fcnb].cnct, tmesh[fcnb].cnct + 4, tmesh[eid].cnct[j]);
			int loc1(it1 - tmesh[fcnb].cnct);
			int id1[2] = {tmesh[fcnb].cnct[(loc1 + 1) % 4], tmesh[fcnb].cnct[(loc1 + 2) % 4]};
			//if(j==0 || j==3)
			//{
			ptmp.coortmp[0] = (6. * cp[id0[0]].coor[0] + 6. * cp[id0[1]].coor[0] + cp[id0[2]].coor[0] + cp[id0[3]].coor[0] + cp[id1[0]].coor[0] + cp[id1[1]].coor[0]) / 16.;
			ptmp.coortmp[1] = (6. * cp[id0[0]].coor[1] + 6. * cp[id0[1]].coor[1] + cp[id0[2]].coor[1] + cp[id0[3]].coor[1] + cp[id1[0]].coor[1] + cp[id1[1]].coor[1]) / 16.;
			ptmp.coortmp[2] = (6. * cp[id0[0]].coor[2] + 6. * cp[id0[1]].coor[2] + cp[id0[2]].coor[2] + cp[id0[3]].coor[2] + cp[id1[0]].coor[2] + cp[id1[1]].coor[2]) / 16.;
			ptmp.wtmp = 1.;
			//}
			//else
			//{
			//	ptmp.coor[0]=cp[tmesh[eid].cnct[0]].coor[0]/16.;
			//	ptmp.coor[1]=cp[tmesh[eid].cnct[0]].coor[1]/16.;
			//	ptmp.coor[2]=cp[tmesh[eid].cnct[0]].coor[2]/16.;
			//}
			cp.push_back(ptmp);
			pid[j + 1] = cp.size() - 1;
			tmedge[tmesh[eid].edge[j]].midpt = pid[j + 1];
			Edge edtmp1, edtmp2, edtmp3;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j + 1];
			edtmp1.len = tmedge[tmesh[eid].edge[j]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[j];
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j + 1];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[j]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[j];
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[(j + 1) % 4]].len / 2.;
			tmedge.push_back(edtmp1);
			edid[3 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].chd[0] = edid[3 * j];
			tmedge.push_back(edtmp2);
			edid[3 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].chd[1] = edid[3 * j + 1];
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].act = 0;
			cp[pid[j + 1]].update = 1;
		}
		else
		{
			pid[j + 1] = tmedge[tmesh[eid].edge[j]].midpt;
			//int id1[2]={tmesh[eid].cnct[(j+2)%4],tmesh[eid].cnct[(j+3)%4]};
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
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[(j + 1) % 4]].len / 2.;
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			int ied(tmedge[tmesh[eid].edge[j]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[3 * j] = ied;
				edid[3 * j + 1] = tmedge[tmesh[eid].edge[j]].chd[1];
			}
			else
			{
				edid[3 * j] = tmedge[tmesh[eid].edge[j]].chd[1];
				edid[3 * j + 1] = ied;
			}
		}
	}

	int nvl(cp[tmesh[eid].cnct[0]].face.size());
	if (cp[tmesh[eid].cnct[0]].update == 0)
	{
		double ccf[3] = {1. - 7. / (4. * nvl), 3. / (2. * nvl * nvl), 1. / (4. * nvl * nvl)};
		cp[tmesh[eid].cnct[0]].coortmp[0] = ccf[0] * cp[tmesh[eid].cnct[0]].coor[0] + ccf[1] * cp[tmesh[eid].cnct[1]].coor[0] + ccf[2] * cp[tmesh[eid].cnct[2]].coor[0];
		cp[tmesh[eid].cnct[0]].coortmp[1] = ccf[0] * cp[tmesh[eid].cnct[0]].coor[1] + ccf[1] * cp[tmesh[eid].cnct[1]].coor[1] + ccf[2] * cp[tmesh[eid].cnct[2]].coor[1];
		cp[tmesh[eid].cnct[0]].coortmp[2] = ccf[0] * cp[tmesh[eid].cnct[0]].coor[2] + ccf[1] * cp[tmesh[eid].cnct[1]].coor[2] + ccf[2] * cp[tmesh[eid].cnct[2]].coor[2];
		cp[tmesh[eid].cnct[0]].update = 2;
		cp[tmesh[eid].cnct[0]].wtmp = 1.;
	}
	else if (cp[tmesh[eid].cnct[0]].update == 2)
	{
		double ccf[2] = {3. / (2. * nvl * nvl), 1. / (4. * nvl * nvl)};
		cp[tmesh[eid].cnct[0]].coortmp[0] += ccf[0] * cp[tmesh[eid].cnct[1]].coor[0] + ccf[1] * cp[tmesh[eid].cnct[2]].coor[0];
		cp[tmesh[eid].cnct[0]].coortmp[1] += ccf[0] * cp[tmesh[eid].cnct[1]].coor[1] + ccf[1] * cp[tmesh[eid].cnct[2]].coor[1];
		cp[tmesh[eid].cnct[0]].coortmp[2] += ccf[0] * cp[tmesh[eid].cnct[1]].coor[2] + ccf[1] * cp[tmesh[eid].cnct[2]].coor[2];
	}
	else
	{
		cout << "Not supported for other udpate types!\n";
		getchar();
	}

	int e_cnct[4][4] = {{tmesh[eid].cnct[0], pid[1], pid[0], pid[4]}, {pid[1], tmesh[eid].cnct[1], pid[2], pid[0]}, {pid[0], pid[2], tmesh[eid].cnct[2], pid[3]}, {pid[4], pid[0], pid[3], tmesh[eid].cnct[3]}};
	int e_edge[4][4] = {{edid[0], edid[2], edid[11], edid[10]}, {edid[1], edid[3], edid[5], edid[2]}, {edid[5], edid[4], edid[6], edid[8]}, {edid[11], edid[8], edid[7], edid[9]}};
	int enewid[4];
	int e_type[4] = {4, 0, 0, 0};
	vector<Element> etmp(4);
	for (int i = 0; i < 4; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = e_type[i];
		etmp[i].lv = tmesh[eid].lv + 1.;
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[j] = e_cnct[i][j];
			etmp[i].edge[j] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::ElementRefine_Unstruct_2(int eid, int dir)
{
	int pid[2], edid[7], pos(dir);
	int cnid[2] = {pos, (pos + 2) % 4};
	for (int j = 0; j < 2; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[cnid[j]], tmesh[eid].cnct[(cnid[j] + 1) % 4]};
		if (tmedge[tmesh[eid].edge[cnid[j]]].act == 1)
		{
			Vertex ptmp;
			cp.push_back(ptmp);
			pid[j] = cp.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].midpt = pid[j];
			Edge edtmp1, edtmp2;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j];
			edtmp1.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[cnid[j]];
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[cnid[j]];
			tmedge.push_back(edtmp1);
			edid[2 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[0] = edid[2 * j];
			tmedge.push_back(edtmp2);
			edid[2 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[1] = edid[2 * j + 1];
			tmedge[tmesh[eid].edge[cnid[j]]].act = 0;
		}
		else
		{
			pid[j] = tmedge[tmesh[eid].edge[cnid[j]]].midpt;
			int ied(tmedge[tmesh[eid].edge[cnid[j]]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[2 * j] = ied;
				edid[2 * j + 1] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
			}
			else
			{
				edid[2 * j] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
				edid[2 * j + 1] = ied;
			}
		}
	}
	Edge edtmp;
	edtmp.act = 1;
	edtmp.pt[0] = pid[0];
	edtmp.pt[1] = pid[1];
	edtmp.len = tmedge[tmesh[eid].edge[(pos + 1) % 4]].len;
	tmedge.push_back(edtmp);
	edid[4] = tmedge.size() - 1;
	edid[5] = tmesh[eid].edge[(pos + 1) % 4];
	edid[6] = tmesh[eid].edge[(pos + 3) % 4];

	int e_cnct[2][4] = {{tmesh[eid].cnct[pos], pid[0], pid[1], tmesh[eid].cnct[(pos + 3) % 4]}, {pid[0], tmesh[eid].cnct[(pos + 1) % 4], tmesh[eid].cnct[(pos + 2) % 4], pid[1]}};
	int e_edge[2][4] = {{edid[0], edid[4], edid[3], edid[6]}, {edid[1], edid[5], edid[2], edid[4]}};
	int enewid[2];
	double chd_org[2][2] = {{0., 0.}, {tmedge[tmesh[eid].edge[0]].len / 2., 0.}};
	if (dir == 1)
	{
		chd_org[1][0] = 0.;
		chd_org[1][1] = tmedge[tmesh[eid].edge[1]].len / 2.;
	}
	vector<Element> etmp(2);
	for (int i = 0; i < 2; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 0;
		etmp[i].lv = tmesh[eid].lv + 0.5;
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[(pos + j) % 4] = e_cnct[i][j];
			etmp[i].edge[(pos + j) % 4] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
		tmesh[eid].chd_o[i][0] = chd_org[i][0];
		tmesh[eid].chd_o[i][1] = chd_org[i][1];
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::ElementRefine_Unstruct_b(int eid)
{
	int pid[2], edid[7], pos(0);
	if (tmedge[tmesh[eid].edge[0]].len == 0.)
		pos = 1;
	int cnid[2] = {pos, (pos + 2) % 4};
	for (int j = 0; j < 2; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[cnid[j]], tmesh[eid].cnct[(cnid[j] + 1) % 4]};
		if (tmedge[tmesh[eid].edge[cnid[j]]].act == 1)
		{
			Vertex ptmp;
			//ptmp.coor[0]=(cp[itmp[0]].coor[0]+cp[itmp[1]].coor[0])/2.;
			//ptmp.coor[1]=(cp[itmp[0]].coor[1]+cp[itmp[1]].coor[1])/2.;
			//ptmp.coor[2]=(cp[itmp[0]].coor[2]+cp[itmp[1]].coor[2])/2.;
			cp.push_back(ptmp);
			pid[j] = cp.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].midpt = pid[j];
			Edge edtmp1, edtmp2;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j];
			edtmp1.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[cnid[j]];
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[cnid[j]];
			tmedge.push_back(edtmp1);
			edid[2 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[0] = edid[2 * j];
			tmedge.push_back(edtmp2);
			edid[2 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[1] = edid[2 * j + 1];
			tmedge[tmesh[eid].edge[cnid[j]]].act = 0;
		}
		else
		{
			pid[j] = tmedge[tmesh[eid].edge[cnid[j]]].midpt;
			int ied(tmedge[tmesh[eid].edge[cnid[j]]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[2 * j] = ied;
				edid[2 * j + 1] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
			}
			else
			{
				edid[2 * j] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
				edid[2 * j + 1] = ied;
			}
		}
	}
	Edge edtmp;
	edtmp.act = 1;
	edtmp.pt[0] = pid[0];
	edtmp.pt[1] = pid[1];
	edtmp.len = tmedge[tmesh[eid].edge[(pos + 1) % 4]].len;
	tmedge.push_back(edtmp);
	edid[4] = tmedge.size() - 1;
	edid[5] = tmesh[eid].edge[(pos + 1) % 4];
	edid[6] = tmesh[eid].edge[(pos + 3) % 4];

	int e_cnct[2][4] = {{tmesh[eid].cnct[pos], pid[0], pid[1], tmesh[eid].cnct[(pos + 3) % 4]}, {pid[0], tmesh[eid].cnct[(pos + 1) % 4], tmesh[eid].cnct[(pos + 2) % 4], pid[1]}};
	int e_edge[2][4] = {{edid[0], edid[4], edid[3], edid[6]}, {edid[1], edid[5], edid[2], edid[4]}};
	int enewid[2];
	double chd_org[2][2] = {{0., 0.}, {tmedge[tmesh[eid].edge[0]].len / 2., 0.}};
	if (pos == 1)
	{
		chd_org[1][0] = 0.;
		chd_org[1][1] = tmedge[tmesh[eid].edge[1]].len / 2.;
	}
	vector<Element> etmp(2);
	for (int i = 0; i < 2; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 2;
		etmp[i].lv = tmesh[eid].lv + 1.;
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[(pos + j) % 4] = e_cnct[i][j];
			etmp[i].edge[(pos + j) % 4] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
		tmesh[eid].chd_o[i][0] = chd_org[i][0];
		tmesh[eid].chd_o[i][1] = chd_org[i][1];
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::Refine_Unstruct(const vector<int> &rfid, const vector<int> &rftype)
{
	//initialize for refinement
	npt_old = cp.size();
	nel_old = tmesh.size();
	for (uint i = 0; i < cp.size(); i++)
	{
		cp[i].update = 0;
		cp[i].aff = 0;
		cp[i].coortmp[0] = 0.;
		cp[i].coortmp[1] = 0.;
		cp[i].coortmp[2] = 0.;
		for (int j = 0; j < 4; j++)
		{
			cp[i].kitvUtmp[j] = 0.;
			cp[i].kitvVtmp[j] = 0.;
		}
		cp[i].truntmp = 0;
		vector<int>().swap(cp[i].tbftmp);
		vector<double>().swap(cp[i].tctmp);
	}
	for (uint i = 0; i < tmesh.size(); i++)
	{
		//tmesh[i].aff=0;
		vector<int>().swap(tmesh[i].IENtmp);
		vector<array<double, 5>>().swap(tmesh[i].patch_kutmp);
		vector<array<double, 5>>().swap(tmesh[i].patch_kvtmp);
	}
	//refine
	for (uint i = 0; i < rfid.size(); i++)
	{
		if (rftype[i] == 0)
		{
			ElementRefine_Unstruct_4(rfid[i]);
		}
		else if (rftype[i] == 1)
		{
			ElementRefine_Unstruct_2(rfid[i], rftype[i] - 1);
		}
		else if (rftype[i] == 2)
		{
			ElementRefine_Unstruct_2(rfid[i], rftype[i] - 1);
		}
		else if (rftype[i] == 3)
		{
			ElementRefine_Unstruct_b(rfid[i]);
		}
		else if (rftype[i] == 4)
		{
			//ElementRefine_Irregular_4(rfid[i]);
			ElementRefine_Invalid_4(rfid[i]);
		}
		else
		{
			cout << "Other types of elements are not supported to be refined!\n";
			getchar();
		}
	}

	UpdateConnect_1();
	FindEdgeTopoDirec_1();
	FindKnotInterval_1();
	UpdateKnotInterval_1(); //might not work
	SetLocalCoorSystem();

	//Find IENtmp
	//FindIEN_2();
	//FindIEN_3();
	FindIEN_4(); //Dpatch

	//calculate control points
	for (uint i = 0; i < rfid.size(); i++)
	{
		if (tmesh[rfid[i]].type == 0 || tmesh[rfid[i]].type == 1 || tmesh[rfid[i]].type == 2)
		{
			UpdatePatchCP_Unstruct_1(rfid[i]);
		}
	}
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].update == 2 || cp[i].update == 3) //EP
		{
			cp[i].coor[0] = cp[i].coortmp[0];
			cp[i].coor[1] = cp[i].coortmp[1];
			cp[i].coor[2] = cp[i].coortmp[2];
		}
	}
	//update basis functions
	for (uint i = 0; i < cp.size(); i++)
	{
		if (i < npt_old)
		{
			if (cp[i].truntmp == 0 && cp[i].aff == 1)
			{
				cp[i].coor[0] = cp[i].coortmp[0];
				cp[i].coor[1] = cp[i].coortmp[1];
				cp[i].coor[2] = cp[i].coortmp[2];
				//for(int j=0; j<4; j++)
				//{
				//	cp[i].kitvU[j]=cp[i].kitvUtmp[j];
				//	cp[i].kitvV[j]=cp[i].kitvVtmp[j];
				//}
			}
			else if (cp[i].truntmp == 1)
			{
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
			//for(int j=0; j<4; j++)
			//{
			//	cp[i].kitvU[j]=cp[i].kitvUtmp[j];
			//	cp[i].kitvV[j]=cp[i].kitvVtmp[j];
			//}
		}
	}

	//for(uint i=0; i<npt_old; i++)
	//{
	//	if(cp[i].update==2)
	//	{
	//		cp[i].coor[0]=cp[i].coortmp[0];
	//		cp[i].coor[1]=cp[i].coortmp[1];
	//		cp[i].coor[2]=cp[i].coortmp[2];
	//		for(int j=0; j<4; j++)
	//		{
	//			cp[i].kitvU[j]=cp[i].kitvUtmp[j];
	//			cp[i].kitvV[j]=cp[i].kitvVtmp[j];
	//		}
	//	}
	//	else if(cp[i].update==3)
	//	{
	//		cp[i].trun=1;
	//		cp[i].tbf.clear();
	//		cp[i].tc.clear();
	//		cp[i].tbf=cp[i].tbftmp;
	//		cp[i].tc=cp[i].tctmp;
	//		cp[i].tbftmp.clear();
	//		cp[i].tctmp.clear();
	//	}
	//}

	//update new elements
	for (uint i = 0; i < rfid.size(); i++)
	{
		if (tmesh[rfid[i]].type == 0 || tmesh[rfid[i]].type == 1 || tmesh[rfid[i]].type == 2)
		{
			for (int j = 0; j < 4; j++)
			{
				int chdid(tmesh[rfid[i]].chd[j]);
				if (chdid != -1)
				{
					tmesh[chdid].IEN.clear();
					tmesh[chdid].patch_ku.clear();
					tmesh[chdid].patch_kv.clear();
					for (uint k = 0; k < tmesh[chdid].IENtmp.size(); k++)
					{
						tmesh[chdid].IEN.push_back(tmesh[chdid].IENtmp[k]);
						if (cp[tmesh[chdid].IENtmp[k]].trun == 0)
						{
							tmesh[chdid].patch_ku.push_back(tmesh[chdid].patch_kutmp[k]);
							tmesh[chdid].patch_kv.push_back(tmesh[chdid].patch_kvtmp[k]);
						}
						else
						{
							vector<int>::iterator it = find(tmesh[rfid[i]].IEN.begin(), tmesh[rfid[i]].IEN.end(), tmesh[chdid].IENtmp[k]);
							if (it != tmesh[rfid[i]].IEN.end())
							{
								int loc(it - tmesh[rfid[i]].IEN.begin());
								array<double, 5> kutmp, kvtmp;
								for (int a = 0; a < 5; a++)
								{
									kutmp[a] = tmesh[rfid[i]].patch_ku[loc][a] - tmesh[rfid[i]].chd_o[j][0];
									kvtmp[a] = tmesh[rfid[i]].patch_kv[loc][a] - tmesh[rfid[i]].chd_o[j][1];
								}
								tmesh[chdid].patch_ku.push_back(kutmp);
								tmesh[chdid].patch_kv.push_back(kvtmp);
							}
							else
							{
								cout << "Cannot find truncated ID in the old set!\n";
								getchar();
							}
						}
					}
					//for(uint k=0; k<tmesh[chdid].IEN.size(); k++)
					//{
					//	cout<<tmesh[chdid].IEN[k]<<" ";
					//}
					//cout<<"\n";
					//getchar();
				}
			}
		}
		else if (tmesh[rfid[i]].type == 4)
		{
			for (int j = 0; j < 4; j++)
			{
				int chdid(tmesh[rfid[i]].chd[j]);
				if (chdid != -1)
				{
					tmesh[chdid].IEN.clear();
					tmesh[chdid].patch_ku.clear();
					tmesh[chdid].patch_kv.clear();
					tmesh[chdid].IEN = tmesh[chdid].IENtmp;
					tmesh[chdid].patch_ku = tmesh[chdid].patch_kutmp;
					tmesh[chdid].patch_kv = tmesh[chdid].patch_kvtmp;
					vector<int>().swap(tmesh[chdid].IENtmp);
					vector<array<double, 5>>().swap(tmesh[chdid].patch_kutmp);
					vector<array<double, 5>>().swap(tmesh[chdid].patch_kvtmp);
				}
			}
		}
	}

	//update old elements
	for (uint i = 0; i < nel_old; i++)
	{
		if (tmesh[i].act == 1 && (tmesh[i].type == 0 || tmesh[i].type == 1 || tmesh[i].type == 2))
		{
			vector<int> ien_old(tmesh[i].IEN);
			vector<array<double, 5>> ku_old(tmesh[i].patch_ku);
			vector<array<double, 5>> kv_old(tmesh[i].patch_kv);
			tmesh[i].IEN.clear();
			tmesh[i].patch_ku.clear();
			tmesh[i].patch_kv.clear();
			for (uint j = 0; j < tmesh[i].IENtmp.size(); j++)
			{
				tmesh[i].IEN.push_back(tmesh[i].IENtmp[j]);
				if (cp[tmesh[i].IENtmp[j]].trun == 0)
				{
					tmesh[i].patch_ku.push_back(tmesh[i].patch_kutmp[j]);
					tmesh[i].patch_kv.push_back(tmesh[i].patch_kvtmp[j]);
				}
				else
				{
					vector<int>::iterator it = find(ien_old.begin(), ien_old.end(), tmesh[i].IENtmp[j]);
					int loc(it - ien_old.begin());
					tmesh[i].patch_ku.push_back(ku_old[loc]);
					tmesh[i].patch_kv.push_back(kv_old[loc]);
				}
			}
			vector<int>().swap(tmesh[i].IENtmp);
			vector<array<double, 5>>().swap(tmesh[i].patch_kutmp);
			vector<array<double, 5>>().swap(tmesh[i].patch_kvtmp);
		}
		else if (tmesh[i].act == 1 && tmesh[i].type == 4)
		{
			vector<int> ien_old(tmesh[i].IEN);
			vector<array<double, 5>> ku_old(tmesh[i].patch_ku);
			vector<array<double, 5>> kv_old(tmesh[i].patch_kv);
			tmesh[i].IEN.clear();
			tmesh[i].patch_ku.clear();
			tmesh[i].patch_kv.clear();
			uint n_1r(2 * cp[tmesh[i].cnct[0]].face.size() + 1);
			for (uint j = 0; j < tmesh[i].IENtmp.size(); j++)
			{
				tmesh[i].IEN.push_back(tmesh[i].IENtmp[j]);
				//if(j>=n_1r)
				//{
				//	if(cp[tmesh[i].IENtmp[j]].trun==0)
				//	{
				//		tmesh[i].patch_ku.push_back(tmesh[i].patch_kutmp[j-n_1r]);
				//		tmesh[i].patch_kv.push_back(tmesh[i].patch_kvtmp[j-n_1r]);
				//	}
				//	else
				//	{
				//		vector<int>::iterator it=find(ien_old.begin(),ien_old.end(),tmesh[i].IENtmp[j]);
				//		int loc(it-ien_old.begin());
				//		tmesh[i].patch_ku.push_back(ku_old[loc-n_1r]);
				//		tmesh[i].patch_kv.push_back(kv_old[loc-n_1r]);
				//	}
				//}
			}
			vector<int>().swap(tmesh[i].IENtmp);
			vector<array<double, 5>>().swap(tmesh[i].patch_kutmp);
			vector<array<double, 5>>().swap(tmesh[i].patch_kvtmp);
		}
	}
}

void TruncatedTspline::Topo_Refine_Unstruct(const vector<int> &rfid, const vector<int> &rftype, vector<int> &rfid_more, vector<int> &rftype_more)
{
	//initialize for refinement
	npt_old = cp.size();
	nel_old = tmesh.size();
	rfid_more.clear();
	rftype_more.clear();
	rfid_more = rfid;
	rftype_more = rftype;
	for (uint i = 0; i < cp.size(); i++)
	{
		cp[i].update = 0; //initialized here!
		cp[i].aff = 0;
		cp[i].coortmp[0] = 0.;
		cp[i].coortmp[1] = 0.;
		cp[i].coortmp[2] = 0.;
		for (int j = 0; j < 4; j++)
		{
			cp[i].kitvUtmp[j] = 0.;
			cp[i].kitvVtmp[j] = 0.;
		}
		cp[i].truntmp = 0;
		vector<int>().swap(cp[i].tbftmp);
		vector<double>().swap(cp[i].tctmp);
	}
	for (uint i = 0; i < tmesh.size(); i++)
	{
		//tmesh[i].aff=0;
		vector<int>().swap(tmesh[i].IENtmp);
		vector<array<double, 5>>().swap(tmesh[i].patch_kutmp);
		vector<array<double, 5>>().swap(tmesh[i].patch_kvtmp);
	}
	//refine
	for (uint i = 0; i < rfid.size(); i++)
	{
		if (rftype[i] == 0)
		{
			ElementRefine_Unstruct_4(rfid[i]);
		}
		else if (rftype[i] == 1)
		{
			ElementRefine_Unstruct_2(rfid[i], rftype[i] - 1);
		}
		else if (rftype[i] == 2)
		{
			ElementRefine_Unstruct_2(rfid[i], rftype[i] - 1);
		}
		else if (rftype[i] == 3)
		{
			ElementRefine_Unstruct_b(rfid[i]);
		}
		else if (rftype[i] == 4)
		{
			ElementRefine_Irregular_4(rfid[i]);
		}
		else if (rftype[i] == 5)
		{
			ElementRefine_Invalid_4(rfid[i]);
		}
		else
		{
			cout << "Other types of elements are not supported to be refined!\n";
			getchar();
		}
	}

	vector<int> rid_b;
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			if (tmesh[i].type == 2)
			{
				for (int j = 0; j < 4; j++)
				{
					if (tmedge[tmesh[i].edge[j]].act == 0)
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
	for (uint i = 0; i < rid_b.size(); i++)
	{
		ElementRefine_Unstruct_b(rid_b[i]);
	}

	//for (uint i = 0; i<cp.size(); i++)
	//{
	//	if (i>=npt_old)
	//	{
	//		cp[i].coor[0] = cp[i].coortmp[0];
	//		cp[i].coor[1] = cp[i].coortmp[1];
	//		cp[i].coor[2] = cp[i].coortmp[2];
	//		cp[i].w = cp[i].wtmp;
	//	}
	//}
}

void TruncatedTspline::Geom_Refine_Unstruct(const vector<int> &rfid, const vector<int> &rftype)
{
	UpdateConnect_1();
	FindEdgeTopoDirec_1();
	FindKnotInterval_1();
	UpdateKnotInterval_1(); //might not work
	SetLocalCoorSystem();
	FindIEN_3();

	//calculate control points
	for (uint i = 0; i < rfid.size(); i++)
	{
		if (tmesh[rfid[i]].type == 0 || tmesh[rfid[i]].type == 1 || tmesh[rfid[i]].type == 2)
		{
			UpdatePatchCP_Unstruct_1(rfid[i]);
		}
	}
	//update basis functions
	for (uint i = 0; i < cp.size(); i++)
	{
		if (i < npt_old)
		{
			if (cp[i].truntmp == 0 && cp[i].aff == 1)
			{
				cp[i].coor[0] = cp[i].coortmp[0];
				cp[i].coor[1] = cp[i].coortmp[1];
				cp[i].coor[2] = cp[i].coortmp[2];
				cp[i].w = cp[i].wtmp;
				//for(int j=0; j<4; j++)
				//{
				//	cp[i].kitvU[j]=cp[i].kitvUtmp[j];
				//	cp[i].kitvV[j]=cp[i].kitvVtmp[j];
				//}
			}
			else if (cp[i].truntmp == 1)
			{
				//if(i==71)
				//{
				//	cout<<"here!\n";
				//	getchar();
				//}
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
			cp[i].w = cp[i].wtmp;
			//for(int j=0; j<4; j++)
			//{
			//	cp[i].kitvU[j]=cp[i].kitvUtmp[j];
			//	cp[i].kitvV[j]=cp[i].kitvVtmp[j];
			//}
		}
	}

	//for(uint i=0; i<npt_old; i++)
	//{
	//	if(cp[i].update==2)
	//	{
	//		cp[i].coor[0]=cp[i].coortmp[0];
	//		cp[i].coor[1]=cp[i].coortmp[1];
	//		cp[i].coor[2]=cp[i].coortmp[2];
	//		for(int j=0; j<4; j++)
	//		{
	//			cp[i].kitvU[j]=cp[i].kitvUtmp[j];
	//			cp[i].kitvV[j]=cp[i].kitvVtmp[j];
	//		}
	//	}
	//	else if(cp[i].update==3)
	//	{
	//		cp[i].trun=1;
	//		cp[i].tbf.clear();
	//		cp[i].tc.clear();
	//		cp[i].tbf=cp[i].tbftmp;
	//		cp[i].tc=cp[i].tctmp;
	//		cp[i].tbftmp.clear();
	//		cp[i].tctmp.clear();
	//	}
	//}

	//update new elements
	for (uint i = 0; i < rfid.size(); i++)
	{
		if (tmesh[rfid[i]].type == 0 || tmesh[rfid[i]].type == 1 || tmesh[rfid[i]].type == 2)
		{
			vector<int> trunIEN;
			vector<array<double, 5>> trun_patch_ku;
			vector<array<double, 5>> trun_patch_kv;
			for (uint j = 0; j < tmesh[rfid[i]].IEN.size(); j++)
			{
				if (cp[tmesh[rfid[i]].IEN[j]].trun == 1)
				{
					trunIEN.push_back(tmesh[rfid[i]].IEN[j]);
					trun_patch_ku.push_back(tmesh[rfid[i]].patch_ku[j]);
					trun_patch_kv.push_back(tmesh[rfid[i]].patch_kv[j]);
				}
			}
			for (int j = 0; j < 4; j++)
			{
				int chdid(tmesh[rfid[i]].chd[j]);
				if (chdid != -1)
				{
					tmesh[chdid].IEN.clear();
					tmesh[chdid].patch_ku.clear();
					tmesh[chdid].patch_kv.clear();
					for (uint k = 0; k < tmesh[chdid].IENtmp.size(); k++)
					{
						tmesh[chdid].IEN.push_back(tmesh[chdid].IENtmp[k]);
						if (cp[tmesh[chdid].IENtmp[k]].trun == 0)
						{
							tmesh[chdid].patch_ku.push_back(tmesh[chdid].patch_kutmp[k]);
							tmesh[chdid].patch_kv.push_back(tmesh[chdid].patch_kvtmp[k]);
						}
						else
						{
							vector<int>::iterator it = find(tmesh[rfid[i]].IEN.begin(), tmesh[rfid[i]].IEN.end(), tmesh[chdid].IENtmp[k]);
							if (it != tmesh[rfid[i]].IEN.end())
							{
								int loc(it - tmesh[rfid[i]].IEN.begin());
								array<double, 5> kutmp, kvtmp;
								for (int a = 0; a < 5; a++)
								{
									kutmp[a] = tmesh[rfid[i]].patch_ku[loc][a] - tmesh[rfid[i]].chd_o[j][0];
									kvtmp[a] = tmesh[rfid[i]].patch_kv[loc][a] - tmesh[rfid[i]].chd_o[j][1];
								}
								tmesh[chdid].patch_ku.push_back(kutmp);
								tmesh[chdid].patch_kv.push_back(kvtmp);
							}
							else
							{
								//cout<<i<<" "<<rfid[i]<<" "<<tmesh[rfid[i]].IEN.size()<<"\n";
								cout << "Cannot find truncated ID in the old set!\n";
								getchar();
							}
						}
					}
					for (uint k = 0; k < trunIEN.size(); k++)
					{
						vector<int>::iterator it = find(tmesh[chdid].IEN.begin(), tmesh[chdid].IEN.end(), trunIEN[k]);
						if (it == tmesh[chdid].IEN.end())
						{
							tmesh[chdid].IEN.push_back(trunIEN[k]);
							array<double, 5> kutmp, kvtmp;
							for (int a = 0; a < 5; a++)
							{
								kutmp[a] = trun_patch_ku[k][a] - tmesh[rfid[i]].chd_o[j][0];
								kvtmp[a] = trun_patch_kv[k][a] - tmesh[rfid[i]].chd_o[j][1];
							}
							tmesh[chdid].patch_ku.push_back(kutmp);
							tmesh[chdid].patch_kv.push_back(kvtmp);
						}
					}
					//for(uint k=0; k<tmesh[chdid].IEN.size(); k++)
					//{
					//	cout<<tmesh[chdid].IEN[k]<<" ";
					//}
					//cout<<"\n";
					//getchar();
				}
			}
		}
		else if (tmesh[rfid[i]].type == 4 || tmesh[rfid[i]].type == 5)
		{
			for (int j = 0; j < 4; j++)
			{
				int chdid(tmesh[rfid[i]].chd[j]);
				if (chdid != -1)
				{
					tmesh[chdid].IEN.clear();
					tmesh[chdid].patch_ku.clear();
					tmesh[chdid].patch_kv.clear();
					tmesh[chdid].IEN = tmesh[chdid].IENtmp;
					tmesh[chdid].patch_ku = tmesh[chdid].patch_kutmp;
					tmesh[chdid].patch_kv = tmesh[chdid].patch_kvtmp;
					vector<int>().swap(tmesh[chdid].IENtmp);
					vector<array<double, 5>>().swap(tmesh[chdid].patch_kutmp);
					vector<array<double, 5>>().swap(tmesh[chdid].patch_kvtmp);
				}
			}
		}
	}

	//update old elements
	for (uint i = 0; i < nel_old; i++)
	{
		if (tmesh[i].act == 1 && (tmesh[i].type == 0 || tmesh[i].type == 1 || tmesh[i].type == 2))
		{
			vector<int> ien_old(tmesh[i].IEN);
			vector<array<double, 5>> ku_old(tmesh[i].patch_ku);
			vector<array<double, 5>> kv_old(tmesh[i].patch_kv);
			tmesh[i].IEN.clear();
			tmesh[i].patch_ku.clear();
			tmesh[i].patch_kv.clear();
			//
			for (uint j = 0; j < ien_old.size(); j++)
			{
				if (cp[ien_old[j]].trun == 1)
				{
					tmesh[i].IEN.push_back(ien_old[j]);
					tmesh[i].patch_ku.push_back(ku_old[j]);
					tmesh[i].patch_kv.push_back(kv_old[j]);
				}
			}
			//
			for (uint j = 0; j < tmesh[i].IENtmp.size(); j++)
			{
				//tmesh[i].IEN.push_back(tmesh[i].IENtmp[j]);
				if (cp[tmesh[i].IENtmp[j]].trun == 0)
				{
					tmesh[i].IEN.push_back(tmesh[i].IENtmp[j]);
					tmesh[i].patch_ku.push_back(tmesh[i].patch_kutmp[j]);
					tmesh[i].patch_kv.push_back(tmesh[i].patch_kvtmp[j]);
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
		}
		else if (tmesh[i].act == 1 && tmesh[i].type == 4)
		{
			vector<int> ien_old(tmesh[i].IEN);
			vector<array<double, 5>> ku_old(tmesh[i].patch_ku);
			vector<array<double, 5>> kv_old(tmesh[i].patch_kv);
			tmesh[i].IEN.clear();
			tmesh[i].patch_ku.clear();
			tmesh[i].patch_kv.clear();
			uint n_1r(2 * cp[tmesh[i].cnct[0]].face.size() + 1);
			for (uint j = 0; j < tmesh[i].IENtmp.size(); j++)
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
		}
	}
}

void TruncatedTspline::Topo_Refine_Unstruct_glb(const vector<int> &rfid, const vector<int> &rftype)
{
	//initialize for refinement
	npt_old = cp.size();
	nel_old = tmesh.size();

	for (uint i = 0; i < cp.size(); i++)
	{
		cp[i].update = 0; //initialized here!
		cp[i].aff = 0;
		cp[i].coortmp[0] = 0.;
		cp[i].coortmp[1] = 0.;
		cp[i].coortmp[2] = 0.;
		cp[i].wtmp = 0.;
		for (int j = 0; j < 3; j++)
		{
			cp[i].coor[j] *= cp[i].w;
		}
		for (int j = 0; j < 4; j++)
		{
			cp[i].kitvUtmp[j] = 0.;
			cp[i].kitvVtmp[j] = 0.;
		}
		cp[i].truntmp = 0;
		vector<int>().swap(cp[i].tbftmp);
		vector<double>().swap(cp[i].tctmp);
	}
	for (uint i = 0; i < tmesh.size(); i++)
	{
		//tmesh[i].aff=0;
		vector<int>().swap(tmesh[i].IENtmp);
		vector<array<double, 5>>().swap(tmesh[i].patch_kutmp);
		vector<array<double, 5>>().swap(tmesh[i].patch_kvtmp);
	}
	//refine
	for (uint i = 0; i < rfid.size(); i++)
	{
		if (tmesh[rfid[i]].type == 0)
		{
			ElementRefine_Unstruct_4(rfid[i]);
		}
		else if (tmesh[rfid[i]].type == 2)
		{
			ElementRefine_Unstruct_b(rfid[i]);
		}
		else if (tmesh[rfid[i]].type == 4)
		{
			ElementRefine_Irregular_4(rfid[i]);
		}
		//else
		//{
		//	cout<<"Not supported in global refinement!\n";
		//	getchar();
		//}
	}
}

void TruncatedTspline::Geom_Refine_Unstruct_glb(const vector<int> &rfid, const vector<int> &rftype)
{
	UpdateConnect_1();
	FindEdgeTopoDirec_1();
	FindKnotInterval_1();
	UpdateKnotInterval_1(); //might not work
	SetLocalCoorSystem();
	FindIEN_3();

	//calculate control points
	for (uint i = 0; i < rfid.size(); i++)
	{
		if (tmesh[rfid[i]].type == 0)
		{
			UpdatePatchCP_Unstruct_2(rfid[i]);
		}
	}
	//update coordinates
	for (uint i = 0; i < cp.size(); i++)
	{
		cp[i].coor[0] = cp[i].coortmp[0] / cp[i].wtmp;
		cp[i].coor[1] = cp[i].coortmp[1] / cp[i].wtmp;
		cp[i].coor[2] = cp[i].coortmp[2] / cp[i].wtmp;
		cp[i].w = cp[i].wtmp;
	}

	//update new elements
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && (tmesh[i].type == 0 || tmesh[i].type == 4))
		{
			tmesh[i].IEN.clear();
			tmesh[i].patch_ku.clear();
			tmesh[i].patch_kv.clear();
			tmesh[i].IEN = tmesh[i].IENtmp;
			tmesh[i].patch_ku = tmesh[i].patch_kutmp;
			tmesh[i].patch_kv = tmesh[i].patch_kvtmp;
			vector<int>().swap(tmesh[i].IENtmp);
			vector<array<double, 5>>().swap(tmesh[i].patch_kutmp);
			vector<array<double, 5>>().swap(tmesh[i].patch_kvtmp);
		}
	}
}

void TruncatedTspline::BezierExtract_Unstruct(vector<BezierElement> &bzmesh)
{
	bzmesh.clear();
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1)
		{
			if (tmesh[eid].type == 0 || tmesh[eid].type == 1)
			{
				BezierElementExtract_Unstruct(eid, bzmesh);
			}
			else if (tmesh[eid].type == 4)
			{
				BezierElementExtract_Unstruct_Irr(eid, bzmesh);
			}
		}
	}
}

void TruncatedTspline::BezierElementExtract_Unstruct(int eid, vector<BezierElement> &bzmesh)
{
	vector<array<double, 4>> be;
	BezierFinder_Unstruct(eid, be);
	uint nbzold = bzmesh.size();
	for (uint i = 0; i < be.size(); i++)
	{
		//BezierUnit_Unstruct(eid,be[i],bzmesh);
		BezierUnit_Unstruct_Trun(eid, be[i], bzmesh);
	}
	for (uint i = nbzold; i < bzmesh.size(); i++)
	{
		bzmesh[i].prt = eid;
	}
}

void TruncatedTspline::BezierFinder_Unstruct(int eid, vector<array<double, 4>> &be)
{
	be.clear();
	array<double, 4> e0 = {0., tmedge[tmesh[eid].edge[0]].len, 0., tmedge[tmesh[eid].edge[3]].len};
	vector<double> ukt, vkt;
	for (uint i = 0; i < tmesh[eid].patch_ku.size(); i++)
	{
		for (int j = 0; j < 5; j++)
		{
			if (tmesh[eid].patch_ku[i][j] > e0[0] && tmesh[eid].patch_ku[i][j] < e0[1])
			{
				vector<double>::iterator it = find(ukt.begin(), ukt.end(), tmesh[eid].patch_ku[i][j]);
				if (it == ukt.end())
					ukt.push_back(tmesh[eid].patch_ku[i][j]);
			}
			if (tmesh[eid].patch_kv[i][j] > e0[2] && tmesh[eid].patch_kv[i][j] < e0[3])
			{
				vector<double>::iterator it = find(vkt.begin(), vkt.end(), tmesh[eid].patch_kv[i][j]);
				if (it == vkt.end())
					vkt.push_back(tmesh[eid].patch_kv[i][j]);
			}
		}
	}
	if (ukt.size() != 0 && vkt.size() == 0)
	{
		sort(ukt.begin(), ukt.end());
		vector<double> ukt1;
		ukt1.push_back(e0[0]);
		for (uint i = 0; i < ukt.size(); i++)
		{
			ukt1.push_back(ukt[i]);
		}
		ukt1.push_back(e0[1]);
		for (uint i = 0; i < ukt1.size() - 1; i++)
		{
			array<double, 4> betmp = {ukt1[i], ukt1[i + 1], e0[2], e0[3]};
			be.push_back(betmp);
		}
	}
	else if (ukt.size() == 0 && vkt.size() != 0)
	{
		sort(vkt.begin(), vkt.end());
		vector<double> vkt1;
		vkt1.push_back(e0[2]);
		for (uint i = 0; i < vkt.size(); i++)
		{
			vkt1.push_back(vkt[i]);
		}
		vkt1.push_back(e0[3]);
		for (uint i = 0; i < vkt1.size() - 1; i++)
		{
			array<double, 4> betmp = {e0[0], e0[1], vkt1[i], vkt1[i + 1]};
			be.push_back(betmp);
		}
	}
	else if (ukt.size() == 0 && vkt.size() == 0)
	{
		be.push_back(e0);
	}
	else
	{
		cerr << "Not available for face intersection!\n";
		getchar();
	}
}

void TruncatedTspline::BezierUnit_Unstruct(int eid, array<double, 4> &kts, vector<BezierElement> &bzmesh)
{
	vector<int> pid = tmesh[eid].IEN;
	BezierElement bzel;
	bzel.cmat.resize(pid.size());
	for (uint i = 0; i < pid.size(); i++)
	{
		for (int j = 0; j < 16; j++)
		{
			bzel.cmat[i][j] = 0.;
		}
	}
	vector<vector<double>> coef(pid.size(), vector<double>(16));
	array<double, 2> ktsU = {kts[0], kts[1]};
	array<double, 2> ktsV = {kts[2], kts[3]};
	double bzku[6] = {ktsU[0], ktsU[0], ktsU[0], ktsU[1], ktsU[1], ktsU[1]};
	double bzkv[6] = {ktsV[0], ktsV[0], ktsV[0], ktsV[1], ktsV[1], ktsV[1]};
	for (uint i = 0; i < pid.size(); i++)
	{
		if (tmesh[eid].patch_ku[i][0] < ktsU[1] && tmesh[eid].patch_ku[i][4] > ktsU[0] && tmesh[eid].patch_kv[i][0] < ktsV[1] && tmesh[eid].patch_kv[i][4] > ktsV[0])
		{
			vector<double> ku(tmesh[eid].patch_ku[i].begin(), tmesh[eid].patch_ku[i].end());
			vector<double> kv(tmesh[eid].patch_kv[i].begin(), tmesh[eid].patch_kv[i].end());
			vector<double> ku1, kv1;
			vector<vector<double>> Tu, Tv;
			BezierInsertKnots(ku, ktsU, ku1);
			BezierInsertKnots(kv, ktsV, kv1);
			TMatrix(ku, ku1, 3, Tu);
			TMatrix(kv, kv1, 3, Tv);
			vector<double>::iterator it1 = search(ku1.begin(), ku1.end(), bzku, bzku + 6);
			vector<double>::iterator it2 = search(kv1.begin(), kv1.end(), bzkv, bzkv + 6);
			int loc1(it1 - ku1.begin() - 1), loc2(it2 - kv1.begin() - 1), count(0);
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					coef[i][count] = Tu[loc1 + k][0] * Tv[loc2 + j][0];
					count++;
				}
			}
		}
	}
	//vector<double> sum(16,0.);
	//for(uint i=0;i<pid.size();i++)
	//{
	//	for(int j=0; j<16; j++)
	//	{
	//		sum[j]+=coef[i][j];
	//	}
	//}
	//for(uint i=0;i<pid.size();i++)
	//{
	//	for(int j=0; j<16; j++)
	//	{
	//		coef[i][j]/=sum[j];
	//	}
	//}
	//cout<<"summation:\n";
	//cout.precision(16);
	//for(int j=0; j<16; j++)
	//{
	//	cout<<sum[j]<<"\n";
	//}
	//getchar();
	for (uint i = 0; i < pid.size(); i++)
	{
		if (tmesh[eid].patch_ku[i][0] < ktsU[1] && tmesh[eid].patch_ku[i][4] > ktsU[0] && tmesh[eid].patch_kv[i][0] < ktsV[1] && tmesh[eid].patch_kv[i][4] > ktsV[0])
		{
			for (int j = 0; j < 16; j++)
			{
				double tmp(coef[i][j]);
				//for(uint k=0;k<cp[pid[i]].tbf.size();k++)
				//{
				//	vector<int>::iterator it=find(pid.begin(),pid.end(),cp[pid[i]].tbf[k]);
				//	if(it!=pid.end())
				//	{
				//		int loc=it-pid.begin();
				//		tmp-=cp[pid[i]].tc[k]*coef[loc][j];
				//	}
				//}
				bzel.pts[j][0] += tmp * cp[pid[i]].coor[0];
				bzel.pts[j][1] += tmp * cp[pid[i]].coor[1];
				bzel.pts[j][2] += tmp * cp[pid[i]].coor[2];
				bzel.cmat[i][j] = tmp;
			}
		}
	}
	//vector<vector<double>> demat;
	//DegreeElevate(demat);
	//bzel.cmat4.resize(tmesh[eid].IEN.size());
	//for(uint i=0; i<tmesh[eid].IEN.size(); i++)
	//{
	//	for(int j=0; j<25; j++)
	//	{
	//		bzel.cmat4[i][j]=0.;
	//	}
	//}
	//for(uint i=0; i<tmesh[eid].IEN.size(); i++)
	//{
	//	for(int j=0; j<16; j++)
	//	{
	//		for(int k=0; k<25; k++)
	//		{
	//			bzel.cmat4[i][k]+=bzel.cmat[i][j]*demat[k][j];
	//		}
	//	}
	//}
	//for(uint i=0;i<pid.size();i++)
	//{
	//	for(int j=0;j<25;j++)
	//	{
	//		bzel.pts4[j][0]+=bzel.cmat4[i][j]*cp[pid[i]].coor[0];
	//		bzel.pts4[j][1]+=bzel.cmat4[i][j]*cp[pid[i]].coor[1];
	//		bzel.pts4[j][2]+=bzel.cmat4[i][j]*cp[pid[i]].coor[2];
	//	}
	//}
	for (uint i = 0; i < pid.size(); i++)
	{
		bzel.IEN.push_back(pid[i]);
	}
	bzmesh.push_back(bzel);
}

void TruncatedTspline::BezierUnit_Unstruct_Irr(int eid, array<double, 4> &kts, vector<BezierElement> &bzmesh)
{
	vector<int> pid = tmesh[eid].IEN;
	BezierElement bzel;
	bzel.order = 4;
	bzel.cmat4.resize(pid.size());
	bzel.cmat.resize(pid.size());
	for (uint i = 0; i < pid.size(); i++)
	{
		for (int j = 0; j < 25; j++)
		{
			bzel.cmat4[i][j] = 0.;
		}
	}
	for (uint i = 0; i < pid.size(); i++)
	{
		for (int j = 0; j < 16; j++)
		{
			bzel.cmat[i][j] = 0.;
		}
	}
	//first 2N+1
	uint nv(cp[tmesh[eid].cnct[0]].face.size());
	vector<vector<double>> bmat, bmat3;
	//SetBezier4TranMatOP(nv,bmat);
	SetBezier4TranMat(nv, bmat);
	SetBezier3TranMat(nv, bmat3);
	int dir(-1);
	double insU(kts[0]), insV(kts[2]);
	//double urang[2]={0.,1.}, vrang[2]={0.,1.};
	double urang[2] = {0., tmedge[tmesh[eid].edge[0]].len}, vrang[2] = {0., tmedge[tmesh[eid].edge[3]].len};
	int pos(0);
	for (int i = 0; i < 2; i++)
	{
		if (kts[i] > urang[0] && kts[i] < urang[1])
		{
			if (i == 0)
				pos = 4;
			dir = 0;
			insU = kts[i];
			break;
		}
		if (kts[i + 2] > vrang[0] && kts[i + 2] < vrang[1])
		{
			if (i == 0)
				pos = 4;
			dir = 1;
			insV = kts[i + 2];
			break;
		}
	}
	vector<double> knv0, knv1;
	if (dir == 0)
	{
		double bzku_0[10] = {urang[0], urang[0], urang[0], urang[0], urang[0], urang[1], urang[1], urang[1], urang[1], urang[1]};
		double bzku_1[14] = {urang[0], urang[0], urang[0], urang[0], urang[0], insU, insU, insU, insU, urang[1], urang[1], urang[1], urang[1], urang[1]};
		knv0.assign(bzku_0, bzku_0 + 10);
		knv1.assign(bzku_1, bzku_1 + 14);
	}
	else if (dir == 1)
	{
		double bzkv_0[10] = {vrang[0], vrang[0], vrang[0], vrang[0], vrang[0], vrang[1], vrang[1], vrang[1], vrang[1], vrang[1]};
		double bzkv_1[14] = {vrang[0], vrang[0], vrang[0], vrang[0], vrang[0], insV, insV, insV, insV, vrang[1], vrang[1], vrang[1], vrang[1], vrang[1]};
		knv0.assign(bzkv_0, bzkv_0 + 10);
		knv1.assign(bzkv_1, bzkv_1 + 14);
	}
	if (dir == 0 || dir == 1)
	{
		vector<vector<double>> tmat0, tmat_u, tmat_v;
		TMatrix(knv0, knv1, 4, tmat0);
		tmat_u.resize(5, vector<double>(5, 0.));
		tmat_v.resize(5, vector<double>(5, 0.));
		if (dir == 0)
		{
			for (uint i = 0; i < 5; i++)
			{
				tmat_v[i][i] = 1.;
				for (uint j = 0; j < 5; j++)
				{
					tmat_u[i][j] = tmat0[pos + i][j];
				}
			}
		}
		else
		{
			for (uint i = 0; i < 5; i++)
			{
				tmat_u[i][i] = 1.;
				for (uint j = 0; j < 5; j++)
				{
					tmat_v[i][j] = tmat0[pos + i][j];
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
		for (uint i = 0; i < 2 * nv + 1; i++)
		{
			for (int j = 0; j < 25; j++)
			{
				for (int k = 0; k < 25; k++)
				{
					int iold(j / 5), jold(j % 5), knew(k / 5), lnew(k % 5);
					bzel.cmat4[i][k] += bmat[i][j] * tmat_u[lnew][jold] * tmat_v[knew][iold];
					//bzel.cmat4[i][k]+=bmat[i][j]*Tmat[j][k];
				}
			}
		}
	}
	else
	{
		for (uint i = 0; i < 2 * nv + 1; i++)
		{
			for (int j = 0; j < 25; j++)
			{
				bzel.cmat4[i][j] = bmat[i][j];
			}
		}
		for (uint i = 0; i < 2 * nv + 1; i++)
		{
			for (int j = 0; j < 16; j++)
			{
				bzel.cmat[i][j] = bmat3[i][j];
			}
		}
	}

	//the others
	vector<vector<double>> coef(tmesh[eid].patch_ku.size(), vector<double>(16, 0.));
	array<double, 2> ktsU = {kts[0], kts[1]};
	array<double, 2> ktsV = {kts[2], kts[3]};
	double bzku[6] = {ktsU[0], ktsU[0], ktsU[0], ktsU[1], ktsU[1], ktsU[1]};
	double bzkv[6] = {ktsV[0], ktsV[0], ktsV[0], ktsV[1], ktsV[1], ktsV[1]};
	for (uint i = 0; i < tmesh[eid].patch_ku.size(); i++)
	{
		if (tmesh[eid].patch_ku[i][0] < ktsU[1] && tmesh[eid].patch_ku[i][4] > ktsU[0] && tmesh[eid].patch_kv[i][0] < ktsV[1] && tmesh[eid].patch_kv[i][4] > ktsV[0])
		{
			vector<double> ku(tmesh[eid].patch_ku[i].begin(), tmesh[eid].patch_ku[i].end());
			vector<double> kv(tmesh[eid].patch_kv[i].begin(), tmesh[eid].patch_kv[i].end());
			vector<double> ku1, kv1;
			vector<vector<double>> Tu, Tv;
			BezierInsertKnots(ku, ktsU, ku1);
			BezierInsertKnots(kv, ktsV, kv1);
			TMatrix(ku, ku1, 3, Tu);
			TMatrix(kv, kv1, 3, Tv);
			vector<double>::iterator it1 = search(ku1.begin(), ku1.end(), bzku, bzku + 6);
			vector<double>::iterator it2 = search(kv1.begin(), kv1.end(), bzkv, bzkv + 6);
			int loc1(it1 - ku1.begin() - 1), loc2(it2 - kv1.begin() - 1), count(0);
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					coef[i][count] = Tu[loc1 + k][0] * Tv[loc2 + j][0];
					count++;
				}
			}
		}
	}
	for (uint i = 0; i < tmesh[eid].patch_ku.size(); i++)
	{
		for (int j = 0; j < 16; j++)
		{
			bzel.cmat[i + (2 * nv + 1)][j] = coef[i][j];
		}
	}
	vector<vector<double>> demat;
	DegreeElevate(demat);
	for (uint i = 0; i < tmesh[eid].patch_ku.size(); i++)
	{
		for (int j = 0; j < 16; j++)
		{
			for (int k = 0; k < 25; k++)
			{
				bzel.cmat4[i + (2 * nv + 1)][k] += coef[i][j] * demat[k][j];
			}
		}
	}
	//vector<double> sum(25,0.);
	//for(uint i=0;i<pid.size();i++)
	//{
	//	for(int j=0; j<25; j++)
	//	{
	//		sum[j]+=bzel.cmat4[i][j];
	//	}
	//}
	//for(uint i=0;i<pid.size();i++)
	//{
	//	for(int j=0; j<25; j++)
	//	{
	//		bzel.cmat4[i][j]/=sum[j];
	//	}
	//}
	//cout<<"summation:\n";
	//cout.precision(16);
	//for(int j=0; j<25; j++)
	//{
	//	cout<<sum[j]<<"\n";
	//}
	//getchar();
	for (uint i = 0; i < pid.size(); i++)
	{
		for (int j = 0; j < 25; j++)
		{
			bzel.pts4[j][0] += bzel.cmat4[i][j] * cp[pid[i]].coor[0];
			bzel.pts4[j][1] += bzel.cmat4[i][j] * cp[pid[i]].coor[1];
			bzel.pts4[j][2] += bzel.cmat4[i][j] * cp[pid[i]].coor[2];
		}
	}
	for (uint i = 0; i < pid.size(); i++)
	{
		for (int j = 0; j < 16; j++)
		{
			bzel.pts[j][0] += bzel.cmat[i][j] * cp[pid[i]].coor[0];
			bzel.pts[j][1] += bzel.cmat[i][j] * cp[pid[i]].coor[1];
			bzel.pts[j][2] += bzel.cmat[i][j] * cp[pid[i]].coor[2];
		}
	}
	for (uint i = 0; i < pid.size(); i++)
	{
		bzel.IEN.push_back(pid[i]);
	}
	bzmesh.push_back(bzel);
}

void TruncatedTspline::BezierUnit_Unstruct_Trun(int eid, array<double, 4> &kts, vector<BezierElement> &bzmesh)
{
	vector<int> pid = tmesh[eid].IEN;
	BezierElement bzel;
	bzel.cmat.resize(pid.size());
	for (uint i = 0; i < pid.size(); i++)
	{
		for (int j = 0; j < 16; j++)
		{
			bzel.cmat[i][j] = 0.;
		}
	}
	vector<vector<double>> coef(pid.size(), vector<double>(16));
	array<double, 2> ktsU = {kts[0], kts[1]};
	array<double, 2> ktsV = {kts[2], kts[3]};
	double bzku[6] = {ktsU[0], ktsU[0], ktsU[0], ktsU[1], ktsU[1], ktsU[1]};
	double bzkv[6] = {ktsV[0], ktsV[0], ktsV[0], ktsV[1], ktsV[1], ktsV[1]};
	for (uint i = 0; i < pid.size(); i++)
	{
		if (tmesh[eid].patch_ku[i][0] < ktsU[1] && tmesh[eid].patch_ku[i][4] > ktsU[0] && tmesh[eid].patch_kv[i][0] < ktsV[1] && tmesh[eid].patch_kv[i][4] > ktsV[0])
		{
			vector<double> ku(tmesh[eid].patch_ku[i].begin(), tmesh[eid].patch_ku[i].end());
			vector<double> kv(tmesh[eid].patch_kv[i].begin(), tmesh[eid].patch_kv[i].end());
			vector<double> ku1, kv1;
			vector<vector<double>> Tu, Tv;
			BezierInsertKnots(ku, ktsU, ku1);
			BezierInsertKnots(kv, ktsV, kv1);
			TMatrix(ku, ku1, 3, Tu);
			TMatrix(kv, kv1, 3, Tv);
			vector<double>::iterator it1 = search(ku1.begin(), ku1.end(), bzku, bzku + 6);
			vector<double>::iterator it2 = search(kv1.begin(), kv1.end(), bzkv, bzkv + 6);
			int loc1(it1 - ku1.begin() - 1), loc2(it2 - kv1.begin() - 1), count(0);
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					coef[i][count] = Tu[loc1 + k][0] * Tv[loc2 + j][0];
					count++;
				}
			}
		}
	}
	for (uint i = 0; i < pid.size(); i++)
	{
		if (tmesh[eid].patch_ku[i][0] < ktsU[1] && tmesh[eid].patch_ku[i][4] > ktsU[0] && tmesh[eid].patch_kv[i][0] < ktsV[1] && tmesh[eid].patch_kv[i][4] > ktsV[0])
		{
			for (int j = 0; j < 16; j++)
			{
				double tmp(coef[i][j]);
				for (uint k = 0; k < cp[pid[i]].tbf.size(); k++)
				{
					vector<int>::iterator it = find(pid.begin(), pid.end(), cp[pid[i]].tbf[k]);
					if (it != pid.end())
					{
						int loc = it - pid.begin();
						tmp -= cp[pid[i]].tc[k] * coef[loc][j];
					}
				}
				bzel.pts[j][0] += tmp * cp[pid[i]].coor[0];
				bzel.pts[j][1] += tmp * cp[pid[i]].coor[1];
				bzel.pts[j][2] += tmp * cp[pid[i]].coor[2];
				bzel.cmat[i][j] = tmp;
			}
		}
	}
	//vector<vector<double>> demat;
	//DegreeElevate(demat);
	//bzel.cmat4.resize(tmesh[eid].IEN.size());
	//for(uint i=0; i<tmesh[eid].IEN.size(); i++)
	//{
	//	for(int j=0; j<25; j++)
	//	{
	//		bzel.cmat4[i][j]=0.;
	//	}
	//}
	//for(uint i=0; i<tmesh[eid].IEN.size(); i++)
	//{
	//	for(int j=0; j<16; j++)
	//	{
	//		for(int k=0; k<25; k++)
	//		{
	//			bzel.cmat4[i][k]+=bzel.cmat[i][j]*demat[k][j];
	//		}
	//	}
	//}
	//for(uint i=0;i<pid.size();i++)
	//{
	//	for(int j=0;j<25;j++)
	//	{
	//		bzel.pts4[j][0]+=bzel.cmat4[i][j]*cp[pid[i]].coor[0];
	//		bzel.pts4[j][1]+=bzel.cmat4[i][j]*cp[pid[i]].coor[1];
	//		bzel.pts4[j][2]+=bzel.cmat4[i][j]*cp[pid[i]].coor[2];
	//	}
	//}
	for (uint i = 0; i < pid.size(); i++)
	{
		bzel.IEN.push_back(pid[i]);
	}
	bzmesh.push_back(bzel);
}

void TruncatedTspline::BezierUnit_Unstruct_Irr_Trun(int eid, array<double, 4> &kts, vector<BezierElement> &bzmesh)
{
	vector<int> pid = tmesh[eid].IEN;
	BezierElement bzel;
	bzel.order = 4;
	bzel.cmat4.resize(pid.size());
	//bzel.cmat.resize(pid.size());
	for (uint i = 0; i < pid.size(); i++)
	{
		for (int j = 0; j < 25; j++)
		{
			bzel.cmat4[i][j] = 0.;
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
	vector<vector<double>> cmat4(pid.size(), vector<double>(25));
	int dir(-1);
	double insU(kts[0]), insV(kts[2]);
	//double urang[2]={0.,1.}, vrang[2]={0.,1.};
	double urang[2] = {0., tmedge[tmesh[eid].edge[0]].len}, vrang[2] = {0., tmedge[tmesh[eid].edge[3]].len};
	int pos(0);
	for (int i = 0; i < 2; i++)
	{
		if (kts[i] > urang[0] && kts[i] < urang[1])
		{
			if (i == 0)
				pos = 4;
			dir = 0;
			insU = kts[i];
			break;
		}
		if (kts[i + 2] > vrang[0] && kts[i + 2] < vrang[1])
		{
			if (i == 0)
				pos = 4;
			dir = 1;
			insV = kts[i + 2];
			break;
		}
	}
	vector<double> knv0, knv1;
	if (dir == 0)
	{
		double bzku_0[10] = {urang[0], urang[0], urang[0], urang[0], urang[0], urang[1], urang[1], urang[1], urang[1], urang[1]};
		double bzku_1[14] = {urang[0], urang[0], urang[0], urang[0], urang[0], insU, insU, insU, insU, urang[1], urang[1], urang[1], urang[1], urang[1]};
		knv0.assign(bzku_0, bzku_0 + 10);
		knv1.assign(bzku_1, bzku_1 + 14);
	}
	else if (dir == 1)
	{
		double bzkv_0[10] = {vrang[0], vrang[0], vrang[0], vrang[0], vrang[0], vrang[1], vrang[1], vrang[1], vrang[1], vrang[1]};
		double bzkv_1[14] = {vrang[0], vrang[0], vrang[0], vrang[0], vrang[0], insV, insV, insV, insV, vrang[1], vrang[1], vrang[1], vrang[1], vrang[1]};
		knv0.assign(bzkv_0, bzkv_0 + 10);
		knv1.assign(bzkv_1, bzkv_1 + 14);
	}
	if (dir == 0 || dir == 1)
	{
		vector<vector<double>> tmat0, tmat_u, tmat_v;
		TMatrix(knv0, knv1, 4, tmat0);
		tmat_u.resize(5, vector<double>(5, 0.));
		tmat_v.resize(5, vector<double>(5, 0.));
		if (dir == 0)
		{
			for (uint i = 0; i < 5; i++)
			{
				tmat_v[i][i] = 1.;
				for (uint j = 0; j < 5; j++)
				{
					tmat_u[i][j] = tmat0[pos + i][j];
				}
			}
		}
		else
		{
			for (uint i = 0; i < 5; i++)
			{
				tmat_u[i][i] = 1.;
				for (uint j = 0; j < 5; j++)
				{
					tmat_v[i][j] = tmat0[pos + i][j];
				}
			}
		}
		for (uint i = 0; i < 2 * nv + 1; i++)
		{
			for (int j = 0; j < 25; j++)
			{
				for (int k = 0; k < 25; k++)
				{
					int iold(j / 5), jold(j % 5), knew(k / 5), lnew(k % 5);
					cmat4[i][k] += tmesh[eid].bemat[i][j] * tmat_u[lnew][jold] * tmat_v[knew][iold];
				}
			}
		}
	}
	else
	{
		for (uint i = 0; i < 2 * nv + 1; i++)
		{
			for (int j = 0; j < 25; j++)
			{
				cmat4[i][j] = tmesh[eid].bemat[i][j];
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
	vector<vector<double>> coef(tmesh[eid].patch_ku.size(), vector<double>(16, 0.));
	array<double, 2> ktsU = {kts[0], kts[1]};
	array<double, 2> ktsV = {kts[2], kts[3]};
	double bzku[6] = {ktsU[0], ktsU[0], ktsU[0], ktsU[1], ktsU[1], ktsU[1]};
	double bzkv[6] = {ktsV[0], ktsV[0], ktsV[0], ktsV[1], ktsV[1], ktsV[1]};
	for (uint i = 0; i < tmesh[eid].patch_ku.size(); i++)
	{
		if (tmesh[eid].patch_ku[i][0] < ktsU[1] && tmesh[eid].patch_ku[i][4] > ktsU[0] && tmesh[eid].patch_kv[i][0] < ktsV[1] && tmesh[eid].patch_kv[i][4] > ktsV[0])
		{
			vector<double> ku(tmesh[eid].patch_ku[i].begin(), tmesh[eid].patch_ku[i].end());
			vector<double> kv(tmesh[eid].patch_kv[i].begin(), tmesh[eid].patch_kv[i].end());
			vector<double> ku1, kv1;
			vector<vector<double>> Tu, Tv;
			BezierInsertKnots(ku, ktsU, ku1);
			BezierInsertKnots(kv, ktsV, kv1);
			TMatrix(ku, ku1, 3, Tu);
			TMatrix(kv, kv1, 3, Tv);
			vector<double>::iterator it1 = search(ku1.begin(), ku1.end(), bzku, bzku + 6);
			vector<double>::iterator it2 = search(kv1.begin(), kv1.end(), bzkv, bzkv + 6);
			int loc1(it1 - ku1.begin() - 1), loc2(it2 - kv1.begin() - 1), count(0);
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					coef[i][count] = Tu[loc1 + k][0] * Tv[loc2 + j][0];
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
	for (uint i = 0; i < tmesh[eid].patch_ku.size(); i++)
	{
		for (int j = 0; j < 16; j++)
		{
			for (int k = 0; k < 25; k++)
			{
				cmat4[i + (2 * nv + 1)][k] += coef[i][j] * demat[k][j];
			}
		}
	}
	//truncation
	for (uint i = 0; i < pid.size(); i++)
	{
		for (int j = 0; j < 25; j++)
		{
			bzel.cmat4[i][j] = cmat4[i][j];
			if (cp[pid[i]].trun == 1)
			{
				for (uint k = 0; k < cp[pid[i]].tbf.size(); k++)
				{
					vector<int>::iterator it = find(pid.begin(), pid.end(), cp[pid[i]].tbf[k]);
					if (it != pid.end())
					{
						int loc = it - pid.begin();
						bzel.cmat4[i][j] -= cp[pid[i]].tc[k] * cmat4[loc][j];
					}
				}
			}
		}
	}
	for (uint i = 0; i < pid.size(); i++)
	{
		for (int j = 0; j < 25; j++)
		{
			bzel.pts4[j][0] += bzel.cmat4[i][j] * cp[pid[i]].coor[0];
			bzel.pts4[j][1] += bzel.cmat4[i][j] * cp[pid[i]].coor[1];
			bzel.pts4[j][2] += bzel.cmat4[i][j] * cp[pid[i]].coor[2];
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
	for (uint i = 0; i < pid.size(); i++)
	{
		bzel.IEN.push_back(pid[i]);
	}
	bzmesh.push_back(bzel);
}

void TruncatedTspline::BezierElementExtract_Unstruct_Irr(int eid, vector<BezierElement> &bzmesh)
{
	vector<array<double, 4>> be;
	BezierFinder_Unstruct(eid, be);
	uint nbzold = bzmesh.size();
	for (uint i = 0; i < be.size(); i++)
	{
		//BezierUnit_Unstruct_Irr(eid,be[i],bzmesh);
		BezierUnit_Unstruct_Irr_Trun(eid, be[i], bzmesh);
	}
	for (uint i = nbzold; i < bzmesh.size(); i++)
	{
		bzmesh[i].prt = eid;
	}
}

void TruncatedTspline::BezierVTK_Unstruct(string fn, vector<BezierElement> &bzmesh)
{
	vector<array<double, 3>> spt;
	vector<double> sval;
	vector<array<double, 3>> norm;
	vector<array<int, 4>> sele;
	vector<array<double, 3>> lpt; //visulize parameter lines
	vector<array<int, 2>> led;	  //line connectivity
	int ns(5), ecount(0);
	vector<double> su(ns), sv(ns);
	for (int i = 0; i < ns; i++)
	{
		su[i] = double(i) / double(ns - 1);
		sv[i] = double(i) / double(ns - 1);
	}

	for (uint e = 0; e < bzmesh.size(); e++)
	{
		int loc(0);
		for (int a = 0; a < ns; a++)
		{
			for (int b = 0; b < ns; b++)
			{
				array<double, 3> pt, nm;
				/*double ptmp[3];
				if(bzmesh[e].order==3)
					bzmesh[e].Para2Phys(su[b],sv[a],ptmp);
				else
					bzmesh[e].Para2Phys4(su[b],sv[a],ptmp);*/
				//array<double,3> ptmp,nmtmp;
				if (bzmesh[e].order == 3)
					bzmesh[e].SurfPointNormal(su[b], sv[a], pt, nm);
				else
					bzmesh[e].SurfPointNormal4(su[b], sv[a], pt, nm);
				//pt.coor[0]=ptmp[0]; pt.coor[1]=ptmp[1]; pt.coor[2]=ptmp[2];
				spt.push_back(pt);
				norm.push_back(nm);
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
				el[0] = ecount * ns * ns + a * ns + b;
				el[1] = ecount * ns * ns + a * ns + b + 1;
				el[2] = ecount * ns * ns + (a + 1) * ns + b + 1;
				el[3] = ecount * ns * ns + (a + 1) * ns + b;
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

	string fname = fn + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i < spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
		}
		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i]<<"\n";
		//}
		fout << "\nPOINT_DATA " << norm.size() << "\nNORMALS Normal FLOAT\n";
		for (uint i = 0; i < norm.size(); i++)
		{
			fout << norm[i][0] << " " << norm[i][1] << " " << norm[i][2] << "\n";
		}
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

void TruncatedTspline::BezierControlMesh_Unstruct(string fn, vector<BezierElement> &bzmesh)
{
	vector<array<double, 3>> spt;
	vector<double> sval;
	vector<array<int, 4>> sele;
	//int ns(2);
	int ecount(0);
	//vector<double> su(ns),sv(ns);
	//for(int i=0;i<ns;i++)
	//{
	//	su[i]=double(i)/double(ns-1);
	//	sv[i]=double(i)/double(ns-1);
	//}

	for (uint e = 0; e < bzmesh.size(); e++)
	{
		int loc(0);
		int ns = bzmesh[e].order + 1;
		int nold = spt.size();
		for (int a = 0; a < ns; a++)
		{
			for (int b = 0; b < ns; b++)
			{
				array<double, 3> pt;
				if (bzmesh[e].order == 3)
				{
					pt[0] = bzmesh[e].pts[loc][0];
					pt[1] = bzmesh[e].pts[loc][1];
					pt[2] = bzmesh[e].pts[loc][2];
				}
				else if (bzmesh[e].order == 4)
				{
					pt[0] = bzmesh[e].pts4[loc][0];
					pt[1] = bzmesh[e].pts4[loc][1];
					pt[2] = bzmesh[e].pts4[loc][2];
				}
				spt.push_back(pt);
				loc++;
			}
		}
		for (int a = 0; a < ns - 1; a++)
		{
			for (int b = 0; b < ns - 1; b++)
			{
				array<int, 4> el;
				el[0] = nold + a * ns + b;
				el[1] = nold + a * ns + b + 1;
				el[2] = nold + (a + 1) * ns + b + 1;
				el[3] = nold + (a + 1) * ns + b;
				sele.push_back(el);
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
		for (uint i = 0; i < spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
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
		cout << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline::run_surf(string fn)
{
	SetLshapeProblem(fn);
	InitialTopology();
	//VisualizeTMesh("Lshape_new/CM_0");
	//VisualizeControlMesh("Lshape_new/CM_0");
	//VisualizeVTK("Lshape_new/surf_0");
	//CollectActives();
	SurfRefine();
	//VisualizeTMesh("Lshape_new/CM_1");
	VisualizeControlMesh(fn + "_CM_0");
	VisualizeVTK(fn + "surf_0");
	//cout<<"DOF: "<<cp.size()<<"\n";
}

void TruncatedTspline::run_surf_1(string fn)
{
	SetProblem_surf(fn);
	InitialTopology();
	//CollectActives();
	SurfRefine_1();
	VisualizeTMesh(fn + "_TM_0");
	//VisualizeControlMesh(fn+"_CM_0");
	VisualizeVTK(fn + "surf_0");
	//cout<<"DOF: "<<cp.size()<<"\n";
}

void TruncatedTspline::run_Laplace(string fn)
{
	SetLshapeProblem(fn);
	InitialTopology();
	//VisualizeTMesh("Lshape_new/CM_0");
	//VisualizeControlMesh("Lshape_new/CM_0");
	//VisualizeVTK("Lshape_new/surf_0");
	//CollectActives();
	//LshapeRefine();
	LshapeRefine_H1();
	//VisualizeTMesh("Lshape_new/CM_1");
	VisualizeControlMesh("Lshape_new/CM_5");
	VisualizeVTK("Lshape_new/surf_5");
	cout << "DOF: " << cp.size() << "\n";
}

void TruncatedTspline::runXP()
{
	//SetProblem2("test7/XP_2_input_pillow.vtk");
	//SetProblem2("test7/XP_0_input.vtk");
	//SetProblem2("test7/squarexp4s_quad.vtk");
	SetProblem2("../io/input/XP_2_input_plane.vtk");
	InitialConnect();
	////VisualizeTMesh("test8/Tmesh_0");
	////ElementRefine_Unstruct_2(25,0);
	////double coortmp[4][3]={{cp[12].coor[0],cp[12].coor[1],cp[12].coor[2]},{cp[1].coor[0],cp[1].coor[1],cp[1].coor[2]},{cp[36].coor[0],cp[36].coor[1],cp[36].coor[2]},{cp[37].coor[0],cp[37].coor[1],cp[37].coor[2]}};
	////cp[12].coor[0]=(5.*coortmp[0][0]+cp[0].coor[0])/6.; cp[12].coor[1]=(5.*coortmp[0][1]+cp[0].coor[1])/6.; cp[12].coor[2]=(5.*coortmp[0][2]+cp[0].coor[2])/6.;
	////cp[1].coor[0]=(5.*coortmp[1][0]+cp[17].coor[0])/6.; cp[1].coor[1]=(5.*coortmp[1][1]+cp[17].coor[1])/6.; cp[1].coor[2]=(5.*coortmp[1][2]+cp[17].coor[2])/6.;
	////cp[36].coor[0]=(5.*coortmp[2][0]+cp[35].coor[0])/6.; cp[36].coor[1]=(5.*coortmp[2][1]+cp[35].coor[1])/6.; cp[36].coor[2]=(5.*coortmp[2][2]+cp[35].coor[2])/6.;
	////cp[37].coor[0]=(5.*coortmp[3][0]+cp[38].coor[0])/6.; cp[37].coor[1]=(5.*coortmp[3][1]+cp[38].coor[1])/6.; cp[37].coor[2]=(5.*coortmp[3][2]+cp[38].coor[2])/6.;
	////UpdateConnect();
	//int rf_id0[2]={25,26};
	//int rf_type0[2]={0,0};
	//vector<int> rf_id(rf_id0,rf_id0+1);
	//vector<int> rf_type(rf_type0,rf_type0+1);
	//Refine_Unstruct(rf_id,rf_type);
	SetBezierMatIrrPatch();
	CollectActives();
	//VisualizeTMesh("test8/Tmesh_3");
	VisualizeSurface("../io/surf1/XP2_1");

	//vector<BezierElement> bzmesh;
	//BezierExtract_Unstruct(bzmesh);
	////BezierControlMesh_Unstruct("test8/bzmesh_2.vtk",bzmesh);
	//BezierVTK_Unstruct("test8/bzmesh_2.vtk",bzmesh);
}

void TruncatedTspline::runXP_Laplace()
{
	SetLshapeProblem_XP("Lshape_new/input_CM_XP");
	//InitialConnect();
	//int rf_id0[2]={25,26};
	//int rf_type0[2]={0,0};
	//vector<int> rf_id(rf_id0,rf_id0+1);
	//vector<int> rf_type(rf_type0,rf_type0+1);
	//Refine_Unstruct(rf_id,rf_type);
	LshapeRefine_XP();
	SetBezierMatIrrPatch();
	CollectActives();
	cout << "DOF: " << cp.size() << "\n";
	//VisualizeTMesh("Lshape_new/Tmesh_XP_1");
	//VisualizeControlMesh("Lshape_new/CM_XP_4");
	VisualizeSurface("Lshape_new/surf_XP_4");

	//vector<BezierElement> bzmesh;
	//BezierExtract_Unstruct(bzmesh);
	//BezierControlMesh_Unstruct("test8/Lshape_bzmesh_CM_2.vtk",bzmesh);
	//BezierVTK_Unstruct("test8/Lshape_bzmesh_2.vtk",bzmesh);
}

void TruncatedTspline::run_PlateHole(string fn)
{
	SetProblem_PlateHole(fn);
	InitialTopology();
	//VisualizeTMesh("Lshape_new/CM_0");
	//VisualizeControlMesh("Lshape_new/CM_0");
	for (int i = 0; i < 3; i++)
		Refine_Global();
	InitialTopology_PlateHole();
	Refine_PlateHole();
	VisualizeControlMesh("plate_hole/CM_adapt_3");
	//VisualizeVTK("plate_hole/surf_3");
	cout << "done with surf\n";
	getchar();
	CollectActives();
	cout << "DOF: " << cp.size() << "\n";
}

void TruncatedTspline::Refine_Global()
{
	npt_old = cp.size();
	nel_old = tmesh.size();

	vector<double> kts_u, kts_v, kts_u1, kts_v1;
	for (uint i = 0; i < 2; i++)
	{
		kts_u.push_back(0.);
		kts_v.push_back(0.);
		kts_u1.push_back(0.);
		kts_v1.push_back(0.);
	}
	for (uint i = 0; i < uanc.size(); i++)
	{
		kts_u.push_back(uanc[i].second);
	}
	for (uint i = 0; i < vanc.size(); i++)
	{
		kts_v.push_back(vanc[i].second);
	}
	for (uint i = 0; i < 2; i++)
	{
		kts_u.push_back(1.);
		kts_v.push_back(1.);
	}
	for (uint i = 0; i < uanc.size() - 1; i++)
	{
		kts_u1.push_back(uanc[i].second);
		if (uanc[i].second < uanc[i + 1].second)
		{
			kts_u1.push_back((uanc[i].second + uanc[i + 1].second) / 2.);
		}
	}
	for (uint i = 0; i < vanc.size() - 1; i++)
	{
		kts_v1.push_back(vanc[i].second);
		if (vanc[i].second < vanc[i + 1].second)
		{
			kts_v1.push_back((vanc[i].second + vanc[i + 1].second) / 2.);
		}
	}
	for (uint i = 0; i < 3; i++)
	{
		kts_u1.push_back(1.);
		kts_v1.push_back(1.);
	}

	vector<vector<double>> Tu, Tv;
	TMatrix(kts_u, kts_u1, 3, Tu);
	TMatrix(kts_v, kts_v1, 3, Tv);

	uint npu(kts_u1.size() - 4), npv(kts_v1.size() - 4), ncp(npu * npv), nfc((npu - 1) * (npv - 1));
	//cout<<npu<<"\n";
	//Tall.resize(ncp,vector<double>(npt_old,0.));
	for (uint i = 0; i < cp.size(); i++)
	{
		cp[i].coortmp[0] = cp[i].coor[0] * cp[i].w;
		cp[i].coortmp[1] = cp[i].coor[1] * cp[i].w;
		cp[i].coortmp[2] = cp[i].coor[2] * cp[i].w;
	}
	vector<array<double, 4>> pts_new(ncp);
	int loc(0), loc1(0);
	for (uint j = 0; j < npv; j++)
	{
		for (uint i = 0; i < npu; i++)
		{
			pts_new[loc][0] = 0.;
			pts_new[loc][1] = 0.;
			pts_new[loc][2] = 0.;
			pts_new[loc][3] = 0.;
			loc1 = 0;
			for (uint j1 = 0; j1 < vanc.size(); j1++)
			{
				for (uint i1 = 0; i1 < uanc.size(); i1++)
				{
					if (Tu[i][i1] != 0. && Tv[j][j1] != 0.)
					{
						double tmp(Tu[i][i1] * Tv[j][j1]);
						for (int dof = 0; dof < 3; dof++)
							pts_new[loc][dof] += tmp * cp[loc1].coortmp[dof];
						pts_new[loc][3] += tmp * cp[loc1].w;
					}
					loc1++;
				}
			}
			pts_new[loc][0] /= pts_new[loc][3];
			pts_new[loc][1] /= pts_new[loc][3];
			pts_new[loc][2] /= pts_new[loc][3];
			loc++;
		}
	}

	vector<Vertex>().swap(cp);
	vector<Element>().swap(tmesh);
	cp.resize(ncp);
	tmesh.resize(nfc);
	loc = 0;
	for (uint j = 0; j < npv; j++)
	{
		for (uint i = 0; i < npu; i++)
		{
			cp[loc].act = 1;
			for (int dof = 0; dof < 3; dof++)
			{
				cp[loc].coor[dof] = pts_new[loc][dof];
			}
			cp[loc].w = pts_new[loc][3];
			for (int k = 0; k < 5; k++)
			{
				cp[loc].knotU[k] = kts_u1[i + k];
				cp[loc].knotV[k] = kts_v1[j + k];
			}
			for (int k = 0; k < 4; k++)
			{
				cp[loc].kitvU[k] = cp[loc].knotU[k + 1] - cp[loc].knotU[k];
				cp[loc].kitvV[k] = cp[loc].knotV[k + 1] - cp[loc].knotV[k];
			}
			cp[loc].pm[0] = cp[loc].knotU[2];
			cp[loc].pm[1] = cp[loc].knotV[2];
			cp[loc].index[0] = i;
			cp[loc].index[1] = j;
			loc++;
		}
	}
	loc = 0;
	for (uint j = 0; j < npv - 1; j++)
	{
		for (uint i = 0; i < npu - 1; i++)
		{
			tmesh[loc].act = 0;
			tmesh[loc].cnct[0] = npu * j + i;
			tmesh[loc].cnct[1] = npu * j + i + 1;
			tmesh[loc].cnct[2] = npu * (j + 1) + i + 1;
			tmesh[loc].cnct[3] = npu * (j + 1) + i;
			double len_u(cp[tmesh[loc].cnct[1]].pm[0] - cp[tmesh[loc].cnct[0]].pm[0]), len_v(cp[tmesh[loc].cnct[3]].pm[1] - cp[tmesh[loc].cnct[0]].pm[1]);
			tmesh[loc].act = 1;
			if (len_u != 0. && len_v != 0.)
			{
				//tmesh[loc].act=1;
				tmesh[loc].type = 0;
				tmesh[loc].IEN.resize(16);
				int tmp[16] = {(j - 1) * npu + i - 1, (j - 1) * npu + i, (j - 1) * npu + i + 1, (j - 1) * npu + i + 2, j * npu + i - 1, j * npu + i, j * npu + i + 1, j * npu + i + 2,
							   (j + 1) * npu + i - 1, (j + 1) * npu + i, (j + 1) * npu + i + 1, (j + 1) * npu + i + 2, (j + 2) * npu + i - 1, (j + 2) * npu + i, (j + 2) * npu + i + 1, (j + 2) * npu + i + 2};
				for (int k = 0; k < 16; k++)
				{
					tmesh[loc].IEN[k] = tmp[k];
				}
			}
			else if (len_u == 0. && len_v == 0.)
			{
				tmesh[loc].type = 3;
			}
			else
			{
				tmesh[loc].type = 2;
			}
			loc++;
		}
	}

	uanc.clear();
	vanc.clear();
	for (uint i = 0; i < npu; i++)
	{
		pair<int, double> utmp(i, kts_u1[i + 2]);
		uanc.push_back(utmp);
	}
	for (uint i = 0; i < npv; i++)
	{
		pair<int, double> vtmp(i, kts_v1[i + 2]);
		vanc.push_back(vtmp);
	}
}

void TruncatedTspline::InitialTopology_PlateHole()
{
	uint i, j;
	tmedge.clear();
	cornerEID.clear();
	for (i = 0; i < tmesh.size(); i++)
	{
		//tmesh[i].act=1;
		//tmesh[i].type=0;
		double lentmp[4];
		for (j = 0; j < 4; j++)
		{
			cp[tmesh[i].cnct[j]].face.push_back(i);
			Edge edtmp;
			edtmp.act = 1;
			edtmp.pt[0] = tmesh[i].cnct[j];
			edtmp.pt[1] = tmesh[i].cnct[(j + 1) % 4];
			if (cp[tmesh[i].cnct[j]].index[0] == cp[tmesh[i].cnct[(j + 1) % 4]].index[0])
			{
				edtmp.len = cp[tmesh[i].cnct[(j + 1) % 4]].pm[1] - cp[tmesh[i].cnct[j]].pm[1];
				if (edtmp.len < 0.)
					edtmp.len = -edtmp.len;
				lentmp[j] = edtmp.len;
			}
			else
			{
				edtmp.len = cp[tmesh[i].cnct[(j + 1) % 4]].pm[0] - cp[tmesh[i].cnct[j]].pm[0];
				if (edtmp.len < 0.)
					edtmp.len = -edtmp.len;
				lentmp[j] = edtmp.len;
			}
			vector<Edge>::iterator it = find(tmedge.begin(), tmedge.end(), edtmp);
			int edid(it - tmedge.begin());
			if (it == tmedge.end())
			{
				tmedge.push_back(edtmp);
			}
			tmesh[i].edge[j] = edid;
			tmedge[edid].face.push_back(i);
		}
		if ((lentmp[0] == 0. && lentmp[1] > 0.) || (lentmp[0] > 0. && lentmp[1] == 0.))
		{
			//tmesh[i].type=2;
		}
		else if (lentmp[0] == 0. && lentmp[1] == 0.)
		{
			//tmesh[i].type=3;
			vector<int>::iterator it = find(cornerEID.begin(), cornerEID.end(), i);
			if (it == cornerEID.end())
				cornerEID.push_back(i);
		}
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		cp[tmedge[i].pt[0]].edge.push_back(i);
		cp[tmedge[i].pt[1]].edge.push_back(i);
	}

	//FindLocalKnotVectors();
	//UpdateKnotVectors();
	//FindIENglb();
}

void TruncatedTspline::NeumannBC_PlateHole(vector<BezierElement> &bzmesh)
{
	double bc_l(-4.), bc_u(4.);
	for (uint i = 0; i < bzmesh.size(); i++)
	{
		int edid(-1);
		if (bzmesh[i].pts[1][0] == bc_l) //left boundary
		{
			edid = 0;
		}
		else if (bzmesh[i].pts[7][0] == bc_l)
		{
			edid = 1;
		}
		else if (bzmesh[i].pts[14][0] == bc_l)
		{
			edid = 2;
		}
		else if (bzmesh[i].pts[8][0] == bc_l)
		{
			edid = 3;
		}
		if (edid != -1)
		{
			vector<int> id;
			for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
			{
				if (cp[bzmesh[i].IEN[j]].coor[0] == bc_l)
				{
					id.push_back(j);
				}
			}
			bzmesh[i].neum_edge.push_back(edid);
			bzmesh[i].neum_ID.push_back(id);
		}
		edid = -1;
		if (bzmesh[i].pts[1][1] == bc_u) //upper boundary
		{
			edid = 0;
		}
		else if (bzmesh[i].pts[7][1] == bc_u)
		{
			edid = 1;
		}
		else if (bzmesh[i].pts[14][1] == bc_u)
		{
			edid = 2;
		}
		else if (bzmesh[i].pts[8][1] == bc_u)
		{
			edid = 3;
		}
		if (edid != -1)
		{
			vector<int> id;
			for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
			{
				if (cp[bzmesh[i].IEN[j]].coor[1] == bc_u)
				{
					id.push_back(j);
					//cout<<bzmesh[i].IEN[j]<<" ";
				}
			}
			//cout<<"\n";
			bzmesh[i].neum_edge.push_back(edid);
			bzmesh[i].neum_ID.push_back(id);
			//cout<<edid<<"\n";
			//cout<<tmesh[bzmesh[i].prt].cnct[edid]<<" "<<tmesh[bzmesh[i].prt].cnct[(edid+1)%4]<<"\n";
			//getchar();
		}
	}
}

void TruncatedTspline::run_PlateHole_XP(string fn)
{
	SetProblem_PlateHole_XP(fn);
	//InitialTopology();
	for (int i = 0; i < 1; i++)
		Refine_Global_XP();
	Refine_PlateHole_XP();
	//VisualizeTMesh("plate_hole_XP/CM_1");
	//VisualizeControlMesh("plate_hole_XP/CM_adapt_3");
	cout << "DOF: " << cp.size() << "\n";
	SetBezierMatIrrPatch();
	CollectActives();

	//VisualizeControlMesh("plate_hole_XP/CM_adapt_2");
	VisualizeSurface("plate_hole_XP/surf_3");
	cout << "surf done!\n";
	getchar();
}

void TruncatedTspline::Refine_Global_XP()
{
	vector<int> rfid, rftype;
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1)
		{
			rfid.push_back(i);
			if (tmesh[i].type == 0)
			{
				rftype.push_back(0);
			}
			else if (tmesh[i].type == 2)
			{
				rftype.push_back(3);
			}
			else if (tmesh[i].type == 4)
			{
				rftype.push_back(4);
			}
			//else
			//{
			//	cout<<"Not supported in global refinement!\n";
			//	getchar();
			//}
		}
	}
	Topo_Refine_Unstruct_glb(rfid, rftype);
	Geom_Refine_Unstruct_glb(rfid, rftype);
}

void TruncatedTspline::CreateMesh_XP_2(string fn, string fn1)
{
	//read seed mesh
	vector<array<double, 3>> coors;
	vector<array<int, 4>> cnct;
	int npt(0), nel(0);
	ifstream fin;
	fin.open(fn);
	if (fin.is_open())
	{
		string tmp;
		for (int i = 0; i < 4; i++)
			getline(fin, tmp);
		fin >> tmp >> npt >> tmp;
		coors.resize(npt);
		for (int i = 0; i < npt; i++)
		{
			fin >> coors[i][0] >> coors[i][1] >> coors[i][2];
		}
		getline(fin, tmp);
		fin >> tmp >> nel >> tmp;
		cnct.resize(nel);
		for (int i = 0; i < nel; i++)
		{
			fin >> tmp >> cnct[i][0] >> cnct[i][1] >> cnct[i][2] >> cnct[i][3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fn << "!\n";
	}

	//subdivide mesh once
	vector<array<int, 4>> cnct1(4 * cnct.size());
	for (int i = 0; i < nel; i++)
	{
		int pid[5];
		for (int j = 0; j < 4; j++)
		{
			array<double, 3> ptmp;
			for (int k = 0; k < 3; k++)
			{
				ptmp[k] = (coors[cnct[i][j]][k] + coors[cnct[i][(j + 1) % 4]][k]) / 2.;
			}
			vector<array<double, 3>>::iterator it = find(coors.begin(), coors.end(), ptmp);
			if (it == coors.end())
			{
				pid[j] = coors.size();
				coors.push_back(ptmp);
			}
			else
			{
				pid[j] = it - coors.begin();
			}
		}
		pid[4] = coors.size();
		array<double, 3> ptmp;
		for (int k = 0; k < 3; k++)
		{
			ptmp[k] = (coors[cnct[i][0]][k] + coors[cnct[i][1]][k] + coors[cnct[i][2]][k] + coors[cnct[i][3]][k]) / 4.;
		}
		coors.push_back(ptmp);
		int etmp[4][4] = {{cnct[i][0], pid[0], pid[4], pid[3]}, {pid[0], cnct[i][1], pid[1], pid[4]}, {pid[4], pid[1], cnct[i][2], pid[2]}, {pid[3], pid[4], pid[2], cnct[i][3]}};
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				cnct1[4 * i + j][k] = etmp[j][k];
			}
		}
	}

	//pillow layers

	//output mesh
	ofstream fout;
	fout.open(fn1);
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << coors.size() << " float\n";
		for (uint i = 0; i < coors.size(); i++)
		{
			//double tmp=coors[i][0]*(4.-coors[i][0])*coors[i][1]*(4.-coors[i][1])/8.;
			//fout<<coors[i][0]<<" "<<coors[i][1]<<" "<<tmp<<"\n";
			fout << coors[i][0] << " " << coors[i][1] << " " << coors[i][2] << "\n";
		}
		fout << "\nCELLS " << cnct1.size() << " " << 5 * cnct1.size() << '\n';
		for (uint i = 0; i < cnct1.size(); i++)
		{
			fout << "4 " << cnct1[i][0] << " " << cnct1[i][1] << " " << cnct1[i][2] << " " << cnct1[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << cnct1.size() << '\n';
		for (uint i = 0; i < cnct1.size(); i++)
		{
			fout << "9\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Cannot open " << fn1 << "!\n";
	}
}

void TruncatedTspline::CreateMesh_XP_2_Pillow(string fn, string fn1)
{
	//read seed mesh
	vector<array<double, 3>> coors;
	vector<array<int, 4>> cnct;
	int npt(0), nel(0);
	ifstream fin;
	fin.open(fn);
	if (fin.is_open())
	{
		string tmp;
		for (int i = 0; i < 4; i++)
			getline(fin, tmp);
		fin >> tmp >> npt >> tmp;
		coors.resize(npt);
		for (int i = 0; i < npt; i++)
		{
			fin >> coors[i][0] >> coors[i][1] >> coors[i][2];
		}
		getline(fin, tmp);
		fin >> tmp >> nel >> tmp;
		cnct.resize(nel);
		for (int i = 0; i < nel; i++)
		{
			fin >> tmp >> cnct[i][0] >> cnct[i][1] >> cnct[i][2] >> cnct[i][3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fn << "!\n";
	}

	//pillow layers
	int bot_id[5] = {0, 12, 1, 17, 2};
	int left_id[7] = {0, 15, 3, 23, 6, 33, 8};
	int right_id[5] = {2, 18, 5, 29, 11};
	int top_id[7] = {8, 32, 9, 26, 10, 30, 11};
	int nlayer(3);
	double len[3] = {0.6, 0.6, 0.2};
	for (int i = 1; i < nlayer; i++)
	{
		len[i] += len[i - 1];
	}
	int count(npt);
	int pnewb[4][5], pnewl[4][7], pnewr[4][5], pnewt[4][7], pbl[4][4], pbr[4][4], ptl[4][4], ptr[4][4];
	//bottom
	for (int j = 0; j < 5; j++)
	{
		pnewb[0][j] = bot_id[j];
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			array<double, 3> ptmp = {coors[bot_id[j]][0], coors[bot_id[j]][1] - len[i], coors[bot_id[j]][2]};
			pnewb[i + 1][j] = count;
			coors.push_back(ptmp);
			count++;
		}
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < 5 - 1; j++)
		{
			array<int, 4> etmp = {pnewb[i + 1][j], pnewb[i + 1][j + 1], pnewb[i][j + 1], pnewb[i][j]};
			cnct.push_back(etmp);
		}
	}
	//left
	for (int j = 0; j < 7; j++)
	{
		pnewl[0][j] = left_id[j];
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			array<double, 3> ptmp = {coors[left_id[j]][0] - len[i], coors[left_id[j]][1], coors[left_id[j]][2]};
			pnewl[i + 1][j] = count;
			coors.push_back(ptmp);
			count++;
		}
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < 7 - 1; j++)
		{
			array<int, 4> etmp = {pnewl[i + 1][j], pnewl[i][j], pnewl[i][j + 1], pnewl[i + 1][j + 1]};
			cnct.push_back(etmp);
		}
	}
	//right
	for (int j = 0; j < 5; j++)
	{
		pnewr[0][j] = right_id[j];
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			array<double, 3> ptmp = {coors[right_id[j]][0] + len[i], coors[right_id[j]][1], coors[right_id[j]][2]};
			pnewr[i + 1][j] = count;
			coors.push_back(ptmp);
			count++;
		}
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < 5 - 1; j++)
		{
			array<int, 4> etmp = {pnewr[i][j], pnewr[i + 1][j], pnewr[i + 1][j + 1], pnewr[i][j + 1]};
			cnct.push_back(etmp);
		}
	}
	//top
	for (int j = 0; j < 7; j++)
	{
		pnewt[0][j] = top_id[j];
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			array<double, 3> ptmp = {coors[top_id[j]][0], coors[top_id[j]][1] + len[i], coors[top_id[j]][2]};
			pnewt[i + 1][j] = count;
			coors.push_back(ptmp);
			count++;
		}
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < 7 - 1; j++)
		{
			array<int, 4> etmp = {pnewt[i][j], pnewt[i][j + 1], pnewt[i + 1][j + 1], pnewt[i + 1][j]};
			cnct.push_back(etmp);
		}
	}
	//bottom left
	for (int j = 0; j < nlayer + 1; j++)
	{
		pbl[0][j] = pnewl[j][0];
		pbl[j][0] = pnewb[j][0];
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < nlayer; j++)
		{
			array<double, 3> ptmp = {coors[left_id[0]][0] - len[j], coors[bot_id[0]][1] - len[i], 0.};
			pbl[i + 1][j + 1] = count;
			coors.push_back(ptmp);
			count++;
		}
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < nlayer; j++)
		{
			array<int, 4> etmp = {pbl[i + 1][j + 1], pbl[i + 1][j], pbl[i][j], pbl[i][j + 1]};
			cnct.push_back(etmp);
		}
	}
	//bottom right
	for (int j = 0; j < nlayer + 1; j++)
	{
		pbr[0][j] = pnewr[j][0];
		pbr[j][0] = pnewb[j][4];
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < nlayer; j++)
		{
			array<double, 3> ptmp = {coors[right_id[0]][0] + len[j], coors[bot_id[0]][1] - len[i], 0.};
			pbr[i + 1][j + 1] = count;
			coors.push_back(ptmp);
			count++;
		}
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < nlayer; j++)
		{
			array<int, 4> etmp = {pbr[i + 1][j], pbr[i + 1][j + 1], pbr[i][j + 1], pbr[i][j]};
			cnct.push_back(etmp);
		}
	}
	//top left
	for (int j = 0; j < nlayer + 1; j++)
	{
		ptl[0][j] = pnewl[j][6];
		ptl[j][0] = pnewt[j][0];
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < nlayer; j++)
		{
			array<double, 3> ptmp = {coors[left_id[0]][0] - len[j], coors[top_id[0]][1] + len[i], 0.};
			ptl[i + 1][j + 1] = count;
			coors.push_back(ptmp);
			count++;
		}
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < nlayer; j++)
		{
			array<int, 4> etmp = {ptl[i][j + 1], ptl[i][j], ptl[i + 1][j], ptl[i + 1][j + 1]};
			cnct.push_back(etmp);
		}
	}
	//top right
	for (int j = 0; j < nlayer + 1; j++)
	{
		ptr[0][j] = pnewr[j][4];
		ptr[j][0] = pnewt[j][6];
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < nlayer; j++)
		{
			array<double, 3> ptmp = {coors[right_id[0]][0] + len[j], coors[top_id[0]][1] + len[i], 0.};
			ptr[i + 1][j + 1] = count;
			coors.push_back(ptmp);
			count++;
		}
	}
	for (int i = 0; i < nlayer; i++)
	{
		for (int j = 0; j < nlayer; j++)
		{
			array<int, 4> etmp = {ptr[i][j], ptr[i][j + 1], ptr[i + 1][j + 1], ptr[i + 1][j]};
			cnct.push_back(etmp);
		}
	}
	//rescale points
	double max_tmp(4. + len[2] * 2);
	for (uint i = 0; i < coors.size(); i++)
	{
		coors[i][0] += len[2];
		coors[i][1] += len[2];
		coors[i][2] = coors[i][0] * (max_tmp - coors[i][0]) * coors[i][1] * (max_tmp - coors[i][1]) / (max_tmp * max_tmp);
	}

	//output mesh
	ofstream fout;
	fout.open(fn1);
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << coors.size() << " float\n";
		for (uint i = 0; i < coors.size(); i++)
		{
			//double tmp=coors[i][0]*(4.-coors[i][0])*coors[i][1]*(4.-coors[i][1])/8.;
			//fout<<coors[i][0]<<" "<<coors[i][1]<<" "<<tmp<<"\n";
			fout << coors[i][0] << " " << coors[i][1] << " " << coors[i][2] << "\n";
		}
		fout << "\nCELLS " << cnct.size() << " " << 5 * cnct.size() << '\n';
		for (uint i = 0; i < cnct.size(); i++)
		{
			fout << "4 " << cnct[i][0] << " " << cnct[i][1] << " " << cnct[i][2] << " " << cnct[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << cnct.size() << '\n';
		for (uint i = 0; i < cnct.size(); i++)
		{
			fout << "9\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Cannot open " << fn1 << "!\n";
	}
}

void TruncatedTspline::SetBezier3TranMat(int N, vector<vector<double>> &bmat)
{
	int nb = 2 * N + 8;
	bmat.resize(nb, vector<double>(16, 0.));
	double a(4. / 9.), bN(4. / (9. * double(N))), cN(1. / (9. * double(N))), b(2. / 9.), c(1. / 9.), d(1. / 18.), e(1. / 36);

	//1 extraordinary node coef
	bmat[0][0] = a;
	for (int i = 0; i < N; i++)
	{
		bmat[2 * i + 1][0] = bN;
		bmat[2 * i + 2][0] = cN;
	}
	//3 regular corner coefs
	int qv[3] = {4, 16, 13};
	int pv[3][9] = {{6, 5, 2 * N + 3, 2 * N + 4, 2 * N + 5, 7, 8, 1, 4}, {5, 4, 2 * N + 7, 2 * N + 6, 2 * N + 2, 2 * N + 3, 2 * N + 4, 6, 1}, {4, 1, 2, 3, 2 * N + 8, 2 * N + 7, 2 * N + 6, 5, 6}};
	if (N == 3)
		pv[0][6] = 2;
	double dv[9] = {a, c, e, c, e, c, e, c, e};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 9; j++)
		{
			bmat[pv[i][j] - 1][qv[i] - 1] = dv[j];
		}
	}
	//4 face coefs
	int qf[4] = {6, 7, 10, 11};
	int pf[4][4] = {{1, 6, 4, 5}, {6, 5, 1, 4}, {4, 1, 5, 6}, {5, 4, 6, 1}};
	double df[4] = {a, b, b, c};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			bmat[pf[i][j] - 1][qf[i] - 1] = df[j];
		}
	}
	//8 edge coefs
	int qe[8] = {2, 3, 5, 9, 8, 12, 14, 15};
	int pe[8][6] = {{1, 6, 4, 8, 5, 7}, {6, 1, 5, 7, 4, 8}, {1, 4, 2, 6, 3, 5}, {4, 1, 3, 5, 2, 6}, {6, 5, 1, 2 * N + 4, 4, 2 * N + 3}, {5, 6, 4, 2 * N + 3, 1, 2 * N + 4}, {4, 5, 1, 2 * N + 7, 6, 2 * N + 6}, {5, 4, 6, 2 * N + 6, 1, 2 * N + 7}};
	if (N == 3)
	{
		pe[0][3] = 2;
		pe[1][5] = 2;
	}
	double de[6] = {a, b, c, c, d, d};
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			bmat[pe[i][j] - 1][qe[i] - 1] = de[j];
		}
	}

	//update first point as Catmull-Clark limit point
	//SingularPatchEval cceigen;
	//cceigen.ReadEigenStruct();
	//vector<double> cc_coef;
	//cceigen.ObtainIVcoef(N,cc_coef);
	//for(uint i=0; i<cc_coef.size(); i++)
	//{
	//	bmat[i][0]=cc_coef[i];
	//}

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

void TruncatedTspline::SetBezier4TranMat(int N, vector<vector<double>> &bmat)
{
	int p(3);
	vector<vector<double>> b3mat;
	SetBezier3TranMat(N, b3mat);
	//degree elevation
	vector<vector<double>> demat(25, vector<double>(16, 0.));
	int loc1(0);
	for (int j = 0; j < 5; j++)
	{
		for (int i = 0; i < 5; i++)
		{
			double a(double(i) / double(p + 1)), b(double(j) / double(p + 1));
			double coef[4] = {(1. - a) * (1. - b), a * (1. - b), (1. - a) * b, a * b};
			int loc0[4] = {4 * j + i, 4 * j + i - 1, 4 * (j - 1) + i, 4 * (j - 1) + i - 1};
			if (i == 0)
			{
				loc0[1] = -1;
				loc0[3] = -1;
			}
			if (j == 0)
			{
				loc0[2] = -1;
				loc0[3] = -1;
			}
			if (i == 4)
			{
				loc0[0] = -1;
				loc0[2] = -1;
			}
			if (j == 4)
			{
				loc0[0] = -1;
				loc0[1] = -1;
			}
			for (int k = 0; k < 4; k++)
			{
				if (loc0[k] != -1)
				{
					demat[loc1][loc0[k]] = coef[k];
				}
			}
			loc1++;
		}
	}
	bmat.resize(2 * N + 8, vector<double>(25, 0.));
	for (uint i = 0; i < bmat.size(); i++)
	{
		for (uint j = 0; j < bmat[i].size(); j++)
		{
			for (int k = 0; k < 16; k++)
			{
				bmat[i][j] += b3mat[i][k] * demat[j][k];
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

void TruncatedTspline::SetBezier4TranMatOP(int N, vector<vector<double>> &bmat)
{
	bmat.clear();
	ifstream fin;
	fin.open("../src/bmat1/bmat_" + to_string(N) + ".txt");
	if (fin.is_open())
	{
		int dim1, dim2;
		fin >> dim1 >> dim2;
		bmat.resize(dim1, vector<double>(dim2));
		for (int i = 0; i < dim1; i++)
		{
			for (int j = 0; j < dim2; j++)
			{
				fin >> bmat[i][j];
			}
		}
	}
	else
	{
		cerr << "Cannot open bmat file!\n";
	}
}

void TruncatedTspline::SurfacePointCal(double u, double v, const vector<int> &IEN, const vector<vector<double>> &bmat, array<double, 3> &pcoor, array<double, 3> &norm)
{
	int nv((IEN.size() - 8) / 2);
	BezierElement be;
	double Nt0[16], dNdt0[16][2];
	be.Basis(u, v, Nt0, dNdt0);
	vector<double> Nt(IEN.size(), 0.);
	pcoor[0] = 0.;
	pcoor[1] = 0.;
	pcoor[2] = 0.;
	double nmtmp[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
	for (uint i = 0; i < Nt.size(); i++)
	{
		double dNdt[2] = {0., 0.};
		for (int j = 0; j < 16; j++)
		{
			Nt[i] += bmat[i][j] * Nt0[j];
			dNdt[0] += bmat[i][j] * dNdt0[j][0];
			dNdt[1] += bmat[i][j] * dNdt0[j][1];
		}
		pcoor[0] += cp[IEN[i]].coor[0] * Nt[i];
		pcoor[1] += cp[IEN[i]].coor[1] * Nt[i];
		pcoor[2] += cp[IEN[i]].coor[2] * Nt[i];
		nmtmp[0][0] += cp[IEN[i]].coor[0] * dNdt[0];
		nmtmp[0][1] += cp[IEN[i]].coor[1] * dNdt[0];
		nmtmp[0][2] += cp[IEN[i]].coor[2] * dNdt[0];
		nmtmp[1][0] += cp[IEN[i]].coor[0] * dNdt[1];
		nmtmp[1][1] += cp[IEN[i]].coor[1] * dNdt[1];
		nmtmp[1][2] += cp[IEN[i]].coor[2] * dNdt[1];
	}
	norm[0] = nmtmp[0][1] * nmtmp[1][2] - nmtmp[0][2] * nmtmp[1][1];
	norm[1] = nmtmp[0][2] * nmtmp[1][0] - nmtmp[0][0] * nmtmp[1][2];
	norm[2] = nmtmp[0][0] * nmtmp[1][1] - nmtmp[0][1] * nmtmp[1][0];
	double len = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
	norm[0] /= len;
	norm[1] /= len;
	norm[2] /= len;
}

void TruncatedTspline::SurfacePointCal1(double u, double v, const vector<int> &IEN, const vector<vector<double>> &bmat, array<double, 3> &pcoor, array<double, 3> &norm)
{
	int nv((IEN.size() - 8) / 2);
	BezierElement be;
	double Nt0[25], dNdt0[25][2];
	be.Basis4(u, v, Nt0, dNdt0);
	vector<double> Nt(IEN.size(), 0.);
	pcoor[0] = 0.;
	pcoor[1] = 0.;
	pcoor[2] = 0.;
	double nmtmp[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
	for (uint i = 0; i < Nt.size(); i++)
	{
		double dNdt[2] = {0., 0.};
		for (int j = 0; j < 25; j++)
		{
			Nt[i] += bmat[i][j] * Nt0[j];
			dNdt[0] += bmat[i][j] * dNdt0[j][0];
			dNdt[1] += bmat[i][j] * dNdt0[j][1];
		}
		pcoor[0] += cp[IEN[i]].coor[0] * Nt[i];
		pcoor[1] += cp[IEN[i]].coor[1] * Nt[i];
		pcoor[2] += cp[IEN[i]].coor[2] * Nt[i];
		nmtmp[0][0] += cp[IEN[i]].coor[0] * dNdt[0];
		nmtmp[0][1] += cp[IEN[i]].coor[1] * dNdt[0];
		nmtmp[0][2] += cp[IEN[i]].coor[2] * dNdt[0];
		nmtmp[1][0] += cp[IEN[i]].coor[0] * dNdt[1];
		nmtmp[1][1] += cp[IEN[i]].coor[1] * dNdt[1];
		nmtmp[1][2] += cp[IEN[i]].coor[2] * dNdt[1];
	}
	norm[0] = nmtmp[0][1] * nmtmp[1][2] - nmtmp[0][2] * nmtmp[1][1];
	norm[1] = nmtmp[0][2] * nmtmp[1][0] - nmtmp[0][0] * nmtmp[1][2];
	norm[2] = nmtmp[0][0] * nmtmp[1][1] - nmtmp[0][1] * nmtmp[1][0];
	double len = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
	norm[0] /= len;
	norm[1] /= len;
	norm[2] /= len;
}

void TruncatedTspline::VisualizeSurface(string fn)
{
	vector<array<double, 3>> spt;
	vector<array<double, 3>> sval;
	vector<array<int, 4>> sele;
	vector<array<double, 3>> lpt; //visulize parameter lines
	vector<array<int, 2>> led;	  //line connectivity
	int ns(11), ecount(0), loc0, loc1, loc2;
	//vector<double> su(ns),sv(ns);
	//for(int i=0; i<ns; i++)
	//{
	//	su[i]=i*1./(ns-1);
	//	sv[i]=i*1./(ns-1);
	//}
	double errmax(0.);

	for (uint e = 0; e < tmesh.size(); e++)
	{
		if (tmesh[e].act == 1 && (tmesh[e].type == 0 || tmesh[e].type == 1 || tmesh[e].type == 4))
		{
			//int loc(0);
			//vector<vector<double>> bmat;
			//SetBezier3TranMat((tmesh[e].IEN.size()-8)/2,bmat);
			//SetBezier4TranMat((tmesh[e].IEN.size()-8)/2,bmat);
			//SetBezier4TranMatOP((tmesh[e].IEN.size()-8)/2,bmat);
			//SolveOptimizeBezierMat((tmesh[e].IEN.size()-8)/2,bmat);
			//for(int i=0; i<25; i++)
			//{
			//	array<double,3> bept={0.,0.,0.};
			//	for(uint j=0; j<bmat.size(); j++)
			//	{
			//		bept[0]+=bmat[j][i]*cp[tmesh[e].IEN[j]].coor[0];
			//		bept[1]+=bmat[j][i]*cp[tmesh[e].IEN[j]].coor[1];
			//		bept[2]+=bmat[j][i]*cp[tmesh[e].IEN[j]].coor[2];
			//	}
			//	spt.push_back(bept);
			//}

			vector<double> su(ns), sv(ns);
			for (int i = 0; i < ns; i++)
			{
				su[i] = i * tmedge[tmesh[e].edge[0]].len / (ns - 1);
				sv[i] = i * tmedge[tmesh[e].edge[3]].len / (ns - 1);
			}

			for (int a = 0; a < ns; a++)
			{
				for (int b = 0; b < ns; b++)
				{
					array<double, 3> pt;
					array<double, 3> nm;
					//SurfacePointCal1(su[b],sv[a],tmesh[e].IEN,bmat,pt,nm);
					//SurfacePointMap(e,su[b],sv[a],pt,nm);
					GeomMap_DPatch(e, su[b], sv[a], pt);
					spt.push_back(pt);
					sval.push_back(nm);
					//double errtmp = fabs(1. - sqrt(pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2]));
					double errtmp = fabs(100. - sqrt(pt[0] * pt[0] + pt[1] * pt[1]));
					if (errtmp > errmax)
						errmax = errtmp;
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
					el[0] = ecount * ns * ns + a * ns + b;
					el[1] = ecount * ns * ns + a * ns + b + 1;
					el[2] = ecount * ns * ns + (a + 1) * ns + b + 1;
					el[3] = ecount * ns * ns + (a + 1) * ns + b;
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
	}

	cout << "err max: " << errmax << "\n";

	string fname = fn + "_geom.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i < spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
		}
		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i]<<"\n";
		//}
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

	string fname1(fn + "_geom-lines.vtk");
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

void TruncatedTspline::SetProblem1(string fn)
{
	//read quad vtk
	string fname(fn), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++)
			getline(fin, stmp); //skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
			//cp[i].coor[2]=cp[i].coor[0]*(4.-cp[i].coor[0])*cp[i].coor[1]*(4.-cp[i].coor[1])/8.;
			cp[i].act = 1;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			tmesh[i].act = 1;
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	int type4[8] = {2, 7, 9, 12, 16, 10, 15, 21};
	int iens[8][18] = {{4, 19, 20, 13, 16, 14, 24, 21, 28, 25, 31, 0, 15, 3, 23, 12, 1, 17}, {4, 25, 31, 19, 20, 13, 16, 14, 24, 21, 28, 2, 17, 1, 12, 18, 5, 29}, {4, 13, 16, 14, 24, 21, 28, 25, 31, 19, 20, 6, 22, 7, 27, 23, 3, 15}, {4, 14, 24, 21, 28, 25, 31, 19, 20, 13, 16, 9, 26, 10, 30, 27, 7, 22}, {4, 21, 28, 25, 31, 19, 20, 13, 16, 14, 24, 11, 29, 5, 18, 30, 10, 26}, {7, 27, 28, 21, 24, 22, 34, 3, 23, 6, 33, 14, 4, 25, -1, -1, -1, -1}, {7, 22, 34, 27, 28, 21, 24, 10, 25, 4, 14, 26, 9, 32, -1, -1, -1, -1}, {7, 21, 24, 22, 34, 27, 28, 8, 32, 9, 26, 33, 6, 23, -1, -1, -1, -1}};
	for (int i = 0; i < 8; i++)
	{
		tmesh[type4[i]].type = 4;
		tmesh[type4[i]].IEN.clear();
		for (int j = 0; j < 18; j++)
		{
			if (iens[i][j] != -1)
				tmesh[type4[i]].IEN.push_back(iens[i][j]);
		}
	}
}

void TruncatedTspline::GlueIrregularPatches4OP(vector<vector<int>> &pats)
{
	vector<int> found;
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].type == 2)
		{
			vector<int>::iterator it = find(found.begin(), found.end(), i);
			if (it == found.end())
			{
				found.push_back(i);
				vector<int> XPid, add0, add1;
				XPid.push_back(i);
				add0.push_back(i);
				while (1)
				{
					for (uint j0 = 0; j0 < add0.size(); j0++)
					{
						for (uint j = 0; j < cp[add0[j0]].face.size(); j++)
						{
							int fid0(cp[XPid.back()].face[j]);
							for (int k = 0; k < 4; k++)
							{
								if (tmedge[tmesh[fid0].edge[k]].act == 1)
								{
									if (tmedge[tmesh[fid0].edge[k]].face.size() == 1)
										break;
									int fid1(tmedge[tmesh[fid0].edge[k]].face[0]);
									if (fid1 == fid0 && tmedge[tmesh[fid0].edge[k]].face.size() == 2)
									{
										fid1 = tmedge[tmesh[fid0].edge[k]].face[1];
									}
									for (int k1 = 0; k1 < 4; k1++)
									{
										if (cp[tmesh[fid1].cnct[k1]].type == 4)
										{
											vector<int>::iterator it1 = find(add1.begin(), add1.end(), tmesh[fid1].cnct[k1]);
											vector<int>::iterator it2 = find(XPid.begin(), XPid.end(), tmesh[fid1].cnct[k1]);
											if (it1 == add1.end() && it2 == XPid.end())
											{
												add1.push_back(tmesh[fid1].cnct[k1]);
												XPid.push_back(tmesh[fid1].cnct[k1]);
												found.push_back(tmesh[fid1].cnct[k1]);
											}
										}
									}
								}
								else
								{
									cerr << "T-junctions around XP are not supported yet!\n";
									getchar();
								}
							}
						}
					}
					if (add1.size() == 0)
						break;
					else
					{
						add0.clear();
						add0 = add1;
						add1.clear();
					}
				}
				vector<int> pat_tmp;
				for (uint j = 0; j < XPid.size(); j++)
				{
					for (uint k = 0; k < cp[XPid[j]].face.size(); k++)
					{
						pat_tmp.push_back(cp[XPid[j]].face[k]);
					}
				}
				pats.push_back(pat_tmp);
			}
		}
	}
}

void TruncatedTspline::CapIndex_Loc2Glb(int N, vector<vector<int>> &ICN)
{
	ICN.resize(N, vector<int>(25));
	for (int i = 0; i < N - 1; i++)
	{
		ICN[i][0] = 0;
		int loc(1);
		for (int j = 1; j < 25; j++)
		{
			int a(j % 5), b(j / 5);
			if (a == 0 && b != 0)
				ICN[i][j] = 20 * (i + 1) + b;
			else
			{
				ICN[i][j] = 20 * i + loc;
				loc++;
			}
		}
	}
	ICN[N - 1][0] = 0;
	int loc(1);
	for (int j = 1; j < 25; j++)
	{
		int a(j % 5), b(j / 5);
		if (a == 0 && b != 0)
			ICN[N - 1][j] = b;
		else
		{
			ICN[N - 1][j] = 20 * (N - 1) + loc;
			loc++;
		}
	}
}

void TruncatedTspline::BuildGeomMat(int N, const vector<int> &bloc, const vector<vector<double>> &bmat, const vector<vector<int>> &ICN, vector<vector<double>> &gmat, vector<double> &gvec)
{
	double PI(3.1415926535);
	double lamda(-2. * cos(2. * PI / N));
	int fix[14] = {3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 16, 17, 21, 22};
	gmat.resize(17 * N + 2, vector<double>(20 * N + 1, 0.));
	gvec.resize(17 * N + 2, 0.);
	int rs(0);
	//gmat[rs][ICN[N-1][1]]=1.; gmat[rs][ICN[0][0]]=-2.-lamda; gmat[rs][ICN[0][1]]=lamda; gmat[rs][ICN[0][5]]=1.;
	//rs++;
	gmat[rs][ICN[0][0]] = -double(N);
	for (int i = 0; i < N; i++)
	{
		gmat[rs][ICN[i][1]] = 1.;
	}
	rs++;
	//gmat[rs][ICN[0][0]]=-double(N);
	//for(int i=0; i<N; i++)
	//{
	//	gmat[rs][ICN[i][6]]=1.;
	//}
	//rs++;
	gmat[rs][0] = 1.;
	gvec[rs] = bmat[bloc[0]][0];
	rs++;
	for (int i = 0; i < N; i++) //loop each face
	{
		for (int j = 0; j < 14; j++) // 14 fixed coefficients
		{
			gmat[rs][ICN[i][fix[j]]] = 1.;
			gvec[rs] = bmat[bloc[i]][fix[j]];
			rs++;
		}
		//6 boundary equations
		int i0((i - 1 + N) % N);
		gmat[rs][ICN[i][0]] = 1.;
		gmat[rs][ICN[i][1]] = -4.;
		gmat[rs][ICN[i][2]] = 6.;
		gmat[rs][ICN[i][3]] = -4.;
		gmat[rs][ICN[i][4]] = 1.;
		rs++;
		//gmat[rs][ICN[i0][1]]=1.; gmat[rs][ICN[i][0]]=-2.-lamda; gmat[rs][ICN[i][1]]=lamda; gmat[rs][ICN[i][5]]=1.;
		//rs++;
		gmat[rs][ICN[i0][6]] = 4.;
		gmat[rs][ICN[i][1]] = -8. - 2. * lamda;
		gmat[rs][ICN[i][0]] = lamda / 2.;
		gmat[rs][ICN[i][3]] = 2. * lamda;
		gmat[rs][ICN[i][4]] = -lamda / 2.;
		gmat[rs][ICN[i][6]] = 4.;
		rs++;
		gmat[rs][ICN[i0][11]] = 4.;
		gmat[rs][ICN[i][2]] = -8.;
		gmat[rs][ICN[i][4]] = 2. * lamda / 3.;
		gmat[rs][ICN[i][3]] = -2. * lamda / 3.;
		gmat[rs][ICN[i][7]] = 4.;
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

void TruncatedTspline::BuildFairMat(int N, const vector<int> &bloc, const vector<vector<double>> &bmat, const vector<vector<int>> &ICN, vector<vector<double>> &fmat, vector<double> &fvec)
{
	fmat.resize(40 * N, vector<double>(20 * N + 1, 0.));
	fvec.resize(40 * N, 0.);
	int count(0);
	for (int i = 0; i < N; i++) //loop each face
	{
		for (int k = 0; k < 5; k++)
		{
			for (int j = 0; j < 4; j++)
			{
				int loc(5 * k + j), loc1(5 * k + j + 1);
				fmat[count][ICN[i][loc]] = 1.;
				fmat[count][ICN[i][loc1]] = -1.;
				fvec[count] = bmat[bloc[i]][loc] - bmat[bloc[i]][loc1];
				count++;
			}
		}
		for (int k = 0; k < 4; k++)
		{
			for (int j = 0; j < 5; j++)
			{
				int loc(5 * k + j), loc1(5 * (k + 1) + j);
				fmat[count][ICN[i][loc]] = 1.;
				fmat[count][ICN[i][loc1]] = -1.;
				fvec[count] = bmat[bloc[i]][loc] - bmat[bloc[i]][loc1];
				count++;
			}
		}
	}
}

void TruncatedTspline::Convert2MatlabData(const vector<vector<double>> &mat, vector<double> &row_id, vector<double> &col_id, vector<double> &coef)
{
	for (uint j = 0; j < mat[0].size(); j++)
	{
		for (uint i = 0; i < mat.size(); i++)
		{
			if (mat[i][j] != 0.)
			{
				row_id.push_back(double(i + 1));
				col_id.push_back(double(j + 1));
				coef.push_back(mat[i][j]);
			}
		}
	}
}

void TruncatedTspline::PreConditionCoef(int N, const vector<int> &bloc, const vector<vector<double>> &bmat, const vector<vector<int>> &ICN, vector<double> &coef)
{
	coef.resize(20 * N + 1, 0.);
	for (uint i = 0; i < ICN.size(); i++)
	{
		for (uint j = 0; j < ICN[i].size(); j++)
		{
			coef[ICN[i][j]] = bmat[bloc[i]][j];
		}
	}
}

void TruncatedTspline::SolveOptimizeBezierMat(int N, vector<vector<double>> &bmatop)
{
	vector<vector<int>> ICN;
	vector<vector<double>> bmat;
	CapIndex_Loc2Glb(N, ICN);
	SetBezier4TranMat(N, bmat);

	bmatop.resize(bmat.size(), vector<double>(bmat[0].size()));
	//ofstream fout;
	//fout.open("bmat1/bmat_5_0.txt");
	for (uint i = 0; i < bmat.size(); i++)
	{
		for (uint j = 0; j < bmat[i].size(); j++)
		{
			bmatop[i][j] = bmat[i][j];
			//fout<<bmat[i][j]<<" ";
		}
		//fout<<"\n";
	}
	//fout.close();
	//getchar();
	vector<vector<int>> G2L(2 * N + 1, vector<int>(N));
	vector<vector<int>> L2G(N, vector<int>(2 * N + 1));
	for (int i = 0; i < N; i++)
	{
		L2G[i][0] = 0;
		G2L[0][i] = 0;
		for (int j = 1; j < 2 * N + 1; j++)
		{
			L2G[i][j] = (-2 * i + j - 1 + 2 * N) % (2 * N) + 1;
			G2L[L2G[i][j]][i] = j;
		}
	}

	for (int bfid = 0; bfid < 3; bfid++) //optimize for one-ring basis functions
	{
		cout << "Basis function: " << bfid << "\n";
		vector<vector<double>> gmat, fmat;
		vector<double> gvec, fvec;
		BuildGeomMat(N, G2L[bfid], bmat, ICN, gmat, gvec);
		BuildFairMat(N, G2L[bfid], bmat, ICN, fmat, fvec);

		vector<double> gid1, gid2, gc, fid1, fid2, fc, precf;
		Convert2MatlabData(gmat, gid1, gid2, gc);
		Convert2MatlabData(fmat, fid1, fid2, fc);
		PreConditionCoef(N, G2L[bfid], bmat, ICN, precf);

		/*MatlabSolver solver;
		cout<<"solving\n";
		vector<double> sol(20*N+1);
		solver.Initilize2(gid1,gid2,gc,gvec,fid1,fid2,fc,fvec,precf);
		solver.Solve3(sol.data());*/

		//for(int i=0; i<25; i++)
		//{
		//	if(bfid==0)
		//	{
		//		bmatop[0][i]=sol[ICN[0][i]];
		//	}
		//	else if(bfid==1)
		//	{
		//		for(int j=0; j<N; j++)
		//		{
		//			int rowid=2*j+1;
		//			bmatop[rowid][i]=sol[ICN[j][i]];
		//			if(bmatop[rowid][i]<1.e-10) bmatop[rowid][i]=0.;
		//		}
		//	}
		//	else if(bfid==2)
		//	{
		//		for(int j=0; j<N; j++)
		//		{
		//			int rowid=2*j+2;
		//			bmatop[rowid][i]=sol[ICN[j][i]];
		//			if(bmatop[rowid][i]<1.e-10) bmatop[rowid][i]=0.;
		//		}
		//	}
		//}
	}

	for (uint i = 0; i < bmatop[0].size(); i++)
	{
		double sum(0.);
		for (uint j = 0; j < bmatop.size(); j++)
		{
			sum += bmatop[j][i];
		}
		for (uint j = 0; j < bmatop.size(); j++)
		{
			bmatop[j][i] /= sum;
		}
	}

	ofstream fout;
	fout.open("bmat1/bmat_" + to_string(N) + ".txt");
	if (fout.is_open())
	{
		fout << bmatop.size() << " " << bmatop[0].size() << "\n";
		for (uint i = 0; i < bmatop.size(); i++)
		{
			for (uint j = 0; j < bmatop[i].size(); j++)
			{
				fout << bmatop[i][j] << " ";
			}
			fout << "\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Cannot open bmat file!\n";
	}
}

void TruncatedTspline::SetProblem_surf(string fn)
{
	//read quad vtk
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++)
			getline(fin, stmp); //skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
			cp[i].act = 1;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	int nx(11), ny(11);
	double kntu[11] = {0., 0., 1., 2., 3., 4., 5., 6., 7., 8., 8.};
	double kntv[11] = {0., 0., 1., 2., 3., 4., 5., 6., 7., 8., 8.};
	int loc(0);
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			cp[loc].index[0] = i;
			cp[loc].index[1] = j;
			cp[loc].pm[0] = kntu[i];
			cp[loc].pm[1] = kntv[j];
			loc++;
		}
	}
	uanc.clear();
	vanc.clear();
	for (int i = 0; i < nx; i++)
	{
		pair<int, double> utmp(i, kntu[i]);
		uanc.push_back(utmp);
	}
	for (int i = 0; i < ny; i++)
	{
		pair<int, double> vtmp(i, kntv[i]);
		vanc.push_back(vtmp);
	}
}

void TruncatedTspline::SetLshapeProblem(string fn)
{
	//read quad vtk
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++)
			getline(fin, stmp); //skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
			cp[i].act = 1;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	double ku[21] = {0., 0., 0., 0., 1., 2., 3., 4., 5., 6., 6., 6., 7., 8., 9., 10., 11., 12., 12., 12., 12};
	double kv[13] = {0., 0., 0., 0., 1., 2., 3., 4., 5., 6., 6., 6., 6.};
	int nu(17), nv(9);

	int loc(0);
	for (int j = 0; j < nv; j++)
	{
		for (int i = 0; i < nu; i++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[loc].knotU[k] = ku[i + k];
				cp[loc].knotV[k] = kv[j + k];
			}
			for (int k = 0; k < 4; k++)
			{
				cp[loc].kitvU[k] = cp[loc].knotU[k + 1] - cp[loc].knotU[k];
				cp[loc].kitvV[k] = cp[loc].knotV[k + 1] - cp[loc].knotV[k];
			}
			cp[loc].pm[0] = cp[loc].knotU[2];
			cp[loc].pm[1] = cp[loc].knotV[2];
			loc++;
		}
	}
	loc = 0;
	for (int j = 0; j < nv - 1; j++)
	{
		for (int i = 0; i < nu - 1; i++)
		{
			if (i != 0 && i != nu - 2 && j != 0 && j != nv - 2)
			{
				tmesh[loc].act = 1;
				tmesh[loc].IEN.resize(16);
				int loc1(0);
				for (int j1 = j - 1; j1 < j + 3; j1++)
				{
					for (int i1 = i - 1; i1 < i + 3; i1++)
					{
						tmesh[loc].IEN[loc1] = j1 * nu + i1;
						loc1++;
					}
				}
			}
			loc++;
		}
	}

	loc = 0;
	for (int j = 0; j < nv; j++)
	{
		for (int i = 0; i < nu; i++)
		{
			cp[loc].index[0] = i;
			cp[loc].index[1] = j;
			loc++;
		}
	}
	uanc.clear();
	vanc.clear();
	for (int i = 0; i < nu; i++)
	{
		pair<int, double> utmp(i, ku[i + 2]);
		uanc.push_back(utmp);
	}
	for (int i = 0; i < nv; i++)
	{
		pair<int, double> vtmp(i, kv[i + 2]);
		vanc.push_back(vtmp);
	}
}

void TruncatedTspline::SetLshapeProblemLocalFeature(string fn)
{
	//read quad vtk
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++)
			getline(fin, stmp); //skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
			cp[i].act = 1;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	//add local features
	vector<Vertex> ptmp(8);
	vector<int> ptmp_id(8);
	for (int i = 0; i < 8; i++)
		ptmp[i].act = 1;
	int pid_old[8][2] = {{8, 7}, {8, 9}, {25, 24}, {25, 26}, {144, 143}, {144, 145}, {127, 126}, {127, 128}};
	for (int i = 0; i < 8; i++)
	{
		ptmp[i].act = 1;
		for (int dof = 0; dof < 3; dof++)
		{
			ptmp[i].coor[dof] = (2. * cp[pid_old[i][0]].coor[dof] + cp[pid_old[i][1]].coor[dof]) / 3.;
		}
		ptmp_id[i] = cp.size();
		cp.push_back(ptmp[i]);
	}
	int pid_old1[2][3] = {{24, 25, 26}, {126, 127, 128}};
	int id_updt[2] = {25, 127};
	for (int i = 0; i < 2; i++)
	{
		for (int dof = 0; dof < 3; dof++)
		{
			cp[id_updt[i]].coor[dof] = (cp[pid_old1[i][0]].coor[dof] + 4. * cp[pid_old1[i][1]].coor[dof] + cp[pid_old1[i][2]].coor[dof]) / 6.;
		}
	}

	//define anchor and index
	double ku[21] = {0., 0., 0., 0., 1., 2., 3., 4., 5., 6., 6., 6., 7., 8., 9., 10., 11., 12., 12., 12., 12};
	double kv[13] = {0., 0., 0., 0., 1., 2., 3., 4., 5., 6., 6., 6., 6.};
	int nu(17), nv(9), loc(0);
	double anc_u[17] = {0., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 14.};
	double anc_v[9] = {0., 0., 1., 2., 3., 4., 5., 6., 6.};
	int index_u[17] = {0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14, 15, 16, 17, 18};
	int index_v[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
	double anc_u1[8] = {7., 7., 7., 7., 7., 7., 7., 7.};
	double anc_v1[8] = {0., 0., 0., 0., 6., 6., 6., 6.};
	int index_u1[8] = {8, 10, 8, 10, 8, 10, 8, 10};
	int index_v1[8] = {0, 0, 1, 1, 8, 8, 7, 7};
	for (int j = 0; j < nv; j++)
	{
		for (int i = 0; i < nu; i++)
		{
			cp[loc].index[0] = i;
			cp[loc].index[1] = j;
			cp[loc].pm[0] = anc_u[i];
			cp[loc].pm[1] = anc_v[j];
			loc++;
		}
	}
	for (int i = 0; i < 8; i++)
	{
		cp[ptmp_id[i]].index[0] = index_u1[i];
		cp[ptmp_id[i]].index[1] = index_v1[i];
		cp[ptmp_id[i]].pm[0] = anc_u1[i];
		cp[ptmp_id[i]].pm[1] = anc_v1[i];
	}

	loc = 0;
	for (int j = 0; j < nv; j++)
	{
		for (int i = 0; i < nu; i++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[loc].knotU[k] = ku[i + k];
				cp[loc].knotV[k] = kv[j + k];
			}
			for (int k = 0; k < 4; k++)
			{
				cp[loc].kitvU[k] = cp[loc].knotU[k + 1] - cp[loc].knotU[k];
				cp[loc].kitvV[k] = cp[loc].knotV[k + 1] - cp[loc].knotV[k];
			}
			cp[loc].pm[0] = cp[loc].knotU[2];
			cp[loc].pm[1] = cp[loc].knotV[2];
			loc++;
		}
	}
	loc = 0;
	for (int j = 0; j < nv - 1; j++)
	{
		for (int i = 0; i < nu - 1; i++)
		{
			if (i != 0 && i != nu - 2 && j != 0 && j != nv - 2)
			{
				tmesh[loc].act = 1;
				tmesh[loc].IEN.resize(16);
				int loc1(0);
				for (int j1 = j - 1; j1 < j + 3; j1++)
				{
					for (int i1 = i - 1; i1 < i + 3; i1++)
					{
						tmesh[loc].IEN[loc1] = j1 * nu + i1;
						loc1++;
					}
				}
			}
			loc++;
		}
	}

	loc = 0;
	for (int j = 0; j < nv; j++)
	{
		for (int i = 0; i < nu; i++)
		{
			cp[loc].index[0] = i;
			cp[loc].index[1] = j;
			loc++;
		}
	}
	uanc.clear();
	vanc.clear();
	for (int i = 0; i < nu; i++)
	{
		pair<int, double> utmp(i, ku[i + 2]);
		uanc.push_back(utmp);
	}
	for (int i = 0; i < nv; i++)
	{
		pair<int, double> vtmp(i, kv[i + 2]);
		vanc.push_back(vtmp);
	}
}

void TruncatedTspline::SetLshapeProblem_XP(string fn)
{
	//read quad vtk
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++)
			getline(fin, stmp); //skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
			cp[i].act = 1;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	//set T-spline info
	InitialConnect_1();
	//VisualizeTMesh("test8/Lshape_Tmesh_1");
	int edge_zero[26] = {15, 17, 19, 26, 28, 53, 55, 62, 64, 157, 159, 166, 168, 212, 214, 221, 243, 247, 249, 255, 275, 279, 349, 355, 375, 379};
	//int el_type3[4]={5,74,96,175};
	for (int i = 0; i < 26; i++)
	{
		tmedge[edge_zero[i]].len = 0.;
	}
	//for(int i=0; i<4; i++)
	//{
	//	tmesh[el_type3[i]].type=3;
	//}
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmedge[tmesh[i].edge[0]].len == 0. || tmedge[tmesh[i].edge[1]].len == 0.)
		{
			tmesh[i].type = 2;
		}
		if (tmedge[tmesh[i].edge[0]].len == 0. && tmedge[tmesh[i].edge[1]].len == 0.)
		{
			tmesh[i].type = 3;
		}
	}
	FindEdgeTopoDirec_1();
	FindKnotInterval_1();
	UpdateKnotInterval_1();
	SetLocalCoorSystem();
	FindIEN_1();
}

//void TruncatedTspline::SetProblem_complex(string fn)//fertility
//{
//	//read quad vtk
//	string fname(fn+".vtk"),stmp;
//	int npts,neles,itmp;
//	ifstream fin;
//	fin.open(fname);
//	if(fin.is_open())
//	{
//		for(int i=0;i<4;i++) getline(fin,stmp);
//		fin>>stmp>>npts>>stmp;
//		cp.resize(npts);
//		for(int i=0;i<npts;i++)
//		{
//			fin>>cp[i].coor[0]>>cp[i].coor[1]>>cp[i].coor[2];
//			cp[i].act=1;
//		}
//		getline(fin,stmp);
//		fin>>stmp>>neles>>itmp;
//		tmesh.resize(neles);
//		for(int i=0;i<neles;i++)
//		{
//			fin>>itmp>>tmesh[i].cnct[0]>>tmesh[i].cnct[1]>>tmesh[i].cnct[2]>>tmesh[i].cnct[3];
//		}
//		fin.close();
//	}
//	else
//	{
//		cerr<<"Cannot open "<<fname<<"!\n";
//	}
//
//	//vector<array<double,3>> pts;
//	//vector<array<int,4>> cnct(tmesh.size());
//	//vector<int> id(cp.size());
//	//double tol(1.e-4);
//	//for(uint i=0; i<cp.size(); i++)
//	//{
//	//	id[i]=-1;
//	//	double dis;
//	//	for(uint j=0; j<pts.size(); j++)
//	//	{
//	//		dis=sqrt((cp[i].coor[0]-pts[j][0])*(cp[i].coor[0]-pts[j][0])+(cp[i].coor[1]-pts[j][1])*(cp[i].coor[1]-pts[j][1])+(cp[i].coor[2]-pts[j][2])*(cp[i].coor[2]-pts[j][2]));
//	//		if(dis<tol)
//	//		{
//	//			id[i]=j; break;
//	//		}
//	//	}
//	//	if(id[i]==-1)
//	//	{
//	//		array<double,3> tmp={cp[i].coor[0],cp[i].coor[1],cp[i].coor[2]};
//	//		id[i]=pts.size();
//	//		pts.push_back(tmp);
//	//	}
//	//}
//	//for(uint i=0; i<tmesh.size(); i++)
//	//{
//	//	for(int j=0; j<4; j++)
//	//		cnct[i][j]=id[tmesh[i].cnct[j]];
//	//}
//
//	//fname=fn+"1.vtk";
//	//ofstream fout;
//	//fout.open(fname.c_str());
//	//if(fout.is_open())
//	//{
//	//	fout<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
//	//	fout<<"POINTS "<<pts.size()<<" float\n";
//	//	for(uint i=0;i<pts.size();i++)
//	//	{
//	//		fout<<pts[i][0]<<" "<<pts[i][1]<<" "<<pts[i][2]<<"\n";
//	//	}
//	//	fout<<"\nCELLS "<<cnct.size()<<" "<<5*cnct.size()<<'\n';
//	//	for(uint i=0;i<cnct.size();i++)
//	//	{
//	//		fout<<"4 "<<cnct[i][0]<<" "<<cnct[i][1]<<" "<<cnct[i][2]<<" "<<cnct[i][3]<<'\n';
//	//	}
//	//	fout<<"\nCELL_TYPES "<<cnct.size()<<'\n';
//	//	for(uint i=0;i<cnct.size();i++)
//	//	{
//	//		fout<<"9\n";
//	//	}
//	//	fout.close();
//	//}
//	//else
//	//{
//	//	cout<<"Cannot open "<<fname<<"!\n";
//	//}
//
//	//for(uint i=0; i<tmesh.size(); i++)
//	//{
//	//	for(int j=0; j<4; j++)
//	//	{
//	//		cp[tmesh[i].cnct[j]].face.push_back(i);
//	//	}
//	//}
//	//uint nv_max(4),nv_min(4);
//	//for(uint i=0; i<cp.size(); i++)
//	//{
//	//	if(cp[i].face.size()>nv_max)
//	//	{
//	//		nv_max=cp[i].face.size();
//	//	}
//	//	if(cp[i].face.size()<nv_min)
//	//	{
//	//		nv_min=cp[i].face.size();
//	//	}
//	//}
//	//cout<<"Max and min valence:\n";
//	//cout<<nv_max<<" "<<nv_min<<"\n";
//	//int ninv(0);
//	for(uint i=0; i<tmesh.size(); i++)
//	{
//		int nxp(0);
//		for(int j=0; j<4; j++)
//		{
//			if(cp[tmesh[i].cnct[j]].face.size()==3 || cp[tmesh[i].cnct[j]].face.size()>4)
//			{
//				nxp++;
//			}
//		}
//		if(nxp>1)
//		{
//			tmesh[i].type=5;
//			//ninv++;
//		}
//	}
//	//cout<<"num of invalid: "<<ninv<<"\n";
//
//	//set T-spline info
//	InitialConnect_1();
//	//Validate_Tmesh();
//	FindEdgeTopoDirec_1();
//	FindKnotInterval_1();
//	UpdateKnotInterval_1();
//	SetLocalCoorSystem();
//	FindIEN_1();
//}

void TruncatedTspline::SetProblem_complex(string fn)
{
	//read quad vtk
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++)
			getline(fin, stmp); //skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
			cp[i].act = 1;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	//for(uint i=0; i<tmesh.size(); i++)
	//{
	//	for(int j=0; j<4; j++)
	//	{
	//		cp[tmesh[i].cnct[j]].face.push_back(i);
	//	}
	//}
	//uint nv_max(4),nv_min(4);
	//for(uint i=0; i<cp.size(); i++)
	//{
	//	if(cp[i].face.size()>nv_max)
	//	{
	//		nv_max=cp[i].face.size();
	//	}
	//	if(cp[i].face.size()<nv_min)
	//	{
	//		nv_min=cp[i].face.size();
	//	}
	//}
	//cout<<"Max and min valence:\n";
	//cout<<nv_max<<" "<<nv_min<<"\n";
	//int ninv(0);
	//for(uint i=0; i<tmesh.size(); i++)
	//{
	//	int nxp(0);
	//	for(int j=0; j<4; j++)
	//	{
	//		if(cp[tmesh[i].cnct[j]].face.size()==3 || cp[tmesh[i].cnct[j]].face.size()>4)
	//		{
	//			nxp++;
	//		}
	//	}
	//	if(nxp>1)
	//	{
	//		tmesh[i].type=5;
	//		ninv++;
	//	}
	//}
	//cout<<"num of invalid: "<<ninv<<"\n";

	//set T-spline info
	InitialConnect_2();
	FindEdgeTopoDirec_1();
	FindKnotInterval_1();
	UpdateKnotInterval_1();
	SetLocalCoorSystem();
	//FindIEN_1();
	FindIEN_3();
	Update_IEN_3();
	//Validate_Tmesh();
	//CollectActives();
}

//void TruncatedTspline::SetBounaryPoints(vector<int>& pid, vector<double>& disp, int& ncp)//bunny
//{
//	ncp=cp.size();
//	int nebc1(9);
//	int ebcid1[]={2332,1351,1350,758,2082,2081,759,2130,793};//disp 100
//	pid.clear();
//	disp.clear();
//	for(int i=0; i<nebc1; i++)
//	{
//		for(int j=0; j<4; j++)
//		{
//			vector<int>::iterator it=find(pid.begin(),pid.end(),tmesh[ebcid1[i]].cnct[j]);
//			if(it==pid.end())
//			{
//				pid.push_back(tmesh[ebcid1[i]].cnct[j]);
//				disp.push_back(100.);
//			}
//		}
//	}
//	int nebc2(16);
//	int ebcid2[]={2356,2746,2498,2497,1277,1174,1175,2714,1276,1173,1172,2093,1741,2183,2184,2092};//disp 0
//	for(int i=0; i<nebc2; i++)
//	{
//		for(int j=0; j<4; j++)
//		{
//			vector<int>::iterator it=find(pid.begin(),pid.end(),tmesh[ebcid2[i]].cnct[j]);
//			if(it==pid.end())
//			{
//				pid.push_back(tmesh[ebcid2[i]].cnct[j]);
//				disp.push_back(0.);
//			}
//		}
//	}
//}

//void TruncatedTspline::SetBounaryPoints(vector<int>& pid, vector<double>& disp, int& ncp)//hand
//{
//	ncp=cp.size();
//	int nebc1(16);
//	int ebcid1[]={1176,1177,1178,1534,1175,1029,1027,1024,1173,1028,1026,1025,1174,1179,3922,3920};//disp 100
//	pid.clear();
//	disp.clear();
//	for(int i=0; i<nebc1; i++)
//	{
//		for(int j=0; j<4; j++)
//		{
//			vector<int>::iterator it=find(pid.begin(),pid.end(),tmesh[ebcid1[i]].cnct[j]);
//			if(it==pid.end())
//			{
//				pid.push_back(tmesh[ebcid1[i]].cnct[j]);
//				disp.push_back(100.);
//			}
//		}
//	}
//	int nebc2(24);
//	int ebcid2[]={3690,3695,3725,3741,3745,3793,3691,3696,3726,3742,3746,3791,3724,3723,3727,3747,3739,3740,3772,3729,3728,3748,3736,3737};//disp 0
//	for(int i=0; i<nebc2; i++)
//	{
//		for(int j=0; j<4; j++)
//		{
//			vector<int>::iterator it=find(pid.begin(),pid.end(),tmesh[ebcid2[i]].cnct[j]);
//			if(it==pid.end())
//			{
//				pid.push_back(tmesh[ebcid2[i]].cnct[j]);
//				disp.push_back(0.);
//			}
//		}
//	}
//}

//void TruncatedTspline::SetBounaryPoints(vector<int>& pid, vector<double>& disp, int& ncp)//fertility
//{
//	ncp=cp.size();
//	int nebc1(9);
//	int ebcid1[]={293,294,295,441,442,274,277,279,292};//disp 100
//	pid.clear();
//	disp.clear();
//	for(int i=0; i<nebc1; i++)
//	{
//		if(tmesh[ebcid1[i]].act==0)
//		{
//			for(int k=0; k<4; k++)
//			{
//				if(tmesh[ebcid1[i]].chd[k]!=-1)
//				{
//					int echd(tmesh[ebcid1[i]].chd[k]);
//					for(int j=0; j<4; j++)
//					{
//						vector<int>::iterator it=find(pid.begin(),pid.end(),tmesh[echd].cnct[j]);
//						if(it==pid.end())
//						{
//							pid.push_back(tmesh[echd].cnct[j]);
//							disp.push_back(100.);
//						}
//					}
//				}
//			}
//		}
//		else
//		{
//			for(int j=0; j<4; j++)
//			{
//				vector<int>::iterator it=find(pid.begin(),pid.end(),tmesh[ebcid1[i]].cnct[j]);
//				if(it==pid.end())
//				{
//					pid.push_back(tmesh[ebcid1[i]].cnct[j]);
//					disp.push_back(100.);
//				}
//			}
//		}
//	}
//	int nebc2(6);
//	int ebcid2[]={341,343,344,345,346,347};//disp 0
//	for(int i=0; i<nebc2; i++)
//	{
//		if(tmesh[ebcid2[i]].act==0)
//		{
//			for(int k=0; k<4; k++)
//			{
//				if(tmesh[ebcid2[i]].chd[k]!=-1)
//				{
//					int echd(tmesh[ebcid2[i]].chd[k]);
//					for(int j=0; j<4; j++)
//					{
//						vector<int>::iterator it=find(pid.begin(),pid.end(),tmesh[echd].cnct[j]);
//						if(it==pid.end())
//						{
//							pid.push_back(tmesh[echd].cnct[j]);
//							disp.push_back(0.);
//						}
//					}
//				}
//			}
//		}
//		else
//		{
//			for(int j=0; j<4; j++)
//			{
//				vector<int>::iterator it=find(pid.begin(),pid.end(),tmesh[ebcid2[i]].cnct[j]);
//				if(it==pid.end())
//				{
//					pid.push_back(tmesh[ebcid2[i]].cnct[j]);
//					disp.push_back(0.);
//				}
//			}
//		}
//	}
//}

//void TruncatedTspline::Validate_Tmesh()//bunny
//{
//	vector<int> rid1;
//	Identify_Invalid_Elements(rid1);
//	int nsub4(2),nsub2(13);
//	int sub4[]={844,2506};//subdivide into 4
//	int sub2[]={843,66,67,223,858,1523,1948,568,1222,2198,1828,307,944};//subdivide into 2
//	int sub2_dir[]={523,198,18,60,261,1815,2331,163,401,2868,1930,1419,2561};//a point to identify direction
//	for(int i=0; i<nsub4; i++)
//	{
//		rid1.push_back(sub4[i]);
//	}
//	vector<int> rid2(sub2,sub2+nsub2);
//	vector<int> rid2_dir(nsub2,0);
//	for(int i=0; i<nsub2; i++)
//	{
//		int* it=find(tmesh[sub2[i]].cnct,tmesh[sub2[i]].cnct+4,sub2_dir[i]);
//		int loc(it-tmesh[sub2[i]].cnct);
//		if(loc==0 || loc==2)
//		{
//			rid2_dir[i]=0;
//		}
//		else if(loc==1 || loc==3)
//		{
//			rid2_dir[i]=1;
//		}
//		else
//		{
//			cerr<<"Cannot find refinement direction!\n";
//			getchar();
//		}
//	}
//
//	//for(uint i=0; i<rid1.size(); i++)
//	//{
//	//	ElementRefine_Unstruct_Topo_Geom_4(rid1[i]);
//	//}
//	//for(uint i=0; i<rid2.size(); i++)
//	//{
//	//	ElementRefine_Unstruct_Topo_Geom_2(rid2[i],rid2_dir[i]);
//	//}
//
//	vector<int> rid,ridtype,rid_more,ridtype_more;
//	for(uint i=0; i<rid1.size(); i++)
//	{
//		rid.push_back(rid1[i]);
//		if(tmesh[rid1[i]].type==4)
//		{
//			ridtype.push_back(4);
//		}
//		else if(tmesh[rid1[i]].type==5)
//		{
//			ridtype.push_back(5);
//		}
//		else
//		{
//			ridtype.push_back(0);
//		}
//	}
//	for(uint i=0; i<rid2.size(); i++)
//	{
//		rid.push_back(rid2[i]);
//		ridtype.push_back(rid2_dir[i]+1);
//	}
//	Topo_Refine_Unstruct(rid,ridtype,rid_more,ridtype_more);
//	Geom_Refine_Unstruct(rid_more,ridtype_more);
//	SetBezierMatIrrPatch();
//}

void TruncatedTspline::SetBounaryPoints(vector<int> &pid, vector<double> &disp, int &ncp) //horse
{
	ncp = cp.size();
	int nebc1(15);
	int ebcid1[] = {136, 10156, 10131, 10130, 32, 33, 13916, 10129, 13876, 13878, 135, 13877, 13917, 200, 13918}; //disp 100
	pid.clear();
	disp.clear();
	for (int i = 0; i < nebc1; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			vector<int>::iterator it = find(pid.begin(), pid.end(), tmesh[ebcid1[i]].cnct[j]);
			if (it == pid.end())
			{
				pid.push_back(tmesh[ebcid1[i]].cnct[j]);
				disp.push_back(100.);
			}
		}
	}
	int nebc2(16);
	int ebcid2[] = {13976, 13978, 5437, 5438, 3380, 12440, 5463, 5466, 3379, 12441, 5464, 5465, 3383, 12678, 6226, 6228}; //disp 0
	for (int i = 0; i < nebc2; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			vector<int>::iterator it = find(pid.begin(), pid.end(), tmesh[ebcid2[i]].cnct[j]);
			if (it == pid.end())
			{
				pid.push_back(tmesh[ebcid2[i]].cnct[j]);
				disp.push_back(0.);
			}
		}
	}
}

//void TruncatedTspline::Validate_Tmesh()//fertility
//{
//	vector<int> rid1;
//	Identify_Invalid_Elements(rid1);
//	int nsub4(1),nsub2(0);
//	int sub4[]={92};//subdivide into 4
//	int sub2[]={0};//subdivide into 2
//	int sub2_dir[]={0};//a point to identify direction
//	for(int i=0; i<nsub4; i++)
//	{
//		rid1.push_back(sub4[i]);
//	}
//	vector<int> rid2(sub2,sub2+nsub2);
//	vector<int> rid2_dir(nsub2,0);
//	for(int i=0; i<nsub2; i++)
//	{
//		int* it=find(tmesh[sub2[i]].cnct,tmesh[sub2[i]].cnct+4,sub2_dir[i]);
//		int loc(it-tmesh[sub2[i]].cnct);
//		if(loc==0 || loc==2)
//		{
//			rid2_dir[i]=0;
//		}
//		else if(loc==1 || loc==3)
//		{
//			rid2_dir[i]=1;
//		}
//		else
//		{
//			cerr<<"Cannot find refinement direction!\n";
//			getchar();
//		}
//	}
//
//	//for(uint i=0; i<rid1.size(); i++)
//	//{
//	//	ElementRefine_Unstruct_Topo_Geom_4(rid1[i]);
//	//}
//	//for(uint i=0; i<rid2.size(); i++)
//	//{
//	//	ElementRefine_Unstruct_Topo_Geom_2(rid2[i],rid2_dir[i]);
//	//}
//
//	vector<int> rid,ridtype,rid_more,ridtype_more;
//	for(uint i=0; i<rid1.size(); i++)
//	{
//		rid.push_back(rid1[i]);
//		if(tmesh[rid1[i]].type==4)
//		{
//			ridtype.push_back(4);
//		}
//		else if(tmesh[rid1[i]].type==5)
//		{
//			ridtype.push_back(5);
//		}
//		else
//		{
//			ridtype.push_back(0);
//		}
//	}
//	for(uint i=0; i<rid2.size(); i++)
//	{
//		rid.push_back(rid2[i]);
//		ridtype.push_back(rid2_dir[i]+1);
//	}
//	Topo_Refine_Unstruct(rid,ridtype,rid_more,ridtype_more);
//	Geom_Refine_Unstruct(rid_more,ridtype_more);
//	SetBezierMatIrrPatch();
//}

//void TruncatedTspline::Validate_Tmesh()//hand
//{
//	vector<int> rid1;
//	Identify_Invalid_Elements(rid1);
//	int nsub4(0),nsub2(20);
//	int sub4[]={0};//subdivide into 4
//	int sub2[]={2937,2936,2811,2813,786,2791,2326,917,3578,1918,2953,62,1911,1903,1898,2938,3511,1197,1431,2959};//subdivide into 2
//	int sub2_dir[]={2334,2336,2332,2333,277,2326,2073,321,3453,997,2346,17,3041,3040,982,2335,3690,447,3093,2356};//a point to identify direction
//	for(int i=0; i<nsub4; i++)
//	{
//		rid1.push_back(sub4[i]);
//	}
//	vector<int> rid2(sub2,sub2+nsub2);
//	vector<int> rid2_dir(nsub2,0);
//	for(int i=0; i<nsub2; i++)
//	{
//		int* it=find(tmesh[sub2[i]].cnct,tmesh[sub2[i]].cnct+4,sub2_dir[i]);
//		int loc(it-tmesh[sub2[i]].cnct);
//		if(loc==0 || loc==2)
//		{
//			rid2_dir[i]=0;
//		}
//		else if(loc==1 || loc==3)
//		{
//			rid2_dir[i]=1;
//		}
//		else
//		{
//			cerr<<"Cannot find refinement direction!\n";
//			getchar();
//		}
//	}
//
//	//for(uint i=0; i<rid1.size(); i++)
//	//{
//	//	ElementRefine_Unstruct_Topo_Geom_4(rid1[i]);
//	//}
//	//for(uint i=0; i<rid2.size(); i++)
//	//{
//	//	ElementRefine_Unstruct_Topo_Geom_2(rid2[i],rid2_dir[i]);
//	//}
//
//	vector<int> rid,ridtype,rid_more,ridtype_more;
//	for(uint i=0; i<rid1.size(); i++)
//	{
//		rid.push_back(rid1[i]);
//		if(tmesh[rid1[i]].type==4)
//		{
//			ridtype.push_back(4);
//		}
//		else if(tmesh[rid1[i]].type==5)
//		{
//			ridtype.push_back(5);
//		}
//		else
//		{
//			ridtype.push_back(0);
//		}
//	}
//	for(uint i=0; i<rid2.size(); i++)
//	{
//		rid.push_back(rid2[i]);
//		ridtype.push_back(rid2_dir[i]+1);
//	}
//	Topo_Refine_Unstruct(rid,ridtype,rid_more,ridtype_more);
//	Geom_Refine_Unstruct(rid_more,ridtype_more);
//	SetBezierMatIrrPatch();
//}

void TruncatedTspline::Validate_Tmesh() //horse
{
	vector<int> rid1;
	Identify_Invalid_Elements(rid1);
	int nsub4(0), nsub2(7);
	int sub4[] = {0};											   //subdivide into 4
	int sub2[] = {2167, 13750, 240, 10021, 6068, 6065, 9672};	   //subdivide into 2
	int sub2_dir[] = {6339, 13298, 5881, 3607, 11470, 2120, 3451}; //a point to identify direction
	for (int i = 0; i < nsub4; i++)
	{
		rid1.push_back(sub4[i]);
	}
	vector<int> rid2(sub2, sub2 + nsub2);
	vector<int> rid2_dir(nsub2, 0);
	for (int i = 0; i < nsub2; i++)
	{
		int *it = find(tmesh[sub2[i]].cnct, tmesh[sub2[i]].cnct + 4, sub2_dir[i]);
		int loc(it - tmesh[sub2[i]].cnct);
		if (loc == 0 || loc == 2)
		{
			rid2_dir[i] = 0;
		}
		else if (loc == 1 || loc == 3)
		{
			rid2_dir[i] = 1;
		}
		else
		{
			cerr << "Cannot find refinement direction!\n";
			getchar();
		}
	}

	//for(uint i=0; i<rid1.size(); i++)
	//{
	//	ElementRefine_Unstruct_Topo_Geom_4(rid1[i]);
	//}
	//for(uint i=0; i<rid2.size(); i++)
	//{
	//	ElementRefine_Unstruct_Topo_Geom_2(rid2[i],rid2_dir[i]);
	//}

	vector<int> rid, ridtype, rid_more, ridtype_more;
	for (uint i = 0; i < rid1.size(); i++)
	{
		rid.push_back(rid1[i]);
		if (tmesh[rid1[i]].type == 4)
		{
			ridtype.push_back(4);
		}
		else if (tmesh[rid1[i]].type == 5)
		{
			ridtype.push_back(5);
		}
		else
		{
			ridtype.push_back(0);
		}
	}
	for (uint i = 0; i < rid2.size(); i++)
	{
		rid.push_back(rid2[i]);
		ridtype.push_back(rid2_dir[i] + 1);
	}
	Topo_Refine_Unstruct(rid, ridtype, rid_more, ridtype_more);
	Geom_Refine_Unstruct(rid_more, ridtype_more);
	SetBezierMatIrrPatch();
}

void TruncatedTspline::SetProblem_PlateHole(string fn)
{
	//read quad vtk
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++)
			getline(fin, stmp); //skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
			cp[i].act = 1;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	double ku[11] = {0., 0., 0., 0., .5, .5, .5, 1., 1., 1., 1.};
	double kv[8] = {0., 0., 0., 0., 1., 1., 1., 1.};
	double nurbs_w[7] = {1., 0.9023689271, 0.8535533906, 0.8535533906, 0.8535533906, 0.9023689271, 1.};
	//double nurbs_w[7]={1.,0.9023689271,0.9023689271,1.,0.9023689271,0.9023689271,1.};
	int nu(7), nv(4);

	for (int i = 0; i < nu; i++)
		cp[i].w = nurbs_w[i];

	int loc(0);
	for (int j = 0; j < nv; j++)
	{
		for (int i = 0; i < nu; i++)
		{
			for (int k = 0; k < 5; k++)
			{
				cp[loc].knotU[k] = ku[i + k];
				cp[loc].knotV[k] = kv[j + k];
			}
			for (int k = 0; k < 4; k++)
			{
				cp[loc].kitvU[k] = cp[loc].knotU[k + 1] - cp[loc].knotU[k];
				cp[loc].kitvV[k] = cp[loc].knotV[k + 1] - cp[loc].knotV[k];
			}
			cp[loc].pm[0] = cp[loc].knotU[2];
			cp[loc].pm[1] = cp[loc].knotV[2];
			loc++;
		}
	}
	loc = 0;
	for (int j = 0; j < nv - 1; j++)
	{
		for (int i = 0; i < nu - 1; i++)
		{
			//if(i!=0 && i!=nu-2 && j!=0 && j!=nv-2)
			//{
			//	tmesh[loc].act=1;
			//	tmesh[loc].IEN.resize(16);
			//	int loc1(0);
			//	for(int j1=j-1; j1<j+3; j1++)
			//	{
			//		for(int i1=i-1; i1<i+3; i1++)
			//		{
			//			tmesh[loc].IEN[loc1]=j1*nu+i1;
			//			loc1++;
			//		}
			//	}
			//}
			double len_u(cp[tmesh[loc].cnct[1]].pm[0] - cp[tmesh[loc].cnct[0]].pm[0]), len_v(cp[tmesh[loc].cnct[3]].pm[1] - cp[tmesh[loc].cnct[0]].pm[1]);
			if (len_u * len_v != 0.)
			{
				tmesh[loc].act = 1;
				tmesh[loc].IEN.resize(16);
				int loc1(0);
				for (int j1 = j - 1; j1 < j + 3; j1++)
				{
					for (int i1 = i - 1; i1 < i + 3; i1++)
					{
						tmesh[loc].IEN[loc1] = j1 * nu + i1;
						loc1++;
					}
				}
			}
			loc++;
		}
	}

	loc = 0;
	for (int j = 0; j < nv; j++)
	{
		for (int i = 0; i < nu; i++)
		{
			cp[loc].index[0] = i;
			cp[loc].index[1] = j;
			loc++;
		}
	}
	uanc.clear();
	vanc.clear();
	for (int i = 0; i < nu; i++)
	{
		pair<int, double> utmp(i, ku[i + 2]);
		uanc.push_back(utmp);
	}
	for (int i = 0; i < nv; i++)
	{
		pair<int, double> vtmp(i, kv[i + 2]);
		vanc.push_back(vtmp);
	}
}

void TruncatedTspline::SetProblem_PlateHole_XP(string fn)
{
	//read quad vtk
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++)
			getline(fin, stmp); //skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
			cp[i].act = 1;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	double nurbs_w[9] = {1., 0.9674563090, 0.9132168241, 0.8644012876, 0.8481294421, 0.8644012876, 0.9132168241, 0.9674563090, 1.};
	for (int i = 0; i < 9; i++)
		cp[i].w = nurbs_w[i];

	InitialConnect_1();
	FindEdgeTopoDirec_1();
	FindKnotInterval_1();
	UpdateKnotInterval_1();
	SetLocalCoorSystem();
	FindIEN_1();
}

//void TruncatedTspline::Validate_Tmesh()//fertility
//{
//	vector<int> rid1;
//	//Identify_Invalid_Elements(rid1);
//	int nsub4(14),nsub2(24);
//	int sub4[]={428,429,430,431,520,521,522,434,435,504,505,506,507,12};//subdivide into 4
//	int sub2[]={151,150,135,134,72,73,74,75,249,248,372,373,374,375,24,25,14,96,98,13,276,277,26,97};//subdivide into 2
//	int sub2_dir[]={320,319,318,317,174,161,202,203,474,473,273,278,291,293,23,11,12,78,79,20,47,48,75,74};//a point to identify direction
//	for(int i=0; i<nsub4; i++)
//	{
//		rid1.push_back(sub4[i]);
//	}
//	Identify_More_Elements(rid1);
//	vector<int> rid2(sub2,sub2+nsub2);
//	vector<int> rid2_dir(nsub2,0);
//	for(int i=0; i<nsub2; i++)
//	{
//		int* it=find(tmesh[sub2[i]].cnct,tmesh[sub2[i]].cnct+4,sub2_dir[i]);
//		int loc(it-tmesh[sub2[i]].cnct);
//		if(loc==0 || loc==2)
//		{
//			rid2_dir[i]=0;
//		}
//		else if(loc==1 || loc==3)
//		{
//			rid2_dir[i]=1;
//		}
//		else
//		{
//			cerr<<"Cannot find refinement direction!\n";
//			cout<<sub2[i]<<"\n";
//			getchar();
//		}
//	}
//
//	//for(uint i=0; i<rid1.size(); i++)
//	//{
//	//	ElementRefine_Unstruct_Topo_Geom_4(rid1[i]);
//	//}
//	//for(uint i=0; i<rid2.size(); i++)
//	//{
//	//	ElementRefine_Unstruct_Topo_Geom_2(rid2[i],rid2_dir[i]);
//	//}
//
//	vector<int> rid,ridtype,rid_more,ridtype_more;
//	for(uint i=0; i<rid1.size(); i++)
//	{
//		rid.push_back(rid1[i]);
//		if(tmesh[rid1[i]].type==4)
//		{
//			ridtype.push_back(4);
//		}
//		else if(tmesh[rid1[i]].type==5)
//		{
//			ridtype.push_back(5);
//		}
//		else
//		{
//			ridtype.push_back(0);
//		}
//	}
//	for(uint i=0; i<rid2.size(); i++)
//	{
//		rid.push_back(rid2[i]);
//		ridtype.push_back(rid2_dir[i]+1);
//	}
//	Topo_Refine_Unstruct(rid,ridtype,rid_more,ridtype_more);
//	Geom_Refine_Unstruct(rid_more,ridtype_more);
//
//	CollectActives();
//	int sub4_1[]={1149,1148,1502};//subdivide into 4
//	//int sub4_1[]={1148};//subdivide into 4
//	int n1(1);
//	vector<int> rfid1(n1);
//	for(int i=0; i<n1; i++)
//	{
//		rfid1[i]=eaid[sub4_1[i]];
//	}
//	vector<int> rfidtype1(n1,0);
//	vector<int> rfid_more1, rfidtype_more1;
//	Topo_Refine_Unstruct(rfid1,rfidtype1,rfid_more1,rfidtype_more1);
//	Geom_Refine_Unstruct(rfid_more1,rfidtype_more1);
//
//	SetBezierMatIrrPatch();
//}

void TruncatedTspline::Identify_Invalid_Elements(vector<int> &rid)
{
	//find two ring elements of XPs
	vector<int> xpid;
	vector<vector<int>> e2r;
	vector<int> loc(cp.size(), -1);
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].face.size() == 3 || cp[i].face.size() > 4)
		{
			loc[i] = count++;
			xpid.push_back(i);
			vector<int> etmp;
			vector<int> p1r;
			for (uint j = 0; j < cp[i].face.size(); j++)
			{
				etmp.push_back(cp[i].face[j]);
				for (uint k = 0; k < 4; k++)
				{
					if (tmesh[cp[i].face[j]].cnct[k] != i)
					{
						vector<int>::iterator it = find(p1r.begin(), p1r.end(), tmesh[cp[i].face[j]].cnct[k]);
						if (it == p1r.end())
						{
							p1r.push_back(tmesh[cp[i].face[j]].cnct[k]);
						}
					}
				}
			}
			for (uint j = 0; j < p1r.size(); j++)
			{
				for (uint k = 0; k < cp[p1r[j]].face.size(); k++)
				{
					vector<int>::iterator it = find(etmp.begin(), etmp.end(), cp[p1r[j]].face[k]);
					if (it == etmp.end())
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
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].type == 5)
		{
			for (int j = 0; j < 4; j++)
			{
				if (cp[tmesh[i].cnct[j]].face.size() == 3 || cp[tmesh[i].cnct[j]].face.size() > 4)
				{
					vector<int>::iterator it = find(xpid0.begin(), xpid0.end(), loc[tmesh[i].cnct[j]]);
					if (it == xpid0.end())
					{
						xpid0.push_back(loc[tmesh[i].cnct[j]]);
					}
				}
			}
		}
	}
	vector<int> rid0;
	for (uint i = 0; i < xpid0.size(); i++)
	{
		for (uint j = 0; j < e2r[xpid0[i]].size(); j++)
		{
			vector<int>::iterator it = find(rid0.begin(), rid0.end(), e2r[xpid0[i]][j]);
			if (it == rid0.end())
			{
				rid0.push_back(e2r[xpid0[i]][j]);
			}
		}
	}
	for (uint i = 0; i < rid0.size(); i++)
	{
		tmesh[rid0[i]].act = 0;
	}
	for (int iter = 0; iter < 3; iter++)
	{
		for (uint i = 0; i < e2r.size(); i++)
		{
			int flag(0);
			for (uint j = 0; j < e2r[i].size(); j++)
			{
				if (tmesh[e2r[i][j]].act == 0)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 1)
			{
				for (uint j = 0; j < e2r[i].size(); j++)
				{
					vector<int>::iterator it = find(rid0.begin(), rid0.end(), e2r[i][j]);
					if (it == rid0.end())
					{
						rid0.push_back(e2r[i][j]);
						tmesh[e2r[i][j]].act = 0;
					}
				}
			}
		}
	}

	rid.clear();
	rid = rid0;
	rid0.clear();
}

void TruncatedTspline::Identify_More_Elements(vector<int> &rid)
{
	//find two ring elements of XPs
	vector<int> xpid;
	vector<vector<int>> e2r;
	vector<int> loc(cp.size(), -1);
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].face.size() == 3 || cp[i].face.size() > 4)
		{
			loc[i] = count++;
			xpid.push_back(i);
			vector<int> etmp;
			vector<int> p1r;
			for (uint j = 0; j < cp[i].face.size(); j++)
			{
				etmp.push_back(cp[i].face[j]);
				for (uint k = 0; k < 4; k++)
				{
					if (tmesh[cp[i].face[j]].cnct[k] != i)
					{
						vector<int>::iterator it = find(p1r.begin(), p1r.end(), tmesh[cp[i].face[j]].cnct[k]);
						if (it == p1r.end())
						{
							p1r.push_back(tmesh[cp[i].face[j]].cnct[k]);
						}
					}
				}
			}
			for (uint j = 0; j < p1r.size(); j++)
			{
				for (uint k = 0; k < cp[p1r[j]].face.size(); k++)
				{
					vector<int>::iterator it = find(etmp.begin(), etmp.end(), cp[p1r[j]].face[k]);
					if (it == etmp.end())
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
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].type == 5)
		{
			for (int j = 0; j < 4; j++)
			{
				if (cp[tmesh[i].cnct[j]].face.size() == 3 || cp[tmesh[i].cnct[j]].face.size() > 4)
				{
					vector<int>::iterator it = find(xpid0.begin(), xpid0.end(), loc[tmesh[i].cnct[j]]);
					if (it == xpid0.end())
					{
						xpid0.push_back(loc[tmesh[i].cnct[j]]);
					}
				}
			}
		}
	}
	vector<int> rid0;
	for (uint i = 0; i < xpid0.size(); i++)
	{
		for (uint j = 0; j < e2r[xpid0[i]].size(); j++)
		{
			vector<int>::iterator it = find(rid0.begin(), rid0.end(), e2r[xpid0[i]][j]);
			if (it == rid0.end())
			{
				rid0.push_back(e2r[xpid0[i]][j]);
			}
		}
	}
	for (uint i = 0; i < rid.size(); i++)
	{
		rid0.push_back(rid[i]);
	}
	for (uint i = 0; i < rid0.size(); i++)
	{
		tmesh[rid0[i]].act = 0;
	}
	for (int iter = 0; iter < 6; iter++)
	{
		for (uint i = 0; i < e2r.size(); i++)
		{
			int flag(0);
			for (uint j = 0; j < e2r[i].size(); j++)
			{
				if (tmesh[e2r[i][j]].act == 0)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 1)
			{
				for (uint j = 0; j < e2r[i].size(); j++)
				{
					vector<int>::iterator it = find(rid0.begin(), rid0.end(), e2r[i][j]);
					if (it == rid0.end())
					{
						rid0.push_back(e2r[i][j]);
						tmesh[e2r[i][j]].act = 0;
					}
				}
			}
		}
	}

	rid.clear();
	rid = rid0;
	rid0.clear();
}

void TruncatedTspline::ElementRefine_Unstruct_Topo_Geom_4(int eid)
{
	int pid[5];
	int edid[12];
	Vertex ptmp1;
	//ptmp1.update=1;
	ptmp1.coor[0] = (cp[tmesh[eid].cnct[0]].coor[0] + cp[tmesh[eid].cnct[1]].coor[0] + cp[tmesh[eid].cnct[2]].coor[0] + cp[tmesh[eid].cnct[3]].coor[0]) / 4.;
	ptmp1.coor[1] = (cp[tmesh[eid].cnct[0]].coor[1] + cp[tmesh[eid].cnct[1]].coor[1] + cp[tmesh[eid].cnct[2]].coor[1] + cp[tmesh[eid].cnct[3]].coor[1]) / 4.;
	ptmp1.coor[2] = (cp[tmesh[eid].cnct[0]].coor[2] + cp[tmesh[eid].cnct[1]].coor[2] + cp[tmesh[eid].cnct[2]].coor[2] + cp[tmesh[eid].cnct[3]].coor[2]) / 4.;
	cp.push_back(ptmp1);
	pid[0] = cp.size() - 1;
	for (int j = 0; j < 4; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[j], tmesh[eid].cnct[(j + 1) % 4]};
		if (tmedge[tmesh[eid].edge[j]].act == 1)
		{
			Vertex ptmp;
			//ptmp.update=1;
			int id0[4] = {tmesh[eid].cnct[j], tmesh[eid].cnct[(j + 1) % 4], tmesh[eid].cnct[(j + 2) % 4], tmesh[eid].cnct[(j + 3) % 4]};
			ptmp.coor[0] = (cp[id0[0]].coor[0] + cp[id0[1]].coor[0]) / 2.;
			ptmp.coor[1] = (cp[id0[0]].coor[1] + cp[id0[1]].coor[1]) / 2.;
			ptmp.coor[2] = (cp[id0[0]].coor[2] + cp[id0[1]].coor[2]) / 2.;
			//ptmp.coortmp[0]=(6.*cp[id0[0]].coor[0]+6.*cp[id0[1]].coor[0]+cp[id0[2]].coor[0]+cp[id0[3]].coor[0])/16.;
			//ptmp.coortmp[1]=(6.*cp[id0[0]].coor[1]+6.*cp[id0[1]].coor[1]+cp[id0[2]].coor[1]+cp[id0[3]].coor[1])/16.;
			//ptmp.coortmp[2]=(6.*cp[id0[0]].coor[2]+6.*cp[id0[1]].coor[2]+cp[id0[2]].coor[2]+cp[id0[3]].coor[2])/16.;
			//{
			//	ptmp.coor[0]=cp[tmesh[eid].cnct[0]].coor[0]/16.;
			//	ptmp.coor[1]=cp[tmesh[eid].cnct[0]].coor[1]/16.;
			//	ptmp.coor[2]=cp[tmesh[eid].cnct[0]].coor[2]/16.;
			//}
			cp.push_back(ptmp);
			pid[j + 1] = cp.size() - 1;
			tmedge[tmesh[eid].edge[j]].midpt = pid[j + 1];
			Edge edtmp1, edtmp2, edtmp3;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j + 1];
			edtmp1.len = tmedge[tmesh[eid].edge[j]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[j];
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j + 1];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[j]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[j];
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[(j + 1) % 4]].len / 2.;
			tmedge.push_back(edtmp1);
			edid[3 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].chd[0] = edid[3 * j];
			tmedge.push_back(edtmp2);
			edid[3 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].chd[1] = edid[3 * j + 1];
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].act = 0;
		}
		else
		{
			pid[j + 1] = tmedge[tmesh[eid].edge[j]].midpt;
			int id1[2] = {tmesh[eid].cnct[(j + 2) % 4], tmesh[eid].cnct[(j + 3) % 4]};
			//cp[pid[j+1]].update=1;
			//if(j==0 || j==3)
			//{
			//cp[pid[j+1]].coortmp[0]+=(cp[id1[0]].coor[0]+cp[id1[1]].coor[0])/16.;
			//cp[pid[j+1]].coortmp[1]+=(cp[id1[0]].coor[1]+cp[id1[1]].coor[1])/16.;
			//cp[pid[j+1]].coortmp[2]+=(cp[id1[0]].coor[2]+cp[id1[1]].coor[2])/16.;
			//}
			//else
			//{
			//	cp[pid[j+1]].coor[0]+=cp[tmesh[eid].cnct[0]].coor[0]/16.;
			//	cp[pid[j+1]].coor[1]+=cp[tmesh[eid].cnct[0]].coor[1]/16.;
			//	cp[pid[j+1]].coor[2]+=cp[tmesh[eid].cnct[0]].coor[2]/16.;
			//}
			Edge edtmp3;
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[(j + 1) % 4]].len / 2.;
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			int ied(tmedge[tmesh[eid].edge[j]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[3 * j] = ied;
				edid[3 * j + 1] = tmedge[tmesh[eid].edge[j]].chd[1];
			}
			else
			{
				edid[3 * j] = tmedge[tmesh[eid].edge[j]].chd[1];
				edid[3 * j + 1] = ied;
			}
		}
	}

	//int nvl(cp[tmesh[eid].cnct[0]].face.size());
	//if(cp[tmesh[eid].cnct[0]].update==0)
	//{
	//	double ccf[3]={1.-7./(4.*nvl),3./(2.*nvl*nvl),1./(4.*nvl*nvl)};
	//	cp[tmesh[eid].cnct[0]].coortmp[0]=ccf[0]*cp[tmesh[eid].cnct[0]].coor[0]+ccf[1]*cp[tmesh[eid].cnct[1]].coor[0]+ccf[2]*cp[tmesh[eid].cnct[2]].coor[0];
	//	cp[tmesh[eid].cnct[0]].coortmp[1]=ccf[0]*cp[tmesh[eid].cnct[0]].coor[1]+ccf[1]*cp[tmesh[eid].cnct[1]].coor[1]+ccf[2]*cp[tmesh[eid].cnct[2]].coor[1];
	//	cp[tmesh[eid].cnct[0]].coortmp[2]=ccf[0]*cp[tmesh[eid].cnct[0]].coor[2]+ccf[1]*cp[tmesh[eid].cnct[1]].coor[2]+ccf[2]*cp[tmesh[eid].cnct[2]].coor[2];
	//	cp[tmesh[eid].cnct[0]].update=2;
	//}
	//else if(cp[tmesh[eid].cnct[0]].update==2)
	//{
	//	double ccf[2]={3./(2.*nvl*nvl),1./(4.*nvl*nvl)};
	//	cp[tmesh[eid].cnct[0]].coortmp[0]+=ccf[0]*cp[tmesh[eid].cnct[1]].coor[0]+ccf[1]*cp[tmesh[eid].cnct[2]].coor[0];
	//	cp[tmesh[eid].cnct[0]].coortmp[1]+=ccf[0]*cp[tmesh[eid].cnct[1]].coor[1]+ccf[1]*cp[tmesh[eid].cnct[2]].coor[1];
	//	cp[tmesh[eid].cnct[0]].coortmp[2]+=ccf[0]*cp[tmesh[eid].cnct[1]].coor[2]+ccf[1]*cp[tmesh[eid].cnct[2]].coor[2];
	//}
	//else
	//{
	//	cout<<"Not supported for other udpate types!\n";
	//	getchar();
	//}

	//for(int j=0; j<4; j++)
	//{
	//	int nvl(cp[tmesh[eid].cnct[j]].face.size());
	//	if(cp[tmesh[eid].cnct[j]].update==0)
	//	{
	//		cp[tmesh[eid].cnct[j]].update=1;
	//		double ccf[3]={1.-7./(4.*nvl),3./(2.*nvl*nvl),1./(4.*nvl*nvl)};
	//		cp[tmesh[eid].cnct[j]].coortmp[0]=ccf[0]*cp[tmesh[eid].cnct[j]].coor[0]+ccf[1]*cp[tmesh[eid].cnct[(j+1)%4]].coor[0]+ccf[2]*cp[tmesh[eid].cnct[(j+2)%4]].coor[0];
	//		cp[tmesh[eid].cnct[j]].coortmp[1]=ccf[0]*cp[tmesh[eid].cnct[j]].coor[1]+ccf[1]*cp[tmesh[eid].cnct[(j+1)%4]].coor[1]+ccf[2]*cp[tmesh[eid].cnct[(j+2)%4]].coor[1];
	//		cp[tmesh[eid].cnct[j]].coortmp[2]=ccf[0]*cp[tmesh[eid].cnct[j]].coor[2]+ccf[1]*cp[tmesh[eid].cnct[(j+1)%4]].coor[2]+ccf[2]*cp[tmesh[eid].cnct[(j+2)%4]].coor[2];
	//	}
	//	else
	//	{
	//		double ccf[2]={3./(2.*nvl*nvl),1./(4.*nvl*nvl)};
	//		cp[tmesh[eid].cnct[j]].coortmp[0]+=ccf[0]*cp[tmesh[eid].cnct[(j+1)%4]].coor[0]+ccf[1]*cp[tmesh[eid].cnct[(j+2)%4]].coor[0];
	//		cp[tmesh[eid].cnct[j]].coortmp[1]+=ccf[0]*cp[tmesh[eid].cnct[(j+1)%4]].coor[1]+ccf[1]*cp[tmesh[eid].cnct[(j+2)%4]].coor[1];
	//		cp[tmesh[eid].cnct[j]].coortmp[2]+=ccf[0]*cp[tmesh[eid].cnct[(j+1)%4]].coor[2]+ccf[1]*cp[tmesh[eid].cnct[(j+2)%4]].coor[2];
	//	}
	//}

	int e_cnct[4][4] = {{tmesh[eid].cnct[0], pid[1], pid[0], pid[4]}, {pid[1], tmesh[eid].cnct[1], pid[2], pid[0]}, {pid[0], pid[2], tmesh[eid].cnct[2], pid[3]}, {pid[4], pid[0], pid[3], tmesh[eid].cnct[3]}};
	int e_edge[4][4] = {{edid[0], edid[2], edid[11], edid[10]}, {edid[1], edid[3], edid[5], edid[2]}, {edid[5], edid[4], edid[6], edid[8]}, {edid[11], edid[8], edid[7], edid[9]}};
	int enewid[4];
	int e_type[4] = {4, 0, 0, 0};
	vector<Element> etmp(4);
	for (int i = 0; i < 4; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = e_type[i];
		etmp[i].lv = tmesh[eid].lv + 1.;
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[j] = e_cnct[i][j];
			etmp[i].edge[j] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::ElementRefine_Unstruct_Topo_Geom_2(int eid, int dir)
{
	int pid[2], edid[7], pos(dir);
	int cnid[2] = {pos, (pos + 2) % 4};
	for (int j = 0; j < 2; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[cnid[j]], tmesh[eid].cnct[(cnid[j] + 1) % 4]};
		if (tmedge[tmesh[eid].edge[cnid[j]]].act == 1)
		{
			Vertex ptmp;
			ptmp.coor[0] = (cp[itmp[0]].coor[0] + cp[itmp[1]].coor[0]) / 2.;
			ptmp.coor[1] = (cp[itmp[0]].coor[1] + cp[itmp[1]].coor[1]) / 2.;
			ptmp.coor[2] = (cp[itmp[0]].coor[2] + cp[itmp[1]].coor[2]) / 2.;
			cp.push_back(ptmp);
			pid[j] = cp.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].midpt = pid[j];
			Edge edtmp1, edtmp2;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j];
			edtmp1.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[cnid[j]];
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[cnid[j]]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[cnid[j]];
			tmedge.push_back(edtmp1);
			edid[2 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[0] = edid[2 * j];
			tmedge.push_back(edtmp2);
			edid[2 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[cnid[j]]].chd[1] = edid[2 * j + 1];
			tmedge[tmesh[eid].edge[cnid[j]]].act = 0;
		}
		else
		{
			pid[j] = tmedge[tmesh[eid].edge[cnid[j]]].midpt;
			int ied(tmedge[tmesh[eid].edge[cnid[j]]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[2 * j] = ied;
				edid[2 * j + 1] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
			}
			else
			{
				edid[2 * j] = tmedge[tmesh[eid].edge[cnid[j]]].chd[1];
				edid[2 * j + 1] = ied;
			}
		}
	}
	Edge edtmp;
	edtmp.act = 1;
	edtmp.pt[0] = pid[0];
	edtmp.pt[1] = pid[1];
	edtmp.len = tmedge[tmesh[eid].edge[(pos + 1) % 4]].len;
	tmedge.push_back(edtmp);
	edid[4] = tmedge.size() - 1;
	edid[5] = tmesh[eid].edge[(pos + 1) % 4];
	edid[6] = tmesh[eid].edge[(pos + 3) % 4];

	int e_cnct[2][4] = {{tmesh[eid].cnct[pos], pid[0], pid[1], tmesh[eid].cnct[(pos + 3) % 4]}, {pid[0], tmesh[eid].cnct[(pos + 1) % 4], tmesh[eid].cnct[(pos + 2) % 4], pid[1]}};
	int e_edge[2][4] = {{edid[0], edid[4], edid[3], edid[6]}, {edid[1], edid[5], edid[2], edid[4]}};
	int enewid[2];
	double chd_org[2][2] = {{0., 0.}, {tmedge[tmesh[eid].edge[0]].len / 2., 0.}};
	if (dir == 1)
	{
		chd_org[1][0] = 0.;
		chd_org[1][1] = tmedge[tmesh[eid].edge[1]].len / 2.;
	}
	vector<Element> etmp(2);
	for (int i = 0; i < 2; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = 0;
		etmp[i].lv = tmesh[eid].lv + 0.5;
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[(pos + j) % 4] = e_cnct[i][j];
			etmp[i].edge[(pos + j) % 4] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
		tmesh[eid].chd_o[i][0] = chd_org[i][0];
		tmesh[eid].chd_o[i][1] = chd_org[i][1];
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::ElementRefine_Invalid_4(int eid)
{
	int pid[5];
	int edid[12];
	Vertex ptmp1;
	ptmp1.update = 1;
	ptmp1.coortmp[0] = (cp[tmesh[eid].cnct[0]].coor[0] + cp[tmesh[eid].cnct[1]].coor[0] + cp[tmesh[eid].cnct[2]].coor[0] + cp[tmesh[eid].cnct[3]].coor[0]) / 4.;
	ptmp1.coortmp[1] = (cp[tmesh[eid].cnct[0]].coor[1] + cp[tmesh[eid].cnct[1]].coor[1] + cp[tmesh[eid].cnct[2]].coor[1] + cp[tmesh[eid].cnct[3]].coor[1]) / 4.;
	ptmp1.coortmp[2] = (cp[tmesh[eid].cnct[0]].coor[2] + cp[tmesh[eid].cnct[1]].coor[2] + cp[tmesh[eid].cnct[2]].coor[2] + cp[tmesh[eid].cnct[3]].coor[2]) / 4.;
	cp.push_back(ptmp1);
	pid[0] = cp.size() - 1;
	for (int j = 0; j < 4; j++)
	{
		int itmp[2] = {tmesh[eid].cnct[j], tmesh[eid].cnct[(j + 1) % 4]};
		if (tmedge[tmesh[eid].edge[j]].act == 1)
		{
			Vertex ptmp;
			int id0[4] = {tmesh[eid].cnct[j], tmesh[eid].cnct[(j + 1) % 4], tmesh[eid].cnct[(j + 2) % 4], tmesh[eid].cnct[(j + 3) % 4]};
			int fcnb(tmedge[tmesh[eid].edge[j]].face[0]);
			if (fcnb == eid)
				fcnb = tmedge[tmesh[eid].edge[j]].face[1];
			int *it1 = find(tmesh[fcnb].cnct, tmesh[fcnb].cnct + 4, tmesh[eid].cnct[j]);
			int loc1(it1 - tmesh[fcnb].cnct);
			int id1[2] = {tmesh[fcnb].cnct[(loc1 + 2) % 4], tmesh[fcnb].cnct[(loc1 + 3) % 4]};
			ptmp.coortmp[0] = (6. * cp[id0[0]].coor[0] + 6. * cp[id0[1]].coor[0] + cp[id0[2]].coor[0] + cp[id0[3]].coor[0] + cp[id1[0]].coor[0] + cp[id1[1]].coor[0]) / 16.;
			ptmp.coortmp[1] = (6. * cp[id0[0]].coor[1] + 6. * cp[id0[1]].coor[1] + cp[id0[2]].coor[1] + cp[id0[3]].coor[1] + cp[id1[0]].coor[1] + cp[id1[1]].coor[1]) / 16.;
			ptmp.coortmp[2] = (6. * cp[id0[0]].coor[2] + 6. * cp[id0[1]].coor[2] + cp[id0[2]].coor[2] + cp[id0[3]].coor[2] + cp[id1[0]].coor[2] + cp[id1[1]].coor[2]) / 16.;
			//if(cp[itmp[0]].type==2 || cp[itmp[1]].type==2)
			//{
			//	ptmp.coortmp[0]=(6.*cp[id0[0]].coor[0]+6.*cp[id0[1]].coor[0]+cp[id0[2]].coor[0]+cp[id0[3]].coor[0])/16.;
			//	ptmp.coortmp[1]=(6.*cp[id0[0]].coor[1]+6.*cp[id0[1]].coor[1]+cp[id0[2]].coor[1]+cp[id0[3]].coor[1])/16.;
			//	ptmp.coortmp[2]=(6.*cp[id0[0]].coor[2]+6.*cp[id0[1]].coor[2]+cp[id0[2]].coor[2]+cp[id0[3]].coor[2])/16.;
			//	ptmp.update=2;
			//}
			//else
			//{
			//	ptmp.coor[0]=cp[tmesh[eid].cnct[0]].coor[0]/16.;
			//	ptmp.coor[1]=cp[tmesh[eid].cnct[0]].coor[1]/16.;
			//	ptmp.coor[2]=cp[tmesh[eid].cnct[0]].coor[2]/16.;
			//}
			cp.push_back(ptmp);
			pid[j + 1] = cp.size() - 1;
			tmedge[tmesh[eid].edge[j]].midpt = pid[j + 1];
			Edge edtmp1, edtmp2, edtmp3;
			edtmp1.act = 1;
			edtmp1.pt[0] = itmp[0];
			edtmp1.pt[1] = pid[j + 1];
			edtmp1.len = tmedge[tmesh[eid].edge[j]].len / 2.;
			edtmp1.prt = tmesh[eid].edge[j];
			edtmp2.act = 1;
			edtmp2.pt[0] = pid[j + 1];
			edtmp2.pt[1] = itmp[1];
			edtmp2.len = tmedge[tmesh[eid].edge[j]].len / 2.;
			edtmp2.prt = tmesh[eid].edge[j];
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[(j + 1) % 4]].len / 2.;
			tmedge.push_back(edtmp1);
			edid[3 * j] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].chd[0] = edid[3 * j];
			tmedge.push_back(edtmp2);
			edid[3 * j + 1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].chd[1] = edid[3 * j + 1];
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[j]].act = 0;
		}
		else
		{
			pid[j + 1] = tmedge[tmesh[eid].edge[j]].midpt;
			int id1[2] = {tmesh[eid].cnct[(j + 2) % 4], tmesh[eid].cnct[(j + 3) % 4]};
			//cp[pid[j+1]].update=1;
			//if(cp[pid[j+1]].update==2)
			//{
			//	cp[pid[j+1]].coortmp[0]+=(cp[id1[0]].coor[0]+cp[id1[1]].coor[0])/16.;
			//	cp[pid[j+1]].coortmp[1]+=(cp[id1[0]].coor[1]+cp[id1[1]].coor[1])/16.;
			//	cp[pid[j+1]].coortmp[2]+=(cp[id1[0]].coor[2]+cp[id1[1]].coor[2])/16.;
			//}
			//else
			//{
			//	cp[pid[j+1]].coor[0]+=cp[tmesh[eid].cnct[0]].coor[0]/16.;
			//	cp[pid[j+1]].coor[1]+=cp[tmesh[eid].cnct[0]].coor[1]/16.;
			//	cp[pid[j+1]].coor[2]+=cp[tmesh[eid].cnct[0]].coor[2]/16.;
			//}
			Edge edtmp3;
			edtmp3.act = 1;
			edtmp3.pt[0] = pid[0];
			edtmp3.pt[1] = pid[j + 1];
			edtmp3.len = tmedge[tmesh[eid].edge[(j + 1) % 4]].len / 2.;
			tmedge.push_back(edtmp3);
			edid[3 * j + 2] = tmedge.size() - 1;
			int ied(tmedge[tmesh[eid].edge[j]].chd[0]);
			if (tmedge[ied].pt[0] == itmp[0] || tmedge[ied].pt[1] == itmp[0])
			{
				edid[3 * j] = ied;
				edid[3 * j + 1] = tmedge[tmesh[eid].edge[j]].chd[1];
			}
			else
			{
				edid[3 * j] = tmedge[tmesh[eid].edge[j]].chd[1];
				edid[3 * j + 1] = ied;
			}
		}
	}

	for (int j = 0; j < 4; j++)
	{
		int nvl(cp[tmesh[eid].cnct[j]].face.size());
		if (nvl == 3 || nvl > 4)
		{
			if (cp[tmesh[eid].cnct[j]].update == 0)
			{
				double ccf[3] = {1. - 7. / (4. * nvl), 3. / (2. * nvl * nvl), 1. / (4. * nvl * nvl)};
				cp[tmesh[eid].cnct[j]].coortmp[0] = ccf[0] * cp[tmesh[eid].cnct[j]].coor[0] + ccf[1] * cp[tmesh[eid].cnct[(j + 1) % 4]].coor[0] + ccf[2] * cp[tmesh[eid].cnct[(j + 2) % 4]].coor[0];
				cp[tmesh[eid].cnct[j]].coortmp[1] = ccf[0] * cp[tmesh[eid].cnct[j]].coor[1] + ccf[1] * cp[tmesh[eid].cnct[(j + 1) % 4]].coor[1] + ccf[2] * cp[tmesh[eid].cnct[(j + 2) % 4]].coor[1];
				cp[tmesh[eid].cnct[j]].coortmp[2] = ccf[0] * cp[tmesh[eid].cnct[j]].coor[2] + ccf[1] * cp[tmesh[eid].cnct[(j + 1) % 4]].coor[2] + ccf[2] * cp[tmesh[eid].cnct[(j + 2) % 4]].coor[2];
				cp[tmesh[eid].cnct[j]].update = 2;
			}
			else if (cp[tmesh[eid].cnct[j]].update == 2)
			{
				double ccf[2] = {3. / (2. * nvl * nvl), 1. / (4. * nvl * nvl)};
				cp[tmesh[eid].cnct[j]].coortmp[0] += ccf[0] * cp[tmesh[eid].cnct[(j + 1) % 4]].coor[0] + ccf[1] * cp[tmesh[eid].cnct[(j + 2) % 4]].coor[0];
				cp[tmesh[eid].cnct[j]].coortmp[1] += ccf[0] * cp[tmesh[eid].cnct[(j + 1) % 4]].coor[1] + ccf[1] * cp[tmesh[eid].cnct[(j + 2) % 4]].coor[1];
				cp[tmesh[eid].cnct[j]].coortmp[2] += ccf[0] * cp[tmesh[eid].cnct[(j + 1) % 4]].coor[2] + ccf[1] * cp[tmesh[eid].cnct[(j + 2) % 4]].coor[2];
			}
			else
			{
				cout << "Not supported for other udpate types!\n";
				getchar();
			}
		}
		else if (nvl == 4)
		{
			if (cp[tmesh[eid].cnct[j]].update == 0)
			{
				int cpid(tmesh[eid].cnct[j]);
				double ccf[3] = {1. - 7. / (4. * nvl), 3. / (2. * nvl * nvl), 1. / (4. * nvl * nvl)};
				cp[cpid].coortmp[0] = ccf[0] * cp[cpid].coor[0];
				cp[cpid].coortmp[1] = ccf[0] * cp[cpid].coor[1];
				cp[cpid].coortmp[2] = ccf[0] * cp[cpid].coor[2];
				for (uint k = 0; k < cp[cpid].edge.size(); k++)
				{
					int pend = tmedge[cp[cpid].edge[k]].pt[0];
					if (pend == cpid)
						pend = tmedge[cp[cpid].edge[k]].pt[1];
					cp[cpid].coortmp[0] += ccf[1] * cp[pend].coor[0];
					cp[cpid].coortmp[1] += ccf[1] * cp[pend].coor[1];
					cp[cpid].coortmp[2] += ccf[1] * cp[pend].coor[2];
				}
				for (uint k = 0; k < cp[cpid].face.size(); k++)
				{
					int fcid(cp[cpid].face[k]);
					int *itfc = find(tmesh[fcid].cnct, tmesh[fcid].cnct + 4, cpid);
					int loctmp = itfc - tmesh[fcid].cnct;
					int pend(tmesh[fcid].cnct[(loctmp + 2) % 4]);
					cp[cpid].coortmp[0] += ccf[2] * cp[pend].coor[0];
					cp[cpid].coortmp[1] += ccf[2] * cp[pend].coor[1];
					cp[cpid].coortmp[2] += ccf[2] * cp[pend].coor[2];
				}
				cp[tmesh[eid].cnct[j]].update = 3;
			}
		}
	}

	int e_cnct[4][4] = {{tmesh[eid].cnct[0], pid[1], pid[0], pid[4]}, {pid[1], tmesh[eid].cnct[1], pid[2], pid[0]}, {pid[0], pid[2], tmesh[eid].cnct[2], pid[3]}, {pid[4], pid[0], pid[3], tmesh[eid].cnct[3]}};
	int e_edge[4][4] = {{edid[0], edid[2], edid[11], edid[10]}, {edid[1], edid[3], edid[5], edid[2]}, {edid[5], edid[4], edid[6], edid[8]}, {edid[11], edid[8], edid[7], edid[9]}};
	int enewid[4];
	int e_type[4] = {4, 0, 0, 0};
	for (int i = 1; i < 4; i++)
	{
		if (cp[tmesh[eid].cnct[i]].face.size() == 3 || cp[tmesh[eid].cnct[i]].face.size() > 4)
		{
			e_type[i] = 4;
			int cnct_tmp[4] = {e_cnct[i][0], e_cnct[i][1], e_cnct[i][2], e_cnct[i][3]};
			int edge_tmp[4] = {e_edge[i][0], e_edge[i][1], e_edge[i][2], e_edge[i][3]};
			for (int j = 0; j < 4; j++)
			{
				e_cnct[i][j] = cnct_tmp[(i + j) % 4];
				e_edge[i][j] = edge_tmp[(i + j) % 4];
			}
		}
	}
	//if(eid==2892)
	//{
	//	cout<<tmesh[eid].cnct[0]<<" "<<tmesh[eid].cnct[1]<<" "<<tmesh[eid].cnct[2]<<" "<<tmesh[eid].cnct[3]<<"\n";
	//	cout<<e_type[0]<<" "<<e_type[1]<<" "<<e_type[2]<<" "<<e_type[3]<<"\n";
	//	getchar();
	//}
	vector<Element> etmp(4);
	for (int i = 0; i < 4; i++)
	{
		etmp[i].act = 1;
		etmp[i].type = e_type[i];
		etmp[i].lv = tmesh[eid].lv + 1.;
		etmp[i].prt = eid;
		for (int j = 0; j < 4; j++)
		{
			etmp[i].cnct[j] = e_cnct[i][j];
			etmp[i].edge[j] = e_edge[i][j];
		}
		tmesh.push_back(etmp[i]);
		enewid[i] = tmesh.size() - 1;
		tmesh[eid].chd[i] = enewid[i];
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::Update_IEN_3()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		vector<int>().swap(tmesh[eid].IEN);
		vector<array<double, 5>>().swap(tmesh[eid].patch_ku);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kv);
		tmesh[eid].IEN = tmesh[eid].IENtmp;
		tmesh[eid].patch_ku = tmesh[eid].patch_kutmp;
		tmesh[eid].patch_kv = tmesh[eid].patch_kvtmp;
		vector<int>().swap(tmesh[eid].IENtmp);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kutmp);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kvtmp);
	}
}

void TruncatedTspline::runXP_complex(string fn)
{
	SetProblem_complex(fn);
	Validate_Tmesh();
	//CollectActives();

	//VisualizeControlMesh("complex/fertility_CM_5");
}

void TruncatedTspline::SurfRefine()
{
	CollectActives();
	int rid_tmp1[18] = {5, 4, 3, 17, 16, 15, 29, 28, 41, 6, 7, 8, 18, 19, 20, 30, 31, 42};
	int nsub1(18);
	vector<int> rid1(nsub1);
	for (int i = 0; i < nsub1; i++)
	{
		rid1[i] = eaid[rid_tmp1[i]];
	}
	Refine_Target(rid1);

	CollectActives();
	int rid_tmp2[18] = {51, 50, 55, 52, 53, 56, 63, 62, 64, 86, 87, 90, 89, 88, 93, 98, 99, 101};
	int nsub2(18);
	vector<int> rid2(nsub2);
	for (int i = 0; i < nsub2; i++)
	{
		rid2[i] = eaid[rid_tmp2[i]];
	}
	Refine_Addition(rid2);
	Refine_Target(rid2);

	CollectActives();
	int rid_tmp3[18] = {109, 108, 113, 110, 111, 114, 121, 120, 122, 144, 145, 148, 147, 146, 151, 156, 157, 159};
	int nsub3(18);
	vector<int> rid3(nsub3);
	for (int i = 0; i < nsub3; i++)
	{
		rid3[i] = eaid[rid_tmp3[i]];
	}
	Refine_Addition(rid3);
	Refine_Target(rid3);

	CollectActives();
	int rid_tmp4[18] = {167, 166, 171, 168, 169, 172, 179, 178, 180, 202, 203, 206, 205, 204, 209, 214, 215, 217};
	int nsub4(18);
	vector<int> rid4(nsub4);
	for (int i = 0; i < nsub4; i++)
	{
		rid4[i] = eaid[rid_tmp4[i]];
	}
	Refine_Addition(rid4);
	Refine_Target(rid4);

	CollectActives();
	int rid_tmp5[20] = {225, 224, 229, 226, 227, 230, 237, 236, 238, 260, 261, 264, 263, 262, 267, 272, 273, 275, 38, 49};
	int nsub5(20);
	vector<int> rid5(nsub5);
	for (int i = 0; i < nsub5; i++)
	{
		rid5[i] = eaid[rid_tmp5[i]];
	}
	Refine_Addition(rid5);
	Refine_Target(rid5);

	CollectActives();
}

void TruncatedTspline::SurfRefine_1()
{
	npt_old = cp.size();
	nel_old = tmesh.size();
	int rfid[4] = {44, 45, 54, 55};
	int rfdir[4] = {0, 1, 1, 0};
	int nrf(4);
	for (int i = 0; i < nrf; i++)
	{
		ElementRefine_Square_2(rfid[i], rfdir[i]);
	}
	UpdateTopology();
	FindLocalKnotVectors();
	GeometryRefine();

	//for (uint i = 0; i < cp.size(); i++)
	//{
	//	if (cp[i].trun == 1)
	//	{
	//		cout << "id: " << i << "\nchild:\n";
	//		for (uint j = 0; j < cp[i].tbf.size(); j++)
	//		{
	//			cout << cp[i].tbf[j] << " ";
	//		}
	//		cout << "\n";
	//		for (uint j = 0; j < cp[i].tc.size(); j++)
	//		{
	//			cout << cp[i].tc[j] << " ";
	//		}
	//		cout << "\n";
	//		getchar();
	//	}
	//}

	//for (uint i = nel_old; i < tmesh.size(); i++)
	//{
	//	cout << "id: " << i << "\nIEN:\n";
	//	for (uint j = 0; j < tmesh[i].IEN.size(); j++)
	//	{
	//		cout << tmesh[i].IEN[j] << " ";
	//	}
	//	cout << "\n";
	//	getchar();
	//}

	//int pid(38);
	//for (int i = 0; i < 5; i++)
	//{
	//	cout << cp[pid].knotU[i] << " ";
	//}
	//cout << "\n";
	//for (int i = 0; i < 5; i++)
	//{
	//	cout << cp[pid].knotV[i] << " ";
	//}
	//cout << "\n";
	//cout << "refine done!\n";
	//getchar();
}

void TruncatedTspline::LshapeRefine()
{
	int rid_tmp1[26] = {20, 21, 22, 36, 37, 38, 52, 53, 54, 68, 69, 70, 97,
						25, 26, 27, 41, 42, 43, 57, 58, 59, 73, 74, 75, 110};
	vector<int> rid1(rid_tmp1, rid_tmp1 + 26);
	Refine_Addition(rid1);
	Refine_Target(rid1);

	int rid_tmp2[22] = {133, 136, 137, 134, 139, 138, 145, 148, 149, 151, 150,
						180, 181, 184, 183, 182, 187, 192, 193, 196, 195, 194};
	vector<int> rid2(rid_tmp2, rid_tmp2 + 22);
	Refine_Addition(rid2);
	Refine_Target(rid2);

	int rid_tmp3[14] = {273, 276, 277, 279, 278, 288, 289,
						312, 313, 316, 315, 314, 324, 325};
	vector<int> rid3(rid_tmp3, rid_tmp3 + 14);
	Refine_Addition(rid3);
	Refine_Target(rid3);

	int rid_tmp4[10] = {396, 397, 399, 398, 405, 416, 417, 419, 418, 428};
	vector<int> rid4(rid_tmp4, rid_tmp4 + 10);
	Refine_Addition(rid4);
	Refine_Target(rid4);

	CollectActives();

	//cout<<tmesh[123].type<<'\n';
	//int ntjc=0;
	//for(int i=0; i<4; i++)
	//{
	//	if(tmedge[tmesh[123].edge[i]].act==0)
	//	{
	//		cout<<i<<'\n';
	//		ntjc++;
	//	}
	//}
	//cout<<ntjc<<'\n';
	//getchar();
}

void TruncatedTspline::LshapeRefine_H1()
{
	CollectActives();
	int rid_tmp1[18] = {5, 4, 3, 17, 16, 15, 29, 28, 41, 6, 7, 8, 18, 19, 20, 30, 31, 42};
	int nsub1(18);
	vector<int> rid1(nsub1);
	for (int i = 0; i < nsub1; i++)
	{
		rid1[i] = eaid[rid_tmp1[i]];
	}
	Refine_Addition(rid1);
	Refine_Target(rid1);

	CollectActives();
	int rid_tmp2[18] = {51, 50, 55, 52, 53, 56, 63, 62, 64, 86, 87, 90, 89, 88, 93, 98, 99, 101};
	int nsub2(18);
	vector<int> rid2(nsub2);
	for (int i = 0; i < nsub2; i++)
	{
		rid2[i] = eaid[rid_tmp2[i]];
	}
	Refine_Addition(rid2);
	Refine_Target(rid2);

	CollectActives();
	int rid_tmp3[18] = {109, 108, 113, 110, 111, 114, 121, 120, 122, 144, 145, 148, 147, 146, 151, 156, 157, 159};
	int nsub3(18);
	vector<int> rid3(nsub3);
	for (int i = 0; i < nsub3; i++)
	{
		rid3[i] = eaid[rid_tmp3[i]];
	}
	Refine_Addition(rid3);
	Refine_Target(rid3);

	CollectActives();
	int rid_tmp4[18] = {167, 166, 171, 168, 169, 172, 179, 178, 180, 202, 203, 206, 205, 204, 209, 214, 215, 217};
	int nsub4(18);
	vector<int> rid4(nsub4);
	for (int i = 0; i < nsub4; i++)
	{
		rid4[i] = eaid[rid_tmp4[i]];
	}
	Refine_Addition(rid4);
	Refine_Target(rid4);

	CollectActives();
	int rid_tmp5[20] = {225, 224, 229, 226, 227, 230, 237, 236, 238, 260, 261, 264, 263, 262, 267, 272, 273, 275, 38, 49};
	int nsub5(20);
	vector<int> rid5(nsub5);
	for (int i = 0; i < nsub5; i++)
	{
		rid5[i] = eaid[rid_tmp5[i]];
	}
	Refine_Addition(rid5);
	Refine_Target(rid5);

	CollectActives();
}

void TruncatedTspline::LshapeRefine_XP()
{
	int rid_tmp1[42] = {7, 2, 3, 8, 13, 12, 11, 14, 15, 42, 43, 41, 40, 38, 39, 16, 19, 18, 17, 23, 20,
						106, 107, 110, 105, 104, 109, 102, 103, 98, 142, 143, 141, 140, 130, 131, 117, 118, 116, 119, 113, 114};
	vector<int> rid1(rid_tmp1, rid_tmp1 + 42), rid1_more;
	vector<int> rftype1(42, 0), rftype1_more;
	for (int i = 0; i < 42; i++)
	{
		if (tmesh[rid_tmp1[i]].type == 4)
		{
			rftype1[i] = 4;
		}
	}
	Topo_Refine_Unstruct(rid1, rftype1, rid1_more, rftype1_more);
	Geom_Refine_Unstruct(rid1_more, rftype1_more);

	CollectActives();
	int rid_tmp2[22] = {79, 78, 83, 80, 81, 84, 91, 90, 95, 92, 93, 194, 195, 190, 197, 196, 193, 182, 183, 178, 185, 184}; //active element id, need to convert
	int rid_tmp2_1[2] = {96, 181};
	int rid_tmp2_dir[2] = {238, 311};
	int nsub2(22), nsub2_1(2);
	for (int i = 0; i < nsub2_1; i++)
	{
		int *it = find(tmesh[eaid[rid_tmp2_1[i]]].cnct, tmesh[eaid[rid_tmp2_1[i]]].cnct + 4, rid_tmp2_dir[i]);
		int loc(it - tmesh[eaid[rid_tmp2_1[i]]].cnct);
		if (loc == 0 || loc == 2)
			rid_tmp2_dir[i] = 0;
		else if (loc == 1 || loc == 3)
			rid_tmp2_dir[i] = 1;
		else
		{
			cout << "Cannot find refinement direction!\n";
			getchar();
		}
	}
	vector<int> rid2(nsub2), rid2_more;
	vector<int> rftype2(nsub2, 0), rftype2_more;
	for (int i = 0; i < nsub2; i++)
	{
		rid2[i] = eaid[rid_tmp2[i]];
		if (tmesh[eaid[rid_tmp2[i]]].type == 4)
		{
			rftype2[i] = 4;
		}
	}
	for (int i = 0; i < nsub2_1; i++)
	{
		rid2.push_back(eaid[rid_tmp2_1[i]]);
		rftype2.push_back(rid_tmp2_dir[i] + 1);
	}
	Topo_Refine_Unstruct(rid2, rftype2, rid2_more, rftype2_more);
	Geom_Refine_Unstruct(rid2_more, rftype2_more);

	CollectActives();
	int rid_tmp3[22] = {223, 222, 227, 224, 225, 228, 235, 234, 239, 236, 237, 266, 267, 270, 269, 268, 273, 278, 279, 282, 281, 280}; //active element id, need to convert
	int rid_tmp3_1[2] = {240, 285};
	int rid_tmp3_dir[2] = {395, 432};
	int nsub3(22), nsub3_1(2);
	for (int i = 0; i < nsub3_1; i++)
	{
		int *it = find(tmesh[eaid[rid_tmp3_1[i]]].cnct, tmesh[eaid[rid_tmp3_1[i]]].cnct + 4, rid_tmp3_dir[i]);
		int loc(it - tmesh[eaid[rid_tmp3_1[i]]].cnct);
		if (loc == 0 || loc == 2)
			rid_tmp3_dir[i] = 0;
		else if (loc == 1 || loc == 3)
			rid_tmp3_dir[i] = 1;
		else
		{
			cout << "Cannot find refinement direction!\n";
			getchar();
		}
	}
	vector<int> rid3(nsub3), rid3_more;
	vector<int> rftype3(nsub3, 0), rftype3_more;
	for (int i = 0; i < nsub3; i++)
	{
		rid3[i] = eaid[rid_tmp3[i]];
		if (tmesh[eaid[rid_tmp3[i]]].type == 4)
		{
			rftype3[i] = 4;
		}
	}
	for (int i = 0; i < nsub3_1; i++)
	{
		rid3.push_back(eaid[rid_tmp3_1[i]]);
		rftype3.push_back(rid_tmp3_dir[i] + 1);
	}
	Topo_Refine_Unstruct(rid3, rftype3, rid3_more, rftype3_more);
	Geom_Refine_Unstruct(rid3_more, rftype3_more);

	CollectActives();
	int rid_tmp4[16] = {291, 290, 295, 292, 293, 296, 303, 302, 334, 335, 338, 337, 336, 341, 346, 347}; //active element id, need to convert
	int rid_tmp4_1[2] = {307, 350};
	int rid_tmp4_dir[2] = {378, 524};
	int nsub4(16), nsub4_1(2);
	for (int i = 0; i < nsub4_1; i++)
	{
		int *it = find(tmesh[eaid[rid_tmp4_1[i]]].cnct, tmesh[eaid[rid_tmp4_1[i]]].cnct + 4, rid_tmp4_dir[i]);
		int loc(it - tmesh[eaid[rid_tmp4_1[i]]].cnct);
		if (loc == 0 || loc == 2)
			rid_tmp4_dir[i] = 0;
		else if (loc == 1 || loc == 3)
			rid_tmp4_dir[i] = 1;
		else
		{
			cout << "Cannot find refinement direction!\n";
			getchar();
		}
	}
	vector<int> rid4(nsub4), rid4_more;
	vector<int> rftype4(nsub4, 0), rftype4_more;
	for (int i = 0; i < nsub4; i++)
	{
		rid4[i] = eaid[rid_tmp4[i]];
		if (tmesh[eaid[rid_tmp4[i]]].type == 4)
		{
			rftype4[i] = 4;
		}
	}
	for (int i = 0; i < nsub4_1; i++)
	{
		rid4.push_back(eaid[rid_tmp4_1[i]]);
		rftype4.push_back(rid_tmp4_dir[i] + 1);
	}
	Topo_Refine_Unstruct(rid4, rftype4, rid4_more, rftype4_more);
	Geom_Refine_Unstruct(rid4_more, rftype4_more);

	//CollectActives();
}

void TruncatedTspline::RepairMesh()
{
	int pid[8] = {58, 411, 67, 72, 57, 66, 465, 530};
	int npid(8);
	vector<array<double, 3>> norm(npid);
	vector<double> mag(npid);
	double scale(1.);
	for (int i = 0; i < npid; i++)
	{
		double nmavg[3] = {0., 0., 0.};
		for (uint j = 0; j < cp[pid[i]].face.size(); j++)
		{
			int eid(cp[pid[i]].face[j]);
			int *it = find(tmesh[eid].cnct, tmesh[eid].cnct + 4, pid[i]);
			if (it != tmesh[eid].cnct + 4)
			{
				int loc(it - tmesh[eid].cnct);
				double nm1[3], nm2[3], nm3[3];
				for (int k = 0; k < 3; k++)
				{
					nm1[k] = cp[tmesh[eid].cnct[(loc + 1) % 4]].coor[k] - cp[tmesh[eid].cnct[loc]].coor[k];
					nm2[k] = cp[tmesh[eid].cnct[(loc + 3) % 4]].coor[k] - cp[tmesh[eid].cnct[loc]].coor[k];
				}
				nm3[0] = nm1[1] * nm2[2] - nm1[2] * nm2[1];
				nm3[1] = nm1[2] * nm2[0] - nm1[0] * nm2[2];
				nm3[2] = nm1[0] * nm2[1] - nm1[1] * nm2[0];
				double dist = sqrt(nm3[0] * nm3[0] + nm3[1] * nm3[1] + nm3[2] * nm3[2]);
				nm3[0] /= dist;
				nm3[1] /= dist;
				nm3[2] /= dist;
				nmavg[0] += nm3[0];
				nmavg[1] += nm3[1];
				nmavg[2] += nm3[2];
			}
			else
			{
				cerr << "Cannot find pid in this element!\n";
				getchar();
			}
		}
		norm[i][0] = nmavg[0] / cp[pid[i]].face.size();
		norm[i][1] = nmavg[1] / cp[pid[i]].face.size();
		norm[i][2] = nmavg[2] / cp[pid[i]].face.size();
		mag[i] = 1.e6;
		for (uint j = 0; j < cp[pid[i]].edge.size(); j++)
		{
			int edid(cp[pid[i]].edge[j]);
			double dir[3] = {cp[tmedge[edid].pt[1]].coor[0] - cp[tmedge[edid].pt[0]].coor[0], cp[tmedge[edid].pt[1]].coor[1] - cp[tmedge[edid].pt[0]].coor[1],
							 cp[tmedge[edid].pt[1]].coor[2] - cp[tmedge[edid].pt[0]].coor[2]};
			double dist = scale * sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
			if (mag[i] > dist)
				mag[i] = dist;
		}
	}
	for (int i = 0; i < npid; i++)
	{
		cp[pid[i]].coor[0] += norm[i][0] * mag[i];
		cp[pid[i]].coor[1] += norm[i][1] * mag[i];
		cp[pid[i]].coor[2] += norm[i][2] * mag[i];
	}
}

void TruncatedTspline::Refine_PlateHole()
{
	CollectActives();
	int rid_tmp1[65] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
						35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 61, 64, 65, 115, 116, 109, 110, 125, 126};
	int nsub1(65);
	vector<int> rid1(nsub1);
	for (int i = 0; i < nsub1; i++)
	{
		rid1[i] = eaid[rid_tmp1[i]];
	}
	Refine_Addition(rid1);
	Refine_Target(rid1);

	CollectActives();
	int rid_tmp2[98] = {54, 55, 58, 59, 62, 63, 66, 67, 70, 71, 74, 75, 78, 79, 82, 83, 86, 87, 90, 91, 94, 95, 98, 99, 102, 103, 106, 107, 110, 111, 114, 115,
						57, 56, 61, 60, 65, 64, 69, 76, 81, 80, 85, 84, 89, 88, 93, 92, 97, 96, 101, 100, 105, 104, 109, 108, 113, 112, 117, 116, 142, 143, 146, 147, 150, 151, 154, 155, 158, 159, //70
						118, 119, 162, 163, 166, 167, 170, 171, 174, 175, 178, 179, 135, 138, 139, 141, 140, 145, 144, 149, 148, 153, 152, 157, 156, 161, 160, 165};
	int nsub2(99);
	vector<int> rid2(nsub2);
	for (int i = 0; i < nsub2; i++)
	{
		rid2[i] = eaid[rid_tmp2[i]];
	}
	Refine_Addition(rid2);
	Refine_Target(rid2);

	CollectActives();
	int rid_tmp3[47] = {221, 222, 225, 226, 229, 230, 233, 234, 237, 238, 241, 242, 245, 265, 266, 269, 270, 273, 274, 277, 277, 281, 282, 285, 286, 289, 290, 293, 294, 297, 298,
						301, 302, 305, 306, 309, 310, 313, 314, 317, 318, 321, 322, 325, 326, 329, 330};
	int nsub3(47);
	vector<int> rid3(nsub3);
	for (int i = 0; i < nsub3; i++)
	{
		rid3[i] = eaid[rid_tmp3[i]];
	}
	Refine_Addition(rid3);
	Refine_Target(rid3);
}

void TruncatedTspline::Refine_PlateHole_XP()
{
	//int rid_tmp1[42]={7,2,3,8,13,12,11,14,15,42,43,41,40,38,39,16,19,18,17,23,20,
	//	106,107,110,105,104,109,102,103,98,142,143,141,140,130,131,117,118,116,119,113,114};
	//vector<int> rid1(rid_tmp1,rid_tmp1+42), rid1_more;
	//vector<int> rftype1(42,0), rftype1_more;
	//for(int i=0; i<42; i++)
	//{
	//	if(tmesh[rid_tmp1[i]].type==4)
	//	{
	//		rftype1[i]=4;
	//	}
	//}
	//Topo_Refine_Unstruct(rid1,rftype1,rid1_more,rftype1_more);
	//Geom_Refine_Unstruct(rid1_more,rftype1_more);

	CollectActives();
	int rid_tmp2[95] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 28, 29, 32, 33, 36, 37, 40, 41, 44, 45,								 //36
						27, 26, 31, 30, 35, 34, 39, 38, 43, 52, 53, 58, 59, 61, 62, 64, 55, 54, 57, 56, 60, 63, 67, 76, 77, 83, 80, 144, 108, 109, 112, 78, 82, 81, 147, 146, 145, 111, 110, //39
						93, 156, 148, 120, 96, 97, 100, 101, 104, 105, 168, 169, 172, 173, 176, 177, 94, 159, 121, 115};																	 //active element id, need to convert
	int rid_tmp2_1[12] = {48, 49, 42, 79, 92, 95, 157, 158, 47, 46, 124, 151};
	int rid_tmp2_dir[12] = {138, 29, 150, 181, 48, 193, 250, 258, 33, 153, 66, 85};
	int nsub2(95), nsub2_1(12);
	for (int i = 0; i < nsub2_1; i++)
	{
		int *it = find(tmesh[eaid[rid_tmp2_1[i]]].cnct, tmesh[eaid[rid_tmp2_1[i]]].cnct + 4, rid_tmp2_dir[i]);
		int loc(it - tmesh[eaid[rid_tmp2_1[i]]].cnct);
		if (loc == 0 || loc == 2)
			rid_tmp2_dir[i] = 0;
		else if (loc == 1 || loc == 3)
			rid_tmp2_dir[i] = 1;
		else
		{
			cout << "Cannot find refinement direction!\n";
			getchar();
		}
	}
	vector<int> rid2(nsub2), rid2_more;
	vector<int> rftype2(nsub2, 0), rftype2_more;
	for (int i = 0; i < nsub2; i++)
	{
		rid2[i] = eaid[rid_tmp2[i]];
		if (tmesh[eaid[rid_tmp2[i]]].type == 4)
		{
			rftype2[i] = 4;
		}
	}
	for (int i = 0; i < nsub2_1; i++)
	{
		rid2.push_back(eaid[rid_tmp2_1[i]]);
		rftype2.push_back(rid_tmp2_dir[i] + 1);
	}
	Topo_Refine_Unstruct(rid2, rftype2, rid2_more, rftype2_more);
	Geom_Refine_Unstruct(rid2_more, rftype2_more);

	CollectActives();
	int rid_tmp3[94] = {73, 74, 77, 78, 89, 90, 106, 109, 110, 121, 122, 125, 126, 137, 138, 141, 142, 153, 154, 157, 158, 91, 96, 95, 108, 107, 112, 111, 124, 123, 128, 127, 140, 139, 144, 143, 156, 155, 160, 159, //40
						85, 86, 81, 82, 88, 87, 84, 83, 104, 103, 100, 99, 120, 119, 116, 115, 136, 135, 132, 131, 174, 177, 178, 181, 182, 185, 186, 189, 190, 193, 194, 197, 198, 201,							   //75
						293, 297, 329, 325, 321, 295, 296, 298, 299, 300, 330, 331, 332, 326, 327, 328, 322, 323, 324, 294};
	int rid_tmp3_1[6] = {93, 92, 105, 101, 152, 173};
	int rid_tmp3_dir[6] = {104, 288, 12, 116, 334, 117};
	int nsub3(94), nsub3_1(6);
	for (int i = 0; i < nsub3_1; i++)
	{
		int *it = find(tmesh[eaid[rid_tmp3_1[i]]].cnct, tmesh[eaid[rid_tmp3_1[i]]].cnct + 4, rid_tmp3_dir[i]);
		int loc(it - tmesh[eaid[rid_tmp3_1[i]]].cnct);
		if (loc == 0 || loc == 2)
			rid_tmp3_dir[i] = 0;
		else if (loc == 1 || loc == 3)
			rid_tmp3_dir[i] = 1;
		else
		{
			cout << "Cannot find refinement direction!\n";
			getchar();
		}
	}
	vector<int> rid3(nsub3), rid3_more;
	vector<int> rftype3(nsub3, 0), rftype3_more;
	for (int i = 0; i < nsub3; i++)
	{
		rid3[i] = eaid[rid_tmp3[i]];
		if (tmesh[eaid[rid_tmp3[i]]].type == 4)
		{
			rftype3[i] = 4;
		}
	}
	for (int i = 0; i < nsub3_1; i++)
	{
		rid3.push_back(eaid[rid_tmp3_1[i]]);
		rftype3.push_back(rid_tmp3_dir[i] + 1);
	}
	Topo_Refine_Unstruct(rid3, rftype3, rid3_more, rftype3_more);
	Geom_Refine_Unstruct(rid3_more, rftype3_more);

	CollectActives();
	int rid_tmp4[20] = {673, 677, 681, 685, 689, 675, 676, 678, 679, 680, 682, 683, 684, 686, 687, 688, 690, 691, 692, 674}; //active element id, need to convert
	int rid_tmp4_1[2] = {307, 350};
	int rid_tmp4_dir[2] = {378, 524};
	int nsub4(20), nsub4_1(0);
	for (int i = 0; i < nsub4_1; i++)
	{
		int *it = find(tmesh[eaid[rid_tmp4_1[i]]].cnct, tmesh[eaid[rid_tmp4_1[i]]].cnct + 4, rid_tmp4_dir[i]);
		int loc(it - tmesh[eaid[rid_tmp4_1[i]]].cnct);
		if (loc == 0 || loc == 2)
			rid_tmp4_dir[i] = 0;
		else if (loc == 1 || loc == 3)
			rid_tmp4_dir[i] = 1;
		else
		{
			cout << "Cannot find refinement direction!\n";
			getchar();
		}
	}
	vector<int> rid4(nsub4), rid4_more;
	vector<int> rftype4(nsub4, 0), rftype4_more;
	for (int i = 0; i < nsub4; i++)
	{
		rid4[i] = eaid[rid_tmp4[i]];
		if (tmesh[eaid[rid_tmp4[i]]].type == 4)
		{
			rftype4[i] = 4;
		}
	}
	for (int i = 0; i < nsub4_1; i++)
	{
		rid4.push_back(eaid[rid_tmp4_1[i]]);
		rftype4.push_back(rid_tmp4_dir[i] + 1);
	}
	Topo_Refine_Unstruct(rid4, rftype4, rid4_more, rftype4_more);
	Geom_Refine_Unstruct(rid4_more, rftype4_more);
}

bool CompareAnchor(const pair<int, double> &anc1, const pair<int, double> &anc2)
{
	//if((anc1.first==9 && anc2.first==10) || (anc1.first==10 && anc2.first==9))
	//{
	//	//double dif(anc1.second-anc2.second);
	//	//cout<<dif<<'\n';
	//	bool tmp(anc1.second < anc2.second);
	//	cout<<tmp<<'\n';
	//	getchar();
	//}
	//double tol(1.e-6);
	//double dif(anc1.second-anc2.second);
	//if(dif<0.) dif=-dif;
	//if(dif<tol) return false;

	//if(anc1.second==anc1.second)
	//{
	//	return anc1.first<anc2.first;
	//}
	//else
	//{
	//	return anc1.second < anc2.second;
	//}

	return (anc1.second < anc2.second) || (anc1.second == anc2.second && anc1.first < anc2.first);
}

void TruncatedTspline::set_PatchTest_HighOrder(string fn)
{
	//read quad vtk
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++)
			getline(fin, stmp); //skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
			//cp[i].coor[0]/=6.8; cp[i].coor[1]/=6.8;
			cp[i].coor[0] /= 12.;
			cp[i].coor[1] /= 12.;
			cp[i].act = 1;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	InitialConnect_1();
	FindEdgeTopoDirec_1();
	FindKnotInterval_1();
	UpdateKnotInterval_1();
	SetLocalCoorSystem();
	FindIEN_1();
}

void TruncatedTspline::run_PatchTest_HighOrder(string fn)
{
	set_PatchTest_HighOrder(fn);
	//VisualizeTMesh("plate_hole_XP/CM_1");
	//VisualizeControlMesh("plate_hole_XP/CM_adapt_3");
	//cout<<"DOF: "<<cp.size()<<"\n";
	SetBezierMatIrrPatch();
	CollectActives();
	//VisualizeSurface("ptest_ho/surf_0");
}

///////////////////////////////////////////////////////////////////////////////////////////
// Xin Li's subdivision idea, hybrid non-uniform recursive subdivision

void TruncatedTspline::RunHybridCC(string fn_in, int nrf)
{
	//string fld("../io/input/");
	//string fn("patch_test_2ep");
	////string fn("patch_test_reg");
	//string fn_in(fld + fn);

	Preprocess(fn_in);
	//VisualizeEdge("../io/test5/square1");
	//VisualizeQuad("../io/test5/square1");
	//cout << "Done outputing CM!\n"; getchar();

	BuildTsplinesDual();
	//VisualizeEdge("../io/test5/square1");
	//cout << "Done building\n"; getchar();

	//int eid(104);
	//tmesh[eid].focus = 1;

	//int nrf(4);
	for (int i = 0; i < nrf; i++)
	{
		GlobalRefineDual();
	}

	//VisualizeEdge("../io/test1/square3");
	//VisualizeQuad("../io/test1/square3");
	//VisualizeGeomDual("../io/test1/square3");
}

void TruncatedTspline::Preprocess(string fn)
{
	ReadInput(fn);
	//Scale2Unit();
	InitialConnect_2();
	GetDual();
	BuildConnectDual();
}

void TruncatedTspline::ReadInput(string fn)
{
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++)
			getline(fin, stmp); //skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
			cp[i].act = 1;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline::Scale2Unit()
{
	double lmax[2] = {0., 0.};
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].coor[0] > lmax[0])
		{
			lmax[0] = cp[i].coor[0];
		}
		if (cp[i].coor[1] > lmax[1])
		{
			lmax[1] = cp[i].coor[1];
		}
	}
	for (uint i = 0; i < cp.size(); i++)
	{
		cp[i].coor[0] /= lmax[0];
		cp[i].coor[1] /= lmax[1];
	}
	//for (uint i = 0; i < cp.size(); i++)
	//{
	//	cp[i].coor[0] *= 4.;
	//	cp[i].coor[1] *= 4.;
	//}
}

void TruncatedTspline::GetDual()
{
	double coefv[3] = {4. / 9., 2. / 9., 1. / 9.};
	double coefp[2] = {2. / 3., 1. / 3.};
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].type == 4 /*|| tmesh[i].type == 5*/) //type==5 not working now
		{
			Vertex ptmp;
			ptmp.coor[0] = coefv[0] * cp[tmesh[i].cnct[0]].coor[0] + coefv[1] * cp[tmesh[i].cnct[1]].coor[0] +
						   coefv[2] * cp[tmesh[i].cnct[2]].coor[0] + coefv[1] * cp[tmesh[i].cnct[3]].coor[0];
			ptmp.coor[1] = coefv[0] * cp[tmesh[i].cnct[0]].coor[1] + coefv[1] * cp[tmesh[i].cnct[1]].coor[1] +
						   coefv[2] * cp[tmesh[i].cnct[2]].coor[1] + coefv[1] * cp[tmesh[i].cnct[3]].coor[1];
			ptmp.coor[2] = coefv[0] * cp[tmesh[i].cnct[0]].coor[2] + coefv[1] * cp[tmesh[i].cnct[1]].coor[2] +
						   coefv[2] * cp[tmesh[i].cnct[2]].coor[2] + coefv[1] * cp[tmesh[i].cnct[3]].coor[2];
			cp.push_back(ptmp);
			tmesh[i].add = cp.size() - 1;

			ptmp.coor[0] = coefp[0] * cp[tmesh[i].cnct[1]].coor[0] + coefp[1] * cp[tmesh[i].cnct[2]].coor[0];
			ptmp.coor[1] = coefp[0] * cp[tmesh[i].cnct[1]].coor[1] + coefp[1] * cp[tmesh[i].cnct[2]].coor[1];
			ptmp.coor[2] = coefp[0] * cp[tmesh[i].cnct[1]].coor[2] + coefp[1] * cp[tmesh[i].cnct[2]].coor[2];
			cp.push_back(ptmp);
			tmedge[tmesh[i].edge[1]].add = cp.size() - 1;
			if (tmedge[tmesh[i].edge[1]].pt[0] == tmesh[i].cnct[1])
				tmedge[tmesh[i].edge[1]].ploc = 0; //near edge.pt[0]
			else
				tmedge[tmesh[i].edge[1]].ploc = 1; //near edge.pt[1]

			ptmp.coor[0] = coefp[1] * cp[tmesh[i].cnct[2]].coor[0] + coefp[0] * cp[tmesh[i].cnct[3]].coor[0];
			ptmp.coor[1] = coefp[1] * cp[tmesh[i].cnct[2]].coor[1] + coefp[0] * cp[tmesh[i].cnct[3]].coor[1];
			ptmp.coor[2] = coefp[1] * cp[tmesh[i].cnct[2]].coor[2] + coefp[0] * cp[tmesh[i].cnct[3]].coor[2];
			cp.push_back(ptmp);
			tmedge[tmesh[i].edge[2]].add = cp.size() - 1;
			if (tmedge[tmesh[i].edge[2]].pt[0] == tmesh[i].cnct[3])
				tmedge[tmesh[i].edge[2]].ploc = 0;
			else
				tmedge[tmesh[i].edge[2]].ploc = 1;
		}
	}
	int flag(1);
	while (flag == 1)
	{
		flag = 0;
		vector<array<int, 3>> edual;
		for (uint i = 0; i < tmesh.size(); i++)
		{
			if (tmesh[i].type != 4 && tmesh[i].type != 5)
			{
				for (int j = 0; j < 4; j++)
				{
					if (tmedge[tmesh[i].edge[j]].add != -1 && tmedge[tmesh[i].edge[(j + 2) % 4]].add == -1)
					{
						array<int, 3> tmp = {i, (j + 2) % 4, (j + 3) % 4};
						int edid(tmesh[i].edge[j]);
						int ptloc(tmedge[tmesh[i].edge[j]].ploc);
						//if (tmesh[i].cnct[j]==tmedge[edid].pt[ptloc])
						//{
						//	tmp[2] = (j + 3) % 4;
						//}
						if (tmesh[i].cnct[(j + 1) % 4] == tmedge[edid].pt[ptloc])
						{
							tmp[2] = (j + 2) % 4;
							//tmp[3] = (j + 3) % 4;
						}
						edual.push_back(tmp);
					}
				}
			}
		}
		if (edual.size() != 0)
			flag = 1;
		for (uint i = 0; i < edual.size(); i++)
		{
			int edid(tmesh[edual[i][0]].edge[edual[i][1]]);
			if (tmedge[edid].add == -1)
			{
				if (tmedge[edid].pt[0] == tmesh[edual[i][0]].cnct[edual[i][2]])
					tmedge[edid].ploc = 0;
				else
					tmedge[edid].ploc = 1;
				int edpt[2] = {tmedge[edid].ploc, 1 - tmedge[edid].ploc};
				Vertex ptmp;
				ptmp.coor[0] = coefp[0] * cp[tmedge[edid].pt[edpt[0]]].coor[0] + coefp[1] * cp[tmedge[edid].pt[edpt[1]]].coor[0];
				ptmp.coor[1] = coefp[0] * cp[tmedge[edid].pt[edpt[0]]].coor[1] + coefp[1] * cp[tmedge[edid].pt[edpt[1]]].coor[1];
				ptmp.coor[2] = coefp[0] * cp[tmedge[edid].pt[edpt[0]]].coor[2] + coefp[1] * cp[tmedge[edid].pt[edpt[1]]].coor[2];
				cp.push_back(ptmp);
				tmedge[edid].add = cp.size() - 1;
			}
		}
	}
}

void TruncatedTspline::BuildConnectDual()
{
	int ned_old(tmedge.size());
	//find all dual edges
	for (uint i = 0; i < tmedge.size(); i++)
	{
		//if ((cp[tmedge[i].pt[0]].face.size() == 3 || cp[tmedge[i].pt[0]].face.size() > 4)
		//	|| (cp[tmedge[i].pt[1]].face.size() == 3 || cp[tmedge[i].pt[1]].face.size() > 4))
		if (cp[tmedge[i].pt[0]].type == 2 || cp[tmedge[i].pt[1]].type == 2)
		{
			tmedge[i].spoke = 1;
			tmedge[i].dual = 1;
		}
		else if (tmedge[i].face.size() == 2)
		{
			int fc[2] = {tmedge[i].face[0], tmedge[i].face[1]};
			int flag[2] = {0, 0};
			for (int j = 0; j < 2; j++)
			{
				int *it = find(tmesh[fc[j]].edge, tmesh[fc[j]].edge + 4, i);
				int loc(it - tmesh[fc[j]].edge);
				if (tmedge[tmesh[fc[j]].edge[(loc + 1) % 4]].add != -1 && tmedge[tmesh[fc[j]].edge[(loc + 3) % 4]].add != -1)
				{
					flag[j] = 1;
				}
			}
			if (flag[0] == 1 && flag[1] == 1) //non-spoke edge
			{
				tmedge[i].dual = 1;
			}
		}
	}
	//find elements influenced by dual edges
	for (uint i = 0; i < tmesh.size(); i++)
	{
		int flag(0);
		for (int j = 0; j < 4; j++)
		{
			if (tmedge[tmesh[i].edge[j]].dual == 1)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 1)
		{
			tmesh[i].dual = 1;
		}
	}
	//add connectivity for "dual" element (type==4), no edge info yet
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].type == 4) //element.dual must be 1
		{
			tmesh[i].cnctd[0] = tmesh[i].add;
			tmesh[i].cnctd[1] = tmedge[tmesh[i].edge[1]].add;
			tmesh[i].cnctd[2] = tmesh[i].cnct[2];
			tmesh[i].cnctd[3] = tmedge[tmesh[i].edge[2]].add;
			//if (tmesh[i].cnctd[0] == -1 || tmesh[i].cnctd[1] == -1 || tmesh[i].cnctd[2] == -1 || tmesh[i].cnctd[3] == -1)
			//{
			//	cout << "Irregular element (type-4) " << i << " has -1 cnctd\n";
			//	getchar();
			//}
		}
		else if (tmesh[i].type != 4 && tmesh[i].type != 5 && tmesh[i].dual == 1)
		{
			int loc(-1);
			for (int j = 0; j < 4; j++)
			{
				if (tmedge[tmesh[i].edge[j]].dual == 1)
				{
					loc = j;
					break;
				}
			}
			if (loc != -1)
			{
				tmesh[i].cnctd[0] = tmedge[tmesh[i].edge[(loc + 3) % 4]].add;
				tmesh[i].cnctd[1] = tmedge[tmesh[i].edge[(loc + 1) % 4]].add;
				tmesh[i].cnctd[2] = tmesh[i].cnct[(loc + 2) % 4];
				tmesh[i].cnctd[3] = tmesh[i].cnct[(loc + 3) % 4];
				//if (tmesh[i].cnctd[0] == -1 || tmesh[i].cnctd[1] == -1 || tmesh[i].cnctd[2] == -1 || tmesh[i].cnctd[3] == -1)
				//{
				//	cout << "Regular element (type-0 or type-1) " << i << " has -1 cnctd\n";
				//	cout << loc << "\n";
				//	cout << "cnct: " << tmesh[i].cnct[0] << " " << tmesh[i].cnct[1] << " " << tmesh[i].cnct[2] << " " << tmesh[i].cnct[3] << "\n";
				//	cout << "cnctd: " << tmesh[i].cnctd[0] << " " << tmesh[i].cnctd[1] << " " << tmesh[i].cnctd[2] << " " << tmesh[i].cnctd[3] << "\n";
				//	getchar();
				//}
			}
			else
			{
				cerr << "loc must be wrong!\n";
				getchar();
			}
		}
	}
	//add new elements, dual faces of EP and edges
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].type == 2) //EP
		{
			//sort the one-ring neighborhood counterclockwise
			vector<int> ring;
			ring.push_back(cp[i].face[0]);
			int enow(cp[i].face[0]), enext(-1), count(0);
			while (enext != cp[i].face[0])
			{
				enext = tmedge[tmesh[enow].edge[3]].face[0];
				if (enext == enow)
					enext = tmedge[tmesh[enow].edge[3]].face[1];
				if (enext != cp[i].face[0])
				{
					ring.push_back(enext);
				}
				else
				{
					break;
				}
				enow = enext;
				count++;
				if (count == cp[i].face.size())
					break;
			}
			Element etmp;
			etmp.act = 1;
			etmp.dual = 1;
			etmp.type = 6; //polygon element
			etmp.poly.resize(ring.size());
			etmp.polyed.resize(ring.size());
			for (uint j = 0; j < ring.size(); j++)
			{
				etmp.poly[j] = tmesh[ring[j]].add;
			}
			tmesh.push_back(etmp);
		}
	}
	for (uint i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].spoke == 1) //dual must be 1
		{
			//find the lower element
			int fcid(tmedge[i].face[0]), fc1(tmedge[i].face[1]);
			int *it = find(tmesh[fcid].edge, tmesh[fcid].edge + 4, i);
			int loc(it - tmesh[fcid].edge);
			if (loc != 3)
			{
				fcid = tmedge[i].face[1];
				fc1 = tmedge[i].face[0];
			}
			Element etmp;
			etmp.act = 1;
			etmp.dual = 1;
			etmp.type = 2;
			etmp.cnctd[0] = tmesh[fcid].add;
			etmp.cnctd[1] = tmedge[tmesh[fcid].edge[2]].add;
			etmp.cnctd[2] = tmedge[tmesh[fc1].edge[1]].add;
			etmp.cnctd[3] = tmesh[fc1].add;
			etmp.edlen[1] = 0.;
			etmp.edlen[3] = 0.;
			tmesh.push_back(etmp);
		}
		else if (tmedge[i].spoke == 0 && tmedge[i].dual == 1)
		{
			//find the lower element
			int fc0(tmedge[i].face[0]), fc1(tmedge[i].face[1]);
			int *it0 = find(tmesh[fc0].edge, tmesh[fc0].edge + 4, i);
			int loc0(it0 - tmesh[fc0].edge);
			int *it1 = find(tmesh[fc1].edge, tmesh[fc1].edge + 4, i);
			int loc1(it1 - tmesh[fc1].edge);
			Element etmp;
			etmp.act = 1;
			etmp.dual = 1;
			etmp.type = 2; //type 2 or 3
			if (tmedge[tmesh[fc0].edge[loc0]].len == 0.)
			{
				etmp.type = 3;
			}
			etmp.cnctd[0] = tmedge[tmesh[fc0].edge[(loc0 + 1) % 4]].add;
			etmp.cnctd[1] = tmedge[tmesh[fc0].edge[(loc0 + 3) % 4]].add;
			etmp.cnctd[2] = tmedge[tmesh[fc1].edge[(loc1 + 1) % 4]].add;
			etmp.cnctd[3] = tmedge[tmesh[fc1].edge[(loc1 + 3) % 4]].add;
			etmp.edlen[0] = tmedge[tmesh[fc0].edge[loc0]].len;
			etmp.edlen[1] = 0.;
			etmp.edlen[2] = tmedge[tmesh[fc0].edge[loc0]].len;
			etmp.edlen[3] = 0.;
			tmesh.push_back(etmp);
		}
	}
	//update connectivity, active and passive
	for (uint i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].dual == 1)
		{
			cp[tmedge[i].pt[0]].act = 0;
			cp[tmedge[i].pt[1]].act = 0;
			tmedge[i].act = 0;
		}
		if (tmedge[i].add != -1)
		{
			tmedge[i].act = 0;
		}
	}
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].dual == 1 && tmesh[i].type != 6)
		{
			for (int j = 0; j < 4; j++)
			{
				tmesh[i].cnct[j] = tmesh[i].cnctd[j];
			}
		}
	}
	//build edges for dual elements
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].dual == 1 && tmesh[i].type != 6)
		{
			for (int j = 0; j < 4; j++)
			{
				Edge edtmp;
				edtmp.pt[0] = tmesh[i].cnct[j];
				edtmp.pt[1] = tmesh[i].cnct[(j + 1) % 4];
				vector<Edge>::iterator it = find(tmedge.begin(), tmedge.end(), edtmp);
				tmesh[i].edge[j] = it - tmedge.begin();
				if (it == tmedge.end())
				{
					tmedge.push_back(edtmp);
				}
				if (tmesh[i].edge[j] >= ned_old)
				{
					tmedge[tmesh[i].edge[j]].len = tmesh[i].edlen[j];
				}
			}
		}
		else if (tmesh[i].dual == 1 && tmesh[i].type == 6)
		{
			for (uint j = 0; j < tmesh[i].poly.size(); j++)
			{
				Edge edtmp;
				edtmp.pt[0] = tmesh[i].poly[j];
				edtmp.pt[1] = tmesh[i].poly[(j + 1) % tmesh[i].poly.size()];
				vector<Edge>::iterator it = find(tmedge.begin(), tmedge.end(), edtmp);
				tmesh[i].polyed[j] = it - tmedge.begin();
				if (it == tmedge.end())
				{
					tmedge.push_back(edtmp);
				}
				tmedge[tmesh[i].polyed[j]].len = 0.;
			}
		}
	}
}

void TruncatedTspline::VisualizeEdge(string fn)
{
	string fname(fn + "_edge.vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		int count(0);
		vector<int> pid_act(cp.size(), -1);
		for (uint i = 0; i < cp.size(); i++)
		{
			if (cp[i].act == 1)
			{
				pid_act[i] = count++;
			}
		}
		fout << "POINTS " << count << " float\n";
		for (uint i = 0; i < cp.size(); i++)
		{
			if (cp[i].act == 1)
			{
				fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
			}
		}
		int ned(0);
		for (uint i = 0; i < tmedge.size(); i++)
		{
			if (tmedge[i].act == 1)
				ned++;
		}
		fout << "\nCELLS " << ned << " " << 3 * ned << '\n';
		for (uint i = 0; i < tmedge.size(); i++)
		{
			if (tmedge[i].act == 1)
				fout << "2 " << pid_act[tmedge[i].pt[0]] << " " << pid_act[tmedge[i].pt[1]] << '\n';
		}
		fout << "\nCELL_TYPES " << ned << '\n';
		for (uint i = 0; i < tmedge.size(); i++)
		{
			if (tmedge[i].act == 1)
				fout << "3\n";
		}
		fout << "\nCELL_DATA " << ned << "\nSCALARS len float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < tmedge.size(); i++)
		{
			if (tmedge[i].act == 1)
			{
				//fout << tmedge[i].len << "\n";
				//fout << tmedge[i].dual << "\n";
				fout << tmedge[i].face.size() << "\n";
				//if (i == 582)
				//{
				//	cout << tmedge[i].face[0] << "\n"; getchar();
				//}
			}
		}
		//fout << "POINT_DATA " << cp.size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<cp.size(); i++)
		//{
		//	fout << cp[i].trun << "\n";
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

void TruncatedTspline::VisualizeQuad(string fn)
{
	//vector<int> IDa(cp.size() + bzcp.size(), -1);
	////vector<int> bzcpaID(bzcp.size(), -1);
	//int count(0);
	//for (uint i = 0; i < cp.size(); i++)
	//{
	//	if (cp[i].act == 1) IDa[i] = count++;
	//}
	//int ncpa(count);
	//for (uint i = 0; i < bzcp.size(); i++)
	//{
	//	IDa[cp.size() + i] = count++;
	//}

	string fname(fn + "_quad.vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << cp.size() + bzcp.size() << " float\n";
		fout.precision(16);
		for (uint i = 0; i < cp.size(); i++)
		{
			fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
		}
		for (uint i = 0; i < bzcp.size(); i++)
		{
			fout << bzcp[i][0] << " " << bzcp[i][1] << " " << bzcp[i][2] << "\n";
		}
		int nquad(0);
		for (uint i = 0; i < tmesh.size(); i++)
		{
			if (tmesh[i].act == 1 && tmesh[i].type != 6)
				nquad++;
			if (tmesh[i].act == 1 && tmesh[i].c1 == 1)
				nquad++;
		}
		fout << "\nCELLS " << nquad << " " << 5 * nquad << '\n';
		for (uint i = 0; i < tmesh.size(); i++)
		{
			if (tmesh[i].act == 1 && tmesh[i].type != 6)
			{
				fout << "4 " << tmesh[i].cnct[0] << " " << tmesh[i].cnct[1] << " " << tmesh[i].cnct[2] << " " << tmesh[i].cnct[3] << "\n";
			}
		}
		for (uint i = 0; i < tmesh.size(); i++)
		{
			if (tmesh[i].act == 1 && tmesh[i].c1 == 1)
			{
				fout << "4 " << cp.size() + tmesh[i].IENc1[0] << " " << cp.size() + tmesh[i].IENc1[1] << " "
					 << cp.size() + tmesh[i].IENc1[2] << " " << cp.size() + tmesh[i].IENc1[3] << "\n";
			}
		}
		fout << "\nCELL_TYPES " << nquad << '\n';
		for (int i = 0; i < nquad; i++)
		{
			fout << "9\n";
		}
		//fout << "\nCELL_DATA " << nquad << "\nSCALARS type float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<tmesh.size(); i++)
		//{
		//	if (tmesh[i].act == 1 && tmesh[i].type != 6)
		//	{
		//		//fout << tmesh[i].dual << "\n";
		//		fout << tmesh[i].c1 << "\n";
		//	}
		//}
		//for (uint i = 0; i<tmesh.size(); i++)
		//{
		//	if (tmesh[i].act == 1 && tmesh[i].c1 == 1)
		//	{
		//		fout << "0\n";
		//	}
		//}
		//fout << "POINT_DATA " << IDa.size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<IDa.size(); i++)
		//{
		//	fout << IDa[i] << "\n";
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

void TruncatedTspline::BuildTsplinesDual()
{
	UpdateConectDual();

	//these functions are more for T-splines
	//FindEdgeTopoDirecDual();
	//FindKnotIntervalDual();//find knottmp
	//UpdateKnotIntervalDual();//assign knottmp to knot
	//SetLocalCoorSystemDual();

	FindIENDual();
	UpdateIENDual();

	//int eid(115);
	//cout << "IEN size: " << tmesh[eid].IEN.size() << "\n";
	//for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	//{
	//	cout << tmesh[eid].IEN[i] << " ";
	//}
	//cout << "\n";
	//for (uint i = 0; i < tmesh[eid].patch_ku.size(); i++)
	//{
	//	cout << i << ": ";
	//	for (int j = 0; j < 5; j++)
	//	{
	//		cout << tmesh[eid].patch_kv[i][j] << " ";
	//	}
	//	cout << "\n";
	//}
	//getchar();

	//int eid(115);
	//cout << tmesh[eid].smat[0] << "\n\n";
	//cout << tmesh[eid].smat[1] << "\n\n";
	//cout << tmesh[eid].smat[2] << "\n\n";
	//cout << tmesh[eid].smat[3] << "\n\n";
	//getchar();

	//int eid(115);
	//double uv[2] = { 0.26,0.26 };
	//vector<double> Nt;
	//vector<array<double, 2>> dNdt;
	//ElementBasis_IrregularDual(eid, uv[0], uv[1], Nt, dNdt);
	//double sum(0.);
	//for (uint i = 0; i < Nt.size(); i++)
	//{
	//	sum += Nt[i];
	//}
	//cout << sum << "\n"; getchar();
}

void TruncatedTspline::UpdateConectDual()
{
	uint i, j, k;
	for (i = 0; i < cp.size(); i++)
	{
		cp[i].face.clear();
		cp[i].edge.clear();
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		tmedge[i].face.clear();
	}
	for (i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].type != 6)
		{
			for (j = 0; j < 4; j++)
			{
				if (cp[tmesh[i].cnct[j]].act == 1)
				{
					cp[tmesh[i].cnct[j]].face.push_back(i);
				}
				if (tmedge[tmesh[i].edge[j]].act == 1)
				{
					tmedge[tmesh[i].edge[j]].face.push_back(i);
				}
			}
		}
		else if (tmesh[i].act == 1 && tmesh[i].type == 6)
		{
			for (j = 0; j < tmesh[i].poly.size(); j++)
			{
				if (cp[tmesh[i].poly[j]].act == 1)
				{
					cp[tmesh[i].poly[j]].face.push_back(i);
				}
			}
			for (j = 0; j < tmesh[i].polyed.size(); j++)
			{
				if (tmedge[tmesh[i].polyed[j]].act == 1)
				{
					tmedge[tmesh[i].polyed[j]].face.push_back(i);
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

void TruncatedTspline::FindEdgeTopoDirecDual()
{
	//T-junctions not considered
	for (uint i = 0; i < tmedge.size(); i++)
	{
		tmedge[i].pn[0][0] = 3;
		tmedge[i].pn[0][1] = -1;
		tmedge[i].pn[1][0] = 3;
		tmedge[i].pn[1][1] = -1; //initialize as end
		if (tmedge[i].act == 1)
		{
			for (int j = 0; j < 2; j++)
			{
				if (cp[tmedge[i].pt[j]].type == 0) //regular (only two types: regular 0 and irregular 2)
				{
					for (uint k = 0; k < cp[tmedge[i].pt[j]].edge.size(); k++)
					{
						if (cp[tmedge[i].pt[j]].edge[k] != i)
						{
							int flag(0);
							for (uint i1 = 0; i1 < tmedge[i].face.size(); i1++)
							{
								if (tmesh[tmedge[i].face[i1]].type != 6)
								{
									for (uint j1 = 0; j1 < 4; j1++)
									{
										int ed(tmesh[tmedge[i].face[i1]].edge[j1]);
										if (tmedge[ed].act == 1 && ed == cp[tmedge[i].pt[j]].edge[k])
										{
											flag = 1;
											break;
										}
									}
								}
								else
								{
									for (uint j1 = 0; j1 < tmesh[tmedge[i].face[i1]].polyed.size(); j1++)
									{
										int ed(tmesh[tmedge[i].face[i1]].polyed[j1]);
										if (tmedge[ed].act == 1 && ed == cp[tmedge[i].pt[j]].edge[k])
										{
											flag = 1;
											break;
										}
									}
								}
							}
							if (flag == 0)
							{
								tmedge[i].pn[j][0] = 0;
								tmedge[i].pn[j][1] = cp[tmedge[i].pt[j]].edge[k];
								break;
							}
						}
					}
				}
				//else if (cp[tmedge[i].pt[j]].type == 1)//T-junctions
				//{
				//	int fid(-1);
				//	int loc(0);
				//	for (uint k = 0; k<cp[tmedge[i].pt[j]].face.size(); k++)
				//	{
				//		int ftmp(cp[tmedge[i].pt[j]].face[k]);
				//		for (int k1 = 0; k1<4; k1++)
				//		{
				//			if (tmedge[tmesh[ftmp].edge[k1]].act == 0 && tmedge[tmesh[ftmp].edge[k1]].midpt == tmedge[i].pt[j])
				//			{
				//				loc = k1; fid = ftmp; break;
				//			}
				//		}
				//		if (fid != -1)
				//		{
				//			break;
				//		}
				//	}
				//	if (fid == -1)
				//	{
				//		cout << "edge id: " << i << "\n";
				//		cout << tmedge[i].pt[0] << " " << tmedge[i].pt[1] << "\n";
				//		cout << tmedge[i].prt << "\n";
				//		//cout<<cp[tmedge[i].pt[0]].face.size()<<" "<<cp[tmedge[i].pt[1]].face.size()<<"\n";
				//		cout << cp[tmedge[i].pt[0]].face[0] << " " << cp[tmedge[i].pt[0]].face[1] << " " << cp[tmedge[i].pt[0]].face[2] << "\n";
				//		cerr << "T-junction cannot be found in any neighboring elements!\n";
				//		getchar();
				//	}
				//	if (tmedge[tmesh[fid].edge[loc]].chd[0] == i)
				//	{
				//		tmedge[i].pn[j][0] = 0;
				//		tmedge[i].pn[j][1] = tmedge[tmesh[fid].edge[loc]].chd[1];
				//	}
				//	else if (tmedge[tmesh[fid].edge[loc]].chd[1] == i)
				//	{
				//		tmedge[i].pn[j][0] = 0;
				//		tmedge[i].pn[j][1] = tmedge[tmesh[fid].edge[loc]].chd[0];
				//	}
				//	else
				//	{
				//		tmedge[i].pn[j][0] = 1;
				//		tmedge[i].pn[j][1] = fid;
				//	}
				//}
				else if (cp[tmedge[i].pt[j]].type == 2) //extraordinary
				{
					//tmedge[i].pn[j][0] = 2;
					cerr << "Not valid: edge next is an EP!\n";
					getchar();
				}
			}
		}
	}
}

void TruncatedTspline::FindKnotIntervalDual()
{
	for (uint i = 0; i < cp.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[i].kitvUtmp[j] = 1.;
			cp[i].kitvVtmp[j] = 1.;
		}
		if (cp[i].act == 1 && cp[i].type != 2) //active regular interior or boundary vertex
		{
			int pos(0);
			cp[i].rfc = -1;
			for (uint j = 0; j < cp[i].face.size(); j++)
			{
				if (tmesh[cp[i].face[j]].type == 0 || tmesh[cp[i].face[j]].type == 4) //type==4 later needs further treatment
				{
					int *it = find(tmesh[cp[i].face[j]].cnct, tmesh[cp[i].face[j]].cnct + 4, i);
					if (it != tmesh[cp[i].face[j]].cnct + 4)
					{
						cp[i].rfc = cp[i].face[j];
						pos = it - tmesh[cp[i].face[j]].cnct;
						break;
					}
				}
			}
			if (cp[i].rfc == -1)
			{
				cerr << "Cannot find correct reference face!\n";
				getchar();
			}
			cp[i].uved[0] = tmesh[cp[i].rfc].edge[pos];
			cp[i].uved[1] = tmesh[cp[i].rfc].edge[(pos + 3) % 4];
			/*if (tmedge[cp[i].uved[0]].act == 0)
			{
				int edtmp(tmedge[cp[i].uved[0]].chd[0]);
				if (tmedge[edtmp].pt[0] == i || tmedge[edtmp].pt[1] == i)
				{
					cp[i].uved[0] = tmedge[cp[i].uved[0]].chd[0];
				}
				else
				{
					cp[i].uved[0] = tmedge[cp[i].uved[0]].chd[1];
				}
			}
			if (tmedge[cp[i].uved[1]].act == 0)
			{
				int edtmp(tmedge[cp[i].uved[1]].chd[0]);
				if (tmedge[edtmp].pt[0] == i || tmedge[edtmp].pt[1] == i)
				{
					cp[i].uved[1] = tmedge[cp[i].uved[1]].chd[0];
				}
				else
				{
					cp[i].uved[1] = tmedge[cp[i].uved[1]].chd[1];
				}
			}*/
			vector<int>::iterator it1 = find(cp[i].edge.begin(), cp[i].edge.end(), cp[i].uved[0]);
			vector<int>::iterator it2 = find(cp[i].edge.begin(), cp[i].edge.end(), cp[i].uved[1]);
			if (it1 == cp[i].edge.end() || it2 == cp[i].edge.end())
			{
				cerr << "Cannot find correct uv edges!\n";
				getchar();
			}
			ShootRay(i, cp[i].uved[0], cp[i].kitvUtmp);
			ShootRay(i, cp[i].uved[1], cp[i].kitvVtmp);
		}
		//else
		//{
		//	for (int j = 0; j<4; j++)
		//	{
		//		cp[i].kitvUtmp[j] = tmedge[cp[i].edge[0]].len; cp[i].kitvVtmp[j] = tmedge[cp[i].edge[0]].len;
		//	}
		//}
	}
}

void TruncatedTspline::ShootRayDual(int pid, int edid, double kv[4])
{
	int loc0(0), loc1(1); //loc0 is the first endpoint, loc1 is the second endpoint
	//positive direction
	if (pid != tmedge[edid].pt[0])
	{
		loc0 = 1;
		loc1 = 0;
	}
	kv[2] = tmedge[edid].len;
	if (tmedge[edid].pn[loc1][0] == 0) //next is edge
	{
		kv[3] = tmedge[tmedge[edid].pn[loc1][1]].len;
	}
	//else if (tmedge[edid].pn[loc1][0] == 1)//next is face
	//{
	//	int fid(tmedge[edid].pn[loc1][1]), pos(0);
	//	for (int i = 0; i<4; i++)
	//	{
	//		if (tmedge[tmesh[fid].edge[i]].act == 0 && tmedge[tmesh[fid].edge[i]].midpt == tmedge[edid].pt[loc1])
	//		{
	//			pos = i; break;
	//		}
	//	}
	//	kv[3] = tmedge[tmesh[fid].edge[(pos + 1) % 4]].len;
	//}
	//else if (tmedge[edid].pn[loc1][0] == 2)//next is XP
	//{
	//	kv[3] = kv[2];
	//}
	else if (tmedge[edid].pn[loc1][0] == 3) //end
	{
		kv[3] = 0.;
	}
	//negative direction
	if (tmedge[edid].pn[loc0][0] == 0) //previous is edge
	{
		int ed0 = tmedge[edid].pn[loc0][1];
		kv[1] = tmedge[ed0].len;
		int a0 = 0, a1 = 1;
		if (tmedge[ed0].pt[0] != pid)
		{
			a0 = 1;
			a1 = 0;
		}
		if (tmedge[ed0].pn[a1][0] == 0) //previous previous is edge
		{
			kv[0] = tmedge[tmedge[ed0].pn[a1][1]].len;
		}
		//else if (tmedge[ed0].pn[a1][0] == 1)
		//{
		//	int pt0(tmedge[ed0].pt[a1]), fid(tmedge[ed0].pn[a1][1]), pos(0);
		//	for (int i = 0; i<4; i++)
		//	{
		//		if (tmedge[tmesh[fid].edge[i]].act == 0 && tmedge[tmesh[fid].edge[i]].midpt == pt0)
		//		{
		//			pos = i; break;
		//		}
		//	}
		//	kv[0] = tmedge[tmesh[fid].edge[(pos + 1) % 4]].len;
		//}
		//else if (tmedge[ed0].pn[a1][0] == 2)
		//{
		//	kv[0] = kv[1];
		//}
		else if (tmedge[ed0].pn[a1][0] == 3) //previous previous is end
		{
			kv[0] = 0.;
		}
	}
	/*else if (tmedge[edid].pn[loc0][0] == 1)
	{
		int fid0(tmedge[edid].pn[loc0][1]), pos(0);
		for (int i = 0; i<4; i++)
		{
			if (tmedge[tmesh[fid0].edge[i]].act == 0 && tmedge[tmesh[fid0].edge[i]].midpt == pid)
			{
				pos = i; break;
			}
		}
		kv[1] = tmedge[tmesh[fid0].edge[(pos + 1) % 4]].len;
		int ed0(tmesh[fid0].edge[(pos + 2) % 4]);
		if (tmedge[ed0].act == 1)
		{
			if (tmedge[ed0].face.size() == 2)
			{
				int fid1(tmedge[ed0].face[0]);
				if (fid1 == fid0) fid1 = tmedge[ed0].face[1];
				int* it = find(tmesh[fid1].edge, tmesh[fid1].edge + 4, ed0);
				int pos1(it - tmesh[fid1].edge);
				kv[0] = tmedge[tmesh[fid1].edge[(pos1 + 1) % 4]].len;
			}
			else
			{
				kv[0] = 0.;
			}
		}
		else
		{
			int pt0(tmedge[ed0].midpt), ed1;
			for (uint i = 0; i<cp[pt0].edge.size(); i++)
			{
				if (cp[pt0].edge[i] != tmedge[ed0].chd[0] && cp[pt0].edge[i] != tmedge[ed0].chd[1])
				{
					ed1 = cp[pt0].edge[i]; break;
				}
			}
			kv[0] = tmedge[ed1].len;
		}
	}
	else if (tmedge[edid].pn[loc0][0] == 2)
	{
		kv[1] = kv[2]; kv[0] = kv[2];
	}*/
	else if (tmedge[edid].pn[loc0][0] == 3) //previous is end
	{
		kv[1] = 0.;
		kv[0] = 0.;
	}
}

void TruncatedTspline::UpdateKnotIntervalDual()
{
	for (uint i = 0; i < cp.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cp[i].kitvU[j] = cp[i].kitvUtmp[j];
			cp[i].kitvV[j] = cp[i].kitvVtmp[j];
		}
	}
}

void TruncatedTspline::SetLocalCoorSystemDual()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		tmesh[eid].node.clear();
		tmesh[eid].lcs.clear();
		if (tmesh[eid].act == 1 && tmesh[eid].type != 6) //active quad element, no matter regular or irregular
		{
			double ul[2] = {tmedge[tmesh[eid].edge[0]].len, tmedge[tmesh[eid].edge[1]].len};
			double uvcoor[8][2] = {{0., 0.}, {ul[0] / 2., 0.}, {ul[0], 0.}, {ul[0], ul[1] / 2.}, {ul[0], ul[1]}, {ul[0] / 2., ul[1]}, {0., ul[1]}, {0., ul[1] / 2.}};
			for (int i = 0; i < 4; i++)
			{
				int uved[2] = {tmesh[eid].edge[i], tmesh[eid].edge[(i + 3) % 4]};
				//for (int j = 0; j<2; j++)
				//{
				//	if (tmedge[uved[j]].act == 0)
				//	{
				//		int edtmp(tmedge[uved[j]].chd[0]);
				//		if (tmedge[edtmp].pt[0] != tmesh[eid].cnct[i] && tmedge[edtmp].pt[1] != tmesh[eid].cnct[i])
				//		{
				//			edtmp = tmedge[uved[j]].chd[1];
				//		}
				//		uved[j] = edtmp;
				//	}
				//}
				tmesh[eid].node.push_back(tmesh[eid].cnct[i]);
				ELCS tmp;
				tmp.u[0] = uvcoor[2 * i][0];
				tmp.u[1] = uvcoor[2 * i][1];
				if (cp[tmesh[eid].cnct[i]].uved[0] == uved[0])
				{
					tmp.rot = i;
				}
				else if (cp[tmesh[eid].cnct[i]].uved[0] == uved[1])
				{
					tmp.rot = (i + 1) % 4;
				}
				else if (cp[tmesh[eid].cnct[i]].uved[1] == uved[0])
				{
					tmp.rot = (i + 3) % 4;
				}
				else
				{
					tmp.rot = (i + 2) % 4;
				}
				tmesh[eid].lcs.push_back(tmp);

				//if (tmedge[tmesh[eid].edge[i]].act == 0)
				//{
				//	int ued(tmedge[tmesh[eid].edge[i]].chd[1]);
				//	if (tmedge[ued].pt[0] == tmesh[eid].cnct[i] || tmedge[ued].pt[1] == tmesh[eid].cnct[i])
				//	{
				//		ued = tmedge[tmesh[eid].edge[i]].chd[0];
				//	}
				//	tmesh[eid].node.push_back(tmedge[tmesh[eid].edge[i]].midpt);
				//	ELCS tmp1;
				//	tmp1.u[0] = uvcoor[2 * i + 1][0]; tmp1.u[1] = uvcoor[2 * i + 1][1];
				//	if (cp[tmedge[tmesh[eid].edge[i]].midpt].uved[1] == ued)
				//	{
				//		tmp1.rot = (i + 3) % 4;
				//	}
				//	else
				//	{
				//		tmp1.rot = (i + 2) % 4;
				//	}
				//	tmesh[eid].lcs.push_back(tmp1);
				//}
			}
		}
	}
}

void TruncatedTspline::FindIENDual()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		tmesh[eid].IENtmp.clear();
		tmesh[eid].patch_kutmp.clear();
		tmesh[eid].patch_kvtmp.clear();
		if (tmesh[eid].act == 1 && tmesh[eid].type == 0)
		{
			FindElementIEN_RegDual(eid);
		}
		else if (tmesh[eid].act == 1 && tmesh[eid].type == 4)
		{
			FindElementIEN_IrrDual(eid);
		}
	}
}

void TruncatedTspline::FindElementIEN_RegDual(int eid)
{
	int vtloc[4][2] = {{5, 0}, {6, 3}, {10, 15}, {9, 12}};
	int edloc[4][2] = {{1, 2}, {7, 11}, {14, 13}, {8, 4}};
	int nbf(16);
	tmesh[eid].IENtmp.resize(nbf);
	tmesh[eid].patch_kutmp.resize(nbf);
	tmesh[eid].patch_kvtmp.resize(nbf);
	//double kint_u[7] = { 0.,0.,0.,tmedge[tmesh[eid].edge[0]].len,0.,0.,0. };
	//double kint_v[7] = { 0.,0.,0.,tmedge[tmesh[eid].edge[3]].len,0.,0.,0. };
	//IEN
	int ednb[4] = {-1, -1, -1, -1};
	for (int i = 0; i < 4; i++)
	{
		ednb[i] = tmedge[tmesh[eid].edge[i]].face[0];
		if (ednb[i] == eid)
			ednb[i] = tmedge[tmesh[eid].edge[i]].face[1];
		int *it = find(tmesh[ednb[i]].cnct, tmesh[ednb[i]].cnct + 4, tmesh[eid].cnct[i]);
		int loc(it - tmesh[ednb[i]].cnct);
		tmesh[eid].IENtmp[edloc[i][0]] = tmesh[ednb[i]].cnct[(loc + 1) % 4];
		tmesh[eid].IENtmp[edloc[i][1]] = tmesh[ednb[i]].cnct[(loc + 2) % 4];
	}
	for (int i = 0; i < 4; i++)
	{
		tmesh[eid].IENtmp[vtloc[i][0]] = tmesh[eid].cnct[i];
		int vtnb(-1);
		for (uint j = 0; j < cp[tmesh[eid].cnct[i]].face.size(); j++)
		{
			int fcid(cp[tmesh[eid].cnct[i]].face[j]);
			if (fcid != eid && fcid != ednb[i] && fcid != ednb[(i + 3) % 4])
			{
				vtnb = fcid;
				break;
			}
		}
		int *it = find(tmesh[vtnb].cnct, tmesh[vtnb].cnct + 4, tmesh[eid].cnct[i]);
		int loc(it - tmesh[vtnb].cnct);
		tmesh[eid].IENtmp[vtloc[i][1]] = tmesh[vtnb].cnct[(loc + 2) % 4];
	}
	//patch kv
	double ku[8], kv[8];
	FindPatchKVDual(eid, ku, kv);
	int loc(0);
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 5; k++)
			{
				tmesh[eid].patch_kutmp[loc][k] = ku[j + k];
				tmesh[eid].patch_kvtmp[loc][k] = kv[i + k];
			}
			loc++;
		}
	}
	tmesh[eid].pku.clear();
	tmesh[eid].pkv.clear();
	tmesh[eid].pku.resize(8);
	tmesh[eid].pkv.resize(8);
	tmesh[eid].pku.assign(ku, ku + 8);
	tmesh[eid].pkv.assign(kv, kv + 8);
}

void TruncatedTspline::FindElementIEN_IrrDual(int eid)
{
	int nv(-1), polyid(-1), ploc;
	for (uint i = 0; i < cp[tmesh[eid].cnct[0]].face.size(); i++)
	{
		if (tmesh[cp[tmesh[eid].cnct[0]].face[i]].type == 6)
		{
			polyid = cp[tmesh[eid].cnct[0]].face[i];
			nv = tmesh[polyid].poly.size();
			vector<int>::iterator it = find(tmesh[polyid].poly.begin(), tmesh[polyid].poly.end(), tmesh[eid].cnct[0]);
			ploc = it - tmesh[polyid].poly.begin();
		}
	}
	if (nv == -1)
	{
		cerr << "EP has regular valence!\n";
		return;
	}
	int nbf(nv + 12);
	tmesh[eid].IENtmp.resize(nbf);
	tmesh[eid].patch_kutmp.resize(nbf); //only the last 7 will be used, and they are associated with B-splines
	tmesh[eid].patch_kvtmp.resize(nbf);
	//IEN
	for (uint i = 0; i < tmesh[polyid].poly.size(); i++)
	{
		tmesh[eid].IENtmp[i] = tmesh[polyid].poly[(ploc + i) % nv];
	}
	for (int i = 1; i < 4; i++)
	{
		tmesh[eid].IENtmp[nv + i] = tmesh[eid].cnct[i];
	}
	int ednb[4] = {-1, -1, -1, -1};
	for (int i = 0; i < 4; i++)
	{
		ednb[i] = tmedge[tmesh[eid].edge[i]].face[0];
		if (ednb[i] == eid)
			ednb[i] = tmedge[tmesh[eid].edge[i]].face[1];
	}
	int iloc[3][3] = {{nv, nv + 5, nv + 6}, {nv + 7, nv + 8, nv + 9}, {nv + 10, nv + 11, nv + 4}};
	for (int i = 1; i < 4; i++)
	{
		for (uint j = 0; j < cp[tmesh[eid].cnct[i]].face.size(); j++)
		{
			int fcid(cp[tmesh[eid].cnct[i]].face[j]);
			if (fcid != eid && fcid != ednb[i] && fcid != ednb[(i + 3) % 4])
			{
				int *it = find(tmesh[fcid].cnct, tmesh[fcid].cnct + 4, tmesh[eid].cnct[i]);
				int loc(it - tmesh[fcid].cnct);
				tmesh[eid].IENtmp[iloc[i - 1][0]] = tmesh[fcid].cnct[(loc + 1) % 4];
				tmesh[eid].IENtmp[iloc[i - 1][1]] = tmesh[fcid].cnct[(loc + 2) % 4];
				tmesh[eid].IENtmp[iloc[i - 1][2]] = tmesh[fcid].cnct[(loc + 3) % 4];
			}
		}
	}
	//patch kv
	double ku[8], kv[8];
	FindPatchKVDual(eid, ku, kv);
	int ist[12][2] = {{2, 0}, {2, 1}, {2, 2}, {1, 2}, {0, 2}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {2, 3}, {1, 3}, {0, 3}};
	for (uint i = 0; i < tmesh[eid].patch_kutmp.size(); i++)
	{
		if (i < nv)
		{
			for (int j = 0; j < 5; j++)
			{
				tmesh[eid].patch_kutmp[i][j] = 0.;
				tmesh[eid].patch_kvtmp[i][j] = 0.;
			}
		}
		else
		{
			for (int j = 0; j < 5; j++)
			{
				tmesh[eid].patch_kutmp[i][j] = ku[ist[i - nv][0] + j];
				tmesh[eid].patch_kvtmp[i][j] = kv[ist[i - nv][1] + j];
			}
		}
	}
	//GetSubMatrixDual(nv, tmesh[eid].smat);
	GetSubMatrixDual_All(nv, tmesh[eid].amat);
	GetSubMatrixDual_Bar(nv, tmesh[eid].abar);
	//GetSubPatchKVDual(tmesh[eid].subku, tmesh[eid].subkv);

	//cout << "eid: " << eid << "\n";
	//double sum(0.);
	////for (uint i = 0; i < tmesh[eid].smat[0].cols(); i++)
	//{
	//	double sum(0.);
	//	for (uint j = 0; j < tmesh[eid].smat[0].rows(); j++)
	//	{
	//		//sum += tmesh[eid].smat[0](j, i);
	//		cout << tmesh[eid].smat[0](j, 6) << " ";
	//	}
	//	//cout << "col: " << i << " " << sum << "\n";
	//}
	//getchar();

	//cout << "eid: " << eid << "\n";
	//cout << tmesh[eid].amat.cols() << "\n";
	//for (uint i = 0; i < tmesh[eid].amat.cols(); i++)
	//{
	//	double sum(0.);
	//	for (uint j = 0; j < tmesh[eid].amat.rows(); j++)
	//	{
	//		sum += tmesh[eid].amat(j, i);
	//		//cout << tmesh[eid].smat[0](j, 6) << " ";
	//		//if (i == 0)
	//		//{
	//		//	cout << tmesh[eid].amat(j, i) << " ";
	//		//}
	//	}
	//	cout << "col: " << i << " " << sum << "\n";
	//}
	//getchar();

	//cout << "eid: " << eid << "\n";
	//cout << tmesh[eid].abar.cols() << "\n";
	//for (uint i = 0; i < tmesh[eid].abar.cols(); i++)
	//{
	//	double sum(0.);
	//	for (uint j = 0; j < tmesh[eid].abar.rows(); j++)
	//	{
	//		sum += tmesh[eid].abar(j, i);
	//		//cout << tmesh[eid].smat[0](j, 6) << " ";
	//		//if (i == 0)
	//		//{
	//		//	cout << tmesh[eid].amat(j, i) << " ";
	//		//}
	//	}
	//	cout << "col: " << i << " " << sum << "\n";
	//}
	//getchar();
}

void TruncatedTspline::FindPatchKVDual(int eid, double ku[8], double kv[8])
{
	double kint_u[7] = {0., 0., 0., tmedge[tmesh[eid].edge[0]].len, 0., 0., 0.};
	double kint_v[7] = {0., 0., 0., tmedge[tmesh[eid].edge[3]].len, 0., 0., 0.};
	//knot intervals, positive direction
	int count(0), enow[2] = {eid, eid}, iu(1), iv(2);
	while (count < 2)
	{
		if (tmedge[tmesh[enow[0]].edge[iu]].face.size() == 2)
		{
			int enext = tmedge[tmesh[enow[0]].edge[iu]].face[0];
			if (enext == enow[0])
				enext = tmedge[tmesh[enow[0]].edge[iu]].face[1];
			int *it = find(tmesh[enext].edge, tmesh[enext].edge + 4, tmesh[enow[0]].edge[iu]);
			int loc(it - tmesh[enext].edge);
			iu = (loc + 2) % 4;
			enow[0] = enext;
		}
		if (tmedge[tmesh[enow[1]].edge[iv]].face.size() == 2)
		{
			int enext = tmedge[tmesh[enow[1]].edge[iv]].face[0];
			if (enext == enow[1])
				enext = tmedge[tmesh[enow[1]].edge[iv]].face[1];
			int *it = find(tmesh[enext].edge, tmesh[enext].edge + 4, tmesh[enow[1]].edge[iv]);
			int loc(it - tmesh[enext].edge);
			iv = (loc + 2) % 4;
			enow[1] = enext;
		}
		kint_u[4 + count] = tmedge[tmesh[enow[0]].edge[(iu + 1) % 4]].len;
		kint_v[4 + count] = tmedge[tmesh[enow[1]].edge[(iv + 1) % 4]].len;
		count++;
	}
	//knot intervals, negative direction
	count = 0;
	enow[0] = eid;
	enow[1] = eid;
	iu = 3;
	iv = 0;
	while (count < 2)
	{
		if (tmedge[tmesh[enow[0]].edge[iu]].face.size() == 2)
		{
			int enext = tmedge[tmesh[enow[0]].edge[iu]].face[0];
			if (enext == enow[0])
				enext = tmedge[tmesh[enow[0]].edge[iu]].face[1];
			int *it = find(tmesh[enext].edge, tmesh[enext].edge + 4, tmesh[enow[0]].edge[iu]);
			int loc(it - tmesh[enext].edge);
			iu = (loc + 2) % 4;
			enow[0] = enext;
		}
		if (tmedge[tmesh[enow[1]].edge[iv]].face.size() == 2)
		{
			int enext = tmedge[tmesh[enow[1]].edge[iv]].face[0];
			if (enext == enow[1])
				enext = tmedge[tmesh[enow[1]].edge[iv]].face[1];
			int *it = find(tmesh[enext].edge, tmesh[enext].edge + 4, tmesh[enow[1]].edge[iv]);
			int loc(it - tmesh[enext].edge);
			iv = (loc + 2) % 4;
			enow[1] = enext;
		}
		kint_u[2 - count] = tmedge[tmesh[enow[0]].edge[(iu + 1) % 4]].len;
		kint_v[2 - count] = tmedge[tmesh[enow[1]].edge[(iv + 1) % 4]].len;
		count++;
	}
	ku[3] = 0.;
	ku[4] = kint_u[3];
	kv[3] = 0.;
	kv[4] = kint_v[3];
	//double ku[8] = { 0.,0.,0.,0.,kint_u[3],0.,0.,0. };
	//double kv[8] = { 0.,0.,0.,0.,kint_v[3],0.,0.,0. };
	for (int i = 0; i < 3; i++)
	{
		ku[5 + i] = ku[4 + i] + kint_u[4 + i];
		ku[2 - i] = ku[3 - i] - kint_u[2 - i];
		kv[5 + i] = kv[4 + i] + kint_v[4 + i];
		kv[2 - i] = kv[3 - i] - kint_v[2 - i];
	}
}

void TruncatedTspline::UpdateIENDual()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		vector<int>().swap(tmesh[eid].IEN);
		vector<array<double, 5>>().swap(tmesh[eid].patch_ku);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kv);
		tmesh[eid].IEN = tmesh[eid].IENtmp;
		tmesh[eid].patch_ku = tmesh[eid].patch_kutmp;
		tmesh[eid].patch_kv = tmesh[eid].patch_kvtmp;
		vector<int>().swap(tmesh[eid].IENtmp);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kutmp);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kvtmp);

		//if (eid == 6)
		//{
		//	cout << "IEN: ";
		//	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
		//	{
		//		cout << tmesh[eid].IEN[i] << " ";
		//	}
		//	cout << "\n\n";
		//	cout << "patch ku:\n";
		//	for (uint i = 0; i < tmesh[eid].patch_ku.size(); i++)
		//	{
		//		cout << "id " << i << ": ";
		//		for (int j = 0; j < 5; j++)
		//		{
		//			cout << tmesh[eid].patch_ku[i][j] << " ";
		//		}
		//		cout << "\n";
		//	}
		//	getchar();
		//}
	}
}

void TruncatedTspline::VisualizeGeomDual(string fn)
{
	vector<array<double, 3>> spt;
	vector<array<double, 3>> sval;
	vector<array<int, 4>> sele;
	vector<array<double, 3>> lpt; //visulize parameter lines
	vector<array<int, 2>> led;	  //line connectivity
	int ns(2), ecount(0), loc0, loc1, loc2;
	vector<double> su(ns), sv(ns);
	for (int i = 0; i < ns; i++)
	{
		su[i] = i * 1. / (ns - 1);
		sv[i] = i * 1. / (ns - 1);
	}

	for (uint e = 0; e < tmesh.size(); e++)
	{
		if (tmesh[e].act == 1 && (tmesh[e].type == 0 || tmesh[e].type == 1 || tmesh[e].type == 4))
		{
			int loc(0);
			//vector<vector<double>> bmat;
			//vector<double> su(ns), sv(ns);
			//for (int i = 0; i<ns; i++)
			//{
			//	su[i] = i*tmedge[tmesh[e].edge[0]].len / (ns - 1);
			//	sv[i] = i*tmedge[tmesh[e].edge[3]].len / (ns - 1);
			//}
			for (int a = 0; a < ns; a++)
			{
				for (int b = 0; b < ns; b++)
				{
					array<double, 3> pt;
					array<double, 3> nm;
					//SurfacePointCal1(su[b],sv[a],tmesh[e].IEN,bmat,pt,nm);
					SurfacePointMap(e, su[b], sv[a], pt, nm);
					spt.push_back(pt);
					sval.push_back(nm);
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
					el[0] = ecount * ns * ns + a * ns + b;
					el[1] = ecount * ns * ns + a * ns + b + 1;
					el[2] = ecount * ns * ns + (a + 1) * ns + b + 1;
					el[3] = ecount * ns * ns + (a + 1) * ns + b;
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
	}

	string fname = fn + "_geom.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i < spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
		}
		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i]<<"\n";
		//}
		//fout << "\nPOINT_DATA " << sval.size() << "\nNORMALS Normal FLOAT\n";
		//for (uint i = 0; i<sval.size(); i++)
		//{
		//	fout << sval[i][0] << " " << sval[i][1] << " " << sval[i][2] << "\n";
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
}

void TruncatedTspline::GetSubMatrixDual(int nv, vector<MatrixXd> &smat)
{
	//assume two irregular elements at least three rings!
	smat.resize(4);
	smat[0] = MatrixXd::Zero(nv + 5, nv + 5);
	smat[1] = MatrixXd::Zero(nv + 5, 16);
	smat[2] = MatrixXd::Zero(nv + 5, 16);
	smat[3] = MatrixXd::Zero(nv + 5, 16);

	double ku0_tmp[8] = {-2., -1., 0., 0., 1., 2., 3., 4.};
	double ku1_tmp[14] = {-2., -1.5, -1., -.5, 0., 0., .5, 1., 1.5, 2., 2.5, 3., 3.5, 4.};
	vector<double> ku0(ku0_tmp, ku0_tmp + 8);
	vector<double> ku1(ku1_tmp, ku1_tmp + 14);
	vector<vector<double>> tmat;
	TMatrix(ku0, ku1, 3, tmat);
	vector<array<int, 2>> ist(nv + 5); //corresponding to old patch
	for (uint i = 0; i < ist.size(); i++)
	{
		ist[i][0] = -1;
		ist[i][1] = -1;
	}
	ist[0][0] = 1;
	ist[0][1] = 1;
	ist[1][0] = 0;
	ist[1][1] = 1;
	int itmp[6][2] = {{1, 0}, {2, 0}, {2, 1}, {2, 2}, {1, 2}, {0, 2}};
	for (int i = 0; i < 6; i++)
	{
		ist[nv - 1 + i][0] = itmp[i][0];
		ist[nv - 1 + i][1] = itmp[i][1];
	}
	int jst[4][2] = {{2, 2}, {3, 2}, {3, 3}, {2, 3}}; //corresponding to three regular sub-patch
	vector<array<int, 2>> kst(nv + 5);				  //corresponding to irregular sub-patch
	for (uint i = 0; i < kst.size(); i++)
	{
		kst[i][0] = -1;
		kst[i][1] = -1;
	}
	kst[0][0] = 3;
	kst[0][1] = 3;
	kst[1][0] = 2;
	kst[1][1] = 3;
	int ktmp[6][2] = {{3, 2}, {4, 2}, {4, 3}, {4, 4}, {3, 4}, {2, 4}};
	for (int i = 0; i < 6; i++)
	{
		kst[nv - 1 + i][0] = ktmp[i][0];
		kst[nv - 1 + i][1] = ktmp[i][1];
	}
	for (int i = 0; i < nv + 5; i++)
	{
		if (ist[i][0] != -1 && ist[i][1] != -1)
		{
			for (int j = 0; j < 16; j++) //three regular sub-patches
			{
				int iloc(j % 4), jloc(j / 4);
				smat[1](i, j) = tmat[jst[1][0] + iloc][ist[i][0]] * tmat[jst[1][1] + jloc][ist[i][1]];
				smat[2](i, j) = tmat[jst[2][0] + iloc][ist[i][0]] * tmat[jst[2][1] + jloc][ist[i][1]];
				smat[3](i, j) = tmat[jst[3][0] + iloc][ist[i][0]] * tmat[jst[3][1] + jloc][ist[i][1]];
			}
			for (int k = 0; k < nv + 5; k++) //smat[0]
			{
				if (kst[k][0] != -1 && kst[k][1] != -1)
				{
					smat[0](i, k) = tmat[kst[k][0]][ist[i][0]] * tmat[kst[k][1]][ist[i][1]];
				}
			}
		}
	}
	//modify smat[1], smat[2], smat[3]
	double coef[3] = {1. / (4. * double(nv)) + 1. / 2., 1. / (4. * double(nv)) + 1. / 8., 1. / (4. * double(nv))};
	vector<double> c0(nv), c_1(nv), c1(nv);
	c0[0] = coef[0];
	c0[1] = coef[1];
	c0[nv - 1] = coef[1];
	for (int i = 2; i < nv - 1; i++)
		c0[i] = coef[2];
	c_1[0] = coef[1];
	c_1[nv - 1] = coef[0];
	c_1[nv - 2] = coef[1];
	for (int i = 1; i < nv - 2; i++)
		c_1[i] = coef[2];
	c1[0] = coef[1];
	c1[1] = coef[0];
	c1[2] = coef[1];
	for (int i = 3; i < nv; i++)
		c1[i] = coef[2];
	for (int i = 0; i < nv; i++)
	{
		smat[1](i, 4) = c0[i];
		smat[1](i, 0) = c_1[i];
		smat[2](i, 0) = c0[i];
		smat[3](i, 1) = c0[i];
		smat[3](i, 0) = c1[i];
	}
	//modify the first nv columns of smat[0]
	for (int i = 0; i < nv; i++)
	{
		smat[0](i, i) = coef[0];
		smat[0]((i + 1) % nv, i) = coef[1];
		smat[0]((i + nv - 1) % nv, i) = coef[1];
		for (int j = 2; j < nv - 1; j++)
		{
			smat[0]((i + j) % nv, i) = coef[2];
		}
	}
}

void TruncatedTspline::GetSubMatrixDual_All(int nv, MatrixXd &amat)
{
	//assume two irregular elements at least three rings!
	int nbf(nv + 12);
	amat = MatrixXd::Zero(nbf, nbf);
	double ku0_tmp[8] = {-2., -1., 0., 0., 1., 2., 3., 4.};
	double ku1_tmp[14] = {-2., -1.5, -1., -.5, 0., 0., .5, 1., 1.5, 2., 2.5, 3., 3.5, 4.};
	vector<double> ku0(ku0_tmp, ku0_tmp + 8);
	vector<double> ku1(ku1_tmp, ku1_tmp + 14);
	vector<vector<double>> tmat;
	TMatrix(ku0, ku1, 3, tmat);
	vector<array<int, 2>> ist(nbf); //corresponding to old patch
	for (uint i = 0; i < ist.size(); i++)
	{
		ist[i][0] = -1;
		ist[i][1] = -1;
	}
	ist[0][0] = 1;
	ist[0][1] = 1;
	ist[1][0] = 0;
	ist[1][1] = 1;
	int itmp[13][2] = {{1, 0}, {2, 0}, {2, 1}, {2, 2}, {1, 2}, {0, 2}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {2, 3}, {1, 3}, {0, 3}};
	for (int i = 0; i < 13; i++)
	{
		ist[nv - 1 + i][0] = itmp[i][0];
		ist[nv - 1 + i][1] = itmp[i][1];
	}
	vector<array<int, 2>> kst(nbf); //corresponding to irregular sub-patch
	//int jst[4][2] = { { 2,2 },{ 3,2 },{ 3,3 },{ 2,3 } };//corresponding to new sub-patch
	for (uint i = 0; i < kst.size(); i++)
	{
		kst[i][0] = -1;
		kst[i][1] = -1;
	}
	kst[0][0] = 3;
	kst[0][1] = 3;
	kst[1][0] = 2;
	kst[1][1] = 3;
	int ktmp[13][2] = {{3, 2}, {4, 2}, {4, 3}, {4, 4}, {3, 4}, {2, 4}, {5, 2}, {5, 3}, {5, 4}, {5, 5}, {4, 5}, {3, 5}, {2, 5}};
	for (int i = 0; i < 13; i++)
	{
		kst[nv - 1 + i][0] = ktmp[i][0];
		kst[nv - 1 + i][1] = ktmp[i][1];
	}
	for (int i = 0; i < nbf; i++) //old
	{
		if (ist[i][0] != -1 && ist[i][1] != -1)
		{
			for (int k = 0; k < nbf; k++) //new
			{
				if (kst[k][0] != -1 && kst[k][1] != -1)
				{
					amat(i, k) = tmat[kst[k][0]][ist[i][0]] * tmat[kst[k][1]][ist[i][1]];
				}
			}
		}
	}

	double PI = 3.14159265358979;
	//modify nv*nv sub-matrix
	//double coef[3] = { 1. / (4.*double(nv)) + 1. / 2., 1. / (4.*double(nv)) + 1. / 8.,1. / (4.*double(nv)) };
	for (int i = 0; i < nv; i++) //old
	{
		int nb[2] = {(i + 1) % nv, (i - 1 + nv) % nv};
		//cout << nb[0] << " " << nb[1] << "\n";
		//cout << i << " " << (i - 1 + nv) % nv << "\n"; getchar();
		for (int j = 0; j < nv; j++) //new
		{
			//cout << j << " ";
			if (j == i)
			{
				//amat(i, j) = coef[0];
				//amat(i, j) = (double(nv) + 5.) / (4.*double(nv));
				amat(i, j) = (lamda * double(nv) + lamda + 2.) / (2. * double(nv));
				//cout << "case 1 " << amat(i, j) << "\n";
			}
			//else if (j == nb[0] || j == nb[1])
			//{
			//	amat(i, j) = coef[1];
			//	//cout << "case 2 " << amat(i, j) << "\n";
			//}
			else
			{
				//amat(i, j) = coef[2];
				//amat(i, j) = (3.+2.*cos(2.*double((j-i))*PI/double(nv)))/ (4.*double(nv));
				amat(i, j) = (2. - lamda + 2. * lamda * cos(2. * double((j - i)) * PI / double(nv))) / (2. * double(nv));
				//cout << "case 3 " << amat(i, j) << "\n";
			}
		}
		//getchar();
	}

	//update to achieve optimal convergence rates
	//for five points, P_i^{1,0} (nv+1), P_i^{1,1} (nv+2), P_i^{nv+3}, P_{i-1}^{0,1} (nv), P_{i+1}^{1,0} (nv+4)
	//where i represents element ID
	//	double lmdbar(1.-lamda);
	//	int target[5] = {nv,nv+1,nv+2,nv+3,nv+4};//new
	//	int pcnct[5][4] = {{nv-1,nv,0,nv+1},{0,nv+1,nv-1,nv},{0,nv+1,nv+3,nv+2},{0,nv+3,1,nv+4},{1,nv+4,0,nv+3}};
	//	double c[7] = {(lmdbar+1)*(lamda+1)/4., lamda*(lamda+1)/4., lmdbar*(lmdbar+1)/4., lmdbar*lamda/4.,
	//				(lmdbar+1)*(lmdbar+1)/4., lamda*(lmdbar+1)/4., lamda*lamda/4.};
	//	double cmat[5][4] = {{c[0],c[1],c[2],c[3]},{c[0],c[1],c[2],c[3]},{c[4],c[5],c[6],c[5]},
	//					  {c[0],c[1],c[2],c[3]},{c[0],c[1],c[2],c[3]}};
	//	for(int i=0; i<nbf; i++)
	//	{
	//		for(int j=0; j<5; j++)
	//		{
	//			amat(i,target[j]) = 0.;
	//		}
	//	}
	//	for(int i=0; i<5; i++)
	//	{
	//		for(int j=0; j<4; j++)
	//		{
	//			amat(pcnct[i][j], target[i]) = cmat[i][j];
	//		}
	//	}
}

void TruncatedTspline::GetSubMatrixDual_Bar(int nv, MatrixXd &abar)
{
	//assume two irregular elements at least three rings!
	int nbf(nv + 12), nbf1(nv + 12 + 9);
	abar = MatrixXd::Zero(nbf, nbf1);
	double ku0_tmp[8] = {-2., -1., 0., 0., 1., 2., 3., 4.};
	double ku1_tmp[14] = {-2., -1.5, -1., -.5, 0., 0., .5, 1., 1.5, 2., 2.5, 3., 3.5, 4.};
	vector<double> ku0(ku0_tmp, ku0_tmp + 8);
	vector<double> ku1(ku1_tmp, ku1_tmp + 14);
	vector<vector<double>> tmat;
	TMatrix(ku0, ku1, 3, tmat);
	vector<array<int, 2>> ist(nbf); //corresponding to old patch
	for (uint i = 0; i < ist.size(); i++)
	{
		ist[i][0] = -1;
		ist[i][1] = -1;
	}
	ist[0][0] = 1;
	ist[0][1] = 1;
	ist[1][0] = 0;
	ist[1][1] = 1;
	int itmp[13][2] = {{1, 0}, {2, 0}, {2, 1}, {2, 2}, {1, 2}, {0, 2}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {2, 3}, {1, 3}, {0, 3}};
	for (int i = 0; i < 13; i++)
	{
		ist[nv - 1 + i][0] = itmp[i][0];
		ist[nv - 1 + i][1] = itmp[i][1];
	}
	vector<array<int, 2>> kst(nbf1); //corresponding to new sub-patch
	for (uint i = 0; i < kst.size(); i++)
	{
		kst[i][0] = -1;
		kst[i][1] = -1;
	}
	kst[0][0] = 3;
	kst[0][1] = 3;
	kst[1][0] = 2;
	kst[1][1] = 3;
	int ktmp[22][2] = {{3, 2}, {4, 2}, {4, 3}, {4, 4}, {3, 4}, {2, 4}, {5, 2}, {5, 3}, {5, 4}, {5, 5}, {4, 5}, {3, 5}, {2, 5}, {6, 2}, {6, 3}, {6, 4}, {6, 5}, {6, 6}, {5, 6}, {4, 6}, {3, 6}, {2, 6}};
	for (int i = 0; i < 22; i++)
	{
		kst[nv - 1 + i][0] = ktmp[i][0];
		kst[nv - 1 + i][1] = ktmp[i][1];
	}
	for (int i = 0; i < nbf; i++) //old
	{
		if (ist[i][0] != -1 && ist[i][1] != -1)
		{
			for (int k = 0; k < nbf1; k++) //new
			{
				if (kst[k][0] != -1 && kst[k][1] != -1)
				{
					abar(i, k) = tmat[kst[k][0]][ist[i][0]] * tmat[kst[k][1]][ist[i][1]];
				}
			}
		}
	}

	double PI = 3.14159265358979;
	//modify nv*nv sub-matrix
	//double coef[3] = { 1. / (4.*double(nv)) + 1. / 2., 1. / (4.*double(nv)) + 1. / 8.,1. / (4.*double(nv)) };
	for (int i = 0; i < nv; i++) //old
	{
		int nb[2] = {(i + 1) % nv, (i - 1 + nv) % nv};
		//cout << nb[0] << " " << nb[1] << "\n";
		//cout << i << " " << (i - 1 + nv) % nv << "\n"; getchar();
		for (int j = 0; j < nv; j++) //new
		{
			//cout << j << " ";
			if (j == i)
			{
				//abar(i, j) = coef[0];
				//abar(i, j) = (double(nv) + 5.) / (4.*double(nv));
				abar(i, j) = (lamda * double(nv) + lamda + 2.) / (2. * double(nv));
				//cout << "case 1 " << amat(i, j) << "\n";
			}
			//else if (j == nb[0] || j == nb[1])
			//{
			//	abar(i, j) = coef[1];
			//	//cout << "case 2 " << amat(i, j) << "\n";
			//}
			else
			{
				//abar(i, j) = coef[2];
				//abar(i, j) = (3. + 2.*cos(2.*double((j - i))*PI / double(nv))) / (4.*double(nv));
				abar(i, j) = (2. - lamda + 2. * lamda * cos(2. * double((j - i)) * PI / double(nv))) / (2. * double(nv));
				//cout << "case 3 " << amat(i, j) << "\n";
			}
		}
		//getchar();
	}

	//update to achieve optimal convergence rates
	//for five points, P_i^{1,0} (nv+1), P_i^{1,1} (nv+2), P_i^{nv+3}, P_{i-1}^{0,1} (nv), P_{i+1}^{1,0} (nv+4)
	//where i represents element ID
	//	double lmdbar(1.-lamda);
	//	int target[5] = {nv,nv+1,nv+2,nv+3,nv+4};//new
	//	int pcnct[5][4] = {{nv-1,nv,0,nv+1},{0,nv+1,nv-1,nv},{0,nv+1,nv+3,nv+2},{0,nv+3,1,nv+4},{1,nv+4,0,nv+3}};
	//	double c[7] = {(lmdbar+1)*(lamda+1)/4., lamda*(lamda+1)/4., lmdbar*(lmdbar+1)/4., lmdbar*lamda/4.,
	//				   (lmdbar+1)*(lmdbar+1)/4., lamda*(lmdbar+1)/4., lamda*lamda/4.};
	//	double cmat[5][4] = {{c[0],c[1],c[2],c[3]},{c[0],c[1],c[2],c[3]},{c[4],c[5],c[6],c[5]},
	//						 {c[0],c[1],c[2],c[3]},{c[0],c[1],c[2],c[3]}};
	//	for(int i=0; i<nbf; i++)
	//	{
	//		for(int j=0; j<5; j++)
	//		{
	//			abar(i,target[j]) = 0.;
	//		}
	//	}
	//	for(int i=0; i<5; i++)
	//	{
	//		for(int j=0; j<4; j++)
	//		{
	//			abar(pcnct[i][j], target[i]) = cmat[i][j];
	//		}
	//	}
}

void TruncatedTspline::GetSubPatchKVDual(vector<array<double, 8>> &subku, vector<array<double, 8>> &subkv)
{
	subku.clear();
	subkv.clear();
	subku.resize(3);
	subkv.resize(3);
	double ku[9] = {-2., -1., 0., 0., 1., 2., 3., 4., 5.};
	int ist[3][2] = {{1, 0}, {1, 1}, {0, 1}};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			subku[i][j] = ku[ist[i][0] + j];
			subkv[i][j] = ku[ist[i][1] + j];
		}
	}
}

void TruncatedTspline::ElementBasisDual(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	if (tmesh[eid].act == 1)
	{
		if (tmesh[eid].type == 0 /*|| tmesh[eid].type == 1*/)
		{
			ElementBasis_Regular(eid, u, v, Nt, dNdt);
		}
		else if (tmesh[eid].type == 4)
		{
			ElementBasis_IrregularDual(eid, u, v, Nt, dNdt);
			//cout << "eid: " << eid << "\n";
			//double sum(0.);
			//for (uint i = 0; i < Nt.size(); i++)
			//{
			//	cout << Nt[i] << " ";
			//	sum += Nt[i];
			//}
			//cout << "\n";
			//cout << "sum: " << sum << "\n";
			//getchar();
		}
	}
}

void TruncatedTspline::ElementBasis_IrregularDual(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	//need to change variabe for Bezier, (u,v) -> (u_b,v_b)
	Nt.clear();
	dNdt.clear();
	Nt.resize(tmesh[eid].IEN.size(), 0.);
	dNdt.resize(tmesh[eid].IEN.size());

	double tol(1.e-3), eps(1.e-8);
	double u1 = max(tol, u);
	double v1 = max(tol, v);
	int n = floor(min(-log2(u1), -log2(v1))) + 1, subid;
	double pow2 = pow(2., double(n - 1));
	u1 *= pow2;
	v1 *= pow2;
	if (v1 < 0.5)
	{
		subid = 0;
		u1 = 2. * u1 - 1.;
		v1 = 2. * v1;
	}
	else if (u1 < 0.5)
	{
		subid = 2;
		u1 = 2. * u1;
		v1 = 2. * v1 - 1.;
	}
	else
	{
		subid = 1;
		u1 = 2. * u1 - 1.;
		v1 = 2. * v1 - 1.;
	}
	//cout << "subid: " << subid << "\n";
	//cout << "uv: " << u1 << " " << v1 << "\n";
	//evaluate first nv+5 and the remaining 7 functions separately
	//the last 7 functions
	vector<double> ku(5, 0.), kv(5, 0.), uval, vval;
	BSplineBasis bu, bv;
	int nv(tmesh[eid].IEN.size() - 12);
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
	double knot[9] = {-2., -1., 0., 0., 1., 2., 3., 4., 5.};
	int ist[3][2] = {{1, 0}, {1, 1}, {0, 1}};
	double shift[3][2] = {{-1., 0.}, {-1., -1.}, {0., -1.}};
	int loc(0);
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
			dNdt1[loc][0] = uval[1] * vval[0]; //problem! need scaling
			dNdt1[loc][1] = uval[0] * vval[1]; //problem!
			loc++;
		}
	}
	int Pk[3][16] = {{nv - 1, nv, nv + 5, nv + 12, 0, nv + 1, nv + 6, nv + 13, nv + 3, nv + 2, nv + 7, nv + 14, nv + 10, nv + 9, nv + 8, nv + 15},
					 {0, nv + 1, nv + 6, nv + 13, nv + 3, nv + 2, nv + 7, nv + 14, nv + 10, nv + 9, nv + 8, nv + 15, nv + 19, nv + 18, nv + 17, nv + 16},
					 {1, 0, nv + 1, nv + 6, nv + 4, nv + 3, nv + 2, nv + 7, nv + 11, nv + 10, nv + 9, nv + 8, nv + 20, nv + 19, nv + 18, nv + 17}};
	//cout << "sum: " << sum << "\n";
	MatrixXd mat(nv + 12, 16);
	for (int i = 0; i < nv + 12; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			mat(i, j) = tmesh[eid].abar(i, Pk[subid][j]);
		}
	}
	//cout << "eid: " << eid << "\n";
	//cout << "IEN size: " << tmesh[eid].IEN.size() << "\n";
	//cout << "nv: " << nv << "\n";
	//cout << "mat: " << mat.rows() << " " << mat.cols() << "\n";
	//cout << "smat[0]: " << tmesh[eid].smat[0].rows() << " " << tmesh[eid].smat[0].rows() << "\n";
	//getchar();
	if (n > 10)
	{
		cout << "n: " << n << "\n";
	}
	for (int i = 0; i < n - 1; i++)
	{
		mat = tmesh[eid].amat * mat;
	}
	for (int i = 0; i < nv + 12; i++)
	{
		Nt[i] = 0.;
		dNdt[i][0] = 0.;
		dNdt[i][1] = 0.;
		for (int j = 0; j < 16; j++)
		{
			Nt[i] += mat(i, j) * Nt1[j];
			dNdt[i][0] += mat(i, j) * dNdt1[j][0];
			dNdt[i][1] += mat(i, j) * dNdt1[j][1];
		}
	}
}

//void TruncatedTspline::ElementBasis_IrregularDual(int eid, double u, double v, vector<double>& Nt, vector<array<double, 2>>& dNdt)
//{
//	//need to change variabe for Bezier, (u,v) -> (u_b,v_b)
//	Nt.clear();
//	dNdt.clear();
//	Nt.resize(tmesh[eid].IEN.size(), 0.);
//	dNdt.resize(tmesh[eid].IEN.size());
//
//	double tol(1.e-3), eps(1.e-8);
//	double u1 = max(tol, u);
//	double v1 = max(tol, v);
//	int n = floor(min(-log2(u1), -log2(v1))) + 1, subid;
//	double pow2 = pow(2., double(n - 1));
//	u1 *= pow2; v1 *= pow2;
//	if (v1<0.5)
//	{
//		subid = 0; u1 = 2.*u1 - 1.; v1 = 2.*v1;
//	}
//	else if (u1<0.5)
//	{
//		subid = 2; u1 = 2.*u1; v1 = 2.*v1 - 1.;
//	}
//	else
//	{
//		subid = 1; u1 = 2.*u1 - 1.; v1 = 2.*v1 - 1.;
//	}
//	//cout << "subid: " << subid << "\n";
//	//cout << "uv: " << u1 << " " << v1 << "\n";
//	//evaluate first nv+5 and the remaining 7 functions separately
//	//the last 7 functions
//	vector<double> ku(5, 0.), kv(5, 0.), uval, vval;
//	BSplineBasis bu, bv;
//	int nv(tmesh[eid].IEN.size() - 12);
//	for (uint i = nv; i < tmesh[eid].IEN.size(); i++)
//	{
//		ku.assign(tmesh[eid].patch_ku[i].begin(), tmesh[eid].patch_ku[i].end());
//		kv.assign(tmesh[eid].patch_kv[i].begin(), tmesh[eid].patch_kv[i].end());
//		bu.Set(3, ku);
//		bv.Set(3, kv);
//		bu.BasisFunction(0, u, 1, uval);
//		bv.BasisFunction(0, v, 1, vval);
//		Nt[i] = uval[0] * vval[0];
//		dNdt[i][0] = uval[1] * vval[0];
//		dNdt[i][1] = uval[0] * vval[1];
//	}
//	//the first nv functions
//	//cout << n << " " << subid << " " << u1 << " " << v1 << "\n";
//	//double sum(0.);
//	double Nt1[16];
//	double dNdt1[16][2];
//	double subku[3][8];
//	double subkv[3][8];
//	double knot[9] = { -2.,-1.,0.,0.,1.,2.,3.,4.,5. };
//	int ist[3][2] = { { 1,0 },{ 1,1 },{ 0,1 } };
//	double shift[3][2] = { { -1.,0. },{ -1.,-1. },{ 0.,-1. } };
//	int loc(0);
//	for (int j = 0; j < 4; j++)
//	{
//		for (int i = 0; i < 4; i++)
//		{
//			for (int k = 0; k < 5; k++)
//			{
//				ku[k] = knot[ist[subid][0] + i + k] + shift[subid][0];
//				kv[k] = knot[ist[subid][1] + j + k] + shift[subid][1];
//			}
//			bu.Set(3, ku);
//			bv.Set(3, kv);
//			bu.BasisFunction(0, u1, 1, uval);
//			bv.BasisFunction(0, v1, 1, vval);
//			Nt1[loc] = uval[0] * vval[0];
//			//cout << Nt1[loc] << "\n";
//			//sum += Nt1[loc];
//			dNdt1[loc][0] = uval[1] * vval[0];//problem!
//			dNdt1[loc][1] = uval[0] * vval[1];//problem!
//			loc++;
//		}
//	}
//	//cout << "sum: " << sum << "\n";
//	MatrixXd mat = tmesh[eid].smat[subid + 1];
//	//cout << "eid: " << eid << "\n";
//	//cout << "IEN size: " << tmesh[eid].IEN.size() << "\n";
//	//cout << "nv: " << nv << "\n";
//	//cout << "mat: " << mat.rows() << " " << mat.cols() << "\n";
//	//cout << "smat[0]: " << tmesh[eid].smat[0].rows() << " " << tmesh[eid].smat[0].rows() << "\n";
//	//getchar();
//	if (n > 10)
//	{
//		cout << "n: " << n << "\n";
//	}
//	for (int i = 0; i < n - 1; i++)
//	{
//		mat = tmesh[eid].smat[0] * mat;
//	}
//	for (int i = 0; i < nv; i++)
//	{
//		Nt[i] = 0.;
//		dNdt[i][0] = 0.; dNdt[i][1] = 0.;
//		for (int j = 0; j < 16; j++)
//		{
//			Nt[i] += mat(i, j)*Nt1[j];
//			dNdt[i][0] += mat(i, j)*dNdt1[j][0];
//			dNdt[i][1] += mat(i, j)*dNdt1[j][1];
//		}
//	}
//}

void TruncatedTspline::GlobalRefineDual()
{
	nel_old = tmesh.size();
	InitializeRefine();
	for (uint eid = 0; eid < nel_old; eid++)
	{
		if (tmesh[eid].act == 1)
		{
			if (tmesh[eid].type == 0)
			{
				ElementRefineReg_Dual(eid);
			}
			else if (tmesh[eid].type == 4)
			{
				ElementRefineIrr_Dual(eid);
			}
			if (tmesh[eid].focus == 1 && (tmesh[eid].type == 0 || tmesh[eid].type == 4))
			{
				for (int i = 0; i < 4; i++)
				{
					if (tmesh[eid].chd[i] != -1)
					{
						tmesh[tmesh[eid].chd[i]].focus = 1;
					}
				}
			}
		}
	}
	//additional refinement of type 2 boundary element
	for (uint eid = 0; eid < nel_old; eid++)
	{
		if (tmesh[eid].act == 1 && tmesh[eid].type == 2)
		{
			ElementRefineBnd_Dual(eid);
		}
	}

	UpdateCP();
	BuildTsplinesDual();
}

void TruncatedTspline::InitializeRefine()
{
	for (uint i = 0; i < cp.size(); i++)
	{
		cp[i].coortmp[0] = 0.;
		cp[i].coortmp[1] = 0.;
		cp[i].coortmp[2] = 0.;
		cp[i].wtmp = 0.;
		cp[i].update = 0;
	}
}

void TruncatedTspline::ElementRefineReg_Dual(int eid)
{
	vector<double> pku1, pkv1;
	int ist[2] = {0, 0};
	for (uint i = 0; i < tmesh[eid].pku.size() - 1; i++)
	{
		pku1.push_back(tmesh[eid].pku[i]);
		if (i == 3)
		{
			ist[0] = pku1.size() - 4;
		}
		if (tmesh[eid].pku[i] < tmesh[eid].pku[i + 1])
		{
			pku1.push_back((tmesh[eid].pku[i] + tmesh[eid].pku[i + 1]) / 2.);
		}
	}
	pku1.push_back(tmesh[eid].pku.back());
	for (uint i = 0; i < tmesh[eid].pkv.size() - 1; i++)
	{
		pkv1.push_back(tmesh[eid].pkv[i]);
		if (i == 3)
		{
			ist[1] = pkv1.size() - 4;
		}
		if (tmesh[eid].pkv[i] < tmesh[eid].pkv[i + 1])
		{
			pkv1.push_back((tmesh[eid].pkv[i] + tmesh[eid].pkv[i + 1]) / 2.);
		}
	}
	pkv1.push_back(tmesh[eid].pkv.back());
	vector<vector<double>> Tu, Tv;
	TMatrix(tmesh[eid].pku, pku1, 3, Tu);
	TMatrix(tmesh[eid].pkv, pkv1, 3, Tv);
	vector<Vertex> pnew(25);
	vector<int> pid(25, 0);
	int loc0(0), loc1(0);
	double ctmp;
	for (int j1 = 0; j1 < 5; j1++)
	{
		for (int i1 = 0; i1 < 5; i1++) //new
		{
			loc0 = 0;
			for (int j0 = 0; j0 < 4; j0++)
			{
				for (int i0 = 0; i0 < 4; i0++) //old
				{
					ctmp = Tu[ist[0] + i1][i0] * Tv[ist[1] + j1][j0];
					if (ctmp != 0.)
					{
						pnew[loc1].coortmp[0] += ctmp * cp[tmesh[eid].IEN[loc0]].coor[0];
						pnew[loc1].coortmp[1] += ctmp * cp[tmesh[eid].IEN[loc0]].coor[1];
						pnew[loc1].coortmp[2] += ctmp * cp[tmesh[eid].IEN[loc0]].coor[2];
					}
					loc0++;
				}
			}
			pnew[loc1].update = 1;
			loc1++;
		}
	}
	int vtloc[4] = {6, 8, 18, 16};
	int edloc[4] = {7, 13, 17, 11};
	int fcloc(12);
	int edbloc[4][3] = {{1, 2, 3}, {9, 14, 19}, {23, 22, 21}, {15, 10, 5}};
	cp.push_back(pnew[fcloc]);
	pid[fcloc] = cp.size() - 1;
	for (int i = 0; i < 4; i++)
	{
		if (cp[tmesh[eid].cnct[i]].update == 0)
		{
			cp[tmesh[eid].cnct[i]].coortmp[0] = pnew[vtloc[i]].coortmp[0];
			cp[tmesh[eid].cnct[i]].coortmp[1] = pnew[vtloc[i]].coortmp[1];
			cp[tmesh[eid].cnct[i]].coortmp[2] = pnew[vtloc[i]].coortmp[2];
			cp[tmesh[eid].cnct[i]].update = 1;
		}
		pid[vtloc[i]] = tmesh[eid].cnct[i];
		//if (tmedge[tmesh[eid].edge[i]].act == 1)
		if (tmedge[tmesh[eid].edge[i]].midpt == -1)
		{
			cp.push_back(pnew[edloc[i]]);
			tmedge[tmesh[eid].edge[i]].midpt = cp.size() - 1;
			//tmedge[tmesh[eid].edge[i]].act = 0;
		}
		pid[edloc[i]] = tmedge[tmesh[eid].edge[i]].midpt;
	}
	//find if there exist neighboring boundary element, add element
	//int nbfc[4][2] = { {-1,0},{ -1,0 },{ -1,0 },{ -1,0 } };
	for (int i = 0; i < 4; i++)
	{
		int fcid(tmedge[tmesh[eid].edge[i]].face[0]);
		if (fcid == eid)
		{
			fcid = tmedge[tmesh[eid].edge[i]].face[1];
		}
		if (tmesh[fcid].type == 2 && tmesh[fcid].act == 1)
		{
			//nbfc[i][0] = fcid;
			int *it = find(tmesh[fcid].cnct, tmesh[fcid].cnct + 4, tmesh[eid].cnct[i]);
			int loc(it - tmesh[fcid].cnct);
			loc = (loc + 1) % 4;
			//nbfc[i][1] = loc;
			if (cp[tmesh[fcid].cnct[loc]].update == 0)
			{
				cp[tmesh[fcid].cnct[loc]].coortmp[0] = pnew[edbloc[i][0]].coortmp[0];
				cp[tmesh[fcid].cnct[loc]].coortmp[1] = pnew[edbloc[i][0]].coortmp[1];
				cp[tmesh[fcid].cnct[loc]].coortmp[2] = pnew[edbloc[i][0]].coortmp[2];
				cp[tmesh[fcid].cnct[loc]].update = 1;
			}
			//if (tmedge[tmesh[fcid].edge[loc]].act == 1)
			if (tmedge[tmesh[fcid].edge[loc]].midpt == -1)
			{
				cp.push_back(pnew[edbloc[i][1]]);
				tmedge[tmesh[fcid].edge[loc]].midpt = cp.size() - 1;
				//tmedge[tmesh[fcid].edge[loc]].act = 0;
			}
			loc = (loc + 1) % 4;
			if (cp[tmesh[fcid].cnct[loc]].update == 0)
			{
				cp[tmesh[fcid].cnct[loc]].coortmp[0] = pnew[edbloc[i][2]].coortmp[0];
				cp[tmesh[fcid].cnct[loc]].coortmp[1] = pnew[edbloc[i][2]].coortmp[1];
				cp[tmesh[fcid].cnct[loc]].coortmp[2] = pnew[edbloc[i][2]].coortmp[2];
				cp[tmesh[fcid].cnct[loc]].update = 1;
			}
		}
	}
	//build new elements and edges
	int fccnct[4][4] = {{6, 7, 12, 11}, {7, 8, 13, 12}, {12, 13, 18, 17}, {11, 12, 17, 16}};
	int edcnct[12][2] = {{6, 7}, {7, 8}, {7, 12}, {8, 13}, {13, 18}, {13, 12}, {18, 17}, {17, 16}, {17, 12}, {16, 11}, {11, 6}, {11, 12}};
	vector<Edge> ednew(12);
	int edid[12];
	for (int i = 0; i < 4; i++)
	{
		int iloc[3] = {3 * i, 3 * i + 1, 3 * i + 2};
		for (int j = 0; j < 3; j++)
		{
			ednew[3 * i + j].pt[0] = pid[edcnct[3 * i + j][0]];
			ednew[3 * i + j].pt[1] = pid[edcnct[3 * i + j][1]];
		}
		//if (tmedge[tmesh[eid].edge[i]].act == 1)
		if (tmedge[tmesh[eid].edge[i]].chd[0] == -1)
		{
			tmedge.push_back(ednew[iloc[0]]);
			tmedge[tmesh[eid].edge[i]].chd[0] = tmedge.size() - 1;
			tmedge.push_back(ednew[iloc[1]]);
			tmedge[tmesh[eid].edge[i]].chd[1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[i]].act = 0;
		}
		edid[iloc[0]] = tmedge[tmesh[eid].edge[i]].chd[0];
		edid[iloc[1]] = tmedge[tmesh[eid].edge[i]].chd[1];
		if (tmedge[edid[iloc[0]]].pt[0] != tmesh[eid].cnct[i] && tmedge[edid[iloc[0]]].pt[1] != tmesh[eid].cnct[i])
		{
			edid[iloc[0]] = tmedge[tmesh[eid].edge[i]].chd[1];
			edid[iloc[1]] = tmedge[tmesh[eid].edge[i]].chd[0];
		}
		tmedge.push_back(ednew[iloc[2]]);
		edid[iloc[2]] = tmedge.size() - 1;
	}
	int fced[4][4] = {{0, 2, 11, 10}, {1, 3, 5, 2}, {5, 4, 6, 8}, {11, 8, 7, 9}};
	vector<Element> enew(4);
	for (int i = 0; i < 4; i++)
	{
		enew[i].act = 1;
		for (int j = 0; j < 4; j++)
		{
			enew[i].cnct[j] = pid[fccnct[i][j]];
			enew[i].edge[j] = edid[fced[i][j]];
		}
		tmesh.push_back(enew[i]);
		tmesh[eid].chd[i] = tmesh.size() - 1;
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::ElementRefineIrr_Dual(int eid)
{
	int nv(tmesh[eid].IEN.size() - 12);
	vector<Vertex> pnew(tmesh[eid].IEN.size());
	vector<int> pid(tmesh[eid].IEN.size());
	for (uint i = 0; i < pnew.size(); i++)
	{
		pnew[i].update = 1;
		for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
		{
			double ctmp(tmesh[eid].amat(j, i));
			if (ctmp > 0.)
			{
				pnew[i].coortmp[0] += ctmp * cp[tmesh[eid].IEN[j]].coor[0];
				pnew[i].coortmp[1] += ctmp * cp[tmesh[eid].IEN[j]].coor[1];
				pnew[i].coortmp[2] += ctmp * cp[tmesh[eid].IEN[j]].coor[2];
			}
		}
	}
	//cout << pnew[0].coortmp[0] << " " << pnew[0].coortmp[1] << "\n";
	//getchar();
	int vtloc[4] = {0, nv + 6, nv + 8, nv + 10};
	int edloc[4] = {nv + 1, nv + 7, nv + 9, nv + 3};
	int fcloc(nv + 2);
	int edbloc[4][3] = {{nv - 1, nv, nv + 5}, {nv + 13, nv + 14, nv + 15}, {nv + 17, nv + 18, nv + 19}, {nv + 11, nv + 4, 1}}; //edbloc[1] and edbloc[2] not useful
	cp.push_back(pnew[fcloc]);
	pid[fcloc] = cp.size() - 1;
	for (int i = 0; i < 4; i++)
	{
		if (cp[tmesh[eid].cnct[i]].update == 0)
		{
			cp[tmesh[eid].cnct[i]].coortmp[0] = pnew[vtloc[i]].coortmp[0];
			cp[tmesh[eid].cnct[i]].coortmp[1] = pnew[vtloc[i]].coortmp[1];
			cp[tmesh[eid].cnct[i]].coortmp[2] = pnew[vtloc[i]].coortmp[2];
			cp[tmesh[eid].cnct[i]].update = 1;
		}
		pid[vtloc[i]] = tmesh[eid].cnct[i];
		//if (tmedge[tmesh[eid].edge[i]].act == 1)
		if (tmedge[tmesh[eid].edge[i]].midpt == -1)
		{
			cp.push_back(pnew[edloc[i]]);
			tmedge[tmesh[eid].edge[i]].midpt = cp.size() - 1;
			//tmedge[tmesh[eid].edge[i]].act = 0;
		}
		pid[edloc[i]] = tmedge[tmesh[eid].edge[i]].midpt;
	}
	//cout << cp[tmesh[eid].cnct[0]].coortmp[0] << " " << cp[tmesh[eid].cnct[0]].coortmp[1] << "\n";
	//getchar();
	//find if there exist neighboring boundary element, add element
	for (int i = 0; i < 4; i++)
	{
		int fcid(tmedge[tmesh[eid].edge[i]].face[0]);
		if (fcid == eid)
		{
			fcid = tmedge[tmesh[eid].edge[i]].face[1];
		}
		if (tmesh[fcid].type == 2 && tmesh[fcid].act == 1)
		{
			int *it = find(tmesh[fcid].cnct, tmesh[fcid].cnct + 4, tmesh[eid].cnct[i]);
			int loc(it - tmesh[fcid].cnct);
			loc = (loc + 1) % 4;
			if (cp[tmesh[fcid].cnct[loc]].update == 0)
			{
				cp[tmesh[fcid].cnct[loc]].coortmp[0] = pnew[edbloc[i][0]].coortmp[0];
				cp[tmesh[fcid].cnct[loc]].coortmp[1] = pnew[edbloc[i][0]].coortmp[1];
				cp[tmesh[fcid].cnct[loc]].coortmp[2] = pnew[edbloc[i][0]].coortmp[2];
				cp[tmesh[fcid].cnct[loc]].update = 1;
			}
			//if (tmedge[tmesh[fcid].edge[loc]].act == 1)
			if (tmedge[tmesh[fcid].edge[loc]].midpt == -1)
			{
				cp.push_back(pnew[edbloc[i][1]]);
				tmedge[tmesh[fcid].edge[loc]].midpt = cp.size() - 1;
				//tmedge[tmesh[fcid].edge[loc]].act = 0;
			}
			loc = (loc + 1) % 4;
			if (cp[tmesh[fcid].cnct[loc]].update == 0)
			{
				cp[tmesh[fcid].cnct[loc]].coortmp[0] = pnew[edbloc[i][2]].coortmp[0];
				cp[tmesh[fcid].cnct[loc]].coortmp[1] = pnew[edbloc[i][2]].coortmp[1];
				cp[tmesh[fcid].cnct[loc]].coortmp[2] = pnew[edbloc[i][2]].coortmp[2];
				cp[tmesh[fcid].cnct[loc]].update = 1;
			}
		}
	}
	//build new elements and edges
	int fccnct[4][4] = {{0, nv + 1, nv + 2, nv + 3}, {nv + 1, nv + 6, nv + 7, nv + 2}, {nv + 2, nv + 7, nv + 8, nv + 9}, {nv + 3, nv + 2, nv + 9, nv + 10}};
	int edcnct[12][2] = {{0, nv + 1}, {nv + 1, nv + 6}, {nv + 1, nv + 2}, {nv + 6, nv + 7}, {nv + 7, nv + 8}, {nv + 7, nv + 2}, {nv + 8, nv + 9}, {nv + 9, nv + 10}, {nv + 9, nv + 2}, {nv + 10, nv + 3}, {nv + 3, 0}, {nv + 3, nv + 2}};
	vector<Edge> ednew(12);
	int edid[12];
	for (int i = 0; i < 4; i++)
	{
		int iloc[3] = {3 * i, 3 * i + 1, 3 * i + 2};
		for (int j = 0; j < 3; j++)
		{
			ednew[3 * i + j].pt[0] = pid[edcnct[3 * i + j][0]];
			ednew[3 * i + j].pt[1] = pid[edcnct[3 * i + j][1]];
		}
		//if (tmedge[tmesh[eid].edge[i]].act == 1)
		if (tmedge[tmesh[eid].edge[i]].chd[0] == -1)
		{
			tmedge.push_back(ednew[iloc[0]]);
			tmedge[tmesh[eid].edge[i]].chd[0] = tmedge.size() - 1;
			tmedge.push_back(ednew[iloc[1]]);
			tmedge[tmesh[eid].edge[i]].chd[1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[i]].act = 0;
		}
		edid[iloc[0]] = tmedge[tmesh[eid].edge[i]].chd[0];
		edid[iloc[1]] = tmedge[tmesh[eid].edge[i]].chd[1];
		if (tmedge[edid[iloc[0]]].pt[0] != tmesh[eid].cnct[i] && tmedge[edid[iloc[0]]].pt[1] != tmesh[eid].cnct[i])
		{
			edid[iloc[0]] = tmedge[tmesh[eid].edge[i]].chd[1];
			edid[iloc[1]] = tmedge[tmesh[eid].edge[i]].chd[0];
		}
		tmedge.push_back(ednew[iloc[2]]);
		edid[iloc[2]] = tmedge.size() - 1;
	}
	int fced[4][4] = {{0, 2, 11, 10}, {1, 3, 5, 2}, {5, 4, 6, 8}, {11, 8, 7, 9}};
	vector<Element> enew(4);
	for (int i = 0; i < 4; i++)
	{
		enew[i].act = 1;
		if (i == 0)
			enew[i].type = 4;
		for (int j = 0; j < 4; j++)
		{
			enew[i].cnct[j] = pid[fccnct[i][j]];
			enew[i].edge[j] = edid[fced[i][j]];
		}
		tmesh.push_back(enew[i]);
		tmesh[eid].chd[i] = tmesh.size() - 1;
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::ElementRefineBnd_Dual(int eid)
{
	int loc(0);
	if (tmedge[tmesh[eid].edge[0]].len < 1.e-8)
	{
		loc = 1;
	}
	int iloc[4] = {loc, (loc + 1) % 4, (loc + 2) % 4, (loc + 3) % 4};
	int pid[6] = {tmesh[eid].cnct[iloc[0]], tmesh[eid].cnct[iloc[1]], tmesh[eid].cnct[iloc[2]], tmesh[eid].cnct[iloc[3]],
				  tmedge[tmesh[eid].edge[iloc[0]]].midpt, tmedge[tmesh[eid].edge[iloc[2]]].midpt};
	//build necessary edges
	int edcnct[5][2] = {{0, 4}, {4, 1}, {2, 5}, {5, 3}, {4, 5}};
	int edid[7];
	int jloc[2] = {iloc[0], iloc[2]};
	for (int i = 0; i < 2; i++)
	{
		int i1[2] = {2 * i, 2 * i + 1};
		if (tmedge[tmesh[eid].edge[jloc[i]]].chd[0] == -1)
		{
			Edge edtmp[2];
			edtmp[0].pt[0] = pid[edcnct[i1[0]][0]];
			edtmp[0].pt[1] = pid[edcnct[i1[0]][1]];
			edtmp[1].pt[0] = pid[edcnct[i1[1]][0]];
			edtmp[1].pt[1] = pid[edcnct[i1[1]][1]];
			tmedge.push_back(edtmp[0]);
			tmedge[tmesh[eid].edge[jloc[i]]].chd[0] = tmedge.size() - 1;
			tmedge.push_back(edtmp[1]);
			tmedge[tmesh[eid].edge[jloc[i]]].chd[1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[jloc[i]]].act = 0;
		}
		edid[i1[0]] = tmedge[tmesh[eid].edge[jloc[i]]].chd[0];
		edid[i1[1]] = tmedge[tmesh[eid].edge[jloc[i]]].chd[1];
		if (tmedge[edid[i1[0]]].pt[0] != tmesh[eid].cnct[jloc[i]] && tmedge[edid[i1[0]]].pt[1] != tmesh[eid].cnct[jloc[i]])
		{
			edid[i1[0]] = tmedge[tmesh[eid].edge[jloc[i]]].chd[1];
			edid[i1[1]] = tmedge[tmesh[eid].edge[jloc[i]]].chd[0];
		}
	}
	Edge ednew0;
	ednew0.pt[0] = pid[edcnct[4][0]];
	ednew0.pt[1] = pid[edcnct[4][1]];
	ednew0.len = 0.;
	tmedge.push_back(ednew0);
	edid[4] = tmedge.size() - 1;
	edid[5] = tmesh[eid].edge[iloc[1]];
	edid[6] = tmesh[eid].edge[iloc[3]];

	int fcvt[2][4] = {{0, 4, 5, 3}, {4, 1, 2, 5}};
	int fced[2][4] = {{0, 4, 3, 6}, {1, 5, 2, 4}};
	Element enew[2];
	for (int i = 0; i < 2; i++)
	{
		enew[i].act = 1;
		enew[i].type = 2;
		for (int j = 0; j < 4; j++)
		{
			enew[i].cnct[iloc[j]] = pid[fcvt[i][j]];
			enew[i].edge[iloc[j]] = edid[fced[i][j]];
		}
		tmesh.push_back(enew[i]);
		tmesh[eid].chd[i] = tmesh.size() - 1;
	}

	tmesh[eid].act = 0;
}

void TruncatedTspline::UpdateCP()
{
	//int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		//cp[i].h2a = -1;
		if (cp[i].update == 1)
		{
			cp[i].coor[0] = cp[i].coortmp[0];
			cp[i].coor[1] = cp[i].coortmp[1];
			cp[i].coor[2] = cp[i].coortmp[2];
		}
		//if (cp[i].act == 1)
		//{
		//	cp[i].h2a = count++;
		//}
	}
}

void TruncatedTspline::BezierExtractDual(vector<BezierElement> &bzmesh)
{
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			cp[i].h2a = count++;
		}
		else
		{
			cp[i].h2a = -1;
		}
	}
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1)
		{
			if (tmesh[eid].type == 0)
			{
				BezierElement bzel;
				bzel.prt = eid;
				bzel.focus = tmesh[eid].focus;
				ElementBezierExtractReg(eid, bzel);
				bzmesh.push_back(bzel);
			}
			else if (tmesh[eid].type == 4)
			{
				BezierElement bzel;
				bzel.prt = eid;
				bzel.focus = tmesh[eid].focus;
				ElementBezierExtractIrr(eid, bzel);
				bzmesh.push_back(bzel);
			}
		}
	}
}

void TruncatedTspline::ElementBezierExtractReg(int eid, BezierElement &bzel)
{
	bzel.type = 0;
	bzel.cmat.resize(16);
	vector<double> pku1, pkv1;
	array<double, 2> bzkt = {0., 1.};
	BezierInsertKnots(tmesh[eid].pku, bzkt, pku1);
	BezierInsertKnots(tmesh[eid].pkv, bzkt, pkv1);
	int ist[2] = {0, 0};
	for (uint i = 0; i < pku1.size() - 1; i++)
	{
		if (fabs(pku1[i] - bzkt[0]) < 1.e-8 && fabs(pku1[i + 1] - bzkt[1]) < 1.e-8)
		{
			ist[0] = i - 3;
		}
	}
	for (uint i = 0; i < pkv1.size() - 1; i++)
	{
		if (fabs(pkv1[i] - bzkt[0]) < 1.e-8 && fabs(pkv1[i + 1] - bzkt[1]) < 1.e-8)
		{
			ist[1] = i - 3;
		}
	}
	//for (uint i = 0; i < tmesh[eid].pku.size(); i++)
	//{
	//	cout << tmesh[eid].pku[i] << " ";
	//}
	//cout << "\n";
	//for (uint i = 0; i < pku1.size(); i++)
	//{
	//	cout << pku1[i] << " ";
	//}
	//cout << "\n";
	//cout << ist[0] << " " << ist[1] << "\n"; getchar();
	vector<vector<double>> Tu, Tv;
	TMatrix(tmesh[eid].pku, pku1, 3, Tu);
	TMatrix(tmesh[eid].pkv, pkv1, 3, Tv);
	int loc0(0), loc1(0);
	double ctmp;
	for (int j1 = 0; j1 < 4; j1++)
	{
		for (int i1 = 0; i1 < 4; i1++) //bezier
		{
			loc0 = 0;
			for (int j0 = 0; j0 < 4; j0++)
			{
				for (int i0 = 0; i0 < 4; i0++) //spline
				{
					ctmp = Tu[ist[0] + i1][i0] * Tv[ist[1] + j1][j0];
					bzel.cmat[loc0][loc1] = ctmp;
					if (ctmp != 0.)
					{
						bzel.pts[loc1][0] += ctmp * cp[tmesh[eid].IEN[loc0]].coor[0];
						bzel.pts[loc1][1] += ctmp * cp[tmesh[eid].IEN[loc0]].coor[1];
						bzel.pts[loc1][2] += ctmp * cp[tmesh[eid].IEN[loc0]].coor[2];
					}
					loc0++;
				}
			}
			loc1++;
		}
	}
	//bzel.IEN = tmesh[eid].IEN;
	bzel.IEN.resize(tmesh[eid].IEN.size());
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		bzel.IEN[i] = cp[tmesh[eid].IEN[i]].h2a;
		if (bzel.IEN[i] == -1)
		{
			cerr << "IEN -1!\n";
			getchar();
		}
	}
}

void TruncatedTspline::ElementBezierExtractIrr(int eid, BezierElement &bzel)
{
	bzel.type = 4;
	bzel.amat = tmesh[eid].amat;
	bzel.abar = tmesh[eid].abar;
	bzel.cp.resize(tmesh[eid].IEN.size());
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		bzel.cp[i][0] = cp[tmesh[eid].IEN[i]].coor[0];
		bzel.cp[i][1] = cp[tmesh[eid].IEN[i]].coor[1];
		bzel.cp[i][2] = cp[tmesh[eid].IEN[i]].coor[2];
	}
	bzel.IEN.resize(tmesh[eid].IEN.size());
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		bzel.IEN[i] = cp[tmesh[eid].IEN[i]].h2a;
		if (bzel.IEN[i] == -1)
		{
			cerr << "IEN -1!\n";
			getchar();
		}
	}
}

void TruncatedTspline::DirichletBC(vector<int> &IDBC, vector<double> &gh)
{
	//unit square, u=0 on BC
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			count++;
		}
	}
	IDBC.resize(count, -1);
	gh.resize(count, 0.);
	int count1(0), loc(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			if (fabs(cp[i].coor[0] - 0.) < 1.e-8 || fabs(cp[i].coor[0] - 1.) < 1.e-8 ||
				fabs(cp[i].coor[1] - 0.) < 1.e-8 || fabs(cp[i].coor[1] - 1.) < 1.e-8)
			{
				gh[loc] = 0.;
			}
			else
			{
				IDBC[loc] = count1++;
			}
			loc++;
		}
	}

	//for (uint i = 0; i < cp.size(); i++)
	//{
	//	if (cp[i].act == 1)
	//	{
	//		if (fabs(cp[i].coor[0] - 0.) < 1.e-8 || fabs(cp[i].coor[0] - 1.) < 1.e-8)
	//		{
	//			gh[loc] = 0.;
	//			gh[loc] = cp[i].coor[0];
	//		}
	//		else
	//		{
	//			IDBC[loc] = count1++;
	//		}
	//		loc++;
	//	}
	//}
}

void TruncatedTspline::DirichletBC_LeastSquare(vector<BezierElement> &bzmesh, vector<int> &IDBC1, vector<double> &gh1)
{
	IDBC1.clear();
	gh1.clear();
	vector<int> bflag(cp.size(), 0);
	for (uint i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].act == 1 && tmedge[i].face.size() == 1)
		{
			bflag[tmedge[i].pt[0]] = 1;
			bflag[tmedge[i].pt[1]] = 1;
		}
	}
	int count(0);
	vector<int> pid(cp.size(), -1);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
			pid[i] = count++;
	}
	vector<int> IDBC(count + bzcp.size(), -1); //assume no bzcp on the boundary
	IDBC1.resize(IDBC.size(), -1);
	vector<double> gh(count + bzcp.size(), 0.);
	count = 0;
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			if (bflag[i] == 1) //solving for boundary points
			{
				//if (count >= 51 && count <= 58)
				//{
				//	cout << i << "\n";
				//}
				IDBC[pid[i]] = count++;
			}
		}
	}
	//getchar();
	count = 0;
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			if (bflag[i] == 0)
			{
				IDBC1[pid[i]] = count++;
			}
		}
	}
	for (uint i0 = 0; i0 < bzmesh.size(); i0++)
	{
		int eid(bzmesh[i0].prt);
		if (eid == -1)
		{
			cerr << "Empty corresponding T-mesh element!\n";
			continue;
		}
		for (int i = 0; i < 4; i++)
		{
			int edid(tmesh[eid].edge[i]);
			int ednb(tmedge[edid].face[0]);
			if (ednb == eid)
				ednb = tmedge[edid].face[1];
			int *it = find(tmesh[ednb].edge, tmesh[ednb].edge + 4, edid);
			int loc(it - tmesh[ednb].edge);
			int edid1(tmesh[ednb].edge[(loc + 2) % 4]);
			if (tmedge[edid1].face.size() == 1)
			{
				bzmesh[i0].bc[i] = 1;
				//cout << eid << " " << i << "\n"; getchar();
			}
		}
	}

	LeastSquare ls;
	ls.SetProblem(IDBC, gh);
	ls.Run_ScalarFitting(bzmesh, "", gh1);
}

//////////////////////////////////////////////

void TruncatedTspline::RunCC(string fn_in, int nrf)
{
	ReadInput(fn_in);
	//Scale2Unit();
	InitialConnect_2();
	//VisualizeEdge("../io/test1/square1");
	//VisualizeQuad("../io/test1/square1");

	BuildCC();

	//int nrf(4);
	for (int i = 0; i < nrf; i++)
	{
		GlobalRefineCC();
	}

	//VisualizeEdge("../io/cctest1/square1");
	//VisualizeQuad("../io/cctest1/square1");
	//VisualizeGeomCC("../io/cctest1/square2");
}

void TruncatedTspline::BuildCC()
{
	UpdateConectDual();
	FindIENCC();
	UpdateIENCC();
}

void TruncatedTspline::FindIENCC()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		tmesh[eid].IENtmp.clear();
		tmesh[eid].patch_kutmp.clear();
		tmesh[eid].patch_kvtmp.clear();
		if (tmesh[eid].act == 1 && tmesh[eid].type == 0)
		{
			FindElementIEN_RegDual(eid);
		}
		else if (tmesh[eid].act == 1 && tmesh[eid].type == 4)
		{
			FindElementIEN_IrrCC(eid);
		}
	}
}

void TruncatedTspline::FindElementIEN_IrrCC(int eid)
{
	int nv(cp[tmesh[eid].cnct[0]].face.size()); //suppose 0 is the only EP
	int nbf(2 * nv + 8);
	tmesh[eid].IENtmp.resize(nbf);
	//tmesh[eid].patch_kutmp.resize(nbf);
	//tmesh[eid].patch_kvtmp.resize(nbf);
	//IEN
	int iloc0[4] = {0, 5, 4, 3};
	for (int i = 0; i < 4; i++)
	{
		tmesh[eid].IENtmp[iloc0[i]] = tmesh[eid].cnct[i];
	}
	int ednb[4] = {-1, -1, -1, -1};
	for (int i = 0; i < 4; i++)
	{
		ednb[i] = tmedge[tmesh[eid].edge[i]].face[0];
		if (ednb[i] == eid)
			ednb[i] = tmedge[tmesh[eid].edge[i]].face[1];
	}
	tmesh[eid].IENtmp[1] = tmesh[ednb[3]].cnct[3];
	int iloc[3][3] = {{6, 2 * nv + 1, 2 * nv + 2}, {2 * nv + 3, 2 * nv + 4, 2 * nv + 5}, {2 * nv + 6, 2 * nv + 7, 2}};
	for (int i = 1; i < 4; i++)
	{
		for (uint j = 0; j < cp[tmesh[eid].cnct[i]].face.size(); j++)
		{
			int fcid(cp[tmesh[eid].cnct[i]].face[j]);
			if (fcid != eid && fcid != ednb[i] && fcid != ednb[(i + 3) % 4])
			{
				int *it = find(tmesh[fcid].cnct, tmesh[fcid].cnct + 4, tmesh[eid].cnct[i]);
				int loc(it - tmesh[fcid].cnct);
				tmesh[eid].IENtmp[iloc[i - 1][0]] = tmesh[fcid].cnct[(loc + 1) % 4];
				tmesh[eid].IENtmp[iloc[i - 1][1]] = tmesh[fcid].cnct[(loc + 2) % 4];
				tmesh[eid].IENtmp[iloc[i - 1][2]] = tmesh[fcid].cnct[(loc + 3) % 4];
			}
		}
	}
	if (cp[tmesh[eid].cnct[0]].face.size() > 3)
	{
		int enow(tmedge[tmesh[ednb[0]].edge[0]].face[0]);
		if (enow == ednb[0])
			enow = tmedge[tmesh[ednb[0]].edge[0]].face[1];
		int count(0), loc(7);
		while (count < cp[tmesh[eid].cnct[0]].face.size() - 3)
		{
			tmesh[eid].IENtmp[loc] = tmesh[enow].cnct[3];
			loc++;
			tmesh[eid].IENtmp[loc] = tmesh[enow].cnct[2];
			loc++;
			int enext(tmedge[tmesh[enow].edge[0]].face[0]);
			if (enext == enow)
				enext = tmedge[tmesh[enow].edge[0]].face[1];
			enow = enext;
			count++;
		}
	}
	GetSubMatrixCC_All(nv, tmesh[eid].amat, tmesh[eid].abar);
}

void TruncatedTspline::GetSubMatrixCC_All(int nv, MatrixXd &amat, MatrixXd &abar)
{
	//assume two irregular elements at least three rings!
	int nv2(2 * nv);
	int nbf(nv2 + 8), nbf1(nv2 + 17);
	amat = MatrixXd::Zero(nbf, nbf);
	abar = MatrixXd::Zero(nbf, nbf1);

	int i7(7);
	if (nv == 3)
		i7 = 1;
	int fcpt_id[8] = {6, nv2 + 8, 2, 4, nv2 + 10, nv2 + 16, nv2 + 14, nv2 + 12};
	int fcpt[8][4] = {{i7, 6, 5, 0}, {6, nv2 + 1, nv2 + 2, 5}, {1, 0, 3, 2}, {0, 5, 4, 3}, {5, nv2 + 2, nv2 + 3, 4}, {2, 3, nv2 + 6, nv2 + 7}, {3, 4, nv2 + 5, nv2 + 6}, {4, nv2 + 3, nv2 + 4, nv2 + 5}};
	double fcc(.25);
	int edpt_id[11] = {5, nv2 + 3, nv2 + 5, 3, 1, nv2 + 1, nv2 + 9, nv2 + 11, nv2 + 13, nv2 + 15, nv2 + 7};
	int edpt[11][6] = {{0, 5, i7, 6, 3, 4}, {5, 4, 0, 3, nv2 + 2, nv2 + 3}, {3, 4, 0, 5, nv2 + 6, nv2 + 5}, {0, 3, 1, 2, 5, 4}, {0, 1, nv2, nv2 - 1, 2, 3}, {5, 6, 0, i7, nv2 + 2, nv2 + 1}, {5, nv2 + 2, 6, nv2 + 1, 4, nv2 + 3}, {4, nv2 + 3, 5, nv2 + 2, nv2 + 5, nv2 + 4}, {4, nv2 + 5, 3, nv2 + 6, nv2 + 3, nv2 + 4}, {3, nv2 + 6, 2, nv2 + 7, 4, nv2 + 5}, {2, 3, 1, 0, nv2 + 7, nv2 + 6}};
	double edc[6] = {3. / 8., 3. / 8., 1. / 16., 1. / 16., 1. / 16., 1. / 16.};
	int vtpt_id[3] = {nv2 + 2, nv2 + 4, nv2 + 6};
	int vtpt[3][9] = {{5, 6, nv2 + 2, 4, 0, i7, nv2 + 1, nv2 + 3, 3}, {4, 5, nv2 + 3, nv2 + 5, 3, 0, nv2 + 2, nv2 + 4, nv2 + 6}, {3, 0, 4, nv2 + 6, 2, 1, 5, nv2 + 5, nv2 + 7}};
	double vtc[9] = {9. / 16., 3. / 32., 3. / 32., 3. / 32., 3. / 32., 1. / 64., 1. / 64., 1. / 64., 1. / 64.};
	//regular region
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			abar(fcpt[i][j], fcpt_id[i]) = fcc;
		}
		if (fcpt_id[i] < nbf)
		{
			for (int j = 0; j < 4; j++)
			{
				amat(fcpt[i][j], fcpt_id[i]) = fcc;
			}
		}
	}
	for (int i = 0; i < 11; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			abar(edpt[i][j], edpt_id[i]) = edc[j];
		}
		if (edpt_id[i] < nbf)
		{
			for (int j = 0; j < 6; j++)
			{
				amat(edpt[i][j], edpt_id[i]) = edc[j];
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 9; j++)
		{
			abar(vtpt[i][j], vtpt_id[i]) = vtc[j];
		}
		if (vtpt_id[i] < nbf)
		{
			for (int j = 0; j < 9; j++)
			{
				amat(vtpt[i][j], vtpt_id[i]) = vtc[j];
			}
		}
	}
	//irregular region
	for (int i = 4; i <= nv; i++) //face points
	{
		int xfcpt_id(2 * i);
		int xfcpt[4] = {0, 2 * i - 1, 2 * i, 2 * i + 1};
		if (i == nv) //the last
		{
			xfcpt[3] = 1;
		}
		for (int j = 0; j < 4; j++)
		{
			abar(xfcpt[j], xfcpt_id) = fcc;
			amat(xfcpt[j], xfcpt_id) = fcc;
		}
	}
	for (int i = 4; i <= nv; i++) //edge points
	{
		int xedpt_id(2 * i - 1);
		int xedpt[6] = {0, 2 * i - 1, 2 * i + 1, 2 * i, 2 * i - 3, 2 * i - 2};
		if (i == nv)
		{
			xedpt[2] = 1;
		}
		for (int j = 0; j < 6; j++)
		{
			abar(xedpt[j], xedpt_id) = edc[j];
			amat(xedpt[j], xedpt_id) = edc[j];
		}
	}
	double xvtc[3] = {1. - 7. / (4. * double(nv)), 3. / (2. * double(nv * nv)), 1. / (4. * double(nv * nv))};
	abar(0, 0) = xvtc[0];
	amat(0, 0) = xvtc[0];
	for (int i = 1; i <= nv; i++)
	{
		abar(2 * i - 1, 0) = xvtc[1];
		abar(2 * i, 0) = xvtc[2];
		amat(2 * i - 1, 0) = xvtc[1];
		amat(2 * i, 0) = xvtc[2];
	}
}

void TruncatedTspline::UpdateIENCC()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		vector<int>().swap(tmesh[eid].IEN);
		vector<array<double, 5>>().swap(tmesh[eid].patch_ku);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kv);
		tmesh[eid].IEN = tmesh[eid].IENtmp;
		tmesh[eid].patch_ku = tmesh[eid].patch_kutmp;
		tmesh[eid].patch_kv = tmesh[eid].patch_kvtmp;
		vector<int>().swap(tmesh[eid].IENtmp);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kutmp);
		vector<array<double, 5>>().swap(tmesh[eid].patch_kvtmp);
	}
}

void TruncatedTspline::VisualizeGeomCC(string fn)
{
	vector<array<double, 3>> spt;
	vector<array<double, 3>> sval;
	vector<array<int, 4>> sele;
	vector<array<double, 3>> lpt; //visulize parameter lines
	vector<array<int, 2>> led;	  //line connectivity
	int ns(2), ecount(0), loc0, loc1, loc2;
	vector<double> su(ns), sv(ns);
	for (int i = 0; i < ns; i++)
	{
		su[i] = i * 1. / (ns - 1);
		sv[i] = i * 1. / (ns - 1);
	}

	for (uint e = 0; e < tmesh.size(); e++)
	{
		if (tmesh[e].act == 1 && (tmesh[e].type == 0 || tmesh[e].type == 4))
		{
			int loc(0);
			for (int a = 0; a < ns; a++)
			{
				for (int b = 0; b < ns; b++)
				{
					array<double, 3> pt;
					array<double, 3> nm;
					SurfacePointMap(e, su[b], sv[a], pt, nm);
					spt.push_back(pt);
					sval.push_back(nm);
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
					el[0] = ecount * ns * ns + a * ns + b;
					el[1] = ecount * ns * ns + a * ns + b + 1;
					el[2] = ecount * ns * ns + (a + 1) * ns + b + 1;
					el[3] = ecount * ns * ns + (a + 1) * ns + b;
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
	}

	string fname = fn + "_geom.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i < spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
		}
		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i]<<"\n";
		//}
		//fout << "\nPOINT_DATA " << sval.size() << "\nNORMALS Normal FLOAT\n";
		//for (uint i = 0; i<sval.size(); i++)
		//{
		//	fout << sval[i][0] << " " << sval[i][1] << " " << sval[i][2] << "\n";
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
}

void TruncatedTspline::ElementBasisCC(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	if (tmesh[eid].act == 1)
	{
		if (tmesh[eid].type == 0)
		{
			ElementBasis_Regular(eid, u, v, Nt, dNdt);
		}
		else if (tmesh[eid].type == 4)
		{
			ElementBasis_IrregularCC(eid, u, v, Nt, dNdt);
		}
	}
}

void TruncatedTspline::ElementBasis_IrregularCC(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	//need to change variabe for Bezier, (u,v) -> (u_b,v_b)
	Nt.clear();
	dNdt.clear();
	Nt.resize(tmesh[eid].IEN.size(), 0.);
	dNdt.resize(tmesh[eid].IEN.size());

	double tol(1.e-3), eps(1.e-8);
	double u1 = max(tol, u);
	double v1 = max(tol, v);
	int n = floor(min(-log2(u1), -log2(v1))) + 1, subid;
	double pow2 = pow(2., double(n - 1));
	u1 *= pow2;
	v1 *= pow2;
	if (v1 < 0.5)
	{
		subid = 0;
		u1 = 2. * u1 - 1.;
		v1 = 2. * v1;
	}
	else if (u1 < 0.5)
	{
		subid = 2;
		u1 = 2. * u1;
		v1 = 2. * v1 - 1.;
	}
	else
	{
		subid = 1;
		u1 = 2. * u1 - 1.;
		v1 = 2. * v1 - 1.;
	}
	//cout << "subid: " << subid << "\n";
	//cout << "uv: " << u1 << " " << v1 << "\n";
	int nv((tmesh[eid].IEN.size() - 8) / 2);
	double Nu[4] = {(1. - 3. * u1 + 3. * u1 * u1 - u1 * u1 * u1) / 6., (4. - 6. * u1 * u1 + 3. * u1 * u1 * u1) / 6.,
					(1. + 3. * u1 + 3. * u1 * u1 - 3. * u1 * u1 * u1) / 6., u1 * u1 * u1 / 6.};
	double Nv[4] = {(1. - 3. * v1 + 3. * v1 * v1 - v1 * v1 * v1) / 6., (4. - 6. * v1 * v1 + 3. * v1 * v1 * v1) / 6.,
					(1. + 3. * v1 + 3. * v1 * v1 - 3. * v1 * v1 * v1) / 6., v1 * v1 * v1 / 6.};
	double dNdu[4] = {0., 0., 0., 0.};
	double dNdv[4] = {0., 0., 0., 0.};
	double Nt1[16];
	double dNdt1[16][2];
	int loc(0);
	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			Nt1[loc] = Nu[i] * Nv[j];
			dNdt1[loc][0] = dNdu[i] * Nv[j];
			dNdt1[loc][1] = Nu[i] * dNdv[j];
			loc++;
		}
	}
	int i7(7), nv2(2 * nv);
	if (nv == 3)
		i7 = 1;
	int Pk[3][16] = {{i7, 6, nv2 + 1, nv2 + 8, 0, 5, nv2 + 2, nv2 + 9, 3, 4, nv2 + 3, nv2 + 10, nv2 + 6, nv2 + 5, nv2 + 4, nv2 + 11},
					 {0, 5, nv2 + 2, nv2 + 9, 3, 4, nv2 + 3, nv2 + 10, nv2 + 6, nv2 + 5, nv2 + 4, nv2 + 11, nv2 + 15, nv2 + 14, nv2 + 13, nv2 + 12},
					 {1, 0, 5, nv2 + 2, 2, 3, 4, nv2 + 3, nv2 + 7, nv2 + 6, nv2 + 5, nv2 + 4, nv2 + 16, nv2 + 15, nv2 + 14, nv2 + 13}};
	//cout << "sum: " << sum << "\n";
	MatrixXd mat(nv2 + 8, 16);
	for (int i = 0; i < nv2 + 8; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			mat(i, j) = tmesh[eid].abar(i, Pk[subid][j]);
		}
	}
	if (n > 10)
	{
		cout << "n: " << n << "\n";
	}
	for (int i = 0; i < n - 1; i++)
	{
		mat = tmesh[eid].amat * mat;
	}
	for (int i = 0; i < nv2 + 8; i++)
	{
		Nt[i] = 0.;
		dNdt[i][0] = 0.;
		dNdt[i][1] = 0.;
		for (int j = 0; j < 16; j++)
		{
			Nt[i] += mat(i, j) * Nt1[j];
			dNdt[i][0] += mat(i, j) * dNdt1[j][0];
			dNdt[i][1] += mat(i, j) * dNdt1[j][1];
		}
	}
}

void TruncatedTspline::GlobalRefineCC()
{
	nel_old = tmesh.size();
	InitializeRefine();
	for (uint eid = 0; eid < nel_old; eid++)
	{
		if (tmesh[eid].act == 1)
		{
			if (tmesh[eid].type == 0)
			{
				ElementRefineReg_Dual(eid);
			}
			else if (tmesh[eid].type == 4)
			{
				ElementRefineIrr_CC(eid);
			}
		}
	}
	//additional refinement of type 2 boundary element
	for (uint eid = 0; eid < nel_old; eid++)
	{
		if (tmesh[eid].act == 1 && tmesh[eid].type == 2)
		{
			ElementRefineBnd_Dual(eid);
		}
	}

	UpdateCP();
	BuildCC();
}

void TruncatedTspline::ElementRefineIrr_CC(int eid)
{
	int nv((tmesh[eid].IEN.size() - 8) / 2);
	vector<Vertex> pnew(tmesh[eid].IEN.size());
	vector<int> pid(tmesh[eid].IEN.size());
	for (uint i = 0; i < pnew.size(); i++)
	{
		pnew[i].update = 1;
		for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
		{
			double ctmp(tmesh[eid].amat(j, i));
			if (ctmp > 0.)
			{
				pnew[i].coortmp[0] += ctmp * cp[tmesh[eid].IEN[j]].coor[0];
				pnew[i].coortmp[1] += ctmp * cp[tmesh[eid].IEN[j]].coor[1];
				pnew[i].coortmp[2] += ctmp * cp[tmesh[eid].IEN[j]].coor[2];
			}
		}
	}
	//cout << pnew[0].coortmp[0] << " " << pnew[0].coortmp[1] << "\n";
	//getchar();
	int nv2(2 * nv);
	int vtloc[4] = {0, nv2 + 2, nv2 + 4, nv2 + 6};
	int edloc[4] = {5, nv2 + 3, nv2 + 5, 3};
	int fcloc(4);
	int edbloc[4][3] = {{nv - 1, nv, nv + 5}, {nv + 13, nv + 14, nv + 15}, {nv + 17, nv + 18, nv + 19}, {nv + 11, nv + 4, 1}}; //edbloc[1] and edbloc[2] not useful
	cp.push_back(pnew[fcloc]);
	pid[fcloc] = cp.size() - 1;
	for (int i = 0; i < 4; i++)
	{
		if (cp[tmesh[eid].cnct[i]].update == 0)
		{
			cp[tmesh[eid].cnct[i]].coortmp[0] = pnew[vtloc[i]].coortmp[0];
			cp[tmesh[eid].cnct[i]].coortmp[1] = pnew[vtloc[i]].coortmp[1];
			cp[tmesh[eid].cnct[i]].coortmp[2] = pnew[vtloc[i]].coortmp[2];
			cp[tmesh[eid].cnct[i]].update = 1;
		}
		pid[vtloc[i]] = tmesh[eid].cnct[i];
		//if (tmedge[tmesh[eid].edge[i]].act == 1)
		if (tmedge[tmesh[eid].edge[i]].midpt == -1)
		{
			cp.push_back(pnew[edloc[i]]);
			tmedge[tmesh[eid].edge[i]].midpt = cp.size() - 1;
			//tmedge[tmesh[eid].edge[i]].act = 0;
		}
		pid[edloc[i]] = tmedge[tmesh[eid].edge[i]].midpt;
	}
	//build new elements and edges
	int fccnct[4][4] = {{0, 5, 4, 3}, {5, nv2 + 2, nv2 + 3, 4}, {4, nv2 + 3, nv2 + 4, nv2 + 5}, {3, 4, nv2 + 5, nv2 + 6}};
	int edcnct[12][2] = {{0, 5}, {5, nv2 + 2}, {4, 5}, {nv2 + 2, nv2 + 3}, {nv2 + 3, nv2 + 4}, {4, nv2 + 3}, {nv2 + 4, nv2 + 5}, {nv2 + 5, nv2 + 6}, {4, nv2 + 5}, {nv2 + 6, 3}, {3, 0}, {4, 3}};
	vector<Edge> ednew(12);
	int edid[12];
	for (int i = 0; i < 4; i++)
	{
		int iloc[3] = {3 * i, 3 * i + 1, 3 * i + 2};
		for (int j = 0; j < 3; j++)
		{
			ednew[iloc[j]].pt[0] = pid[edcnct[iloc[j]][0]];
			ednew[iloc[j]].pt[1] = pid[edcnct[iloc[j]][1]];
		}
		//if (tmedge[tmesh[eid].edge[i]].act == 1)
		if (tmedge[tmesh[eid].edge[i]].chd[0] == -1)
		{
			tmedge.push_back(ednew[iloc[0]]);
			tmedge[tmesh[eid].edge[i]].chd[0] = tmedge.size() - 1;
			tmedge.push_back(ednew[iloc[1]]);
			tmedge[tmesh[eid].edge[i]].chd[1] = tmedge.size() - 1;
			tmedge[tmesh[eid].edge[i]].act = 0;
		}
		edid[iloc[0]] = tmedge[tmesh[eid].edge[i]].chd[0];
		edid[iloc[1]] = tmedge[tmesh[eid].edge[i]].chd[1];
		if (tmedge[edid[iloc[0]]].pt[0] != tmesh[eid].cnct[i] && tmedge[edid[iloc[0]]].pt[1] != tmesh[eid].cnct[i])
		{
			edid[iloc[0]] = tmedge[tmesh[eid].edge[i]].chd[1];
			edid[iloc[1]] = tmedge[tmesh[eid].edge[i]].chd[0];
		}
		tmedge.push_back(ednew[iloc[2]]);
		edid[iloc[2]] = tmedge.size() - 1;
	}
	int fced[4][4] = {{0, 2, 11, 10}, {1, 3, 5, 2}, {5, 4, 6, 8}, {11, 8, 7, 9}};
	vector<Element> enew(4);
	for (int i = 0; i < 4; i++)
	{
		enew[i].act = 1;
		if (i == 0)
			enew[i].type = 4;
		for (int j = 0; j < 4; j++)
		{
			enew[i].cnct[j] = pid[fccnct[i][j]];
			enew[i].edge[j] = edid[fced[i][j]];
		}
		tmesh.push_back(enew[i]);
		tmesh[eid].chd[i] = tmesh.size() - 1;
	}

	tmesh[eid].act = 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////

void TruncatedTspline::Run_UTSP(string fn_in, string fn_out)
{
	InputFromVTK_Quad(fn_in);
	//RescaleCoor();
	InitialConnect_UT();

	//int rf_id0[9] = { 28,64,100,34,29,65,70,101,106 };
	//int rf_type0[9]={0,0,0,3,3,3,3,3,3};
	//vector<int> rf_id(rf_id0,rf_id0+9);
	//vector<int> rf_type(rf_type0,rf_type0+9);
	//Refine_Unstruct(rf_id,rf_type);

	//SetDPatch_Irr();//design space
	//CollectActives();

	InitializeProjectFacePoints(SMOOTH_TYPE, SMOOTH_BETA); //analysis space, initial step

	int nref(0);
	for (int i = 0; i < nref; i++)
	{
		vector<int> rftype, rfid, c1id;
		FindC1Element(c1id);
		FindGlobalRefineID(rfid, rftype);
		Refine_Unstruct(rfid, rftype);
		Refine_C1Element(c1id);
		if (i == nref - 1)
		{
			AddProjectFacePoints(SMOOTH_TYPE, SMOOTH_BETA);
			SetTrunMat();
		}
	}

	//AddProjectFacePoints(SMOOTH_TYPE, SMOOTH_BETA);
	//SetTrunMat();
	CollectActives_DPatchAnalysis();

	//int nfit(1);
	//for (int i = 0; i < nfit; i++)
	//{
	//	vector<BezierElement> bzmesh;
	//	BezierExtract_DPatch(bzmesh);
	//	//OutputBezierMesh(bzmesh, "../io/surf2/hemisphere0"); getchar();
	//	LeastSquare ls;
	//	vector<array<double, 3>> cpnew;
	//	ls.SetProblem(cp.size());
	//	ls.Run(bzmesh, "", cpnew);
	//	UpdateCP_Fitting(cpnew);
	//}

	vector<BezierElement> bzmesh;
	////BezierExtract_DPatch(bzmesh);
	BezierExtract_DPatch_AS(bzmesh);
	////ReadBezierMat("../io/hemisphere/h0_bezex_deepesh.txt",bzmesh);//to compare with Deepesh's matrix
	OutputBezierMesh(bzmesh, fn_out);

	//VisualizeSurface(fn_out);
	VisualizeQuad(fn_out);

	//WriteBezier_Shell(bzmesh, fn_out);

	//OutputBezierPoint22("../io/surf2/square_XP2_1_bzpt");
}

void TruncatedTspline::Run_UTSP_1(string fn_in, string fn_out, vector<BezierElement> &bzmesh, vector<int> &IDBC, vector<double> &gh)
{
	InputFromVTK_Quad(fn_in);
	//RescaleCoor();
	InitialConnect_UT();

	//VisualizeEdge("../io/centralpillar/cenpil"); cout << "done\n"; getchar();

	//SetDPatch_Irr();//design space
	//CollectActives();

	InitializeProjectFacePoints(SMOOTH_TYPE, SMOOTH_BETA); //analysis space, initial step

	int nref(0);
	for (int i = 0; i < nref; i++)
	{
		vector<int> rftype, rfid, c1id;
		FindC1Element(c1id);
		FindGlobalRefineID(rfid, rftype);
		Refine_Unstruct(rfid, rftype);
		Refine_C1Element(c1id);
		//EnlargeOneRing();
		//if (i == nref - 1)
		//{
		//	AddProjectFacePoints(SMOOTH_TYPE, SMOOTH_BETA);
		//	SetTrunMat();
		//}
	}
	//AddProjectFacePoints(SMOOTH_TYPE, SMOOTH_BETA);
	//SetTrunMat();

	//HemisphereLocalRefine1_6();//local refinement based on the mesh after gloabl refinement once, input hemisphere_2

	//HemisphereLocalRefine2_6();//local refinement based on the mesh after gloabl refinement twice, input hemishpere_2

	//HemisphereLocalRefine("../io/hemisphere_adp4/lev", 3);//local refinement
	//HemisphereLocalRefine("../io/cylinder_reg/lev", 4);//local refinement

	//AddProjectFacePoints(SMOOTH_TYPE, SMOOTH_BETA);
	//SetTrunMat();
	CollectActives_DPatchAnalysis(); //set paid

	//int nfit(1);
	//for (int i = 0; i < nfit; i++)//design sapce
	//{
	//	vector<BezierElement> bzmesh;
	//	BezierExtract_DPatch(bzmesh);
	//	//OutputBezierMesh(bzmesh, "../io/surf2/hemisphere0"); getchar();
	//	LeastSquare ls;
	//	vector<array<double, 3>> cpnew;
	//	ls.SetProblem(cp.size());
	//	ls.Run(bzmesh, "", cpnew);
	//	UpdateCP_Fitting(cpnew);
	//}

	//int nfit(1);
	//for (int i = 0; i < nfit; i++)//analysis sapce
	//{
	//	vector<BezierElement> bzmesh1, bzmesh2;
	//	BezierExtract_DPatch_AS(bzmesh1);
	//	BezierMeshUpdateIndex(bzmesh1, bzmesh2);//need to use paid
	//	//OutputBezierMesh(bzmesh, "../io/surf2/hemisphere0"); getchar();
	//	vector<int> IDBC1;
	//	vector<double> gh1;
	//	BC_Hemisphere(IDBC1, gh1);
	//	LeastSquare ls;
	//	vector<array<double, 3>> cpnew;
	//	ls.SetProblem(IDBC1, gh1);
	//	ls.Run(bzmesh2, "", cpnew);
	//	UpdateCP_Fitting(cpnew);
	//}

	vector<BezierElement> bzmesh0;
	////BezierExtract_DPatch(bzmesh);
	BezierExtract_DPatch_AS(bzmesh0);
	////ReadBezierMat("../io/hemisphere/h0_bezex_deepesh.txt",bzmesh);//to compare with Deepesh's matrix
	OutputBezierMesh(bzmesh0, fn_out);

	//BezierMeshUpdateIndex(bzmesh0, bzmesh);
	//ReadBezier("../io/hemisphere/hemisphere_analysis0_bezex.txt", bzmesh);
	//BC_Square(IDBC, gh);

	//string fn_out("../io/UTSP_test0/square2");
	VisualizeSurface(fn_out);
	//VisualizeQuad(fn_out);
	VisualizeCM_DPatch(fn_out);
	VisualizeControlMesh(fn_out);

	WriteBezier_Shell(bzmesh0, fn_out);
	//WriteBezier_LSDYNA_C12(fn_out, bzmesh0);

	//OutputBezierPoint22("../io/surf2/square_XP2_1_bzpt");
}

void TruncatedTspline::Run_UTSP_AS(string fn_in, string fn_out, int ngref, int nlref, int fit_flag, string fn_lref, vector<BezierElement> &bzmesh, vector<int> &IDBC, vector<double> &gh)
{
	InputFromVTK_Quad(fn_in);
	RescaleCoor();
	InitialConnect_UT();

	//HemisphereLocalRefine(fn_lref, nlref);//local refinement

	//D-patch first and then global refinement
	InitializeProjectFacePoints(SMOOTH_TYPE, SMOOTH_BETA); //analysis space, initial step

	int nref(ngref); //global refinement
	for (int i = 0; i < nref; i++)
	{
		vector<int> rftype, rfid, c1id;
		FindC1Element(c1id);
		FindGlobalRefineID(rfid, rftype);
		Refine_Unstruct(rfid, rftype);
		Refine_C1Element(c1id);
		if (i == nref - 1)
		{
			AddProjectFacePoints(SMOOTH_TYPE, SMOOTH_BETA);
			SetTrunMat();
		}
	}

	HemisphereLocalRefine(fn_lref, nlref);
	////global refinement first and then D-patch
	//int nref(1);//global refinement
	//for (int i = 0; i < nref; i++)
	//{
	//	vector<int> rftype, rfid, c1id;
	//	FindGlobalRefineID(rfid, rftype);
	//	Refine_Unstruct(rfid, rftype);
	//	if (i == nref - 1)
	//	{
	//		InitializeProjectFacePoints(SMOOTH_TYPE, SMOOTH_BETA);
	//	}
	//}

	CollectActives_DPatchAnalysis(); //set paid

	int nfit(fit_flag);
	for (int i = 0; i < nfit; i++) //analysis sapce
	{
		vector<BezierElement> bzmesh1, bzmesh2;
		BezierExtract_DPatch_AS(bzmesh1);
		BezierMeshUpdateIndex(bzmesh1, bzmesh2); //need to use paid
		//OutputBezierMesh(bzmesh, "../io/surf2/hemisphere0"); getchar();
		vector<int> IDBC1;
		vector<double> gh1;
		BC_Hemisphere(IDBC1, gh1);
		LeastSquare ls;
		vector<array<double, 3>> cpnew;
		ls.SetProblem(IDBC1, gh1);
		ls.Run(bzmesh2, "", cpnew);
		UpdateCP_Fitting(cpnew);
	}

	//vector<BezierElement> bzmesh0;
	BezierExtract_DPatch_AS(bzmesh);
	OutputBezierMesh(bzmesh, fn_out);

	//VisualizeSurface(fn_out);
	VisualizeCM_DPatch(fn_out);
	VisualizeControlMesh(fn_out, 0);

	//WriteBezier_Shell(bzmesh, fn_out);
	WriteBezier_LSDYNA_C12(fn_out, bzmesh);
	//WriteBezier_Angran(bzmesh, fn_out);
}

void TruncatedTspline::Run_UTSP_Design(string fn_in, string fn_out)
{
	cout << "Constructing design space...\n"; //EP not allowed in 3-ring of boundary. T-junction not allowed in 4-ring of an EP. Consistent with previous version

	bool Rectangle2UnitSquare(false); //rescale a rectangle domain to a unit square domain if needed
	bool LocalRefine(false);		  //local refinement performed on the input mesh
	bool GeometricFitting(false);

	bool OutputControlMesh(false);	///Angran add, to output Control Mesh after geometric fitting
	bool OutputShellHugo(true);		//output shell Bezier elements in Hugo's format
	bool OutputShellLSDYNA(false);	//output shell Bezier elements in LS-DYNA's format
	bool OutputGeometry(false);		//output geometric mapping result in the physical domain, not available
	bool BezierPointsCompare(true); //used to check if the geometry is the same before and after refinement, by checking Bezier points.

	InputFromVTK_Quad(fn_in);
	if (Rectangle2UnitSquare)
		RescaleCoor();
	InitialConnect_UT();

	if (LocalRefine)
	{
		vector<int> rfid(1, 89);  //element id to be refined
		vector<int> rftype(1, 1); //refinement type
		Refine_Unstruct(rfid, rftype);
	}

	SetDPatch_Irr(); //design space
	CollectActives();

	//int nfit(1);
	if (GeometricFitting)
	{
		vector<BezierElement> bzmesh;
		BezierExtract_DPatch(bzmesh);
		//OutputBezierMesh(bzmesh, "../io/surf2/hemisphere0"); getchar();
		LeastSquare ls;
		vector<array<double, 3>> cpnew;
		ls.SetProblem(cp.size());
		ls.Run(bzmesh, "", cpnew);
		UpdateCP_Fitting(cpnew);
	}

	vector<BezierElement> bzmesh;
	BezierExtract_DPatch(bzmesh);
	OutputBezierMesh(bzmesh, fn_out);

	if (BezierPointsCompare)
	{
		OutputBezierMeshAll(bzmesh, fn_out);
	}

	//	if(OutputGeometry)
	//	{
	//		VisualizeSurface(fn_out);
	//	}

	if (OutputControlMesh)
	{
		VisualizeControlMesh(fn_out);
	}

	if (OutputShellHugo)
	{
		WriteBezier_Shell(bzmesh, fn_out);
	}
}

void TruncatedTspline::Run_UTSP_AS_Scale(string fn_in, string fn_out, int ngref, int nlref, int fit_flag, string fn_lref, vector<BezierElement> &bzmesh,
										 vector<int> &IDBC, vector<double> &gh)
{
	cout << "Constructing analysis space with scaling to enforce partition of unit...\n";

	bool Rectangle2UnitSquare(false);	  //rescale a rectangle domain to a unit square domain if needed
	bool LocalRefineOnInputMesh(false);	  //local refinement performed on the input mesh
	bool LocalRefineOnRefinedMesh(false); //local refinement performed on the globally refined mesh
	//bool GeometricFitting(false);
	bool GeometricFitting(fit_flag);
	int numGlobalRefine(ngref);

	bool OutputShellHugo(true);	   //output shell Bezier elements in Hugo's format
	bool OutputShellLSDYNA(false); //output shell Bezier elements in LS-DYNA's format
	bool OutputGeometry(false);	   //output geometric mapping result in the physical domain
	bool OutputControlMesh(true);  //output control mesh to local element id for local refinement

	bool BezierPointsCompare(false); //used to check if the geometry is the same before and after refinement, by checking Bezier points.
	int CompareLevel[2] = {0, 1};	 //Compare Bezier points from level 0 with level 1, or any two consecutive levels, refined coarse Bezier points should match with fine Bezier points

	InputFromVTK_Quad(fn_in); //EP can't be a corner of an element in the 2-ring of boundary, but can be in the 3rd-ring, which is new!!
	if (Rectangle2UnitSquare)
		RescaleCoor();
	InitialConnect_UT();

	if (LocalRefineOnInputMesh)
	{
		//local refinement here first before calling InitializeProjectedPoints_Scale(), i.e., construct C1 elements
		//T-junctions not allowed in the 3-ring of an EP (including the boundary of 3-ring), but can be in the 4th-ring, which is new!!
		//vector<int> rfid(1,89); //rfid[1]=47; rfid[2]=52; rfid[3]=83;//element id to be refined
		//vector<int> rftype(1,1);//refinement type
		//Refine_Unstruct(rfid,rftype);
		//HemisphereLocalRefine("../io/cross_pipe_new/lev", nlref);
		HemisphereLocalRefine(fn_lref, nlref);
		//HemisphereLocalRefine("../io/hemisphere_24by24/lev", 3);//local refinement
		//HemisphereLocalRefine("../io/figure2/lev", 1);//local refinement
	}

	//D-patch first and then global refinement
	InitializeProjectFacePoints_Scale(SMOOTH_TYPE, SMOOTH_BETA); //analysis space, initial step
	if (!LocalRefineOnInputMesh)								 //global refinement is not available if a local refinement has been performed
	{
		for (int i = 0; i < numGlobalRefine; i++)
		{
			vector<int> rftype, rfid, c1id;
			FindC1Element(c1id);
			FindGlobalRefineID(rfid, rftype);
			Refine_Unstruct(rfid, rftype);
			Refine_C1Element(c1id);
			EnlargeOneRing();
			if (SMOOTH_TYPE == 1)
				SMOOTH_BETA *= 0.5;
		}
	}

	if (LocalRefineOnRefinedMesh) //local refinement after globally refined mesh
	{
		//HemisphereLocalRefine("../io/hemisphere_6by6_new/lev", 3);//local refinement
		//HemisphereLocalRefine("../io/cross_pipe_new/refine3_lev", 1);
		//HemisphereLocalRefine("../io/bpillar_paper/refine2_lev", 1);//local refinement
		HemisphereLocalRefine(fn_lref, nlref);
	}
	if (!LocalRefineOnInputMesh && numGlobalRefine > 0)
	{
		FindSplineC1_Scale();
		FindSplineC2_Scale();
	}
	CollectActives_DPatchAnalysis(); //set paid

	if (GeometricFitting) //geometric fitting when the geometry has an analytical form
	{
		vector<BezierElement> bzmesh1, bzmesh2;
		BezierExtract_DPatch_AS_Scale(bzmesh1);
		BezierMeshUpdateIndex(bzmesh1, bzmesh2); //need to use paid
		//OutputBezierMesh(bzmesh, "../io/surf2/hemisphere0"); getchar();
		vector<int> IDBC1;
		vector<double> gh1;
		BC_Hemisphere(IDBC1, gh1);
		LeastSquare ls;
		vector<array<double, 3>> cpnew;
		ls.SetProblem(IDBC1, gh1);
		ls.Run(bzmesh2, "", cpnew);
		UpdateCP_Fitting(cpnew);
	}

	if (BezierPointsCompare)
	{
		if (numGlobalRefine == CompareLevel[0])
		{
			BezierExtract_DPatch_AS_Scale_Refine(bzmesh);
			OutputBezierMeshAll(bzmesh, fn_out);
		}
		else if (numGlobalRefine == CompareLevel[1])
		{
			BezierExtract_DPatch_AS_Scale(bzmesh);
			OutputBezierMeshAll(bzmesh, fn_out);
		}
	}
	else
	{
		BezierExtract_DPatch_AS_Scale(bzmesh);
		OutputBezierMesh(bzmesh, fn_out);
		//		OutputBezierMeshAll(bzmesh, fn_out);
	}

	if (OutputGeometry)
	{
		VisualizeSurface_Scale(fn_out);
	}

	if (OutputShellHugo)
	{
		VisualizeCM_DPatch(fn_out);
		WriteBezier_Shell(bzmesh, fn_out);
	}

	if (OutputControlMesh)
	{
		VisualizeControlMesh(fn_out);
	}

	if (OutputShellLSDYNA)
	{
		WriteBezier_LSDYNA_C12(fn_out, bzmesh);
	}
}

void TruncatedTspline::Run_UTSP_AS_Scale_Honda(string fn_in, string fn_out, int ngref, int nlref, string fn_lref, vector<BezierElement> &bzmesh,
											   vector<int> &IDBC, vector<double> &gh)
{
	cout << "Constructing analysis space with scaling to enforce partition of unit...\n";

	bool Rectangle2UnitSquare(false);	  //rescale a rectangle domain to a unit square domain if needed
	bool LocalRefineOnInputMesh(true);	  //local refinement performed on the input mesh
	bool LocalRefineOnRefinedMesh(false); //local refinement performed on the globally refined mesh
	bool GeometricFitting(false);
	int numGlobalRefine(ngref);

	bool OutputShellHugo(false);   //output shell Bezier elements in Hugo's format
	bool OutputShellLSDYNA(true);  //output shell Bezier elements in LS-DYNA's format
	bool OutputGeometry(true);	   //output geometric mapping result in the physical domain
	bool OutputControlMesh(false); //output control mesh to local element id for local refinement

	bool BezierPointsCompare(false); //used to check if the geometry is the same before and after refinement, by checking Bezier points.
	int CompareLevel[2] = {0, 1};	 //Compare Bezier points from level 0 with level 1, or any two consecutive levels, refined coarse Bezier points should match with fine Bezier points

	InputFromVTK_Quad(fn_in); //EP can't be a corner of an element in the 2-ring of boundary, but can be in the 3rd-ring, which is new!!
	if (Rectangle2UnitSquare)
		RescaleCoor();
	InitialConnect_UT();

	if (LocalRefineOnInputMesh)
	{
		//local refinement here first before calling InitializeProjectedPoints_Scale(), i.e., construct C1 elements
		//T-junctions not allowed in the 3-ring of an EP (including the boundary of 3-ring), but can be in the 4th-ring, which is new!!
		//vector<int> rfid(1,89); //rfid[1]=47; rfid[2]=52; rfid[3]=83;//element id to be refined
		//vector<int> rftype(1,1);//refinement type
		//Refine_Unstruct(rfid,rftype);
		//HemisphereLocalRefine("../io/cross_pipe_new/lev", nlref);
		//HemisphereLocalRefine(fn_lref, nlref);
		//HemisphereLocalRefine("../io/hemisphere_24by24/lev", 3);//local refinement
		//HemisphereLocalRefine("../io/figure2/lev", 1);//local refinement
		HemisphereLocalRefine_Honda(fn_lref, nlref);
	}

	//D-patch first and then global refinement
	InitializeProjectFacePoints_Scale(SMOOTH_TYPE, SMOOTH_BETA); //analysis space, initial step
	if (!LocalRefineOnInputMesh)								 //global refinement is not available if a local refinement has been performed
	{
		for (int i = 0; i < numGlobalRefine; i++)
		{
			vector<int> rftype, rfid, c1id;
			FindC1Element(c1id);
			FindGlobalRefineID(rfid, rftype);
			Refine_Unstruct(rfid, rftype);
			Refine_C1Element(c1id);
			EnlargeOneRing();
			if (SMOOTH_TYPE == 1)
				SMOOTH_BETA *= 0.5;
		}
	}

	if (LocalRefineOnRefinedMesh) //local refinement after globally refined mesh
	{
		//HemisphereLocalRefine("../io/hemisphere_6by6_new/lev", 3);//local refinement
		//HemisphereLocalRefine("../io/cross_pipe_new/refine3_lev", 1);
		//HemisphereLocalRefine("../io/bpillar_paper/refine2_lev", 1);//local refinement
		//HemisphereLocalRefine(fn_lref, nlref);
	}
	if (!LocalRefineOnInputMesh && numGlobalRefine > 0)
	{
		FindSplineC1_Scale();
		FindSplineC2_Scale();
	}
	CollectActives_DPatchAnalysis(); //set paid

	if (GeometricFitting) //geometric fitting when the geometry has an analytical form
	{
		vector<BezierElement> bzmesh1, bzmesh2;
		BezierExtract_DPatch_AS_Scale(bzmesh1);
		BezierMeshUpdateIndex(bzmesh1, bzmesh2); //need to use paid
		//OutputBezierMesh(bzmesh, "../io/surf2/hemisphere0"); getchar();
		vector<int> IDBC1;
		vector<double> gh1;
		BC_Hemisphere(IDBC1, gh1);
		LeastSquare ls;
		vector<array<double, 3>> cpnew;
		ls.SetProblem(IDBC1, gh1);
		ls.Run(bzmesh2, "", cpnew);
		UpdateCP_Fitting(cpnew);
	}

	if (BezierPointsCompare)
	{
		if (numGlobalRefine == CompareLevel[0])
		{
			BezierExtract_DPatch_AS_Scale_Refine(bzmesh);
			OutputBezierMeshAll(bzmesh, fn_out);
		}
		else if (numGlobalRefine == CompareLevel[1])
		{
			BezierExtract_DPatch_AS_Scale(bzmesh);
			OutputBezierMeshAll(bzmesh, fn_out);
		}
	}
	else
	{
		BezierExtract_DPatch_AS_Scale(bzmesh);
		//OutputBezierMesh(bzmesh, fn_out);
		//		OutputBezierMeshAll(bzmesh, fn_out);
	}

	if (OutputGeometry)
	{
		VisualizeSurface_Scale(fn_out);
	}

	if (OutputShellHugo)
	{
		VisualizeCM_DPatch(fn_out);
		WriteBezier_Shell(bzmesh, fn_out);
	}

	if (OutputControlMesh)
	{
		VisualizeControlMesh(fn_out);
	}

	if (OutputShellLSDYNA)
	{
		WriteBezier_LSDYNA_C12(fn_out, bzmesh);
	}
}

void TruncatedTspline::HemisphereLocalRefine_Honda(string fpre, int nlev)
{
	for (int it = 0; it < nlev; it++)
	{
		string fn4 = fpre + to_string(it + 1) + ".txt";
		string fn2 = fpre + to_string(it + 1) + "_bi.txt";
		vector<int> rf4, rf2, rf2_cn;
		ReadRef4ID_Honda(fn4, rf4);
		ReadRef2ID(fn2, rf2, rf2_cn);
		HemisphereLocalRefine(rf4, rf2, rf2_cn);
	}
}

void TruncatedTspline::ReadRef4ID_Honda(string fn, vector<int> &rfid)
{
	rfid.clear();
	ifstream fin;
	fin.open(fn);
	if (fin.is_open())
	{
		string stmp;
		vector<int> itmp_vec;
		int itmp;
		while (getline(fin, stmp))
		{
			stringstream ss(stmp);
			while (ss >> itmp)
				itmp_vec.push_back(itmp);
			itmp = itmp_vec.back();
			if (tmesh[itmp].act == 1 && tmesh[itmp].type == 0)
			{
				rfid.push_back(itmp);
			}
		}
		fin.close();
	}
	else
	{
		cerr << "Can't open " << fn << "!\n";
	}
}

void TruncatedTspline::Run_UTSP_AS_Angran(string fn_in, string fn_out, int ngref, int nlref, int fit_flag, string fn_lref, vector<BezierElement> &bzmesh, vector<int> &IDBC, vector<double> &gh)
{
	//InputFromVTK_Quad(fn_in);
	InputFromVTK_Quad_Angran(fn_in);
	RescaleCoor();
	InitialConnect_UT();

	//HemisphereLocalRefine(fn_lref, nlref);//local refinement

	//D-patch first and then global refinement
	InitializeProjectFacePoints(SMOOTH_TYPE, SMOOTH_BETA); //analysis space, initial step

	int nref(ngref); //global refinement
	for (int i = 0; i < nref; i++)
	{
		vector<int> rftype, rfid, c1id;
		FindC1Element(c1id);
		FindGlobalRefineID(rfid, rftype);
		Refine_Unstruct(rfid, rftype);
		Refine_C1Element(c1id);
		if (i == nref - 1)
		{
			AddProjectFacePoints(SMOOTH_TYPE, SMOOTH_BETA);
			SetTrunMat();
		}
	}

	HemisphereLocalRefine(fn_lref, nlref);
	////global refinement first and then D-patch
	//int nref(1);//global refinement
	//for (int i = 0; i < nref; i++)
	//{
	//	vector<int> rftype, rfid, c1id;
	//	FindGlobalRefineID(rfid, rftype);
	//	Refine_Unstruct(rfid, rftype);
	//	if (i == nref - 1)
	//	{
	//		InitializeProjectFacePoints(SMOOTH_TYPE, SMOOTH_BETA);
	//	}
	//}

	CollectActives_DPatchAnalysis(); //set paid

	int nfit(fit_flag);
	for (int i = 0; i < nfit; i++) //analysis sapce
	{
		vector<BezierElement> bzmesh1, bzmesh2;
		BezierExtract_DPatch_AS(bzmesh1);
		BezierMeshUpdateIndex(bzmesh1, bzmesh2); //need to use paid
		//OutputBezierMesh(bzmesh, "../io/surf2/hemisphere0"); getchar();
		vector<int> IDBC1;
		vector<double> gh1;
		BC_Hemisphere(IDBC1, gh1);
		LeastSquare ls;
		vector<array<double, 3>> cpnew;
		ls.SetProblem(IDBC1, gh1);
		ls.Run(bzmesh2, "", cpnew);
		UpdateCP_Fitting(cpnew);
	}

	//vector<BezierElement> bzmesh0;
	BezierExtract_DPatch_AS(bzmesh);
	//OutputBezierMesh(bzmesh, fn_out);

	//VisualizeSurface(fn_out);
	VisualizeCM_DPatch(fn_out);
	//VisualizeControlMesh(fn_out, 0);

	//WriteBezier_Shell(bzmesh, fn_out);
	//WriteBezier_LSDYNA_C12(fn_out, bzmesh);
	//WriteBezier_Angran(bzmesh, fn_out);

	WriteBezier_Angran(bzmesh, fn_out);
}

void TruncatedTspline::InputFromVTK_Quad(string fn)
{
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++)
			getline(fin, stmp); //skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
			//cp[i].act = 1;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			//tmesh[i].act = 1;
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline::InputFromVTK_Quad_Angran(string fn)
{
	string fname(fn + "/controlmesh.vtk"), stmp, stmp_line;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		while (getline(fin, stmp))
		{
			if (fin)
			{
				string upper(stmp);

				std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
				if (upper.find("POINTS") != string::npos)
				{
					stringstream sstmp(upper);
					sstmp >> stmp_line >> npts >> stmp_line;
					cout << "Start reading node # of points: " << npts << endl;
					cp.resize(npts);
					for (int i = 0; i < npts; i++)
					{
						fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
					}
					//cout << "Node number: " << cp.size() << endl;
				}
				else if (upper.find("CELLS") != string::npos)
				{
					stringstream sstmp(upper);
					sstmp >> stmp_line >> neles >> stmp_line;
					cout << "Start reading cells # of elements: " << neles << endl;
					tmesh.resize(neles);
					for (int i = 0; i < neles; i++)
					{
						fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3];
					}
					//cout << "Node number: " << cp.size() << endl;
				}
				else if (upper.find("POINT_DATA") != string::npos)
				{
					stringstream sstmp(upper);
					sstmp >> stmp_line >> npts;
					cout << "Start reading point labels # of labels: " << npts << endl;
					for (int i = 0; i < 2; i++)
					{
						getline(fin, stmp); //skip lines
					}
					for (int i = 0; i < npts; i++)
					{
						fin >> cp[i].label;
						//cout <<  cp[i].label << endl;
					}
					//cout << "Node number: " << cp.size() << endl;
				}
				else if (upper.find("VECTORS") != string::npos)
				{
					cout << "Start reading initial velocity" << endl;
					for (int i = 0; i < cp.size(); i++)
					{
						fin >> cp[i].velocity[0] >> cp[i].velocity[1] >> cp[i].velocity[2];
						//cout <<  cp[i].label << endl;
					}
				}
				// 2.) Look for the "*Element" section
				continue;
			}
			//if (fin.eof())
			//{
			//	cout << "Read to the end!" << endl;
			//	break;
			//}
		}
		fin.close();
		cout << "File is closed!" << endl;
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline::RescaleCoor()
{
	double xymax[2] = {-1.e5, -1.e5};
	double xymin[2] = {1.e5, 1.e5};
	for (uint i = 0; i < cp.size(); i++)
	{
		xymax[0] = cp[i].coor[0] > xymax[0] ? cp[i].coor[0] : xymax[0];
		xymax[1] = cp[i].coor[1] > xymax[1] ? cp[i].coor[1] : xymax[1];
		xymin[0] = cp[i].coor[0] < xymin[0] ? cp[i].coor[0] : xymin[0];
		xymin[1] = cp[i].coor[1] < xymin[1] ? cp[i].coor[1] : xymin[1];
	}
	double dim[2] = {xymax[0] - xymin[0], xymax[1] - xymin[1]};
	cout << dim[0] << " " << dim[1];
	//getchar();
	//for (uint i = 0; i < cp.size(); i++)
	//{
	//	cp[i].coor[0] /= dim[0];
	//	cp[i].coor[1] /= dim[1];
	//	//cp[i].coor[2] = 10.*cp[i].coor[0]* cp[i].coor[1]*(1. - cp[i].coor[0])*(1. - cp[i].coor[1]);
	//	//cp[i].coor[2] = (1. - (cp[i].coor[0] - xymin[0]) / dim[0])*(1. - (cp[i].coor[1] - xymin[1]) / dim[1]);
	//}
}

void TruncatedTspline::InitialConnect_UT()
{
	uint i, j, k;
	for (i = 0; i < tmesh.size(); i++)
	{
		tmesh[i].act = 1;
		for (j = 0; j < 4; j++)
		{
			cp[tmesh[i].cnct[j]].face.push_back(i);
		}
	}
	for (i = 0; i < cp.size(); i++)
	{
		cp[i].act = 1;
		if (cp[i].face.size() == 3 || cp[i].face.size() > 4) //boundary extraordinary points not considered
		{
			cp[i].type = 2;
		}
	}
	for (i = 0; i < tmesh.size(); i++)
	{
		int pos(-1), count(0);
		for (j = 0; j < 4; j++)
		{
			if (cp[tmesh[i].cnct[j]].type == 2)
			{
				pos = j;
				count++;
			}
		}
		if (count != 0)
		{
			if (count == 1)
				tmesh[i].type = 4;
			else
				tmesh[i].type = 5;
			int cnctnew[4] = {tmesh[i].cnct[pos], tmesh[i].cnct[(pos + 1) % 4], tmesh[i].cnct[(pos + 2) % 4], tmesh[i].cnct[(pos + 3) % 4]};
			for (j = 0; j < 4; j++)
			{
				tmesh[i].cnct[j] = cnctnew[j];
			}
		}
	}
	for (i = 0; i < tmesh.size(); i++)
	{
		for (j = 0; j < 4; j++)
		{
			Edge edtmp;
			edtmp.act = 1;
			edtmp.pt[0] = tmesh[i].cnct[j];
			edtmp.pt[1] = tmesh[i].cnct[(j + 1) % 4];
			edtmp.len = 1.;
			int flag(-1);
			for (k = 0; k < cp[tmesh[i].cnct[j]].face.size(); k++)
			{
				int fcid(cp[tmesh[i].cnct[j]].face[k]);
				if (fcid != i && tmesh[fcid].edge[0] != -1)
				{
					for (int k1 = 0; k1 < 4; k1++)
					{
						if ((tmedge[tmesh[fcid].edge[k1]].pt[0] == edtmp.pt[0] && tmedge[tmesh[fcid].edge[k1]].pt[1] == edtmp.pt[1]) ||
							(tmedge[tmesh[fcid].edge[k1]].pt[0] == edtmp.pt[1] && tmedge[tmesh[fcid].edge[k1]].pt[1] == edtmp.pt[0]))
						{
							flag = tmesh[fcid].edge[k1];
							break;
						}
					}
				}
				if (flag != -1)
					break;
			}
			if (flag == -1)
			{
				tmedge.push_back(edtmp);
				tmesh[i].edge[j] = tmedge.size() - 1;
			}
			else
			{
				tmesh[i].edge[j] = flag;
			}
		}
		for (j = 0; j < 4; j++)
		{
			tmedge[tmesh[i].edge[j]].face.push_back(i);
		}
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		cp[tmedge[i].pt[0]].edge.push_back(i);
		cp[tmedge[i].pt[1]].edge.push_back(i);
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].face.size() == 1)
		{
			int eid(tmedge[i].face[0]);
			if (tmesh[eid].type == 0)
				tmesh[eid].type = 2;
			else if (tmesh[eid].type == 2)
				tmesh[eid].type = 3;
			int *it = find(tmesh[eid].edge, tmesh[eid].edge + 4, i);
			int loc(it - tmesh[eid].edge);
			tmedge[tmesh[eid].edge[(loc + 1) % 4]].len = 0.;
			tmedge[tmesh[eid].edge[(loc + 3) % 4]].len = 0.;
		}
	}

	FindEdgeTopoDirec_1();
	FindKnotInterval_1();
	UpdateKnotInterval_1();
	SetLocalCoorSystem();
	//	FindIEN_3();
	FindIEN_4();
	Update_IEN_3();
}

void TruncatedTspline::SetDPatch_Irr()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		tmesh[eid].update = 0;
	}
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1 && tmesh[eid].type == 4)
		{
			for (uint i = 0; i < tmesh[eid].bemat.size(); i++)
			{
				tmesh[eid].bemat[i].clear();
			}
			tmesh[eid].bemat.clear();
			uint nv(cp[tmesh[eid].cnct[0]].face.size());
			SetBezierMat_DPatch22(SMOOTH_TYPE, nv, SMOOTH_BETA, tmesh[eid].bemat22);
		}
	}
}

void TruncatedTspline::ProjectVector(int type, int N, int eloc, double beta, vector<double> &pv)
{
	pv.clear();
	pv.resize(9, 0.);
	const double PI(3.141592653589793);
	double phi(2. * PI / double(N));
	double re = cos(-phi / 2.) - beta * sin(phi) * sin(-phi / 2.);
	double im = sin(-phi / 2.) + beta * sin(phi) * cos(-phi / 2.);
	double psi = atan2(im, re);
	if (type == 0)
	{
		pv[0] = 0.;
		pv[3] = 0.;
		pv[6] = 0.;
		pv[1] = 1. / (2. * double(N));
		pv[2] = pv[1];
		pv[4] = (1. + cos(double(eloc) * phi)) / (2. * double(N));
		pv[8] = pv[4];
		pv[5] = (1. + cos(2. * psi + double(eloc) * phi)) / (2. * double(N));
		pv[7] = (1. + cos(2. * psi - double(eloc) * phi)) / (2. * double(N));
	}
	else if (type == 1)
	{
		pv[0] = 1. / double(3. * N);
		pv[1] = pv[0];
		pv[2] = pv[0];
		pv[3] = pv[0];
		pv[6] = pv[0];
		pv[4] = (1. + 3. * cos(double(eloc) * phi)) / (3. * double(N));
		pv[8] = pv[4];
		pv[5] = (1. + 3. * cos(2. * psi + double(eloc) * phi)) / (3. * double(N));
		pv[7] = (1. + 3. * cos(2. * psi - double(eloc) * phi)) / (3. * double(N));
	}
	else
	{
		cout << "Other types of projections not available!\n";
	}
}

void TruncatedTspline::ProjectMatrix(int type, int N, double beta, vector<vector<double>> &pmat)
{
	//dimension of pmat: 3N by 3N, 9 submatrices
	for (uint i = 0; i < pmat.size(); i++)
	{
		pmat[i].clear();
	}
	pmat.clear();
	pmat.resize(3 * N, vector<double>(3 * N, 0.));
	vector<vector<double>> pv_all;
	for (int j = 0; j < N; j++)
	{
		vector<double> pv;
		ProjectVector(type, N, j, beta, pv);
		pv_all.push_back(pv);
	}
	for (int j = 0; j < N; j++)
	{
		for (int k = 0; k < N; k++)
		{
			int loc(0);
			for (int a = 0; a < 3; a++)
			{
				for (int b = 0; b < 3; b++) //0 to 8 (1 to 9)
				{
					pmat[a * N + j][b * N + k] = pv_all[(j - k + N) % N][loc]; //loc=3*a+b
					loc++;
				}
			}
		}
	}
}

void TruncatedTspline::VertexFaceTranMat(int N, vector<vector<double>> &bmat)
{
	for (uint i = 0; i < bmat.size(); i++)
		bmat[i].clear();
	bmat.clear();
	bmat.resize(2 * N + 1, vector<double>(3 * N, 0.));
	double a(4. / 9.), b(2. / 9.), c(1. / 9.);
	vector<array<int, 4>> fcloc(N);
	int fc0[3][4] = {{0, 5, 4, 3}, {0, 3, 2, 1}, {0, 7, 6, 5}};
	if (N == 3)
		fc0[2][1] = 1;
	for (int i = 0; i < 4; i++)
	{
		fcloc[0][i] = fc0[0][i];
		fcloc[1][i] = fc0[1][i];
		fcloc[N - 1][i] = fc0[2][i];
	}
	if (N >= 5)
	{
		fcloc[2][0] = 0;
		fcloc[2][1] = 1;
		fcloc[2][2] = 2 * N;
		fcloc[2][3] = 2 * N - 1;
		fcloc[N - 2][0] = 0;
		fcloc[N - 2][1] = 9;
		fcloc[N - 2][2] = 8;
		fcloc[N - 2][3] = 7;
	}
	for (int i = 3; i < N - 2; i++)
	{
		fcloc[i][0] = 0;
		fcloc[i][1] = 2 * N - 1 - 2 * (i - 3);
		fcloc[i][2] = 2 * N - 1 - 2 * (i - 3) - 1;
		fcloc[i][3] = 2 * N - 1 - 2 * (i - 3) - 2;
	}

	double df[3][4] = {{a, b, b, c}, {b, a, b, c}, {b, c, b, a}};
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				bmat[fcloc[i][k]][j * N + i] = df[j][k];
			}
		}
	}
}

void TruncatedTspline::SetBezierMat_DPatch(int type, int N, double beta, vector<vector<double>> &bemat)
{
	for (uint i = 0; i < bemat.size(); i++)
	{
		bemat[i].clear();
	}
	bemat.clear();
	int nb(2 * N + 8), n1r(2 * N + 1);
	bemat.resize(nb, vector<double>(16, 0.));
	vector<vector<double>> m1, m2, mb0;
	SetBezier3TranMat(N, mb0);		  //only used for outermost Bezier points
	VertexFaceTranMat(N, m1);		  //m1: (2N+1) by 3N
	ProjectMatrix(type, N, beta, m2); //m2: 3N by 3N

	vector<vector<double>> m12(n1r, vector<double>(3 * N, 0.));
	for (int i = 0; i < n1r; i++)
	{
		for (int j = 0; j < 3 * N; j++)
		{
			for (int k = 0; k < 3 * N; k++)
			{
				m12[i][k] += m1[i][j] * m2[k][j];
			}
		}
	}

	//vector<array<int, 4>> fcloc(N);
	//int fc0[3][4] = { { 0,5,4,3 },{ 0,3,2,1 },{ 0,7,6,5 } };
	//if (N == 3) fc0[2][1] = 1;
	//for (int i = 0; i < 4; i++)
	//{
	//	fcloc[0][i] = fc0[0][i];
	//	fcloc[1][i] = fc0[1][i];
	//	fcloc[N - 1][i] = fc0[2][i];
	//}
	//if (N >= 5)
	//{
	//	fcloc[2][0] = 0; fcloc[2][1] = 1; fcloc[2][2] = 2 * N; fcloc[2][3] = 2 * N - 1;
	//	fcloc[N - 2][0] = 0; fcloc[N - 2][1] = 9; fcloc[N - 2][2] = 8; fcloc[N - 2][3] = 7;
	//}
	//for (int i = 3; i < N - 2; i++)
	//{
	//	fcloc[i][0] = 0;
	//	fcloc[i][1] = 2 * N - 1 - 2 * (i - 3);
	//	fcloc[i][2] = 2 * N - 1 - 2 * (i - 3) - 1;
	//	fcloc[i][3] = 2 * N - 1 - 2 * (i - 3) - 2;
	//}

	/*cout << "m12:\n";
	for (int i = 0; i < N; i++)
	{
		cout << "i-th element: " << i << "\n";
		cout << m12[0][2*N+i] << " ";
		int ist= (1 - 2 * i + 2 * N) % (2 * N);
		for (int j = 1; j < n1r; j++)
		{
			int iloc = ist + j - 1;
			if (iloc > 2 * N)
			{
				iloc = iloc % (2 * N);
			}
			cout << m12[iloc][2*N+i] << " ";
		}
		cout << "\n";
		getchar();
	}*/

	//face
	for (int j = 0; j < n1r; j++)
	{
		//bemat[fcloc[0][j]][5] = m12[fcloc[0][j]][0];
		//bemat[fcloc[0][j]][6] = m12[fcloc[0][j]][N];
		//bemat[fcloc[0][j]][9] = m12[fcloc[0][j]][2 * N];
		bemat[j][5] = m12[j][0];
		bemat[j][6] = m12[j][N];
		bemat[j][9] = m12[j][2 * N];

		bemat[j][4] = (m12[j][0] + m12[j][1]) / 2.;
		bemat[j][8] = (m12[j][2 * N] + m12[j][N + 1]) / 2.;
		bemat[j][1] = (m12[j][0] + m12[j][N - 1]) / 2.;
		bemat[j][2] = (m12[j][N] + m12[j][3 * N - 1]) / 2.;

		for (int i = 0; i < N; i++)
		{
			bemat[j][0] += m12[j][i] / double(N);
		}
	}
	//cout << "bemat face:\n";
	//for (int i = 0; i < bemat[0].size(); i++)
	//{
	//	if (i == 5 || i == 6 || i == 9)
	//	{
	//		double sum(0.);
	//		for (int j = 0; j < bemat.size(); j++)
	//		{
	//			sum += bemat[j][i];
	//		}
	//		if (fabs(sum - 1.) > 1.e-6)
	//		{
	//			cout << sum << "\n";
	//			for (int j = 0; j < bemat.size(); j++)
	//			{
	//				cout << bemat[j][i] << " ";
	//			}
	//			getchar();
	//		}
	//		//double sum(0.);
	//		//for (int j = 0; j < m12.size(); j++)
	//		//{
	//		//	sum += m12[j][i];
	//		//}
	//		//if (fabs(sum - 1.) > 1.e-6)
	//		//{
	//		//	cout << sum << "\n";
	//		//	for (int j = 0; j < m12.size(); j++)
	//		//	{
	//		//		cout << m12[j][i] << " ";
	//		//	}
	//		//	getchar();
	//		//}
	//	}
	//}
	//getchar();

	//edge
	//for (int j = 0; j < 4; j++)
	//{
	//	bemat[fcloc[0][j]][4] += m12[fcloc[0][j]][0] / 2.;
	//	bemat[fcloc[1][j]][4] += m12[fcloc[1][j]][1] / 2.;
	//	bemat[fcloc[0][j]][8] += m12[fcloc[0][j]][2 * N] / 2.;
	//	bemat[fcloc[1][j]][8] += m12[fcloc[1][j]][N + 1] / 2.;
	//	bemat[fcloc[0][j]][1] += m12[fcloc[0][j]][0] / 2.;
	//	bemat[fcloc[N - 1][j]][1] += m12[fcloc[N - 1][j]][N - 1] / 2.;
	//	bemat[fcloc[0][j]][2] += m12[fcloc[0][j]][N] / 2.;
	//	bemat[fcloc[N - 1][j]][2] += m12[fcloc[N - 1][j]][3 * N - 1] / 2.;
	//}
	//EP corner
	//for (int i = 0; i < N; i++)
	//{
	//	for (int j = 0; j < 4; j++)
	//	{
	//		bemat[fcloc[i][j]][0] += m12[fcloc[i][j]][i] / double(N);
	//	}
	//}

	int id_same[8] = {3, 7, 10, 11, 12, 13, 14, 15};
	for (int i = 0; i < nb; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			bemat[i][id_same[j]] = mb0[i][id_same[j]]; //problematic
		}
	}

	//for (int j = 0; j < 16; j++)
	//{
	//	cout << "Bezier id: " << j << "\n";
	//	double sum(0.);
	//	for (int i = 0; i < nb; i++)
	//	{
	//		sum += bemat[i][j];
	//	}
	//	cout << sum << "\n";
	//	if (fabs(sum - 1.) > 1.e-8)
	//	{
	//		for (int i = 0; i < nb; i++)
	//		{
	//			cout << bemat[i][j] << " ";
	//			//if (i % 8 == 0 && i != 0)
	//			//	cout << "\n";
	//		}
	//	}
	//	getchar();
	//}
}

void TruncatedTspline::SetBezierMat_DPatch22(int type, int N, double beta, vector<vector<vector<double>>> &bemat)
{
	bemat.clear();
	bemat.resize(4);
	int nb(2 * N + 8), n1r(2 * N + 1);
	//bemat.resize(nb, vector<double>(16, 0.));
	vector<vector<double>> mb0;
	SetBezier3TranMat(N, mb0); //mb0: nb by 16
	double smat1d[7][4] = {{1., 0., 0., 0.}, {.5, .5, 0., 0.}, {.25, .5, .25, 0.}, {.125, .375, .375, .125}, {0., .25, .5, .25}, {0., 0., .5, .5}, {0., 0., 0., 1.}};
	vector<vector<double>> smat2d(49, vector<double>(16, 0.));
	int loc(0);
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			//cout << "new id: " << loc << "\n";
			int loc0(0);
			for (int i0 = 0; i0 < 4; i0++)
			{
				for (int j0 = 0; j0 < 4; j0++)
				{
					smat2d[loc][loc0] = smat1d[j][j0] * smat1d[i][i0];
					//cout << smat2d[loc][loc0] << " ";
					loc0++;
				}
			}
			//cout << "\n";
			//getchar();
			loc++;
		}
	}
	vector<vector<double>> mb1(nb, vector<double>(49, 0.));
	for (int i = 0; i < nb; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			for (int k = 0; k < 49; k++)
			{
				mb1[i][k] += mb0[i][j] * smat2d[k][j];
			}
		}
	}

	vector<vector<double>> m1(n1r, vector<double>(3 * N, 0.));
	for (int i = 0; i < N; i++)
	{
		m1[0][i] = mb1[0][8];
		m1[0][N + i] = mb1[0][9];
		m1[0][2 * N + i] = mb1[0][15];
		int ist = (1 - 2 * i + 2 * N) % (2 * N);
		//cout << "local eid: " << i << "\n";
		for (int j = 1; j < n1r; j++)
		{
			int iloc = ist + j - 1;
			if (iloc > 2 * N)
			{
				iloc = iloc % (2 * N);
			}
			m1[iloc][i] = mb1[j][8];
			m1[iloc][N + i] = mb1[j][9];
			m1[iloc][2 * N + i] = mb1[j][15];
			//cout << m1[iloc][i] << "(" << iloc << ") ";
		}
		//cout << "\n";
		//getchar();
	}

	vector<vector<double>> m2;
	ProjectMatrix(type, N, beta, m2); //3N by 3N

	vector<vector<double>> m12(n1r, vector<double>(3 * N, 0.));
	for (int i = 0; i < n1r; i++)
	{
		for (int j = 0; j < 3 * N; j++)
		{
			for (int k = 0; k < 3 * N; k++)
			{
				m12[i][k] += m1[i][j] * m2[k][j];
			}
		}
	}

	for (int j = 0; j < n1r; j++)
	{
		mb1[j][8] = m12[j][0];
		mb1[j][9] = m12[j][N];
		mb1[j][15] = m12[j][2 * N];

		mb1[j][7] = (m12[j][0] + m12[j][1]) / 2.;
		mb1[j][14] = (m12[j][2 * N] + m12[j][N + 1]) / 2.;
		mb1[j][1] = (m12[j][0] + m12[j][N - 1]) / 2.;
		mb1[j][2] = (m12[j][N] + m12[j][3 * N - 1]) / 2.;

		mb1[j][10] = (mb1[j][9] + mb1[j][11]) / 2.;
		mb1[j][22] = (mb1[j][15] + mb1[j][29]) / 2.;

		//mb1[j][3] = (m12[j][3 * N - 1] + mb1[j][29] + mb1[j][9] + mb1[j][11]) / 4.;
		//mb1[j][21] = (m12[j][N + 1] + mb1[j][15] + mb1[j][11] + mb1[j][29]) / 4.;

		mb1[j][0] = 0.;
		for (int i = 0; i < N; i++)
		{
			mb1[j][0] += m12[j][i] / double(N);
		}
	}
	mb1[0][3] = (m12[0][3 * N - 1] + mb1[0][29] + mb1[0][9] + mb1[0][11]) / 4.;
	mb1[0][21] = (m12[0][N + 1] + mb1[0][15] + mb1[0][11] + mb1[0][29]) / 4.;
	for (int j = 1; j < n1r; j++)
	{
		int iloc = j + 2;
		if (iloc > 2 * N)
		{
			iloc = iloc % (2 * N);
		}
		mb1[j][21] = (m12[j][N + 1] + mb1[j][15] + mb1[iloc][11] + mb1[j][29]) / 4.;
		iloc = j - 2;
		if (iloc <= 0)
		{
			iloc = iloc + 2 * N;
		}
		mb1[j][3] = (m12[j][3 * N - 1] + mb1[iloc][29] + mb1[j][9] + mb1[j][11]) / 4.;
	}

	//for (int i = 0; i < 49; i++)
	//{
	//	double sum(0.);
	//	for (int j = 0; j < nb; j++)
	//	{
	//		sum += mb1[j][i];
	//	}
	//	if (fabs(sum - 1.) > 1.e-6)
	//	{
	//		cout << "bezie id: " << i << "\n";
	//		cout << sum << "\n";
	//		for (int j = 0; j < nb; j++)
	//		{
	//			cout << mb1[j][i] << " ";
	//		}
	//		cout << "\n";
	//		getchar();
	//	}
	//}

	int bzcnct[4][16] = {{0, 1, 2, 3, 7, 8, 9, 10, 14, 15, 16, 17, 21, 22, 23, 24}, {3, 4, 5, 6, 10, 11, 12, 13, 17, 18, 19, 20, 24, 25, 26, 27}, {21, 22, 23, 24, 28, 29, 30, 31, 35, 36, 37, 38, 42, 43, 44, 45}, {24, 25, 26, 27, 31, 32, 33, 34, 38, 39, 40, 41, 45, 46, 47, 48}};
	for (int i = 0; i < 4; i++)
	{
		bemat[i].resize(nb, vector<double>(16, 0.));
		//bemat[i] = MatrixXd::Zero(nb,16);
		for (int j = 0; j < nb; j++)
		{
			for (int k = 0; k < 16; k++)
			{
				bemat[i][j][k] = mb1[j][bzcnct[i][k]];
			}
		}
	}

	//for (int i = 0; i < 4; i++)
	//{
	//	for (int j = 0; j < 16; j++)
	//	{
	//		double sum(0.);
	//		for (int k = 0; k < nb; k++)
	//		{
	//			sum += bemat[i](k, j);
	//		}
	//		if (fabs(sum - 1.) > 1.e-6)
	//		{
	//			cout << "subeid: " << i << "\n";
	//			cout << "bzid: " << j << "\n";
	//			cout << sum << "\n";
	//			for (int k = 0; k < nb; k++)
	//			{
	//				cout << bemat[i](k, j) << " ";
	//			}
	//			cout << "\n";
	//			getchar();
	//		}
	//	}
	//}
}

void TruncatedTspline::ElementBasis_Irregular_DPatch(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	//need to change variabe for Bezier, (u,v) -> (u_b,v_b)
	Nt.clear();
	dNdt.clear();
	Nt.resize(tmesh[eid].IEN.size(), 0.);
	dNdt.resize(tmesh[eid].IEN.size());
	vector<double> Nt1(tmesh[eid].IEN.size(), 0.);
	vector<array<double, 2>> dNdt1(tmesh[eid].IEN.size());
	uint nv(cp[tmesh[eid].cnct[0]].face.size());
	//vector<vector<double>> bmat(tmesh[eid].bemat.size(), vector<double>(tmesh[eid].bemat[0].size()));
	//for (uint i = 0; i<tmesh[eid].bemat.size(); i++)
	//{
	//	for (uint j = 0; j<tmesh[eid].bemat[i].size(); j++)
	//	{
	//		bmat[i][j] = tmesh[eid].bemat[i][j];
	//	}
	//}
	BezierElement be;
	double Nt0[16], dNdt0[16][2];
	double u_b(u / tmedge[tmesh[eid].edge[0]].len), v_b(v / tmedge[tmesh[eid].edge[3]].len);
	be.Basis(u_b, v_b, Nt0, dNdt0);
	//cout << "bemat size: " << tmesh[eid].bemat.size() << " " << tmesh[eid].bemat[0].size() << "\n";
	//getchar();
	for (uint i = 0; i < 2 * nv + 8; i++)
	{
		dNdt1[i][0] = 0.;
		dNdt1[i][1] = 0.;
		for (int j = 0; j < 16; j++)
		{
			//cout << tmesh[eid].bemat[i][j] << " ";
			//if (j % 8 == 0 && j != 0) cout << "\n";
			Nt1[i] += tmesh[eid].bemat[i][j] * Nt0[j];
			dNdt1[i][0] += tmesh[eid].bemat[i][j] * dNdt0[j][0];
			dNdt1[i][1] += tmesh[eid].bemat[i][j] * dNdt0[j][1];
		}
		//cout << "\n";
		//getchar();
		dNdt1[i][0] /= tmedge[tmesh[eid].edge[0]].len;
		dNdt1[i][1] /= tmedge[tmesh[eid].edge[3]].len;
	}
	/*vector<double> ku(5, 0.), kv(5, 0.), uval, vval;
	BSplineBasis bu, bv;
	for (uint i = 2 * nv + 1; i<tmesh[eid].IEN.size(); i++)
	{
		ku.assign(tmesh[eid].patch_ku[i - (2 * nv + 1)].begin(), tmesh[eid].patch_ku[i - (2 * nv + 1)].end());
		kv.assign(tmesh[eid].patch_kv[i - (2 * nv + 1)].begin(), tmesh[eid].patch_kv[i - (2 * nv + 1)].end());
		bu.Set(3, ku);
		bv.Set(3, kv);
		bu.BasisFunction(0, u, 1, uval);
		bv.BasisFunction(0, v, 1, vval);
		Nt1[i] = uval[0] * vval[0];
		dNdt1[i][0] = uval[1] * vval[0];
		dNdt1[i][1] = uval[0] * vval[1];
	}*/
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		Nt[i] = Nt1[i];
		dNdt[i][0] = dNdt1[i][0];
		dNdt[i][1] = dNdt1[i][1];
		//if (cp[tmesh[eid].IEN[i]].trun == 1)
		//{
		//	int pid(tmesh[eid].IEN[i]);
		//	for (uint j = 0; j<cp[pid].tbf.size(); j++)
		//	{
		//		vector<int>::iterator it = find(tmesh[eid].IEN.begin(), tmesh[eid].IEN.end(), cp[pid].tbf[j]);
		//		if (it != tmesh[eid].IEN.end())
		//		{
		//			int loc(it - tmesh[eid].IEN.begin());
		//			Nt[i] -= cp[pid].tc[j] * Nt1[loc];
		//			dNdt[i][0] -= cp[pid].tc[j] * dNdt1[loc][0];
		//			dNdt[i][1] -= cp[pid].tc[j] * dNdt1[loc][1];
		//		}
		//	}
		//}
	}
}

void TruncatedTspline::ElementBasis_Irregular_DPatch22(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	//need to change variabe for Bezier, (u,v) -> (u_b,v_b)
	Nt.clear();
	dNdt.clear();
	Nt.resize(tmesh[eid].IEN.size(), 0.);
	dNdt.resize(tmesh[eid].IEN.size());
	vector<double> Nt1(tmesh[eid].IEN.size(), 0.);
	vector<array<double, 2>> dNdt1(tmesh[eid].IEN.size());
	uint nv(cp[tmesh[eid].cnct[0]].face.size());
	//vector<vector<double>> bmat(tmesh[eid].bemat.size(), vector<double>(tmesh[eid].bemat[0].size()));
	//for (uint i = 0; i<tmesh[eid].bemat.size(); i++)
	//{
	//	for (uint j = 0; j<tmesh[eid].bemat[i].size(); j++)
	//	{
	//		bmat[i][j] = tmesh[eid].bemat[i][j];
	//	}
	//}
	BezierElement be;
	double Nt0[16], dNdt0[16][2];
	double edm[2] = {tmedge[tmesh[eid].edge[0]].len, tmedge[tmesh[eid].edge[3]].len};
	double u_b(u / edm[0]), v_b(v / edm[1]); //problematic
	double uhalf(.5 * edm[0]), vhalf(.5 * edm[1]);
	int subeid(0);
	if (u_b >= uhalf && v_b < vhalf)
	{
		subeid = 1;
		u_b = 2. * (u_b - uhalf);
		v_b = 2. * v_b;
	}
	else if (u_b < uhalf && v_b >= vhalf)
	{
		subeid = 2;
		u_b = 2. * u_b;
		v_b = 2. * (v_b - vhalf);
	}
	else if (u_b >= uhalf && v_b >= vhalf)
	{
		subeid = 3;
		u_b = 2. * (u_b - uhalf);
		v_b = 2. * (v_b - vhalf);
	}
	else
	{
		subeid = 0;
		u_b = 2. * u_b;
		v_b = 2. * v_b;
	}
	//cout << u << " " << v << "\n";
	//cout << u_b << " " << v_b << "\n";
	//cout << subeid << "\n";
	//getchar();
	be.Basis(u_b, v_b, Nt0, dNdt0);
	//cout << "bemat size: " << tmesh[eid].bemat.size() << " " << tmesh[eid].bemat[0].size() << "\n";
	//getchar();
	//double sum(0.);
	for (uint i = 0; i < 2 * nv + 8; i++)
	{
		dNdt1[i][0] = 0.;
		dNdt1[i][1] = 0.;
		for (int j = 0; j < 16; j++)
		{
			Nt1[i] += tmesh[eid].bemat22[subeid][i][j] * Nt0[j];
			dNdt1[i][0] += tmesh[eid].bemat22[subeid][i][j] * dNdt0[j][0];
			dNdt1[i][1] += tmesh[eid].bemat22[subeid][i][j] * dNdt0[j][1];
		}
		//for (int j = 0; j<16; j++)
		//{
		//	Nt1[i] += tmesh[eid].bemat[i][j] * Nt0[j];
		//	dNdt1[i][0] += tmesh[eid].bemat[i][j] * dNdt0[j][0];
		//	dNdt1[i][1] += tmesh[eid].bemat[i][j] * dNdt0[j][1];
		//}
		//sum += Nt1[i];
		//cout << "\n";
		//getchar();
		//dNdt1[i][0] /= tmedge[tmesh[eid].edge[0]].len;
		//dNdt1[i][1] /= tmedge[tmesh[eid].edge[3]].len;
		dNdt1[i][0] *= uhalf;
		dNdt1[i][1] *= vhalf;
	}
	/*vector<double> ku(5, 0.), kv(5, 0.), uval, vval;
	BSplineBasis bu, bv;
	for (uint i = 2 * nv + 1; i<tmesh[eid].IEN.size(); i++)
	{
		ku.assign(tmesh[eid].patch_ku[i - (2 * nv + 1)].begin(), tmesh[eid].patch_ku[i - (2 * nv + 1)].end());
		kv.assign(tmesh[eid].patch_kv[i - (2 * nv + 1)].begin(), tmesh[eid].patch_kv[i - (2 * nv + 1)].end());
		bu.Set(3, ku);
		bv.Set(3, kv);
		bu.BasisFunction(0, u, 1, uval);
		bv.BasisFunction(0, v, 1, vval);
		Nt1[i] = uval[0] * vval[0];
		dNdt1[i][0] = uval[1] * vval[0];
		dNdt1[i][1] = uval[0] * vval[1];
		sum += Nt1[i];
	}*/
	//cout << sum << "\n";
	//getchar();
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		Nt[i] = Nt1[i];
		dNdt[i][0] = dNdt1[i][0];
		dNdt[i][1] = dNdt1[i][1];
		//if (cp[tmesh[eid].IEN[i]].trun == 1)
		//{
		//	int pid(tmesh[eid].IEN[i]);
		//	for (uint j = 0; j<cp[pid].tbf.size(); j++)
		//	{
		//		vector<int>::iterator it = find(tmesh[eid].IEN.begin(), tmesh[eid].IEN.end(), cp[pid].tbf[j]);
		//		if (it != tmesh[eid].IEN.end())
		//		{
		//			int loc(it - tmesh[eid].IEN.begin());
		//			Nt[i] -= cp[pid].tc[j] * Nt1[loc];
		//			dNdt[i][0] -= cp[pid].tc[j] * dNdt1[loc][0];
		//			dNdt[i][1] -= cp[pid].tc[j] * dNdt1[loc][1];
		//		}
		//	}
		//}
	}
}

void TruncatedTspline::OutputBezierPoint22(string fn)
{
	vector<array<double, 3>> pts0(16), pts1(16), pts2(16), pts3(16);
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].type == 4 && eid == 116)
		{
			for (int i = 0; i < 16; i++)
			{
				pts0[i][0] = 0.;
				pts0[i][1] = 0.;
				pts0[i][2] = 0.;
				pts1[i][0] = 0.;
				pts1[i][1] = 0.;
				pts1[i][2] = 0.;
				pts2[i][0] = 0.;
				pts2[i][1] = 0.;
				pts2[i][2] = 0.;
				pts3[i][0] = 0.;
				pts3[i][1] = 0.;
				pts3[i][2] = 0.;
				for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
				{
					pts0[i][0] += tmesh[eid].bemat22[0][j][i] * cp[tmesh[eid].IEN[j]].coor[0];
					pts0[i][1] += tmesh[eid].bemat22[0][j][i] * cp[tmesh[eid].IEN[j]].coor[1];
					pts0[i][2] += tmesh[eid].bemat22[0][j][i] * cp[tmesh[eid].IEN[j]].coor[2];

					pts1[i][0] += tmesh[eid].bemat22[1][j][i] * cp[tmesh[eid].IEN[j]].coor[0];
					pts1[i][1] += tmesh[eid].bemat22[1][j][i] * cp[tmesh[eid].IEN[j]].coor[1];
					pts1[i][2] += tmesh[eid].bemat22[1][j][i] * cp[tmesh[eid].IEN[j]].coor[2];

					pts2[i][0] += tmesh[eid].bemat22[2][j][i] * cp[tmesh[eid].IEN[j]].coor[0];
					pts2[i][1] += tmesh[eid].bemat22[2][j][i] * cp[tmesh[eid].IEN[j]].coor[1];
					pts2[i][2] += tmesh[eid].bemat22[2][j][i] * cp[tmesh[eid].IEN[j]].coor[2];

					pts3[i][0] += tmesh[eid].bemat22[3][j][i] * cp[tmesh[eid].IEN[j]].coor[0];
					pts3[i][1] += tmesh[eid].bemat22[3][j][i] * cp[tmesh[eid].IEN[j]].coor[1];
					pts3[i][2] += tmesh[eid].bemat22[3][j][i] * cp[tmesh[eid].IEN[j]].coor[2];
				}
			}
			//break;
		}
	}

	string fname = fn + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		int npt = pts0.size() + pts1.size() + pts2.size() + pts3.size();
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << npt << " float\n";
		for (uint i = 0; i < pts0.size(); i++)
		{
			fout << pts0[i][0] << " " << pts0[i][1] << " " << pts0[i][2] << "\n";
		}
		for (uint i = 0; i < pts1.size(); i++)
		{
			fout << pts1[i][0] << " " << pts1[i][1] << " " << pts1[i][2] << "\n";
		}
		for (uint i = 0; i < pts2.size(); i++)
		{
			fout << pts2[i][0] << " " << pts2[i][1] << " " << pts2[i][2] << "\n";
		}
		for (uint i = 0; i < pts3.size(); i++)
		{
			fout << pts3[i][0] << " " << pts3[i][1] << " " << pts3[i][2] << "\n";
		}
		fout << "\nCELLS " << npt << " " << 2 * npt << '\n';
		for (int i = 0; i < npt; i++)
		{
			fout << "1 " << i << '\n';
		}
		fout << "\nCELL_TYPES " << npt << '\n';
		for (int i = 0; i < npt; i++)
		{
			fout << "1\n";
		}
		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i]<<"\n";
		//}
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

void TruncatedTspline::OutputBezierMesh(const vector<BezierElement> &bzmesh, string fn)
{
	int bzloc[4] = {0, 3, 15, 12};
	string fname = fn + "_bzmesh.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << 4 * bzmesh.size() << " float\n";
		for (uint i = 0; i < bzmesh.size(); i++)
		{
			for (int j = 0; j < 4; j++)
			{
				fout << bzmesh[i].pts[bzloc[j]][0] << " " << bzmesh[i].pts[bzloc[j]][1] << " " << bzmesh[i].pts[bzloc[j]][2] << "\n";
			}
		}
		fout << "\nCELLS " << bzmesh.size() << " " << 5 * bzmesh.size() << '\n';
		for (int i = 0; i < bzmesh.size(); i++)
		{
			fout << "4 " << 4 * i << " " << 4 * i + 1 << " " << 4 * i + 2 << " " << 4 * i + 3 << '\n';
		}
		fout << "\nCELL_TYPES " << bzmesh.size() << '\n';
		for (int i = 0; i < bzmesh.size(); i++)
		{
			fout << "9\n";
		}

		/*fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << 16 * bzmesh.size() << " float\n";
		for (uint i = 0; i < bzmesh.size(); i++)
		{
			for (int j = 0; j < 16; j++)
			{
				fout << bzmesh[i].pts[j][0] << " " << bzmesh[i].pts[j][1] << " " << bzmesh[i].pts[j][2] << "\n";
			}
		}
		fout << "\nCELLS " << 9 * bzmesh.size() << " " << 45 * bzmesh.size() << '\n';
		for (int i = 0; i<bzmesh.size(); i++)
		{
			int ist(16 * i);
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					fout << "4 " << ist+4*j+k << " " << ist + 4 * j + k+1 << " " << ist + 4 * (j+1) + k + 1 << " " << ist + 4 * (j + 1) + k << '\n';
				}
			}
		}
		fout << "\nCELL_TYPES " << 9*bzmesh.size() << '\n';
		for (int i = 0; i<9*bzmesh.size(); i++)
		{
			fout << "9\n";
		}*/

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline::BezierExtract_DPatch(vector<BezierElement> &bzmesh)
{
	cout << "Bezier extracting...\n";
	bzmesh.clear();
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1)
		{
			if (tmesh[eid].type == 0 || tmesh[eid].type == 1)
			{
				BezierElementExtract_Unstruct(eid, bzmesh);
			}
			else if (tmesh[eid].type == 4)
			{
				BezierUnit_DPatch22(eid, bzmesh);
			}
		}
	}
}

void TruncatedTspline::BezierExtract_DPatch_AS(vector<BezierElement> &bzmesh) //Analysis space
{
	cout << "Bezier extracting...\n";
	bzmesh.clear();
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1)
		{
			if (tmesh[eid].c1 == 0 && (tmesh[eid].type == 0 || tmesh[eid].type == 1))
			{
				BezierElementExtract_Unstruct(eid, bzmesh);
			}
			else if (tmesh[eid].c1 == 1 && tmesh[eid].type == 4)
			{
				BezierUnit_DPatch_Irr22(eid, bzmesh);
			}
			else if (tmesh[eid].c1 == 1 && tmesh[eid].type != 4)
			{
				BezierUnit_DPatch_Trs(eid, bzmesh);
			}
			else if (tmesh[eid].c1 == 2)
			{
				BezierUnit_DPatch_Trs(eid, bzmesh);
			}
		}
	}
}

void TruncatedTspline::BezierUnit_DPatch22(int eid, vector<BezierElement> &bzmesh)
{
	vector<int> pid = tmesh[eid].IEN;
	for (int is = 0; is < 4; is++)
	{
		BezierElement bzel;
		bzel.prt = eid;
		bzel.cmat.resize(pid.size());
		for (uint j = 0; j < pid.size(); j++)
		{
			for (int k = 0; k < 16; k++)
			{
				bzel.cmat[j][k] = tmesh[eid].bemat22[is][j][k];
			}
		}
		for (uint i = 0; i < pid.size(); i++)
		{
			for (int j = 0; j < 16; j++)
			{
				bzel.pts[j][0] += bzel.cmat[i][j] * cp[pid[i]].coor[0];
				bzel.pts[j][1] += bzel.cmat[i][j] * cp[pid[i]].coor[1];
				bzel.pts[j][2] += bzel.cmat[i][j] * cp[pid[i]].coor[2];
			}
		}
		for (uint i = 0; i < pid.size(); i++)
		{
			bzel.IEN.push_back(pid[i]);
		}
		bzmesh.push_back(bzel);
	}
}

void TruncatedTspline::BezierUnit_DPatch_Irr22(int eid, vector<BezierElement> &bzmesh)
{
	vector<int> pid = tmesh[eid].IEN;
	for (uint i = 0; i < tmesh[eid].IENc1.size(); i++)
	{
		pid.push_back(cp.size() + tmesh[eid].IENc1[i]);
	}
	for (int is = 0; is < 4; is++)
	{
		BezierElement bzel;
		bzel.prt = eid;
		bzel.cmat.resize(pid.size());
		bzel.IEN = pid;
		for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
		{
			for (int k = 0; k < 16; k++)
			{
				bzel.cmat[j][k] = tmesh[eid].bemat22[is][j][k];
				bzel.pts[k][0] += bzel.cmat[j][k] * cp[pid[j]].coor[0];
				bzel.pts[k][1] += bzel.cmat[j][k] * cp[pid[j]].coor[1];
				bzel.pts[k][2] += bzel.cmat[j][k] * cp[pid[j]].coor[2];
			}
		}
		int ist(tmesh[eid].IEN.size());
		//if (tmesh[eid].IENc1.size() != tmesh[eid].c1mat22[is].size())
		//{
		//	cout << "eid: " << eid << "\n";
		//	cout << "nv: " << cp[tmesh[eid].cnct[0]].face.size() << "\n";
		//	cout << tmesh[eid].c1mat22[is].size() << " " << tmesh[eid].IENc1.size() << "\n"; getchar();
		//}
		for (uint j = 0; j < tmesh[eid].IENc1.size(); j++)
		{
			int jloc(ist + j);
			for (int k = 0; k < 16; k++)
			{
				bzel.cmat[jloc][k] = tmesh[eid].c1mat22[is][j][k];
				bzel.pts[k][0] += bzel.cmat[jloc][k] * bzcp[pid[jloc] - cp.size()][0];
				bzel.pts[k][1] += bzel.cmat[jloc][k] * bzcp[pid[jloc] - cp.size()][1];
				bzel.pts[k][2] += bzel.cmat[jloc][k] * bzcp[pid[jloc] - cp.size()][2];
			}
		}
		bzmesh.push_back(bzel);

		//if (is == 0)
		//{
		//	/*cout << "bzid-3:\n";
		//	int bzid(3);
		//	cout << "c2:\n";
		//	for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
		//	{
		//		cout << bzel.cmat[j][bzid] << "(" << j << ") ";
		//	}
		//	cout << "\n";
		//	cout << "c1:\n";
		//	for (uint j = 0; j < tmesh[eid].IENc1.size(); j++)
		//	{
		//		cout << bzel.cmat[ist+j][bzid] << "(" << j << ") ";
		//	}
		//	cout << "\n";
		//	cout << "\nbzid-12:\n";
		//	bzid = 12;
		//	cout << "c2:\n";
		//	for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
		//	{
		//		cout << bzel.cmat[j][bzid] << "(" << j << ") ";
		//	}
		//	cout << "\n";
		//	cout << "c1:\n";
		//	for (uint j = 0; j < tmesh[eid].IENc1.size(); j++)
		//	{
		//		cout << bzel.cmat[ist + j][bzid] << "(" << j << ") ";
		//	}
		//	cout << "\n";*/

		//	cout << "eid: " << eid << "\n";
		//	cout << "c2:\n";
		//	for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
		//	{
		//		cout << tmesh[eid].IEN[j] << " ";
		//	}
		//	cout << "\n";
		//	cout << "c1:\n";
		//	for (uint j = 0; j < tmesh[eid].IENc1.size(); j++)
		//	{
		//		cout << tmesh[eid].IENc1[j]+cp.size() << " ";
		//	}
		//	cout << "\n";

		//	getchar();
		//}
	}
}

void TruncatedTspline::BezierUnit_DPatch_Trs(int eid, vector<BezierElement> &bzmesh)
{
	vector<int> pid = tmesh[eid].IEN;
	for (uint i = 0; i < tmesh[eid].IENc1.size(); i++)
	{
		pid.push_back(cp.size() + tmesh[eid].IENc1[i]);
	}
	BezierElement bzel;
	bzel.prt = eid;
	bzel.cmat.resize(pid.size());
	bzel.IEN = pid;
	for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
	{
		for (int k = 0; k < 16; k++)
		{
			bzel.cmat[j][k] = tmesh[eid].bemat[j][k];
			bzel.pts[k][0] += bzel.cmat[j][k] * cp[pid[j]].coor[0];
			bzel.pts[k][1] += bzel.cmat[j][k] * cp[pid[j]].coor[1];
			bzel.pts[k][2] += bzel.cmat[j][k] * cp[pid[j]].coor[2];
		}
	}
	int ist(tmesh[eid].IEN.size());
	for (uint j = 0; j < tmesh[eid].IENc1.size(); j++)
	{
		int jloc(ist + j);
		for (int k = 0; k < 16; k++)
		{
			bzel.cmat[jloc][k] = tmesh[eid].c1mat[j][k];
			bzel.pts[k][0] += bzel.cmat[jloc][k] * bzcp[pid[jloc] - cp.size()][0];
			bzel.pts[k][1] += bzel.cmat[jloc][k] * bzcp[pid[jloc] - cp.size()][1];
			bzel.pts[k][2] += bzel.cmat[jloc][k] * bzcp[pid[jloc] - cp.size()][2];
		}
	}
	bzmesh.push_back(bzel);
}

void TruncatedTspline::UpdateCP_Fitting(const vector<array<double, 3>> &cpnew)
{
	int loc(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			cp[i].coor[0] = cpnew[loc][0];
			cp[i].coor[1] = cpnew[loc][1];
			//cp[i].coor[2] = cpnew[loc][2];
			loc++;
		}
	}
	for (uint i = 0; i < bzcp.size(); i++)
	{
		bzcp[i][0] = cpnew[loc][0];
		bzcp[i][1] = cpnew[loc][1];
		//bzcp[i][2] = cpnew[loc][2];
		loc++;
	}
}

void TruncatedTspline::WriteBezier_Shell(const vector<BezierElement> &bzmesh, string fn)
{
	vector<int> bcflag(cp.size(), 0);
	vector<int> bcnb(cp.size(), -1);
	for (uint i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].act == 1 && tmedge[i].face.size() == 1)
		{
			bcflag[tmedge[i].pt[0]] = 1;
			bcflag[tmedge[i].pt[1]] = 1;
		}
	}
	for (uint i = 0; i < cp.size(); i++)
	{
		if (bcflag[i] == 1)
		{
			int edid(-1);
			for (uint j = 0; j < cp[i].edge.size(); j++)
			{
				int edid1(cp[i].edge[j]);
				if (tmedge[edid1].face.size() != 1)
				{
					edid = edid1;
					break;
				}
			}
			if (edid != -1)
			{
				bcnb[i] = tmedge[edid].pt[0];
				if (bcnb[i] == i)
					bcnb[i] = tmedge[edid].pt[1];
			}
		}
	}

	////actives, the same as paid
	//vector<int> IDa(cp.size() + bzcp.size(), -1);
	////vector<int> bzcpaID(bzcp.size(), -1);
	//int count(0);
	//for (uint i = 0; i < cp.size(); i++)
	//{
	//	if (cp[i].act == 1) IDa[i] = count++;
	//}
	//int ncpa(count);
	//for (uint i = 0; i < bzcp.size(); i++)
	//{
	//	IDa[cp.size() + i] = count++;
	//}

	cout << "Writing file...\n";
	string fname = fn + "_bezex.txt";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "nsd " << 3 << "\n";
		fout << "degree " << 3 << " " << 3 << "\n";
		fout << "funcs " << npta << "\n";
		fout.precision(16);
		fout << fixed;
		double w(1.);
		for (uint i = 0; i < cp.size(); i++)
		{
			if (cp[i].act == 1)
			{
				fout << paid[i] << " " << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << " " << w << " ";
				if (bcnb[i] == -1)
					fout << "-1\n";
				else
					fout << paid[bcnb[i]] << "\n";
			}
		}
		for (uint i = 0; i < bzcp.size(); i++)
		{
			fout << paid[cp.size() + i] << " " << bzcp[i][0] << " " << bzcp[i][1] << " " << bzcp[i][2] << " " << w << " " << -1 << "\n";
		}
		fout << "elems " << bzmesh.size() << "\n";
		for (uint i = 0; i < bzmesh.size(); i++)
		{
			if (i != 0 && i % 500 == 0)
			{
				cout << i << " ";
			}
			vector<int> loc;
			vector<int> nnzv;
			//if (bzmesh[i].order == 3)
			{
				for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
				{
					if (bzmesh[i].IEN[j] >= 0 && paid[bzmesh[i].IEN[j]] != -1)
					{
						int nnz(0);
						for (uint k = 0; k < bzmesh[i].cmat[j].size(); k++)
						{
							if (bzmesh[i].cmat[j][k] != 0.)
								nnz++;
						}
						if (nnz != 0)
						{
							loc.push_back(j);
							nnzv.push_back(nnz);
						}
					}
				}
			}

			fout << i << " " << loc.size() << "\n";
			for (uint j = 0; j < loc.size(); j++)
			{
				fout << paid[bzmesh[i].IEN[loc[j]]] << " ";
			}
			fout << "\n";
			//if (bzmesh[i].order == 3)
			{
				for (uint j = 0; j < loc.size(); j++)
				{
					for (uint k = 0; k < bzmesh[i].cmat[loc[j]].size(); k++)
						fout << bzmesh[i].cmat[loc[j]][k] << " ";
					fout << "\n";
				}
			}
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	//boundary
	fname.clear();
	fname = fn + "_boundarycp.txt";
	string fname2 = fn + "_boundarycp_inner_layer.txt";
	OutputBoundaryCPID_SimpleSquare(fname, fname2);
	//OutputBoundaryCPID_Hemisphere(fname);
	//OutputBoundaryCPID_Crosspipe(fname);
	//OutputBoundaryCPID_Cylinder(fname);

	//corner
	fname.clear();
	fname = fn + "_corner.txt";
	vector<array<double, 3>> corn;
	GetCornerInfo_Hemisphere(corn);
	//GetCornerInfo_Cylinder(corn);
	OutputCorner(bzmesh, corn, 1.e-3, fname);

	cout << "End of writing!\n";
}

void TruncatedTspline::OutputBoundaryCPID_SimpleSquare(string fn_bc, string fn_bcin)
{
	//for a quarter of hemisphere
	vector<vector<int>> bcpid(8);
	double tol(1.e-4);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (fabs(cp[i].coor[0]) < tol)
		{
			bcpid[0].push_back(i);
		}
		if (fabs(1.0 - cp[i].coor[0]) < tol)
		{
			bcpid[1].push_back(i);
		}
		if (fabs(cp[i].coor[0]) < tol)
		{
			bcpid[2].push_back(i);
		}
		if (fabs(1.0 - cp[i].coor[1]) < tol)
		{
			bcpid[3].push_back(i);
		}
		if (fabs(cp[i].coor[0] - 0.01) < tol)
		{
			bcpid[4].push_back(i);
		}
		if (fabs(0.99 - cp[i].coor[0]) < tol)
		{
			bcpid[5].push_back(i);
		}
		if (fabs(cp[i].coor[0] - 0.01) < tol)
		{
			bcpid[6].push_back(i);
		}
		if (fabs(0.99 - cp[i].coor[1]) < tol)
		{
			bcpid[7].push_back(i);
		}
	}

	ofstream fout;
	fout.open(fn_bc.c_str());
	if (fout.is_open())
	{
		fout << "Left " << bcpid[0].size() << "\n";
		for (uint i = 0; i < bcpid[0].size(); i++)
		{
			fout << paid[bcpid[0][i]] << " ";
		}
		fout << "\n";
		fout << "Right " << bcpid[1].size() << "\n";
		for (uint i = 0; i < bcpid[1].size(); i++)
		{
			fout << paid[bcpid[1][i]] << " ";
		}
		fout << "\n";
		fout << "Bottom " << bcpid[2].size() << "\n";
		for (uint i = 0; i < bcpid[2].size(); i++)
		{
			fout << paid[bcpid[2][i]] << " ";
		}
		fout << "\n";
		fout << "Top " << bcpid[3].size() << "\n";
		for (uint i = 0; i < bcpid[3].size(); i++)
		{
			fout << paid[bcpid[3][i]] << " ";
		}
		fout << "\n";
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fn_bc << "!\n";
	}

	fout.open(fn_bcin.c_str());
	if (fout.is_open())
	{
		fout << "Left " << bcpid[4].size() << "\n";
		for (uint i = 0; i < bcpid[4].size(); i++)
		{
			fout << paid[bcpid[4][i]] << " ";
		}
		fout << "\n";
		fout << "Right " << bcpid[5].size() << "\n";
		for (uint i = 0; i < bcpid[5].size(); i++)
		{
			fout << paid[bcpid[5][i]] << " ";
		}
		fout << "\n";
		fout << "Bottom " << bcpid[6].size() << "\n";
		for (uint i = 0; i < bcpid[6].size(); i++)
		{
			fout << paid[bcpid[6][i]] << " ";
		}
		fout << "\n";
		fout << "Top " << bcpid[7].size() << "\n";
		for (uint i = 0; i < bcpid[7].size(); i++)
		{
			fout << paid[bcpid[7][i]] << " ";
		}
		fout << "\n";
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fn_bcin << "!\n";
	}
}

void TruncatedTspline::OutputBoundaryCPID_Hemisphere(string fn)
{
	//for a quarter of hemisphere
	vector<vector<int>> bcpid(3);
	double tol(1.e-4);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (fabs(cp[i].coor[1]) < tol)
		{
			bcpid[0].push_back(i);
		}
		if (fabs(cp[i].coor[0]) < tol)
		{
			bcpid[1].push_back(i);
		}
		if (fabs(cp[i].coor[2]) < tol)
		{
			bcpid[2].push_back(i);
		}
	}

	ofstream fout;
	fout.open(fn.c_str());
	if (fout.is_open())
	{
		fout << "Left " << bcpid[0].size() << "\n";
		for (uint i = 0; i < bcpid[0].size(); i++)
		{
			fout << paid[bcpid[0][i]] << " ";
		}
		fout << "\n";
		fout << "Right " << bcpid[1].size() << "\n";
		for (uint i = 0; i < bcpid[1].size(); i++)
		{
			fout << paid[bcpid[1][i]] << " ";
		}
		fout << "\n";
		fout << "Bottom " << bcpid[2].size() << "\n";
		for (uint i = 0; i < bcpid[2].size(); i++)
		{
			fout << paid[bcpid[2][i]] << " ";
		}
		fout << "\n";
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fn << "!\n";
	}
}

void TruncatedTspline::OutputBoundaryCPID_Crosspipe(string fn)
{
	//for a quarter of cylinder
	vector<vector<int>> bcpid(4);
	double tol(1.e-4);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (fabs(cp[i].coor[0]) < tol && fabs(cp[i].coor[1] - (-375.0)) < tol) //left: y==-0.0519419
		{
			bcpid[0].push_back(i);
		}
		if (fabs(cp[i].coor[0]) < tol && fabs(cp[i].coor[1] - (375.0)) < tol) //right: y==0.0519419
		{
			bcpid[1].push_back(i);
		}
		if (fabs(cp[i].coor[2] - 1200.0) < tol) //top: z==0.25146
		{
			bcpid[2].push_back(i);
		}
		if (fabs(cp[i].coor[2] - 0.0) < tol) //bottom: z==0.00254
		{
			bcpid[3].push_back(i);
		}
	}

	ofstream fout;
	fout.open(fn.c_str());
	if (fout.is_open())
	{
		fout << "Left " << bcpid[0].size() << "\n";
		for (uint i = 0; i < bcpid[0].size(); i++)
		{
			fout << paid[bcpid[0][i]] << " ";
		}
		fout << "\n";
		fout << "Right " << bcpid[1].size() << "\n";
		for (uint i = 0; i < bcpid[1].size(); i++)
		{
			fout << paid[bcpid[1][i]] << " ";
		}
		fout << "\n";
		fout << "Top " << bcpid[2].size() << "\n";
		for (uint i = 0; i < bcpid[2].size(); i++)
		{
			fout << paid[bcpid[2][i]] << " ";
		}
		fout << "\n";
		fout << "Bottom " << bcpid[3].size() << "\n";
		for (uint i = 0; i < bcpid[3].size(); i++)
		{
			fout << paid[bcpid[3][i]] << " ";
		}
		fout << "\n";
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fn << "!\n";
	}
}

void TruncatedTspline::OutputBoundaryCPID_Cylinder(string fn)
{
	//for a quarter of cylinder
	vector<vector<int>> bcpid(4);
	double tol(1.e-4);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (fabs(cp[i].coor[1]) < tol) //left: y==0
		{
			bcpid[0].push_back(i);
		}
		if (fabs(cp[i].coor[0]) < tol) //right: x==0
		{
			bcpid[1].push_back(i);
		}
		if (fabs(cp[i].coor[2] - 100.) < tol) //top: z==100
		{
			bcpid[2].push_back(i);
		}
		if (fabs(cp[i].coor[2]) < tol) //bottom: z==0
		{
			bcpid[3].push_back(i);
		}
	}

	ofstream fout;
	fout.open(fn.c_str());
	if (fout.is_open())
	{
		fout << "Left " << bcpid[0].size() << "\n";
		for (uint i = 0; i < bcpid[0].size(); i++)
		{
			fout << paid[bcpid[0][i]] << " ";
		}
		fout << "\n";
		fout << "Right " << bcpid[1].size() << "\n";
		for (uint i = 0; i < bcpid[1].size(); i++)
		{
			fout << paid[bcpid[1][i]] << " ";
		}
		fout << "\n";
		fout << "Top " << bcpid[2].size() << "\n";
		for (uint i = 0; i < bcpid[2].size(); i++)
		{
			fout << paid[bcpid[2][i]] << " ";
		}
		fout << "\n";
		fout << "Bottom " << bcpid[3].size() << "\n";
		for (uint i = 0; i < bcpid[3].size(); i++)
		{
			fout << paid[bcpid[3][i]] << " ";
		}
		fout << "\n";
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fn << "!\n";
	}
}

void TruncatedTspline::OutputCorner(const vector<BezierElement> &bzmesh, const vector<array<double, 3>> &corn, double tol, string fn)
{
	int bzcn[4] = {0, 3, 15, 12};
	double bzpar[4][2] = {{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}};
	vector<int> cnid;
	vector<array<double, 2>> cnu;

	for (uint ic = 0; ic < corn.size(); ic++)
	{
		for (uint i = 0; i < bzmesh.size(); i++) //hemisphere first corner
		{
			int flag(0);
			for (int j = 0; j < 4; j++)
			{
				if (fabs(bzmesh[i].pts[bzcn[j]][0] - corn[ic][0]) < tol &&
					fabs(bzmesh[i].pts[bzcn[j]][1] - corn[ic][1]) < tol &&
					fabs(bzmesh[i].pts[bzcn[j]][2] - corn[ic][2]) < tol)
				{
					cnid.push_back(i);
					array<double, 2> partmp = {bzpar[j][0], bzpar[j][1]};
					cnu.push_back(partmp);
					flag = 1;
					break;
				}
			}
			if (flag == 1)
				break;
		}
	}

	ofstream fout;
	fout.open(fn.c_str());
	if (fout.is_open())
	{
		for (uint i = 0; i < cnid.size(); i++)
		{
			fout << cnid[i] << "\t" << cnu[i][0] << "\t" << cnu[i][1] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fn << "!\n";
	}
}

void TruncatedTspline::GetCornerInfo_Hemisphere(vector<array<double, 3>> &corn)
{
	corn.clear();
	corn.resize(2);
	double cntmp[2][3] = {{10., 0., 0.}, {0., 10., 0.}};
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			corn[i][j] = cntmp[i][j];
		}
	}
}

void TruncatedTspline::GetCornerInfo_Cylinder(vector<array<double, 3>> &corn)
{
	corn.clear();
	corn.resize(2);
	double cntmp[2][3] = {{100., 0., 100.}, {0., 100., 100.}};
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			corn[i][j] = cntmp[i][j];
		}
	}
}

void TruncatedTspline::ReadBezierMat(string fn, vector<BezierElement> &bzmesh)
{
	ifstream fin;
	fin.open(fn.c_str());
	if (fin.is_open())
	{
		string stmp;
		int itmp, n, nIEN;
		fin >> stmp >> itmp;
		fin >> stmp >> itmp >> itmp;
		fin >> stmp >> n;
		for (int i = 0; i < n + 1; i++)
		{
			getline(fin, stmp);
		}
		fin >> stmp >> n;
		//cout << stmp << "\n";
		//cout << n << "\n"; getchar();
		bzmesh.resize(n);
		for (int i = 0; i < n; i++)
		{
			fin >> itmp >> nIEN;
			bzmesh[i].IEN.resize(nIEN);
			bzmesh[i].cmat.resize(nIEN);
			for (int j = 0; j < nIEN; j++)
			{
				fin >> bzmesh[i].IEN[j];
			}
			for (int j = 0; j < nIEN; j++)
			{
				for (int k = 0; k < 16; k++)
				{
					fin >> bzmesh[i].cmat[j][k];
					bzmesh[i].pts[k][0] += bzmesh[i].cmat[j][k] * cp[bzmesh[i].IEN[j]].coor[0];
					bzmesh[i].pts[k][1] += bzmesh[i].cmat[j][k] * cp[bzmesh[i].IEN[j]].coor[1];
					bzmesh[i].pts[k][2] += bzmesh[i].cmat[j][k] * cp[bzmesh[i].IEN[j]].coor[2];
				}
			}
		}
		fin.close();
	}
	else
	{
		cerr << "Can't open " << fn << "!\n";
	}
}

void TruncatedTspline::InitializeProjectFacePoints(int proj_type, double proj_beta)
{
	//Intialization step, afterwards it will be included in the refinement step
	bzcp.clear();
	int bzloc[4] = {5, 6, 10, 9};
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].type == 4)
		{
			tmesh[i].c1 = 1; //irregular element to be c1 element
		}
	}
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].c1 == 0)
		{
			int flag(0);
			for (int j = 0; j < 4; j++)
			{
				for (uint k = 0; k < cp[tmesh[i].cnct[j]].face.size(); k++)
				{
					if (tmesh[cp[tmesh[i].cnct[j]].face[k]].c1 == 1)
					{
						flag = 1;
						break;
					}
				}
				if (flag == 1)
					break;
			}
			if (flag == 1)
				tmesh[i].c1 = 2; //transition
		}
	}
	for (uint i = 0; i < cp.size(); i++)
	{
		int flag(0);
		for (uint j = 0; j < cp[i].face.size(); j++)
		{
			if (tmesh[cp[i].face[j]].c1 != 1)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
		{
			cp[i].act = 0;
		}
		else
		{
			cp[i].act = 1;
		}
	}

	double a(4. / 9.), b(2. / 9.), c(1. / 9.);
	double a0(2. / 3.), b0(1. / 3.);
	double bzcf[4] = {a, b, c, b};
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1 && tmesh[eid].c1 == 1)
		{
			//		    array<array<double,3>,4> fcpt;
			//		    fcpt[0][0] = a * cp[tmesh[eid].cnct[0]].coor[0] + b * cp[tmesh[eid].cnct[1]].coor[0] +
			//		                 c * cp[tmesh[eid].cnct[2]].coor[0] + b * cp[tmesh[eid].cnct[3]].coor[0];
			//            fcpt[0][1] = a * cp[tmesh[eid].cnct[0]].coor[1] + b * cp[tmesh[eid].cnct[1]].coor[1] +
			//                         c * cp[tmesh[eid].cnct[2]].coor[1] + b * cp[tmesh[eid].cnct[3]].coor[1];
			//            fcpt[0][2] = a * cp[tmesh[eid].cnct[0]].coor[2] + b * cp[tmesh[eid].cnct[1]].coor[2] +
			//                         c * cp[tmesh[eid].cnct[2]].coor[2] + b * cp[tmesh[eid].cnct[3]].coor[2];
			//            fcpt[1][0] = a0 * cp[tmesh[eid].cnct[0]].coor[0] + b0 * cp[tmesh[eid].cnct[3]].coor[0];
			//            fcpt[1][1] = a0 * cp[tmesh[eid].cnct[0]].coor[1] + b0 * cp[tmesh[eid].cnct[3]].coor[1];
			//            fcpt[1][2] = a0 * cp[tmesh[eid].cnct[0]].coor[2] + b0 * cp[tmesh[eid].cnct[3]].coor[2];
			//            fcpt[2][0] = cp[tmesh[eid].cnct[0]].coor[0];
			//            fcpt[2][1] = cp[tmesh[eid].cnct[0]].coor[1];
			//            fcpt[2][2] = cp[tmesh[eid].cnct[0]].coor[2];
			//            fcpt[3][0] = a0 * cp[tmesh[eid].cnct[0]].coor[0] + b0 * cp[tmesh[eid].cnct[1]].coor[0];
			//            fcpt[3][1] = a0 * cp[tmesh[eid].cnct[0]].coor[1] + b0 * cp[tmesh[eid].cnct[1]].coor[1];
			//            fcpt[3][2] = a0 * cp[tmesh[eid].cnct[0]].coor[2] + b0 * cp[tmesh[eid].cnct[1]].coor[2];
			//            for(int i=0; i<4; i++)
			//            {
			//                bzcp.push_back(fcpt[i]);
			//                tmesh[eid].IENc1.push_back(bzcp.size() - 1);
			//            }

			for (int i = 0; i < 4; i++)
			{
				array<double, 3> fcpt = {0., 0., 0.};
				for (int j = 0; j < 4; j++)
				{
					fcpt[0] += bzcf[j] * cp[tmesh[eid].cnct[(i + j) % 4]].coor[0];
					fcpt[1] += bzcf[j] * cp[tmesh[eid].cnct[(i + j) % 4]].coor[1];
					fcpt[2] += bzcf[j] * cp[tmesh[eid].cnct[(i + j) % 4]].coor[2];
				}
				bzcp.push_back(fcpt);
				tmesh[eid].IENc1.push_back(bzcp.size() - 1); //c1mat will be set in FindSplineC1_Irr()
			}
		}
	}

	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1)
		{
			if (tmesh[eid].c1 == 1)
			{
				FindSplineC1_Irr(eid);
			}
			//else if (tmesh[eid].c1 == 2)
			//{
			//	FindSplineC1_Trs(eid);
			//}
			//else if (tmesh[eid].c1 == 11)
			//{
			//	//FindSplineC1_Trs(eid);
			//}
		}
	}

	//for (uint eid = 0; eid < tmesh.size(); eid++)
	//{
	//	if (tmesh[eid].act == 1)
	//	{
	//		if (tmesh[eid].c1 == 2)
	//		{
	//			FindSplineC1_Trs(eid);
	//		}
	//	}
	//}

	SetTrunMat();

	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].c1 == 1)
		{
			if (tmesh[i].type == 4)
			{
				vector<BezierElement> bzmesh;
				BezierUnit_DPatch_Irr22(i, bzmesh); //bzmesh.size()==4
				tmesh[i].bzpts.resize(16 * 4);
				for (int j = 0; j < 4; j++)
				{
					for (int k = 0; k < 16; k++)
					{
						tmesh[i].bzpts[j * 16 + k][0] = bzmesh[j].pts[k][0];
						tmesh[i].bzpts[j * 16 + k][1] = bzmesh[j].pts[k][1];
						tmesh[i].bzpts[j * 16 + k][2] = bzmesh[j].pts[k][2];
					}
				}
			}
			else
			{
				vector<BezierElement> bzmesh;
				BezierUnit_DPatch_Trs(i, bzmesh); //bzmesh.size()==1
				tmesh[i].bzpts.resize(16);
				for (int k = 0; k < 16; k++)
				{
					tmesh[i].bzpts[k][0] = bzmesh[0].pts[k][0];
					tmesh[i].bzpts[k][1] = bzmesh[0].pts[k][1];
					tmesh[i].bzpts[k][2] = bzmesh[0].pts[k][2];
				}
			}
		}
	}
}

void TruncatedTspline::InitializeProjectFacePoints_Scale(int proj_type, double proj_beta)
{
	//Intialization step, afterwards it will be included in the refinement step
	bzcp.clear();
	wc1.clear();
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].type == 4)
		{
			tmesh[i].c1 = 1; //irregular element to be c1 element
		}
	}
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].c1 == 0)
		{
			int flag(0);
			for (int j = 0; j < 4; j++)
			{
				for (uint k = 0; k < cp[tmesh[i].cnct[j]].face.size(); k++)
				{
					if (tmesh[cp[tmesh[i].cnct[j]].face[k]].c1 == 1)
					{
						flag = 1;
						break;
					}
				}
				if (flag == 1)
					break;
			}
			if (flag == 1)
				tmesh[i].c1 = 2; //transition
		}
	}
	for (uint i = 0; i < cp.size(); i++)
	{
		int flag(0);
		for (uint j = 0; j < cp[i].face.size(); j++)
		{
			if (tmesh[cp[i].face[j]].c1 != 1)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
		{
			cp[i].act = 0;
		}
		else
		{
			cp[i].act = 1;
		}
	}

	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1 && tmesh[eid].c1 == 1)
		{
			double ki[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
			ki[0][1] = tmedge[tmesh[eid].edge[0]].len;
			ki[1][1] = tmedge[tmesh[eid].edge[3]].len;
			double kitmp[4];
			for (int i = 0; i < 4; i++)
			{
				int ednb(tmedge[tmesh[eid].edge[i]].face[0]);
				if (ednb == eid)
					ednb = tmedge[tmesh[eid].edge[i]].face[1];
				int *it = find(tmesh[ednb].edge, tmesh[ednb].edge + 4, tmesh[eid].edge[i]);
				int loc(it - tmesh[ednb].edge);
				kitmp[i] = tmedge[tmesh[ednb].edge[(loc + 1) % 4]].len;
			}
			ki[0][0] = kitmp[3];
			ki[0][2] = kitmp[1];
			ki[1][0] = kitmp[0];
			ki[1][2] = kitmp[2];
			double kilen[2] = {ki[0][0] + ki[0][1] + ki[0][2], ki[1][0] + ki[1][1] + ki[1][2]};
			double au[2] = {(ki[0][1] + ki[0][2]) / kilen[0], ki[0][0] / kilen[0]};
			double av[2] = {(ki[1][1] + ki[1][2]) / kilen[1], ki[1][0] / kilen[1]};
			double c[4] = {au[0] * av[0], au[1] * av[0], au[1] * av[1], au[0] * av[1]};
			array<array<double, 3>, 4> fcpt;
			fcpt[0][0] = c[0] * cp[tmesh[eid].cnct[0]].coor[0] + c[1] * cp[tmesh[eid].cnct[1]].coor[0] +
						 c[2] * cp[tmesh[eid].cnct[2]].coor[0] + c[3] * cp[tmesh[eid].cnct[3]].coor[0];
			fcpt[0][1] = c[0] * cp[tmesh[eid].cnct[0]].coor[1] + c[1] * cp[tmesh[eid].cnct[1]].coor[1] +
						 c[2] * cp[tmesh[eid].cnct[2]].coor[1] + c[3] * cp[tmesh[eid].cnct[3]].coor[1];
			fcpt[0][2] = c[0] * cp[tmesh[eid].cnct[0]].coor[2] + c[1] * cp[tmesh[eid].cnct[1]].coor[2] +
						 c[2] * cp[tmesh[eid].cnct[2]].coor[2] + c[3] * cp[tmesh[eid].cnct[3]].coor[2];
			fcpt[1][0] = av[0] * cp[tmesh[eid].cnct[0]].coor[0] + av[1] * cp[tmesh[eid].cnct[3]].coor[0];
			fcpt[1][1] = av[0] * cp[tmesh[eid].cnct[0]].coor[1] + av[1] * cp[tmesh[eid].cnct[3]].coor[1];
			fcpt[1][2] = av[0] * cp[tmesh[eid].cnct[0]].coor[2] + av[1] * cp[tmesh[eid].cnct[3]].coor[2];
			fcpt[2][0] = cp[tmesh[eid].cnct[0]].coor[0];
			fcpt[2][1] = cp[tmesh[eid].cnct[0]].coor[1];
			fcpt[2][2] = cp[tmesh[eid].cnct[0]].coor[2];
			fcpt[3][0] = au[0] * cp[tmesh[eid].cnct[0]].coor[0] + au[1] * cp[tmesh[eid].cnct[1]].coor[0];
			fcpt[3][1] = au[0] * cp[tmesh[eid].cnct[0]].coor[1] + au[1] * cp[tmesh[eid].cnct[1]].coor[1];
			fcpt[3][2] = au[0] * cp[tmesh[eid].cnct[0]].coor[2] + au[1] * cp[tmesh[eid].cnct[1]].coor[2];
			double wctmp[4] = {1., ki[0][2] / kilen[0], ki[0][2] * ki[1][2] / (kilen[0] * kilen[1]), ki[1][2] / kilen[1]};
			for (int i = 0; i < 4; i++)
			{
				bzcp.push_back(fcpt[i]);
				wc1.push_back(wctmp[i]);
				tmesh[eid].IENc1.push_back(bzcp.size() - 1);
			}
		}
	}

	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1)
		{
			if (tmesh[eid].c1 == 1)
			{
				FindSplineC1_Irr(eid);
			}
			else if (tmesh[eid].c1 == 2)
			{
				FindSplineC1_Scale_Trs(eid);
			}
		}
	}

	FindSplineC2_Scale();

	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].c1 == 1)
		{
			if (tmesh[i].type == 4)
			{
				vector<BezierElement> bzmesh;
				BezierUnit_DPatch_Irr22_Scale(i, bzmesh); //bzmesh.size()==4
				tmesh[i].bzpts.resize(16 * 4);
				for (int j = 0; j < 4; j++)
				{
					for (int k = 0; k < 16; k++)
					{
						tmesh[i].bzpts[j * 16 + k][0] = bzmesh[j].pts[k][0];
						tmesh[i].bzpts[j * 16 + k][1] = bzmesh[j].pts[k][1];
						tmesh[i].bzpts[j * 16 + k][2] = bzmesh[j].pts[k][2];
					}
				}
			}
			else
			{
				vector<BezierElement> bzmesh;
				BezierUnit_DPatch_Trs_Scale(i, bzmesh); //bzmesh.size()==1
				tmesh[i].bzpts.resize(16);
				for (int k = 0; k < 16; k++)
				{
					tmesh[i].bzpts[k][0] = bzmesh[0].pts[k][0];
					tmesh[i].bzpts[k][1] = bzmesh[0].pts[k][1];
					tmesh[i].bzpts[k][2] = bzmesh[0].pts[k][2];
				}
			}
		}
	}
}

void TruncatedTspline::AddProjectFacePoints(int proj_type, double proj_beta)
{
	//bzcp.clear();
	int bzloc[4] = {5, 6, 10, 9};
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].c1 == 0)
		{
			int flag(0);
			for (int j = 0; j < 4; j++)
			{
				for (uint k = 0; k < cp[tmesh[i].cnct[j]].face.size(); k++)
				{
					if (tmesh[cp[tmesh[i].cnct[j]].face[k]].c1 == 1)
					{
						flag = 1;
						break;
					}
				}
				if (flag == 1)
					break;
			}
			if (flag == 1)
				tmesh[i].c1 = 2; //transition
		}
	}
	for (uint i = 0; i < cp.size(); i++)
	{
		int flag(0);
		for (uint j = 0; j < cp[i].face.size(); j++)
		{
			if (tmesh[cp[i].face[j]].c1 != 1)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
		{
			cp[i].act = 0;
		}
		else
		{
			cp[i].act = 1;
		}
	}

	//for (uint eid = 0; eid < tmesh.size(); eid++)
	//{
	//	if (tmesh[eid].act == 1 && tmesh[eid].c1 == 1)
	//	{
	//		//int nv(cp[tmesh[eid].cnct[0]].face.size());
	//		//vector<vector<double>> bemat;
	//		//SetBezierMat_DPatch(proj_type, nv, proj_beta, bemat);
	//		//for (int i = 0; i < 4; i++)
	//		//{
	//		//	array<double, 3> fcpt = { 0.,0.,0. };
	//		//	for (uint j = 0; j < 2 * nv + 1; j++)
	//		//	{
	//		//		fcpt[0] += bemat[j][bzloc[i]] * cp[tmesh[eid].IEN[j]].coor[0];
	//		//		fcpt[1] += bemat[j][bzloc[i]] * cp[tmesh[eid].IEN[j]].coor[1];
	//		//		fcpt[2] += bemat[j][bzloc[i]] * cp[tmesh[eid].IEN[j]].coor[2];
	//		//	}
	//		//	bzcp.push_back(fcpt);
	//		//	tmesh[eid].IENc1.push_back(bzcp.size() - 1);//c1mat will be set in FindSplineC1_Irr()
	//		//}
	//	}
	//	//else if (tmesh[eid].act == 1 && tmesh[eid].c1 == 11)//irregular ones after refinement
	//	//{
	//	//	//vector<vector<double>> bemat;
	//	//	//for (int i = 0; i < 4; i++)
	//	//	//{
	//	//	//	array<double, 3> fcpt = { 0.,0.,0. };
	//	//	//	for (uint j = 0; j < 2 * nv + 1; j++)
	//	//	//	{
	//	//	//		fcpt[0] += bemat[j][bzloc[i]] * cp[tmesh[eid].IEN[j]].coor[0];
	//	//	//		fcpt[1] += bemat[j][bzloc[i]] * cp[tmesh[eid].IEN[j]].coor[1];
	//	//	//		fcpt[2] += bemat[j][bzloc[i]] * cp[tmesh[eid].IEN[j]].coor[2];
	//	//	//	}
	//	//	//	bzcp.push_back(fcpt);
	//	//	//	tmesh[eid].IENc1.push_back(bzcp.size() - 1);
	//	//	//}
	//	//}
	//}

	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1)
		{
			if (tmesh[eid].c1 == 1 && tmesh[eid].type == 4)
			{
				FindSplineC1_Irr(eid);
			}
			else if (tmesh[eid].c1 == 1 && tmesh[eid].type != 4)
			{
				FindSplineC1_C1Element(eid);
			}
			//else if (tmesh[eid].c1 == 2)
			//{
			//	FindSplineC1_Trs(eid);
			//}
			//else if (tmesh[eid].c1 == 11)
			//{
			//	//FindSplineC1_Trs(eid);
			//}
		}
	}

	//for (uint eid = 0; eid < tmesh.size(); eid++)
	//{
	//	if (tmesh[eid].act == 1)
	//	{
	//		if (tmesh[eid].c1 == 2)
	//		{
	//			FindSplineC1_Trs(eid);
	//		}
	//	}
	//}
}

void TruncatedTspline::FindSplineC1_Irr(int eid)
{
	int nv(cp[tmesh[eid].cnct[0]].face.size());
	int ednb[4], cnnb[3];
	for (int i = 0; i < 4; i++)
	{
		ednb[i] = tmedge[tmesh[eid].edge[i]].face[0];
		if (ednb[i] == eid)
			ednb[i] = tmedge[tmesh[eid].edge[i]].face[1];
	}
	for (int i = 1; i < 4; i++)
	{
		cnnb[i - 1] = -1;
		for (uint j = 0; j < cp[tmesh[eid].cnct[i]].face.size(); j++)
		{
			int fcid(cp[tmesh[eid].cnct[i]].face[j]);
			if (fcid != eid && fcid != ednb[i] && fcid != ednb[(i + 3) % 4])
			{
				cnnb[i - 1] = fcid;
				break;
			}
		}
	}

	int enext(ednb[3]);
	for (int i = 1; i < nv; i++)
	{
		tmesh[eid].IENc1.push_back(tmesh[enext].IENc1[0]);
		tmesh[eid].IENc1.push_back(tmesh[enext].IENc1[1]);
		tmesh[eid].IENc1.push_back(tmesh[enext].IENc1[3]);
		int enow(enext);
		enext = tmedge[tmesh[enow].edge[3]].face[0];
		if (enext == enow)
			enext = tmedge[tmesh[enow].edge[3]].face[1];
	}

	SetFaceBezierMat_DPatch22(SMOOTH_TYPE, nv, SMOOTH_BETA, tmesh[eid].c1mat22);

	//tmesh[eid].IENc1.push_back(tmesh[ednb[0]].IENc1[0]);//IENc1[4]
	//tmesh[eid].IENc1.push_back(tmesh[ednb[0]].IENc1[3]);//IENc1[5]
	//tmesh[eid].IENc1.push_back(tmesh[ednb[3]].IENc1[0]);//IENc1[6]
	//tmesh[eid].IENc1.push_back(tmesh[ednb[3]].IENc1[1]);//IENc1[7]
	//if (nv >= 5)
	//{
	//	int enext(tmedge[tmesh[ednb[3]].edge[3]].face[0]);
	//	if (enext == ednb[3]) enext = tmedge[tmesh[ednb[3]].edge[3]].face[1];
	//	for (int i = 0; i < nv - 3; i++)
	//	{
	//		tmesh[eid].IENc1.push_back(tmesh[enext].IENc1[0]);
	//		int enow(enext);
	//		enext = tmedge[tmesh[enow].edge[3]].face[0];
	//		if (enext == enow) enext = tmedge[tmesh[enow].edge[3]].face[1];
	//	}
	//}
	//SetFaceBezierMat_DPatch(SMOOTH_TYPE, nv, SMOOTH_BETA, tmesh[eid].c1mat);

	//other c1 functions from edge or corner neighbors, whcih applies to the case after refinement
	MatrixXd smat = MatrixXd::Zero(16, 49);
	double smat1d[7][4] = {{1., 0., 0., 0.}, {.5, .5, 0., 0.}, {.25, .5, .25, 0.}, {.125, .375, .375, .125}, {0., .25, .5, .25}, {0., 0., .5, .5}, {0., 0., 0., 1.}};
	int loc(0);
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			int loc0(0);
			for (int i0 = 0; i0 < 4; i0++)
			{
				for (int j0 = 0; j0 < 4; j0++)
				{
					smat(loc0, loc) = smat1d[j][j0] * smat1d[i][i0];
					loc0++;
				}
			}
			loc++;
		}
	}
	int bzcnct[4][16] = {{0, 1, 2, 3, 7, 8, 9, 10, 14, 15, 16, 17, 21, 22, 23, 24}, {3, 4, 5, 6, 10, 11, 12, 13, 17, 18, 19, 20, 24, 25, 26, 27}, {21, 22, 23, 24, 28, 29, 30, 31, 35, 36, 37, 38, 42, 43, 44, 45}, {24, 25, 26, 27, 31, 32, 33, 34, 38, 39, 40, 41, 45, 46, 47, 48}};

	int ed_eid_bz[4][2][2] = {{{1, 0}, {2, 3}}, {{7, 3}, {11, 15}}, {{14, 15}, {13, 12}}, {{8, 12}, {4, 0}}};
	int ed_nb_bz[4][2][2] = {{{4, 0}, {8, 12}}, {{2, 3}, {1, 0}}, {{11, 15}, {7, 3}}, {{13, 12}, {14, 15}}};
	for (int i = 1; i < 3; i++)
	{
		if (tmesh[ednb[i]].c1 == 1 && tmesh[ednb[i]].type != 4)
		{
			int *it1 = find(tmesh[ednb[i]].cnct, tmesh[ednb[i]].cnct + 4, tmesh[eid].cnct[i]);
			int iloc = it1 - tmesh[ednb[i]].cnct;
			vector<double> ctmp1(16, 0.);
			int j = 0;
			tmesh[eid].IENc1.push_back(tmesh[ednb[i]].IENc1[iloc]);
			ctmp1[ed_eid_bz[i][j][0]] = 0.5;
			ctmp1[ed_eid_bz[i][j][1]] = 0.25;
			vector<double> ctmp2(16, 0.);
			iloc = (iloc + 3) % 4;
			j = 1;
			tmesh[eid].IENc1.push_back(tmesh[ednb[i]].IENc1[iloc]);
			ctmp2[ed_eid_bz[i][j][0]] = 0.5;
			ctmp2[ed_eid_bz[i][j][1]] = 0.25;
			vector<double> c1(49, 0.), c2(49, 0.);
			for (j = 0; j < 16; j++)
			{
				for (int k = 0; k < 49; k++)
				{
					c1[k] += ctmp1[j] * smat(j, k);
					c2[k] += ctmp2[j] * smat(j, k);
				}
			}
			for (j = 0; j < 4; j++)
			{
				for (int k = 0; k < 16; k++)
				{
					ctmp1[k] = c1[bzcnct[j][k]];
					ctmp2[k] = c2[bzcnct[j][k]];
				}
				tmesh[eid].c1mat22[j].push_back(ctmp1);
				tmesh[eid].c1mat22[j].push_back(ctmp2);
			}
		}
	}

	int cn_eid_bz[3] = {3, 15, 12};
	int cn_nb_bz[3] = {3, 15, 12};
	for (int i = 0; i < 3; i++)
	{
		if (tmesh[cnnb[i]].c1 == 1 && tmesh[cnnb[i]].type != 4)
		{
			int *it1 = find(tmesh[cnnb[i]].cnct, tmesh[cnnb[i]].cnct + 4, tmesh[eid].cnct[i + 1]);
			int iloc = it1 - tmesh[cnnb[i]].cnct;
			vector<double> ctmp(16, 0.);
			tmesh[eid].IENc1.push_back(tmesh[cnnb[i]].IENc1[iloc]);
			ctmp[cn_eid_bz[i]] = 0.25;
			vector<double> c1(49, 0.);
			for (int j = 0; j < 16; j++)
			{
				for (int k = 0; k < 49; k++)
				{
					c1[k] += ctmp[j] * smat(j, k);
				}
			}
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 16; k++)
				{
					ctmp[k] = c1[bzcnct[j][k]];
				}
				tmesh[eid].c1mat22[j].push_back(ctmp);
			}
		}
	}
}

void TruncatedTspline::FindSplineC1_C1Element(int eid)
{
	//c1==1 && type!=4
	//	int nv(cp[tmesh[eid].cnct[0]].face.size());
	int ednb[4], cnnb[4];
	for (int i = 0; i < 4; i++)
	{
		ednb[i] = tmedge[tmesh[eid].edge[i]].face[0];
		if (ednb[i] == eid)
			ednb[i] = tmedge[tmesh[eid].edge[i]].face[1];
	}
	for (int i = 0; i < 4; i++)
	{
		cnnb[i] = -1;
		for (uint j = 0; j < cp[tmesh[eid].cnct[i]].face.size(); j++)
		{
			int fcid(cp[tmesh[eid].cnct[i]].face[j]);
			if (fcid != eid && fcid != ednb[i] && fcid != ednb[(i + 3) % 4])
			{
				cnnb[i] = fcid;
				break;
			}
		}
	}

	//first add the ordinates of the four face points
	int bzfc[4][4] = {{5, 1, 4, 0}, {6, 2, 7, 3}, {10, 11, 14, 15}, {9, 8, 13, 12}};
	double cfc[4] = {1., 0.5, 0.5, 0.25};
	for (int i = 0; i < 4; i++)
	{
		vector<double> ctmp(16, 0.);
		for (int j = 0; j < 4; j++)
		{
			ctmp[bzfc[i][j]] = cfc[j];
		}
		tmesh[eid].c1mat.push_back(ctmp);
	}

	//other c1 functions from edge or corner neighbors, whcih applies to the case after refinement
	/*int ed_eid_bz[4][2][2] = { { { 1,0 },{ 2,3 } },{ { 7,3 },{ 11,15 } },{ { 14,15 },{ 13,12 } },{ { 8,12 },{ 4,0 } } };
	int ed_nb_bz[4][2][2] = { { { 4,0 },{ 8,12 } },{ { 2,3 },{ 1,0 } },{ { 11,15 },{ 7,3 } },{ { 13,12 },{ 14,15 } } };
	for (int i = 0; i < 4; i++)
	{
		if (tmesh[ednb[i]].c1 == 1)
		{
			int* it1 = find(tmesh[ednb[i]].cnct, tmesh[ednb[i]].cnct + 4, tmesh[eid].cnct[i]);
			int iloc = it1 - tmesh[ednb[i]].cnct;
			vector<double> ctmp1(16, 0.);
			int j = 0;
			tmesh[eid].IENc1.push_back(tmesh[ednb[i]].IENc1[iloc]);
			ctmp1[ed_eid_bz[i][j][0]] = 0.5;
			ctmp1[ed_eid_bz[i][j][1]] = 0.25;
			tmesh[eid].c1mat.push_back(ctmp1);
			vector<double> ctmp2(16, 0.);
			iloc = (iloc + 3) % 4;
			j = 1;
			tmesh[eid].IENc1.push_back(tmesh[ednb[i]].IENc1[iloc]);
			ctmp2[ed_eid_bz[i][j][0]] = 0.5;
			ctmp2[ed_eid_bz[i][j][1]] = 0.25;
			tmesh[eid].c1mat.push_back(ctmp2);
		}
	}

	int cn_eid_bz[4] = { 0,3,15,12 };
	int cn_nb_bz[4] = { 0,3,15,12 };
	for (int i = 0; i < 4; i++)
	{
		if (tmesh[cnnb[i]].c1 == 1)
		{
			int* it1 = find(tmesh[cnnb[i]].cnct, tmesh[cnnb[i]].cnct + 4, tmesh[eid].cnct[i]);
			int iloc = it1 - tmesh[cnnb[i]].cnct;
			vector<double> ctmp(16, 0.);
			tmesh[eid].IENc1.push_back(tmesh[cnnb[i]].IENc1[iloc]);
			ctmp[cn_eid_bz[i]] = 0.25;
			tmesh[eid].c1mat.push_back(ctmp);
		}
	}*/
}

void TruncatedTspline::FindSplineC1_Trs(int eid)
{
	int ednb[4][2], cnnb[4][2];
	for (int i = 0; i < 4; i++)
	{
		ednb[i][0] = tmedge[tmesh[eid].edge[i]].face[0];
		if (ednb[i][0] == eid)
			ednb[i][0] = tmedge[tmesh[eid].edge[i]].face[1];
		int *it = find(tmesh[ednb[i][0]].cnct, tmesh[ednb[i][0]].cnct + 4, tmesh[eid].cnct[i]);
		ednb[i][1] = it - tmesh[ednb[i][0]].cnct;
	}
	for (int i = 0; i < 4; i++)
	{
		int pos(-1);
		for (uint j = 0; j < cp[tmesh[eid].cnct[i]].face.size(); j++)
		{
			int fcid(cp[tmesh[eid].cnct[i]].face[j]);
			if (fcid != eid && fcid != ednb[i][0] && fcid != ednb[(i + 3) % 4][0])
			{
				pos = fcid;
				break;
			}
		}
		if (pos != -1)
		{
			cnnb[i][0] = pos;
			int *it = find(tmesh[cnnb[i][0]].cnct, tmesh[cnnb[i][0]].cnct + 4, tmesh[eid].cnct[i]);
			cnnb[i][1] = it - tmesh[cnnb[i][0]].cnct;
		}
		else
		{
			cerr << "Error! Can't find corner neighbor!\n";
			getchar();
		}
	}

	int ed_eid_bz[4][2][2] = {{{1, 0}, {2, 3}}, {{7, 3}, {11, 15}}, {{14, 15}, {13, 12}}, {{8, 12}, {4, 0}}};
	int ed_nb_bz[4][2][2] = {{{4, 0}, {8, 12}}, {{2, 3}, {1, 0}}, {{11, 15}, {7, 3}}, {{13, 12}, {14, 15}}};
	for (int i = 0; i < 4; i++)
	{
		if (tmesh[ednb[i][0]].c1 == 1 || tmesh[ednb[i][0]].c1 == 11)
		{
			vector<double> ctmp1(16, 0.);
			int iloc(ednb[i][1]);
			int j = 0;
			tmesh[eid].IENc1.push_back(tmesh[ednb[i][0]].IENc1[iloc]);
			ctmp1[ed_eid_bz[i][j][0]] = tmesh[ednb[i][0]].c1mat[iloc][ed_nb_bz[ednb[i][1]][j][0]];
			ctmp1[ed_eid_bz[i][j][1]] = tmesh[ednb[i][0]].c1mat[iloc][ed_nb_bz[ednb[i][1]][j][1]];
			tmesh[eid].c1mat.push_back(ctmp1);
			vector<double> ctmp2(16, 0.);
			iloc = (ednb[i][1] + 3) % 4;
			j = 1;
			tmesh[eid].IENc1.push_back(tmesh[ednb[i][0]].IENc1[iloc]);
			ctmp2[ed_eid_bz[i][j][0]] = tmesh[ednb[i][0]].c1mat[iloc][ed_nb_bz[ednb[i][1]][j][0]];
			ctmp2[ed_eid_bz[i][j][1]] = tmesh[ednb[i][0]].c1mat[iloc][ed_nb_bz[ednb[i][1]][j][1]];
			tmesh[eid].c1mat.push_back(ctmp2);
		}
	}

	int cn_eid_bz[4] = {0, 3, 15, 12};
	int cn_nb_bz[4] = {0, 3, 15, 12};
	for (int i = 0; i < 4; i++)
	{
		if (tmesh[cnnb[i][0]].c1 == 1 || tmesh[cnnb[i][0]].c1 == 11)
		{
			vector<double> ctmp(16, 0.);
			int iloc(cnnb[i][1]);
			tmesh[eid].IENc1.push_back(tmesh[cnnb[i][0]].IENc1[iloc]);
			ctmp[cn_eid_bz[i]] = tmesh[cnnb[i][0]].c1mat[iloc][cn_nb_bz[cnnb[i][1]]];
			tmesh[eid].c1mat.push_back(ctmp);
		}
	}
}

void TruncatedTspline::SetFaceBezierMat_DPatch(int type, int N, double beta, vector<vector<double>> &bemat)
{
	//problamtic
	bemat.clear();
	int nc1(N + 5);
	vector<vector<double>> proj_mat_all, proj_mat(nc1, vector<double>(nc1, 0.));
	ProjectMatrix(type, N, beta, proj_mat_all); //3N by 3N, from arbitray to projected
	int indx0[8] = {0, N, -1, 2 * N, N - 1, 3 * N - 1, 1, N + 1};
	vector<int> indx(indx0, indx0 + 8);
	for (int i = 2; i < N - 1; i++)
	{
		indx.push_back(i);
	}
	for (int i = 0; i < nc1; i++)
	{
		for (int j = 0; j < nc1; j++)
		{
			if (indx[i] != -1 && indx[j] != -1)
			{
				proj_mat[i][j] = proj_mat_all[indx[i]][indx[j]];
			}
		}
	}
	proj_mat[2][2] = 1.;

	//cout << "proj_mat\n";
	//for (uint i = 0; i < proj_mat.size(); i++)
	//{
	//	cout << i << "\n";
	//	double sum(0.);
	//	for (uint j = 0; j < proj_mat[i].size(); j++)
	//	{
	//		cout << proj_mat[i][j] << "(" << j << ") ";
	//		sum += proj_mat[i][j];
	//	}
	//	cout << "\n";
	//	cout << "sum: " << sum << "\n";
	//}
	////for (uint i = 0; i < proj_mat_all.size(); i++)
	////{
	////	cout << proj_mat_all[0][i] << " " << proj_mat_all[1][i] << " " << proj_mat_all[2][i] << "\n";
	////}
	//getchar();

	int loc4[4][4] = {{5, 1, 4, 0}, {6, 2, 7, 3}, {10, 11, 14, 15}, {9, 8, 13, 12}};
	double coef4[5] = {1., .5, .5, .25, 1. / double(N)};
	int loc2[4][2] = {{1, 0}, {2, 3}, {4, 0}, {8, 12}};
	double coef2[3] = {.5, .25, 1. / double(N)};
	vector<vector<double>> avg(nc1, vector<double>(16, 0.));
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			avg[i][loc4[i][j]] = coef4[j];
		}
	}
	avg[0][0] = coef4[4];
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			avg[i + 4][loc2[i][j]] = coef2[j];
		}
	}
	avg[4][0] = coef2[2];
	avg[6][0] = coef2[2];
	for (int i = 0; i < N - 3; i++)
	{
		avg[8 + i][0] = coef2[2];
	}

	bemat.resize(nc1, vector<double>(16, 0.));
	for (int i = 0; i < nc1; i++) //before project
	{
		for (int j = 0; j < nc1; j++) //after project
		{
			for (int k = 0; k < 16; k++)
			{
				bemat[i][k] += proj_mat[j][i] * avg[j][k];
			}
		}
	}
}

void TruncatedTspline::SetFaceBezierMat_DPatch22(int type, int N, double beta, vector<vector<vector<double>>> &c1mat)
{
	c1mat.clear();
	c1mat.resize(4);
	MatrixXd amat = MatrixXd::Zero(N + 5, 16);
	for (int i = 0; i < N; i++)
		amat(i, 0) = 1. / double(N);
	int fc1[4] = {0, N, N + 1, N + 4};
	int bz1[4][4] = {{5, 1, 4, 0}, {6, 2, 7, 3}, {9, 8, 13, 12}, {10, 11, 14, 15}};
	double coef1[4] = {1., 0.5, 0.5, 0.25};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (bz1[i][j] != 0)
				amat(fc1[i], bz1[i][j]) = coef1[j];
		}
	}
	int fc2[4] = {N - 1, N + 3, 1, N + 2};
	int bz2[4][2] = {{1, 0}, {2, 3}, {4, 0}, {8, 12}};
	double coef2[2] = {0.5, 0.25};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			if (bz2[i][j] != 0)
				amat(fc2[i], bz2[i][j]) = coef2[j];
		}
	}

	MatrixXd smat = MatrixXd::Zero(16, 49);
	double smat1d[7][4] = {{1., 0., 0., 0.}, {.5, .5, 0., 0.}, {.25, .5, .25, 0.}, {.125, .375, .375, .125}, {0., .25, .5, .25}, {0., 0., .5, .5}, {0., 0., 0., 1.}};
	int loc(0);
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			int loc0(0);
			for (int i0 = 0; i0 < 4; i0++)
			{
				for (int j0 = 0; j0 < 4; j0++)
				{
					smat(loc0, loc) = smat1d[j][j0] * smat1d[i][i0];
					loc0++;
				}
			}
			loc++;
		}
	}
	MatrixXd me = amat * smat; //(N+5) by 49

	MatrixXd m1 = MatrixXd::Zero(3 * N, 3 * N); //(i,j) i: macro c1 functions; and j: micro c1
	vector<vector<double>> pmat;
	ProjectMatrix(type, N, beta, pmat); //3N by 3N, from arbitray to projected
	for (int eid = 0; eid < N; eid++)	//counterclock-wise
	{
		int iglb;
		for (int i = 0; i < N; i++) //first N macro c1
		{
			iglb = (eid + i) % N;
			m1(iglb, eid) = me(i, 8);
			m1(iglb, N + eid) = me(i, 9);
			m1(iglb, 2 * N + eid) = me(i, 15);
		}
		iglb = N + eid;
		m1(iglb, eid) = me(N, 8);
		m1(iglb, N + eid) = me(N, 9);
		m1(iglb, 2 * N + eid) = me(N, 15);
		iglb = 2 * N + eid;
		m1(iglb, eid) = me(N + 1, 8);
		m1(iglb, N + eid) = me(N + 1, 9);
		m1(iglb, 2 * N + eid) = me(N + 1, 15);
		iglb = N + (1 + eid) % N;
		m1(iglb, eid) = me(N + 2, 8);
		m1(iglb, N + eid) = me(N + 2, 9);
		m1(iglb, 2 * N + eid) = me(N + 2, 15);
		iglb = 2 * N + (N - 1 + eid) % N;
		m1(iglb, eid) = me(N + 3, 8);
		m1(iglb, N + eid) = me(N + 3, 9);
		m1(iglb, 2 * N + eid) = me(N + 3, 15);
	}
	MatrixXd msmth = MatrixXd::Zero(3 * N, 3 * N);
	for (int i = 0; i < 3 * N; i++) //macro c1 functions
	{
		for (int j = 0; j < 3 * N; j++) //before projection
		{
			for (int k = 0; k < 3 * N; k++) //after projection
			{
				msmth(i, k) += m1(i, j) * pmat[k][j];
			}
		}
	}
	MatrixXd mfinal = MatrixXd::Zero(3 * N + 1, 49);
	vector<int> indx(N + 5);
	for (int i = 0; i < N; i++)
		indx[i] = i;
	indx[N] = N;
	indx[N + 1] = 2 * N;
	indx[N + 2] = N + 1;
	indx[N + 3] = 3 * N - 1;
	indx[N + 4] = 3 * N;
	for (int i = 0; i < N + 5; i++)
	{
		for (int j = 0; j < 49; j++)
		{
			mfinal(indx[i], j) = me(i, j);
		}
	}
	for (int i = 0; i < 3 * N; i++)
	{
		mfinal(i, 8) = msmth(i, 0);
		mfinal(i, 9) = msmth(i, N);
		mfinal(i, 15) = msmth(i, 2 * N);
		mfinal(i, 10) = (mfinal(i, 9) + mfinal(i, 11)) / 2.;
		mfinal(i, 22) = (mfinal(i, 15) + mfinal(i, 29)) / 2.;
		mfinal(i, 1) = (mfinal(i, 8) + msmth(i, N - 1)) / 2.;
		mfinal(i, 2) = (mfinal(i, 9) + msmth(i, 3 * N - 1)) / 2.;
		mfinal(i, 7) = (mfinal(i, 8) + msmth(i, 1)) / 2.;
		mfinal(i, 14) = (mfinal(i, 15) + msmth(i, N + 1)) / 2.;
		int iloc = (i / N) * N + (i + 1) % N;
		//cout << "i/iloc: " << i << "/" << iloc << "\n"; getchar();
		mfinal(i, 3) = (mfinal(i, 9) + mfinal(i, 11) + msmth(i, 3 * N - 1) + mfinal(iloc, 29)) / 4.;
		iloc = (i / N) * N + (i + N - 1) % N;
		mfinal(i, 21) = (mfinal(i, 15) + mfinal(i, 29) + msmth(i, N + 1) + mfinal(iloc, 11)) / 4.;
		mfinal(i, 0) = 0.;
		for (int j = 0; j < N; j++)
		{
			mfinal(i, 0) += msmth(i, j);
		}
		mfinal(i, 0) /= double(N);
	}

	//cout << "mfinal:\n";
	//for (int i = 0; i < 3*N+1; i++)
	//{
	//	if (i == 3 * N)
	//	{
	//		cout << "c1 ID: " << i << "\n";
	//		for (int j = 0; j < 49; j++)
	//		{
	//			cout << mfinal(i, j) << "(" << j << ") ";
	//		}
	//		cout << "\n";
	//		getchar();
	//	}
	//}
	//getchar();

	vector<int> indx1(3 * N + 1);
	indx1[0] = 0;
	indx1[N] = 1;
	indx1[2 * N] = 3;
	indx1[3 * N] = 2;
	for (int i = 1; i < N; i++)
	{
		indx1[i] = 3 * i + 1;
		indx1[N + i] = 3 * i + 2;
		indx1[2 * N + i] = 3 * i + 3;
	}

	//for (uint i = 0; i < indx1.size(); i++)
	//{
	//	cout << "proj/IENc1: " << i << "/"<< indx1[i] << "\n";
	//}
	//getchar();

	int bzcnct[4][16] = {{0, 1, 2, 3, 7, 8, 9, 10, 14, 15, 16, 17, 21, 22, 23, 24}, {3, 4, 5, 6, 10, 11, 12, 13, 17, 18, 19, 20, 24, 25, 26, 27}, {21, 22, 23, 24, 28, 29, 30, 31, 35, 36, 37, 38, 42, 43, 44, 45}, {24, 25, 26, 27, 31, 32, 33, 34, 38, 39, 40, 41, 45, 46, 47, 48}};
	for (int i = 0; i < 4; i++)
	{
		c1mat[i].resize(3 * N + 1, vector<double>(16, 0.));
		for (int j = 0; j < 3 * N + 1; j++)
		{
			for (int k = 0; k < 16; k++)
			{
				c1mat[i][indx1[j]][k] = mfinal(j, bzcnct[i][k]);
			}
		}
	}
}

void TruncatedTspline::SetTrunMat()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1)
		{
			if (tmesh[eid].c1 == 1 && tmesh[eid].type == 4)
			{
				//SetTrunMat_Irr(eid);
				SetTrunMat_Irr22(eid);
			}
			else if (tmesh[eid].c1 == 1 && tmesh[eid].type != 4)
			{
				SetTrunMat_C1Element(eid);
			}
			else if (tmesh[eid].c1 == 2)
			{
				SetTrunMat_Trs(eid);
			}
		}
	}
}

void TruncatedTspline::SetTrunMat_Irr(int eid)
{
	//at most one EP in an irregular element, not considering post-refinement case
	for (uint i = 0; i < tmesh[eid].bemat.size(); i++)
	{
		tmesh[eid].bemat[i].clear();
	}
	tmesh[eid].bemat.clear();
	int nv(cp[tmesh[eid].cnct[0]].face.size());
	tmesh[eid].bemat.resize(tmesh[eid].IEN.size(), vector<double>(16, 0.));

	int fc[7][4] = {{5, 6, 2 * nv + 1, 2 * nv + 2}, {5, 2 * nv + 2, 2 * nv + 3, 4}, {4, 5, 2 * nv + 2, 2 * nv + 3}, {4, 2 * nv + 3, 2 * nv + 4, 2 * nv + 5}, {4, 2 * nv + 5, 2 * nv + 6, 3}, {3, 4, 2 * nv + 5, 2 * nv + 6}, {3, 2 * nv + 6, 2 * nv + 7, 2}};
	double coef[4] = {4. / 9., 2. / 9., 1. / 9., 2. / 9.};
	int bzid[7][2] = {{-1, 3}, {7, 3}, {11, 15}, {-1, 15}, {14, 15}, {13, 12}, {-1, 12}};

	/*int enb[5];
	enb[1] = tmedge[tmesh[eid].edge[1]].face[0];
	if (enb[1] == eid) enb[1] = tmedge[tmesh[eid].edge[1]].face[1];
	enb[3] = tmedge[tmesh[eid].edge[2]].face[0];
	if (enb[3] == eid) enb[3] = tmedge[tmesh[eid].edge[2]].face[1];
	int* it = find(tmesh[enb[1]].cnct, tmesh[enb[1]].cnct + 4, tmesh[eid].cnct[1]);
	int loc = it - tmesh[enb[1]].cnct;
	enb[0] = tmedge[tmesh[enb[1]].edge[loc]].face[0];
	if (enb[0] == enb[1]) enb[0] = tmedge[tmesh[enb[1]].edge[loc]].face[1];
	loc = (loc + 2) % 4;
	enb[2] = tmedge[tmesh[enb[1]].edge[loc]].face[0];
	if (enb[2] == enb[1]) enb[2] = tmedge[tmesh[enb[1]].edge[loc]].face[1];
	it = find(tmesh[enb[3]].cnct, tmesh[enb[3]].cnct + 4, tmesh[eid].cnct[3]);
	loc = it - tmesh[enb[3]].cnct;
	loc = (loc + 3) % 4;
	enb[4] = tmedge[tmesh[enb[3]].edge[loc]].face[0];
	if (enb[4] == enb[3]) enb[4] = tmedge[tmesh[enb[3]].edge[loc]].face[1];
	int enb1[7] = { enb[0],enb[1],enb[1],enb[2],enb[3],enb[3],enb[4] };
	*/

	int ednb[4], cnnb[3];
	for (int i = 0; i < 4; i++)
	{
		ednb[i] = tmedge[tmesh[eid].edge[i]].face[0];
		if (ednb[i] == eid)
			ednb[i] = tmedge[tmesh[eid].edge[i]].face[1];
	}
	for (int i = 1; i < 4; i++)
	{
		cnnb[i - 1] = -1;
		for (uint j = 0; j < cp[tmesh[eid].cnct[i]].face.size(); j++)
		{
			int fcid(cp[tmesh[eid].cnct[i]].face[j]);
			if (fcid != eid && fcid != ednb[i] && fcid != ednb[(i + 3) % 4])
			{
				cnnb[i - 1] = fcid;
				break;
			}
		}
	}
	int enb[7] = {cnnb[0], ednb[1], ednb[1], cnnb[1], ednb[2], ednb[2], cnnb[2]};
	vector<int> trun(7, 0);
	for (int i = 0; i < 7; i++)
	{
		if (tmesh[enb[i]].c1 == 1 || tmesh[enb[i]].c1 == 11)
		{
			trun[i] = 1;
		}
	}

	for (int i = 0; i < 7; i++)
	{
		if (bzid[i][0] != -1 && trun[i] == 0)
		{
			for (int j = 0; j < 4; j++)
			{
				tmesh[eid].bemat[fc[i][j]][bzid[i][0]] += coef[j] / 2.;
			}
		}
		if (trun[i] == 0)
		{
			for (int j = 0; j < 4; j++)
			{
				tmesh[eid].bemat[fc[i][j]][bzid[i][1]] += coef[j] / 4.;
			}
		}
	}
}

void TruncatedTspline::SetTrunMat_Irr22(int eid)
{
	//set truncated C2 functions
	//at most one EP in an irregular element, not considering post-refinement case
	tmesh[eid].bemat22.clear();
	int nv(cp[tmesh[eid].cnct[0]].face.size());
	tmesh[eid].bemat22.resize(4);

	int fc[7][4] = {{5, 6, 2 * nv + 4, 2 * nv + 3}, {5, 2 * nv + 3, 2 * nv + 2, 4}, {4, 5, 2 * nv + 3, 2 * nv + 2}, {4, 2 * nv + 2, 2 * nv + 1, 2 * nv + 5}, {4, 2 * nv + 5, 2 * nv + 6, 3}, {3, 4, 2 * nv + 5, 2 * nv + 6}, {3, 2 * nv + 6, 2 * nv + 7, 2}};
	double coef[4] = {4. / 9., 2. / 9., 1. / 9., 2. / 9.}; //this is restricted, can be made more general
	int bzid[7][2] = {{-1, 3}, {7, 3}, {11, 15}, {-1, 15}, {14, 15}, {13, 12}, {-1, 12}};

	int ednb[4], cnnb[3];
	for (int i = 0; i < 4; i++)
	{
		ednb[i] = tmedge[tmesh[eid].edge[i]].face[0];
		if (ednb[i] == eid)
			ednb[i] = tmedge[tmesh[eid].edge[i]].face[1];
	}
	for (int i = 1; i < 4; i++)
	{
		cnnb[i - 1] = -1;
		for (uint j = 0; j < cp[tmesh[eid].cnct[i]].face.size(); j++)
		{
			int fcid(cp[tmesh[eid].cnct[i]].face[j]);
			if (fcid != eid && fcid != ednb[i] && fcid != ednb[(i + 3) % 4])
			{
				cnnb[i - 1] = fcid;
				break;
			}
		}
	}
	int enb[7] = {cnnb[0], ednb[1], ednb[1], cnnb[1], ednb[2], ednb[2], cnnb[2]};
	vector<int> trun(7, 0);
	for (int i = 0; i < 7; i++)
	{
		if (tmesh[enb[i]].c1 == 1)
		{
			trun[i] = 1;
		}
	}

	MatrixXd bemat0 = MatrixXd::Zero(tmesh[eid].IEN.size(), 16);
	for (int i = 0; i < 7; i++)
	{
		if (bzid[i][0] != -1 && trun[i] == 0)
		{
			for (int j = 0; j < 4; j++)
			{
				bemat0(fc[i][j], bzid[i][0]) += coef[j] / 2.;
			}
		}
		if (trun[i] == 0)
		{
			for (int j = 0; j < 4; j++)
			{
				bemat0(fc[i][j], bzid[i][1]) += coef[j] / 4.;
			}
		}
	}
	MatrixXd smat = MatrixXd::Zero(16, 49);
	double smat1d[7][4] = {{1., 0., 0., 0.}, {.5, .5, 0., 0.}, {.25, .5, .25, 0.}, {.125, .375, .375, .125}, {0., .25, .5, .25}, {0., 0., .5, .5}, {0., 0., 0., 1.}};
	int loc(0);
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			int loc0(0);
			for (int i0 = 0; i0 < 4; i0++)
			{
				for (int j0 = 0; j0 < 4; j0++)
				{
					smat(loc0, loc) = smat1d[j][j0] * smat1d[i][i0];
					loc0++;
				}
			}
			loc++;
		}
	}
	MatrixXd bemat = bemat0 * smat;

	/*cout << "bemat:\n";
	for (int i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		cout << "c2 ID: " << i << "\n";
		for (int j = 0; j < 49; j++)
		{
			cout << bemat(i, j) << "(" << j << ") ";
		}
		cout << "\n";
		getchar();
	}
	getchar();*/

	int bzcnct[4][16] = {{0, 1, 2, 3, 7, 8, 9, 10, 14, 15, 16, 17, 21, 22, 23, 24}, {3, 4, 5, 6, 10, 11, 12, 13, 17, 18, 19, 20, 24, 25, 26, 27}, {21, 22, 23, 24, 28, 29, 30, 31, 35, 36, 37, 38, 42, 43, 44, 45}, {24, 25, 26, 27, 31, 32, 33, 34, 38, 39, 40, 41, 45, 46, 47, 48}};
	for (int i = 0; i < 4; i++)
	{
		tmesh[eid].bemat22[i].resize(tmesh[eid].IEN.size(), vector<double>(16, 0.));
		for (int j = 0; j < tmesh[eid].IEN.size(); j++)
		{
			for (int k = 0; k < 16; k++)
			{
				tmesh[eid].bemat22[i][j][k] = bemat(j, bzcnct[i][k]);
			}
		}
	}
}

void TruncatedTspline::SetTrunMat_C1Element(int eid)
{
	//a transition element should also be a B-spline element without influence of T-junctions

	for (uint i = 0; i < tmesh[eid].bemat.size(); i++)
	{
		tmesh[eid].bemat[i].clear();
	}
	tmesh[eid].bemat.clear();

	//tmesh[eid].IENc1.clear();
	//for (uint i = 0; i < tmesh[eid].c1mat.size(); i++)
	//{
	//	tmesh[eid].c1mat[i].clear();
	//}
	//tmesh[eid].c1mat.clear();

	int ednb[4], cnnb[4];
	for (int i = 0; i < 4; i++)
	{
		ednb[i] = tmedge[tmesh[eid].edge[i]].face[0];
		if (ednb[i] == eid)
			ednb[i] = tmedge[tmesh[eid].edge[i]].face[1];
	}
	for (int i = 0; i < 4; i++)
	{
		cnnb[i] = -1;
		for (uint j = 0; j < cp[tmesh[eid].cnct[i]].face.size(); j++)
		{
			int fcid(cp[tmesh[eid].cnct[i]].face[j]);
			if (fcid != eid && fcid != ednb[i] && fcid != ednb[(i + 3) % 4])
			{
				cnnb[i] = fcid;
				break;
			}
		}
		if (cnnb[i] == -1)
		{
			cerr << "Can't find corner neighbor!\n";
			getchar();
		}
	}
	//knot vectors of C1 B-splines
	double ktsu[8] = {tmesh[eid].patch_ku[0][1], tmesh[eid].patch_ku[0][1], tmesh[eid].patch_ku[0][2], tmesh[eid].patch_ku[0][2],
					  tmesh[eid].patch_ku[0][3], tmesh[eid].patch_ku[0][3], tmesh[eid].patch_ku[0][4], tmesh[eid].patch_ku[0][4]};
	double ktsv[8] = {tmesh[eid].patch_kv[0][1], tmesh[eid].patch_kv[0][1], tmesh[eid].patch_kv[0][2], tmesh[eid].patch_kv[0][2],
					  tmesh[eid].patch_kv[0][3], tmesh[eid].patch_kv[0][3], tmesh[eid].patch_kv[0][4], tmesh[eid].patch_kv[0][4]};
	vector<vector<double>> c21mat(tmesh[eid].IEN.size(), vector<double>(16, 0.));
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		int loc(0);
		for (int j = 0; j < 4; j++) //v
		{
			for (int k = 0; k < 4; k++) //u
			{
				if (tmesh[eid].patch_ku[i][0] <= ktsu[k] && tmesh[eid].patch_ku[i][4] >= ktsu[k + 4] &&
					tmesh[eid].patch_kv[i][0] <= ktsv[j] && tmesh[eid].patch_kv[i][4] >= ktsv[j + 4])
				{
					vector<double> ku1, kv1;
					InsertKnotsC1(tmesh[eid].patch_ku[i], ku1);
					InsertKnotsC1(tmesh[eid].patch_kv[i], kv1);
					vector<double>::iterator itu = search(ku1.begin(), ku1.end(), ktsu + k, ktsu + k + 5);
					vector<double>::iterator itv = search(kv1.begin(), kv1.end(), ktsv + j, ktsv + j + 5);
					if (itu != ku1.end() && itv != kv1.end())
					{
						vector<double> ku(tmesh[eid].patch_ku[i].begin(), tmesh[eid].patch_ku[i].end());
						vector<double> kv(tmesh[eid].patch_kv[i].begin(), tmesh[eid].patch_kv[i].end());
						vector<vector<double>> Tu, Tv;
						TMatrix(ku, ku1, 3, Tu);
						TMatrix(kv, kv1, 3, Tv);
						int uloc = itu - ku1.begin();
						int vloc = itv - kv1.begin();
						c21mat[i][loc] = Tu[uloc][0] * Tv[vloc][0];
					}
				}
				loc++;
			}
		}
	}

	int cntrun[4] = {0, 3, 15, 12};
	int edtrun[4][2] = {{1, 2}, {7, 11}, {13, 14}, {4, 8}};
	int cnc0[4] = {0, 3, 15, 12};
	int cnc1[4][4] = {{0, 1, 4, 5}, {2, 3, 6, 7}, {10, 11, 14, 15}, {8, 9, 12, 13}};
	int edc0[8] = {1, 2, 7, 11, 14, 13, 8, 4};
	int edc1[8][2] = {{1, 5}, {2, 6}, {6, 7}, {10, 11}, {10, 14}, {9, 13}, {8, 9}, {4, 5}};
	int fcc0[4] = {5, 6, 10, 9};
	int fcc1[4] = {5, 6, 10, 9};

	if (tmesh[eid].c1 == 1)
	{
		for (int i = 0; i < 4; i++)
		{
			for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
			{
				c21mat[j][fcc1[i]] = 0.;
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		if (tmesh[cnnb[i]].c1 == 1)
		{
			for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
			{
				c21mat[j][cntrun[i]] = 0.;
			}
			int *it = find(tmesh[cnnb[i]].cnct, tmesh[cnnb[i]].cnct + 4, tmesh[eid].cnct[i]);
			int c1loc = it - tmesh[cnnb[i]].cnct;
			tmesh[eid].IENc1.push_back(tmesh[cnnb[i]].IENc1[c1loc]);
			vector<double> c1tmp(16, 0.);
			c1tmp[cnc0[i]] = 0.25;
			tmesh[eid].c1mat.push_back(c1tmp);
		}
		if (tmesh[ednb[i]].c1 == 1)
		{
			for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
			{
				c21mat[j][edtrun[i][0]] = 0.;
				c21mat[j][edtrun[i][1]] = 0.;
			}
			int *it = find(tmesh[ednb[i]].cnct, tmesh[ednb[i]].cnct + 4, tmesh[eid].cnct[i]);
			int c1loc = it - tmesh[ednb[i]].cnct;
			tmesh[eid].IENc1.push_back(tmesh[ednb[i]].IENc1[c1loc]);
			tmesh[eid].IENc1.push_back(tmesh[ednb[i]].IENc1[(c1loc + 3) % 4]);
			vector<double> c1tmp1(16, 0.), c1tmp2(16, 0.);
			c1tmp1[cnc0[i]] = 0.25;
			c1tmp1[edc0[2 * i]] = 0.5;
			tmesh[eid].c1mat.push_back(c1tmp1);
			c1tmp2[cnc0[(i + 1) % 4]] = 0.25;
			c1tmp2[edc0[2 * i + 1]] = 0.5;
			tmesh[eid].c1mat.push_back(c1tmp2);
		}
	}

	tmesh[eid].bemat.resize(tmesh[eid].IEN.size(), vector<double>(16, 0.));
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		int j, k;
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < 4; k++)
			{
				tmesh[eid].bemat[i][cnc0[j]] += c21mat[i][cnc1[j][k]] / 4.;
			}
		}
		for (j = 0; j < 8; j++)
		{
			for (k = 0; k < 2; k++)
			{
				tmesh[eid].bemat[i][edc0[j]] += c21mat[i][edc1[j][k]] / 2.;
			}
		}
		for (j = 0; j < 4; j++)
		{
			tmesh[eid].bemat[i][fcc0[j]] = c21mat[i][fcc1[j]];
		}
	}
}

void TruncatedTspline::SetTrunMat_Trs(int eid)
{
	//a transition element should also be a B-spline element without influence of T-junctions

	for (uint i = 0; i < tmesh[eid].bemat.size(); i++)
	{
		tmesh[eid].bemat[i].clear();
	}
	tmesh[eid].bemat.clear();

	tmesh[eid].IENc1.clear();
	for (uint i = 0; i < tmesh[eid].c1mat.size(); i++)
	{
		tmesh[eid].c1mat[i].clear();
	}
	tmesh[eid].c1mat.clear();

	int ednb[4], cnnb[4];
	for (int i = 0; i < 4; i++)
	{
		ednb[i] = tmedge[tmesh[eid].edge[i]].face[0];
		if (ednb[i] == eid)
			ednb[i] = tmedge[tmesh[eid].edge[i]].face[1];
	}
	for (int i = 0; i < 4; i++)
	{
		cnnb[i] = -1;
		for (uint j = 0; j < cp[tmesh[eid].cnct[i]].face.size(); j++)
		{
			int fcid(cp[tmesh[eid].cnct[i]].face[j]);
			if (fcid != eid && fcid != ednb[i] && fcid != ednb[(i + 3) % 4])
			{
				cnnb[i] = fcid;
				break;
			}
		}
		if (cnnb[i] == -1)
		{
			cerr << "Can't find corner neighbor!\n";
			getchar();
		}
	}
	//knot vectors of C1 B-splines
	double ktsu[8] = {tmesh[eid].patch_ku[0][1], tmesh[eid].patch_ku[0][1], tmesh[eid].patch_ku[0][2], tmesh[eid].patch_ku[0][2],
					  tmesh[eid].patch_ku[0][3], tmesh[eid].patch_ku[0][3], tmesh[eid].patch_ku[0][4], tmesh[eid].patch_ku[0][4]};
	double ktsv[8] = {tmesh[eid].patch_kv[0][1], tmesh[eid].patch_kv[0][1], tmesh[eid].patch_kv[0][2], tmesh[eid].patch_kv[0][2],
					  tmesh[eid].patch_kv[0][3], tmesh[eid].patch_kv[0][3], tmesh[eid].patch_kv[0][4], tmesh[eid].patch_kv[0][4]};
	vector<vector<double>> c21mat(tmesh[eid].IEN.size(), vector<double>(16, 0.));
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		int loc(0);
		for (int j = 0; j < 4; j++) //v
		{
			for (int k = 0; k < 4; k++) //u
			{
				if (tmesh[eid].patch_ku[i][0] <= ktsu[k] && tmesh[eid].patch_ku[i][4] >= ktsu[k + 4] &&
					tmesh[eid].patch_kv[i][0] <= ktsv[j] && tmesh[eid].patch_kv[i][4] >= ktsv[j + 4])
				{
					vector<double> ku1, kv1;
					InsertKnotsC1(tmesh[eid].patch_ku[i], ku1);
					InsertKnotsC1(tmesh[eid].patch_kv[i], kv1);
					vector<double>::iterator itu = search(ku1.begin(), ku1.end(), ktsu + k, ktsu + k + 5);
					vector<double>::iterator itv = search(kv1.begin(), kv1.end(), ktsv + j, ktsv + j + 5);
					if (itu != ku1.end() && itv != kv1.end())
					{
						vector<double> ku(tmesh[eid].patch_ku[i].begin(), tmesh[eid].patch_ku[i].end());
						vector<double> kv(tmesh[eid].patch_kv[i].begin(), tmesh[eid].patch_kv[i].end());
						vector<vector<double>> Tu, Tv;
						TMatrix(ku, ku1, 3, Tu);
						TMatrix(kv, kv1, 3, Tv);
						int uloc = itu - ku1.begin();
						int vloc = itv - kv1.begin();
						c21mat[i][loc] = Tu[uloc][0] * Tv[vloc][0];
					}
				}
				loc++;
			}
		}
	}

	int cntrun[4] = {0, 3, 15, 12};
	int edtrun[4][2] = {{1, 2}, {7, 11}, {13, 14}, {4, 8}};
	int cnc0[4] = {0, 3, 15, 12};
	int cnc1[4][4] = {{0, 1, 4, 5}, {2, 3, 6, 7}, {10, 11, 14, 15}, {8, 9, 12, 13}};
	int edc0[8] = {1, 2, 7, 11, 14, 13, 8, 4};
	int edc1[8][2] = {{1, 5}, {2, 6}, {6, 7}, {10, 11}, {10, 14}, {9, 13}, {8, 9}, {4, 5}};
	int fcc0[4] = {5, 6, 10, 9};
	int fcc1[4] = {5, 6, 10, 9};

	for (int i = 0; i < 4; i++)
	{
		if (tmesh[cnnb[i]].c1 == 1)
		{
			for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
			{
				c21mat[j][cntrun[i]] = 0.;
			}
			int *it = find(tmesh[cnnb[i]].cnct, tmesh[cnnb[i]].cnct + 4, tmesh[eid].cnct[i]);
			int c1loc = it - tmesh[cnnb[i]].cnct;
			tmesh[eid].IENc1.push_back(tmesh[cnnb[i]].IENc1[c1loc]);
			vector<double> c1tmp(16, 0.);
			c1tmp[cnc0[i]] = 0.25;
			tmesh[eid].c1mat.push_back(c1tmp);
		}
		if (tmesh[ednb[i]].c1 == 1)
		{
			for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
			{
				c21mat[j][edtrun[i][0]] = 0.;
				c21mat[j][edtrun[i][1]] = 0.;
			}
			int *it = find(tmesh[ednb[i]].cnct, tmesh[ednb[i]].cnct + 4, tmesh[eid].cnct[i]);
			int c1loc = it - tmesh[ednb[i]].cnct;
			tmesh[eid].IENc1.push_back(tmesh[ednb[i]].IENc1[c1loc]);
			tmesh[eid].IENc1.push_back(tmesh[ednb[i]].IENc1[(c1loc + 3) % 4]);
			vector<double> c1tmp1(16, 0.), c1tmp2(16, 0.);
			c1tmp1[cnc0[i]] = 0.25;
			c1tmp1[edc0[2 * i]] = 0.5;
			tmesh[eid].c1mat.push_back(c1tmp1);
			c1tmp2[cnc0[(i + 1) % 4]] = 0.25;
			c1tmp2[edc0[2 * i + 1]] = 0.5;
			tmesh[eid].c1mat.push_back(c1tmp2);
		}
	}

	tmesh[eid].bemat.resize(tmesh[eid].IEN.size(), vector<double>(16, 0.));
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		int j, k;
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < 4; k++)
			{
				tmesh[eid].bemat[i][cnc0[j]] += c21mat[i][cnc1[j][k]] / 4.;
			}
		}
		for (j = 0; j < 8; j++)
		{
			for (k = 0; k < 2; k++)
			{
				tmesh[eid].bemat[i][edc0[j]] += c21mat[i][edc1[j][k]] / 2.;
			}
		}
		for (j = 0; j < 4; j++)
		{
			tmesh[eid].bemat[i][fcc0[j]] = c21mat[i][fcc1[j]];
		}
	}
}

void TruncatedTspline::CollectActives_DPatchAnalysis()
{
	paid.clear();
	paid.resize(cp.size() + bzcp.size());
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			paid[i] = count;
			count++;
		}
		else
		{
			paid[i] = -1;
		}
	}
	for (uint i = 0; i < bzcp.size(); i++)
	{
		paid[i + cp.size()] = count;
		count++;
	}
	npta = count;

	eaid.clear();
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].type != 2 && tmesh[i].type != 3)
		{
			eaid.push_back(i);
		}
	}
}

void TruncatedTspline::GeomMap_DPatch(int eid, double u, double v, array<double, 3> &pt)
{
	vector<double> Nt;
	vector<array<double, 2>> dNdt;
	ElementBasis_DPatch(eid, u, v, Nt, dNdt);
	pt[0] = 0.;
	pt[1] = 0.;
	pt[2] = 0.;
	double sum(0.);
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		pt[0] += cp[tmesh[eid].IEN[i]].coor[0] * Nt[i];
		pt[1] += cp[tmesh[eid].IEN[i]].coor[1] * Nt[i];
		pt[2] += cp[tmesh[eid].IEN[i]].coor[2] * Nt[i];
		sum += Nt[i];
	}
	int ist(tmesh[eid].IEN.size());
	for (uint i = 0; i < tmesh[eid].IENc1.size(); i++)
	{
		pt[0] += bzcp[tmesh[eid].IENc1[i]][0] * Nt[ist + i];
		pt[1] += bzcp[tmesh[eid].IENc1[i]][1] * Nt[ist + i];
		pt[2] += bzcp[tmesh[eid].IENc1[i]][2] * Nt[ist + i];
		sum += Nt[ist + i];
	}

	//if (tmesh[eid].c1 == 1 && tmesh[eid].type != 4 /*&& fabs(sum - 1.) > 1.e-6*/)
	//{
	//	cout << "eid: " << eid << "\n";
	//	//cout << "uv: " << u << " " << v << "\n";
	//	cout << "sum: " << sum << "\n";
	//	cout << "c2: ";
	//	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	//	{
	//		cout << Nt[i] << " ";
	//	}
	//	cout << "\n";
	//	cout << "c1: ";
	//	for (uint i = 0; i < tmesh[eid].IENc1.size(); i++)
	//	{
	//		cout << Nt[ist+i] << " ";
	//	}
	//	cout << "\n";
	//	//cout << "c1mat 0:\n";
	//	//cout << tmesh[eid].c1mat22[0] << "\n";
	//	//cout << "c1mat 1:\n";
	//	//cout << tmesh[eid].c1mat22[1] << "\n";
	//	//cout << "c1mat 2:\n";
	//	//cout << tmesh[eid].c1mat22[2] << "\n";
	//	//cout << "c1mat 3:\n";
	//	//cout << tmesh[eid].c1mat22[3] << "\n";
	//	//cout << "bemat 0:\n";
	//	//cout << tmesh[eid].bemat22[0] << "\n";
	//	//cout << "bemat 1:\n";
	//	//cout << tmesh[eid].bemat22[1] << "\n";
	//	//cout << "bemat 2:\n";
	//	//cout << tmesh[eid].bemat22[2] << "\n";
	//	//cout << "bemat 3:\n";
	//	//cout << tmesh[eid].bemat22[3] << "\n";
	//	getchar();
	//}
}

void TruncatedTspline::ElementBasis_DPatch(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	if (tmesh[eid].act == 1)
	{
		if (tmesh[eid].c1 == 0 && (tmesh[eid].type == 0 || tmesh[eid].type == 1))
		{
			ElementBasis_Regular(eid, u, v, Nt, dNdt);
		}
		else if (tmesh[eid].c1 == 1 && tmesh[eid].type == 4)
		{
			//ElementBasis_Irregular_DPatch(eid, u, v, Nt, dNdt);
			//ElementBasis_Irregular_DPatch22(eid, u, v, Nt, dNdt);
			ElementBasis_DPatch_Irr22(eid, u, v, Nt, dNdt);
		}
		else if (tmesh[eid].c1 == 1 && tmesh[eid].type != 4)
		{
			ElementBasis_DPatch_Trs(eid, u, v, Nt, dNdt);
		}
		else if (tmesh[eid].c1 == 2)
		{
			ElementBasis_DPatch_Trs(eid, u, v, Nt, dNdt);
		}
	}
}

void TruncatedTspline::ElementBasis_DPatch_Irr22(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	//need to change variabe for Bezier, (u,v) -> (u_b,v_b)
	Nt.clear();
	dNdt.clear();
	Nt.resize(tmesh[eid].IEN.size() + tmesh[eid].IENc1.size(), 0.);
	dNdt.resize(tmesh[eid].IEN.size() + tmesh[eid].IENc1.size());
	//uint nv(cp[tmesh[eid].cnct[0]].face.size());
	BezierElement be;
	double Nt0[16], dNdt0[16][2];
	//double u_b(u / tmedge[tmesh[eid].edge[0]].len), v_b(v / tmedge[tmesh[eid].edge[3]].len);
	//be.Basis(u_b, v_b, Nt0, dNdt0);
	double edm[2] = {tmedge[tmesh[eid].edge[0]].len, tmedge[tmesh[eid].edge[3]].len};
	double u_b(u / edm[0]), v_b(v / edm[1]);
	double uhalf(.5), vhalf(.5);
	int subeid(0);
	if (u_b >= uhalf && v_b < vhalf)
	{
		subeid = 1;
		u_b = 2. * (u_b - uhalf);
		v_b = 2. * v_b;
	}
	else if (u_b < uhalf && v_b >= vhalf)
	{
		subeid = 2;
		u_b = 2. * u_b;
		v_b = 2. * (v_b - vhalf);
	}
	else if (u_b >= uhalf && v_b >= vhalf)
	{
		subeid = 3;
		u_b = 2. * (u_b - uhalf);
		v_b = 2. * (v_b - vhalf);
	}
	else
	{
		subeid = 0;
		u_b = 2. * u_b;
		v_b = 2. * v_b;
	}
	be.Basis(u_b, v_b, Nt0, dNdt0);

	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		dNdt[i][0] = 0.;
		dNdt[i][1] = 0.;
		for (int j = 0; j < 16; j++)
		{
			Nt[i] += tmesh[eid].bemat22[subeid][i][j] * Nt0[j];
			dNdt[i][0] += tmesh[eid].bemat22[subeid][i][j] * dNdt0[j][0];
			dNdt[i][1] += tmesh[eid].bemat22[subeid][i][j] * dNdt0[j][1];
		}
		dNdt[i][0] *= 2. / edm[0];
		dNdt[i][1] *= 2. / edm[1];
	}
	int ist(tmesh[eid].IEN.size());
	for (uint i = 0; i < tmesh[eid].IENc1.size(); i++)
	{
		int iloc(ist + i);
		dNdt[iloc][0] = 0.;
		dNdt[iloc][1] = 0.;
		for (int j = 0; j < 16; j++)
		{
			Nt[iloc] += tmesh[eid].c1mat22[subeid][i][j] * Nt0[j];
			dNdt[iloc][0] += tmesh[eid].c1mat22[subeid][i][j] * dNdt0[j][0];
			dNdt[iloc][1] += tmesh[eid].c1mat22[subeid][i][j] * dNdt0[j][1];
		}
		dNdt[iloc][0] *= 2. / edm[0];
		dNdt[iloc][1] *= 2. / edm[1];
	}
}

void TruncatedTspline::ElementBasis_DPatch_Trs(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	//need to change variabe for Bezier, (u,v) -> (u_b,v_b)
	Nt.clear();
	dNdt.clear();
	Nt.resize(tmesh[eid].IEN.size() + tmesh[eid].IENc1.size(), 0.);
	dNdt.resize(tmesh[eid].IEN.size() + tmesh[eid].IENc1.size());
	//uint nv(cp[tmesh[eid].cnct[0]].face.size());
	BezierElement be;
	double Nt0[16], dNdt0[16][2];
	double u_b(u / tmedge[tmesh[eid].edge[0]].len), v_b(v / tmedge[tmesh[eid].edge[3]].len);
	be.Basis(u_b, v_b, Nt0, dNdt0);
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		dNdt[i][0] = 0.;
		dNdt[i][1] = 0.;
		for (int j = 0; j < 16; j++)
		{
			Nt[i] += tmesh[eid].bemat[i][j] * Nt0[j];
			dNdt[i][0] += tmesh[eid].bemat[i][j] * dNdt0[j][0];
			dNdt[i][1] += tmesh[eid].bemat[i][j] * dNdt0[j][1];
		}
		dNdt[i][0] /= tmedge[tmesh[eid].edge[0]].len;
		dNdt[i][1] /= tmedge[tmesh[eid].edge[3]].len;
	}
	int ist(tmesh[eid].IEN.size());
	for (uint i = 0; i < tmesh[eid].IENc1.size(); i++)
	{
		int iloc(ist + i);
		dNdt[iloc][0] = 0.;
		dNdt[iloc][1] = 0.;
		for (int j = 0; j < 16; j++)
		{
			Nt[iloc] += tmesh[eid].c1mat[i][j] * Nt0[j];
			dNdt[iloc][0] += tmesh[eid].c1mat[i][j] * dNdt0[j][0];
			dNdt[iloc][1] += tmesh[eid].c1mat[i][j] * dNdt0[j][1];
		}
		dNdt[iloc][0] /= tmedge[tmesh[eid].edge[0]].len;
		dNdt[iloc][1] /= tmedge[tmesh[eid].edge[3]].len;
	}
}

void TruncatedTspline::FindGlobalRefineID(vector<int> &rfid, vector<int> &rftype)
{
	rfid.clear();
	rftype.clear();
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1)
		{
			if (tmesh[eid].type == 0) //note that there will no type 1
			{
				rfid.push_back(eid);
				rftype.push_back(0);
			}
			//else if (tmesh[eid].type == 1)
			//{
			//	rfid.push_back(eid);
			//	rftype.push_back(1);
			//}
			else if (tmesh[eid].type == 2)
			{
				rfid.push_back(eid);
				rftype.push_back(3);
			}
			else if (tmesh[eid].type == 4)
			{
				rfid.push_back(eid);
				rftype.push_back(4);
			}
		}
	}
}

void TruncatedTspline::FindC1Element(vector<int> &c1id)
{
	c1id.clear();
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1 && tmesh[eid].c1 == 1)
		{
			c1id.push_back(eid);
		}
	}
}

void TruncatedTspline::Refine_C1Element(const vector<int> &c1id)
{
	//need to add face points
	bzcp.clear();
	int eb[4] = {0, 1, 3, 2};
	int bzfc[4] = {5, 6, 10, 9};
	int bzcnct[4][16] = {{0, 1, 2, 3, 7, 8, 9, 10, 14, 15, 16, 17, 21, 22, 23, 24}, {3, 4, 5, 6, 10, 11, 12, 13, 17, 18, 19, 20, 24, 25, 26, 27}, {21, 22, 23, 24, 28, 29, 30, 31, 35, 36, 37, 38, 42, 43, 44, 45}, {24, 25, 26, 27, 31, 32, 33, 34, 38, 39, 40, 41, 45, 46, 47, 48}};

	for (uint i = 0; i < c1id.size(); i++)
	{
		int eid(c1id[i]);
		if (tmesh[eid].type == 4 && tmesh[eid].bzpts.size() == 16 * 4)
		{
			for (int j = 0; j < 4; j++)
			{
				if (tmesh[eid].chd[j] != -1) //must be true
				{
					int chid(tmesh[eid].chd[j]);
					tmesh[chid].c1 = tmesh[eid].c1; //==1
					tmesh[chid].bzpts.clear();
					tmesh[chid].bzpts.resize(16);
					for (int k = 0; k < 16; k++)
					{
						tmesh[chid].bzpts[k] = tmesh[eid].bzpts[16 * eb[j] + k];
					}
					tmesh[chid].IENc1.clear();
					for (int k = 0; k < 4; k++)
					{
						bzcp.push_back(tmesh[chid].bzpts[bzfc[k]]);
						tmesh[chid].IENc1.push_back(bzcp.size() - 1);
					}
				}
			}
		}
		else
		{
			//refinement first
			vector<array<double, 3>> pts1;
			//cout << tmesh[eid].bzpts.size() << "\n";
			BezierRefineBi3(tmesh[eid].bzpts, pts1);
			for (int j = 0; j < 4; j++)
			{
				if (tmesh[eid].chd[j] != -1) //must be true
				{
					int chid(tmesh[eid].chd[j]);
					tmesh[chid].c1 = tmesh[eid].c1; //==1
					tmesh[chid].bzpts.clear();
					tmesh[chid].bzpts.resize(16);
					for (int k = 0; k < 16; k++)
					{
						tmesh[chid].bzpts[k] = pts1[bzcnct[eb[j]][k]];
					}
					tmesh[chid].IENc1.clear();
					for (int k = 0; k < 4; k++)
					{
						bzcp.push_back(tmesh[chid].bzpts[bzfc[k]]);
						tmesh[chid].IENc1.push_back(bzcp.size() - 1);
					}
				}
			}
		}
	}
}

void TruncatedTspline::BezierMeshUpdateIndex(const vector<BezierElement> &bzmesh0, vector<BezierElement> &bzmesh)
{
	bzmesh.resize(bzmesh0.size());
	for (uint eid = 0; eid < bzmesh0.size(); eid++)
	{
		bzmesh[eid].prt = bzmesh0[eid].prt;
		for (uint i = 0; i < bzmesh0[eid].IEN.size(); i++)
		{
			if (paid[bzmesh0[eid].IEN[i]] != -1)
			{
				bzmesh[eid].IEN.push_back(paid[bzmesh0[eid].IEN[i]]);
				bzmesh[eid].cmat.push_back(bzmesh0[eid].cmat[i]);
			}
		}
		for (int i = 0; i < 16; i++)
		{
			bzmesh[eid].pts[i][0] = bzmesh0[eid].pts[i][0];
			bzmesh[eid].pts[i][1] = bzmesh0[eid].pts[i][1];
			bzmesh[eid].pts[i][2] = bzmesh0[eid].pts[i][2];
		}
	}
}

void TruncatedTspline::BC_Square(vector<int> &IDBC, vector<double> &gh)
{
	IDBC.clear();
	gh.clear();
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			count++;
		}
	}
	IDBC.resize(count + bzcp.size(), -1);
	gh.resize(IDBC.size(), 0.);
	int loc(0);
	count = 0;
	double tol(1.e-3);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			if (fabs(cp[i].coor[0]) < tol || fabs(cp[i].coor[0] - 1.) < tol ||
				fabs(cp[i].coor[1]) < tol || fabs(cp[i].coor[1] - 1.) < tol)
			{
				IDBC[loc] = -1;
				gh[loc] = 0.;
			}
			else
			{
				IDBC[loc] = count++;
			}
			loc++;
		}
	}
	for (uint i = 0; i < bzcp.size(); i++)
	{
		IDBC[loc] = count++;
		loc++;
	}
}

void TruncatedTspline::BC_Hemisphere(vector<int> &IDBC, vector<double> &gh)
{
	IDBC.clear();
	gh.clear();
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			count++;
		}
	}
	IDBC.resize(count + bzcp.size(), -1);
	gh.resize(IDBC.size(), 0.);
	int loc(0);
	count = 0;
	double tol(1.e-3);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			IDBC[loc] = count++;
			loc++;
		}
	}
	for (uint i = 0; i < bzcp.size(); i++)
	{
		IDBC[loc] = count++;
		loc++;
	}
}

void TruncatedTspline::ReadBezier(string fn, vector<BezierElement> &bzmesh)
{
	vector<array<double, 3>> pts;
	ifstream fin;
	fin.open(fn);
	if (fin.is_open())
	{
		string stmp;
		int itmp, ncp, nel, loc, nIEN;
		double dtmp;
		fin >> stmp >> itmp;
		fin >> stmp >> itmp >> itmp;
		fin >> stmp >> ncp;
		cout << stmp << " " << ncp << "\n";
		pts.resize(ncp);
		for (int i = 0; i < ncp; i++)
		{
			fin >> loc >> pts[i][0] >> pts[i][1] >> pts[i][2] >> dtmp >> itmp;
		}
		fin >> stmp >> nel;
		cout << "elems " << stmp << "\n";
		bzmesh.resize(nel);
		for (int i = 0; i < nel; i++)
		{
			fin >> loc >> nIEN;
			bzmesh[i].IEN.resize(nIEN);
			bzmesh[i].cmat.resize(nIEN);
			for (int j = 0; j < nIEN; j++)
			{
				fin >> bzmesh[i].IEN[j];
			}
			for (int j = 0; j < nIEN; j++)
			{
				for (int k = 0; k < 16; k++)
				{
					fin >> bzmesh[i].cmat[j][k];
				}
			}
			for (int k = 0; k < 16; k++)
			{
				bzmesh[i].pts[k][0] = 0.;
				bzmesh[i].pts[k][1] = 0.;
				bzmesh[i].pts[k][2] = 0.;
				for (int j = 0; j < nIEN; j++)
				{
					bzmesh[i].pts[k][0] += bzmesh[i].cmat[j][k] * pts[bzmesh[i].IEN[j]][0];
					bzmesh[i].pts[k][1] += bzmesh[i].cmat[j][k] * pts[bzmesh[i].IEN[j]][1];
					bzmesh[i].pts[k][2] += bzmesh[i].cmat[j][k] * pts[bzmesh[i].IEN[j]][2];
				}
			}
		}

		fin.close();
	}
	else
	{
		cerr << "Can't open " << fn << "!\n";
	}
}

void TruncatedTspline::VisualizeCM_DPatch(string fn)
{
	vector<int> IDa(cp.size() + bzcp.size(), -1);
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
			IDa[i] = count++;
	}
	int ncpa(count);
	for (uint i = 0; i < bzcp.size(); i++)
	{
		IDa[cp.size() + i] = count++;
	}

	string fname(fn + "/controlmesh_label.vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << ncpa + bzcp.size() << " float\n";
		fout.precision(16);
		for (uint i = 0; i < cp.size(); i++)
		{
			if (cp[i].act == 1)
			{
				fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
			}
		}
		for (uint i = 0; i < bzcp.size(); i++)
		{
			fout << bzcp[i][0] << " " << bzcp[i][1] << " " << bzcp[i][2] << "\n";
		}
		vector<int> quad, qdc1;
		for (uint i = 0; i < tmesh.size(); i++)
		{
			if (tmesh[i].act == 1)
			{
				int flag(0);
				for (int j = 0; j < 4; j++)
				{
					if (cp[tmesh[i].cnct[j]].act == 0)
					{
						flag = 1;
						break;
					}
				}
				if (flag == 0)
				{
					quad.push_back(i);
				}
				if (tmesh[i].c1 == 1)
					qdc1.push_back(i);
			}
		}
		int nquad(quad.size() + qdc1.size());
		fout << "\nCELLS " << nquad << " " << 5 * nquad << '\n';
		for (uint i = 0; i < quad.size(); i++)
		{
			int i1(quad[i]);
			fout << "4 " << IDa[tmesh[i1].cnct[0]] << " " << IDa[tmesh[i1].cnct[1]] << " "
				 << IDa[tmesh[i1].cnct[2]] << " " << IDa[tmesh[i1].cnct[3]] << "\n";
		}
		for (uint i = 0; i < qdc1.size(); i++)
		{
			int i1(qdc1[i]);
			fout << "4 " << IDa[cp.size() + tmesh[i1].IENc1[0]] << " " << IDa[cp.size() + tmesh[i1].IENc1[1]] << " "
				 << IDa[cp.size() + tmesh[i1].IENc1[2]] << " " << IDa[cp.size() + tmesh[i1].IENc1[3]] << "\n";
		}
		fout << "\nCELL_TYPES " << nquad << '\n';
		for (int i = 0; i < nquad; i++)
		{
			fout << "9\n";
		}
		//fout << "\nCELL_DATA " << nquad << "\nSCALARS type float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<tmesh.size(); i++)
		//{
		//	if (tmesh[i].act == 1 && tmesh[i].type != 6)
		//	{
		//		//fout << tmesh[i].dual << "\n";
		//		fout << tmesh[i].c1 << "\n";
		//	}
		//}
		//for (uint i = 0; i<tmesh.size(); i++)
		//{
		//	if (tmesh[i].act == 1 && tmesh[i].c1 == 1)
		//	{
		//		fout << "0\n";
		//	}
		//}
		fout << "POINT_DATA " << ncpa + bzcp.size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < cp.size(); i++)
		{
			if (cp[i].act == 1)
			{
				fout << cp[i].label << "\n";
				//fout<<"1\n";
			}
		}
		for (uint i = 0; i < bzcp.size(); i++)
		{
			fout << -1 << "\n";
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
		cout << "Cannot open " << fname << "!\n";
	}

	string fname4(fn + "/initial_velocityfield.txt");
	//ofstream fout;
	fout.open(fname4.c_str());
	if (fout.is_open())
	{
		for (uint i = 0; i < cp.size(); i++)
		{
			if (cp[i].act == 1)
			{
				fout << cp[i].velocity[0] << " " << cp[i].velocity[1] << " " << cp[i].velocity[2] << "\n";
			}
		}
		for (uint i = 0; i < bzcp.size(); i++)
		{
			fout << "0 0 0"
				 << "\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname4 << '\n';
	}
}

void TruncatedTspline::SetBezierMatIrr_NonUniform(int eid, vector<vector<double>> &bemat)
{
	//outermost 7 ordered arbitrarily, works with FindIEN4()
	//some coefficients not calculated because later they will be set as zero
	for (uint i = 0; i < bemat.size(); i++)
		bemat[i].clear();
	bemat.clear();
	int nv(cp[tmesh[eid].cnct[0]].face.size());
	int nbf(2 * nv + 8), n1r(2 * nv + 1), nout(7);
	bemat.resize(nbf, vector<double>(16, 0.));
	double ki[2][7] = {{1., 0., 0., 0., 0., 0., 1.}, {1., 0., 0., 0., 0., 0., 1.}};
	ki[0][3] = tmedge[tmesh[eid].edge[0]].len;
	ki[1][3] = tmedge[tmesh[eid].edge[3]].len;
	double kitmp[4][2] = {{0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}};
	for (int i = 0; i < 4; i++)
	{
		int nb1(tmedge[tmesh[eid].edge[i]].face[0]);
		if (nb1 == eid)
			nb1 = tmedge[tmesh[eid].edge[i]].face[1];
		int *it1 = find(tmesh[nb1].edge, tmesh[nb1].edge + 4, tmesh[eid].edge[i]);
		int loc1(it1 - tmesh[nb1].edge);
		kitmp[i][0] = tmedge[tmesh[nb1].edge[(loc1 + 1) % 4]].len;
		loc1 = (loc1 + 2) % 4;
		if (tmedge[tmesh[nb1].edge[loc1]].face.size() == 2)
		{
			int nb2(tmedge[tmesh[nb1].edge[loc1]].face[0]);
			if (nb2 == nb1)
				nb2 = tmedge[tmesh[nb1].edge[loc1]].face[1];
			int *it2 = find(tmesh[nb2].edge, tmesh[nb2].edge + 4, tmesh[nb1].edge[loc1]);
			int loc2(it2 - tmesh[nb2].edge);
			kitmp[i][1] = tmedge[tmesh[nb2].edge[(loc2 + 1) % 4]].len;
		}
		else
		{
			kitmp[i][1] = 0.;
		}
	}
	ki[0][4] = kitmp[1][0];
	ki[0][5] = kitmp[1][1];
	ki[0][2] = kitmp[3][0];
	ki[0][1] = kitmp[3][1];
	ki[1][4] = kitmp[2][0];
	ki[1][5] = kitmp[2][1];
	ki[1][2] = kitmp[0][0];
	ki[1][1] = kitmp[0][1];
	array<vector<double>, 2> kv;
	kv[0].resize(8);
	kv[1].resize(8);
	kv[0][3] = 0.;
	kv[0][4] = ki[0][3];
	kv[1][3] = 0.;
	kv[1][4] = ki[1][3];
	for (int i = 0; i < 3; i++)
	{
		kv[0][i + 5] = kv[0][i + 4] + ki[0][i + 4];
		kv[0][2 - i] = kv[0][3 - i] - ki[0][2 - i];
		kv[1][i + 5] = kv[1][i + 4] + ki[1][i + 4];
		kv[1][2 - i] = kv[1][3 - i] - ki[1][2 - i];
	}
	double kvbtmp[2][12] = {{kv[0][0], kv[0][1], kv[0][2], kv[0][3], kv[0][3], kv[0][3], kv[0][4], kv[0][4], kv[0][4], kv[0][5], kv[0][6], kv[0][7]},
							{kv[1][0], kv[1][1], kv[1][2], kv[1][3], kv[1][3], kv[1][3], kv[1][4], kv[1][4], kv[1][4], kv[1][5], kv[1][6], kv[1][7]}};
	array<vector<double>, 2> kvb;
	kvb[0].assign(kvbtmp[0], kvbtmp[0] + 12);
	kvb[1].assign(kvbtmp[1], kvbtmp[1] + 12);
	vector<vector<double>> umat, vmat;
	TMatrix(kv[0], kvb[0], 3, umat);
	TMatrix(kv[1], kvb[1], 3, vmat);
	int st[2] = {2, 2}, loc0(0);
	vector<vector<double>> cmat(16, vector<double>(16, 0.));
	for (int i0 = 0; i0 < 4; i0++)
	{
		for (int j0 = 0; j0 < 4; j0++)
		{
			int loc1(0);
			for (int i1 = 0; i1 < 4; i1++)
			{
				for (int j1 = 0; j1 < 4; j1++)
				{
					cmat[loc0][loc1] = umat[st[0] + j1][j0] * vmat[st[1] + i1][i0];
					loc1++;
				}
			}
			loc0++;
		}
	}

	int id[15][2] = {{1, 4}, {2, 8}, {3, 9}, {4, 10}, {5, 6}, {6, 2}, {7, 1}, {2 * nv + 1, 15}, {2 * nv + 2, 11}, {2 * nv + 3, 7}, {2 * nv + 4, 3}, {2 * nv + 5, 14}, {2 * nv + 6, 13}, {2 * nv + 7, 12}};
	if (nv == 3)
		id[6][0] = 1;
	for (int i = 0; i < 15; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			bemat[id[i][0]][j] += cmat[id[i][1]][j];
		}
	}
	//
	//    for(int i=0; i<2; i++)
	//    {
	//        ki_out[i].clear();
	//        ki_out[i].resize(5);
	//        for(int j=0; j<5; j++)
	//        {
	//            ki_out[i][j]=ki[i][j+1];
	//        }
	//    }
}

void TruncatedTspline::FindSplineC2_Scale()
{
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1)
		{
			if (tmesh[eid].c1 == 1 && tmesh[eid].type == 4)
			{
				FindSplineC2_Scale_Irr22(eid);
			}
			else if (tmesh[eid].c1 == 1 && tmesh[eid].type != 4)
			{
				FindSplineC2_Scale_Trs(eid);
			}
			else if (tmesh[eid].c1 == 2)
			{
				FindSplineC2_Scale_Trs(eid);
			}
		}
	}
}

void TruncatedTspline::FindSplineC2_Scale_Irr22(int eid)
{
	//set truncated C2 functions
	//at most one EP in an irregular element, not considering post-refinement case
	int enpta(0);
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		if (cp[tmesh[eid].IEN[i]].act == 1)
			enpta++;
	}
	if (enpta == 0)
		return;

	tmesh[eid].bemat22.clear();
	tmesh[eid].IENc2.clear();
	int nv(cp[tmesh[eid].cnct[0]].face.size());
	tmesh[eid].bemat22.resize(4);

	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		if (cp[tmesh[eid].IEN[i]].act == 1)
		{
			tmesh[eid].IENc2.push_back(tmesh[eid].IEN[i]);
		}
	}

	vector<int> b2ed(nv);
	vector<array<int, 2>> b1ed(nv);
	vector<int> b2cn(nv);
	vector<array<int, 2>> b1cn(nv);
	for (int i = 0; i < nv; i++)
	{
		b2ed[i] = 2 * i + 1;
		b2cn[i] = 2 * i + 2;
	}
	b1ed[0][0] = 6;
	b1ed[0][1] = 8;
	b1ed[1][0] = 3;
	b1ed[1][1] = 5;
	b1ed[2][0] = 1;
	b1ed[2][1] = 3 * nv;
	b1cn[0][0] = 5;
	b1cn[0][1] = 6;
	b1cn[1][0] = 1;
	b1cn[1][1] = 3;
	b1cn[2][0] = 3 * nv - 1;
	b1cn[2][1] = 3 * nv;
	for (int i = 2; i < nv - 1; i++)
	{
		int i1 = (2 * nv - 1 - 2 * (i - 2) - 1) / 2;
		b1ed[i1][0] = 3 * i + 3;
		b1ed[i1][1] = 3 * i + 5;
		i1 = (2 * nv - 2 * (i - 2)) / 2 - 1;
		b1cn[i1][0] = 3 * i + 2;
		b1cn[i1][1] = 3 * i + 3;
	}
	double coef[2] = {4. / 9., 2. / 9.};

	array<MatrixXd, 4> bemat22;
	for (int is = 0; is < 4; is++)
		bemat22[is] = MatrixXd::Zero(2 * nv + 1, 16);
	for (int ip = 0; ip < nv; ip++)
	{
		for (int is = 0; is < 4; is++)
		{
			for (int i = 0; i < 16; i++)
			{
				bemat22[is](b2ed[ip], i) = coef[0] * (tmesh[eid].c1mat22[is][b1ed[ip][0]][i] + tmesh[eid].c1mat22[is][b1ed[ip][1]][i]);
				bemat22[is](b2cn[ip], i) = coef[1] * (tmesh[eid].c1mat22[is][b1cn[ip][0]][i] + tmesh[eid].c1mat22[is][b1cn[ip][1]][i]);
			}
		}
	}
	for (int is = 0; is < 4; is++)
	{
		for (int i = 0; i < 16; i++)
		{
			bemat22[is](3, i) += coef[1] * tmesh[eid].c1mat22[is][2][i];
			bemat22[is](4, i) += coef[0] * tmesh[eid].c1mat22[is][2][i];
			bemat22[is](5, i) += coef[1] * tmesh[eid].c1mat22[is][2][i];
		}
	}

	vector<vector<double>> bemat0;
	SetBezierMatIrr_NonUniform(eid, bemat0); //only consider outermost 7
											 //	int cs1[3]={2,4,6};
											 //	int bzero1[3][4]={{0,4,-1,-1},{0,1,4,5},{0,1,-1,-1}};
											 //	int cs2[4]={1,3,5,7};
											 //	if(nv==3) cs2[3]=1;
											 //	int bzero2[4][8]={{0,4,8,12,-1,-1,-1,-1},{0,1,2,3,4,5,6,7},{0,4,8,12,1,5,9,13},{0,1,2,3,-1,-1,-1,-1}};
											 //	for(int i=0; i<3; i++)
											 //	{
											 //		for(int j=0; j<4; j++)
											 //		{
											 //			if(bzero1[i][j] != -1)
											 //			{
											 //				bemat0[cs1[i]][bzero1[i][j]]=0.;
											 //			}
											 //		}
											 //	}
											 //	for(int i=0; i<4; i++)
											 //	{
											 //		for(int j=0; j<8; j++)
											 //		{
											 //			if(bzero2[i][j] != -1)
											 //			{
											 //				bemat0[cs2[i]][bzero2[i][j]]=0.;
											 //			}
											 //		}
											 //	}

	MatrixXd bemat1 = MatrixXd::Zero(7, 16);
	for (int i = 0; i < 7; i++)
	{
		int i1(2 * nv + 1 + i);
		for (int j = 0; j < 16; j++)
			bemat1(i, j) = bemat0[i1][j];
	}
	MatrixXd mtmp0 = MatrixXd::Zero(5, 16);

	mtmp0(0, 12) = bemat0[2][12] - coef[1] * .25;
	int bzid3[4] = {12, 13, 14, 15};
	double rep3[4] = {2. * coef[0] * 0.25, coef[0] * 0.5, coef[1] * 0.5, coef[1] * 0.25};
	int bzid4[7] = {12, 13, 14, 15, 11, 7, 3};
	double rep4[7] = {coef[1] * 0.25, coef[1] * 0.5, coef[0] * 0.5, coef[0] * 0.25, coef[0] * 0.5, coef[1] * 0.5, coef[1] * 0.25};
	int bzid5[4] = {3, 7, 11, 15};
	double rep5[4] = {2. * coef[0] * 0.25, coef[0] * 0.5, coef[1] * 0.5, coef[1] * 0.25};
	for (int j = 0; j < 4; j++)
	{
		mtmp0(1, bzid3[j]) = bemat0[3][bzid3[j]] - rep3[j];
		mtmp0(3, bzid5[j]) = bemat0[5][bzid5[j]] - rep5[j];
	}
	for (int j = 0; j < 7; j++)
	{
		mtmp0(2, bzid4[j]) = bemat0[4][bzid4[j]] - rep4[j];
	}
	mtmp0(4, 3) = bemat0[6][3] - coef[1] * .25;

	MatrixXd smat = MatrixXd::Zero(16, 49);
	double smat1d[7][4] = {{1., 0., 0., 0.}, {.5, .5, 0., 0.}, {.25, .5, .25, 0.}, {.125, .375, .375, .125}, {0., .25, .5, .25}, {0., 0., .5, .5}, {0., 0., 0., 1.}};
	int loc(0);
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			int loc0(0);
			for (int i0 = 0; i0 < 4; i0++)
			{
				for (int j0 = 0; j0 < 4; j0++)
				{
					smat(loc0, loc) = smat1d[j][j0] * smat1d[i][i0];
					loc0++;
				}
			}
			loc++;
		}
	}
	MatrixXd bemat = bemat1 * smat;
	MatrixXd mtmp = mtmp0 * smat;

	int bzcnct[4][16] = {{0, 1, 2, 3, 7, 8, 9, 10, 14, 15, 16, 17, 21, 22, 23, 24}, {3, 4, 5, 6, 10, 11, 12, 13, 17, 18, 19, 20, 24, 25, 26, 27}, {21, 22, 23, 24, 28, 29, 30, 31, 35, 36, 37, 38, 42, 43, 44, 45}, {24, 25, 26, 27, 31, 32, 33, 34, 38, 39, 40, 41, 45, 46, 47, 48}};
	for (int is = 0; is < 4; is++)
	{
		for (int id = 0; id < 5; id++)
		{
			for (int i = 0; i < 16; i++)
			{
				bemat22[is](id + 2, i) += mtmp(id, bzcnct[is][i]);
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		tmesh[eid].bemat22[i].resize(enpta, vector<double>(16, 0.));
		loc = 0;
		//for (int j = 0; j < tmesh[eid].IEN.size(); j++)
		for (int j = 0; j < 2 * nv + 1; j++)
		{
			if (cp[tmesh[eid].IEN[j]].act == 1)
			{
				for (int k = 0; k < 16; k++)
				{
					tmesh[eid].bemat22[i][loc][k] = bemat22[i](j, k);
				}
				loc++;
			}
		}
		for (int j = 0; j < 7; j++)
		{
			if (cp[tmesh[eid].IEN[2 * nv + 1 + j]].act == 1)
			{
				for (int k = 0; k < 16; k++)
				{
					tmesh[eid].bemat22[i][loc][k] = bemat(j, bzcnct[i][k]);
				}
				loc++;
			}
		}
	}
}

void TruncatedTspline::FindSplineC2_Scale_Trs(int eid)
{
	//c1 non-irregular or transition element
	int enpta(0);
	for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
	{
		if (cp[tmesh[eid].IEN[i]].act == 1)
			enpta++;
	}
	if (enpta == 0)
		return;

	tmesh[eid].bemat.clear();
	tmesh[eid].IENc2.clear();
	for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
	{
		if (cp[tmesh[eid].IEN[j]].act == 1)
		{
			tmesh[eid].IENc2.push_back(tmesh[eid].IEN[j]);
		}
	}

	vector<BezierElement> bzmesh;
	BezierElementExtract_Unstruct(eid, bzmesh);

	int ednb[4][2] = {{-1, -1}, {-1, -1}, {-1, -1}, {-1, -1}};
	int cnnb[4][2] = {{-1, -1}, {-1, -1}, {-1, -1}, {-1, -1}};
	vector<int> ien(16, -1);
	vector<int> ienloc(16, -1);
	for (int i = 0; i < 4; i++)
	{
		ednb[i][0] = tmedge[tmesh[eid].edge[i]].face[0];
		if (ednb[i][0] == eid)
			ednb[i][0] = tmedge[tmesh[eid].edge[i]].face[1];
		int *it = find(tmesh[ednb[i][0]].cnct, tmesh[ednb[i][0]].cnct + 4, tmesh[eid].cnct[i]);
		ednb[i][1] = it - tmesh[ednb[i][0]].cnct;
	}
	for (int i = 0; i < 4; i++)
	{
		for (uint j = 0; j < cp[tmesh[eid].cnct[i]].face.size(); j++)
		{
			int fcid(cp[tmesh[eid].cnct[i]].face[j]);
			if (fcid != eid && fcid != ednb[i][0] && fcid != ednb[(i + 3) % 4][0])
			{
				cnnb[i][0] = fcid;
				int *it = find(tmesh[cnnb[i][0]].cnct, tmesh[cnnb[i][0]].cnct + 4, tmesh[eid].cnct[i]);
				cnnb[i][1] = it - tmesh[cnnb[i][0]].cnct;
				break;
			}
		}
	}
	int imap[4][4] = {{5, 4, 0, 1}, {6, 2, 3, 7}, {10, 11, 15, 14}, {9, 13, 12, 8}};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			ien[imap[i][j]] = tmesh[cnnb[i][0]].cnct[(cnnb[i][1] + j) % 4];
		}
	}
	for (uint i = 0; i < ien.size(); i++)
	{
		if (cp[ien[i]].act == 1)
		{
			vector<int>::iterator it = find(tmesh[eid].IEN.begin(), tmesh[eid].IEN.end(), ien[i]);
			ienloc[i] = it - tmesh[eid].IEN.begin();
		}
	}

	//	int sfc[4]={5,6,10,9};
	int sed[8] = {1, 2, 7, 11, 14, 13, 8, 4};
	int scn[4] = {0, 3, 15, 12};
	int crsfc1[4][4] = {{10, 11, 14, 15}, {8, 9, 12, 13}, {0, 1, 4, 5}, {2, 3, 6, 7}};
	int crsfc2[4][2][8] = {{{2, 6, 10, 14, 3, 7, 11, 15}, {8, 9, 10, 11, 12, 13, 14, 15}}, {{8, 9, 10, 11, 12, 13, 14, 15}, {0, 4, 8, 12, 1, 5, 9, 13}}, {{0, 4, 8, 12, 1, 5, 9, 13}, {0, 1, 2, 3, 4, 5, 6, 7}}, {{0, 1, 2, 3, 4, 5, 6, 7}, {2, 6, 10, 14, 3, 7, 11, 15}}};
	int crsed1[8][2] = {{2, 3}, {1, 0}, {11, 15}, {7, 3}, {13, 12}, {14, 15}, {4, 0}, {8, 12}};
	int crsed2[4][4] = {{0, 1, 2, 3}, {3, 7, 11, 15}, {12, 13, 14, 15}, {0, 4, 8, 12}};
	int crscn[4] = {0, 3, 15, 12};

	if (tmesh[eid].c1 == 1)
	{
		for (int i = 0; i < 4; i++)
		{
			if (cp[tmesh[eid].cnct[i]].act == 1)
			{
				if (tmesh[ednb[i][0]].c1 == 1)
				{
					for (int j = 0; j < 8; j++)
					{
						bzmesh[0].cmat[i][crsfc2[i][0][j]] = 0.;
					}
				}
				else if (tmesh[ednb[(i + 3) % 4][0]].c1 == 1)
				{
					for (int j = 0; j < 8; j++)
					{
						bzmesh[0].cmat[i][crsfc2[i][1][j]] = 0.;
					}
				}
				else
				{
					for (int j = 0; j < 4; j++)
					{
						bzmesh[0].cmat[i][crsfc1[i][j]] = 0.;
					}
				}
			}
		}
	}
	for (int i = 0; i < 4; i++)
	{
		if (tmesh[ednb[i][0]].c1 == 1)
		{
			if (cp[ien[sed[2 * i]]].act == 1)
			{
				if (tmesh[cnnb[i][0]].c1 == 1)
				{
					for (int j = 0; j < 4; j++)
					{
						bzmesh[0].cmat[ienloc[sed[2 * i]]][crsed2[i][j]] = 0.;
					}
				}
				else
				{
					bzmesh[0].cmat[ienloc[sed[2 * i]]][crsed1[2 * i][0]] = 0.;
					bzmesh[0].cmat[ienloc[sed[2 * i]]][crsed1[2 * i][1]] = 0.;
				}
			}
			if (cp[ien[sed[2 * i + 1]]].act == 1)
			{
				if (tmesh[cnnb[(i + 1) % 4][0]].c1 == 1)
				{
					for (int j = 0; j < 4; j++)
					{
						bzmesh[0].cmat[ienloc[sed[2 * i + 1]]][crsed2[i][j]] = 0.;
					}
				}
				else
				{
					bzmesh[0].cmat[ienloc[sed[2 * i + 1]]][crsed1[2 * i + 1][0]] = 0.;
					bzmesh[0].cmat[ienloc[sed[2 * i + 1]]][crsed1[2 * i + 1][1]] = 0.;
				}
			}
		}
	}
	for (int i = 0; i < 4; i++)
	{
		if (tmesh[cnnb[i][0]].c1 == 1)
		{
			if (cp[ien[scn[i]]].act == 1)
			{
				bzmesh[0].cmat[ienloc[scn[i]]][crscn[i]] = 0.;
			}
		}
	}

	tmesh[eid].bemat.resize(enpta, vector<double>(16, 0.));
	int loc = 0;
	for (uint j = 0; j < tmesh[eid].IEN.size(); j++)
	{
		if (cp[tmesh[eid].IEN[j]].act == 1)
		{
			for (int k = 0; k < 16; k++)
			{
				tmesh[eid].bemat[loc][k] = bzmesh[0].cmat[j][k];
			}
			loc++;
		}
	}
}

void TruncatedTspline::FindSplineC1_Scale_Trs(int eid)
{
	int ednb[4], cnnb[4];
	for (int i = 0; i < 4; i++)
	{
		ednb[i] = tmedge[tmesh[eid].edge[i]].face[0];
		if (ednb[i] == eid)
			ednb[i] = tmedge[tmesh[eid].edge[i]].face[1];
	}
	for (int i = 0; i < 4; i++)
	{
		cnnb[i] = -1;
		for (uint j = 0; j < cp[tmesh[eid].cnct[i]].face.size(); j++)
		{
			int fcid(cp[tmesh[eid].cnct[i]].face[j]);
			if (fcid != eid && fcid != ednb[i] && fcid != ednb[(i + 3) % 4])
			{
				cnnb[i] = fcid;
				break;
			}
		}
		if (cnnb[i] == -1)
		{
			cerr << "Can't find corner neighbor!\n";
			getchar();
		}
	}

	int cnc0[4] = {0, 3, 15, 12};
	//	int cnc1[4][4] = { { 0,1,4,5 },{ 2,3,6,7 },{ 10,11,14,15 },{ 8,9,12,13 } };
	int edc0[8] = {1, 2, 7, 11, 14, 13, 8, 4};
	//	int edc1[8][2] = { { 1,5 },{ 2,6 },{ 6,7 },{ 10,11 },{ 10,14 },{ 9,13 },{ 8,9 },{ 4,5 } };
	int fcc0[4] = {5, 6, 10, 9};
	//	int fcc1[4] = { 5,6,10,9 };
	int fcbz[4][4] = {{5, 1, 4, 0}, {6, 2, 7, 3}, {10, 14, 11, 15}, {9, 13, 8, 12}};
	int edbz[8][2] = {{1, 0}, {2, 3}, {7, 3}, {11, 15}, {14, 15}, {13, 12}, {8, 12}, {4, 0}};

	double ki[2][3];
	ki[0][1] = tmedge[tmesh[eid].edge[0]].len;
	ki[1][1] = tmedge[tmesh[eid].edge[1]].len;
	double kitmp[4];
	for (int i = 0; i < 4; i++)
	{
		int *it = find(tmesh[ednb[i]].edge, tmesh[ednb[i]].edge + 4, tmesh[eid].edge[i]);
		int loc(it - tmesh[ednb[i]].edge);
		kitmp[i] = tmedge[tmesh[ednb[i]].edge[(loc + 1) % 4]].len;
	}
	ki[0][0] = kitmp[3];
	ki[0][2] = kitmp[1];
	ki[1][0] = kitmp[0];
	ki[1][2] = kitmp[2];
	double cu[4] = {ki[0][1] / (ki[0][0] + ki[0][1]), ki[0][0] / (ki[0][0] + ki[0][1]),
					ki[0][2] / (ki[0][1] + ki[0][2]), ki[0][1] / (ki[0][1] + ki[0][2])};
	double cv[4] = {ki[1][1] / (ki[1][0] + ki[1][1]), ki[1][0] / (ki[1][0] + ki[1][1]),
					ki[1][2] / (ki[1][1] + ki[1][2]), ki[1][1] / (ki[1][1] + ki[1][2])};
	double cfc[4][4] = {{1., cv[1], cu[1], cu[1] * cv[1]}, {1., cv[1], cu[2], cu[2] * cv[1]}, {1., cv[2], cu[2], cu[2] * cv[2]}, {1., cv[2], cu[1], cu[1] * cv[2]}};
	double ced[8][2] = {{cv[0], cu[1] * cv[0]}, {cv[0], cu[2] * cv[0]}, {cu[3], cu[3] * cv[1]}, {cu[3], cu[3] * cv[2]}, {cv[3], cu[2] * cv[3]}, {cv[3], cu[1] * cv[3]}, {cu[0], cu[0] * cv[2]}, {cu[0], cu[0] * cv[1]}};
	double ccn[4] = {cu[0] * cv[0], cu[3] * cv[0], cu[3] * cv[3], cu[0] * cv[3]};

	if (tmesh[eid].c1 == 1 && tmesh[eid].IENc1.size() == 4 && tmesh[eid].c1mat.size() == 0)
	{
		for (int i = 0; i < 4; i++)
		{
			vector<double> c1tmp(16, 0.);
			for (int j = 0; j < 4; j++)
			{
				c1tmp[fcbz[i][j]] = cfc[i][j];
			}
			tmesh[eid].c1mat.push_back(c1tmp);
		}
	}

	for (int i = 0; i < 4; i++)
	{
		if (tmesh[cnnb[i]].c1 == 1)
		{
			int *it = find(tmesh[cnnb[i]].cnct, tmesh[cnnb[i]].cnct + 4, tmesh[eid].cnct[i]);
			int c1loc = it - tmesh[cnnb[i]].cnct;
			tmesh[eid].IENc1.push_back(tmesh[cnnb[i]].IENc1[c1loc]);
			vector<double> c1tmp(16, 0.);
			c1tmp[cnc0[i]] = ccn[i]; //0.25
			tmesh[eid].c1mat.push_back(c1tmp);
		}
		if (tmesh[ednb[i]].c1 == 1)
		{
			int *it = find(tmesh[ednb[i]].cnct, tmesh[ednb[i]].cnct + 4, tmesh[eid].cnct[i]);
			int c1loc = it - tmesh[ednb[i]].cnct;
			tmesh[eid].IENc1.push_back(tmesh[ednb[i]].IENc1[c1loc]);
			tmesh[eid].IENc1.push_back(tmesh[ednb[i]].IENc1[(c1loc + 3) % 4]);
			vector<double> c1tmp1(16, 0.), c1tmp2(16, 0.);
			c1tmp1[edbz[2 * i][0]] = ced[2 * i][0];
			c1tmp1[edbz[2 * i][1]] = ced[2 * i][1];
			tmesh[eid].c1mat.push_back(c1tmp1);
			c1tmp2[edbz[2 * i + 1][0]] = ced[2 * i + 1][0];
			c1tmp2[edbz[2 * i + 1][1]] = ced[2 * i + 1][1];
			tmesh[eid].c1mat.push_back(c1tmp2);
		}
	}
}

void TruncatedTspline::BezierUnit_DPatch_Irr22_Scale(int eid, vector<BezierElement> &bzmesh)
{
	vector<int> pid = tmesh[eid].IENc2;
	for (uint i = 0; i < tmesh[eid].IENc1.size(); i++)
	{
		pid.push_back(cp.size() + tmesh[eid].IENc1[i]);
	}
	for (int is = 0; is < 4; is++)
	{
		BezierElement bzel;
		bzel.prt = eid;
		bzel.cmat.resize(pid.size());
		bzel.IEN = pid;
		for (uint j = 0; j < tmesh[eid].IENc2.size(); j++)
		{
			for (int k = 0; k < 16; k++)
			{
				bzel.cmat[j][k] = tmesh[eid].bemat22[is][j][k];
				bzel.pts[k][0] += bzel.cmat[j][k] * cp[pid[j]].coor[0];
				bzel.pts[k][1] += bzel.cmat[j][k] * cp[pid[j]].coor[1];
				bzel.pts[k][2] += bzel.cmat[j][k] * cp[pid[j]].coor[2];
			}
		}
		int ist(tmesh[eid].IENc2.size());
		for (uint j = 0; j < tmesh[eid].IENc1.size(); j++)
		{
			int jloc(ist + j);
			for (int k = 0; k < 16; k++)
			{
				bzel.cmat[jloc][k] = tmesh[eid].c1mat22[is][j][k] * wc1[tmesh[eid].IENc1[j]];
				bzel.pts[k][0] += bzel.cmat[jloc][k] * bzcp[tmesh[eid].IENc1[j]][0];
				bzel.pts[k][1] += bzel.cmat[jloc][k] * bzcp[tmesh[eid].IENc1[j]][1];
				bzel.pts[k][2] += bzel.cmat[jloc][k] * bzcp[tmesh[eid].IENc1[j]][2];
			}
		}
		bzmesh.push_back(bzel);
	}
}

void TruncatedTspline::BezierUnit_DPatch_Trs_Scale(int eid, vector<BezierElement> &bzmesh)
{
	vector<int> pid = tmesh[eid].IENc2;
	for (uint i = 0; i < tmesh[eid].IENc1.size(); i++)
	{
		pid.push_back(cp.size() + tmesh[eid].IENc1[i]);
	}
	BezierElement bzel;
	bzel.prt = eid;
	bzel.cmat.resize(pid.size());
	bzel.IEN = pid;
	for (uint j = 0; j < tmesh[eid].IENc2.size(); j++)
	{
		for (int k = 0; k < 16; k++)
		{
			bzel.cmat[j][k] = tmesh[eid].bemat[j][k];
			bzel.pts[k][0] += bzel.cmat[j][k] * cp[pid[j]].coor[0];
			bzel.pts[k][1] += bzel.cmat[j][k] * cp[pid[j]].coor[1];
			bzel.pts[k][2] += bzel.cmat[j][k] * cp[pid[j]].coor[2];
		}
	}
	int ist(tmesh[eid].IENc2.size());
	for (uint j = 0; j < tmesh[eid].IENc1.size(); j++)
	{
		int jloc(ist + j);
		for (int k = 0; k < 16; k++)
		{
			bzel.cmat[jloc][k] = tmesh[eid].c1mat[j][k] * wc1[tmesh[eid].IENc1[j]];
			bzel.pts[k][0] += bzel.cmat[jloc][k] * bzcp[tmesh[eid].IENc1[j]][0];
			bzel.pts[k][1] += bzel.cmat[jloc][k] * bzcp[tmesh[eid].IENc1[j]][1];
			bzel.pts[k][2] += bzel.cmat[jloc][k] * bzcp[tmesh[eid].IENc1[j]][2];
		}
	}
	bzmesh.push_back(bzel);
}

void TruncatedTspline::BezierExtract_DPatch_AS_Scale(vector<BezierElement> &bzmesh) //Analysis space
{
	cout << "Bezier extracting...\n";
	bzmesh.clear();
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1)
		{
			if (tmesh[eid].c1 == 0 && (tmesh[eid].type == 0 || tmesh[eid].type == 1))
			{
				BezierElementExtract_Unstruct(eid, bzmesh);
			}
			else if (tmesh[eid].c1 == 1 && tmesh[eid].type == 4)
			{
				BezierUnit_DPatch_Irr22_Scale(eid, bzmesh);
			}
			else if (tmesh[eid].c1 == 1 && tmesh[eid].type != 4)
			{
				BezierUnit_DPatch_Trs_Scale(eid, bzmesh);
			}
			else if (tmesh[eid].c1 == 2)
			{
				BezierUnit_DPatch_Trs_Scale(eid, bzmesh);
			}
		}
	}
}

void TruncatedTspline::VisualizeSurface_Scale(string fn, int ns)
{
	vector<array<double, 3>> spt;
	vector<array<double, 3>> sval;
	vector<array<int, 4>> sele;
	vector<array<double, 3>> lpt; //visulize parameter lines
	vector<array<int, 2>> led;	  //line connectivity
	int ecount(0), loc0, loc1, loc2;

	for (uint e = 0; e < tmesh.size(); e++)
	{
		if (tmesh[e].act == 1 && (tmesh[e].type == 0 || tmesh[e].type == 1 || tmesh[e].type == 4))
		{
			vector<double> su(ns), sv(ns);
			for (int i = 0; i < ns; i++)
			{
				su[i] = i * tmedge[tmesh[e].edge[0]].len / (ns - 1);
				sv[i] = i * tmedge[tmesh[e].edge[3]].len / (ns - 1);
			}

			for (int a = 0; a < ns; a++)
			{
				for (int b = 0; b < ns; b++)
				{
					array<double, 3> pt;
					GeomMap_DPatch_Scale(e, su[b], sv[a], pt);
					spt.push_back(pt);
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
					el[0] = ecount * ns * ns + a * ns + b;
					el[1] = ecount * ns * ns + a * ns + b + 1;
					el[2] = ecount * ns * ns + (a + 1) * ns + b + 1;
					el[3] = ecount * ns * ns + (a + 1) * ns + b;
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
	}

	string fname = fn + "_geom.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i < spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
		}
		//fout<<"\nPOINT_DATA "<<sval.size()<<"\nSCALARS sum FLOAT\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<sval.size();i++)
		//{
		//	fout<<sval[i]<<"\n";
		//}
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

	string fname1(fn + "_geom-lines.vtk");
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

void TruncatedTspline::GeomMap_DPatch_Scale(int eid, double u, double v, array<double, 3> &pt)
{
	vector<double> Nt;
	vector<array<double, 2>> dNdt;
	ElementBasis_DPatch_Scale(eid, u, v, Nt, dNdt);
	pt[0] = 0.;
	pt[1] = 0.;
	pt[2] = 0.;
	double sum(0.);
	if (tmesh[eid].c1 == 1 || tmesh[eid].c1 == 2)
	{
		for (uint i = 0; i < tmesh[eid].IENc2.size(); i++)
		{
			pt[0] += cp[tmesh[eid].IENc2[i]].coor[0] * Nt[i];
			pt[1] += cp[tmesh[eid].IENc2[i]].coor[1] * Nt[i];
			pt[2] += cp[tmesh[eid].IENc2[i]].coor[2] * Nt[i];
			sum += Nt[i];
		}
		int ist(tmesh[eid].IENc2.size());
		for (uint i = 0; i < tmesh[eid].IENc1.size(); i++)
		{
			pt[0] += bzcp[tmesh[eid].IENc1[i]][0] * Nt[ist + i];
			pt[1] += bzcp[tmesh[eid].IENc1[i]][1] * Nt[ist + i];
			pt[2] += bzcp[tmesh[eid].IENc1[i]][2] * Nt[ist + i];
			sum += Nt[ist + i];
		}
	}
	else
	{
		for (uint i = 0; i < tmesh[eid].IEN.size(); i++)
		{
			pt[0] += cp[tmesh[eid].IEN[i]].coor[0] * Nt[i];
			pt[1] += cp[tmesh[eid].IEN[i]].coor[1] * Nt[i];
			pt[2] += cp[tmesh[eid].IEN[i]].coor[2] * Nt[i];
			sum += Nt[i];
		}
	}
}

void TruncatedTspline::ElementBasis_DPatch_Scale(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	if (tmesh[eid].act == 1)
	{
		if (tmesh[eid].c1 == 0 && (tmesh[eid].type == 0 || tmesh[eid].type == 1))
		{
			ElementBasis_Regular(eid, u, v, Nt, dNdt);
		}
		else if (tmesh[eid].c1 == 1 && tmesh[eid].type == 4)
		{
			ElementBasis_DPatch_Irr22_Scale(eid, u, v, Nt, dNdt);
		}
		else if (tmesh[eid].c1 == 1 && tmesh[eid].type != 4)
		{
			ElementBasis_DPatch_Trs_Scale(eid, u, v, Nt, dNdt);
		}
		else if (tmesh[eid].c1 == 2)
		{
			ElementBasis_DPatch_Trs_Scale(eid, u, v, Nt, dNdt);
		}
	}
}

void TruncatedTspline::ElementBasis_DPatch_Irr22_Scale(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	//need to change variabe for Bezier, (u,v) -> (u_b,v_b)
	Nt.clear();
	dNdt.clear();
	Nt.resize(tmesh[eid].IENc2.size() + tmesh[eid].IENc1.size(), 0.);
	dNdt.resize(tmesh[eid].IENc2.size() + tmesh[eid].IENc1.size());
	//uint nv(cp[tmesh[eid].cnct[0]].face.size());
	BezierElement be;
	double Nt0[16], dNdt0[16][2];
	//double u_b(u / tmedge[tmesh[eid].edge[0]].len), v_b(v / tmedge[tmesh[eid].edge[3]].len);
	//be.Basis(u_b, v_b, Nt0, dNdt0);
	double edm[2] = {tmedge[tmesh[eid].edge[0]].len, tmedge[tmesh[eid].edge[3]].len};
	double u_b(u / edm[0]), v_b(v / edm[1]);
	double uhalf(.5), vhalf(.5);
	int subeid(0);
	if (u_b >= uhalf && v_b < vhalf)
	{
		subeid = 1;
		u_b = 2. * (u_b - uhalf);
		v_b = 2. * v_b;
	}
	else if (u_b < uhalf && v_b >= vhalf)
	{
		subeid = 2;
		u_b = 2. * u_b;
		v_b = 2. * (v_b - vhalf);
	}
	else if (u_b >= uhalf && v_b >= vhalf)
	{
		subeid = 3;
		u_b = 2. * (u_b - uhalf);
		v_b = 2. * (v_b - vhalf);
	}
	else
	{
		subeid = 0;
		u_b = 2. * u_b;
		v_b = 2. * v_b;
	}
	be.Basis(u_b, v_b, Nt0, dNdt0);

	double ctmp;
	for (uint i = 0; i < tmesh[eid].IENc2.size(); i++)
	{
		dNdt[i][0] = 0.;
		dNdt[i][1] = 0.;
		for (int j = 0; j < 16; j++)
		{
			ctmp = tmesh[eid].bemat22[subeid][i][j];
			Nt[i] += ctmp * Nt0[j];
			dNdt[i][0] += ctmp * dNdt0[j][0];
			dNdt[i][1] += ctmp * dNdt0[j][1];
		}
		dNdt[i][0] *= 2. / edm[0];
		dNdt[i][1] *= 2. / edm[1];
	}
	int ist(tmesh[eid].IENc2.size());
	for (uint i = 0; i < tmesh[eid].IENc1.size(); i++)
	{
		int iloc(ist + i);
		dNdt[iloc][0] = 0.;
		dNdt[iloc][1] = 0.;
		for (int j = 0; j < 16; j++)
		{
			ctmp = tmesh[eid].c1mat22[subeid][i][j] * wc1[tmesh[eid].IENc1[i]];
			Nt[iloc] += ctmp * Nt0[j];
			dNdt[iloc][0] += ctmp * dNdt0[j][0];
			dNdt[iloc][1] += ctmp * dNdt0[j][1];
		}
		dNdt[iloc][0] *= 2. / edm[0];
		dNdt[iloc][1] *= 2. / edm[1];
	}
}

void TruncatedTspline::ElementBasis_DPatch_Trs_Scale(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt)
{
	//need to change variabe for Bezier, (u,v) -> (u_b,v_b)
	Nt.clear();
	dNdt.clear();
	Nt.resize(tmesh[eid].IENc2.size() + tmesh[eid].IENc1.size(), 0.);
	dNdt.resize(tmesh[eid].IENc2.size() + tmesh[eid].IENc1.size());
	//uint nv(cp[tmesh[eid].cnct[0]].face.size());
	BezierElement be;
	double Nt0[16], dNdt0[16][2];
	double u_b(u / tmedge[tmesh[eid].edge[0]].len), v_b(v / tmedge[tmesh[eid].edge[3]].len);
	be.Basis(u_b, v_b, Nt0, dNdt0);
	double ctmp;
	for (uint i = 0; i < tmesh[eid].IENc2.size(); i++)
	{
		dNdt[i][0] = 0.;
		dNdt[i][1] = 0.;
		for (int j = 0; j < 16; j++)
		{
			ctmp = tmesh[eid].bemat[i][j];
			Nt[i] += ctmp * Nt0[j];
			dNdt[i][0] += ctmp * dNdt0[j][0];
			dNdt[i][1] += ctmp * dNdt0[j][1];
		}
		dNdt[i][0] /= tmedge[tmesh[eid].edge[0]].len;
		dNdt[i][1] /= tmedge[tmesh[eid].edge[3]].len;
	}
	int ist(tmesh[eid].IENc2.size());
	for (uint i = 0; i < tmesh[eid].IENc1.size(); i++)
	{
		int iloc(ist + i);
		dNdt[iloc][0] = 0.;
		dNdt[iloc][1] = 0.;
		for (int j = 0; j < 16; j++)
		{
			ctmp = tmesh[eid].c1mat[i][j] * wc1[tmesh[eid].IENc1[i]];
			Nt[iloc] += ctmp * Nt0[j];
			dNdt[iloc][0] += ctmp * dNdt0[j][0];
			dNdt[iloc][1] += ctmp * dNdt0[j][1];
		}
		dNdt[iloc][0] /= tmedge[tmesh[eid].edge[0]].len;
		dNdt[iloc][1] /= tmedge[tmesh[eid].edge[3]].len;
	}
}

void TruncatedTspline::EnlargeOneRing()
{
	vector<int> elist;
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1 && tmesh[eid].c1 == 0)
		{
			int flag(0);
			for (int i = 0; i < 4; i++)
			{
				for (int fcid : cp[tmesh[eid].cnct[i]].face)
				{
					if (tmesh[fcid].c1 == 1)
					{
						flag = 1;
						break;
					}
				}
				if (flag == 1)
					break;
			}
			if (flag == 1)
				elist.push_back(eid);
		}
	}

	int bzfc[4] = {5, 6, 10, 9};
	for (int eid : elist)
	{
		tmesh[eid].c1 = 1; //need to get all its Bezier points and add C1 ones to bzcp
		vector<BezierElement> bzmesh;
		BezierElementExtract_Unstruct(eid, bzmesh);
		tmesh[eid].bzpts.clear();
		tmesh[eid].bzpts.resize(16);
		for (int k = 0; k < 16; k++)
		{
			tmesh[eid].bzpts[k][0] = bzmesh[0].pts[k][0];
			tmesh[eid].bzpts[k][1] = bzmesh[0].pts[k][1];
			tmesh[eid].bzpts[k][2] = bzmesh[0].pts[k][2];
		}
		tmesh[eid].IENc1.clear();
		for (int k = 0; k < 4; k++)
		{
			bzcp.push_back(tmesh[eid].bzpts[bzfc[k]]);
			tmesh[eid].IENc1.push_back(bzcp.size() - 1);
		}
	}
}

void TruncatedTspline::FindSplineC1_Scale()
{
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].c1 == 0)
		{
			int flag(0);
			for (int j = 0; j < 4; j++)
			{
				for (uint k = 0; k < cp[tmesh[i].cnct[j]].face.size(); k++)
				{
					if (tmesh[cp[tmesh[i].cnct[j]].face[k]].c1 == 1)
					{
						flag = 1;
						break;
					}
				}
				if (flag == 1)
					break;
			}
			if (flag == 1)
				tmesh[i].c1 = 2; //transition
		}
	}
	for (uint i = 0; i < cp.size(); i++)
	{
		int flag(0);
		for (uint j = 0; j < cp[i].face.size(); j++)
		{
			if (tmesh[cp[i].face[j]].c1 != 1)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
		{
			cp[i].act = 0;
		}
		else
		{
			cp[i].act = 1;
		}
	}

	//set scaling coefficients
	wc1.clear();
	wc1.resize(bzcp.size(), 1.);
	double coef[2] = {1. / 3., 1. / 9.};
	double a(4. / 9.), b(2. / 9.);
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1 && tmesh[eid].c1 == 1)
		{
			int flag[4] = {0, 0, 0, 0};
			for (int i = 0; i < 4; i++)
			{
				for (int fcid : cp[tmesh[eid].cnct[i]].face)
				{
					if (tmesh[fcid].c1 == 2)
					{
						flag[i]++;
					}
				}
			}
			for (int i = 0; i < 4; i++)
			{
				if (flag[i] == 3)
				{
					int loc[3] = {i, (i + 1) % 4, (i + 3) % 4};
					int id[3] = {tmesh[eid].IENc1[loc[0]], tmesh[eid].IENc1[loc[1]], tmesh[eid].IENc1[loc[2]]};
					int cn[3] = {tmesh[eid].cnct[loc[0]], tmesh[eid].cnct[loc[1]], tmesh[eid].cnct[loc[2]]};
					wc1[id[0]] = coef[1];
					wc1[id[1]] = coef[0];
					wc1[id[2]] = coef[0];
					bzcp[id[0]][0] -= (a * cp[cn[0]].coor[0] + b * cp[cn[1]].coor[0] + b * cp[cn[2]].coor[0]);
					bzcp[id[0]][1] -= (a * cp[cn[0]].coor[1] + b * cp[cn[1]].coor[1] + b * cp[cn[2]].coor[1]);
					bzcp[id[0]][2] -= (a * cp[cn[0]].coor[2] + b * cp[cn[1]].coor[2] + b * cp[cn[2]].coor[2]);
					bzcp[id[1]][0] -= (a * cp[cn[1]].coor[0] + b * cp[cn[0]].coor[0]);
					bzcp[id[1]][1] -= (a * cp[cn[1]].coor[1] + b * cp[cn[0]].coor[1]);
					bzcp[id[1]][2] -= (a * cp[cn[1]].coor[2] + b * cp[cn[0]].coor[2]);
					bzcp[id[2]][0] -= (a * cp[cn[2]].coor[0] + b * cp[cn[0]].coor[0]);
					bzcp[id[2]][1] -= (a * cp[cn[2]].coor[1] + b * cp[cn[0]].coor[1]);
					bzcp[id[2]][2] -= (a * cp[cn[2]].coor[2] + b * cp[cn[0]].coor[2]);
					for (int j = 0; j < 3; j++)
					{
						bzcp[id[0]][j] /= coef[1];
						bzcp[id[1]][j] /= coef[0];
						bzcp[id[2]][j] /= coef[0];
					}
				}
				else if (flag[i] == 2 && flag[(i + 1) % 4] == 2)
				{
					int loc[2] = {i, (i + 1) % 4};
					int id[2] = {tmesh[eid].IENc1[loc[0]], tmesh[eid].IENc1[loc[1]]};
					int cn[2] = {tmesh[eid].cnct[loc[0]], tmesh[eid].cnct[loc[1]]};
					wc1[id[0]] = coef[0];
					wc1[id[1]] = coef[0];
					bzcp[id[0]][0] -= (a * cp[cn[0]].coor[0] + b * cp[cn[1]].coor[0]);
					bzcp[id[0]][1] -= (a * cp[cn[0]].coor[1] + b * cp[cn[1]].coor[1]);
					bzcp[id[0]][2] -= (a * cp[cn[0]].coor[2] + b * cp[cn[1]].coor[2]);
					bzcp[id[1]][0] -= (a * cp[cn[1]].coor[0] + b * cp[cn[0]].coor[0]);
					bzcp[id[1]][1] -= (a * cp[cn[1]].coor[1] + b * cp[cn[0]].coor[1]);
					bzcp[id[1]][2] -= (a * cp[cn[1]].coor[2] + b * cp[cn[0]].coor[2]);
					for (int j = 0; j < 3; j++)
					{
						bzcp[id[0]][j] /= coef[0];
						bzcp[id[1]][j] /= coef[0];
					}
				}
			}
		}
	}

	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1)
		{
			if (tmesh[eid].c1 == 1 && tmesh[eid].type == 4)
			{
				FindSplineC1_Irr(eid);
			}
			else if (tmesh[eid].c1 == 1 && tmesh[eid].type != 4)
			{
				FindSplineC1_Scale_Trs(eid);
			}
			else if (tmesh[eid].c1 == 2)
			{
				FindSplineC1_Scale_Trs(eid);
			}
		}
	}
}

void TruncatedTspline::OutputBezierMeshAll(const vector<BezierElement> &bzmesh, string fn)
{
	string fname = fn + "_bzmesh_all.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << 16 * bzmesh.size() << " float\n";
		for (uint i = 0; i < bzmesh.size(); i++)
		{
			for (int j = 0; j < 16; j++)
			{
				fout << bzmesh[i].pts[j][0] << " " << bzmesh[i].pts[j][1] << " " << bzmesh[i].pts[j][2] << "\n";
			}
		}
		fout << "\nCELLS " << 9 * bzmesh.size() << " " << 45 * bzmesh.size() << '\n';
		for (int eid = 0; eid < bzmesh.size(); eid++)
		{
			int ist(16 * eid);
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					fout << "4 " << ist + 4 * i + j << " " << ist + 4 * i + j + 1 << " " << ist + 4 * (i + 1) + j + 1 << " " << ist + 4 * (i + 1) + j << '\n';
				}
			}
		}
		fout << "\nCELL_TYPES " << 9 * bzmesh.size() << '\n';
		for (int i = 0; i < bzmesh.size(); i++)
		{
			for (int i = 0; i < 9; i++)
				fout << "9\n";
		}

		/*fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
        fout << "POINTS " << 16 * bzmesh.size() << " float\n";
        for (uint i = 0; i < bzmesh.size(); i++)
        {
            for (int j = 0; j < 16; j++)
            {
                fout << bzmesh[i].pts[j][0] << " " << bzmesh[i].pts[j][1] << " " << bzmesh[i].pts[j][2] << "\n";
            }
        }
        fout << "\nCELLS " << 9 * bzmesh.size() << " " << 45 * bzmesh.size() << '\n';
        for (int i = 0; i<bzmesh.size(); i++)
        {
            int ist(16 * i);
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    fout << "4 " << ist+4*j+k << " " << ist + 4 * j + k+1 << " " << ist + 4 * (j+1) + k + 1 << " " << ist + 4 * (j + 1) + k << '\n';
                }
            }
        }
        fout << "\nCELL_TYPES " << 9*bzmesh.size() << '\n';
        for (int i = 0; i<9*bzmesh.size(); i++)
        {
            fout << "9\n";
        }*/

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void TruncatedTspline::BezierExtract_DPatch_AS_Scale_Refine(vector<BezierElement> &bzmesh) //compare refinement
{
	cout << "Bezier extracting...\n";
	bzmesh.clear();
	int bzcnct[4][16] = {{0, 1, 2, 3, 7, 8, 9, 10, 14, 15, 16, 17, 21, 22, 23, 24}, {3, 4, 5, 6, 10, 11, 12, 13, 17, 18, 19, 20, 24, 25, 26, 27}, {21, 22, 23, 24, 28, 29, 30, 31, 35, 36, 37, 38, 42, 43, 44, 45}, {24, 25, 26, 27, 31, 32, 33, 34, 38, 39, 40, 41, 45, 46, 47, 48}};
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		if (tmesh[eid].act == 1)
		{
			vector<BezierElement> bzmesh0;
			if (tmesh[eid].c1 == 0 && (tmesh[eid].type == 0 || tmesh[eid].type == 1))
			{
				BezierElementExtract_Unstruct(eid, bzmesh0);
			}
			else if (tmesh[eid].c1 == 1 && tmesh[eid].type == 4)
			{
				BezierUnit_DPatch_Irr22_Scale(eid, bzmesh0);
			}
			else if (tmesh[eid].c1 == 1 && tmesh[eid].type != 4)
			{
				BezierUnit_DPatch_Trs_Scale(eid, bzmesh0);
			}
			else if (tmesh[eid].c1 == 2)
			{
				BezierUnit_DPatch_Trs_Scale(eid, bzmesh0);
			}

			if (bzmesh0.size() != 0)
			{
				vector<array<double, 3>> pts0(16), pts1;
				for (int i = 0; i < 16; i++)
				{
					pts0[i][0] = bzmesh0[0].pts[i][0];
					pts0[i][1] = bzmesh0[0].pts[i][1];
					pts0[i][2] = bzmesh0[0].pts[i][2];
				}
				BezierRefineBi3(pts0, pts1);
				for (int is = 0; is < 4; is++)
				{
					BezierElement bzel;
					for (int i = 0; i < 16; i++)
					{
						bzel.pts[i][0] = pts1[bzcnct[is][i]][0];
						bzel.pts[i][1] = pts1[bzcnct[is][i]][1];
						bzel.pts[i][2] = pts1[bzcnct[is][i]][2];
					}
					bzmesh.push_back(bzel);
				}
			}
			if (bzmesh0.size() == 4)
			{
				for (int i = 1; i < 4; i++)
					bzmesh.push_back(bzmesh0[i]);
			}
		}
	}
}

void TruncatedTspline::HemisphereLocalRefine()
{
	//refine on input
	int ref4[8] = {21, 22, 27, 28, 57, 58, 63, 64};
	int ref2[2] = {15, 56};
	int ref2_cn[2] = {25, 66};
	int nref4(8);
	int nref2(2);

	vector<int> rfid, rftype;
	for (int i = 0; i < nref4; i++)
	{
		rfid.push_back(ref4[i]);
		rftype.push_back(0);
		int eid(ref4[i]);
		for (int j = 0; j < 4; j++)
		{
			int enb(tmedge[tmesh[eid].edge[j]].face[0]);
			if (enb == eid)
				enb = tmedge[tmesh[eid].edge[j]].face[1];
			if (tmesh[enb].act == 1 && tmesh[enb].type == 2)
			{
				vector<int>::iterator it = find(rfid.begin(), rfid.end(), enb);
				if (it == rfid.end())
				{
					rfid.push_back(enb);
					rftype.push_back(3);
				}
			}
		}
	}
	for (int i = 0; i < nref2; i++)
	{
		rfid.push_back(ref2[i]);
		int type(1);
		int *it = find(tmesh[ref2[i]].cnct, tmesh[ref2[i]].cnct + 4, ref2_cn[i]);
		int pos = it - tmesh[ref2[i]].cnct;
		if (pos == 4)
		{
			cerr << "Can't find corner when refining type 2!\n";
			getchar();
		}
		if (pos == 1 || pos == 3)
			type = 2;
		rftype.push_back(type);

		int eid(ref2[i]);
		for (int j = 0; j < 4; j++)
		{
			int enb(tmedge[tmesh[eid].edge[j]].face[0]);
			if (enb == eid)
				enb = tmedge[tmesh[eid].edge[j]].face[1];
			if (tmesh[enb].act == 1 && tmesh[enb].type == 2)
			{
				int type1(1);
				if (j == 1 || j == 3)
					type1 = 2;
				vector<int>::iterator it1 = find(rfid.begin(), rfid.end(), enb);
				if (it1 == rfid.end() && type == type1)
				{
					rfid.push_back(enb);
					rftype.push_back(3);
				}
			}
		}
	}

	Refine_Unstruct(rfid, rftype);
}

void TruncatedTspline::HemisphereLocalRefine1()
{
	//2nd time refinement
	int ref4[8] = {126, 127, 125, 124, 149, 150, 148, 151};
	int ref2[2] = {115, 143};
	int ref2_cn[2] = {135, 74};
	int nref4(8);
	int nref2(2);

	vector<int> rfid, rftype;
	for (int i = 0; i < nref4; i++)
	{
		rfid.push_back(ref4[i]);
		rftype.push_back(0);
		int eid(ref4[i]);
		for (int j = 0; j < 4; j++)
		{
			int enb(tmedge[tmesh[eid].edge[j]].face[0]);
			if (enb == eid)
				enb = tmedge[tmesh[eid].edge[j]].face[1];
			if (tmesh[enb].act == 1 && tmesh[enb].type == 2)
			{
				vector<int>::iterator it = find(rfid.begin(), rfid.end(), enb);
				if (it == rfid.end())
				{
					rfid.push_back(enb);
					rftype.push_back(3);
				}
			}
		}
	}
	for (int i = 0; i < nref2; i++)
	{
		rfid.push_back(ref2[i]);
		int type(1);
		int *it = find(tmesh[ref2[i]].cnct, tmesh[ref2[i]].cnct + 4, ref2_cn[i]);
		int pos = it - tmesh[ref2[i]].cnct;
		if (pos == 4)
		{
			cerr << "Can't find corner when refining type 2!\n";
			getchar();
		}
		if (pos == 1 || pos == 3)
			type = 2;
		rftype.push_back(type);

		int eid(ref2[i]);
		for (int j = 0; j < 4; j++)
		{
			int enb(tmedge[tmesh[eid].edge[j]].face[0]);
			if (enb == eid)
				enb = tmedge[tmesh[eid].edge[j]].face[1];
			if (tmesh[enb].act == 1 && tmesh[enb].type == 2)
			{
				int type1(1);
				if (j == 1 || j == 3)
					type1 = 2;
				vector<int>::iterator it1 = find(rfid.begin(), rfid.end(), enb);
				if (it1 == rfid.end() && type == type1)
				{
					rfid.push_back(enb);
					rftype.push_back(3);
				}
			}
		}
	}

	Refine_Unstruct(rfid, rftype);
}

void TruncatedTspline::HemisphereLocalRefine2()
{
	//2nd time refinement
	int ref4[8] = {162, 163, 161, 160, 191, 192, 190, 193};
	int ref2[2] = {177, 203};
	int ref2_cn[2] = {170, 162};
	int nref4(8);
	int nref2(2);

	vector<int> rfid, rftype;
	for (int i = 0; i < nref4; i++)
	{
		rfid.push_back(ref4[i]);
		rftype.push_back(0);
		int eid(ref4[i]);
		for (int j = 0; j < 4; j++)
		{
			int enb(tmedge[tmesh[eid].edge[j]].face[0]);
			if (enb == eid)
				enb = tmedge[tmesh[eid].edge[j]].face[1];
			if (tmesh[enb].act == 1 && tmesh[enb].type == 2)
			{
				vector<int>::iterator it = find(rfid.begin(), rfid.end(), enb);
				if (it == rfid.end())
				{
					rfid.push_back(enb);
					rftype.push_back(3);
				}
			}
		}
	}
	for (int i = 0; i < nref2; i++)
	{
		rfid.push_back(ref2[i]);
		int type(1);
		int *it = find(tmesh[ref2[i]].cnct, tmesh[ref2[i]].cnct + 4, ref2_cn[i]);
		int pos = it - tmesh[ref2[i]].cnct;
		if (pos == 4)
		{
			cerr << "Can't find corner when refining type 2!\n";
			getchar();
		}
		if (pos == 1 || pos == 3)
			type = 2;
		rftype.push_back(type);

		int eid(ref2[i]);
		for (int j = 0; j < 4; j++)
		{
			int enb(tmedge[tmesh[eid].edge[j]].face[0]);
			if (enb == eid)
				enb = tmedge[tmesh[eid].edge[j]].face[1];
			if (tmesh[enb].act == 1 && tmesh[enb].type == 2)
			{
				int type1(1);
				if (j == 1 || j == 3)
					type1 = 2;
				vector<int>::iterator it1 = find(rfid.begin(), rfid.end(), enb);
				if (it1 == rfid.end() && type == type1)
				{
					rfid.push_back(enb);
					rftype.push_back(3);
				}
			}
		}
	}

	Refine_Unstruct(rfid, rftype);
}

void TruncatedTspline::HemisphereLocalRefine(const vector<int> &ref4, const vector<int> &ref2, const vector<int> &ref2_cn)
{
	int nref4(ref4.size());
	int nref2(ref2.size());

	vector<int> rfid, rftype;
	for (int i = 0; i < nref4; i++)
	{
		rfid.push_back(ref4[i]);
		rftype.push_back(0);
		int eid(ref4[i]);
		for (int j = 0; j < 4; j++)
		{
			int enb(tmedge[tmesh[eid].edge[j]].face[0]);
			if (enb == eid)
				enb = tmedge[tmesh[eid].edge[j]].face[1];
			if (tmesh[enb].act == 1 && tmesh[enb].type == 2)
			{
				vector<int>::iterator it = find(rfid.begin(), rfid.end(), enb);
				if (it == rfid.end())
				{
					rfid.push_back(enb);
					rftype.push_back(3);
				}
			}
		}
	}
	for (int i = 0; i < nref2; i++)
	{
		rfid.push_back(ref2[i]);
		int type(1);
		int *it = find(tmesh[ref2[i]].cnct, tmesh[ref2[i]].cnct + 4, ref2_cn[i]);
		int pos = it - tmesh[ref2[i]].cnct;
		if (pos == 4)
		{
			cerr << "Can't find corner when refining type 2!\n";
			getchar();
		}
		if (pos == 1 || pos == 3)
			type = 2;
		rftype.push_back(type);

		int eid(ref2[i]);
		for (int j = 0; j < 4; j++)
		{
			int enb(tmedge[tmesh[eid].edge[j]].face[0]);
			if (enb == eid)
				enb = tmedge[tmesh[eid].edge[j]].face[1];
			if (tmesh[enb].act == 1 && tmesh[enb].type == 2)
			{
				int type1(1);
				if (j == 1 || j == 3)
					type1 = 2;
				vector<int>::iterator it1 = find(rfid.begin(), rfid.end(), enb);
				if (it1 == rfid.end() && type == type1)
				{
					rfid.push_back(enb);
					rftype.push_back(3);
				}
			}
		}
	}

	Refine_Unstruct(rfid, rftype);
}

void TruncatedTspline::HemisphereLocalRefine1_6()
{
	string fpre("../io/hemisphere_adp1/lev");
	//refine on input
	string fn1(fpre + "1.txt");
	vector<int> rf4_1;
	ReadRef4ID(fn1, rf4_1);
	//int ref4_1[8] = { 454,455,453,452,213,214,212,215 };
	int ref2_1[2] = {161, 283};
	int ref2_1_cn[2] = {174, 263};
	//int nref4_1(8);
	int nref2_1(2);
	//vector<int> rf4_1(ref4_1, ref4_1 + nref4_1);
	vector<int> rf2_1(ref2_1, ref2_1 + nref2_1);
	vector<int> rf2_1_cn(ref2_1_cn, ref2_1_cn + nref2_1);
	HemisphereLocalRefine(rf4_1, rf2_1, rf2_1_cn);

	string fn2(fpre + "2.txt");
	vector<int> rf4_2;
	ReadRef4ID(fn2, rf4_2);
	//int ref4_2[8] = { 470,471,469,468,499,500,498,501 };
	int ref2_2[2] = {507, 625};
	int ref2_2_cn[2] = {434, 527};
	//int nref4_2(8);
	int nref2_2(2);
	//vector<int> rf4_2(ref4_2, ref4_2 + nref4_2);
	vector<int> rf2_2(ref2_2, ref2_2 + nref2_2);
	vector<int> rf2_2_cn(ref2_2_cn, ref2_2_cn + nref2_2);
	HemisphereLocalRefine(rf4_2, rf2_2, rf2_2_cn);

	string fn3(fpre + "3.txt");
	vector<int> rf4_3;
	ReadRef4ID(fn3, rf4_3);
	//int ref4_3[8] = { 522,523,521,520,551,552,550,553 };
	int ref2_3[2] = {825, 947};
	int ref2_3_cn[2] = {643, 776};
	//int nref4_3(8);
	int nref2_3(2);
	//vector<int> rf4_3(ref4_3, ref4_3 + nref4_3);
	vector<int> rf2_3(ref2_3, ref2_3 + nref2_3);
	vector<int> rf2_3_cn(ref2_3_cn, ref2_3_cn + nref2_3);
	HemisphereLocalRefine(rf4_3, rf2_3, rf2_3_cn);

	string fn4(fpre + "4.txt");
	vector<int> rf4_4;
	ReadRef4ID(fn4, rf4_4);
	//int ref4_4[8] = { 574,575,573,572,603,604,602,605 };
	int ref2_4[2] = {1103, 967};
	int ref2_4_cn[2] = {902, 790};
	//int nref4_4(8);
	int nref2_4(2);
	//vector<int> rf4_4(ref4_4, ref4_4 + nref4_4);
	vector<int> rf2_4(ref2_4, ref2_4 + nref2_4);
	vector<int> rf2_4_cn(ref2_4_cn, ref2_4_cn + nref2_4);
	HemisphereLocalRefine(rf4_4, rf2_4, rf2_4_cn);

	string fn5(fpre + "5.txt");
	vector<int> rf4_5;
	ReadRef4ID(fn5, rf4_5);
	//int ref4_5[8] = { 626,627,625,624,655,656,654,657 };
	int ref2_5[2] = {1231, 1353};
	int ref2_5_cn[2] = {1003, 1097};
	//int nref4_5(8);
	int nref2_5(2);
	//vector<int> rf4_5(ref4_5, ref4_5 + nref4_5);
	vector<int> rf2_5(ref2_5, ref2_5 + nref2_5);
	vector<int> rf2_5_cn(ref2_5_cn, ref2_5_cn + nref2_5);
	HemisphereLocalRefine(rf4_5, rf2_5, rf2_5_cn);
}

void TruncatedTspline::HemisphereLocalRefine2_6()
{
	string fpre("../io/hemisphere_adp2/lev");
	//refine on input
	string fn1(fpre + "1.txt");
	vector<int> rf4_1;
	ReadRef4ID(fn1, rf4_1);
	int ref2_1[2] = {675, 1121};
	int ref2_1_cn[2] = {566, 263};
	int nref2_1(2);
	vector<int> rf2_1(ref2_1, ref2_1 + nref2_1);
	vector<int> rf2_1_cn(ref2_1_cn, ref2_1_cn + nref2_1);
	HemisphereLocalRefine(rf4_1, rf2_1, rf2_1_cn);

	string fn2(fpre + "2.txt");
	vector<int> rf4_2;
	ReadRef4ID(fn2, rf4_2);
	int ref2_2[2] = {1823, 2309};
	int ref2_2_cn[2] = {1422, 973};
	int nref2_2(2);
	vector<int> rf2_2(ref2_2, ref2_2 + nref2_2);
	vector<int> rf2_2_cn(ref2_2_cn, ref2_2_cn + nref2_2);
	HemisphereLocalRefine(rf4_2, rf2_2, rf2_2_cn);

	string fn3(fpre + "3.txt");
	vector<int> rf4_3;
	ReadRef4ID(fn3, rf4_3);
	int ref2_3[2] = {2863, 3197};
	int ref2_3_cn[2] = {2231, 1889};
	int nref2_3(2);
	vector<int> rf2_3(ref2_3, ref2_3 + nref2_3);
	vector<int> rf2_3_cn(ref2_3_cn, ref2_3_cn + nref2_3);
	HemisphereLocalRefine(rf4_3, rf2_3, rf2_3_cn);

	string fn4(fpre + "4.txt");
	vector<int> rf4_4;
	ReadRef4ID(fn4, rf4_4);
	int ref2_4[2] = {3679, 4125};
	int ref2_4_cn[2] = {2857, 2541};
	int nref2_4(2);
	vector<int> rf2_4(ref2_4, ref2_4 + nref2_4);
	vector<int> rf2_4_cn(ref2_4_cn, ref2_4_cn + nref2_4);
	HemisphereLocalRefine(rf4_4, rf2_4, rf2_4_cn);

	string fn5(fpre + "5.txt");
	vector<int> rf4_5;
	ReadRef4ID(fn5, rf4_5);
	int ref2_5[2] = {4563, 5009};
	int ref2_5_cn[2] = {3539, 3276};
	int nref2_5(2);
	vector<int> rf2_5(ref2_5, ref2_5 + nref2_5);
	vector<int> rf2_5_cn(ref2_5_cn, ref2_5_cn + nref2_5);
	HemisphereLocalRefine(rf4_5, rf2_5, rf2_5_cn);
}

void TruncatedTspline::HemisphereLocalRefine(string fpre, int nlev)
{
	for (int it = 0; it < nlev; it++)
	{
		string fn4 = fpre + to_string(it + 1) + ".txt";
		string fn2 = fpre + to_string(it + 1) + "_bi.txt";
		vector<int> rf4, rf2, rf2_cn;
		ReadRef4ID(fn4, rf4);
		ReadRef2ID(fn2, rf2, rf2_cn);
		HemisphereLocalRefine(rf4, rf2, rf2_cn);
	}

	/*string fn1(fpre + "1.txt");
	vector<int> rf4_1;
	ReadRef4ID(fn1, rf4_1);
	int ref2_1[2] = { 675,1121 };
	int ref2_1_cn[2] = { 566,263 };
	int nref2_1(2);
	vector<int> rf2_1(ref2_1, ref2_1 + nref2_1);
	vector<int> rf2_1_cn(ref2_1_cn, ref2_1_cn + nref2_1);
	HemisphereLocalRefine(rf4_1, rf2_1, rf2_1_cn);

	string fn2(fpre + "2.txt");
	vector<int> rf4_2;
	ReadRef4ID(fn2, rf4_2);
	int ref2_2[2] = { 1823,2309 };
	int ref2_2_cn[2] = { 1422,973 };
	int nref2_2(2);
	vector<int> rf2_2(ref2_2, ref2_2 + nref2_2);
	vector<int> rf2_2_cn(ref2_2_cn, ref2_2_cn + nref2_2);
	HemisphereLocalRefine(rf4_2, rf2_2, rf2_2_cn);

	string fn3(fpre + "3.txt");
	vector<int> rf4_3;
	ReadRef4ID(fn3, rf4_3);
	int ref2_3[2] = { 2863,3197 };
	int ref2_3_cn[2] = { 2231,1889 };
	int nref2_3(2);
	vector<int> rf2_3(ref2_3, ref2_3 + nref2_3);
	vector<int> rf2_3_cn(ref2_3_cn, ref2_3_cn + nref2_3);
	HemisphereLocalRefine(rf4_3, rf2_3, rf2_3_cn);

	string fn4(fpre + "4.txt");
	vector<int> rf4_4;
	ReadRef4ID(fn4, rf4_4);
	int ref2_4[2] = { 3679,4125 };
	int ref2_4_cn[2] = { 2857,2541 };
	int nref2_4(2);
	vector<int> rf2_4(ref2_4, ref2_4 + nref2_4);
	vector<int> rf2_4_cn(ref2_4_cn, ref2_4_cn + nref2_4);
	HemisphereLocalRefine(rf4_4, rf2_4, rf2_4_cn);

	string fn5(fpre + "5.txt");
	vector<int> rf4_5;
	ReadRef4ID(fn5, rf4_5);
	int ref2_5[2] = { 4563,5009 };
	int ref2_5_cn[2] = { 3539,3276 };
	int nref2_5(2);
	vector<int> rf2_5(ref2_5, ref2_5 + nref2_5);
	vector<int> rf2_5_cn(ref2_5_cn, ref2_5_cn + nref2_5);
	HemisphereLocalRefine(rf4_5, rf2_5, rf2_5_cn);*/
}

void TruncatedTspline::ReadRef4ID(string fn, vector<int> &rfid)
{
	rfid.clear();
	ifstream fin;
	fin.open(fn);
	if (fin.is_open())
	{
		string stmp;
		int itmp;
		while (getline(fin, stmp))
		{
			stringstream ss(stmp);
			ss >> itmp;
			if (tmesh[itmp].act == 1 && tmesh[itmp].type == 0)
			{
				rfid.push_back(itmp);
			}
		}
		fin.close();
	}
	else
	{
		cerr << "Can't open " << fn << "!\n";
	}
}

void TruncatedTspline::ReadRef2ID(string fn, vector<int> &rfid, vector<int> &rfid_cn)
{
	rfid.clear();
	rfid_cn.clear();

	ifstream fin;
	fin.open(fn);
	if (fin.is_open())
	{
		string stmp;
		int itmp, cntmp;
		while (getline(fin, stmp))
		{
			stringstream ss(stmp);
			ss >> itmp >> cntmp;
			if (tmesh[itmp].act == 1 && (tmesh[itmp].type == 0 || tmesh[itmp].type == 1))
			{
				rfid.push_back(itmp);
				rfid_cn.push_back(cntmp);
			}
		}
		fin.close();
	}
	else
	{
		cerr << "Can't open " << fn << "!\n";
	}
}

void TruncatedTspline::SmoothPlanar(string fn_in, string fn_out, int nstep)
{
	InputFromVTK_Quad(fn_in);
	//RescaleCoor();
	InitialConnect_UT();

	vector<int> bflag(cp.size(), 0);
	for (uint i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].face.size() == 1)
		{
			bflag[tmedge[i].pt[0]] = 1;
			bflag[tmedge[i].pt[1]] = 1;
		}
	}

	for (int it = 0; it < nstep; it++)
	{
		for (uint i = 0; i < cp.size(); i++)
		{
			if (bflag[i] == 0 && (i == 183 || i == 21 || i == 178))
			{
				double tmp[3] = {0., 0., 0.};
				for (uint j = 0; j < cp[i].edge.size(); j++)
				{
					int pid(tmedge[cp[i].edge[j]].pt[0]);
					if (pid == i)
						pid = tmedge[cp[i].edge[j]].pt[1];
					tmp[0] += cp[pid].coor[0];
					tmp[1] += cp[pid].coor[1];
					tmp[2] += cp[pid].coor[2];
				}
				cp[i].coor[0] = tmp[0] / cp[i].edge.size();
				cp[i].coor[1] = tmp[1] / cp[i].edge.size();
				cp[i].coor[2] = tmp[2] / cp[i].edge.size();
			}
		}
	}

	VisualizeControlMesh(fn_out, 1);
}

void TruncatedTspline::DirichletLS_Scalar(vector<BezierElement> &bzmesh, vector<int> &IDBC, vector<double> &gh)
{
	//assume bzmesh[i].IEN[j] can be a passive one
}

void TruncatedTspline::Extrude(string fn, const vector<BezierElement> &bzmesh, double t, int nel, int dir, int p)
{
	//create an open knot vector with nel non-zero intervals
	vector<double> kv, kvb;
	for (int i = 0; i < p; i++)
	{
		kv.push_back(0.);  // spline knot vector
		kvb.push_back(0.); // bezier decomposition knot vector
	}
	for (int i = 0; i < nel; i++)
	{
		double tmp(double(i) / double(nel));
		kv.push_back(tmp);
		kvb.push_back(tmp);
		if (i > 0)
		{
			for (int j = 1; j < p; j++)
				kvb.push_back(tmp);
		}
	}
	for (int i = 0; i < p + 1; i++)
	{
		kv.push_back(1.);
		kvb.push_back(1.);
	}
	vector<vector<double>> matw;
	TMatrix(kv, kvb, p, matw);
	//for (uint i = 0; i < matw.size(); i++)
	//{
	//	for (uint j = 0; j < matw[i].size(); j++)
	//	{
	//		cout << matw[i][j] << " ";
	//	}
	//	cout << "\n";
	//}
	//getchar();

	int npw(kv.size() - p - 1); // number of control points
	cout << "npw = " << npw << endl;
	vector<double> ww(npw, 0.);
	for (int i = 0; i < npw; i++)
	{
		for (int j = 0; j < p; j++)
		{
			ww[i] += kv[i + j + 1];
		}
		ww[i] /= double(p);
	}

	//mid-surface control mesh
	vector<array<double, 3>> ptmid;
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			array<double, 3> tmp = {cp[i].coor[0], cp[i].coor[1], cp[i].coor[2]};
			ptmid.push_back(tmp);
		}
	}
	for (uint i = 0; i < bzcp.size(); i++)
		ptmid.push_back(bzcp[i]);
	vector<array<double, 3>> nm(ptmid.size());
	int loc(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			array<double, 3> tmp = {0., 0., 0.};
			int nfc(0);
			for (uint j = 0; j < cp[i].face.size(); j++)
			{
				int fcid(cp[i].face[j]);
				int *it = find(tmesh[fcid].cnct, tmesh[fcid].cnct + 4, i);
				int pos0 = it - tmesh[fcid].cnct;
				int pos[2] = {tmesh[fcid].cnct[(pos0 + 1) % 4], tmesh[fcid].cnct[(pos0 + 3) % 4]};
				array<double, 3> v1 = {cp[pos[0]].coor[0] - cp[i].coor[0], cp[pos[0]].coor[1] - cp[i].coor[1], cp[pos[0]].coor[2] - cp[i].coor[2]};
				array<double, 3> v2 = {cp[pos[1]].coor[0] - cp[i].coor[0], cp[pos[1]].coor[1] - cp[i].coor[1], cp[pos[1]].coor[2] - cp[i].coor[2]};
				array<double, 3> nmtmp = {v1[1] * v2[2] - v2[1] * v1[2], -v1[0] * v2[2] + v2[0] * v1[2], v1[0] * v2[1] - v2[0] * v1[1]};
				double dst = sqrt(nmtmp[0] * nmtmp[0] + nmtmp[1] * nmtmp[1] + nmtmp[2] * nmtmp[2]);
				if (dst < 1.e-12)
				{
					cerr << "Zero magnitute of normal!\n";
				}
				else
				{
					nmtmp[0] /= dst;
					nmtmp[1] /= dst;
					nmtmp[2] /= dst;
					tmp[0] += nmtmp[0];
					tmp[1] += nmtmp[1];
					tmp[2] += nmtmp[2];
					nfc++;
				}
			}
			if (nfc != 0)
			{
				nm[loc][0] = tmp[0] / double(nfc);
				nm[loc][1] = tmp[1] / double(nfc);
				nm[loc][2] = tmp[2] / double(nfc);
				double dst = sqrt(nm[loc][0] * nm[loc][0] + nm[loc][1] * nm[loc][1] + nm[loc][2] * nm[loc][2]);
				nm[loc][0] /= dst;
				nm[loc][1] /= dst;
				nm[loc][2] /= dst;
			}
			else
			{
				cerr << "Can't compute normal from neighboring faces! Default normal assumed!\n";
				nm[loc][0] = 0.;
				nm[loc][1] = 0.;
				nm[loc][2] = 1.;
			}
			loc++;
		}
	}
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].c1 == 1)
		{
			int pos[3] = {tmesh[i].cnct[0], tmesh[i].cnct[1], tmesh[i].cnct[3]};
			array<double, 3> v1 = {cp[pos[1]].coor[0] - cp[pos[0]].coor[0], cp[pos[1]].coor[1] - cp[pos[0]].coor[1], cp[pos[1]].coor[2] - cp[pos[0]].coor[2]};
			array<double, 3> v2 = {cp[pos[2]].coor[0] - cp[pos[0]].coor[0], cp[pos[2]].coor[1] - cp[pos[0]].coor[1], cp[pos[2]].coor[2] - cp[pos[0]].coor[2]};
			array<double, 3> nmtmp = {v1[1] * v2[2] - v2[1] * v1[2], -v1[0] * v2[2] + v2[0] * v1[2], v1[0] * v2[1] - v2[0] * v1[1]};
			double dst = sqrt(nmtmp[0] * nmtmp[0] + nmtmp[1] * nmtmp[1] + nmtmp[2] * nmtmp[2]);
			if (dst < 1.e-12)
			{
				cerr << "Zero magnitute of normal!\n";
				nmtmp[0] = 0.;
				nmtmp[1] = 0.;
				nmtmp[2] = 1.;
			}
			else
			{
				nmtmp[0] /= dst;
				nmtmp[1] /= dst;
				nmtmp[2] /= dst;
			}
			for (int j = 0; j < 4; j++)
			{
				int iloc(loc + tmesh[i].IENc1[j]);
				nm[iloc][0] = nmtmp[0];
				nm[iloc][1] = nmtmp[1];
				nm[iloc][2] = nmtmp[2];
			}
		}
	}
	//upper and lower layer of control mesh
	vector<array<double, 3>> pts(npw * ptmid.size());
	double thalf = 0.5 * t;
	if (dir == 1) //extruding along +w direction
	{
		for (uint i = 0; i < ptmid.size(); i++)
		{
			int i1(i + (npw - 1) * ptmid.size());
			pts[i][0] = ptmid[i][0];
			pts[i][1] = ptmid[i][1];
			pts[i][2] = ptmid[i][2];
			pts[i1][0] = ptmid[i][0] + t * nm[i][0];
			pts[i1][1] = ptmid[i][1] + t * nm[i][1];
			pts[i1][2] = ptmid[i][2] + t * nm[i][2];
		}
	}
	else if (dir == -1) //extruding along -w direction
	{
		for (uint i = 0; i < ptmid.size(); i++)
		{
			int i1(i + (npw - 1) * ptmid.size());
			pts[i][0] = ptmid[i][0] - t * nm[i][0];
			pts[i][1] = ptmid[i][1] - t * nm[i][1];
			pts[i][2] = ptmid[i][2] - t * nm[i][2];
			pts[i1][0] = ptmid[i][0];
			pts[i1][1] = ptmid[i][1];
			pts[i1][2] = ptmid[i][2];
		}
	}
	else
	{
		for (uint i = 0; i < ptmid.size(); i++) //input as midsurface
		{
			int i1(i + (npw - 1) * ptmid.size());
			pts[i][0] = ptmid[i][0] - thalf * nm[i][0];
			pts[i][1] = ptmid[i][1] - thalf * nm[i][1];
			pts[i][2] = ptmid[i][2] - thalf * nm[i][2];
			pts[i1][0] = ptmid[i][0] + thalf * nm[i][0];
			pts[i1][1] = ptmid[i][1] + thalf * nm[i][1];
			pts[i1][2] = ptmid[i][2] + thalf * nm[i][2];
		}
	}
	//middle layers
	for (uint iw = 1; iw < ww.size() - 1; iw++)
	{
		for (uint i = 0; i < ptmid.size(); i++)
		{
			int i1(i + iw * ptmid.size());
			int i2(i + (npw - 1) * ptmid.size());
			pts[i1][0] = (1. - ww[iw]) * pts[i][0] + ww[iw] * pts[i2][0];
			pts[i1][1] = (1. - ww[iw]) * pts[i][1] + ww[iw] * pts[i2][1];
			pts[i1][2] = (1. - ww[iw]) * pts[i][2] + ww[iw] * pts[i2][2];
		}
	}

	CollectActives_DPatchAnalysis(); //set paid

	vector<vector<int>> IENa(bzmesh.size());
	vector<vector<int>> loca(bzmesh.size());
	for (uint eid = 0; eid < bzmesh.size(); eid++)
	{
		for (uint i = 0; i < bzmesh[eid].IEN.size(); i++)
		{
			if (bzmesh[eid].IEN[i] >= 0 && paid[bzmesh[eid].IEN[i]] != -1)
			{
				int nnz(0);
				for (uint j = 0; j < bzmesh[eid].cmat[i].size(); j++)
				{
					if (bzmesh[eid].cmat[i][j] != 0.)
						nnz++;
				}
				if (nnz != 0)
				{
					IENa[eid].push_back(paid[bzmesh[eid].IEN[i]]);
					loca[eid].push_back(i);
				}
			}
		}
	}

	//cout << nel << " " << bzmesh.size() << "\n"; getchar();

	cout << "Writing file...\n";

	string fname = fn + "_extrude_BEXT.txt";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "B E X T\ntype hex\n";
		fout << "nodeN " << pts.size() << "\n";
		fout << "elemN " << nel * bzmesh.size() << "\n";
		fout.precision(8);
		int width(10);
		int npew(p + 1);
		int nnztol((bzmesh[0].order + 1) * (bzmesh[0].order + 1));
		if (nnztol < (bzmesh[0].order + 1) * npew)
			nnztol = (bzmesh[0].order + 1) * npew;
		for (uint i = 0; i < pts.size(); i++)
		{
			fout << "gnode " << setw(width) << pts[i][0] << " " << setw(width) << pts[i][1] << " " << setw(width) << pts[i][2] << setw(width) << " 1\n";
		}
		for (int iw = 0; iw < nel; iw++)
		{
			cout << "\nlayer: " << iw << "\n";
			for (uint eid = 0; eid < bzmesh.size(); eid++)
			{
				if (eid != 0 && eid % 500 == 0)
					cout << eid << " ";
				int npe3d(npew * IENa[eid].size()), nb3d((bzmesh[eid].order + 1) * (bzmesh[eid].order + 1) * npew);
				fout << "belem " << npe3d << " " << bzmesh[eid].order << " " << bzmesh[eid].order << " " << p << "\n";
				for (int i = 0; i < npew; i++)
				{
					for (uint j = 0; j < IENa[eid].size(); j++)
					{
						fout << setw(width) << (iw + i) * ptmid.size() + IENa[eid][j] << " ";
					}
				}
				fout << "\n";
				vector<vector<double>> mat3d(npe3d, vector<double>(nb3d, 0.));
				int counts(0);
				vector<int> nnzv(npe3d);
				for (int i = 0; i < npew; i++)
				{
					for (uint j = 0; j < IENa[eid].size(); j++)
					{
						int countb(0), nnz(0);
						for (int ki = 0; ki < npew; ki++)
						{
							double wtmp(matw[iw * p + ki][iw + i]);
							for (int kj = 0; kj < 16; kj++)
							{
								mat3d[counts][countb] = bzmesh[eid].cmat[loca[eid][j]][kj] * wtmp;
								if (mat3d[counts][countb] != 0.)
									nnz++;
								countb++;
							}
						}
						nnzv[counts] = nnz;
						counts++;
					}
				}
				for (uint i = 0; i < mat3d.size(); i++)
				{
					if (nnzv[i] > nnztol) //dense
					{
						fout << setw(width) << "d ";
						for (uint j = 0; j < mat3d[i].size(); j++)
							fout << setw(width) << mat3d[i][j] << " ";
						fout << "\n";
					}
					else //sparse
					{
						fout << setw(width) << "s " << setw(width) << nnzv[i] << " ";
						for (uint j = 0; j < mat3d[i].size(); j++)
						{
							if (mat3d[i][j] != 0.)
							{
								fout << setw(width) << j << " " << setw(width) << mat3d[i][j] << " ";
							}
						}
						fout << "\n";
					}
				}
			}
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	cout << "Done writing!\n";
}

void TruncatedTspline::Extrude_GEM_layer(string fn, const vector<BezierElement> &bzmesh, double t, int nel, int dir, int p)
{
	//create an open knot vector with nel non-zero intervals
	vector<double> kv, kvb;
	for (int i = 0; i < p; i++)
	{
		kv.push_back(0.);  // spline knot vector
		kvb.push_back(0.); // bezier decomposition knot vector
	}
	for (int i = 0; i < nel; i++)
	{
		double tmp(double(i) / double(nel));
		kv.push_back(tmp);
		kvb.push_back(tmp);
		if (i > 0)
		{
			for (int j = 1; j < p; j++)
				kvb.push_back(tmp);
		}
	}
	for (int i = 0; i < p + 1; i++)
	{
		kv.push_back(1.);
		kvb.push_back(1.);
	}
	vector<vector<double>> matw;
	TMatrix(kv, kvb, p, matw);

	int npw(kv.size() - p - 1); // number of control points
	cout << "npw = " << npw << endl;
	vector<double> ww(npw, 0.);
	for (int i = 0; i < npw; i++)
	{
		for (int j = 0; j < p; j++)
		{
			ww[i] += kv[i + j + 1];
		}
		ww[i] /= double(p);
	}

	//mid-surface control mesh
	vector<array<double, 3>> ptmid;
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			array<double, 3> tmp = {cp[i].coor[0], cp[i].coor[1], cp[i].coor[2]};
			ptmid.push_back(tmp);
		}
	}
	for (uint i = 0; i < bzcp.size(); i++)
		ptmid.push_back(bzcp[i]);
	vector<array<double, 3>> nm(ptmid.size());
	int loc(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			array<double, 3> tmp = {0., 0., 0.};
			int nfc(0);
			for (uint j = 0; j < cp[i].face.size(); j++)
			{
				int fcid(cp[i].face[j]);
				int *it = find(tmesh[fcid].cnct, tmesh[fcid].cnct + 4, i);
				int pos0 = it - tmesh[fcid].cnct;
				int pos[2] = {tmesh[fcid].cnct[(pos0 + 1) % 4], tmesh[fcid].cnct[(pos0 + 3) % 4]};
				array<double, 3> v1 = {cp[pos[0]].coor[0] - cp[i].coor[0], cp[pos[0]].coor[1] - cp[i].coor[1], cp[pos[0]].coor[2] - cp[i].coor[2]};
				array<double, 3> v2 = {cp[pos[1]].coor[0] - cp[i].coor[0], cp[pos[1]].coor[1] - cp[i].coor[1], cp[pos[1]].coor[2] - cp[i].coor[2]};
				array<double, 3> nmtmp = {v1[1] * v2[2] - v2[1] * v1[2], -v1[0] * v2[2] + v2[0] * v1[2], v1[0] * v2[1] - v2[0] * v1[1]};
				double dst = sqrt(nmtmp[0] * nmtmp[0] + nmtmp[1] * nmtmp[1] + nmtmp[2] * nmtmp[2]);
				if (dst < 1.e-12)
				{
					cerr << "Zero magnitute of normal!\n";
				}
				else
				{
					nmtmp[0] /= dst;
					nmtmp[1] /= dst;
					nmtmp[2] /= dst;
					tmp[0] += nmtmp[0];
					tmp[1] += nmtmp[1];
					tmp[2] += nmtmp[2];
					nfc++;
				}
			}
			if (nfc != 0)
			{
				nm[loc][0] = tmp[0] / double(nfc);
				nm[loc][1] = tmp[1] / double(nfc);
				nm[loc][2] = tmp[2] / double(nfc);
				double dst = sqrt(nm[loc][0] * nm[loc][0] + nm[loc][1] * nm[loc][1] + nm[loc][2] * nm[loc][2]);
				nm[loc][0] /= dst;
				nm[loc][1] /= dst;
				nm[loc][2] /= dst;
			}
			else
			{
				cerr << "Can't compute normal from neighboring faces! Default normal assumed!\n";
				nm[loc][0] = 0.;
				nm[loc][1] = 0.;
				nm[loc][2] = 1.;
			}
			loc++;
		}
	}
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].c1 == 1)
		{
			int pos[3] = {tmesh[i].cnct[0], tmesh[i].cnct[1], tmesh[i].cnct[3]};
			array<double, 3> v1 = {cp[pos[1]].coor[0] - cp[pos[0]].coor[0], cp[pos[1]].coor[1] - cp[pos[0]].coor[1], cp[pos[1]].coor[2] - cp[pos[0]].coor[2]};
			array<double, 3> v2 = {cp[pos[2]].coor[0] - cp[pos[0]].coor[0], cp[pos[2]].coor[1] - cp[pos[0]].coor[1], cp[pos[2]].coor[2] - cp[pos[0]].coor[2]};
			array<double, 3> nmtmp = {v1[1] * v2[2] - v2[1] * v1[2], -v1[0] * v2[2] + v2[0] * v1[2], v1[0] * v2[1] - v2[0] * v1[1]};
			double dst = sqrt(nmtmp[0] * nmtmp[0] + nmtmp[1] * nmtmp[1] + nmtmp[2] * nmtmp[2]);
			if (dst < 1.e-12)
			{
				cerr << "Zero magnitute of normal!\n";
				nmtmp[0] = 0.;
				nmtmp[1] = 0.;
				nmtmp[2] = 1.;
			}
			else
			{
				nmtmp[0] /= dst;
				nmtmp[1] /= dst;
				nmtmp[2] /= dst;
			}
			for (int j = 0; j < 4; j++)
			{
				int iloc(loc + tmesh[i].IENc1[j]);
				nm[iloc][0] = nmtmp[0];
				nm[iloc][1] = nmtmp[1];
				nm[iloc][2] = nmtmp[2];
			}
		}
	}
	//upper and lower layer of control mesh
	vector<array<double, 3>> pts(npw * ptmid.size());
	double thalf = 0.5 * t;
	if (dir == 1) //extruding along +w direction
	{
		for (uint i = 0; i < ptmid.size(); i++)
		{
			int i1(i + (npw - 1) * ptmid.size());
			pts[i][0] = ptmid[i][0];
			pts[i][1] = ptmid[i][1];
			pts[i][2] = ptmid[i][2];
			pts[i1][0] = ptmid[i][0] + t * nm[i][0];
			pts[i1][1] = ptmid[i][1] + t * nm[i][1];
			pts[i1][2] = ptmid[i][2] + t * nm[i][2];
		}
	}
	else if (dir == -1) //extruding along -w direction
	{
		for (uint i = 0; i < ptmid.size(); i++)
		{
			int i1(i + (npw - 1) * ptmid.size());
			pts[i][0] = ptmid[i][0] - t * nm[i][0];
			pts[i][1] = ptmid[i][1] - t * nm[i][1];
			pts[i][2] = ptmid[i][2] - t * nm[i][2];
			pts[i1][0] = ptmid[i][0];
			pts[i1][1] = ptmid[i][1];
			pts[i1][2] = ptmid[i][2];
		}
	}
	else
	{
		for (uint i = 0; i < ptmid.size(); i++) //input as midsurface
		{
			int i1(i + (npw - 1) * ptmid.size());
			pts[i][0] = ptmid[i][0] - thalf * nm[i][0];
			pts[i][1] = ptmid[i][1] - thalf * nm[i][1];
			pts[i][2] = ptmid[i][2] - thalf * nm[i][2];
			pts[i1][0] = ptmid[i][0] + thalf * nm[i][0];
			pts[i1][1] = ptmid[i][1] + thalf * nm[i][1];
			pts[i1][2] = ptmid[i][2] + thalf * nm[i][2];
		}
	}
	//middle layers
	for (uint iw = 1; iw < ww.size() - 1; iw++)
	{
		for (uint i = 0; i < ptmid.size(); i++)
		{
			int i1(i + iw * ptmid.size());
			int i2(i + (npw - 1) * ptmid.size());
			pts[i1][0] = (1. - ww[iw]) * pts[i][0] + ww[iw] * pts[i2][0];
			pts[i1][1] = (1. - ww[iw]) * pts[i][1] + ww[iw] * pts[i2][1];
			pts[i1][2] = (1. - ww[iw]) * pts[i][2] + ww[iw] * pts[i2][2];
		}
	}

	CollectActives_DPatchAnalysis(); //set paid

	vector<vector<int>> IENa(bzmesh.size());
	vector<vector<int>> loca(bzmesh.size());
	for (uint eid = 0; eid < bzmesh.size(); eid++)
	{
		for (uint i = 0; i < bzmesh[eid].IEN.size(); i++)
		{
			if (bzmesh[eid].IEN[i] >= 0 && paid[bzmesh[eid].IEN[i]] != -1)
			{
				int nnz(0);
				for (uint j = 0; j < bzmesh[eid].cmat[i].size(); j++)
				{
					if (bzmesh[eid].cmat[i][j] != 0.)
						nnz++;
				}
				if (nnz != 0)
				{
					IENa[eid].push_back(paid[bzmesh[eid].IEN[i]]);
					loca[eid].push_back(i);
				}
			}
		}
	}

	//cout << nel << " " << bzmesh.size() << "\n"; getchar();

	cout << "Writing file...\n";

	string fname = fn + "_extrude_BEXT.txt";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "B E X T\ntype hex\n";
		fout << "nodeN " << pts.size() << "\n";
		fout << "elemN " << nel * bzmesh.size() << "\n";
		fout.precision(8);
		int width(10);
		int npew(p + 1);
		int nnztol((bzmesh[0].order + 1) * (bzmesh[0].order + 1));
		if (nnztol < (bzmesh[0].order + 1) * npew)
			nnztol = (bzmesh[0].order + 1) * npew;
		for (uint i = 0; i < pts.size(); i++)
		{
			fout << "gnode " << setw(width) << pts[i][0] << " " << setw(width) << pts[i][1] << " " << setw(width) << pts[i][2] << setw(width) << " 1\n";
		}
		for (int iw = 0; iw < nel; iw++)
		{
			cout << "\nlayer: " << iw << "\n";
			for (uint eid = 0; eid < bzmesh.size(); eid++)
			{
				if (eid != 0 && eid % 500 == 0)
					cout << eid << " ";
				int npe3d(npew * IENa[eid].size()), nb3d((bzmesh[eid].order + 1) * (bzmesh[eid].order + 1) * npew);
				fout << "belem " << npe3d << " " << bzmesh[eid].order << " " << bzmesh[eid].order << " " << p << "\n";
				for (int i = 0; i < npew; i++)
				{
					for (uint j = 0; j < IENa[eid].size(); j++)
					{
						fout << setw(width) << (iw + i) * ptmid.size() + IENa[eid][j] << " ";
					}
				}
				fout << "\n";
				vector<vector<double>> mat3d(npe3d, vector<double>(nb3d, 0.));
				int counts(0);
				vector<int> nnzv(npe3d);
				for (int i = 0; i < npew; i++)
				{
					for (uint j = 0; j < IENa[eid].size(); j++)
					{
						int countb(0), nnz(0);
						for (int ki = 0; ki < npew; ki++)
						{
							double wtmp(matw[iw * p + ki][iw + i]);
							for (int kj = 0; kj < 16; kj++)
							{
								mat3d[counts][countb] = bzmesh[eid].cmat[loca[eid][j]][kj] * wtmp;
								if (mat3d[counts][countb] != 0.)
									nnz++;
								countb++;
							}
						}
						nnzv[counts] = nnz;
						counts++;
					}
				}
				for (uint i = 0; i < mat3d.size(); i++)
				{
					if (nnzv[i] > nnztol) //dense
					{
						fout << setw(width) << "d ";
						for (uint j = 0; j < mat3d[i].size(); j++)
							fout << setw(width) << mat3d[i][j] << " ";
						fout << "\n";
					}
					else //sparse
					{
						fout << setw(width) << "s " << setw(width) << nnzv[i] << " ";
						for (uint j = 0; j < mat3d[i].size(); j++)
						{
							if (mat3d[i][j] != 0.)
							{
								fout << setw(width) << j << " " << setw(width) << mat3d[i][j] << " ";
							}
						}
						fout << "\n";
					}
				}
			}
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	vector<int> IDa(cp.size() + bzcp.size(), -1);
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
			IDa[i] = count++;
	}
	int ncpa(count);
	for (uint i = 0; i < bzcp.size(); i++)
	{
		IDa[cp.size() + i] = count++;
	}

	fname = fn + "_extrude_CM_volume.vtk";
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		fout.precision(8);
		for (int i = 0; i < pts.size(); i++)
		{
			fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
		}

		vector<int> quad, qdc1;
		for (uint i = 0; i < tmesh.size(); i++)
		{
			if (tmesh[i].act == 1)
			{
				int flag(0);
				for (int j = 0; j < 4; j++)
				{
					if (cp[tmesh[i].cnct[j]].act == 0)
					{
						flag = 1;
						break;
					}
				}
				if (flag == 0)
				{
					quad.push_back(i);
				}
				if (tmesh[i].c1 == 1)
					qdc1.push_back(i);
			}
		}
		int nquad(quad.size() + qdc1.size());
		int npew(kv.size() - p - 2);
		//int npew(p);
		fout << "\nCELLS " << npew * nquad << " " << 9 * npew * nquad << '\n';
		for (int i = 0; i < npew; i++)
		{
			int shift1 = i * (pts.size()) / (npew + 1);
			int shift2 = (i + 1) * (pts.size()) / (npew + 1);
			for (uint e = 0; e < quad.size(); e++)
			{
				int i1(quad[e]);
				fout << "8 " << IDa[tmesh[i1].cnct[0]] + shift1 << " " << IDa[tmesh[i1].cnct[1]] + shift1 << " "
					 << IDa[tmesh[i1].cnct[2]] + shift1 << " " << IDa[tmesh[i1].cnct[3]] + shift1 << " "
					 << IDa[tmesh[i1].cnct[0]] + shift2 << " " << IDa[tmesh[i1].cnct[1]] + shift2 << " "
					 << IDa[tmesh[i1].cnct[2]] + shift2 << " " << IDa[tmesh[i1].cnct[3]] + shift2 << "\n";
			}
			for (uint e = 0; e < qdc1.size(); e++)
			{
				int i1(qdc1[e]);
				fout << "8 " << IDa[cp.size() + tmesh[i1].IENc1[0]] + shift1 << " " << IDa[cp.size() + tmesh[i1].IENc1[1]] + shift1 << " "
					 << IDa[cp.size() + tmesh[i1].IENc1[2]] + shift1 << " " << IDa[cp.size() + tmesh[i1].IENc1[3]] + shift1 << " "
					 << IDa[cp.size() + tmesh[i1].IENc1[0]] + shift2 << " " << IDa[cp.size() + tmesh[i1].IENc1[1]] + shift2 << " "
					 << IDa[cp.size() + tmesh[i1].IENc1[2]] + shift2 << " " << IDa[cp.size() + tmesh[i1].IENc1[3]] + shift2 << "\n";
			}
		}

		fout << "\nCELL_TYPES " << npew * nquad << '\n';
		for (int i = 0; i < npew * nquad; i++)
		{
			fout << "12\n";
		}

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	fname = fn + "_extrude_CM_layer.vtk";
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		fout.precision(8);
		for (int i = 0; i < pts.size(); i++)
		{
			fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
		}
		vector<int> quad, qdc1;
		for (uint i = 0; i < tmesh.size(); i++)
		{
			if (tmesh[i].act == 1)
			{
				int flag(0);
				for (int j = 0; j < 4; j++)
				{
					if (cp[tmesh[i].cnct[j]].act == 0)
					{
						flag = 1;
						break;
					}
				}
				if (flag == 0)
				{
					quad.push_back(i);
				}
				if (tmesh[i].c1 == 1)
					qdc1.push_back(i);
			}
		}
		int nquad(quad.size() + qdc1.size());
		int npew(kv.size() - p - 1);
		//int npew(p + 1);
		fout << "\nCELLS " << npew * nquad << " " << 5 * npew * nquad << '\n';
		for (int i = 0; i < npew; i++)
		{
			int shift = i * pts.size() / npew;
			for (uint e = 0; e < quad.size(); e++)
			{
				int i1(quad[e]);
				fout << "4 " << IDa[tmesh[i1].cnct[0]] + shift << " " << IDa[tmesh[i1].cnct[1]] + shift << " "
					 << IDa[tmesh[i1].cnct[2]] + shift << " " << IDa[tmesh[i1].cnct[3]] + shift << "\n";
			}
			for (uint e = 0; e < qdc1.size(); e++)
			{
				int i1(qdc1[e]);
				fout << "4 " << IDa[cp.size() + tmesh[i1].IENc1[0]] + shift << " " << IDa[cp.size() + tmesh[i1].IENc1[1]] + shift << " "
					 << IDa[cp.size() + tmesh[i1].IENc1[2]] + shift << " " << IDa[cp.size() + tmesh[i1].IENc1[3]] + shift << "\n";
			}
		}

		fout << "\nCELL_TYPES " << npew * nquad << '\n';
		for (int i = 0; i < npew * nquad; i++)
		{
			fout << "9\n";
		}

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	// Need further debugging
	fname = fn + "_extrude_CM_NewInterface.vtk";
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << ptmid.size() << " float\n";
		fout.precision(8);
		for (int i = 0; i < ptmid.size(); i++)
		{
			int i1(i + (npw - 1) * ptmid.size());
			if (dir == 1)
				fout << pts[i1][0] << " " << pts[i1][1] << " " << pts[i1][2] << "\n";
			if (dir == -1)
				fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
			if (dir == 0)
				fout << ptmid[i][0] << " " << ptmid[i][1] << " " << ptmid[i][2] << "\n";
		}
		//int npew(kv.size() - p - 1);
		//int npew(p + 1);
		//fout << "\nCELLS " << tmesh.size() << " " << 5 * tmesh.size() << '\n';

		int shift = 0;
		vector<int> quad, qdc1;
		for (uint i = 0; i < tmesh.size(); i++)
		{
			if (tmesh[i].act == 1)
			{
				int flag(0);
				for (int j = 0; j < 4; j++)
				{
					if (cp[tmesh[i].cnct[j]].act == 0)
					{
						flag = 1;
						break;
					}
				}
				if (flag == 0)
				{
					quad.push_back(i);
				}
				if (tmesh[i].c1 == 1)
					qdc1.push_back(i);
			}
		}
		int nquad(quad.size() + qdc1.size());
		fout << "\nCELLS " << nquad << " " << 5 * nquad << '\n';
		for (uint i = 0; i < quad.size(); i++)
		{
			int i1(quad[i]);
			fout << "4 " << IDa[tmesh[i1].cnct[0]] << " " << IDa[tmesh[i1].cnct[1]] << " "
				 << IDa[tmesh[i1].cnct[2]] << " " << IDa[tmesh[i1].cnct[3]] << "\n";
		}
		for (uint i = 0; i < qdc1.size(); i++)
		{
			int i1(qdc1[i]);
			fout << "4 " << IDa[cp.size() + tmesh[i1].IENc1[0]] << " " << IDa[cp.size() + tmesh[i1].IENc1[1]] << " "
				 << IDa[cp.size() + tmesh[i1].IENc1[2]] << " " << IDa[cp.size() + tmesh[i1].IENc1[3]] << "\n";
		}

		fout << "\nCELL_TYPES " << tmesh.size() << '\n';
		for (int i = 0; i < tmesh.size(); i++)
		{
			fout << "9\n";
		}

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	//fname = fn + "_extrude_CM_volume.vtk";
	//fout.open(fname.c_str());
	//if (fout.is_open())
	//{
	//	fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	//	fout << "POINTS " << pts.size() << " float\n";
	//	fout.precision(8);
	//	for (int i = 0; i < pts.size(); i++)
	//	{
	//		fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
	//	}
	//	int npew(kv.size() - p - 2);
	//	//int npew(p);
	//	fout << "\nCELLS " << npew * tmesh.size() << " " << 9 * npew * tmesh.size() << '\n';
	//	for (int i = 0; i < npew; i++)
	//	{
	//		int shift1 = i * pts.size() / (npew+1);
	//		int shift2 = (i+1) * pts.size() / (npew + 1);
	//		for (int e = 0; e < tmesh.size(); e++)
	//		{
	//			fout << "8 " << tmesh[e].cnct[0] + shift1 << ' ' << tmesh[e].cnct[1] + shift1 << ' ' << tmesh[e].cnct[2] + shift1 << ' ' << tmesh[e].cnct[3] + shift1 <<' '
	//				<< tmesh[e].cnct[0] + shift2 << ' ' << tmesh[e].cnct[1] + shift2 << ' ' << tmesh[e].cnct[2] + shift2 << ' ' << tmesh[e].cnct[3] + shift2 << '\n';
	//		}
	//	}

	//	fout << "\nCELL_TYPES " << npew * tmesh.size() << '\n';
	//	for (int i = 0; i < npew * tmesh.size(); i++)
	//	{
	//		fout << "12\n";
	//	}

	//	fout.close();
	//}
	//else
	//{
	//	cout << "Cannot open " << fname << "!\n";
	//}

	//fname = fn + "_extrude_CM_layer.vtk";
	//fout.open(fname.c_str());
	//if (fout.is_open())
	//{
	//	fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	//	fout << "POINTS " << pts.size() << " float\n";
	//	fout.precision(8);
	//	for (int i = 0; i < pts.size(); i++)
	//	{
	//		fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
	//	}
	//	int npew(kv.size() - p - 1);
	//	//int npew(p + 1);
	//	fout << "\nCELLS " << npew * tmesh.size() << " " << 5 * npew * tmesh.size() << '\n';
	//	for (int i = 0; i < npew; i++)
	//	{
	//		int shift = i * pts.size() / npew;
	//		for (int e = 0; e < tmesh.size(); e++)
	//		{
	//			fout << "4 " << tmesh[e].cnct[0] + shift << ' ' << tmesh[e].cnct[1] + shift << ' ' << tmesh[e].cnct[2] + shift << ' ' << tmesh[e].cnct[3] + shift << '\n';
	//		}
	//	}

	//	fout << "\nCELL_TYPES " << npew * tmesh.size() << '\n';
	//	for (int i = 0; i < npew * tmesh.size(); i++)
	//	{
	//		fout << "9\n";
	//	}

	//	fout.close();
	//}
	//else
	//{
	//	cout << "Cannot open " << fname << "!\n";
	//}

	//fname = fn + "_extrude_CM_NewInterface.vtk";
	//fout.open(fname.c_str());
	//if (fout.is_open())
	//{
	//	fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	//	fout << "POINTS " << ptmid.size() << " float\n";
	//	fout.precision(8);
	//	for (int i = 0; i < ptmid.size(); i++)
	//	{
	//		int i1(i + (npw - 1)*ptmid.size());
	//		if(dir == 1)
	//			fout << pts[i1][0] << " " << pts[i1][1] << " " << pts[i1][2] << "\n";
	//		if(dir == -1)
	//			fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
	//	}
	//	int npew(kv.size() - p - 1);
	//	//int npew(p + 1);
	//	fout << "\nCELLS " << tmesh.size() << " " << 5 * tmesh.size() << '\n';

	//	int shift = 0;
	//	for (int e = 0; e < tmesh.size(); e++)
	//	{
	//		fout << "4 " << tmesh[e].cnct[0] + shift << ' ' << tmesh[e].cnct[1] + shift << ' ' << tmesh[e].cnct[2] + shift << ' ' << tmesh[e].cnct[3] + shift << '\n';
	//	}

	//	fout << "\nCELL_TYPES " << tmesh.size() << '\n';
	//	for (int i = 0; i < tmesh.size(); i++)
	//	{
	//		fout << "9\n";
	//	}

	//	fout.close();
	//}
	//else
	//{
	//	cout << "Cannot open " << fname << "!\n";
	//}

	fname = fn + "_extrude_BEXT_layer.txt";
	fout.open(fname.c_str());
	vector<int> pid(cp.size() + bzcp.size(), -1);
	count = 0;
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			pid[i] = count++;
		}
	}
	for (uint i = 0; i < bzcp.size(); i++)
	{
		pid[i + cp.size()] = count++;
	}
	if (fout.is_open())
	{
		fout << "B E X T\ntype plane\n";
		fout << "nodeN " << pts.size() << "\n";
		fout << "elemN " << bzmesh.size() * (p + 1) << "\n";
		int width(16);
		for (uint i = 0; i < pts.size(); i++)
		{
			fout << "gnode " << setw(width) << pts[i][0] << " " << setw(width) << pts[i][1] << " " << setw(width) << pts[i][2] << setw(width) << " 1\n";
		}
		int npew(p + 1);
		for (uint jj = 0; jj < npew; jj++)
		{
			int shift = jj * pts.size() / npew;
			for (uint i = 0; i < bzmesh.size(); i++)
			{
				if (i != 0 && i % 500 == 0)
				{
					cout << i << " ";
				}
				vector<int> loc;
				vector<int> nnzv;
				if (bzmesh[i].order == 3)
				{
					for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
					{
						if (bzmesh[i].IEN[j] >= 0 && pid[bzmesh[i].IEN[j]] != -1)
						{
							int nnz(0);
							for (uint k = 0; k < bzmesh[i].cmat[j].size(); k++)
							{
								if (bzmesh[i].cmat[j][k] != 0.)
									nnz++;
							}
							if (nnz != 0)
							{
								loc.push_back(j);
								nnzv.push_back(nnz);
							}
						}
					}
				}
				else if (bzmesh[i].order == 4)
				{
					for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
					{
						if (bzmesh[i].IEN[j] >= 0 && pid[bzmesh[i].IEN[j]] != -1)
						{
							int nnz(0);
							for (uint k = 0; k < bzmesh[i].cmat4[j].size(); k++)
							{
								if (bzmesh[i].cmat4[j][k] != 0.)
									nnz++;
							}
							if (nnz != 0)
							{
								loc.push_back(j);
								nnzv.push_back(nnz);
							}
						}
					}
				}

				fout << "belem " << loc.size() << " " << bzmesh[i].order << " " << bzmesh[i].order << "\n";
				width = 10;
				for (uint j = 0; j < loc.size(); j++)
				{
					fout << setw(width) << pid[bzmesh[i].IEN[loc[j]]] + shift << " ";
				}
				fout << "\n";
				if (bzmesh[i].order == 3)
				{
					for (uint j = 0; j < loc.size(); j++)
					{
						for (uint k = 0; k < bzmesh[i].cmat[loc[j]].size(); k++)
							fout << setw(width) << bzmesh[i].cmat[loc[j]][k] << " ";
						fout << "\n";
					}
				}
				else if (bzmesh[i].order == 4)
				{
					for (uint j = 0; j < loc.size(); j++)
					{
						for (uint k = 0; k < bzmesh[i].cmat4[loc[j]].size(); k++)
							fout << setw(width) << bzmesh[i].cmat4[loc[j]][k] << " ";
						fout << "\n";
					}
				}
			}

			//for (uint j = 0; j < loc.size(); j++)
			//{
			//	//if (nnzv[j] > 20)//dense
			//	//{
			//	//	fout << setw(width) << "d ";
			//	//	for (uint k = 0; k < bzmesh[i].cmat[loc[j]].size(); k++) fout << setw(width) << bzmesh[i].cmat[loc[j]][k] << " ";
			//	//	fout << "\n";
			//	//}
			//	//else//sparse
			//	//{
			//	//	fout << setw(width) << "s " << setw(width) << nnzv[j] << " ";
			//	//	for (uint k = 0; k < bzmesh[i].cmat[loc[j]].size(); k++)
			//	//	{
			//	//		if (bzmesh[i].cmat[loc[j]][k] != 0.) fout << setw(width) << k << " " << setw(width) << bzmesh[i].cmat[loc[j]][k] << " ";
			//	//	}
			//	//	fout << "\n";
			//	//}
			//}
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	cout << "Done writing!\n";
}

void TruncatedTspline::Extrude_GEM_Shell(string fn, const vector<BezierElement> &bzmesh, double t, int nel, int dir, int p)
{
	//create an open knot vector with nel non-zero intervals
	vector<double> kv, kvb;
	for (int i = 0; i < p; i++)
	{
		kv.push_back(0.);  // spline knot vector
		kvb.push_back(0.); // bezier decomposition knot vector
	}
	for (int i = 0; i < nel; i++)
	{
		double tmp(double(i) / double(nel));
		kv.push_back(tmp);
		kvb.push_back(tmp);
		if (i > 0)
		{
			for (int j = 1; j < p; j++)
				kvb.push_back(tmp);
		}
	}
	for (int i = 0; i < p + 1; i++)
	{
		kv.push_back(1.);
		kvb.push_back(1.);
	}
	vector<vector<double>> matw;
	TMatrix(kv, kvb, p, matw);

	int npw(kv.size() - p - 1); // number of control points
	cout << "npw = " << npw << endl;
	vector<double> ww(npw, 0.);
	for (int i = 0; i < npw; i++)
	{
		for (int j = 0; j < p; j++)
		{
			ww[i] += kv[i + j + 1];
		}
		ww[i] /= double(p);
	}

	//mid-surface control mesh
	vector<array<double, 3>> ptmid;
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			array<double, 3> tmp = {cp[i].coor[0], cp[i].coor[1], cp[i].coor[2]};
			ptmid.push_back(tmp);
		}
	}
	for (uint i = 0; i < bzcp.size(); i++)
		ptmid.push_back(bzcp[i]);
	vector<array<double, 3>> nm(ptmid.size());
	int loc(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			array<double, 3> tmp = {0., 0., 0.};
			int nfc(0);
			for (uint j = 0; j < cp[i].face.size(); j++)
			{
				int fcid(cp[i].face[j]);
				int *it = find(tmesh[fcid].cnct, tmesh[fcid].cnct + 4, i);
				int pos0 = it - tmesh[fcid].cnct;
				int pos[2] = {tmesh[fcid].cnct[(pos0 + 1) % 4], tmesh[fcid].cnct[(pos0 + 3) % 4]};
				array<double, 3> v1 = {cp[pos[0]].coor[0] - cp[i].coor[0], cp[pos[0]].coor[1] - cp[i].coor[1], cp[pos[0]].coor[2] - cp[i].coor[2]};
				array<double, 3> v2 = {cp[pos[1]].coor[0] - cp[i].coor[0], cp[pos[1]].coor[1] - cp[i].coor[1], cp[pos[1]].coor[2] - cp[i].coor[2]};
				array<double, 3> nmtmp = {v1[1] * v2[2] - v2[1] * v1[2], -v1[0] * v2[2] + v2[0] * v1[2], v1[0] * v2[1] - v2[0] * v1[1]};
				double dst = sqrt(nmtmp[0] * nmtmp[0] + nmtmp[1] * nmtmp[1] + nmtmp[2] * nmtmp[2]);
				if (dst < 1.e-12)
				{
					cerr << "Zero magnitute of normal!\n";
				}
				else
				{
					nmtmp[0] /= dst;
					nmtmp[1] /= dst;
					nmtmp[2] /= dst;
					tmp[0] += nmtmp[0];
					tmp[1] += nmtmp[1];
					tmp[2] += nmtmp[2];
					nfc++;
				}
			}
			if (nfc != 0)
			{
				nm[loc][0] = tmp[0] / double(nfc);
				nm[loc][1] = tmp[1] / double(nfc);
				nm[loc][2] = tmp[2] / double(nfc);
				double dst = sqrt(nm[loc][0] * nm[loc][0] + nm[loc][1] * nm[loc][1] + nm[loc][2] * nm[loc][2]);
				nm[loc][0] /= dst;
				nm[loc][1] /= dst;
				nm[loc][2] /= dst;
			}
			else
			{
				cerr << "Can't compute normal from neighboring faces! Default normal assumed!\n";
				nm[loc][0] = 0.;
				nm[loc][1] = 0.;
				nm[loc][2] = 1.;
			}
			loc++;
		}
	}
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].act == 1 && tmesh[i].c1 == 1)
		{
			int pos[3] = {tmesh[i].cnct[0], tmesh[i].cnct[1], tmesh[i].cnct[3]};
			array<double, 3> v1 = {cp[pos[1]].coor[0] - cp[pos[0]].coor[0], cp[pos[1]].coor[1] - cp[pos[0]].coor[1], cp[pos[1]].coor[2] - cp[pos[0]].coor[2]};
			array<double, 3> v2 = {cp[pos[2]].coor[0] - cp[pos[0]].coor[0], cp[pos[2]].coor[1] - cp[pos[0]].coor[1], cp[pos[2]].coor[2] - cp[pos[0]].coor[2]};
			array<double, 3> nmtmp = {v1[1] * v2[2] - v2[1] * v1[2], -v1[0] * v2[2] + v2[0] * v1[2], v1[0] * v2[1] - v2[0] * v1[1]};
			double dst = sqrt(nmtmp[0] * nmtmp[0] + nmtmp[1] * nmtmp[1] + nmtmp[2] * nmtmp[2]);
			if (dst < 1.e-12)
			{
				cerr << "Zero magnitute of normal!\n";
				nmtmp[0] = 0.;
				nmtmp[1] = 0.;
				nmtmp[2] = 1.;
			}
			else
			{
				nmtmp[0] /= dst;
				nmtmp[1] /= dst;
				nmtmp[2] /= dst;
			}
			for (int j = 0; j < 4; j++)
			{
				int iloc(loc + tmesh[i].IENc1[j]);
				nm[iloc][0] = nmtmp[0];
				nm[iloc][1] = nmtmp[1];
				nm[iloc][2] = nmtmp[2];
			}
		}
	}
	//upper and lower layer of control mesh
	vector<array<double, 3>> pts(npw * ptmid.size());
	double thalf = 0.5 * t;
	if (dir == 1) //extruding along +w direction
	{
		for (uint i = 0; i < ptmid.size(); i++)
		{
			int i1(i + (npw - 1) * ptmid.size());
			pts[i][0] = ptmid[i][0];
			pts[i][1] = ptmid[i][1];
			pts[i][2] = ptmid[i][2];
			pts[i1][0] = ptmid[i][0] + t * nm[i][0];
			pts[i1][1] = ptmid[i][1] + t * nm[i][1];
			pts[i1][2] = ptmid[i][2] + t * nm[i][2];
		}
	}
	else if (dir == -1) //extruding along -w direction
	{
		for (uint i = 0; i < ptmid.size(); i++)
		{
			int i1(i + (npw - 1) * ptmid.size());
			pts[i][0] = ptmid[i][0] - t * nm[i][0];
			pts[i][1] = ptmid[i][1] - t * nm[i][1];
			pts[i][2] = ptmid[i][2] - t * nm[i][2];
			pts[i1][0] = ptmid[i][0];
			pts[i1][1] = ptmid[i][1];
			pts[i1][2] = ptmid[i][2];
		}
	}
	else
	{
		for (uint i = 0; i < ptmid.size(); i++) //input as midsurface
		{
			int i1(i + (npw - 1) * ptmid.size());
			pts[i][0] = ptmid[i][0] - thalf * nm[i][0];
			pts[i][1] = ptmid[i][1] - thalf * nm[i][1];
			pts[i][2] = ptmid[i][2] - thalf * nm[i][2];
			pts[i1][0] = ptmid[i][0] + thalf * nm[i][0];
			pts[i1][1] = ptmid[i][1] + thalf * nm[i][1];
			pts[i1][2] = ptmid[i][2] + thalf * nm[i][2];
		}
	}
	//middle layers
	for (uint iw = 1; iw < ww.size() - 1; iw++)
	{
		for (uint i = 0; i < ptmid.size(); i++)
		{
			int i1(i + iw * ptmid.size());
			int i2(i + (npw - 1) * ptmid.size());
			pts[i1][0] = (1. - ww[iw]) * pts[i][0] + ww[iw] * pts[i2][0];
			pts[i1][1] = (1. - ww[iw]) * pts[i][1] + ww[iw] * pts[i2][1];
			pts[i1][2] = (1. - ww[iw]) * pts[i][2] + ww[iw] * pts[i2][2];
		}
	}

	CollectActives_DPatchAnalysis(); //set paid

	vector<vector<int>> IENa(bzmesh.size());
	vector<vector<int>> loca(bzmesh.size());
	for (uint eid = 0; eid < bzmesh.size(); eid++)
	{
		for (uint i = 0; i < bzmesh[eid].IEN.size(); i++)
		{
			if (bzmesh[eid].IEN[i] >= 0 && paid[bzmesh[eid].IEN[i]] != -1)
			{
				int nnz(0);
				for (uint j = 0; j < bzmesh[eid].cmat[i].size(); j++)
				{
					if (bzmesh[eid].cmat[i][j] != 0.)
						nnz++;
				}
				if (nnz != 0)
				{
					IENa[eid].push_back(paid[bzmesh[eid].IEN[i]]);
					loca[eid].push_back(i);
				}
			}
		}
	}

	//cout << nel << " " << bzmesh.size() << "\n"; getchar();

	cout << "Writing file...\n";
	string fname;
	ofstream fout;

	vector<int> IDa(cp.size() + bzcp.size(), -1);
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
			IDa[i] = count++;
	}
	int ncpa(count);
	for (uint i = 0; i < bzcp.size(); i++)
	{
		IDa[cp.size() + i] = count++;
	}

	/// Connectivity information
	fname = fn + "_BotAndTop_CM_connectivity.vtk";
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		fout.precision(8);
		for (int i = 0; i < pts.size(); i++)
		{
			fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
		}
		vector<int> quad, qdc1;
		for (uint i = 0; i < tmesh.size(); i++)
		{
			if (tmesh[i].act == 1)
			{
				int flag(0);
				for (int j = 0; j < 4; j++)
				{
					if (cp[tmesh[i].cnct[j]].act == 0)
					{
						flag = 1;
						break;
					}
				}
				if (flag == 0)
				{
					quad.push_back(i);
				}
				if (tmesh[i].c1 == 1)
					qdc1.push_back(i);
			}
		}
		int nquad(quad.size() + qdc1.size());
		int npew(kv.size() - p - 1);
		//int npew(p + 1);
		fout << "\nCELLS " << npew * nquad << " " << 5 * npew * nquad << '\n';
		for (int i = 0; i < npew; i++)
		{
			int shift = i * pts.size() / npew;
			for (uint e = 0; e < quad.size(); e++)
			{
				int i1(quad[e]);
				fout << "4 " << IDa[tmesh[i1].cnct[0]] + shift << " " << IDa[tmesh[i1].cnct[1]] + shift << " "
					 << IDa[tmesh[i1].cnct[2]] + shift << " " << IDa[tmesh[i1].cnct[3]] + shift << "\n";
			}
			for (uint e = 0; e < qdc1.size(); e++)
			{
				int i1(qdc1[e]);
				fout << "4 " << IDa[cp.size() + tmesh[i1].IENc1[0]] + shift << " " << IDa[cp.size() + tmesh[i1].IENc1[1]] + shift << " "
					 << IDa[cp.size() + tmesh[i1].IENc1[2]] + shift << " " << IDa[cp.size() + tmesh[i1].IENc1[3]] + shift << "\n";
			}
		}

		fout << "\nCELL_TYPES " << npew * nquad << '\n';
		for (int i = 0; i < npew * nquad; i++)
		{
			fout << "9\n";
		}

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	vector<int> pid(cp.size() + bzcp.size(), -1);
	count = 0;
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			pid[i] = count++;
		}
	}
	for (uint i = 0; i < bzcp.size(); i++)
	{
		pid[i + cp.size()] = count++;
	}

	/// Reorder bezier elements
	map<int, vector<int>> bzmesh_dict;
	int tmpIEN;
	for (uint i = 0; i < IENa.size(); i++)
	{
		tmpIEN = IENa[i].size();
		if (bzmesh_dict.find(tmpIEN) == bzmesh_dict.end())
		{
			vector<int> tmp(1, i);
			bzmesh_dict[tmpIEN] = tmp;
			//cout << "tmpIEN " << tmpIEN << "inserted\n";
		}
		else
		{
			bzmesh_dict[tmpIEN].push_back(i);
		}
	}
	map<int, vector<int>>::iterator it = bzmesh_dict.begin();
	while (it != bzmesh_dict.end())
	{
		cout << "Type: " << it->first << endl;
		it++;
	}

	string fname_inp_bot = fn + "_Bot.inp";
	string fname_inp_top = fn + "_Top.inp";
	string fname_tspline_bot = fn + "_Bot.tsplines";
	string fname_tspline_top = fn + "_Top.tsplines";
	ofstream fout_inp[2];
	ofstream fout_tspline[2];
	fout_inp[0].open(fname_inp_bot.c_str());
	fout_inp[1].open(fname_inp_top.c_str());
	fout_tspline[0].open(fname_tspline_bot.c_str());
	fout_tspline[1].open(fname_tspline_top.c_str());

#pragma region output_inp
	if (fout_inp[0].is_open() && fout_inp[1].is_open())
	{
		fout_inp[0] << "*Heading"
					<< "\n";
		fout_inp[0] << "*Part, name=Part-1"
					<< "\n";

		fout_inp[1] << "*Heading"
					<< "\n";
		fout_inp[1] << "*Part, name=Part-1"
					<< "\n";

		fout_inp[0] << "*Node"
					<< "\n";
		fout_inp[1] << "*Node"
					<< "\n";

		int width(8);
		for (uint i = 0; i < pts.size() / 2; i++)
		{
			fout_inp[0] << i + 1 << "," << pts[i][0] << "," << pts[i][1] << "," << pts[i][2] << "\n";
			fout_inp[1] << i + 1 + pts.size() / 2 << "," << pts[i + pts.size() / 2][0] << "," << pts[i + pts.size() / 2][1] << "," << pts[i + pts.size() / 2][2] << "\n";
		}
		int npew(p + 1);
		for (uint jj = 0; jj < npew; jj++)
		{
			int shift = jj * pts.size() / npew;
			//int shift_ele = jj * bzmesh.size() / npew;
			//cout << shift_ele;
			map<int, vector<int>>::iterator it = bzmesh_dict.begin();
			while (it != bzmesh_dict.end())
			{
				fout_inp[jj] << "*USER ELEMENT, NODES=" << it->first << ", TYPE=U" << it->first << ", PROP=3, COORD=3,VARIABLES=216\n";
				fout_inp[jj] << "1\n";
				fout_inp[jj] << "*UEL PROPERTY, ELSET=uel" << it->first << "\n";
				fout_inp[jj] << "1,2,3\n";
				fout_inp[jj] << "*Element, type=U" << it->first << ", ELSET=uel" << it->first << "\n";
				vector<int> list_ele = it->second;
				for (int idx = 0; idx < list_ele.size(); idx++)
				{
					int count = 1;
					fout_inp[jj] << list_ele[idx] + 1 + bzmesh.size() * jj << ",";
					for (int pidx = 0; pidx < IENa[list_ele[idx]].size() - 1; pidx++)
					{
						fout_inp[jj] << IENa[list_ele[idx]][pidx] + 1 + shift << ",";
						count++;
						if (count % 16 == 0)
							fout_inp[jj] << "\n";
					}
					fout_inp[jj] << IENa[list_ele[idx]][IENa[list_ele[idx]].size() - 1] + 1 + shift << "\n";
				}

				it++;
			}
		}

		fout_inp[0] << "*End Part\n";
		fout_inp[0] << "**\n";
		fout_inp[0] << "**\n";
		fout_inp[0] << "** ASSEMBLY\n";
		fout_inp[0] << "**\n";
		fout_inp[0] << "Assembly, name=Assembly\n";
		fout_inp[0] << "**\n";
		fout_inp[0] << "*Instance, name=Part-1-1, part=Part-1\n";
		fout_inp[0] << "*End Instance\n";
		fout_inp[0] << "**\n";

		fout_inp[1] << "*End Part\n";
		fout_inp[1] << "**\n";
		fout_inp[1] << "**\n";
		fout_inp[1] << "** ASSEMBLY\n";
		fout_inp[1] << "**\n";
		fout_inp[1] << "Assembly, name=Assembly\n";
		fout_inp[1] << "**\n";
		fout_inp[1] << "*Instance, name=Part-1-1, part=Part-1\n";
		fout_inp[1] << "*End Instance\n";
		fout_inp[1] << "**\n";

		fout_inp[0].close();
		fout_inp[1].close();
	}
	else
	{
		cout << "Cannot open " << fname_inp_bot << "!\n";
	}
#pragma endregion

#pragma region output_tspline
	if (fout_tspline[0].is_open() && fout_tspline[1].is_open())
	{
		fout_tspline[0] << "**this is the T-spline input file for bot surface\n";
		fout_tspline[0] << "*thickness\n";
		fout_tspline[0] << t << "\n";
		fout_tspline[0] << "*number of layers\n";
		fout_tspline[0] << "1\n";
		fout_tspline[0] << "*fiber angle\n";
		fout_tspline[0] << "0.0\n";
		fout_tspline[0] << "*knotW\n";
		fout_tspline[0] << "0,1\n";
		fout_tspline[0] << "*elRangeW\n";
		fout_tspline[0] << "0,1\n";

		fout_tspline[1] << "**this is the T-spline input file for top surface\n";
		fout_tspline[1] << "*thickness\n";
		fout_tspline[1] << t << "\n";
		fout_tspline[1] << "*number of layers\n";
		fout_tspline[1] << "1\n";
		fout_tspline[1] << "*fiber angle\n";
		fout_tspline[1] << "0.0\n";
		fout_tspline[1] << "*knotW\n";
		fout_tspline[1] << "0,1\n";
		fout_tspline[1] << "*elRangeW\n";
		fout_tspline[1] << "0,1\n";

		fout_tspline[0] << "*coordinates\n";
		fout_tspline[1] << "*coordinates\n";
		int width(8);
		for (uint i = 0; i < pts.size() / 2; i++)
		{
			fout_tspline[0] << pts[i][0] << "\t" << pts[i][1] << "\t" << pts[i][2] << "\n";
			fout_tspline[1] << pts[i + pts.size() / 2][0] << "\t" << pts[i + pts.size() / 2][1] << "\t" << pts[i + pts.size() / 2][2] << "\n";
		}

		fout_tspline[0] << "*weights\n";
		fout_tspline[1] << "*weights\n";
		for (uint i = 0; i < pts.size() / 2; i++)
		{
			fout_tspline[0] << "1\n";
			fout_tspline[1] << "1\n";
		}

		fout_tspline[0] << "*bexts\n";
		fout_tspline[1] << "*bexts\n";
		int npew(p + 1);
		for (uint jj = 0; jj < npew; jj++)
		{
			int shift = jj * pts.size() / npew;
			for (uint i = 0; i < bzmesh.size(); i++)
			{
				if (i != 0 && i % 500 == 0)
				{
					cout << i << " ";
				}
				vector<int> loc;
				vector<int> nnzv;
				if (bzmesh[i].order == 3)
				{
					for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
					{
						if (bzmesh[i].IEN[j] >= 0 && pid[bzmesh[i].IEN[j]] != -1)
						{
							int nnz(0);
							for (uint k = 0; k < bzmesh[i].cmat[j].size(); k++)
							{
								if (bzmesh[i].cmat[j][k] != 0.)
									nnz++;
							}
							if (nnz != 0)
							{
								loc.push_back(j);
								nnzv.push_back(nnz);
							}
						}
					}
				}
				else if (bzmesh[i].order == 4)
				{
					for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
					{
						if (bzmesh[i].IEN[j] >= 0 && pid[bzmesh[i].IEN[j]] != -1)
						{
							int nnz(0);
							for (uint k = 0; k < bzmesh[i].cmat4[j].size(); k++)
							{
								if (bzmesh[i].cmat4[j][k] != 0.)
									nnz++;
							}
							if (nnz != 0)
							{
								loc.push_back(j);
								nnzv.push_back(nnz);
							}
						}
					}
				}

				fout_tspline[jj] << "belem " << loc.size() << " " << bzmesh[i].order << " " << bzmesh[i].order << "\n";
				width = 8;
				for (uint j = 0; j < loc.size(); j++)
				{
					fout_tspline[jj] << setw(width) << pid[bzmesh[i].IEN[loc[j]]] + shift << "\t";
				}
				fout_tspline[jj] << "\n";
				if (bzmesh[i].order == 3)
				{
					for (uint j = 0; j < loc.size(); j++)
					{
						for (uint k = 0; k < bzmesh[i].cmat[loc[j]].size(); k++)
							fout_tspline[jj] << setw(width) << bzmesh[i].cmat[loc[j]][k] << "\t";
						fout_tspline[jj] << "\n";
					}
				}
				else if (bzmesh[i].order == 4)
				{
					for (uint j = 0; j < loc.size(); j++)
					{
						for (uint k = 0; k < bzmesh[i].cmat4[loc[j]].size(); k++)
							fout_tspline[jj] << setw(width) << bzmesh[i].cmat4[loc[j]][k] << "\t";
						fout_tspline[jj] << "\n";
					}
				}
			}
		}
		fout_tspline[0].close();
		fout_tspline[1].close();
	}
	else
	{
		cout << "Cannot open tspline file!\n";
	}

#pragma endregion

	//fname = fn + "_extrude_BEXT_layer.txt";
	//fout.open(fname.c_str());
	////vector<int> pid(cp.size() + bzcp.size(), -1);
	////count = 0;
	////for (uint i = 0; i < cp.size(); i++)
	////{
	////	if (cp[i].act == 1)
	////	{
	////		pid[i] = count++;
	////	}
	////}
	////for (uint i = 0; i < bzcp.size(); i++)
	////{
	////	pid[i + cp.size()] = count++;
	////}
	//if (fout.is_open())
	//{
	//	fout << "B E X T\ntype plane\n";
	//	fout << "nodeN " << pts.size() << "\n";
	//	fout << "elemN " << bzmesh.size() * (p + 1) << "\n";
	//	int width(16);
	//	for (uint i = 0; i < pts.size(); i++)
	//	{
	//		fout << "gnode " << setw(width) << pts[i][0] << " " << setw(width) << pts[i][1] << " " << setw(width) << pts[i][2] << setw(width) << " 1\n";
	//	}
	//	int npew(p + 1);
	//	for (uint jj = 0; jj < npew; jj++)
	//	{
	//		int shift = jj * pts.size() / npew;
	//		for (uint i = 0; i < bzmesh.size(); i++)
	//		{
	//			if (i != 0 && i % 500 == 0)
	//			{
	//				cout << i << " ";
	//			}
	//			vector<int> loc;
	//			vector<int> nnzv;
	//			if (bzmesh[i].order == 3)
	//			{
	//				for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
	//				{
	//					if (bzmesh[i].IEN[j] >= 0 && pid[bzmesh[i].IEN[j]] != -1)
	//					{
	//						int nnz(0);
	//						for (uint k = 0; k < bzmesh[i].cmat[j].size(); k++)
	//						{
	//							if (bzmesh[i].cmat[j][k] != 0.) nnz++;
	//						}
	//						if (nnz != 0)
	//						{
	//							loc.push_back(j);
	//							nnzv.push_back(nnz);
	//						}
	//					}
	//				}
	//			}
	//			else if (bzmesh[i].order == 4)
	//			{
	//				for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
	//				{
	//					if (bzmesh[i].IEN[j] >= 0 && pid[bzmesh[i].IEN[j]] != -1)
	//					{
	//						int nnz(0);
	//						for (uint k = 0; k < bzmesh[i].cmat4[j].size(); k++)
	//						{
	//							if (bzmesh[i].cmat4[j][k] != 0.) nnz++;
	//						}
	//						if (nnz != 0)
	//						{
	//							loc.push_back(j);
	//							nnzv.push_back(nnz);
	//						}
	//					}
	//				}
	//			}

	//			fout << "belem " << loc.size() << " " << bzmesh[i].order << " " << bzmesh[i].order << "\n";
	//			width = 10;
	//			for (uint j = 0; j < loc.size(); j++)
	//			{
	//				fout << setw(width) << pid[bzmesh[i].IEN[loc[j]]] + shift << " ";
	//			}
	//			fout << "\n";
	//			if (bzmesh[i].order == 3)
	//			{
	//				for (uint j = 0; j < loc.size(); j++)
	//				{
	//					for (uint k = 0; k < bzmesh[i].cmat[loc[j]].size(); k++) fout << setw(width) << bzmesh[i].cmat[loc[j]][k] << " ";
	//					fout << "\n";
	//				}
	//			}
	//			else if (bzmesh[i].order == 4)
	//			{
	//				for (uint j = 0; j < loc.size(); j++)
	//				{
	//					for (uint k = 0; k < bzmesh[i].cmat4[loc[j]].size(); k++) fout << setw(width) << bzmesh[i].cmat4[loc[j]][k] << " ";
	//					fout << "\n";
	//				}
	//			}
	//		}
	//	}
	//	fout.close();
	//}
	//else
	//{
	//	cout << "Cannot open " << fname << "!\n";
	//}

	cout << "Done writing!\n";
}

void TruncatedTspline::WriteBezier_LSDYNA(string fn, vector<BezierElement> &bzmesh)
{
	cout << "Writing file...\n";
	string fname = fn + ".txt";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "B E X T\ntype plane\n";
		fout << "nodeN " << cp.size() << "\n";
		fout << "elemN " << bzmesh.size() << "\n";
		int width(16);
		for (uint i = 0; i < cp.size(); i++)
		{
			fout << "gnode " << setw(width) << cp[i].coor[0] << " " << setw(width) << cp[i].coor[1] << " " << setw(width) << cp[i].coor[2] << setw(width) << " 1\n";
		}
		for (uint i = 0; i < bzmesh.size(); i++)
		{
			if (i != 0 && i % 500 == 0)
			{
				cout << i << " ";
			}
			vector<int> loc;
			vector<int> nnzv;
			if (bzmesh[i].order == 3)
			{
				for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
				{
					if (bzmesh[i].IEN[j] >= 0)
					{
						int nnz(0);
						for (uint k = 0; k < bzmesh[i].cmat[j].size(); k++)
						{
							if (bzmesh[i].cmat[j][k] != 0.)
								nnz++;
						}
						if (nnz != 0)
						{
							loc.push_back(j);
							nnzv.push_back(nnz);
						}
					}
				}
			}
			else if (bzmesh[i].order == 4)
			{
				for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
				{
					if (bzmesh[i].IEN[j] >= 0)
					{
						int nnz(0);
						for (uint k = 0; k < bzmesh[i].cmat4[j].size(); k++)
						{
							if (bzmesh[i].cmat4[j][k] != 0.)
								nnz++;
						}
						if (nnz != 0)
						{
							loc.push_back(j);
							nnzv.push_back(nnz);
						}
					}
				}
			}

			fout << "belem " << loc.size() << " " << bzmesh[i].order << " " << bzmesh[i].order << "\n";
			width = 10;
			for (uint j = 0; j < loc.size(); j++)
			{
				fout << setw(width) << bzmesh[i].IEN[loc[j]] << " ";
			}
			fout << "\n";
			if (bzmesh[i].order == 3)
			{
				for (uint j = 0; j < loc.size(); j++)
				{
					for (uint k = 0; k < bzmesh[i].cmat[loc[j]].size(); k++)
						fout << setw(width) << bzmesh[i].cmat[loc[j]][k] << " ";
					fout << "\n";
				}
			}
			else if (bzmesh[i].order == 4)
			{
				for (uint j = 0; j < loc.size(); j++)
				{
					for (uint k = 0; k < bzmesh[i].cmat4[loc[j]].size(); k++)
						fout << setw(width) << bzmesh[i].cmat4[loc[j]][k] << " ";
					fout << "\n";
				}
			}
			//for (uint j = 0; j < loc.size(); j++)
			//{
			//	//if (nnzv[j] > 20)//dense
			//	//{
			//	//	fout << setw(width) << "d ";
			//	//	for (uint k = 0; k < bzmesh[i].cmat[loc[j]].size(); k++) fout << setw(width) << bzmesh[i].cmat[loc[j]][k] << " ";
			//	//	fout << "\n";
			//	//}
			//	//else//sparse
			//	//{
			//	//	fout << setw(width) << "s " << setw(width) << nnzv[j] << " ";
			//	//	for (uint k = 0; k < bzmesh[i].cmat[loc[j]].size(); k++)
			//	//	{
			//	//		if (bzmesh[i].cmat[loc[j]][k] != 0.) fout << setw(width) << k << " " << setw(width) << bzmesh[i].cmat[loc[j]][k] << " ";
			//	//	}
			//	//	fout << "\n";
			//	//}
			//}
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
	cout << "End of writing!\n";
}

void TruncatedTspline::WriteBezier_LSDYNA_C12(string fn, vector<BezierElement> &bzmesh)
{
	vector<int> pid(cp.size() + bzcp.size(), -1);
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			pid[i] = count++;
		}
	}
	for (uint i = 0; i < bzcp.size(); i++)
	{
		pid[i + cp.size()] = count++;
	}

	cout << "Writing file...\n";
	string fname = fn + "_BEXT.txt";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "B E X T\ntype plane\n";
		fout << "nodeN " << count << "\n";
		fout << "elemN " << bzmesh.size() << "\n";
		int width(16);
		for (uint i = 0; i < cp.size(); i++)
		{
			if (cp[i].act == 1)
			{
				fout << "gnode " << setw(width) << cp[i].coor[0] << " " << setw(width) << cp[i].coor[1] << " " << setw(width) << cp[i].coor[2] << setw(width) << " 1\n";
			}
		}
		for (uint i = 0; i < bzcp.size(); i++)
		{
			fout << "gnode " << setw(width) << bzcp[i][0] << " " << setw(width) << bzcp[i][1] << " " << setw(width) << bzcp[i][2] << setw(width) << " 1\n";
		}
		for (uint i = 0; i < bzmesh.size(); i++)
		{
			if (i != 0 && i % 500 == 0)
			{
				cout << i << " ";
			}
			vector<int> loc;
			vector<int> nnzv;
			if (bzmesh[i].order == 3)
			{
				for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
				{
					if (bzmesh[i].IEN[j] >= 0 && pid[bzmesh[i].IEN[j]] != -1)
					{
						int nnz(0);
						for (uint k = 0; k < bzmesh[i].cmat[j].size(); k++)
						{
							if (bzmesh[i].cmat[j][k] != 0.)
								nnz++;
						}
						if (nnz != 0)
						{
							loc.push_back(j);
							nnzv.push_back(nnz);
						}
					}
				}
			}
			else if (bzmesh[i].order == 4)
			{
				for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
				{
					if (bzmesh[i].IEN[j] >= 0 && pid[bzmesh[i].IEN[j]] != -1)
					{
						int nnz(0);
						for (uint k = 0; k < bzmesh[i].cmat4[j].size(); k++)
						{
							if (bzmesh[i].cmat4[j][k] != 0.)
								nnz++;
						}
						if (nnz != 0)
						{
							loc.push_back(j);
							nnzv.push_back(nnz);
						}
					}
				}
			}

			fout << "belem " << loc.size() << " " << bzmesh[i].order << " " << bzmesh[i].order << "\n";
			width = 10;
			for (uint j = 0; j < loc.size(); j++)
			{
				fout << setw(width) << pid[bzmesh[i].IEN[loc[j]]] << " ";
			}
			fout << "\n";
			if (bzmesh[i].order == 3)
			{
				for (uint j = 0; j < loc.size(); j++)
				{
					for (uint k = 0; k < bzmesh[i].cmat[loc[j]].size(); k++)
						fout << setw(width) << bzmesh[i].cmat[loc[j]][k] << " ";
					fout << "\n";
				}
			}
			else if (bzmesh[i].order == 4)
			{
				for (uint j = 0; j < loc.size(); j++)
				{
					for (uint k = 0; k < bzmesh[i].cmat4[loc[j]].size(); k++)
						fout << setw(width) << bzmesh[i].cmat4[loc[j]][k] << " ";
					fout << "\n";
				}
			}
			//for (uint j = 0; j < loc.size(); j++)
			//{
			//	//if (nnzv[j] > 20)//dense
			//	//{
			//	//	fout << setw(width) << "d ";
			//	//	for (uint k = 0; k < bzmesh[i].cmat[loc[j]].size(); k++) fout << setw(width) << bzmesh[i].cmat[loc[j]][k] << " ";
			//	//	fout << "\n";
			//	//}
			//	//else//sparse
			//	//{
			//	//	fout << setw(width) << "s " << setw(width) << nnzv[j] << " ";
			//	//	for (uint k = 0; k < bzmesh[i].cmat[loc[j]].size(); k++)
			//	//	{
			//	//		if (bzmesh[i].cmat[loc[j]][k] != 0.) fout << setw(width) << k << " " << setw(width) << bzmesh[i].cmat[loc[j]][k] << " ";
			//	//	}
			//	//	fout << "\n";
			//	//}
			//}
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
	cout << "End of writing!\n";
}

void TruncatedTspline::WriteBezier_Angran(vector<BezierElement> &bzmesh, string fn)
{
	int bzloc[4] = {0, 3, 15, 12};

	vector<int> paid(cp.size() + bzcp.size(), -1);
	int count(0);
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].act == 1)
		{
			paid[i] = count++;
		}
	}
	for (uint i = 0; i < bzcp.size(); i++)
	{
		paid[i + cp.size()] = count++;
	}

	string fname = fn + "/bzmesh.vtk";
	ofstream fout;
	//fout.open(fname.c_str());
	//if (fout.is_open())
	//{
	//	fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	//	fout << "POINTS " << 4 * bzmesh.size() << " float\n";
	//	for (uint i = 0; i < bzmesh.size(); i++)
	//	{
	//		for (int j = 0; j < 4; j++)
	//		{
	//			fout << bzmesh[i].pts[bzloc[j]][0] << " " << bzmesh[i].pts[bzloc[j]][1] << " " << bzmesh[i].pts[bzloc[j]][2] << "\n";
	//		}
	//	}
	//	fout << "\nCELLS " << bzmesh.size() << " " << 5 * bzmesh.size() << '\n';
	//	for (int i = 0; i < bzmesh.size(); i++)
	//	{
	//		fout << "4 " << 4 * i << " " << 4 * i + 1 << " " << 4 * i + 2 << " " << 4 * i + 3 << '\n';
	//	}
	//	fout << "\nCELL_TYPES " << bzmesh.size() << '\n';
	//	for (int i = 0; i < bzmesh.size(); i++)
	//	{
	//		fout << "9\n";
	//	}
	//}
	string fname3(fn + "/bzmeshinfo.txt");
	//ofstream fout;
	fout.open(fname3.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() << "\n";
		for (uint i = 0; i < bzmesh.size(); i++)
		{
			if (i != 0 && i % 500 == 0)
			{
				cout << i << " ";
			}
			vector<int> loc;
			vector<int> nnzv;
			{
				for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
				{
					if (bzmesh[i].IEN[j] >= 0 && paid[bzmesh[i].IEN[j]] != -1)
					{
						int nnz(0);
						for (uint k = 0; k < bzmesh[i].cmat[j].size(); k++)
						{
							if (bzmesh[i].cmat[j][k] != 0.)
								nnz++;
						}
						if (nnz != 0)
						{
							loc.push_back(j);
							nnzv.push_back(nnz);
						}
					}
				}
			}

			for (uint j = 0; j < loc.size(); j++)
			{
				fout << paid[bzmesh[i].IEN[loc[j]]] + 1;
				if (j == loc.size() - 1)
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

	string fname1(fn + "/cmat.txt");
	//ofstream fout;
	fout.open(fname1.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() << "\n";
		for (uint i = 0; i < bzmesh.size(); i++)
		{
			if (i != 0 && i % 500 == 0)
			{
				cout << i << " ";
			}
			vector<int> loc;
			vector<int> nnzv;
			{
				for (uint j = 0; j < bzmesh[i].IEN.size(); j++)
				{
					if (bzmesh[i].IEN[j] >= 0 && paid[bzmesh[i].IEN[j]] != -1)
					{
						int nnz(0);
						for (uint k = 0; k < bzmesh[i].cmat[j].size(); k++)
						{
							if (bzmesh[i].cmat[j][k] != 0.)
								nnz++;
						}
						if (nnz != 0)
						{
							loc.push_back(j);
							nnzv.push_back(nnz);
						}
					}
				}
			}

			fout << i << " " << loc.size() << " " << bzmesh[i].type << "\n";
			for (uint j = 0; j < loc.size(); j++)
			{
				fout << paid[bzmesh[i].IEN[loc[j]]];
				if (j == loc.size() - 1)
				{
					fout << "\n";
				}
				else
				{
					fout << " ";
				}
			}
			for (uint j = 0; j < loc.size(); j++)
			{
				for (uint k = 0; k < bzmesh[i].cmat[loc[j]].size(); k++)
				{
					fout << bzmesh[i].cmat[loc[j]][k] << " ";
					if (k == bzmesh[i].cmat[loc[j]].size() - 1)
					{
						fout << "\n";
					}
					else
					{
						fout << " ";
					}
				}
			}

			//for (uint j = 0; j < bzmesh[i].cmat.size(); j++)
			//{
			//	for (uint k = 0; k < bzmesh[i].cmat[j].size(); k++)
			//	{
			//		fout << bzmesh[i].cmat[j][k];
			//		if (k == bzmesh[i].cmat[j].size() - 1)
			//		{
			//			fout << "\n";
			//		}
			//		else
			//		{
			//			fout << " ";
			//		}
			//	}
			//}
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname1 << '\n';
	}

	string fname2(fn + "/bzpt.txt");
	//ofstream fout;
	fout.open(fname2.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() * 16 << "\n";
		for (uint i = 0; i < bzmesh.size(); i++)
		{
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

void TruncatedTspline::ReadBEXT(string fn, vector<BezierElement> &bzmesh)
{
	ifstream fin;
	fin.open(fn);
	if (fin.is_open())
	{
		string stmp, stmp1;
		int npt, nel, pcount(0), ecount(0);
		double dtmp;
		cp.clear();
		bzmesh.clear();
		while (getline(fin, stmp))
		{
			stringstream ss(stmp);
			if (stmp.compare(0, 5, "nodeN") == 0)
			{
				ss >> stmp1 >> npt;
				cp.resize(npt);
			}
			else if (stmp.compare(0, 5, "elemN") == 0)
			{
				ss >> stmp1 >> nel;
				bzmesh.resize(nel);
			}
			else if (stmp.compare(0, 5, "gnode") == 0)
			{
				ss >> stmp1 >> cp[pcount].coor[0] >> cp[pcount].coor[1] >> cp[pcount].coor[2] >> cp[pcount].w;
				pcount++;
			}
			else if (stmp.compare(0, 5, "belem") == 0)
			{
				int nIEN, ord1, ord2, nbn;
				ss >> stmp1 >> nIEN >> ord1 >> ord2;
				nbn = (ord1 + 1) * (ord2 + 1);
				bzmesh[ecount].IEN.resize(nIEN);
				bzmesh[ecount].cmat.resize(nIEN);
				for (int i = 0; i < nIEN; i++)
				{
					fin >> bzmesh[ecount].IEN[i];
				}
				for (int i = 0; i < nIEN; i++)
				{
					for (int j = 0; j < nbn; j++)
					{
						fin >> bzmesh[ecount].cmat[i][j];
						bzmesh[ecount].pts[j][0] += bzmesh[ecount].cmat[i][j] * cp[bzmesh[ecount].IEN[i]].coor[0];
						bzmesh[ecount].pts[j][1] += bzmesh[ecount].cmat[i][j] * cp[bzmesh[ecount].IEN[i]].coor[1];
						bzmesh[ecount].pts[j][2] += bzmesh[ecount].cmat[i][j] * cp[bzmesh[ecount].IEN[i]].coor[2];
					}
				}
				ecount++;
			}
		}
		fin.close();
	}
	else
	{
		cerr << "Can't open " << fn << "!\n";
	}
}

void TruncatedTspline::GeomMapBEXT(double u, double v, BezierElement &bzel, array<double, 3> &pt)
{
	double Bt[16], dBdt[16][2];
	bzel.Basis(u, v, Bt, dBdt);
	vector<double> Nt(bzel.IEN.size());
	vector<array<double, 2>> dNdt(bzel.IEN.size());
	pt[0] = 0.;
	pt[1] = 0.;
	pt[2] = 0.;
	for (uint i = 0; i < bzel.cmat.size(); i++)
	{
		Nt[i] = 0;
		dNdt[i][0] = 0.;
		dNdt[i][1] = 0.;
		for (uint j = 0; j < bzel.cmat[i].size(); j++)
		{
			Nt[i] += bzel.cmat[i][j] * Bt[j];
		}
		pt[0] += Nt[i] * cp[bzel.IEN[i]].coor[0];
		pt[1] += Nt[i] * cp[bzel.IEN[i]].coor[1];
		pt[2] += Nt[i] * cp[bzel.IEN[i]].coor[2];
	}
}

void TruncatedTspline::OutputBEXTGeom(string fn, vector<BezierElement> &bzmesh)
{
	vector<array<double, 3>> spt;
	vector<array<double, 3>> sval;
	vector<array<int, 4>> sele;
	vector<array<double, 3>> lpt; //visulize parameter lines
	vector<array<int, 2>> led;	  //line connectivity
	int ns(11), ecount(0);
	vector<double> su(ns), sv(ns);
	for (int i = 0; i < ns; i++)
	{
		su[i] = i * 1. / (ns - 1);
		sv[i] = i * 1. / (ns - 1);
	}

	for (uint e = 0; e < bzmesh.size(); e++)
	{
		for (int a = 0; a < ns; a++)
		{
			for (int b = 0; b < ns; b++)
			{
				array<double, 3> pt;
				GeomMapBEXT(su[b], sv[a], bzmesh[e], pt);
				spt.push_back(pt);
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
				el[0] = ecount * ns * ns + a * ns + b;
				el[1] = ecount * ns * ns + a * ns + b + 1;
				el[2] = ecount * ns * ns + (a + 1) * ns + b + 1;
				el[3] = ecount * ns * ns + (a + 1) * ns + b;
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

	string fname = fn + "_geom.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (uint i = 0; i < spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (uint i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	string fname1(fn + "_geom-lines.vtk");
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

void TruncatedTspline::ReadBEXT3D(string fn, vector<BezierElement> &bzmesh)
{
	//note that bzel.cp and bzel.amat is used instead of bzel.pts and bzel.cmat
	ifstream fin;
	fin.open(fn);
	if (fin.is_open())
	{
		string stmp, stmp1;
		int npt, nel, pcount(0), ecount(0);
		double dtmp;
		cp.clear();
		bzmesh.clear();
		while (getline(fin, stmp))
		{
			stringstream ss(stmp);
			if (stmp.compare(0, 5, "nodeN") == 0)
			{
				ss >> stmp1 >> npt;
				cp.resize(npt);
			}
			else if (stmp.compare(0, 5, "elemN") == 0)
			{
				ss >> stmp1 >> nel;
				bzmesh.resize(nel);
			}
			else if (stmp.compare(0, 5, "gnode") == 0)
			{
				ss >> stmp1 >> cp[pcount].coor[0] >> cp[pcount].coor[1] >> cp[pcount].coor[2] >> cp[pcount].w;
				pcount++;
			}
			else if (stmp.compare(0, 5, "belem") == 0)
			{
				int nIEN, ord1, ord2, ord3, nbn;
				ss >> stmp1 >> nIEN >> ord1 >> ord2 >> ord3;
				nbn = (ord1 + 1) * (ord2 + 1) * (ord3 + 1);
				//cout << nIEN << " " << nbn << "\n"; getchar();
				bzmesh[ecount].IEN.resize(nIEN);
				bzmesh[ecount].amat = MatrixXd::Zero(nIEN, nbn);
				bzmesh[ecount].cp.resize(nbn);
				for (int i = 0; i < nbn; i++)
				{
					bzmesh[ecount].cp[i][0] = 0.;
					bzmesh[ecount].cp[i][1] = 0.;
					bzmesh[ecount].cp[i][2] = 0.;
				}
				for (int i = 0; i < nIEN; i++)
				{
					fin >> bzmesh[ecount].IEN[i];
				}
				for (int i = 0; i < nIEN; i++)
				{
					string mflag;
					fin >> mflag;
					if (mflag.compare("s") == 0)
					{
						int nnz, loc;
						fin >> nnz;
						for (int j = 0; j < nnz; j++)
						{
							fin >> loc;
							fin >> bzmesh[ecount].amat(i, loc);
						}
					}
					else if (mflag.compare("d") == 0)
					{
						for (int j = 0; j < nbn; j++)
						{
							fin >> bzmesh[ecount].amat(i, j);
						}
					}
					else
					{
						cerr << "matrix flag read in error: " << mflag << "!\n";
						getchar();
					}
					for (int j = 0; j < nbn; j++)
					{
						if (bzmesh[ecount].amat(i, j) != 0.)
						{
							//cout << i <<" "<<j<<" "<< bzmesh[ecount].amat(i, j) << " "; getchar();
							bzmesh[ecount].cp[j][0] += bzmesh[ecount].amat(i, j) * cp[bzmesh[ecount].IEN[i]].coor[0];
							bzmesh[ecount].cp[j][1] += bzmesh[ecount].amat(i, j) * cp[bzmesh[ecount].IEN[i]].coor[1];
							bzmesh[ecount].cp[j][2] += bzmesh[ecount].amat(i, j) * cp[bzmesh[ecount].IEN[i]].coor[2];
						}
					}
				}
				ecount++;
			}
		}
		fin.close();
	}
	else
	{
		cerr << "Can't open " << fn << "!\n";
	}
}

void TruncatedTspline::OutputBezierMesh_BEXT3D(string fn, int p, vector<BezierElement> &bzmesh)
{
	int bzcn[8] = {0, 3, 15, 12, 48, 51, 63, 60}; //
	//int bzcn[8] = { 0,3,15,12,16,19,31,28 };

	if (p != 3)
		for (int i = 0; i < 4; i++)
			bzcn[i + 4] = bzcn[i] + 16 * p;

	string fname = fn + "_bzmesh3D.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nVolume\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << 8 * bzmesh.size() << " float\n";
		for (uint i = 0; i < bzmesh.size(); i++)
		{
			for (int j = 0; j < 8; j++)
			{
				fout << bzmesh[i].cp[bzcn[j]][0] << " " << bzmesh[i].cp[bzcn[j]][1] << " " << bzmesh[i].cp[bzcn[j]][2] << "\n";
			}
		}
		fout << "\nCELLS " << bzmesh.size() << " " << 9 * bzmesh.size() << '\n';
		for (uint i = 0; i < bzmesh.size(); i++)
		{
			fout << "8 ";
			for (int j = 0; j < 8; j++)
			{
				fout << 8 * i + j << " ";
			}
			fout << "\n";
		}
		fout << "\nCELL_TYPES " << bzmesh.size() << '\n';
		for (uint i = 0; i < bzmesh.size(); i++)
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

//void TruncatedTspline::Extrude(string fn, const vector<BezierElement>& bzmesh, double t, int nel, int dir, int p)
//{
//	//create an open knot vector with nel non-zero intervals
//	vector<double> kv, kvb;
//	for (int i = 0; i < p; i++)
//	{
//		kv.push_back(0.);
//		kvb.push_back(0.);
//	}
//	for (int i = 0; i < nel; i++)
//	{
//		double tmp(double(i) / double(nel));
//		kv.push_back(tmp);
//		kvb.push_back(tmp);
//		if (i > 0)
//		{
//			for (int j = 1; j < p; j++) kvb.push_back(tmp);
//		}
//	}
//	for (int i = 0; i < p + 1; i++)
//	{
//		kv.push_back(1.);
//		kvb.push_back(1.);
//	}
//	vector<vector<double>> matw;
//	TMatrix(kv, kvb, p, matw);
//	//for (uint i = 0; i < matw.size(); i++)
//	//{
//	//	for (uint j = 0; j < matw[i].size(); j++)
//	//	{
//	//		cout << matw[i][j] << " ";
//	//	}
//	//	cout << "\n";
//	//}
//	//getchar();
//
//	int npw(kv.size() - p - 1);
//	vector<double> ww(npw, 0.);
//	for (int i = 0; i < npw; i++)
//	{
//		for (int j = 0; j < p; j++)
//		{
//			ww[i] += kv[i + j + 1];
//		}
//		ww[i] /= double(p);
//	}
//
//	//mid-surface control mesh
//	vector<array<double, 3>> ptmid;
//	for (uint i = 0; i < cp.size(); i++)
//	{
//		if (cp[i].act == 1)
//		{
//			array<double, 3> tmp = { cp[i].coor[0],cp[i].coor[1], cp[i].coor[2] };
//			ptmid.push_back(tmp);
//		}
//	}
//	for (uint i = 0; i < bzcp.size(); i++) ptmid.push_back(bzcp[i]);
//	vector<array<double, 3>> nm(ptmid.size());
//	int loc(0);
//	for (uint i = 0; i < cp.size(); i++)
//	{
//		if (cp[i].act == 1)
//		{
//			array<double, 3> tmp = { 0.,0.,0. };
//			int nfc(0);
//			for (uint j = 0; j < cp[i].face.size(); j++)
//			{
//				int fcid(cp[i].face[j]);
//				int* it = find(tmesh[fcid].cnct, tmesh[fcid].cnct + 4, i);
//				int pos0 = it - tmesh[fcid].cnct;
//				int pos[2] = { tmesh[fcid].cnct[(pos0 + 1) % 4],tmesh[fcid].cnct[(pos0 + 3) % 4] };
//				array<double, 3> v1 = { cp[pos[0]].coor[0] - cp[i].coor[0],cp[pos[0]].coor[1] - cp[i].coor[1],cp[pos[0]].coor[2] - cp[i].coor[2] };
//				array<double, 3> v2 = { cp[pos[1]].coor[0] - cp[i].coor[0],cp[pos[1]].coor[1] - cp[i].coor[1],cp[pos[1]].coor[2] - cp[i].coor[2] };
//				array<double, 3> nmtmp = { v1[1] * v2[2] - v2[1] * v1[2],-v1[0] * v2[2] + v2[0] * v1[2],v1[0] * v2[1] - v2[0] * v1[1] };
//				double dst = sqrt(nmtmp[0] * nmtmp[0] + nmtmp[1] * nmtmp[1] + nmtmp[2] * nmtmp[2]);
//				if (dst < 1.e-12)
//				{
//					cerr << "Zero magnitute of normal!\n";
//				}
//				else
//				{
//					nmtmp[0] /= dst; nmtmp[1] /= dst; nmtmp[2] /= dst;
//					tmp[0] += nmtmp[0]; tmp[1] += nmtmp[1]; tmp[2] += nmtmp[2];
//					nfc++;
//				}
//			}
//			if (nfc != 0)
//			{
//				nm[loc][0] = tmp[0] / double(nfc); nm[loc][1] = tmp[1] / double(nfc); nm[loc][2] = tmp[2] / double(nfc);
//				double dst = sqrt(nm[loc][0] * nm[loc][0] + nm[loc][1] * nm[loc][1] + nm[loc][2] * nm[loc][2]);
//				nm[loc][0] /= dst; nm[loc][1] /= dst; nm[loc][2] /= dst;
//			}
//			else
//			{
//				cerr << "Can't compute normal from neighboring faces! Default normal assumed!\n";
//				nm[loc][0] = 0.; nm[loc][1] = 0.; nm[loc][2] = 1.;
//			}
//			loc++;
//		}
//	}
//	for (uint i = 0; i < tmesh.size(); i++)
//	{
//		if (tmesh[i].act == 1 && tmesh[i].c1 == 1)
//		{
//			int pos[3] = { tmesh[i].cnct[0],tmesh[i].cnct[1],tmesh[i].cnct[3] };
//			array<double, 3> v1 = { cp[pos[1]].coor[0] - cp[pos[0]].coor[0],cp[pos[1]].coor[1] - cp[pos[0]].coor[1],cp[pos[1]].coor[2] - cp[pos[0]].coor[2] };
//			array<double, 3> v2 = { cp[pos[2]].coor[0] - cp[pos[0]].coor[0],cp[pos[2]].coor[1] - cp[pos[0]].coor[1],cp[pos[2]].coor[2] - cp[pos[0]].coor[2] };
//			array<double, 3> nmtmp = { v1[1] * v2[2] - v2[1] * v1[2],-v1[0] * v2[2] + v2[0] * v1[2],v1[0] * v2[1] - v2[0] * v1[1] };
//			double dst = sqrt(nmtmp[0] * nmtmp[0] + nmtmp[1] * nmtmp[1] + nmtmp[2] * nmtmp[2]);
//			if (dst < 1.e-12)
//			{
//				cerr << "Zero magnitute of normal!\n";
//				nmtmp[0] = 0.; nmtmp[1] = 0.; nmtmp[2] = 1.;
//			}
//			else
//			{
//				nmtmp[0] /= dst; nmtmp[1] /= dst; nmtmp[2] /= dst;
//			}
//			for (int j = 0; j < 4; j++)
//			{
//				int iloc(loc + tmesh[i].IENc1[j]);
//				nm[iloc][0] = nmtmp[0]; nm[iloc][1] = nmtmp[1]; nm[iloc][2] = nmtmp[2];
//			}
//		}
//	}
//	//upper and lower layer of control mesh
//	vector<array<double, 3>> pts(npw*ptmid.size());
//	double thalf = 0.5*t;
//	if (dir == 1)//extruding along +w direction
//	{
//		for (uint i = 0; i < ptmid.size(); i++)
//		{
//			int i1(i + (npw - 1)*ptmid.size());
//			pts[i][0] = ptmid[i][0];
//			pts[i][1] = ptmid[i][1];
//			pts[i][2] = ptmid[i][2];
//			pts[i1][0] = ptmid[i][0] + t * nm[i][0];
//			pts[i1][1] = ptmid[i][1] + t * nm[i][1];
//			pts[i1][2] = ptmid[i][2] + t * nm[i][2];
//		}
//	}
//	else if (dir == -1)//extruding along -w direction
//	{
//		for (uint i = 0; i < ptmid.size(); i++)
//		{
//			int i1(i + (npw - 1)*ptmid.size());
//			pts[i][0] = ptmid[i][0] - t * nm[i][0];
//			pts[i][1] = ptmid[i][1] - t * nm[i][1];
//			pts[i][2] = ptmid[i][2] - t * nm[i][2];
//			pts[i1][0] = ptmid[i][0];
//			pts[i1][1] = ptmid[i][1];
//			pts[i1][2] = ptmid[i][2];
//		}
//	}
//	else
//	{
//		for (uint i = 0; i < ptmid.size(); i++)//input as midsurface
//		{
//			int i1(i + (npw - 1)*ptmid.size());
//			pts[i][0] = ptmid[i][0] - thalf * nm[i][0];
//			pts[i][1] = ptmid[i][1] - thalf * nm[i][1];
//			pts[i][2] = ptmid[i][2] - thalf * nm[i][2];
//			pts[i1][0] = ptmid[i][0] + thalf * nm[i][0];
//			pts[i1][1] = ptmid[i][1] + thalf * nm[i][1];
//			pts[i1][2] = ptmid[i][2] + thalf * nm[i][2];
//		}
//	}
//	//middle layers
//	for (uint iw = 1; iw < ww.size() - 1; iw++)
//	{
//		for (uint i = 0; i < ptmid.size(); i++)
//		{
//			int i1(i + iw * ptmid.size());
//			int i2(i + (npw - 1)*ptmid.size());
//			pts[i1][0] = (1. - ww[iw])*pts[i][0] + ww[iw] * pts[i2][0];
//			pts[i1][1] = (1. - ww[iw])*pts[i][1] + ww[iw] * pts[i2][1];
//			pts[i1][2] = (1. - ww[iw])*pts[i][2] + ww[iw] * pts[i2][2];
//		}
//	}
//
//	CollectActives_DPatchAnalysis();//set paid
//
//	vector<vector<int>> IENa(bzmesh.size());
//	vector<vector<int>> loca(bzmesh.size());
//	for (uint eid = 0; eid < bzmesh.size(); eid++)
//	{
//		for (uint i = 0; i < bzmesh[eid].IEN.size(); i++)
//		{
//			if (bzmesh[eid].IEN[i] >= 0 && paid[bzmesh[eid].IEN[i]] != -1)
//			{
//				int nnz(0);
//				for (uint j = 0; j < bzmesh[eid].cmat[i].size(); j++)
//				{
//					if (bzmesh[eid].cmat[i][j] != 0.) nnz++;
//				}
//				if (nnz != 0)
//				{
//					IENa[eid].push_back(paid[bzmesh[eid].IEN[i]]);
//					loca[eid].push_back(i);
//				}
//			}
//		}
//	}
//
//	//cout << nel << " " << bzmesh.size() << "\n"; getchar();
//
//	cout << "Writing file...\n";
//	string fname = fn + "_extrude_BEXT.txt";
//	ofstream fout;
//	fout.open(fname.c_str());
//	if (fout.is_open())
//	{
//		fout << "B E X T\ntype volume\n";
//		fout << "nodeN " << pts.size() << "\n";
//		fout << "elemN " << nel * bzmesh.size() << "\n";
//		fout.precision(8);
//		int width(10);
//		int npew(p + 1);
//		int nnztol((bzmesh[0].order + 1)*(bzmesh[0].order + 1));
//		if (nnztol < (bzmesh[0].order + 1)*npew) nnztol = (bzmesh[0].order + 1)*npew;
//		for (uint i = 0; i < pts.size(); i++)
//		{
//			fout << "gnode " << setw(width) << pts[i][0] << " " << setw(width) << pts[i][1] << " " << setw(width) << pts[i][2] << setw(width) << " 1\n";
//		}
//		for (int iw = 0; iw < nel; iw++)
//		{
//			cout << "\nlayer: " << iw << "\n";
//			for (uint eid = 0; eid < bzmesh.size(); eid++)
//			{
//				if (eid != 0 && eid % 500 == 0) cout << eid << " ";
//				int npe3d(npew*IENa[eid].size()), nb3d((bzmesh[eid].order + 1)*(bzmesh[eid].order + 1)*npew);
//				fout << "belem " << npe3d << " " << bzmesh[eid].order << " " << bzmesh[eid].order << " " << p << "\n";
//				for (int i = 0; i < npew; i++)
//				{
//					for (uint j = 0; j < IENa[eid].size(); j++)
//					{
//						fout << setw(width) << (iw + i)*ptmid.size() + IENa[eid][j] << " ";
//					}
//				}
//				fout << "\n";
//				vector<vector<double>> mat3d(npe3d, vector<double>(nb3d, 0.));
//				int counts(0);
//				vector<int> nnzv(npe3d);
//				for (int i = 0; i < npew; i++)
//				{
//					for (uint j = 0; j < IENa[eid].size(); j++)
//					{
//						int countb(0), nnz(0);
//						for (int ki = 0; ki < npew; ki++)
//						{
//							double wtmp(matw[iw*p + ki][iw + i]);
//							for (int kj = 0; kj < 16; kj++)
//							{
//								mat3d[counts][countb] = bzmesh[eid].cmat[loca[eid][j]][kj] * wtmp;
//								if (mat3d[counts][countb] != 0.) nnz++;
//								countb++;
//							}
//						}
//						nnzv[counts] = nnz;
//						counts++;
//					}
//				}
//				for (uint i = 0; i < mat3d.size(); i++)
//				{
//					if (nnzv[i] > nnztol)//dense
//					{
//						fout << setw(width) << "d ";
//						for (uint j = 0; j < mat3d[i].size(); j++) fout << setw(width) << mat3d[i][j] << " ";
//						fout << "\n";
//					}
//					else//sparse
//					{
//						fout << setw(width) << "s " << setw(width) << nnzv[i] << " ";
//						for (uint j = 0; j < mat3d[i].size(); j++)
//						{
//							if (mat3d[i][j] != 0.)
//							{
//								fout << setw(width) << j << " " << setw(width) << mat3d[i][j] << " ";
//							}
//						}
//						fout << "\n";
//					}
//				}
//			}
//		}
//		fout.close();
//	}
//	else
//	{
//		cout << "Cannot open " << fname << "!\n";
//	}
//	cout << "Done writing!\n";
//}

//void TruncatedTspline::ReadBEXT3D(string fn, vector<BezierElement>& bzmesh)
//{
//	//note that bzel.cp and bzel.amat is used instead of bzel.pts and bzel.cmat
//	ifstream fin;
//	fin.open(fn);
//	if (fin.is_open())
//	{
//		string stmp, stmp1;
//		int npt, nel, pcount(0), ecount(0);
//		double dtmp;
//		cp.clear();
//		bzmesh.clear();
//		while (getline(fin, stmp))
//		{
//			stringstream ss(stmp);
//			if (stmp.compare(0, 5, "nodeN") == 0)
//			{
//				ss >> stmp1 >> npt;
//				cp.resize(npt);
//			}
//			else if (stmp.compare(0, 5, "elemN") == 0)
//			{
//				ss >> stmp1 >> nel;
//				bzmesh.resize(nel);
//			}
//			else if (stmp.compare(0, 5, "gnode") == 0)
//			{
//				ss >> stmp1 >> cp[pcount].coor[0] >> cp[pcount].coor[1] >> cp[pcount].coor[2] >> cp[pcount].w;
//				pcount++;
//			}
//			else if (stmp.compare(0, 5, "belem") == 0)
//			{
//				int nIEN, ord1, ord2, ord3, nbn;
//				ss >> stmp1 >> nIEN >> ord1 >> ord2 >> ord3;
//				nbn = (ord1 + 1)*(ord2 + 1)*(ord3 + 1);
//				//cout << nIEN << " " << nbn << "\n"; getchar();
//				bzmesh[ecount].IEN.resize(nIEN);
//				bzmesh[ecount].amat = MatrixXd::Zero(nIEN, nbn);
//				bzmesh[ecount].cp.resize(nbn);
//				for (int i = 0; i < nbn; i++)
//				{
//					bzmesh[ecount].cp[i][0] = 0.;
//					bzmesh[ecount].cp[i][1] = 0.;
//					bzmesh[ecount].cp[i][2] = 0.;
//				}
//				for (int i = 0; i < nIEN; i++)
//				{
//					fin >> bzmesh[ecount].IEN[i];
//				}
//				for (int i = 0; i < nIEN; i++)
//				{
//					string mflag;
//					fin >> mflag;
//					if (mflag.compare("s") == 0)
//					{
//						int nnz, loc;
//						fin >> nnz;
//						for (int j = 0; j < nnz; j++)
//						{
//							fin >> loc;
//							fin >> bzmesh[ecount].amat(i, loc);
//						}
//					}
//					else if (mflag.compare("d") == 0)
//					{
//						for (int j = 0; j < nbn; j++)
//						{
//							fin >> bzmesh[ecount].amat(i, j);
//						}
//					}
//					else
//					{
//						cerr << "matrix flag read in error: " << mflag << "!\n";
//						getchar();
//					}
//					for (int j = 0; j < nbn; j++)
//					{
//						if (bzmesh[ecount].amat(i, j) != 0.)
//						{
//							//cout << i <<" "<<j<<" "<< bzmesh[ecount].amat(i, j) << " "; getchar();
//							bzmesh[ecount].cp[j][0] += bzmesh[ecount].amat(i, j) * cp[bzmesh[ecount].IEN[i]].coor[0];
//							bzmesh[ecount].cp[j][1] += bzmesh[ecount].amat(i, j) * cp[bzmesh[ecount].IEN[i]].coor[1];
//							bzmesh[ecount].cp[j][2] += bzmesh[ecount].amat(i, j) * cp[bzmesh[ecount].IEN[i]].coor[2];
//						}
//					}
//				}
//				ecount++;
//			}
//		}
//		fin.close();
//	}
//	else
//	{
//		cerr << "Can't open " << fn << "!\n";
//	}
//}
//
//void TruncatedTspline::OutputBezierMesh_BEXT3D(string fn, vector<BezierElement>& bzmesh)
//{
//	int bzcn[8] = { 0,3,15,12,48,51,63,60 };
//	string fname = fn + "_bzmesh3D.vtk";
//	ofstream fout;
//	fout.open(fname.c_str());
//	if (fout.is_open())
//	{
//		fout << "# vtk DataFile Version 2.0\nVolume\nASCII\nDATASET UNSTRUCTURED_GRID\n";
//		fout << "POINTS " << 8 * bzmesh.size() << " float\n";
//		for (uint i = 0; i < bzmesh.size(); i++)
//		{
//			for (int j = 0; j < 8; j++)
//			{
//				fout << bzmesh[i].cp[bzcn[j]][0] << " " << bzmesh[i].cp[bzcn[j]][1] << " " << bzmesh[i].cp[bzcn[j]][2] << "\n";
//			}
//		}
//		fout << "\nCELLS " << bzmesh.size() << " " << 9 * bzmesh.size() << '\n';
//		for (uint i = 0; i < bzmesh.size(); i++)
//		{
//			fout << "8 ";
//			for (int j = 0; j < 8; j++)
//			{
//				fout << 8 * i + j << " ";
//			}
//			fout << "\n";
//		}
//		fout << "\nCELL_TYPES " << bzmesh.size() << '\n';
//		for (uint i = 0; i < bzmesh.size(); i++)
//		{
//			fout << "12\n";
//		}
//		fout.close();
//	}
//	else
//	{
//		cout << "Cannot open " << fname << "!\n";
//	}
//}
