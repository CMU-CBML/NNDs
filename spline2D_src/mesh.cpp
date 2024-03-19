#include "mesh.h"
#include <iomanip>
#include <sstream>
#include <cmath>

using namespace std;

typedef unsigned int uint;

void CreateHemiSphere(string fn)
{
	double r(1.);
	double corner[3][3] = { {r,0.,0.},{ 0.,r,0. },{ 0.,0.,r } };
	double edmid[3][3];
	double center[3] = { 0.,0.,0. };
	for (int i = 0; i < 3; i++)
	{
		edmid[i][0] = (corner[i][0] + corner[(i + 1) % 3][0]) / 2.;
		edmid[i][1] = (corner[i][1] + corner[(i + 1) % 3][1]) / 2.;
		edmid[i][2] = (corner[i][2] + corner[(i + 1) % 3][2]) / 2.;
		center[0] += corner[i][0] / 3.;
		center[1] += corner[i][1] / 3.;
		center[2] += corner[i][2] / 3.;
	}
	vector<array<double, 3>> pts0;
	for (int i = 0; i < 3; i++)
	{
		array<double, 3> tmp = { corner[i][0],corner[i][1], corner[i][2] };
		pts0.push_back(tmp);
	}
	for (int i = 0; i < 3; i++)
	{
		double len = sqrt(edmid[i][0] * edmid[i][0] + edmid[i][1] * edmid[i][1] + edmid[i][2] * edmid[i][2]);
		len = r / len;
		edmid[i][0] *= len; edmid[i][1] *= len; edmid[i][2] *= len;
		array<double, 3> tmp = { edmid[i][0],edmid[i][1], edmid[i][2] };
		pts0.push_back(tmp);
	}
	double lenc = sqrt(center[0] * center[0] + center[1] * center[1] + center[2] * center[2]);
	lenc = r / lenc;
	center[0] *= lenc; center[1] *= lenc; center[2] *= lenc;
	array<double, 3> tmpc = { center[0],center[1], center[2] };
	pts0.push_back(tmpc);

	int cnct0_tmp[3][4] = { {6,5,0,3},{6,3,1,4},{6,4,2,5} };
	vector<array<int, 4>> cnct0;
	for (int i = 0; i < 3; i++)
	{
		array<int, 4> tmp = { cnct0_tmp[i][0],cnct0_tmp[i][1],cnct0_tmp[i][2],cnct0_tmp[i][3] };
		cnct0.push_back(tmp);
	}

	int p(3);
	int nel(7);//8-1
	vector<double> kv;
	for (int i = 0; i < p; i++) kv.push_back(0.);
	for (int i = 0; i < nel; i++)
	{
		kv.push_back(double(i)/double(nel));
	}
	for (int i = 0; i < p + 1; i++) kv.push_back(1.);
	vector<double> w(kv.size() - p - 1, 0.);
	for (int i = 0; i < w.size(); i++)
	{
		for (int j = 0; j < p; j++)
		{
			w[i] += kv[i + j + 1];
		}
		w[i] /= double(p);
	}

	vector<array<double, 3>> pts;
	vector<array<int, 4>> cnct;
	for (int eid = 0; eid < cnct0.size(); eid++)
	{
		for (int i = 0; i < w.size(); i++)//v
		{
			for (int j = 0; j < w.size(); j++)//u
			{
				double wtmp[4] = { (1. - w[i]) * (1. - w[j]),(1. - w[i]) * w[j],w[i] * w[j],w[i] * (1. - w[j]) };
				array<double, 3> ptmp = { 0.,0.,0. };
				for (int k = 0; k < 4; k++)
				{
					ptmp[0] += wtmp[k] * pts0[cnct0[eid][k]][0];
					ptmp[1] += wtmp[k] * pts0[cnct0[eid][k]][1];
					ptmp[2] += wtmp[k] * pts0[cnct0[eid][k]][2];
				}
				double len = sqrt(ptmp[0] * ptmp[0] + ptmp[1] * ptmp[1] + ptmp[2] * ptmp[2]);
				len = r / len;
				ptmp[0] *= len; ptmp[1] *= len; ptmp[2] *= len;
				pts.push_back(ptmp);
			}
		}
		int ist = eid*w.size()*w.size();
		for (int i = 0; i < w.size() - 1; i++)//v
		{
			for (int j = 0; j < w.size() - 1; j++)//u
			{
				//array<int, 4> ctmp = { ist + i*w.size() + j,ist + i*w.size() + j + 1,ist + (i + 1)*w.size() + j + 1,ist + (i + 1)*w.size() + j };
				array<int, 4> ctmp;
				ctmp[0]=ist + i*w.size() + j;
				ctmp[1]=ist + i*w.size() + j + 1;
				ctmp[2]=ist + (i + 1)*w.size() + j + 1;
				ctmp[3]=ist + (i + 1)*w.size() + j;
				cnct.push_back(ctmp);
			}
		}
	}

	vector<int> pid(pts.size(),-1);
	vector<int> pflag(pts.size(), 0);
	int count(0);
	for (int i = 0; i < pts.size(); i++)
	{
		if (pflag[i] == 0)
		{
			pid[i] = count++;
			for (int j = i+1; j < pts.size(); j++)
			{
				double dif[3] = { pts[j][0] - pts[i][0],pts[j][1] - pts[i][1],pts[j][2] - pts[i][2] };
				double dist = sqrt(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]);
				if (dist < 1.e-6)
				{
					pid[j] = pid[i];
					pflag[j] = 1;
				}
			}
		}
	}

	string fname(fn + ".vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << count << " float\n";
		//fout << "POINTS " << pts.size() << " float\n";
		fout << setprecision(15);
		for (int i = 0; i < pts.size(); i++)
		{
			if (pflag[i] == 0)
			{
				fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
			}
		}
		fout << "\nCELLS " << cnct.size() << " " << 5 * cnct.size() << '\n';
		for (int i = 0; i<cnct.size(); i++)
		{
			fout << "4 " << pid[cnct[i][0]] << " " << pid[cnct[i][1]] << " " << pid[cnct[i][2]] << " " << pid[cnct[i][3]] << '\n';
			//fout << "4 " << cnct[i][0] << " " << cnct[i][1] << " " << cnct[i][2] << " " << cnct[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << cnct.size() << '\n';
		for (int i = 0; i < cnct.size(); i++)
		{
			fout << "9\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname << "!\n";
	}
}

void CreateCylinder(double radius, double len, int nel_u, int nel_v, string fn)
{
	const double PI(3.14159265359);
	const double PI_half(PI / 2.);
	double corner[4][3] = { { radius,0.,0. },{ 0.,radius,0. },{ radius,0.,len },{0.,radius,len} };
	//get weights
	int p(3);
	vector<double> ku, kv;
	for (int i = 0; i < p; i++)
	{
		ku.push_back(0.);
		kv.push_back(0.);
	}
	for (int i = 0; i < nel_u; i++) ku.push_back(double(i) / double(nel_u));
	for (int i = 0; i < nel_v; i++) kv.push_back(double(i) / double(nel_v));
	for (int i = 0; i < p + 1; i++)
	{
		ku.push_back(1.);
		kv.push_back(1.);
	}
	vector<double> wu(ku.size() - p - 1, 0.), wv(kv.size() - p - 1, 0.);
	for (int i = 0; i < wu.size(); i++)
	{
		for (int j = 0; j < p; j++)
		{
			wu[i] += ku[i + j + 1];
		}
		wu[i] /= double(p);
	}
	for (int i = 0; i < wv.size(); i++)
	{
		for (int j = 0; j < p; j++)
		{
			wv[i] += kv[i + j + 1];
		}
		wv[i] /= double(p);
	}

	int nptu(wu.size()), nptv(wv.size());
	vector<array<double, 3>> pts;
	vector<array<int, 4>> cnct;
	for (int i = 0; i < nptv; i++)
	{
		for (int j = 0; j < nptu; j++)
		{
			double wtmp[4] = { (1. - wv[i]) * (1. - wu[j]),(1. - wv[i]) * wu[j],wv[i] * (1. - wu[j]),wv[i] * wu[j] };
			array<double, 3> ptmp = { 0.,0.,0. };
			ptmp[0] = radius*cos(wu[j] * PI_half);
			ptmp[1] = radius*sin(wu[j] * PI_half);
			for (int k = 0; k < 4; k++)
			{
				//ptmp[0] += wtmp[k] * corner[k][0];
				//ptmp[1] += wtmp[k] * corner[k][1];
				ptmp[2] += wtmp[k] * corner[k][2];
			}
			//double dst = sqrt(ptmp[0] * ptmp[0] + ptmp[1] * ptmp[1]);
			//dst = radius / dst;
			//ptmp[0] *= dst; ptmp[1] *= dst;
			pts.push_back(ptmp);
		}
	}
	for (int i = 0; i < nptv - 1; i++)
	{
		for (int j = 0; j < nptu - 1; j++)
		{
			array<int, 4> ctmp = { i*nptu + j,i*nptu + j + 1,(i + 1)*nptu + j + 1,(i + 1)*nptu + j };
			cnct.push_back(ctmp);
		}
	}

	string fname(fn + ".vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		fout << setprecision(15);
		for (int i = 0; i < pts.size(); i++)
		{
			fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
		}
		fout << "\nCELLS " << cnct.size() << " " << 5 * cnct.size() << '\n';
		for (int i = 0; i<cnct.size(); i++)
		{
			fout << "4 " << cnct[i][0] << " " << cnct[i][1] << " " << cnct[i][2] << " " << cnct[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << cnct.size() << '\n';
		for (int i = 0; i < cnct.size(); i++)
		{
			fout << "9\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname << "!\n";
	}
}

void CreateCylinder_Irr(string fn_in, string fn_out, double radius, double len)
{
	vector<array<double, 3>> pts0;
	vector<array<int, 4>> cnct;
	ReadVtk_Quad(fn_in, pts0, cnct);
	double corner[4][3] = { { radius,0.,0. },{ 0.,radius,0. },{ radius,0.,len },{ 0.,radius,len } };
	//int cnid[4] = { 115,124,133,142 };//correspond to pts
	double urg[2] = { 1.e8,-1.e8 }, vrg[2] = { 1.e8,-1.e8 };
	for (uint i = 0; i < pts0.size(); i++)
	{
		if (pts0[i][0] < urg[0]) urg[0] = pts0[i][0];
		if (pts0[i][0] > urg[1]) urg[1] = pts0[i][0];
		if (pts0[i][1] < vrg[0]) vrg[0] = pts0[i][1];
		if (pts0[i][1] > vrg[1]) vrg[1] = pts0[i][1];
	}
	double udim(urg[1] - urg[0]), vdim(vrg[1] - vrg[0]);

	vector<array<double, 3>> pts;
	for (uint i = 0; i < pts0.size(); i++)
	{
		double wu = (pts0[i][0] - urg[0]) / udim;
		double wv = (pts0[i][1] - vrg[0]) / vdim;
		double wtmp[4] = { (1. - wv) * (1. - wu),(1. - wv) * wu,wv * (1. - wu),wv * wu };
		array<double, 3> ptmp = { 0.,0.,0. };
		for (int k = 0; k < 4; k++)
		{
			ptmp[0] += wtmp[k] * corner[k][0];
			ptmp[1] += wtmp[k] * corner[k][1];
			ptmp[2] += wtmp[k] * corner[k][2];
		}
		double dst = sqrt(ptmp[0] * ptmp[0] + ptmp[1] * ptmp[1]);
		dst = radius / dst;
		ptmp[0] *= dst; ptmp[1] *= dst;
		pts.push_back(ptmp);
	}

	WriteVtk_Quad(fn_out, pts, cnct);
}

void ReadAbaqusInp_Quad(string fn, vector<array<double, 3>>& pts, vector<array<int, 4>>& cnct,
	vector<int>& indxp, vector<int>& indxe)
{
	ifstream fin;
	fin.open(fn);
	if (fin.is_open())
	{
		string stmp;
		char ch;
		int nflag(0), eflag(0);
		int ncount(0), ecount(0);
		while (getline(fin,stmp))
		{
			//read *NODE
			if (stmp.compare(0, 5, "*NODE")==0 || stmp.compare(0, 5, "*Node") == 0)
			{
				nflag = 1;
				eflag = 0;
				if(!getline(fin, stmp)) break;
			}
			else if (stmp.compare(0, 8, "*ELEMENT") == 0 || stmp.compare(0, 8, "*Element") == 0)
			{
				eflag = 1;
				nflag = 0;
				if (!getline(fin, stmp)) break;
			}
			else if (stmp.compare(0,1,"*") == 0)
			{
				nflag = 0;
				eflag = 0;
			}
			if (nflag)
			{
				stringstream ss(stmp);
				array<double, 3> ptmp;
				int itmp;
				ss >> itmp >> ch >> ptmp[0] >> ch >> ptmp[1] >> ch >> ptmp[2];
				pts.push_back(ptmp);
				indxp.push_back(itmp);
			}
			if (eflag)
			{
				stringstream ss(stmp);
				array<int, 4> etmp;
				int itmp;
				ss >> itmp >> ch >> etmp[0] >> ch >> etmp[1] >> ch >> etmp[2] >> ch >> etmp[3];
				cnct.push_back(etmp);
				indxe.push_back(itmp);
			}
		}
		fin.close();
	}
	else
	{
		cerr << "Can't open " << fn << "!\n";
	}
}

void WriteVtk_Quad(string fn, const vector<array<double, 3>>& pts, const vector<array<int, 4>>& cnct)
{
	unsigned int i;
	ofstream fout;
	fout.open(fn);
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		fout.precision(16);
		for (i = 0; i<pts.size(); i++)
		{
			fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
		}
		fout << "\nCELLS " << cnct.size() << " " << 5 * cnct.size() << '\n';
		for (i = 0; i<cnct.size(); i++)
		{
			fout << "4 " << cnct[i][0] << ' ' << cnct[i][1] << ' ' << cnct[i][2] << ' ' << cnct[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << cnct.size() << '\n';
		for (i = 0; i<cnct.size(); i++)
		{
			fout << "9\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fn << "!\n";
	}
}

void ReadVtk_Quad(string fn, vector<array<double, 3>>& pts, vector<array<int, 4>>& cnct)
{
	string stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fn);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		pts.resize(npts);
		for (int i = 0; i<npts; i++)
		{
			fin >> pts[i][0] >> pts[i][1] >> pts[i][2];
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		cnct.resize(neles);
		for (int i = 0; i<neles; i++)
		{
			fin >> itmp >> cnct[i][0] >> cnct[i][1] >> cnct[i][2] >> cnct[i][3];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fn << "!\n";
	}
}

void Inp2Vtk_Quad(string fn_in, string fn_out)
{
	vector<array<double, 3>> pts;
	vector<array<int, 4>> cnct;
	vector<int> indxp, indxe;
	ReadAbaqusInp_Quad(fn_in, pts, cnct, indxp, indxe);
	int pmax(0), pmin(1000000000);
	unsigned int i;
	for (i = 0; i < indxp.size(); i++)
	{
		if (pmax < indxp[i]) pmax = indxp[i];
		if (pmin > indxp[i]) pmin = indxp[i];
	}
	int nresv(pmax - pmin + 1);
	vector<int> aflag(nresv, 0), pid(nresv, -1);
	for (i = 0; i < indxe.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			aflag[cnct[i][j] - pmin] = 1;
		}
	}
	vector<array<double, 3>> pts1;
	for (i = 0; i < indxp.size(); i++)
	{
		if (aflag[indxp[i] - pmin] == 1)
		{
			pts1.push_back(pts[i]);
		}
	}
	int count(0);
	for (i = 0; i < pid.size(); i++)
	{
		if (aflag[i] == 1)
		{
			pid[i] = count++;
		}
	}
	for (i = 0; i < indxe.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cnct[i][j] = pid[cnct[i][j] - pmin];
		}
	}

	WriteVtk_Quad(fn_out, pts1, cnct);
}

void RemoveDuplicatePoints(string fn_in, string fn_out)
{
	vector<array<double, 3>> pts, pts1;
	vector<array<int, 4>> cnct;
	ReadVtk_Quad(fn_in, pts, cnct);
	double tol(1.e5);
	for (uint i = 0; i < cnct.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			int loc1(cnct[i][j]), loc2(cnct[i][(j + 1) % 4]);
			double dist = (pts[loc1][0] - pts[loc2][0])*(pts[loc1][0] - pts[loc2][0]) +
				(pts[loc1][1] - pts[loc2][1])*(pts[loc1][1] - pts[loc2][1]) +
				(pts[loc1][2] - pts[loc2][2])*(pts[loc1][2] - pts[loc2][2]);
			dist = sqrt(dist);
			if (dist < tol) tol = dist;
		}
	}
	tol *= 0.2;
	vector<int> rept(pts.size(), -1), pid(pts.size(), -1);
	for (uint i = 0; i < pts.size(); i++)
	{
		if (rept[i] == -1)
		{
			for (uint j = i + 1; j < pts.size(); j++)
			{
				double dist = (pts[i][0] - pts[j][0])*(pts[i][0] - pts[j][0]) +
					(pts[i][1] - pts[j][1])*(pts[i][1] - pts[j][1]) +
					(pts[i][2] - pts[j][2])*(pts[i][2] - pts[j][2]);
				dist = sqrt(dist);
				if (dist < tol)
				{
					rept[j] = i;
				}
			}
		}
	}
	int count(0);
	for (uint i = 0; i < pts.size(); i++)
	{
		if (rept[i] == -1)
		{
			pid[i] = count++;
			pts1.push_back(pts[i]);
		}
		else
		{
			pid[i] = pid[rept[i]];
		}
	}
	for (uint i = 0; i < cnct.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cnct[i][j] = pid[cnct[i][j]];
		}
	}

	WriteVtk_Quad(fn_out, pts1, cnct);
}