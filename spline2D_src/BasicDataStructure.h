#ifndef BASIC_DATA_STRUCTURE_H
#define BASIC_DATA_STRUCTURE_H

#include <vector>
#include <array>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//const int phys_dim(3);

/////////////////////////////////

class ELCS
{
public:
	int rot;//0, 1, 2, 3
	double u[2];
};

class Edge
{
public:
	Edge();
	int act;
	int pt[2];//two ends
	vector<int> face;//faces sharing this edge
	double len;//parametric length
	//int pned[2];//previous and next edge (2nd end), -1 for end, -2 for face, -3 for XP
	//int pnfc[2];//previous face and next face (T-jucntion case)
	int pn[2][2];//previous edge (or face) and next edge (or face), 1st flag (0 for edge, 1 for face, 2 for XP, 3 for end), 2nd ID
	int prt;
	int chd[2];
	int midpt;//middle point that splits the edge into two
	bool operator==(const Edge& ed);

	int add;
	//double pdual[3];
	int ploc;//correspond to which end
	int dual;
	int spoke;
};

class Vertex
{
public:
	double coor[3];
	double coortmp[3];//not useful
	double w;//NURBS weight
	double wtmp;
	double pm[2];//parametric coordinates, not useful
	int index[2];
	int uved[2];
	//PointIPP ipp;
	int update;//0 not updated, 1 updated, 2 partially updated (EPs), 3 later update tbf and tc
	double knotU[5];//temporally store local knot vector in its own hierarchy, only for boundary elements
	double knotV[5];//modify later
	double kitvU[4];//knot interval in local definition
	double kitvV[4];
	double min_itv[2];
	vector<int> uchd;
	vector<int> vchd;
	double kitvUtmp[4];
	double kitvVtmp[4];
	vector<int> vn;//store element number that connects to this vertex
	vector<int> sup;
	//vector<int> r1;
	//vector<int> r2;
	//int type;//extraordinary(2), boundary(1), regular(0), not necessary
	int act;//not useful
	int type;//0 for regular (including corner and boundary), 1 for T-junction, 2 for extraordinary

	int label; // ! Added by Angran, define boundary condition
	double velocity[3]; // ! Added by Angran, generate initial velocity field

	//double ctmp[3];//store a temporary coordinates of corner truncation, used to calculate new knot insertion
	int lev;
	int h2a;
	//int id;//0 regular, 1 truncated, 2 refined
	//int vln;//valence number
	//int ntr;//number of two ring elements
	//double err;
	int trun;
	int truntmp;
	int aff;//affected by knot insertion
	double kutmp[5];
	double kvtmp[5];
	//int cal;//0 for not calculated, 1 for calculated
	//double sc;//scaling coefficient for type 1
	vector<int> tbf;//truncated basis functions fot type 2
	vector<double> tc;//truncation coefficients for type 2
	vector<int> tbftmp;
	vector<double> tctmp;
	//vector<int> chd;
	//vector<double> chdc;

	vector<array<int,2>> child;//last dimention [lev,position]
	vector<double> c;//25 (5*5) children of cubic case

	vector<int> edge;//edges that connect to this vertex
	vector<int> face;//faces that share this vertex
	int rfc;//direction references

	Vertex();
	Vertex operator+(const Vertex& v);
	Vertex operator-(const Vertex& v);
	Vertex operator*(double a);
	Vertex operator/(double a);
	bool operator==(const Vertex& v);
};

class Element
{
public:
	int cnct[4];
	//vector<vector<int>> edge;
	int edge[4];
	vector<vector<int>> edge_act;
	int act;
	//int aff;
	int type;//0 for square, 1 for rectangular, 2 for zero area (boundary), 3 for corner boundary, 4 for extraordinary, 5 for invalid
	//int square;//1 for square, 0 for rectangular, not useful
	int lev;
	double lv;//element level, allow 0.5
	int h2a;
	//int nxp;//number of extraordinary vertices in the element
	int prt;
	int chd[4];
	double chd_o[4][2];
	int trun;//0 is not truncated and 1 is
	int ref;//0 for element without refinement, 1 for element to be refined, 2 for element with two T-jucntions (next), 3 for 2 T-junctions (opposite), 4 for 3 T-jucntions, 5 for 4 T-jucntions
	int nT;//# of T-junctions
	int Tjunc[4];//-1 for no T-junctions, other for T-junction ID, at most one T-jucntion on one edge, not useful
	//int nbed[4];//-1 for none, other for edge connecting T-junction
	int pn[4][2];//0 and 1 for egdge 3 and 1 (u-direction); 2 and 3 for edge 0 and 2 (v-direction); 1st flag and 2nd edge or face ID
	//int edge_nb[4][2];//first: -1 for none (boundary), other for element ID; second: locate which edge
	//vector<vector<int>> ednb;
	//int corn_nb[4];//(without XPs) -1 for none, other for element ID
	vector<int> newpID;//5 new points by subdivision
	//vector<array<int,2>> edge;//tmp vector, no use in fact
	vector<array<double,10>> kv;
	vector<int> IEN;//order in the desired way, in each hierarchy!
	vector<array<double,2>> pmcoor;
	vector<int> UV;//0: no rotate; 1: 90; 2: 180; 3: 270
	vector<array<int,2>> IENh;
	vector<int> nb;//one ring neighbor
	vector<int> nbrot;//rotation of corresponding edge neighbor
	vector<array<double,2>> nbog;
	//vector<Matrix3d> rotM;
	//vector<int> edge_neibors;//tmp vector to find vertx_neibors
	//vector<int> corner_neibors;//tmp vector to find vertx_neibors
	//vector<int> corner_nex;//store corner neighbors for extraordinary node in some order
	array<int,2> parent;
	vector<array<int,2>> child;
	//array<double,4> uhc;
	//double err;
	vector<int> node;
	vector<ELCS> lcs;//local coordinate system for each node
	vector<array<double,5>> patch_ku;
	vector<array<double,5>> patch_kv;
	vector<vector<double>> bemat;//for irregular element

	vector<int> IENtmp;
	vector<array<double,5>> patch_kutmp;
	vector<array<double,5>> patch_kvtmp;

	vector<int> poly;//when element type=6
	vector<int> polyed;

	int add;
	int dual;
	int focus;//0 not output error, and 1 output
	//double pdual[3];//assume only one EP per element
	int cnctd[4];//add connectivity for dual points, tmp use, later update to cnct
	double edlen[4];
	vector<MatrixXd> smat;
	MatrixXd amat;
	MatrixXd abar;
	vector<double> pku;
	vector<double> pkv;
	//vector<array<double, 8>> subku;
	//vector<array<double, 8>> subkv;

	int update;
	int c1;//0 for c2 element, 1 for c1 element, 2 for transition
	vector<vector<vector<double>>> bemat22;//corresponds to IEN, C2 functions, needs truncation
	vector<int> IENc1;//analysis space
	vector<vector<double>> c1mat;
	vector<vector<vector<double>>> c1mat22;
	vector<int> IENc2;//analysis space, only store active ones
	vector<array<double, 3>> bzpts;

	Element();
	bool operator==(const Element& e);
	//Element& operator=(const Element& e);
};

class BezierElement
{
public:
	BezierElement();
	int order;
	int prt;
	double pts[16][3];
//	vector<array<double,3>> pts;
	double pts4[25][3];
//	vector<array<double,3>> pts4;
	double w[16];
	double w4[25];
	vector<int> IEN;
	vector<double> w_nurbs;
	vector<array<double,16>> cmat;
	vector<array<double,25>> cmat4;
	vector<vector<int>> neum_ID;//store local of IEN
	vector<int> neum_edge;

	int type;
	int focus;
	MatrixXd amat;//used for subdivision basis functions
	MatrixXd abar;
	vector<array<double, 3>> cp;

	int bc[4];

	void Basis(double u, double v, double Nx[16], double dNdx[16][2]) const;
	void Basis4(double u, double v, double Nx[25], double dNdx[25][2]) const;
	void Para2Phys(double u, double v, double pt[3]);
	void Para2Phys4(double u, double v, double pt[3]);
	void SurfPointNormal(double u, double v, array<double,3>& pt, array<double,3>& nm) const;
	void SurfPointNormal4(double u, double v, array<double,3>& pt, array<double,3>& nm) const;
};

class PointIPP
{
public:
	int index[2];//index space
	double pm[2];//parameter space
	//double coor[3];//physical coordinates
	bool operator==(const PointIPP& pt);
};

//bool ComparePointIPPu(PointIPP ipp1, PointIPP ipp2)
//{
//	return (ipp1.pm[0] < ipp2.pm[0]);
//}
//
//bool ComparePointIPPv(PointIPP ipp1, PointIPP ipp2)
//{
//	return (ipp1.pm[1] < ipp2.pm[1]);
//}

class RegularPatchBasis
{
public:
	RegularPatchBasis(){for(int i=0;i<4;i++) {val[i]=0.; Dval[i]=0.; D2val[i]=0.;}};
	double val[4];
	double Dval[4];
	double D2val[4];
	void Evaluate(double x,double a=1.,double b=0.);
	void Clear();
};

#endif