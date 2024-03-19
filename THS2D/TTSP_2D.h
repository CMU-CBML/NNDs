#ifndef TTSP_2D_H
#define TTSP_2D_H

#include <vector>
#include <utility>
#include <string>
//#include "T_mesh.h"
#include "BasicDataStructure.h"
using namespace std;

class TruncatedTspline_2D
{
public:
	TruncatedTspline_2D();
	void VisualizeControlMesh(string fn);
	double PartitionOfUnity(int eid,double u,double v);
	void CollectActives();
	void StrongBalanceCheck(const vector<int>& rid, vector<int>& ridsb);
	void FaceIntersectCheck(vector<int>& ridtp);//one T-junction requirement
	void TjuncExtentCheck();
	void TjuncExtentCheck_1(vector<int>& ridtjx);
	void TjuncExtentCheck_2(vector<int>& ridtjx);
	void StrongBalanceRefine(const vector<int>& rid);
	void TargetRefine(const vector<int>& rid);
	void OneTjunctionRefine(const vector<int>& rid);

	//void ShootRay_Edge(int edid, int pid, double kv[4]);//take an edge as reference
	//void ShootRay_Face(int fcid, int pid, double kv[4]);//take a face as reference
	void VisualizeTMesh(string fn);
	bool CheckSubKnotVector(const double ku1[5], const double kv1[5], const double ku2[5], const double kv2[5]);//check whether v1 is a subvector of v2
	bool CheckSubKnotVector(const array<double,5>& ku1, const array<double,5>& kv1, const array<double,5>& ku2, const array<double,5>& kv2);

	void InitialConnect();
	void UpdateConnect();
	void FindEdgeTopoDirec();
	void FindKnotInterval();
	void UpdateKnotInterval();
	void ShootRay(int pid, int edid, double kv[4], int& trun_flag);//start from point pid with the edge edid
	void FindIEN_Unstruct();
	void FindIEN_Invalid();//include invalid elements
	void Update_IEN();
	void FindNextRing(const vector<int>& pr0, const vector<int>& er0, vector<int>& pr1, vector<int>& er1, vector<int>& pr1_pref, vector<int>& pr1_eref);
	void SetLocalCoorSystem();
	void FindRotateAndUVCoor(int pref,int rot_ref,const array<double,2>& uv_ref,int eid,int pid,int& rot,array<double,2>& uv);
	void FindLocalKnotVector(int id,int rot,const array<double,2>& uv,array<double,5>& ku,array<double,5>& kv);
	bool CheckSupport(const array<double,2>& u,const array<double,2>& v,const array<double,5>& ku,const array<double,5>& kv);
	void UpdatePatchCP_Unstruct(int eid);

	void Truncation();
	//void FindChildren();

	void ElementBasis(int eid, double u, double v, vector<double>& Nt, vector<array<double,2>>& dNdt);
	void ElementBasis_Regular(int eid, double u, double v, vector<double>& Nt, vector<array<double,2>>& dNdt);
	void ElementBasis_Irregular(int eid, double u, double v, vector<double>& Nt, vector<array<double,2>>& dNdt);
	void SetBezierMatIrrPatch();
	void FindPatchKnotVector_Irr(int eid, vector<array<double,5>>& patch_ku, vector<array<double,5>>& patch_kv);
	void SurfacePointMap(int eid, double u, double v, array<double,3>& pt, array<double,3>& norm);
	void ElementRefine_Unstruct_4(int eid);
	void ElementRefine_Irregular_4(int eid);
	void ElementRefine_Invalid_4(int eid);
	void ElementRefine_Unstruct_2(int eid, int dir);
	void ElementRefine_Unstruct_b(int eid);
	//void ElementRefine_Unstruct_Topo_Geom_4(int eid);
	//void ElementRefine_Unstruct_Topo_Geom_2(int eid, int dir);
	void Topo_Refine_Unstruct(const vector<int>& rfid, const vector<int>& rftype, vector<int>& rfid_more, vector<int>& rftype_more);
	void Geom_Refine_Unstruct(const vector<int>& rfid, const vector<int>& rftype);

	void BezierExtract_Unstruct(vector<BezierElement2D>& bzmesh);
	void BezierElementExtract_Unstruct(int eid,vector<BezierElement2D>& bzmesh);
	void BezierElementExtract_Unstruct_Irr(int eid,vector<BezierElement2D>& bzmesh);
	void BezierUnit_Unstruct(int eid,array<double,4>& kts,vector<BezierElement2D>& bzmesh);
	void BezierUnit_Unstruct_Irr(int eid,array<double,4>& kts,vector<BezierElement2D>& bzmesh);
	void BezierUnit_Unstruct_Trun(int eid,array<double,4>& kts,vector<BezierElement2D>& bzmesh);
	void BezierUnit_Unstruct_Irr_Trun(int eid,array<double,4>& kts,vector<BezierElement2D>& bzmesh);
	void BezierFinder_Unstruct(int eid,vector<array<double,4>>& be);
	void BezierVTK_Unstruct(string fn,vector<BezierElement2D>& bzmesh);
	void BezierControlMesh_Unstruct(string fn,vector<BezierElement2D>& bzmesh);

	void run_surf_XP(string fn, vector<int> ids);
	void OutputMesh(const vector<BezierElement2D>& bzmesh, string fn);
	void runXP_Laplace();
	void runXP_complex(string fn);

	void SetBezier3TranMat(int N, vector<vector<double>>& bmat);
	void SetBezier4TranMat(int N, vector<vector<double>>& bmat);
	void SetBezier4TranMatOP(int N, vector<vector<double>>& bmat);
	void VisualizeSurface(string fn);
	void CapIndex_Loc2Glb(int N, vector<vector<int>>& ICN);
	void BuildGeomMat(int N, const vector<int>& bloc, const vector<vector<double>>& bmat, const vector<vector<int>>& ICN, vector<vector<double>>& gmat, vector<double>& gvec);
	void BuildFairMat(int N, const vector<int>& bloc, const vector<vector<double>>& bmat, const vector<vector<int>>& ICN, vector<vector<double>>& fmat, vector<double>& fvec);
	void PreConditionCoef(int N, const vector<int>& bloc, const vector<vector<double>>& bmat, const vector<vector<int>>& ICN, vector<double>& coef);
	void Convert2MatlabData(const vector<vector<double>>& mat, vector<double>& row_id, vector<double>& col_id, vector<double>& coef);
	void SolveOptimizeBezierMat(int N, vector<vector<double>>& bmatop);
	void Topology_Irregular();

	void SetProblem_surf_XP(string fn);
	void SetLshapeProblem_XP(string fn);
	void SetProblem_complex(string fn);

	//void SurfRefine();
	//void LshapeRefine();
	//void LshapeRefine_H1();
	//void LshapeRefine_XP();

	void Identify_Invalid_Elements(vector<int>& rid);
	void Identify_More_Elements(vector<int>& rid);
	void SetBounaryPoints(vector<int>& pid, vector<double>& disp, int& ncp);

	void Refine_Surf_Select_0(vector<int>& rfid, vector<int>& rftype);
	// void Refine_Surf_Select_0(vector<int>& rfid, vector<int>& rftype);
	void Refine_Surf_Test_0(vector<int> ids);

	//void Topo_Refine_Struct(const vector<int>& rfid, const vector<int>& rftype, vector<int>& rfid_more, vector<int>& rftype_more);
	//void Geom_Refine_Struct(const vector<int>& rfid, const vector<int>& rftype);
	//void ElementRefine_Struct_4(int eid);//subdivide into 4
	//void ElementRefine_Struct_3(int eid, int dir);//subdivide into 3
	//void ElementRefine_Struct_2(int eid, int dir);//subdivide into 2
	//void ElementRefine_Struct_b(int eid);//subdivide into 2, boundary element

	//vector<int> paid;//active control points
	vector<int> eaid;//active elements

private:
	vector<Vertex2D> cp;//control points
	vector<Element2D> tmesh;//elements in T-mesh
	vector<Edge> tmedge;

	//vector<pair<int,double>> uanc;
	//vector<pair<int,double>> vanc;
	//vector<vector<int>> uedge;
	//vector<vector<int>> vedge;

	unsigned int npt_old;
	unsigned int nel_old;
};

#endif
