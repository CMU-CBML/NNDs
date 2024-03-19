#ifndef TTSP_3D_H
#define TTSP_3D_H

#include <vector>
#include <utility>
#include <string>
//#include "T_mesh.h"
#include "BasicDataStructure.h"
using namespace std;

class TruncatedTspline_3D
{
public:
	TruncatedTspline_3D();
	void CreateUniformCube(string fn);
	void VisualizeControlMesh(string fn);
	double PartionOfUnity(int eid, const array<double, 3>& u);
	void CollectActives();

	void VisualizeTMesh(string fn);
	void VisualizeFaceMesh(string fn);
	bool CheckSubKnotVector(const double ku1[5], const double kv1[5], const double ku2[5], const double kv2[5]);//check whether v1 is a subvector of v2
	bool CheckSubKnotVector(const array<double,5>& ku1, const array<double,5>& kv1, const array<double,5>& ku2, const array<double,5>& kv2);

	void InitialConnect();
	void InitialRotate();
	void getElementRotate(vector<Matrix3d>& mat);
	void getCornerRotate(int hxid, int uvw_loc[3], int uvw_ref[3], Matrix3d& mat);
	void UpdateConnect();
	void FindEdgeTopoDirec();
	void EdgeConnect(int ed0, int pid0, int& type1, int& ed1);
	//void FindFaceTopoDirec();
	void FaceConnect(int fc0, int ed0, int& type1, int& fc1);
	void setReference();
	void FindKnotInterval();
	//void UpdateKnotInterval();
	void ShootRay(int pid, int edid, double kv[4], int& trun_flag);//start from point pid with the edge edid
	void SetLocalCoorSystem();
	void Find_Neighbor_Rot(int hxid, int pid, Matrix3d& rot);
	void getElementRotate_Unit(int loc, Matrix3d& rot);
	//void FindIEN_Unstruct();
	void FindIEN_PatchKV();
	void Update_IEN();
	void FindNextRing(const vector<int>& pr0, const vector<int>& er0, vector<int>& pr1, vector<int>& er1, vector<int>& pr1_pref, vector<int>& pr1_eref);
	void TranslateLCS(int pref,Matrix4d& lcs_ref,int eid,int pid,Matrix4d& lcs);
	void getLCS_inverse(const Matrix4d& mat_in, Matrix4d& mat_out);
	void FindLocalKnotVector(int pid, Matrix4d lcs, array<double, 5>& ku, array<double, 5>& kv, array<double, 5>& kw);
	bool CheckSupport(const array<double, 2>& u, const array<double, 2>& v, const array<double, 2>& w, const array<double, 5>& ku, const array<double, 5>& kv, const array<double, 5>& kw);
	void SetBezierMatIrrPatch(int eid);
	void AllBezierPatch();
	void Para2Physical(int eid, const array<double, 3>& u, array<double, 3>& pt);
	void UpdatePatchCP_Unstruct(int eid);

	void Truncation();
	bool CheckSubKnotVector(const array<double, 5>& ku1, const array<double, 5>& kv1, const array<double, 5>& kw1,
		const array<double, 5>& ku2, const array<double, 5>& kv2, const array<double, 5>& kw2);
	//void FindChildren();

	void ElementBasis(int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);
	//void ElementBasis_Bezier(int eid, const array<double,3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);
	void ElementBasis_Regular(int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);
	void ElementBasis_Irregular(int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);
	void FindPatchKnotVector_Irr(int eid, vector<array<double,5>>& patch_ku, vector<array<double,5>>& patch_kv);
	void SurfacePointMap(int eid, double u, double v, array<double,3>& pt, array<double,3>& norm);

	void ElementSubdivide_8(int eid);
	void ElementSubdivide_5(int eid, int dir);
	void ElementSubdivide_4(int eid, int dir);
	void ElementSubdivide_2(int eid, int dir);//add a yz plane if dir=0, xz for dir=1, xy for dir=2

	void FaceSubdivision_4(int fcid);
	void FaceSubdivision_2(int fcid, int dir);
	void FaceSubdivision_24(int fcid, int dir);//subdivide a face that has been split
	void EdgeSubdivision_2(int edid);

	void SolidFaceDirection(int eid, int fcid, int& dir, int& pos);
	void EdgeIndex_in_Face(Face3D& fc, const vector<int>& edid);
	void EdgeFaceIndex_in_Solid(Element3D& hex, const vector<int>& edid, const vector<int>& fcid);
	void Index_Direction(int eid, int vloc[8], int edloc[12], int fcloc[6], int dir);

	void Construct_NewEdge(int pt[2], int lev, double len);
	void Construct_NewFace(int cnct[4], int edid[4], double lv, double len);
	void Construct_BaseFace_Subdv(int fcid, int dir, int& ed_base, int fc_base[2]);
	void Construct_BaseFace_Bisct(int fcid, int dir);

	void IdentifyTarget(vector<int>& target);
	void IdentifyAddition(vector<int>& target);
	bool IdentifyAddition();
	void IdentifyPropagate(vector<int>& id_in, vector<int>& id_out);
	void RefineTopology();
	void RefineGeometry();
	void RefineGeometry_v1();

	void CalPatchCP_Regular(int eid);
	bool CheckFullRefine(const vector<array<double, 3>>& spt, const array<double, 5>& motu, const array<double, 5>& motv, const array<double, 5>& motw,
		const vector<array<double, 5>>& chdu, const vector<array<double, 5>>& chdv, const vector<array<double, 5>>& chdw, const vector<double>& coef);
	//void CheckTruncation();
	void UpdateGeometry(vector<int>& rfid);

	bool Identify_FaceFace(int eid);
	bool Identify_FaceEdge(int eid);
	bool Template_FaceFace();
	bool Template_FaceEdge();

	void RefineTest_1();
	void RefineTest_2();

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

	//void run_surf_XP(string fn);
	//void runXP_Laplace();
	void run(string fn);
	void runh(string fn);

	void SetBezier3TranMat(int N, vector<vector<double>>& bmat);
	void VisualizeSolidVTK(string fn);

	void SetProblem_surf_XP(string fn);
	void SetLshapeProblem_XP(string fn);
	void SetProblem(string fn);
	void SetSharpFeature();
	void SetDomain();

	//void SurfRefine();
	//void LshapeRefine();
	//void LshapeRefine_H1();
	//void LshapeRefine_XP();

	void Identify_Invalid_Elements(vector<int>& rid);
	void Identify_More_Elements(vector<int>& rid);
	void SetBounaryPoints(vector<int>& pid, vector<double>& disp, int& ncp);

	void Refine_Surf_Select_0(vector<int>& rfid, vector<int>& rftype);
	void Refine_Surf_Test_0();

	void ElementDomain();
	void PatchRefine_Regular(int lev, int eid);
	void BisectKnotInterval(const array<double, 5>& kv_in, vector<double>& kv_out);
	void PatchRefine_Irregular(int lev, int eid);
	void CatmullClark(int lev, int eid);
	void PatchRefine_Boundary(int lev, int eid);
	void ConstructBezierBasis(int lev, int eid);
	void ConstructBezierBasis_1(int lev, int eid);
	void ConstructBezierBasis_Boundary(int lev, int eid);
	void ConstructBezierBasis_Feature(int lev, int eid);
	void ConstructConnect(int lev);
	void ConstructFaceEdge(int lev, int eid);
	void Refine_Ghost(const vector<array<int, 2>>& rfid);

	void Selection(int lev);
	void SetSupport(int lev);
	void Select();
	void Truncate(int lev);
	void Basis_Regular(int lev, int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);
	void Basis_Irregular(int lev, int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);
	void Identify_Pseudo(vector<array<int,2>>& rfid);
	void Identify_Pseudo_1(vector<array<int, 2>>& rfid);
	void Identify_Test_1(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void Identify_Test_2(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void Identify_Test_3(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void Identify_Test_4(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);//boundary
	void Refine(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	double BasisSum(int lev, int eid, const array<double,3>& u);

	void OutputCM(int lev, string fn);
	void OutputFace(int lev, string fn);
	void OutputEdge(int lev, string fn);
	void OutputGeom(int lev, string fn);
	void GeomMap(int lev, int eid, const array<double,3>& u, array<double,3>& pt);
	void GeomMap_Bezier(int lev, int eid, const array<double, 3>& u, array<double, 3>& pt);
	void OutputGeom_All(string fn);

	void GeomMap_Lev(int lev, int eid, const array<double, 3>& u, array<double, 3>& pt);
	double BasisSum_Lev(int lev, int eid, const array<double, 3>& u);
	void Basis_Regular_Lev(int lev, int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);
	void Basis_Irregular_Lev(int lev, int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);

	void AllBezierLev(int lev);
	void AnalysisInterface_Elastic(vector<BezierElement3D>& bzmesh, vector<int>& DrchBC, vector<double>& gh);
	void AnalysisInterface_Poisson(vector<BezierElement3D>& bzmesh, vector<int>& DrchBC, vector<double>& gh);//only for cube domain [0,1]^3
	void AnalysisInterface_Poisson_1(vector<BezierElement3D>& bzmesh, vector<int>& DrchBC, vector<double>& gh);//for arbitrary shapes
	void AnalysisInterface_Laplace(const vector<array<int, 2>>& pbc, const vector<double>& pdisp, vector<BezierElement3D>& bzmesh, vector<int>& DrchBC, vector<double>& gh);
	void AnalysisInterface_LeastSquare(vector<BezierElement3D>& bzmesh, vector<int>& DrchBC, vector<double>& gh);
	double SpecifyDirichBC(double x[3]);
	double SpecifyDirichBC_2(double x[3]);
	double SpecifyDirichBC_3(double x[3]);
	double SpecifyDirichBC_4(double x[3]);
	double SpecifyDirichBC_5(double x[3]);
	double SpecifyDirichBC_6(double x[3]);
	double SpecifyDirichBC_7(double x[3]);
	double SpecifyDirichBC_8(double x[3]);
	void Identify_Poisson(const vector<array<double, 2>>& ehid, const vector<double>& err, vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void Identify_Poisson_1(const vector<array<double, 2>>& ehid, const vector<double>& err, vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);//this is new
	void Identify_Laplace(const vector<array<double, 2>>& ehid, const vector<double>& err, vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void Identify_LeastSquare(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);//sphere
	void Identify_LeastSquare_Line(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);//line
	void SetInitialBC(vector<array<int, 2>>& ebc, vector<double>& edisp, vector<array<int, 2>>& pbc, vector<double>& pdisp);
	void SetBC(vector<array<int, 2>>& ebc, vector<double>& edisp, vector<array<int, 2>>& pbc, vector<double>& pdisp);
	void FittingBC(vector<int>& DrchBC, vector<double>& gh);

	void BezierPoints_ini();
	void BezierPoints_Refine(int lev, int pos, const vector<MatrixXd>& bsmat);//input is the father element
	void BezierSubdivMatrix(vector<MatrixXd>& bsmat);
	void OutputRefineID(string fn, const vector<array<int,2>>& rfid, const vector<array<int,2>>& gst);
	void InputRefineID(string fn, vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void GetRemoveRegion(double xy[3][2]);
	void OutputRemoveCM(string fn, double xy[3][2]);
	void CreateRemoveView(string fn_in, string fn_out);

	void OutputBasis(int lev, int pid, string fn);
	void PillowCube(string fn);
	void VisualizeBezier(const vector<BezierElement3D>& bzmesh, string fn);
	void InputCheck(string fn);
	void MeshRepair(string fn_in, string fn_out);
	void ReportXP();

	void GlobalRefine(int niter);//this is just using functions of local refinement, not efficient
	double MaxElementSize(const vector<BezierElement3D>& bzmesh);

	//true global refinement, can also be included in local refinement
	void InitializeMesh(string fn);
	void SetSharpFeature_1();//use tmesh, tmface, tmedge
	void Global_Subdivide();
	void Global_Subdivide_Simple();
	void BuildBasisFunction();
	void OutputCM(string fn);
	void OutputFace(string fn);
	void OutputEdge(string fn);
	void run_Global(int nref, string fn);

	void SetDomainRange(double rg[3][2], double nm[3], double& a);

	//reimplement local refinement for efficiency
	void Refine_eff(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void PatchRefine_eff(int lev, int eid);

	//for pipeline paper
	void PipelineDataProcess(string fn_in, string fn_out);
	void PipelineBezierExtract(vector<BezierElement3D>& bzmesh);

	//B-splines for test
	void CreateBsplines(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void CreateBsplines(int nex[3]);
	void Bsplines_Refine();//global refinement
	void Bspline_BezierExtract(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void Bspline_BezierExtract_fit(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void FittingBC_Bsplines(const vector<BezierElement3D>& bzall, vector<int>& DrchBC, vector<double>& gh);
	void SetDomainRange_Bsplines(double rg[3][2], double nm[3], double& a);

	//vector<int> paid;//active control points
	vector<int> haid;//active hex elements
	vector<int> faid;//active faces
	vector<int> eaid;//active edges

private:
	vector<Vertex3D> cp;//control points
	vector<Element3D> tmesh;//elements in T-mesh
	vector<Edge3D> tmedge;
	vector<Face3D> tmface;
	unsigned int npt_old;
	unsigned int nel_old;
	vector<vector<Vertex3D>> hcp;
	vector<vector<Element3D>> hmesh;
	vector<vector<Edge3D>> hedge;
	vector<vector<Face3D>> hface;
	vector<vector<double>> kvec;

	//used for solution
	double dmrg[3][2];
	double nmpl[3];
	double acoef;
	double dmlen[3];
};

#endif