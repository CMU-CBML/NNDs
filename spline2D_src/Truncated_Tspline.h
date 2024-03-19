#ifndef TRUNCATED_TSPLINE_H
#define TRUNCATED_TSPLINE_H

#include <vector>
#include <utility>
#include <string>
//#include "T_mesh.h"
#include "BasicDataStructure.h"
using namespace std;

class TruncatedTspline
{
public:
	TruncatedTspline();
	void ReadTMesh(string fn);
	void FindPatchKnotVector(int eid, vector<double> &ku, vector<double> &kv);
	void Refine(const vector<int> &rf, int id = 0); //0 for both directions, 1 for u, 2 for v
	void ElementSubdiv(int eid);
	void BasisEvaluate(int eid, double u, double v, vector<double> &Nt);
	void Parmt2Phys(int eid, double u, double v, Vertex &pt);
	void Parmt2Phys_v0(int eid, double u, double v, Vertex &pt);
	void VisualizeVTK(string fn);
	void VisualizeControlMesh(string fn, int outflag = 0);
	double PartitionOfUnity(int eid, double u, double v);
	void PseudoProblem();
	void Test();
	void CreateInputMesh(string fn, int neu, int nev);
	void SetPseudoProblem(string fn);
	void ElementSubdivTest(int eid);
	void IdentifyTest(vector<int> &rid);
	void CollectActives();
	void BezierExtract(vector<BezierElement> &bzmesh);
	void BezierElementExtract(int eid, vector<BezierElement> &bzmesh);
	void BezierUnit(vector<int> &pid, array<double, 4> &kts, vector<BezierElement> &bzmesh);
	void BezierUnit_glb(vector<int> &pid, array<double, 4> &kts, vector<BezierElement> &bzmesh);
	void BezierFinder(int eid, vector<array<double, 4>> &be);
	void BezierFinder_1(int eid, vector<array<double, 4>> &be);
	void BezierVTK(string fn, vector<BezierElement> &bzmesh);
	bool CheckSameKnotVector(const Vertex &p1, const Vertex &p2);
	void RefineTest();
	void RefineTest1();
	void RefineTest2();
	void RefineTest3();
	void RefineTest3_1();
	void RefineTest4();
	void RefineTest5();
	void RefineTest6();
	void RefineTest7();
	void RefineTest8();
	void RefineTest9();
	void RefineTest10();
	void RefineTest10_1();
	void Refine1();
	void Refine2();
	void Refine3();
	void Truncation();
	void ElementRefine(int eid);
	void ElementRefine1(int eid);
	void UpdateIEN();
	void TopologyRefine(const vector<int> &rid);
	void TopologyRefine_1(const vector<int> &rid);
	void TopologyRefine_2(const vector<int> &rid);
	void ElementTopologyRefine(int eid);
	void InsertionCheck(int);
	//void TopologyRefine_2(const vector<int>& rid);
	void StrongBalanceCheck(const vector<int> &rid, vector<int> &ridsb);
	void ElementRefine_Square_4(int eid);		   //subdivide into 4
	void ElementRefine_Square_3(int eid);		   //subdivide into 3
	void ElementRefine_Square_2(int eid);		   //subdivide into 2
	void ElementRefine_Square_3(int eid, int dir); //subdivide into 3
	void ElementRefine_Square_2(int eid, int dir); //subdivide into 2
	void ElementRefine_Rectangular(int eid);	   //subdivide into 2
	void ElementRefine_Boundary(int eid);		   //subdivide into 2
	void ElementRefine_Singular(int eid);		   //subdivide into 4
	void OneTjunctionCheck(vector<int> &ridtp);	   //one T-junction requirement
	void TjuncExtentCheck();
	void TjuncExtentCheck_1(vector<int> &ridtjx);
	void TjuncExtentCheck_2(vector<int> &ridtjx);
	void UpdateTopology();
	void InitialTopology(); //assume input as a coarse quad mesh
	void StrongBalanceRefine(const vector<int> &rid);
	void TargetRefine(const vector<int> &rid);
	//void TargetRefine(const vector<int>& rid, const vector<int>& dir_id);
	void OneTjunctionRefine(const vector<int> &rid);
	void TjuncExtentRefine(const vector<int> &ridtjx);
	void InitializeTopologyDirect();
	void FindTopologyDirect();
	void FindKnotInterval();
	void ShootRay_Edge(int edid, int pid, double kv[4]); //take an edge as reference
	void ShootRay_Face(int fcid, int pid, double kv[4]); //take a face as reference
	void FindIEN();
	int FindLocalUV(int eid, int pid); //UV of pid wrt eid
	void RelativeElementDirect();
	void ElementKnotVectors(int eid);

	void GeometryRefine();
	void UpdateControlPoints();
	void UpdateControlPoints_1();
	void UpdateControlPoints_2();
	void UpdateControlPoints_3();
	void UpdateControlPoints_4();
	void UpdateKnotVectors();
	void CheckTruncation();
	void CheckTruncation_1();
	double EvaluateTrunBF(int pid);
	double EvaluateTrunBF1(int pid);
	double EvaluateTrunBF2(int pid);
	void Refine_Addition(const vector<int> &rid);
	void Refine_Target(const vector<int> &rid);

	void AssignIndex_NewFacePoint(int fcid, Vertex &pt);
	void AssignIndex_NewEdgePoint(int edid, Vertex &pt);
	void SortEdge();
	void FindLocalKnotVectors();
	void FindIENglb();

	void TopologyRefineTest();
	void TopologyRefineTest1();
	void TopologyRefineTest2();
	void TopologyRefineTest3();
	void TopologyRefineTest4();
	void IdentifyTest1(vector<int> &rid);
	void IdentifyTest2(vector<int> &rid);
	void VisualizeTMesh(string fn);
	void OutputNumbers();

	void ElementSubdivide_4(int eid); //subdivide into 4
	void ElementSubdivide_2(int eid); //subdivide into 4
	void UpdateControlPoints_v0();
	void UpdateIEN_v0();
	void FindEdge_v0();
	void RefineTest_v0_1();

	bool CheckSubKnotVector(const double ku1[5], const double kv1[5], const double ku2[5], const double kv2[5]); //check whether v1 is a subvector of v2
	bool CheckSubKnotVector(const array<double, 5> &ku1, const array<double, 5> &kv1, const array<double, 5> &ku2, const array<double, 5> &kv2);
	int CheckKnotInsertion(const double ku1[5], const double kv1[5], const double ku2[5], const double kv2[5]); //-1 no insertion, 0 u, 1 v, 2 uv
	bool NewKnots(const double kv1[5], const double kv2[5]);													//kv1 is a result of knot insertion of kv2
	void PlotBasisFunctions(string fn, vector<int> &pid);
	void PlotTrunBasisFunctions(string fn, vector<int> &pid);
	void FunctionMagnitude(int eid, int pid, double u, double v, Vertex &pt);
	void TrunFunctionMagnitude(int eid, int pid, double u, double v, Vertex &pt);

	void InitialConnect();
	void InitialConnect_1();
	void InitialConnect_2();
	void UpdateConnect();
	void UpdateConnect_1();
	void FindEdgeTopoDirec();
	void FindEdgeTopoDirec_1();
	void FindKnotInterval_1();
	void UpdateKnotInterval_1();
	void ShootRay(int pid, int edid, double kv[4]); //start from point pid with the edge edid
	//double FindNextInterval_Edge(int edid, int end, int& flag, int& id, int& dir);//flag=0 means next is edge, flag=1 means next is face, flag=2 XP, 3 for boundary
	//double FindNextInterval_Face(int fcid, int dir, int& flag, int& id, int& dir);
	void FindIEN_1();
	void FindIEN_2();
	void FindIEN_3(); //include invalid elements
	void Update_IEN_3();
	//void AdditionIEN_TrunPatch();
	void FindNextRing(const vector<int> &pr0, const vector<int> &er0, vector<int> &pr1, vector<int> &er1, vector<int> &pr1_pref, vector<int> &pr1_eref);
	void SetLocalCoorSystem();
	void FindRotateAndUVCoor(int pref, int rot_ref, const array<double, 2> &uv_ref, int eid, int pid, int &rot, array<double, 2> &uv);
	void FindLocalKnotVector_1(int id, int rot, const array<double, 2> &uv, array<double, 5> &ku, array<double, 5> &kv);
	bool CheckSupport(const array<double, 2> &u, const array<double, 2> &v, const array<double, 5> &ku, const array<double, 5> &kv);
	void FindPatchKnotVector_1();
	void UpdatePatchCP_Unstruct(int eid);
	void UpdatePatchCP_Unstruct_1(int eid);
	void UpdatePatchCP_Unstruct_2(int eid);
	bool CheckFullChildren(const array<double, 5> &motu, const array<double, 5> &motv, const vector<array<double, 5>> &chdu, const vector<array<double, 5>> &chdv, const vector<double> &coef);
	bool CheckFullChildren_1(const vector<array<double, 2>> &spt, const array<double, 5> &motu, const array<double, 5> &motv, const vector<array<double, 5>> &chdu, const vector<array<double, 5>> &chdv, const vector<double> &coef);
	double CheckFullChildren_2(const vector<array<double, 2>> &spt, const array<double, 5> &motu, const array<double, 5> &motv, const vector<array<double, 5>> &chdu, const vector<array<double, 5>> &chdv, const vector<double> &coef);
	void SetProblem2(string fn);
	void SetProblem_complex(string fn);
	void ElementBasis(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);
	void ElementBasis_Regular(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);
	void ElementBasis_Irregular(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);
	void SetBezierMatIrrPatch();
	void FindPatchKnotVector_Irr(int eid, vector<array<double, 5>> &patch_ku, vector<array<double, 5>> &patch_kv);
	void SurfacePointMap(int eid, double u, double v, array<double, 3> &pt, array<double, 3> &norm);
	void ElementRefine_Unstruct_4(int eid);
	void ElementRefine_Irregular_4(int eid);
	void ElementRefine_Invalid_4(int eid);
	void ElementRefine_Unstruct_2(int eid, int dir);
	void ElementRefine_Unstruct_b(int eid);
	void ElementRefine_Unstruct_Topo_Geom_4(int eid);
	void ElementRefine_Unstruct_Topo_Geom_2(int eid, int dir);
	void Refine_Unstruct(const vector<int> &rfid, const vector<int> &rftype); //type 0 for regular element into 4, 1 for regular element into 2 in dir=0, 2 for in dir=1, 3 for boundary element, 4 for irregular element into 4
	void Topo_Refine_Unstruct(const vector<int> &rfid, const vector<int> &rftype, vector<int> &rfid_more, vector<int> &rftype_more);
	void Geom_Refine_Unstruct(const vector<int> &rfid, const vector<int> &rftype);
	void Topo_Refine_Unstruct_glb(const vector<int> &rfid, const vector<int> &rftype);
	void Geom_Refine_Unstruct_glb(const vector<int> &rfid, const vector<int> &rftype);
	void Validate_Tmesh();
	void BezierExtract_Unstruct(vector<BezierElement> &bzmesh);
	void BezierElementExtract_Unstruct(int eid, vector<BezierElement> &bzmesh);
	void BezierElementExtract_Unstruct_Irr(int eid, vector<BezierElement> &bzmesh);
	void BezierUnit_Unstruct(int eid, array<double, 4> &kts, vector<BezierElement> &bzmesh);
	void BezierUnit_Unstruct_Irr(int eid, array<double, 4> &kts, vector<BezierElement> &bzmesh);
	void BezierUnit_Unstruct_Trun(int eid, array<double, 4> &kts, vector<BezierElement> &bzmesh);
	void BezierUnit_Unstruct_Irr_Trun(int eid, array<double, 4> &kts, vector<BezierElement> &bzmesh);
	void BezierFinder_Unstruct(int eid, vector<array<double, 4>> &be);
	void BezierVTK_Unstruct(string fn, vector<BezierElement> &bzmesh);
	void BezierControlMesh_Unstruct(string fn, vector<BezierElement> &bzmesh);
	void run_surf(string fn);
	void run_surf_1(string fn);
	void run_Laplace(string fn);
	void runXP();
	void runXP_Laplace();
	void runXP_complex(string fn);
	void run_PlateHole(string fn);
	void Refine_Global();
	void NeumannBC_PlateHole(vector<BezierElement> &bzmesh);

	void run_PlateHole_XP(string fn);
	void Refine_Global_XP();
	void CircleError(const vector<double> &kv, const vector<array<double, 4>> &ptsw, vector<double> &err);

	void CreateMesh_XP_2(string fn0, string fn1);
	void CreateMesh_XP_2_Pillow(string fn0, string fn1);
	void SetBezier3TranMat(int N, vector<vector<double>> &bmat);
	void SetBezier4TranMat(int N, vector<vector<double>> &bmat);
	void SetBezier4TranMatOP(int N, vector<vector<double>> &bmat);
	void SurfacePointCal(double u, double v, const vector<int> &IEN, const vector<vector<double>> &bmat, array<double, 3> &pcoor, array<double, 3> &norm);
	void SurfacePointCal1(double u, double v, const vector<int> &IEN, const vector<vector<double>> &bmat, array<double, 3> &pcoor, array<double, 3> &norm);
	void VisualizeSurface(string fn);
	void SetProblem1(string fn);
	void GlueIrregularPatches4OP(vector<vector<int>> &pats);
	void CapIndex_Loc2Glb(int N, vector<vector<int>> &ICN);
	void BuildGeomMat(int N, const vector<int> &bloc, const vector<vector<double>> &bmat, const vector<vector<int>> &ICN, vector<vector<double>> &gmat, vector<double> &gvec);
	void BuildFairMat(int N, const vector<int> &bloc, const vector<vector<double>> &bmat, const vector<vector<int>> &ICN, vector<vector<double>> &fmat, vector<double> &fvec);
	void PreConditionCoef(int N, const vector<int> &bloc, const vector<vector<double>> &bmat, const vector<vector<int>> &ICN, vector<double> &coef);
	void Convert2MatlabData(const vector<vector<double>> &mat, vector<double> &row_id, vector<double> &col_id, vector<double> &coef);
	void SolveOptimizeBezierMat(int N, vector<vector<double>> &bmatop);
	void Topology_Irregular();

	void SetProblem_surf(string fn);
	void SetLshapeProblem(string fn);
	void SetLshapeProblemLocalFeature(string fn);
	void SetLshapeProblem_XP(string fn);
	void SurfRefine();
	void SurfRefine_1();
	void LshapeRefine();
	void LshapeRefine_H1();
	void LshapeRefine_XP();

	void SetProblem_PlateHole(string fn);
	void SetProblem_PlateHole_XP(string fn);
	void Refine_PlateHole();
	void Refine_PlateHole_XP();
	void InitialTopology_PlateHole();

	void Identify_Invalid_Elements(vector<int> &rid);
	void Identify_More_Elements(vector<int> &rid);
	void SetBounaryPoints(vector<int> &pid, vector<double> &disp, int &ncp);

	void Topo_Refine_Struct(const vector<int> &rfid, const vector<int> &rftype, vector<int> &rfid_more, vector<int> &rftype_more);
	void Geom_Refine_Struct(const vector<int> &rfid, const vector<int> &rftype);
	void ElementRefine_Struct_4(int eid);		   //subdivide into 4
	void ElementRefine_Struct_3(int eid, int dir); //subdivide into 3
	void ElementRefine_Struct_2(int eid, int dir); //subdivide into 2
	void ElementRefine_Struct_b(int eid);		   //subdivide into 2, boundary element

	void RepairMesh();

	void run_PatchTest_HighOrder(string fn);
	void set_PatchTest_HighOrder(string fn);

	//Xin Li's subdivision idea
	void RunHybridCC(string fn, int nrf);
	void Preprocess(string fn);
	void ReadInput(string fn);
	void Scale2Unit();
	void GetDual();
	void BuildConnectDual();
	void VisualizeEdge(string fn);
	void VisualizeQuad(string fn);
	void BuildTsplinesDual();
	void UpdateConectDual();
	void FindEdgeTopoDirecDual();
	void FindKnotIntervalDual();
	void ShootRayDual(int pid, int edid, double kv[4]);
	void UpdateKnotIntervalDual();
	void SetLocalCoorSystemDual();
	void FindIENDual();
	void FindElementIEN_RegDual(int eid);
	void FindElementIEN_IrrDual(int eid);
	void FindPatchKVDual(int eid, double ku[8], double kv[8]);
	void UpdateIENDual();
	void VisualizeGeomDual(string fn);
	void GetSubMatrixDual(int nv, vector<MatrixXd> &smat);
	void GetSubMatrixDual_All(int nv, MatrixXd &amat);
	void GetSubMatrixDual_Bar(int nv, MatrixXd &amat);
	void GetSubPatchKVDual(vector<array<double, 8>> &subku, vector<array<double, 8>> &subkv);
	void ElementBasisDual(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);
	void ElementBasis_IrregularDual(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);
	void GlobalRefineDual();
	void InitializeRefine();
	void ElementRefineReg_Dual(int eid);
	void ElementRefineIrr_Dual(int eid);
	void ElementRefineBnd_Dual(int eid);
	void UpdateCP();
	void BezierExtractDual(vector<BezierElement> &bzmesh);
	void ElementBezierExtractReg(int eid, BezierElement &bzel);
	void ElementBezierExtractIrr(int eid, BezierElement &bzel);
	void DirichletBC(vector<int> &IDBC, vector<double> &gh);
	void DirichletBC_LeastSquare(vector<BezierElement> &bzmesh, vector<int> &IDBC, vector<double> &gh);

	//Catmull-Clark
	void RunCC(string fn_in, int nrf);
	void BuildCC();
	void FindIENCC();
	void FindElementIEN_IrrCC(int eid);
	void GetSubMatrixCC_All(int nv, MatrixXd &amat, MatrixXd &abar);
	void UpdateIENCC();
	void GlobalRefineCC();
	void ElementRefineIrr_CC(int eid);
	void VisualizeGeomCC(string fn);
	void ElementBasisCC(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);
	void ElementBasis_IrregularCC(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);

	//Uunstructured T-splines with degenerated C1
	void Run_UTSP(string fn_in, string fn_out);
	void Run_UTSP_1(string fn_in, string fn_out, vector<BezierElement> &bzmesh, vector<int> &IDBC, vector<double> &gh);
	void Run_UTSP_AS(string fn_in, string fn_out, int ngref, int nlref, int fit_flag, string fn_lref, vector<BezierElement> &bzmesh, vector<int> &IDBC, vector<double> &gh);
	void Run_UTSP_Design(string fn_in, string fn_out);
	void Run_UTSP_AS_Scale(string fn_in, string fn_out, int ngref, int nlref, int fit_flag, string fn_lref, vector<BezierElement> &bzmesh,
						   vector<int> &IDBC, vector<double> &gh);

	void Run_UTSP_AS_Scale_Honda(string fn_in, string fn_out, int ngref, int nlref, string fn_lref, vector<BezierElement> &bzmesh,
								 vector<int> &IDBC, vector<double> &gh);
	void HemisphereLocalRefine_Honda(string fpre, int nlev);
	void ReadRef4ID_Honda(string fn, vector<int> &rfid);
	void Run_UTSP_AS_Honda(string fn_in, string fn_out, int ngref, int nlref, int fit_flag, string fn_lref,
						   vector<BezierElement> &bzmesh, vector<int> &IDBC, vector<double> &gh);
	void Run_UTSP_AS_GEM(string fn_in, string fn_out, int ngref, int nlref, int fit_flag, string fn_lref,
						 vector<BezierElement> &bzmesh, vector<int> &IDBC, vector<double> &gh);
	void Run_UTSP_AS_Angran(string fn_in, string fn_out, int ngref, int nlref, int fit_flag, string fn_lref,
							vector<BezierElement> &bzmesh, vector<int> &IDBC, vector<double> &gh);

	void InputFromVTK_Quad_Angran(string fn);
	void InputFromVTK_Quad(string fn);
	void InputFromINP_Quad(string fn);
	void RescaleCoor();
	void InitialConnect_UT();
	void FindIEN_4();
	void SetDPatch_Irr();
	void ProjectVector(int type, int N, int eloc, double beta, vector<double> &pv); //type, 0 for non-negative, 1 for idempotent
	void ProjectMatrix(int type, int N, double beta, vector<vector<double>> &pmat);
	void VertexFaceTranMat(int N, vector<vector<double>> &bmat);
	void SetBezierMat_DPatch(int type, int N, double beta, vector<vector<double>> &bemat);
	void SetBezierMat_DPatch22(int type, int N, double beta, vector<vector<vector<double>>> &bemat);
	void ElementBasis_Irregular_DPatch(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);
	void ElementBasis_Irregular_DPatch22(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);
	void OutputBezierPoint22(string fn);
	void OutputBezierMesh(const vector<BezierElement> &bzmesh, string fn);
	void BezierExtract_DPatch(vector<BezierElement> &bzmesh);
	void BezierExtract_DPatch_AS(vector<BezierElement> &bzmesh);
	void BezierUnit_DPatch22(int eid, vector<BezierElement> &bzmesh);
	void BezierUnit_DPatch_Irr22(int eid, vector<BezierElement> &bzmesh);
	void BezierUnit_DPatch_Trs(int eid, vector<BezierElement> &bzmesh);
	void UpdateCP_Fitting(const vector<array<double, 3>> &cpnew);
	void ReadBezierMat(string fn, vector<BezierElement> &bzmesh);
	void InitializeProjectFacePoints(int proj_type, double proj_beta);
	void InitializeProjectFacePoints_Scale(int proj_type, double proj_beta);
	void AddProjectFacePoints(int proj_type, double proj_beta);
	void FindSplineC1_Irr(int eid);
	void FindSplineC1_C1Element(int eid);
	void FindSplineC1_Trs(int eid);
	void SetFaceBezierMat_DPatch(int type, int N, double beta, vector<vector<double>> &bemat);
	void SetFaceBezierMat_DPatch22(int type, int N, double beta, vector<vector<vector<double>>> &c1mat);
	void SetTrunMat();
	void SetTrunMat_Irr(int eid);
	void SetTrunMat_Irr22(int eid);
	void SetTrunMat_C1Element(int eid);
	void SetTrunMat_Trs(int eid);
	void CollectActives_DPatchAnalysis();
	void GeomMap_DPatch(int eid, double u, double v, array<double, 3> &pt);
	void ElementBasis_DPatch(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);
	void ElementBasis_DPatch_Irr22(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);
	void ElementBasis_DPatch_Trs(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);
	void FindGlobalRefineID(vector<int> &rfid, vector<int> &rftype);
	void FindC1Element(vector<int> &c1id);
	void Refine_C1Element(const vector<int> &c1id);
	void EnlargeOneRing();
	void BezierMeshUpdateIndex(const vector<BezierElement> &bzmesh0, vector<BezierElement> &bzmesh);
	void BC_Square(vector<int> &IDBC, vector<double> &gh);
	void BC_Hemisphere(vector<int> &IDBC, vector<double> &gh);
	void ReadBezier(string fn, vector<BezierElement> &bzmesh);
	void VisualizeCM_DPatch(string fn);

	void SetBezierMatIrr_NonUniform(int eid, vector<vector<double>> &bemat);
	//	void SetBezierMatTrs(int eid, vector<vector<double>>& bemat);
	void FindSplineC2_Scale();
	void FindSplineC2_Scale_Irr22(int eid);
	void FindSplineC2_Scale_Trs(int eid);
	void FindSplineC1_Scale();
	void FindSplineC1_Scale_Trs(int eid);
	void BezierUnit_DPatch_Irr22_Scale(int eid, vector<BezierElement> &bzmesh);
	void BezierUnit_DPatch_Trs_Scale(int eid, vector<BezierElement> &bzmesh);
	void BezierExtract_DPatch_AS_Scale(vector<BezierElement> &bzmesh);
	void VisualizeSurface_Scale(string fn, int ns = 5);
	void GeomMap_DPatch_Scale(int eid, double u, double v, array<double, 3> &pt);
	void ElementBasis_DPatch_Scale(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);
	void ElementBasis_DPatch_Irr22_Scale(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);
	void ElementBasis_DPatch_Trs_Scale(int eid, double u, double v, vector<double> &Nt, vector<array<double, 2>> &dNdt);
	void OutputBezierMeshAll(const vector<BezierElement> &bzmesh, string fn);
	void BezierExtract_DPatch_AS_Scale_Refine(vector<BezierElement> &bzmesh);

	void HemisphereLocalRefine();
	void HemisphereLocalRefine1();
	void HemisphereLocalRefine2();

	void HemisphereLocalRefine(const vector<int> &ref4, const vector<int> &ref2, const vector<int> &ref2_cn);
	void HemisphereLocalRefine1_6();
	void HemisphereLocalRefine2_6();
	void HemisphereLocalRefine(string fn, int nlev);
	void ReadRef4ID(string fn, vector<int> &rfid);
	void ReadRef2ID(string fn, vector<int> &rfid, vector<int> &rfid_cn);

	void WriteBezier_Shell(const vector<BezierElement> &bzmesh, string fn);
	void OutputBoundaryCPID_Hemisphere(string fn);
	void OutputBoundaryCPID_SimpleSquare(string fn_bc, string fn_bcin);
	void OutputBoundaryCPID_Crosspipe(string fn);
	void OutputBoundaryCPID_Cylinder(string fn);
	void OutputCorner(const vector<BezierElement> &bzmesh, const vector<array<double, 3>> &corn, double tol, string fn);
	void GetCornerInfo_Hemisphere(vector<array<double, 3>> &corn);
	void GetCornerInfo_Cylinder(vector<array<double, 3>> &corn);

	void SmoothPlanar(string fn_in, string fn_out, int nstep);

	void WriteBezier_LSDYNA(string fn, vector<BezierElement> &bzmesh);
	void WriteBezier_LSDYNA_C12(string fn, vector<BezierElement> &bzmesh);
	void ReadBEXT(string fn, vector<BezierElement> &bzmesh);
	void GeomMapBEXT(double u, double v, BezierElement &bzel, array<double, 3> &pt);
	void OutputBEXTGeom(string fn, vector<BezierElement> &bzmesh);

	void WriteBezier_Angran(vector<BezierElement> &bzmesh, string fn);
	void Extrude(string fn, const vector<BezierElement> &bzmesh, double t, int nel = 1, int dir = 0, int p = 3);		   //input surface is the midsurface
	void Extrude_GEM_layer(string fn, const vector<BezierElement> &bzmesh, double t, int nel = 1, int dir = 0, int p = 3); //Output layer information for GEM composite shell
	void Extrude_GEM_Shell(string fn, const vector<BezierElement> &bzmesh, double t, int nel = 1, int dir = 0, int p = 1); //Output layer information for GEM composite shell

	void ReadBEXT3D(string fn, vector<BezierElement> &bzmesh);
	void OutputBezierMesh_BEXT3D(string fn, int p, vector<BezierElement> &bzmesh);

	void DirichletLS_Scalar(vector<BezierElement> &bzmesh, vector<int> &IDBC, vector<double> &gh);

	//vector<BezierElement> bzmesh;//bezier mesh
	vector<int> paid;  //active control points
	vector<int> eaid;  //active elements
	vector<Vertex> cp; //control points
	int npta;

private:
	//vector<Vertex> cp;//control points
	vector<Element> tmesh; //elements in T-mesh
	vector<Edge> tmedge;
	//vector<vector<int>> phid;//id in different levels
	vector<int> cornerEID; //corner element ID

	vector<pair<int, double>> uanc;
	vector<pair<int, double>> vanc;
	vector<vector<int>> uedge;
	vector<vector<int>> vedge;

	unsigned int npt_old;
	unsigned int nel_old;
	//int npt_new;

	vector<array<double, 3>> bzcp;
	vector<array<double, 3>> bzcptmp;
	vector<double> wc1;
};

bool CompareAnchor(const pair<int, double> &anc1, const pair<int, double> &anc2);

#endif