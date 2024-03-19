#ifndef _MATLAB_SOLVER_H_
#define _MATLAB_SOLVER_H_
 
#include <string.h>
//#include "engine.h"
//#include "matrix.h"
#include <sstream>
#include <vector>
#include <iostream>

using namespace std;

class MatlabSolver
{
public:
	MatlabSolver()
	{ 
		
		ix = NULL;
		iy = NULL;
		dxy = NULL;
		bxy = NULL;
		rxy = NULL;
	}
	~MatlabSolver()
	{ 
		
	}

	void Initilize(vector<double> &idx, vector<double> &idy, vector<double> &dataxy, vector<double> &bdata, int lenx) {
		ep = engOpen("");
		int len = idx.size(), lenb = bdata.size();
		int sz = 8 * len, sz1 = 8 * lenb;
		ix = mxCreateNumericMatrix(1, len, mxDOUBLE_CLASS, mxREAL);
		iy = mxCreateNumericMatrix(1, len, mxDOUBLE_CLASS, mxREAL);
		dxy = mxCreateNumericMatrix(1, len, mxDOUBLE_CLASS, mxREAL);
		bxy = mxCreateNumericMatrix(1, lenb, mxDOUBLE_CLASS, mxREAL);
		rxy = mxCreateNumericMatrix(1, lenx, mxDOUBLE_CLASS, mxREAL);

		memcpy((void *)mxGetPr(ix), (void *)idx.data(), sz);
		memcpy((void *)mxGetPr(iy), (void *)idy.data(), sz);
		memcpy((void *)mxGetPr(dxy), (void *)dataxy.data(), sz);
		memcpy((void *)mxGetPr(bxy), (void *)bdata.data(), sz1);

		lene = lenb;
		lenr = lenx;
		
		return;
	}
	void Solve(double *rlt) {

		engPutVariable(ep, "ix", ix);
		engPutVariable(ep, "iy", iy);
		engPutVariable(ep, "dxy", dxy);
		engPutVariable(ep, "b", bxy);
		ss << "s = sparse(ix,iy,dxy,";
		ss << lene << "," << lenr << ");";
		engEvalString(ep, ss.str().c_str());
		//engEvalString(ep, "s = full(s);");
		engEvalString(ep, "b = b';");
		engEvalString(ep, "[sz1,sz2]=size(b);");
		engEvalString(ep, "lb = zeros(sz1,1);");
		engEvalString(ep, "ub = 0.5*ones(sz1,1);");
		engEvalString(ep, "x = lsqlin(s, b, [], [], [], [], lb, ub);");
		rxy = engGetVariable(ep, "x");
		memcpy(rlt, (void *)mxGetPr(rxy), lenr * 8);

		mxDestroyArray(ix);
		mxDestroyArray(iy);
		mxDestroyArray(dxy);
		mxDestroyArray(bxy);
		mxDestroyArray(rxy);
		
		engClose(ep);
		return;
	}
	void Initilize1(vector<double> &idx, vector<double> &idy, vector<double> &dataxy, vector<double> &bdata, 
					vector<double> &fidx, vector<double> &fidy, vector<double> &fdataxy, vector<double> &fbdata,  int flenx) {
		ep = engOpen("");
		int len = idx.size(), lenb = bdata.size();
		int sz = 8 * len, sz1 = 8 * lenb;
		ix = mxCreateNumericMatrix(1, len, mxDOUBLE_CLASS, mxREAL);
		iy = mxCreateNumericMatrix(1, len, mxDOUBLE_CLASS, mxREAL);
		dxy = mxCreateNumericMatrix(1, len, mxDOUBLE_CLASS, mxREAL);
		bxy = mxCreateNumericMatrix(1, lenb, mxDOUBLE_CLASS, mxREAL);

		memcpy((void *)mxGetPr(ix), (void *)idx.data(), sz);
		memcpy((void *)mxGetPr(iy), (void *)idy.data(), sz);
		memcpy((void *)mxGetPr(dxy), (void *)dataxy.data(), sz);
		memcpy((void *)mxGetPr(bxy), (void *)bdata.data(), sz1);
		 
		lene = lenb;
		lenr = bdata.size();

		int flen = fidx.size(), flenb = fbdata.size();
		int fsz = 8 * flen, fsz1 = 8 * flenb;
		fix = mxCreateNumericMatrix(1, flen, mxDOUBLE_CLASS, mxREAL);
		fiy = mxCreateNumericMatrix(1, flen, mxDOUBLE_CLASS, mxREAL);
		fdxy = mxCreateNumericMatrix(1, flen, mxDOUBLE_CLASS, mxREAL);
		fbxy = mxCreateNumericMatrix(1, flenb, mxDOUBLE_CLASS, mxREAL);
		frxy = mxCreateNumericMatrix(1, flenx, mxDOUBLE_CLASS, mxREAL);

		memcpy((void *)mxGetPr(fix), (void *)fidx.data(), fsz);
		memcpy((void *)mxGetPr(fiy), (void *)fidy.data(), fsz);
		memcpy((void *)mxGetPr(fdxy), (void *)fdataxy.data(), fsz);
		memcpy((void *)mxGetPr(fbxy), (void *)fbdata.data(), fsz1);

		flene = flenb;
		flenr = flenx;

		return;
	}
	void Solve1(double *rlt ) {

		engPutVariable(ep, "ix", ix);
		engPutVariable(ep, "iy", iy);
		engPutVariable(ep, "dxy", dxy);
		engPutVariable(ep, "b", bxy);
		ss << "s = sparse(ix,iy,dxy,";
		ss << lene << "," << lenr << ");";
		engEvalString(ep, ss.str().c_str());
		engEvalString(ep, "b = b';");

		engPutVariable(ep, "fix", fix);
		engPutVariable(ep, "fiy", fiy);
		engPutVariable(ep, "fdxy", fdxy);
		engPutVariable(ep, "fb", fbxy);

		ss.clear();
		ss.str("");
		ss << "fs = sparse(fix,fiy,fdxy,";
		ss << flene << "," << flenr << ");";
		engEvalString(ep, ss.str().c_str());
		//engEvalString(ep, "fs = full(fs);");
		engEvalString(ep, "fb = fb';");
		engEvalString(ep, "[sz1,sz2]=size(fb);");
		engEvalString(ep, "lb = zeros(sz1,1);");
		engEvalString(ep, "ub = 0.5*ones(sz1,1);");
		engEvalString(ep, "ix = ix';");
		engEvalString(ep, "options = optimoptions('lsqlin','Algorithm','active-set');");
		//engEvalString(ep, "x = lsqlin(fs, fb, [], [], s, b, lb, ub);");
		//engEvalString(ep, "x = lsqlin(fs, fb, [], [], [],[], lb, ub);");
		engEvalString(ep, "x = lsqlin(fs, fb, [], [], [], [],lb, ub, ix, options);");
		//engEvalString(ep, "x = lsqnonneg(fs, fb);");
		frxy = engGetVariable(ep, "x");
		memcpy(rlt, (void *)mxGetPr(frxy), flenr * 8);


		mxDestroyArray(ix);
		mxDestroyArray(iy);
		mxDestroyArray(dxy);
		mxDestroyArray(bxy);

		mxDestroyArray(fix);
		mxDestroyArray(fiy);
		mxDestroyArray(fdxy);
		mxDestroyArray(fbxy);
		mxDestroyArray(frxy); 

		engClose(ep);
		return;
	}
	 
	void Initilize2(vector<double> &idx, vector<double> &idy, vector<double> &dataxy, vector<double> &bdata, 
		vector<double> &fidx, vector<double> &fidy, vector<double> &fdataxy, vector<double> &fbdata, vector<double> &pre_value) {
		ep = engOpen("");
		engEvalString(ep, "clear;");
		int len = idx.size(), lenb = bdata.size();
		//for(int i=0; i<len; i++)
		//{
		//	cout<<idx[i]<<" "<<idy[i]<<" "<<dataxy[i]<<"\n";
		//	getchar();
		//}
		int sz = 8 * len, sz1 = 8 * lenb;
		ix = mxCreateNumericMatrix(1, len, mxDOUBLE_CLASS, mxREAL);
		iy = mxCreateNumericMatrix(1, len, mxDOUBLE_CLASS, mxREAL);
		dxy = mxCreateNumericMatrix(1, len, mxDOUBLE_CLASS, mxREAL);
		bxy = mxCreateNumericMatrix(1, lenb, mxDOUBLE_CLASS, mxREAL);
		rxy = mxCreateNumericMatrix(1, pre_value.size(), mxDOUBLE_CLASS, mxREAL);
		

		memcpy((void *)mxGetPr(ix), (void *)idx.data(), sz);
		memcpy((void *)mxGetPr(iy), (void *)idy.data(), sz);
		memcpy((void *)mxGetPr(dxy), (void *)dataxy.data(), sz);
		memcpy((void *)mxGetPr(bxy), (void *)bdata.data(), sz1);

		lene = lenb;
		lenr = pre_value.size();
		flenr = pre_value.size();
		flene = fbdata.size();

		int flen = fidx.size(), flenb = fbdata.size();
		int fsz = 8 * flen, fsz1 = 8 * flenb;
		fix = mxCreateNumericMatrix(1, flen, mxDOUBLE_CLASS, mxREAL);
		fiy = mxCreateNumericMatrix(1, flen, mxDOUBLE_CLASS, mxREAL);
		fdxy = mxCreateNumericMatrix(1, flen, mxDOUBLE_CLASS, mxREAL);
		fbxy = mxCreateNumericMatrix(1, flenb, mxDOUBLE_CLASS, mxREAL);
		frxy = mxCreateNumericMatrix(1, flenr, mxDOUBLE_CLASS, mxREAL);
		pre_xy = mxCreateNumericMatrix(1, flenr, mxDOUBLE_CLASS, mxREAL);

		memcpy((void *)mxGetPr(fix), (void *)fidx.data(), fsz);
		memcpy((void *)mxGetPr(fiy), (void *)fidy.data(), fsz);
		memcpy((void *)mxGetPr(fdxy), (void *)fdataxy.data(), fsz);
		memcpy((void *)mxGetPr(fbxy), (void *)fbdata.data(), fsz1);
		memcpy((void *)mxGetPr(pre_xy), (void *)pre_value.data(), flenr * 8);

		
		

		return;
	}

	void Initilize3(vector<double> &idx, vector<double> &idy, vector<double> &dataxy, vector<double> &bdata,
		vector<double> &idx1, vector<double> &idy1, vector<double> &dataxy1, vector<double> &bdata1,
		vector<double> &fidx, vector<double> &fidy, vector<double> &fdataxy, vector<double> &fbdata, vector<double> &pre_value) {
		ep = engOpen("");
		engEvalString(ep, "clear;");
		int len = idx.size(), lenb = bdata.size();
		int sz = 8 * len, sz1 = 8 * lenb;
		ix = mxCreateNumericMatrix(1, len, mxDOUBLE_CLASS, mxREAL);
		iy = mxCreateNumericMatrix(1, len, mxDOUBLE_CLASS, mxREAL);
		dxy = mxCreateNumericMatrix(1, len, mxDOUBLE_CLASS, mxREAL);
		bxy = mxCreateNumericMatrix(1, lenb, mxDOUBLE_CLASS, mxREAL);
		rxy = mxCreateNumericMatrix(1, pre_value.size(), mxDOUBLE_CLASS, mxREAL);

		memcpy((void *)mxGetPr(ix), (void *)idx.data(), sz);
		memcpy((void *)mxGetPr(iy), (void *)idy.data(), sz);
		memcpy((void *)mxGetPr(dxy), (void *)dataxy.data(), sz);
		memcpy((void *)mxGetPr(bxy), (void *)bdata.data(), sz1);


		int len1 = idx1.size(), lenb1 = bdata1.size();
		int sz10 = 8 * len1, sz11 = 8 * lenb1;
		ix1 = mxCreateNumericMatrix(1, len1, mxDOUBLE_CLASS, mxREAL);
		iy1 = mxCreateNumericMatrix(1, len1, mxDOUBLE_CLASS, mxREAL);
		dxy1 = mxCreateNumericMatrix(1, len1, mxDOUBLE_CLASS, mxREAL);
		bxy1 = mxCreateNumericMatrix(1, lenb1, mxDOUBLE_CLASS, mxREAL);

		memcpy((void *)mxGetPr(ix1), (void *)idx1.data(), sz10);
		memcpy((void *)mxGetPr(iy1), (void *)idy1.data(), sz10);
		memcpy((void *)mxGetPr(dxy1), (void *)dataxy1.data(), sz10);
		memcpy((void *)mxGetPr(bxy1), (void *)bdata1.data(), sz11);
		

		lene1 = lenb1;
		lenr1 = pre_value.size();

		lene = lenb;
		lenr = pre_value.size();
		flenr = pre_value.size();
		flene = fbdata.size();

		int flen = fidx.size(), flenb = fbdata.size();
		int fsz = 8 * flen, fsz1 = 8 * flenb;
		fix = mxCreateNumericMatrix(1, flen, mxDOUBLE_CLASS, mxREAL);
		fiy = mxCreateNumericMatrix(1, flen, mxDOUBLE_CLASS, mxREAL);
		fdxy = mxCreateNumericMatrix(1, flen, mxDOUBLE_CLASS, mxREAL);
		fbxy = mxCreateNumericMatrix(1, flenb, mxDOUBLE_CLASS, mxREAL);
		frxy = mxCreateNumericMatrix(1, flenr, mxDOUBLE_CLASS, mxREAL);
		pre_xy = mxCreateNumericMatrix(1, flenr, mxDOUBLE_CLASS, mxREAL);

		memcpy((void *)mxGetPr(fix), (void *)fidx.data(), fsz);
		memcpy((void *)mxGetPr(fiy), (void *)fidy.data(), fsz);
		memcpy((void *)mxGetPr(fdxy), (void *)fdataxy.data(), fsz);
		memcpy((void *)mxGetPr(fbxy), (void *)fbdata.data(), fsz1);
		memcpy((void *)mxGetPr(pre_xy), (void *)pre_value.data(), flenr * 8);




		return;
	}

	void Solve2(double *rlt) {
		 
		engEvalString(ep, "clear;");
		engPutVariable(ep, "ix", ix);
		engPutVariable(ep, "iy", iy);
		engPutVariable(ep, "dxy", dxy);
		engPutVariable(ep, "b", bxy);
		ss.clear();
		ss.str("");
		ss << "s = sparse(ix,iy,dxy,";
		ss << lene << "," << lenr << ");";
		engEvalString(ep, ss.str().c_str());
		engEvalString(ep, "b = b';");

		engPutVariable(ep, "fix", fix);
		engPutVariable(ep, "fiy", fiy);
		engPutVariable(ep, "fdxy", fdxy);
		engPutVariable(ep, "fb", fbxy);
		engPutVariable(ep, "x0", pre_xy);

		ss.clear();
		ss.str("");
		ss << "fs = sparse(fix,fiy,fdxy,";
		ss << flene << "," << flenr << ");";
		engEvalString(ep, ss.str().c_str());
		//engEvalString(ep, "fs = full(fs);");
		engEvalString(ep, "fb = fb';");
		engEvalString(ep, "[sz1,sz2]=size(fb);");
		engEvalString(ep, "lb = zeros(sz1,1);");
		engEvalString(ep, "ub = 0.5*ones(sz1,1);");
		engEvalString(ep, "x0 = x0';");
		engEvalString(ep, "options = optimoptions('lsqlin','Algorithm','active-set');");
		//engEvalString(ep, "x = lsqlin(fs, fb, [], [], s, b, lb, ub);");
		//engEvalString(ep, "x = lsqlin(fs, fb, [], [], [],[], lb, ub);");
		engEvalString(ep, "x = lsqlin(fs, fb, [], [], s, b,lb, ub, x0, options);");
		//engEvalString(ep, "x = lsqnonneg(fs, fb);");
		frxy = engGetVariable(ep, "x");
		memcpy(rlt, (void *)mxGetPr(frxy), flenr * 8);


		mxDestroyArray(ix);
		mxDestroyArray(iy);
		mxDestroyArray(dxy);
		mxDestroyArray(bxy);
		mxDestroyArray(rxy);

		mxDestroyArray(fix);
		mxDestroyArray(fiy);
		mxDestroyArray(fdxy);
		mxDestroyArray(fbxy);
		mxDestroyArray(frxy);
		mxDestroyArray(pre_xy);

		engClose(ep);
		return;
	}
	void Solve3(double *rlt) {

		engEvalString(ep, "clear;");
		engPutVariable(ep, "ix", ix);
		engPutVariable(ep, "iy", iy);
		engPutVariable(ep, "dxy", dxy);
		engPutVariable(ep, "b", bxy);
		ss.clear();
		ss.str("");
		ss << "s = sparse(ix,iy,dxy,";
		ss << lene << "," << lenr << ");";
		engEvalString(ep, ss.str().c_str());
		engEvalString(ep, "b = b';");

		engPutVariable(ep, "fix", fix);
		engPutVariable(ep, "fiy", fiy);
		engPutVariable(ep, "fdxy", fdxy);
		engPutVariable(ep, "fb", fbxy);
		engPutVariable(ep, "x0", pre_xy);

		ss.clear();
		ss.str("");
		ss << "fs = sparse(fix,fiy,fdxy,";
		ss << flene << "," << flenr << ");";
		engEvalString(ep, ss.str().c_str());
		//engEvalString(ep, "fs = full(fs);");
		engEvalString(ep, "fb = fb';");
		engEvalString(ep, "x0 = x0';");
		//engEvalString(ep, "[sz1,sz2]=size(fb);");
		engEvalString(ep, "[sz1,sz2]=size(x0);");
		engEvalString(ep, "lb = zeros(sz1,1);");
		engEvalString(ep, "ub = 0.5*ones(sz1,1);");
		//engEvalString(ep, "x0 = x0';");
		engEvalString(ep, "options = optimoptions('lsqlin','Algorithm','active-set');");
		//engEvalString(ep, "x = lsqlin(fs, fb, [], [], s, b, lb, ub);");
		//engEvalString(ep, "x = lsqlin(fs, fb, [], [], [],[], lb, ub);");
		cout<<"manually solving...\n";
		//getchar();
		//engEvalString(ep, "x = lsqlin(fs, fb, [], [], s, b,lb, ub, x0, options);");
		engEvalString(ep, "x = lsqlin(fs, fb, [], [], s, b);");
		//getchar();
		//engEvalString(ep, "x = lsqlin(s, b);");
		//engEvalString(ep, "x = lsqlin(s, b);");
		//engEvalString(ep, "x = lsqlin(fs, fb, [], [], s, b,lb, ub, x0);");
		//getchar();
		//engEvalString(ep, "x = lsqlin(fs, fb, [], [], s, b, lb, ub, x0);");
		//engEvalString(ep, "x = lsqlin(fs, fb, [], [], [], [],[], [], x0);");
		//engEvalString(ep, "x = lsqnonneg(fs, fb);");
		getchar();
		frxy = engGetVariable(ep, "x");
		memcpy(rlt, (void *)mxGetPr(frxy), flenr * 8);


		mxDestroyArray(ix);
		mxDestroyArray(iy);
		mxDestroyArray(dxy);
		mxDestroyArray(bxy);
		mxDestroyArray(rxy);

		mxDestroyArray(fix);
		mxDestroyArray(fiy);
		mxDestroyArray(fdxy);
		mxDestroyArray(fbxy);
		mxDestroyArray(frxy);
		mxDestroyArray(pre_xy);

		engClose(ep);
		return;
	}
	void Solve31(double *rlt) {

		engEvalString(ep, "clear;");
		engPutVariable(ep, "ix", ix);
		engPutVariable(ep, "iy", iy);
		engPutVariable(ep, "dxy", dxy);
		engPutVariable(ep, "b", bxy);
		ss.clear();
		ss.str("");
		ss << "s = sparse(ix,iy,dxy,";
		ss << lene << "," << lenr << ");";
		engEvalString(ep, ss.str().c_str());
		engEvalString(ep, "b = b';");
		engPutVariable(ep, "x0", pre_xy);

		engEvalString(ep, "[sz1,sz2]=size(b);");
		engEvalString(ep, "lb = zeros(sz1,1);");
		engEvalString(ep, "ub = 0.4*ones(sz1,1);");
		engEvalString(ep, "x0 = x0';");
		engEvalString(ep, "options = optimoptions('lsqlin','Algorithm','active-set');");
		//engEvalString(ep, "x = lsqlin(fs, fb, [], [], s, b, lb, ub);");
		//engEvalString(ep, "x = lsqlin(fs, fb, [], [], [],[], lb, ub);");
		engEvalString(ep, "x = lsqlin(s, b, [], [], [], [],lb, ub, x0, options);");
		//engEvalString(ep, "x = lsqnonneg(fs, fb);");
		frxy = engGetVariable(ep, "x");
		memcpy(rlt, (void *)mxGetPr(frxy), flenr * 8);


		mxDestroyArray(ix);
		mxDestroyArray(iy);
		mxDestroyArray(dxy);
		mxDestroyArray(bxy);
		mxDestroyArray(rxy);

		mxDestroyArray(fix);
		mxDestroyArray(fiy);
		mxDestroyArray(fdxy);
		mxDestroyArray(fbxy);
		mxDestroyArray(frxy);
		mxDestroyArray(pre_xy);

		engClose(ep);
		return;
	}
	void Solve4(double *rlt) {
		 
		engEvalString(ep, "clear;");


		engEvalString(ep, "clear;");
		engPutVariable(ep, "ix1", ix1);
		engPutVariable(ep, "iy1", iy1);
		engPutVariable(ep, "dxy1", dxy1);
		engPutVariable(ep, "b1", bxy1);
		ss.clear();
		ss.str("");
		ss << "s1 = sparse(ix1,iy1,dxy1,";
		ss << lene1 << "," << lenr << ");";
		engEvalString(ep, ss.str().c_str());
		engEvalString(ep, "b1 = b1';");


		engPutVariable(ep, "ix", ix);
		engPutVariable(ep, "iy", iy);
		engPutVariable(ep, "dxy", dxy);
		engPutVariable(ep, "b", bxy);
		ss.clear();
		ss.str("");
		ss << "s = sparse(ix,iy,dxy,";
		ss << lene << "," << lenr << ");";
		engEvalString(ep, ss.str().c_str());
		engEvalString(ep, "b = b';");

		engPutVariable(ep, "fix", fix);
		engPutVariable(ep, "fiy", fiy);
		engPutVariable(ep, "fdxy", fdxy);
		engPutVariable(ep, "fb", fbxy);
		engPutVariable(ep, "x0", pre_xy);

		ss.clear();
		ss.str("");
		ss << "fs = sparse(fix,fiy,fdxy,";
		ss << flene << "," << flenr << ");";
		engEvalString(ep, ss.str().c_str());
		//engEvalString(ep, "fs = full(fs);");
		engEvalString(ep, "fb = fb';");
		engEvalString(ep, "[sz1,sz2]=size(fb);");
		engEvalString(ep, "lb = zeros(sz1,1);");
		engEvalString(ep, "ub = 0.4*ones(sz1,1);");
		engEvalString(ep, "x0 = x0';");
		engEvalString(ep, "options = optimoptions('lsqlin','Algorithm','active-set');");
		//engEvalString(ep, "x = lsqlin(fs, fb, [], [], s, b, lb, ub);");
		//engEvalString(ep, "x = lsqlin(fs, fb, s1, b1, s, b, lb, ub);");
		engEvalString(ep, "x = lsqlin(fs, fb, s1, b1, s, b, lb, ub, x0, options);");
		//engEvalString(ep, "x = lsqnonneg(fs, fb);");
		frxy = engGetVariable(ep, "x");
		memcpy(rlt, (void *)mxGetPr(frxy), flenr * 8);


		mxDestroyArray(ix);
		mxDestroyArray(iy);
		mxDestroyArray(dxy);
		mxDestroyArray(bxy);
		mxDestroyArray(rxy);

		mxDestroyArray(fix);
		mxDestroyArray(fiy);
		mxDestroyArray(fdxy);
		mxDestroyArray(fbxy);
		mxDestroyArray(frxy);
		mxDestroyArray(pre_xy);

		engClose(ep);
		return;
	}

	void Reset() {
		 
		ix1 = NULL;
		iy1 = NULL;
		dxy1 = NULL;
		bxy1 = NULL;


		ix = NULL;
		iy = NULL;
		dxy = NULL;
		bxy = NULL;
		rxy = NULL;

		fix = NULL;
		fiy = NULL;
		fdxy = NULL;
		fbxy = NULL;
		frxy = NULL;
		pre_xy = NULL;

	}
protected:
private:

	Engine *ep;

	mxArray *ix;
	mxArray *iy;
	mxArray *dxy;
	mxArray *bxy;

	mxArray *ix1;
	mxArray *iy1;
	mxArray *dxy1;
	mxArray *bxy1;

	mxArray *fix;
	mxArray *fiy;
	mxArray *fdxy;
	mxArray *fbxy;
	mxArray *pre_xy;

	mxArray *rxy;
	mxArray *frxy;
	std::stringstream ss;
	int  lenr1;
	int  lene1;

	int  lenr;
	int  lene;

	int  flenr;
	int  flene;
};

#endif