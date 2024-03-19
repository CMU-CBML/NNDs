#include "kernel.h"
#include <ctime>

// void kernel::run()
// {
// 	int niter(20);
// 	unsigned int i;
// 	TruncatedTspline_3D tt3;
// 	//tt3.SetProblem("../io/hex_input/cube_dense");
// 	//tt3.SetProblem("../io/hex_input/cube4");
// 	tt3.SetProblem("../io/hex_input/cube_coarse_0");
// 	string fld("../io/hex_out_2/"), fn("test16_");
// 	vector<int> dof_list(niter,0);
// 	vector<double> err_list(niter,0.);
// 	for (int itr = 0; itr < niter; itr++)
// 	{
// 		vector<BezierElement3D> bzmesh;
// 		vector<int> IDBC;
// 		vector<double> gh, err;
// 		//cout << "interface...\n";
// 		tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
// 		//tt3.AnalysisInterface_Laplace(bzmesh, IDBC, gh);
// 		//cout << "interface done!\n";
// 		//getchar();

// 		//cout << "npt: " << IDBC.size() << "\n";
// 		//getchar();

// 		Laplace lap;
// 		lap.SetProblem(IDBC, gh);
// 		stringstream ss;
// 		ss << itr;
// 		//cout << "before simulation...\n";
// 		//getchar();
// 		lap.Run(bzmesh, fld+fn+ss.str(), err);
// 		//cout << "simulation done!\n";
// 		//getchar();

// 		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

// 		double errL2(0.);
// 		for (i = 0; i < err.size(); i++) errL2 += err[i];
// 		errL2 = sqrt(errL2);
// 		dof_list[itr] = IDBC.size();
// 		err_list[itr] = errL2;
// 		//cout << "DOF: " << IDBC.size() << "\n";
// 		//cout << "L2-norm error: " << errL2 << "\n";
// 		//getchar();

// 		if (itr < niter - 1)
// 		{
// 			cout << itr << " refining...\n";
// 			//distribute error
// 			vector<array<double, 2>> eh(bzmesh.size());
// 			for (i = 0; i < bzmesh.size(); i++)
// 			{
// 				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
// 			}
// 			vector<array<int, 2>> rfid, gst;
// 			//tt3.Identify_Poisson(eh, err, rfid, gst);
// 			tt3.Identify_Poisson_1(eh, err, rfid, gst);
// 			//tt3.Identify_Laplace(eh, err, rfid, gst);
// 			tt3.Refine(rfid, gst);
// 			cout << "Refining done\n";
// 			tt3.OutputGeom_All(fld+fn + ss.str()+"_geom");
// 			//cout << "Output Geom done!\n";
// 			//getchar();
// 		}
// 	}

// 	//output error
// 	output_err(fld + fn + "err", dof_list, err_list);
// }

// void kernel::run_complex()
// {
// 	int niter(2);
// 	unsigned int i;
// 	TruncatedTspline_3D tt3;

// 	tt3.SetProblem("../io/hex_input/statue");

// 	string fld("../io/complex2/");

// 	string fn("statue1_");

// 	double remove[3][2];
// 	tt3.GetRemoveRegion(remove);

// 	vector<int> dof_list(niter, 0);
// 	vector<double> err_list(niter, 0.);
// 	vector<array<int, 2>> ebc, pbc;
// 	vector<double> edisp, pdisp;
// 	tt3.SetInitialBC(ebc, edisp, pbc, pdisp);
// 	for (int itr = 0; itr < niter; itr++)
// 	{
// 		vector<BezierElement3D> bzmesh;
// 		vector<int> IDBC;
// 		vector<double> gh, err;
// 		tt3.SetBC(ebc, edisp, pbc, pdisp);
// 		//cout << pbc.size() << "\n";
// 		//cout << "interface...\n";
// 		//tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
// 		//tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
// 		tt3.AnalysisInterface_Laplace(pbc,pdisp,bzmesh, IDBC, gh);
// 		//cout << "interface done!\n";
// 		//getchar();

// 		//cout << "npt: " << IDBC.size() << "\n";
// 		//getchar();

// 		Laplace lap;
// 		lap.SetProblem(IDBC, gh);
// 		stringstream ss;
// 		ss << itr;
// 		//cout << "before simulation...\n";
// 		//getchar();
// 		lap.GetRemoveRegion(remove);
// 		lap.Run(bzmesh, fld + fn + ss.str(), err);
// 		//cout << "simulation done!\n";
// 		//getchar();

// 		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

// 		double errL2(0.);
// 		for (i = 0; i < err.size(); i++) errL2 += err[i];
// 		errL2 = sqrt(errL2);
// 		dof_list[itr] = IDBC.size();
// 		err_list[itr] = errL2;
// 		//cout << "DOF: " << IDBC.size() << "\n";
// 		//cout << "L2-norm error: " << errL2 << "\n";
// 		//getchar();

// 		output_err(fld + fn+ss.str() + "_err", dof_list, err_list);

// 		if (itr < niter - 1)
// 		{
// 			cout << itr << " refining...\n";
// 			//distribute error
// 			vector<array<double, 2>> eh(bzmesh.size());
// 			for (i = 0; i < bzmesh.size(); i++)
// 			{
// 				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
// 			}
// 			vector<array<int, 2>> rfid, gst;
// 			//tt3.Identify_Poisson_1(eh, err, rfid, gst);
// 			tt3.Identify_Laplace(eh, err, rfid, gst);
// 			//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
// 			//tt3.InputRefineID("../io/benchmark1/cube1_" + ss.str(), rfid, gst);
// 			tt3.Refine(rfid, gst);
// 			cout << "Refining done\n";
// 			//tt3.OutputGeom_All(fld + fn + ss.str() + "_geom");
// 			//cout << "Output Geom done!\n";
// 			//getchar();
// 		}
// 	}

// 	//output error
// 	//output_err(fld + fn + "err", dof_list, err_list);
// }

void kernel::OutputMesh(const vector<BezierElement3D>& bzmesh, string fn)
{
	int cn[8] = { 0, 3, 15, 12, 48, 51, 63, 60 };
	string fname = fn + "bzmesh.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nBezier mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << 8 * bzmesh.size() << " float\n";
		for (int i = 0; i<bzmesh.size(); i++)
		{
			for (int j = 0; j < 8; j++)
			{
				fout << bzmesh[i].pts[cn[j]][0] << " " << bzmesh[i].pts[cn[j]][1] << " " << bzmesh[i].pts[cn[j]][2] << "\n";
			}
		}
		fout << "\nCELLS " << bzmesh.size() << " " << 9 * bzmesh.size() << '\n';
		for (int i = 0; i<bzmesh.size(); i++)
		{
			fout << "8 " << 8 * i << " " << 8 * i + 1 << " " << 8 * i + 2 << " " << 8 * i + 3
				<< " " << 8 * i + 4 << " " << 8 * i + 5 << " " << 8 * i + 6 << " " << 8 * i + 7 << '\n';
		}
		fout << "\nCELL_TYPES " << bzmesh.size() << '\n';
		for (int i = 0; i<bzmesh.size(); i++)
		{
			fout << "12\n";
		}
		//fout << "POINT_DATA " << sdisp.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		//for (int i = 0; i<sdisp.size(); i++)
		//{
		//	fout << sdisp[i] << "\n";
		//}
		fout << "\nCELL_DATA " << bzmesh.size() << "\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		for (int i = 0; i < bzmesh.size(); i++)
		{
			fout << bzmesh[i].type << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
	string fname3(fn + "bzmeshinfo.txt");
	//ofstream fout;
	fout.open(fname3.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() << "\n";
		for (int i = 0; i<bzmesh.size(); i++)
		{
			for (int l = 0; l < bzmesh[i].IEN.size(); l++)
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
	fout.open(fname1.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() << "\n";
		for (int i = 0; i<bzmesh.size(); i++)
		{
			fout << i << " " << bzmesh[i].IEN.size() << " " << bzmesh[i].type << "\n";
			for (int l = 0; l < bzmesh[i].IEN.size(); l++)
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
			for (int j = 0; j < bzmesh[i].cmat.size(); j++)
			{
				for (int k = 0; k < bzmesh[i].cmat[j].size(); k++)
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
	fout.open(fname2.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() * 64 << "\n";
		for (int i = 0; i<bzmesh.size(); i++)
		{
			for (int j = 0; j < 64; j++)
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

void kernel::run_complex_fit()
{
	int niter(2);
	double thresh(0.25);//cube
	unsigned int i;
	double xy[3][2], nm[3], a(50.);
	TruncatedTspline_3D tt3;

	clock_t begin = clock();

	tt3.SetProblem("../io/hex_input/cube5");
	// tt3.SetProblem("./io/plate_input/input_CM");

	string fld("../io/kk_test/");

	string fn("base_");

	tt3.SetDomainRange(xy, nm, a);

	vector<int> dof_list(niter, 0);
	vector<double> err_list(niter, 0.);

	int itr;
	double errL2(1.e6);
	for (itr = 0; itr < niter; itr++)
	{
		vector<BezierElement3D> bzmesh;
		vector<int> IDBC;
		vector<double> gh, err;

		tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);

 		Laplace lap;
		lap.SetProblem(IDBC, gh);
		stringstream ss;
		ss << itr;
		lap.GetEqParameter(xy, nm, a);
		lap.Run(bzmesh, fld + fn + ss.str(), err);

		std::cout << "+++++++++++++++++++++" << std::endl;
		std::cout << err.size() << std::endl;
		std::cout << "+++++++++++++++++++++" << std::endl;

		errL2 = 0.;
		for (i = 0; i < err.size(); i++) errL2 += err[i];
		errL2 = sqrt(errL2);
		dof_list[itr] = IDBC.size();
		err_list[itr] = errL2;

		output_err(fld + fn + ss.str() + "_err", dof_list, err_list);

		// if (errL2 > thresh)
		if (itr > 0)
		{
			cout << itr << " refining...\n";
			//distribute error
			vector<array<double, 2>> eh(bzmesh.size());
			for (i = 0; i < bzmesh.size(); i++)
			{
				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
			}
			
			vector<array<int, 2>> rfid, gst;
			// tt3.Identify_Poisson_1(eh, err, rfid, gst);
			tt3.Identify_Laplace(eh, err, rfid, gst);			
			tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
			tt3.InputRefineID("../io/kk_test/base_" + ss.str(), rfid, gst);
			tt3.Refine(rfid, gst);
			cout << "Refining done\n";

			tt3.OutputGeom_All(fld + fn + ss.str() + "_geom");
			cout << "Output Geom done!\n";
		}
		OutputMesh(bzmesh, "../ioTHS3D/");

	}

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "\nElapsed time: " << elapsed_secs << "\n";

	//output error
	//output_err(fld + fn + "err", dof_list, err_list);
}

// void kernel::run_complex_glb()
// {
// 	int niter(1);
// 	unsigned int i;
// 	double xy[3][2], nm[3], a(50.);
// 	TruncatedTspline_3D tt3;
// 	//tt3.InitializeMesh("../io/hex_input/cube_coarse_0");
// 	tt3.InitializeMesh("../io/hex_input/rod");
// 	//tt3.InitializeMesh("../io/hex_input/hook");
// 	//tt3.InitializeMesh("../io/hex_input/cad1");
// 	//tt3.InitializeMesh("../io/hex_input/statue");
// 	//tt3.InitializeMesh("../io/hex_input/gear");

// 	//tt3.SetProblem("../io/hex_input/cube_dense");
// 	//tt3.SetProblem("../io/hex_input/cube3");
// 	//tt3.SetProblem("../io/hex_input/gear");

// 	string fld("../io/global3/");
// 	//string fld("../io/benchmark3/");

// 	//string fn("cube0_4");
// 	string fn("rod1_2");
// 	//string fn("hook1_");
// 	//string fn("base0_");
// 	//string fn("head1_");
// 	//string fn("gear1_");

// 	tt3.run_Global(2, fld + fn);
// 	cout << "Done global refinement!\n";
// 	//getchar();

// 	//double remove[3][2];
// 	//tt3.GetRemoveRegion(remove);
// 	tt3.SetDomainRange(xy, nm, a);

// 	//tt3.OutputRemoveCM(fld+fn,remove);
// 	//cout << "done\n";
// 	//getchar();

// 	vector<int> dof_list(niter, 0);
// 	vector<double> h_max(niter, 0.);
// 	vector<double> err_list(niter, 0.);
// 	for (int itr = 0; itr < niter; itr++)
// 	{
// 		vector<BezierElement3D> bzmesh;
// 		vector<int> IDBC;
// 		vector<double> gh, err;
// 		tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);

// 		Laplace lap;
// 		lap.SetProblem(IDBC, gh);
// 		stringstream ss;
// 		ss << itr;
// 		//lap.GetRemoveRegion(remove);
// 		lap.GetEqParameter(xy, nm, a);
// 		lap.Run(bzmesh, fld + fn, err);

// 		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());
// 		//cout << "Done visualize Bezier!\n";
// 		//getchar();

// 		double errL2(0.);
// 		for (i = 0; i < err.size(); i++) errL2 += err[i];
// 		errL2 = sqrt(errL2);
// 		dof_list[itr] = IDBC.size();
// 		h_max[itr] = tt3.MaxElementSize(bzmesh);
// 		err_list[itr] = errL2;
// 		//cout << "DOF: " << IDBC.size() << "\n";
// 		//cout << "L2-norm error: " << errL2 << "\n";
// 		//getchar();

// 		output_err(fld + fn + "_err", dof_list, err_list);
// 		output_err(fld + fn + "_err_h", h_max, err_list);
// 	}
// }

// void kernel::run_Bspline()
// {
// 	int niter(1);
// 	unsigned int i;
// 	string fld("../io/BSP/");
// 	//string fn("ZBSP3_");
// 	string fn("ZBSP3_");
// 	TruncatedTspline_3D tt3;
// 	vector<BezierElement3D> bzmesh;
// 	vector<int> IDBC;
// 	vector<double> gh, err;
// 	int nsp(4);
// 	int nex[3] = { nsp, nsp, nsp };
// 	//tt3.CreateBsplines(bzmesh, IDBC, gh);
// 	tt3.CreateBsplines(nex);
// 	int nref(3);
// 	for (int i = 0; i < nref; i++)
// 	{
// 		tt3.Bsplines_Refine();
// 	}
// 	tt3.Bspline_BezierExtract(bzmesh, IDBC, gh);

// 	//tt3.VisualizeBezier(bzmesh, fld + fn);
// 	//cout << "Done visualize Bezier!\n";
// 	//getchar();

// 	Laplace lap;
// 	lap.SetProblem(IDBC, gh);
// 	lap.Run(bzmesh, fld + fn, err);

// 	vector<int> dof_list(niter, 0);
// 	vector<double> h_max(niter, 0.);
// 	vector<double> err_list(niter, 0.);
// 	double errL2(0.);
// 	for (i = 0; i < err.size(); i++) errL2 += err[i];
// 	errL2 = sqrt(errL2);
// 	dof_list[0] = IDBC.size();
// 	h_max[0] = tt3.MaxElementSize(bzmesh);
// 	err_list[0] = errL2;

// 	output_err(fld + fn + "_err", dof_list, err_list);
// 	output_err(fld + fn + "_err_h", h_max, err_list);
// }

// void kernel::run_Bspline_fit()
// {
// 	int niter(1);
// 	unsigned int i;
// 	double xy[3][2], nm[3], a(50.);
// 	string fld("../io/BSP/");
// 	string fn("BSP2_fit4");
// 	TruncatedTspline_3D tt3;
// 	vector<BezierElement3D> bzmesh;
// 	vector<int> IDBC;
// 	vector<double> gh, err;
// 	int nsp(2);
// 	int nex[3] = { nsp, nsp, nsp };
// 	//tt3.CreateBsplines(bzmesh, IDBC, gh);
// 	tt3.CreateBsplines(nex);
// 	int nref(4);
// 	for (int i = 0; i < nref; i++)
// 	{
// 		tt3.Bsplines_Refine();
// 	}
// 	tt3.SetDomainRange_Bsplines(xy, nm, a);//for solution 7
// 	tt3.Bspline_BezierExtract_fit(bzmesh, IDBC, gh);
// 	tt3.FittingBC_Bsplines(bzmesh, IDBC, gh);

// 	//tt3.VisualizeBezier(bzmesh, fld + fn);
// 	//cout << "Done visualize Bezier!\n";
// 	//getchar();

// 	Laplace lap;
// 	lap.SetProblem(IDBC, gh);
// 	lap.GetEqParameter(xy, nm, a);//for solution 7
// 	lap.Run(bzmesh, fld + fn, err);

// 	vector<int> dof_list(niter, 0);
// 	vector<double> h_max(niter, 0.);
// 	vector<double> err_list(niter, 0.);
// 	double errL2(0.);
// 	for (i = 0; i < err.size(); i++) errL2 += err[i];
// 	errL2 = sqrt(errL2);
// 	dof_list[0] = IDBC.size();
// 	h_max[0] = tt3.MaxElementSize(bzmesh);
// 	err_list[0] = errL2;

// 	output_err(fld + fn + "_err", dof_list, err_list);
// 	output_err(fld + fn + "_err_h", h_max, err_list);
// }

// void kernel::run_benchmark()
// {
// 	int niter(4);
// 	unsigned int i;
// 	TruncatedTspline_3D tt3;
// 	//tt3.SetProblem("../io/hex_input/cube_dense");
// 	//tt3.SetProblem("../io/hex_input/cube9");
// 	//tt3.SetProblem("../io/hex_input/cube_coarse_0");
// 	tt3.SetProblem("../io/hex_input/cube_coarse_2");
// 	string fld("../io/benchmark4/");
// 	string fn("cube_THS1_");
// 	//string fn("cube_HS1_");

// 	//double remove[3][2];
// 	//tt3.GetRemoveRegion(remove);
// 	//tt3.OutputRemoveCM(fld+fn,remove);
// 	//cout << "done\n";
// 	//getchar();

// 	vector<int> dof_list(niter, 0);
// 	vector<double> err_list(niter, 0.);
// 	vector<array<int, 2>> ebc, pbc;
// 	vector<double> edisp, pdisp;
// 	//tt3.SetInitialBC(ebc, edisp, pbc, pdisp);
// 	for (int itr = 0; itr < niter; itr++)
// 	{
// 		vector<BezierElement3D> bzmesh;
// 		vector<int> IDBC;
// 		vector<double> gh, err;
// 		//tt3.SetBC(ebc, edisp, pbc, pdisp);
// 		//cout << pbc.size() << "\n";
// 		//cout << "interface...\n";
// 		tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
// 		//tt3.AnalysisInterface_Laplace(pbc,pdisp,bzmesh, IDBC, gh);
// 		//cout << "interface done!\n";
// 		//getchar();

// 		//cout << "npt: " << IDBC.size() << "\n";
// 		//getchar();

// 		Laplace lap;
// 		lap.SetProblem(IDBC, gh);
// 		stringstream ss;
// 		ss << itr;
// 		//cout << "before simulation...\n";
// 		//getchar();
// 		//lap.GetRemoveRegion(remove);
// 		lap.Run(bzmesh, fld + fn + ss.str(), err);
// 		//cout << "simulation done!\n";
// 		//getchar();

// 		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

// 		double errL2(0.);
// 		for (i = 0; i < err.size(); i++) errL2 += err[i];
// 		errL2 = sqrt(errL2);
// 		dof_list[itr] = IDBC.size();
// 		err_list[itr] = errL2;
// 		//cout << "DOF: " << IDBC.size() << "\n";
// 		//cout << "L2-norm error: " << errL2 << "\n";
// 		//getchar();

// 		output_err(fld + fn + ss.str() + "_err", dof_list, err_list);

// 		if (itr < niter - 1)
// 		{
// 			cout << itr << " refining...\n";
// 			//distribute error
// 			vector<array<double, 2>> eh(bzmesh.size());
// 			for (i = 0; i < bzmesh.size(); i++)
// 			{
// 				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
// 			}
// 			vector<array<int, 2>> rfid, gst;
// 			//tt3.Identify_Poisson_1(eh, err, rfid, gst);
// 			//tt3.Identify_Laplace(eh, err, rfid, gst);
// 			//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
// 			tt3.InputRefineID("../io/benchmark3/cube_THS1_" + ss.str(), rfid, gst);
// 			tt3.Refine(rfid, gst);
// 			cout << "Refining done\n";
// 			//tt3.OutputGeom_All(fld + fn + ss.str() + "_geom");
// 			//cout << "Output Geom done!\n";
// 			//getchar();
// 		}
// 	}

// 	//output error
// 	//output_err(fld + fn + "err", dof_list, err_list);
// }

// void kernel::run_leastsquare()
// {
// 	int niter(2);
// 	unsigned int i;
// 	TruncatedTspline_3D tt3;
// 	//tt3.SetProblem("../io/hex_input/cube_dense");
// 	//tt3.SetProblem("../io/hex_input/cube4");
// 	tt3.SetProblem("../io/hex_input/cube_coarse_2");
// 	string fld("../io/least_square1/");
// 	//string fn("cube_THS2_2");
// 	string fn("cube_HS2_2");

// 	//tt3.GlobalRefine(2);

// 	for (int itr = 0; itr < niter; itr++)
// 	{
// 		cout << itr << " refining...\n";
// 		vector<array<int, 2>> rfid, gst;
// 		//tt3.Identify_Poisson_1(eh, err, rfid, gst);
// 		//tt3.Identify_Laplace(eh, err, rfid, gst);
// 		//tt3.Identify_LeastSquare(rfid, gst);
// 		tt3.Identify_LeastSquare_Line(rfid, gst);
// 		//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
// 		//tt3.InputRefineID("../io/benchmark3/cube_THS1_" + ss.str(), rfid, gst);
// 		tt3.Refine(rfid, gst);
// 		stringstream ss;
// 		ss << itr;
// 		cout << "Refining done\n";
// 		//tt3.OutputGeom_All(fld + fn + ss.str());
// 		//cout << "Output Geom done!\n";
// 		//getchar();
// 	}

// 	//cout << "Refining done\n";
// 	//tt3.OutputGeom_All(fld + fn + "_geom");
// 	//cout << "geom done!\n";
// 	//getchar();

// 	vector<BezierElement3D> bzmesh;
// 	vector<int> IDBC;
// 	vector<double> gh, err;
// 	//tt3.AnalysisInterface_LeastSquare(bzmesh, IDBC, gh);
// 	tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
// 	LeastSquare lap;
// 	lap.SetProblem(IDBC, gh);
// 	lap.Run(bzmesh, fld + fn, err);

// 	//output error
// 	//output_err(fld + fn + "err", dof_list, err_list);
// }

// void kernel::run_leastsquare_1()
// {
// 	int niter(4);
// 	unsigned int i;
// 	TruncatedTspline_3D tt3;
// 	//tt3.SetProblem("../io/hex_input/cube_dense");
// 	//tt3.SetProblem("../io/hex_input/cube9");
// 	//tt3.SetProblem("../io/hex_input/cube_coarse_0");
// 	tt3.SetProblem("../io/hex_input/cube_coarse_2");
// 	string fld("../io/least_square1/");
// 	string fn("cube_THS1_");
// 	//string fn("cube_HS1_");

// 	//double remove[3][2];
// 	//tt3.GetRemoveRegion(remove);
// 	//tt3.OutputRemoveCM(fld+fn,remove);
// 	//cout << "done\n";
// 	//getchar();

// 	vector<int> dof_list(niter, 0);
// 	vector<double> err_list(niter, 0.);
// 	vector<array<int, 2>> ebc, pbc;
// 	vector<double> edisp, pdisp;
// 	//tt3.SetInitialBC(ebc, edisp, pbc, pdisp);
// 	for (int itr = 0; itr < niter; itr++)
// 	{
// 		vector<BezierElement3D> bzmesh;
// 		vector<int> IDBC;
// 		vector<double> gh, err;
// 		//tt3.SetBC(ebc, edisp, pbc, pdisp);
// 		//cout << pbc.size() << "\n";
// 		//cout << "interface...\n";
// 		//tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
// 		//tt3.AnalysisInterface_Laplace(pbc,pdisp,bzmesh, IDBC, gh);
// 		tt3.AnalysisInterface_LeastSquare(bzmesh, IDBC, gh);
// 		//cout << "interface done!\n";
// 		//getchar();

// 		//cout << "npt: " << IDBC.size() << "\n";
// 		//getchar();

// 		Laplace lap;
// 		lap.SetProblem(IDBC, gh);
// 		stringstream ss;
// 		ss << itr;
// 		//cout << "before simulation...\n";
// 		//getchar();
// 		//lap.GetRemoveRegion(remove);
// 		lap.Run(bzmesh, fld + fn + ss.str(), err);
// 		//cout << "simulation done!\n";
// 		//getchar();

// 		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

// 		double errL2(0.);
// 		for (i = 0; i < err.size(); i++) errL2 += err[i];
// 		errL2 = sqrt(errL2);
// 		dof_list[itr] = IDBC.size();
// 		err_list[itr] = errL2;
// 		//cout << "DOF: " << IDBC.size() << "\n";
// 		//cout << "L2-norm error: " << errL2 << "\n";
// 		//getchar();

// 		output_err(fld + fn + ss.str() + "_err", dof_list, err_list);

// 		if (itr < niter - 1)
// 		{
// 			cout << itr << " refining...\n";
// 			//distribute error
// 			vector<array<double, 2>> eh(bzmesh.size());
// 			for (i = 0; i < bzmesh.size(); i++)
// 			{
// 				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
// 			}
// 			vector<array<int, 2>> rfid, gst;
// 			//tt3.Identify_Poisson_1(eh, err, rfid, gst);
// 			//tt3.Identify_Laplace(eh, err, rfid, gst);
// 			//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
// 			tt3.InputRefineID("../io/benchmark3/cube_THS1_" + ss.str(), rfid, gst);
// 			tt3.Refine(rfid, gst);
// 			cout << "Refining done\n";
// 			//tt3.OutputGeom_All(fld + fn + ss.str() + "_geom");
// 			//cout << "Output Geom done!\n";
// 			//getchar();
// 		}
// 	}

// 	//output error
// 	//output_err(fld + fn + "err", dof_list, err_list);
// }

// void kernel::run_pipeline()
// {
// 	int niter(5);
// 	unsigned int i;
// 	TruncatedTspline_3D tt3;
// 	tt3.SetProblem("../io/pipeline/neuron");
// 	string fld("../io/pipeline/");
// 	string fn("neuron1");

// 	//double remove[3][2];
// 	//tt3.GetRemoveRegion(remove);
// 	//tt3.OutputRemoveCM(fld+fn,remove);
// 	//cout << "done\n";
// 	//getchar();

// 	vector<BezierElement3D> bzmesh;
// 	tt3.PipelineBezierExtract(bzmesh);

// 	Laplace lap;
// 	lap.PipelineTmp(bzmesh,fld+fn);
// }

void kernel::output_err(string fn, const vector<int>& dof, const vector<double>& err)
{
	string fname = fn + ".txt";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		for (i = 0; i < dof.size(); i++)
		{
			fout << dof[i] << " " << err[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void kernel::output_err(string fn, const vector<double>& dof, const vector<double>& err)
{
	string fname = fn + ".txt";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		for (i = 0; i < dof.size(); i++)
		{
			fout << dof[i] << " " << err[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}