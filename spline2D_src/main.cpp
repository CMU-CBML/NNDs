//#include "HB-2D.h"
//#include "SingularEval.h"
//#include "BSplineBasis.h"
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include "KnotInsertion.h"
#include "Truncated_Tspline.h"
#include "Laplace.h"
#include "LeastSquare.h"

#include "mesh.h"
//#include <sstream>
#include <iomanip>
#include "cxxopts.hpp"
#include <filesystem>

using namespace std;

void Commandline_Angran(int argc, char *argv[]);

int main(int argc, char *argv[])
{

	TruncatedTspline tt;
	try
	{
		cxxopts::Options options(argv[0], "CMU Surface Spline Generation");
		options
			.positional_help("[optional args]")
			.show_positional_help();

		int flag_fit = 0;
		int ngref = 0;
		int nlref = 0;

		vector<BezierElement> bzmesh;
		vector<int> IDBC;
		vector<double> gh;
		//tt.Run_UTSP_1(fn_in, fn_out, bzmesh, IDBC, gh);

		string fn_in;
		string fn_out;
		string fn_lref;

		options.add_options("General")("h,help", "Print help")("g,globalref", "Set the level of global refinement (Default value: 0)", cxxopts::value<int>(ngref))("i,input", "Input path for the mesh", cxxopts::value<std::string>(fn_in))
		/*("o,output","Output path", cxxopts::value<std::string>(fn_out))*/
		// ("E,extrusion", "Extrusion mode", cxxopts::value<bool>(flag_extrude))("t,thick", "Extrusion thickness", cxxopts::value<double>(thick))("d,dir", "Extrusion direction: 0: Extrude from mid surface to top and bottom directions; 1: Extrude from bottom to top; -1: Extrude from top to bottom (Default value: 0)", cxxopts::value<int>(dir))("n,nel", "Number of element in thickness (Default value: 1)", cxxopts::value<int>(nel))("p,poly", "The degree of element in thickness (Default value: 3)", cxxopts::value<int>(p))

#ifdef CXXOPTS_USE_UNICODE
			("unicode", u8"A help option with non-ascii: ��. Here the size of the"
						" string should be correct")
#endif
			;

		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			cout << options.help({"General"}) << endl;
			exit(0);
		}
		if (result.count("input"))
		{
			fn_out = fn_in;
			tt.Run_UTSP_AS_Angran(fn_in, fn_out, ngref, nlref, flag_fit, fn_lref, bzmesh, IDBC, gh);

			Laplace lap;
			lap.VisualizeVTK(fn_out + "/bz", bzmesh);

			exit(0);
		}
	}
	catch (const cxxopts::OptionException &e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}

	return 0;
}