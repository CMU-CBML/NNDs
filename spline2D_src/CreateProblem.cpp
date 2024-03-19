#include "CreateProblem.h"
#include <fstream>
#include <iostream>

void CreateProblem::NURBSPatch(string fn)
{
	ofstream fout;
	fout.open(fn);
	if(fout.is_open())
	{
		fout.close();
	}
	else
	{
		cerr<<"Cannot open "<<fn<<"!\n";
	}
}