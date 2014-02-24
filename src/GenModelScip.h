/***************************************************************************
 *  GenModelScip.h
 *  A generic model creation interface for different solver
 *
 *  October 5 11:32 2007
 *  Copyright  2007  Mathieu Bouchard
 *  mathbouchard@gmail.com
 ****************************************************************************/

#ifndef GenModelScip_H
#define GenModelScip_H

#if defined WIN64 || defined WIN32
	#ifndef snprintf
		#define snprintf sprintf_s
	#endif
#endif

#include "GenModel.h"
#include <scip/scip.h>
#include <scip/lpi.h>
#include <scip/scipdefplugins.h>

using namespace std;

class ScipData
{
public:
	ScipData();
	~ScipData();
	long Reset();
	long Delete();
    long ClearStructure();
	SCIP* scip;
    SCIP_VAR** vars;
    SCIP_CONS** cons;
    SCIP_VARTYPE* type;
    long nc;
    long qvar_index;
    long qcons_index;
	long nr;
    //char** cname;
	//char** rname;
};

class GenModelScip : public GenModel
{
public:
    GenModelScip();
    ~GenModelScip();
	long Init(string name);
    long CreateModel(string filename, int type=0, string dn="");
	long CreateModel();
	long AddSolverCol(vector<int>& ind, vector<double>& val, double obj, double lb, double ub, string name, char type = 'C');
	long AddSolverRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name);
	long AddCol(int* newi, double* newcol, int nz, double obj, double lb, double ub, const char* name, char type = 'C');
	long AddCut(int* cols, double* vals, int nz, double rhs, char sense, const char* name);
	long ChangeBulkBounds(int count, int * ind, char * type, double * vals);
	long ChangeBulkObjectives(int count, int * ind, double * vals);
    long ChangeBulkNz(int count, int* rind, int* cind, double* vals);
    long WriteProblemToLpFile(string filename);
    long WriteSolutionToFile(string filename);
    long SwitchToMip();
    long SwitchToLp();
	long DeleteMipStarts();
	long Solve();
	long SetSol();
	long Clean();
	double GetMIPRelativeGap();
    long SetDirectParam(string whichparam, genmodel_param value, string type, string message);
    long SetParam(string param, string whichparam, string type, string message, bool implemented = true);
};

#endif // GenModelScip_H
