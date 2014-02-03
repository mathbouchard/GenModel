/***************************************************************************
 *  GenModel.h
 *  A generic model creation interface for different solver
 *
 *  October 5 11:32 2007
 *  Copyright  2007  Mathieu Bouchard
 *  mathbouchard@gmail.com
 ****************************************************************************/

#ifndef GenModelCplex_H
#define GenModelCplex_H

#if defined WIN64 || defined WIN32
	#ifndef snprintf
		#define snprintf sprintf_s
	#endif
#endif

#include "GenModel.h"
#include <ilcplex/cplex.h>

using namespace std;

class CplexData
{
public:
	CplexData();
	~CplexData();
	long Reset();
	long Delete();
    long ClearStructure();
	CPXENVptr		env;
	CPXLPptr		lp;
	int* mat_c;
	int* mat_r;
	double* mat_v;
	double* lrhs;
	double* urhs;
	char* sense;
	double* ub;
	double* lb;
	char* type;
	double* obj;
	double* x;
	double* dual;
	double* rcost;
	double* slack;
	char** cname;
	char** rname;
	long nc;
	long nr;
	long onc;
	long onr;
};

class GenModelCplex : public GenModel
{
public:
    GenModelCplex();
    ~GenModelCplex();
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
    long SwitchToMip();
    long SwitchToLp();
	long DeleteMipStarts();
	long Solve();
	long SetSol();
	long Clean();
	double GetMIPRelativeGap();
};

#endif // GenModelCplex_H
