/***************************************************************************
 *  GenModel.h
 *  A generic model creation interface for different solver
 *
 *  October 5 11:32 2009
 *  Copyright  2007  Mathieu Bouchard
 *  mathbouchard@gmail.com
 ****************************************************************************/

#ifndef GenModel_H
#define GenModel_H

#if defined WIN64 || defined WIN32
	#ifndef snprintf
		#define snprintf sprintf_s
	#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <time.h>
//#include <sys/time.h>
#include <stdarg.h>
#include <float.h>

using namespace std;

class ModVars
{
public:
	ModVars() {n=0; defub=DBL_MAX; deflb=-DBL_MAX;}
	long AddVar(string nn, double o, double l, double u, char t);
	long AddVars(string nn, long size, double o, double l, double u, char t);
	long SetQpCoef(long i, long j, double val);
	long Print();
	vector<string> name;
	vector<double> obj;
	vector<char> type;
	map<string, long> offset;
	vector<double> ub;
	vector<double> lb;
	vector<double> sol;
	vector<double> rc;
	vector<double> qobj;
	vector<long> qi;
	vector<long> qj;
	double defub;
	double deflb;

	unsigned long n;
};

class ModConsts
{
public:
	ModConsts() {nz=0; lrhs=0.0; sense = 'E'; urhs=0.0;}
	long AddNz(long c, double v);
	string name;
	vector<long> cols;
	vector<double> coefs;
	double dual;
	double slack;
	double lrhs;
	char sense;
	double urhs;
	long id;
	unsigned long nz;
};

class GenModel
{
public:
	GenModel();
	virtual ~GenModel()  { ClearStructure(); };
	long AddIndexToCoef(string coef, string index);
	long AddCoef(string coef);
	long CoefIndex(string coef, int nargs, ...);
	long AddConst(string cname);
	long AddConst(string cname, double rhs, char sense);
	long AddVar(string nn, double o, double l, double u, char t);
	long AddVars(string nn, long size, double o, double l, double u, char t);
	long AddModelCol(vector<int>& ind, vector<double>& val, double obj, double lb, double ub, string name, char type = 'C');
	long AddModelRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name);
	long SetQpCoef(long i, long j, double val);
	long AddNz(long row, long col, double val);
	long AddNzToLast(long col, double val);
	long SetNumbers();
	long ClearStructure();
	long PrintSol();
	long PrintSol(string v);
	long PrintSolNz();
	long PrintSolNz(string v);
	long PrintObjVal();
	long SetLongParam(string param, long val);
	long SetDblParam(string param, double val);
	long SetBoolParam(string param, bool val);
	long SetStrParam(string param, string val);
	virtual long Init(string name) = 0;
	virtual long CreateModel() = 0;
	virtual long Solve() = 0;
	virtual long SetSol() = 0;
    virtual long ChangeBulkBounds(int count, int * ind, char * type, double * vals);
	virtual long WriteProblemToLpFile(string filename);
	virtual long ChangeBulkObjectives(int count, int * ind, double * vals);
	virtual long DeleteMipStarts();
	virtual double GetMIPRelativeGap();
	
	vector<ModConsts> consts;
	map<string,long> ci;
	unsigned long nc;
	unsigned long nr;
	unsigned long nz;
	ModVars vars;

	double objval;
	int solstat;
	bool feasible;
	bool dualfeasible;
	bool hassolution;
	void* solverdata;
	map<string, long> longParam;
	map<string, double> dblParam;
	map<string, bool> boolParam;
	map<string, string> strParam;
};

#endif // GenModel_H
