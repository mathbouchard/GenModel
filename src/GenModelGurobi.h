/*
 * GenModelGurobi.h
 *
 *  Created on: 2011-04-12
 *      Author: mbouchard
 */

#ifndef GENMODELGUROBI_H_
#define GENMODELGUROBI_H_

#if defined WIN64 || defined WIN32
	#ifndef snprintf
		#define snprintf sprintf_s
	#endif
#endif

#include "GenModel.h"
#include "gurobi_c.h"
#include "gurobi_c++.h"

using namespace std;

class GurobiData
{
public:
	GurobiData();
	~GurobiData();
	long Reset();
	long Delete();
	GRBEnv*	env;
	GRBModel* model;
	GRBVar*	grb_v;
	double* ub;
	double* lb;
	long* equiv;
	char* type;
	double* obj;
	string* cname;
	long nc;
	long nr;
};

class GenModelGurobi : public GenModel
{
public:
	~GenModelGurobi() {if (solverdata != NULL) delete static_cast<GurobiData*>(solverdata);}
	long Init(string name);
	long CreateModel();
    long CreateModel(string filename, int type, string dn);
	long AddSolverRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name);
	long AddSolverCol(vector<int>& ind, vector<double>& val, double obj, double lb, double ub, string name, char type = 'C');
	long AddCut(int* cols, double* vals, int nz, double rhs, char sense, const char* name);
	long AddCol(int* newi, double* newcol, int nz, double obj, double lb, double ub, const char* name, char type = 'C');
	long Solve();
	long SetSol();
	long Clean();
    long SetDirectParam(int whichparam, genmodel_param value, string type, string message);
    long SetParam(string param, int whichparam, string type, string message, bool implemented = true);
};


#endif /* GENMODELGUROBI_H_ */
