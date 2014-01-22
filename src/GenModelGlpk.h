/*
 * GenModelGlpk.h
 *
 *  Created on: 2011-03-18
 *      Author: mbouchard
 */

#ifndef GENMODELGLPK_H_
#define GENMODELGLPK_H_

#if defined WIN64 || defined WIN32
	#ifndef snprintf
		#define snprintf sprintf_s
	#endif
#endif

#include "GenModel.h"
#include "glpk.h"

using namespace std;

class GlpkData
{
public:
	GlpkData();
	~GlpkData();
	long Reset();
	long Delete();
	glp_prob*  model;
	glp_iptcp interior_param;
	glp_smcp simplex_param;
	glp_iocp mip_param;
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
	char** cname;
	char** rname;
	long nc;
	long nr;
};

class GenModelGlpk : public GenModel
{
public:
	~GenModelGlpk() {if (solverdata != NULL) delete static_cast<GlpkData*>(solverdata);}
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
};


#endif /* GENMODELGLPK_H_ */
