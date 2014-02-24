/*
 * GenModelOsi.h
 *
 *  Created on: 2012-10-01
 *      Author: mbouchard
 */

#ifndef GENMODELOSI_H_
#define GENMODELOSI_H_

#if defined WIN64 || defined WIN32
	#ifndef snprintf
		#define snprintf sprintf_s
	#endif
#endif

#include "GenModel.h"
#include "OsiClpSolverInterface.hpp"
#include "OsiCpxSolverInterface.hpp"
#include "OsiSpxSolverInterface.hpp"
#include "OsiGrbSolverInterface.hpp"
#include "OsiGlpkSolverInterface.hpp"
#include "OsiVolSolverInterface.hpp"

using namespace std;

class OsiData
{
public:
	OsiData();
	~OsiData();
	long Reset();
	long Delete();
    
	OsiSolverInterface* model;
    CoinBigIndex* mat_beg;
	int* mat_r;
	double* mat_v;
	double* lrhs;
	double* urhs;
	double* ub;
	double* lb;
	int* typei;
    int* typec;
	double* obj;
	char** cname;
	char** rname;

    int solvertype;
	long nc;
	long nr;
	long nz;
};

class GenModelOsi : public GenModel
{
public:
	~GenModelOsi() {if (solverdata != NULL) delete static_cast<OsiData*>(solverdata);}
	long Init(string name, int type=0);
	long CreateModel();
    long CreateModel(string filename, int type=0, string dn="");
	long AddSolverRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name);
	long AddSolverCol(vector<int>& ind, vector<double>& val, double obj, double lb, double ub, string name, char type = 'C');
	long AddCut(int* cols, double* vals, int nz, double rhs, char sense, const char* name);
	long AddCol(int* newi, double* newcol, int nz, double obj, double lb, double ub, const char* name, char type = 'C');
	long Solve();
	long SetSol();
	long Clean();
};


#endif /* GENMODELOSI_H_ */
