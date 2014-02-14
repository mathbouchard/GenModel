/*
 * GenModelCoin.h
 *
 *  Created on: 2011-03-18
 *      Author: mbouchard
 */

#ifndef GENMODELCOIN_H_
#define GENMODELCOIN_H_

#if defined WIN64 || defined WIN32
	#ifndef snprintf
		#define snprintf sprintf_s
	#endif
#endif

#include "GenModel.h"
#include "ClpSimplex.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#include "CoinBuild.hpp"
#include "CoinModel.hpp"
#include "ClpInterior.hpp"

using namespace std;

class CoinData
{
public:
	CoinData();
	~CoinData();
	long Reset();
	long Delete();
	ClpSimplex*  model;

	CoinBigIndex* mat_beg;
	int* mat_r;
	double* mat_v;
	CoinBigIndex* Q_beg;
	int* Q_r;
	double* Q_v;
	double* lrhs;
	double* urhs;
	double* ub;
	double* lb;
	char* type;
	double* obj;
	char** cname;
	char** rname;

	int* vnewi;
	double* vnewcol;

	long nc;
	long nr;
	long nq;
};

class GenModelCoin : public GenModel
{
public:
	~GenModelCoin() {if (solverdata != NULL) delete static_cast<CoinData*>(solverdata);}
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


#endif /* GENMODELCOIN_H_ */
