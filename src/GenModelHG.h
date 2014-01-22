/***************************************************************************
 *  GenModelHG.h
 *  A generic model creation interface for different solver
 *
 *  October 5 11:32 2007
 *  Copyright  2007  Mathieu Bouchard
 *  mathbouchard@gmail.com
 ****************************************************************************/

#ifndef GenModelHG_H
#define GenModelHG_H

#if defined WIN64 || defined WIN32
	#ifndef snprintf
		#define snprintf sprintf_s
	#endif
#endif

#include "GenModel.h"
#include "GenModelCplex.h"
#include "HyperGraph.h"
#include <ilcplex/cplex.h>
#include <math.h>

using namespace std;

class GenModelHG;

class PrimalProblem
{
public:
	PrimalProblem() {}
	~PrimalProblem() {}
    
	int Solve(vector<double>& dual, vector<int>& lb2he, vector<int>& ub2he, vector<vector<pair<int,double> > >& exc, vector<vector<pair<int,double> > >& exv);
	int Init(GenModelHG& gm, vector<vector<pair<int,double> > > exc, vector<double>& lb, vector<string>& lbname, vector<double>& ub, vector<string>& ubname);
	int AddVar(vector<double>& sol, double val, char type);
	int AddColAlt(vector<double>& sol, double val, char type);
	int CreateProblem();
    
    int lbindex;
    int ubindex;
    int cindex;
    vector<int> cons2extrac;
    
	GenModelCplex gm;
	string algo;
};

class HGData
{
public:
	HGData();
	~HGData();
    
    int AddOrigin(GenModelHG& gm);
    int CalcLayers();
    int SolveHgMp();
    int SolveLayer(uint i);
    long Delete();
    long Reset();
    
    string solvername;
	char tmp[4096];
	//PreciseTimer timer;
	//PreciseTimer dualtimer;
	//PreciseTimer primtimer;
	HyperGraph hg;
	
	vector<uint> insol;
	vector<double> inval;
    
	vector<int> he2lb;
    vector<int> he2ub;
	vector<int> lb2he;
    vector<int> ub2he;
    vector<double> lb;
	vector<string> lbname;
	vector<double> ub;
	vector<string> ubname;
	vector<vector<double> > currsol;
	vector<vector<double> > colrep;
	vector<double> currobj;
	vector<double> realobj;
    vector<bool> mippath;
    vector<char> mipedges;
    
	int extradeb;
    int coldeb;
    int norigin_index;
    int curr_layer;
    int incol;
    int numcol;
    double objmult;
    
	vector<vector<uint> > nlayers;
    
    vector<vector<pair<int,double> > > extrac;
    vector<vector<pair<int,double> > > extrav;
    
    
	vector<double> extcost;
    
	PrimalProblem prim;
    
	int probid;
    
};

class GenModelHG : public GenModel
{
public:
	~GenModelHG() {if (solverdata != NULL) delete static_cast<HGData*>(solverdata);}
	long Init(string name);
    long CreateModel(string filename, int type=0, string dn="");
	long CreateModel();
	long AddSolverCol(vector<int>& ind, vector<double>& val, double obj, double lb, double ub, string name, char type = 'C');
	long AddSolverRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name);
	long AddCol(int* newi, double* newcol, int nz, double obj, double lb, double ub, const char* name, char type = 'C');
	long AddCut(int* cols, double* vals, int nz, double rhs, char sense, const char* name);
	long Solve();
	long SetSol();
	long Clean();
};

#endif // GenModelHG_H
