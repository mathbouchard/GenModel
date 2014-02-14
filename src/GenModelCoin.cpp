/*
 * GenModelCoin.cpp
 *
 *  Created on: 2011-03-18
 *      Author: mbouchard
 */

#include "GenModelCoin.h"
#include "Coin_C_defines.h"
#include "ClpSimplex.hpp"
#include "ClpInterior.hpp"
#include "ClpCholeskyBase.hpp"
//#include "ProblemReader.h"
#include <limits>

using namespace std;

long GenModelCoin::Solve()
{
	CoinData* d = static_cast<CoinData*>(solverdata);

	if(strParam.count("algo") > 0 && strParam["algo"] == "interior")
	{
		ClpInterior intmod;
		if(boolParam.count("qp") > 0 && boolParam["qp"])
		{
			ClpCholeskyBase * cholesky = new ClpCholeskyBase();
			cholesky->setKKT(true);
			intmod.borrowModel(*(d->model));
			intmod.setCholesky(cholesky);
			intmod.primalDual();
			intmod.returnModel(*(d->model));
		}
		else
		{
			printf("interior\n");
			ClpCholeskyBase * cholesky = new ClpCholeskyBase();
			intmod.borrowModel(*(d->model));
			intmod.setCholesky(cholesky);
			intmod.primalDual();
			intmod.returnModel(*(d->model));
		}
	}
	else if(strParam.count("algo") > 0 && strParam["algo"] == "dual")
		d->model->dual();
	else
		d->model->primal();

	return 0;
}

long GenModelCoin::SetSol()
{
	vars.sol.clear();
	vars.sol.resize(vars.n,0);
	vars.rc.clear();
	vars.rc.resize(vars.n,0);
	CoinData* d = (CoinData*)solverdata;

	double* tx;
	double* tdual;
	double* trcost;
	double* tslack;

	tx = d->model->primalColumnSolution();
	tdual = d->model->dualRowSolution();
	trcost = d->model->dualColumnSolution();
	tslack = d->model->primalRowSolution();
	solstat = d->model->status();

	objval=d->model->getObjValue();

	for(long i = 0; i < long(nc); i++)
	{
		vars.sol[i] = tx[i];
		vars.rc[i] = trcost[i];
	}
	for(long i = 0; i < long(nr); i++)
	{
		consts[i].dual = tdual[i];
		consts[i].slack = tslack[i];
	}

	return 0;
}

long GenModelCoin::CreateModel()
{
	CoinData* d = (CoinData*)solverdata;

	d->nc = nc;
	d->nr = nr;
	d->nq = 0;

	if(boolParam.count("maximize") > 0 && boolParam["maximize"])
		d->model->setOptimizationDirection(-1.0);

	d->lrhs = new double[nr];
	d->urhs = new double[nr];
	d->ub = new double[nc];
	d->lb = new double[nc];
	d->obj = new double[nc];
	d->type = new char[nc];
	d->mat_r = new int[nz];
	d->mat_beg = new CoinBigIndex[nc+1];
	d->mat_v = new double[nz];
	d->cname = new char*[nc];
	d->rname = new char*[nr];
	memset(d->mat_beg, 0, (nc+1)*sizeof(CoinBigIndex));
	memset(d->Q_beg, 0, (nc+1)*sizeof(CoinBigIndex));


	vector<vector<pair<int,double> > > tvect;
	tvect.resize(nc);
	for(unsigned long i = 0; i < nr; i++)
	{
		d->rname[i] = new char[consts[i].name.length()+1];
		snprintf(d->rname[i], consts[i].name.length()+1, "%s", consts[i].name.c_str());
		for(unsigned long j = 0; j < consts[i].nz; j++)
		{
			d->mat_beg[consts[i].cols[j]]++;
			tvect[consts[i].cols[j]].push_back(pair<int,double>(i,consts[i].coefs[j]));
		}

		if(consts[i].lrhs == numeric_limits<double>::infinity())
			d->lrhs[i] = COIN_DBL_MAX;
		else if(consts[i].lrhs == -numeric_limits<double>::infinity())
			d->lrhs[i] = -COIN_DBL_MAX;
		else
			d->lrhs[i] = consts[i].lrhs;
		if(consts[i].urhs == numeric_limits<double>::infinity())
			d->urhs[i] = COIN_DBL_MAX;
		else if(consts[i].urhs == -numeric_limits<double>::infinity())
			d->urhs[i] = -COIN_DBL_MAX;
		else
			d->urhs[i] = consts[i].urhs;
		if(consts[i].sense == 'G')
		{
			//d->lrhs[i] = d->lrhs[i];
			d->urhs[i] = COIN_DBL_MAX;
		}
		else if(consts[i].sense == 'L')
		{
			d->urhs[i] = d->lrhs[i];
			d->lrhs[i] = -COIN_DBL_MAX;
		}
		else if(consts[i].sense == 'E')
		{
			d->urhs[i] = d->lrhs[i];
			//d->lrhs[i] = d->lrhs[i];
		}
	}
	int begcsum = 0;
	nz=0;
	for(unsigned long i = 0; i < nc; i++)
	{
		int temp = begcsum;
		begcsum+=d->mat_beg[i];
		d->mat_beg[i]=temp;
		for(unsigned int k = 0; k < (unsigned int)(tvect[i].size()); k++)
		{
			d->mat_r[nz] = tvect[i][k].first;
			d->mat_v[nz] = tvect[i][k].second;
			nz++;
		}

		d->cname[i] = new char[vars.name[i].length()+1];
		snprintf(d->cname[i], vars.name[i].length()+1, "%s", vars.name[i].c_str());
		d->obj[i] = vars.obj[i];
		if(vars.ub[i] == numeric_limits<double>::infinity())
			d->ub[i] = COIN_DBL_MAX;
		else if(vars.ub[i] == -numeric_limits<double>::infinity())
			d->ub[i] = -COIN_DBL_MAX;
		else
			d->ub[i] = vars.ub[i];
		if(vars.lb[i] == numeric_limits<double>::infinity())
			d->lb[i] = COIN_DBL_MAX;
		else if(vars.lb[i] == -numeric_limits<double>::infinity())
			d->lb[i] = -COIN_DBL_MAX;
		else
			d->lb[i] = vars.lb[i];
		d->type[i] = vars.type[i];
	}

	d->mat_beg[nc]=begcsum;

	d->model->loadProblem(nc,nr,d->mat_beg,d->mat_r,d->mat_v, d->lb,d->ub,d->obj,d->lrhs, d->urhs);

	vector<long>::iterator iti;
	vector<long>::iterator itj = vars.qj.begin();
	vector<double>::iterator itv = vars.qobj.begin();
	for(iti = vars.qi.begin(); iti != vars.qi.end(); iti++, itj++, itv++)
	{
		boolParam["qp"] = true;
		d->Q_r[d->nq] = *iti;
		d->Q_beg[*itj]++;
		d->Q_v[d->nq] = *itv;
		d->nq++;
	}
	int Qcsum = 0;
	for(unsigned long i = 0; i < nc; i++)
	{
		int temp = Qcsum;
		Qcsum += d->Q_beg[i];
		d->Q_beg[i] = temp;
	}
	d->Q_beg[nc] = Qcsum;

	if(!vars.qi.empty())
		d->model->loadQuadraticObjective(nc,d->Q_beg,d->Q_r,d->Q_v);

	return 0;
}

long GenModelCoin::CreateModel(string filename, int type, string dn)
{
    //ReadFromFile(static_cast<GenModel*>(this), filename, type);
    SetNumbers();
    CreateModel();
    
    return 0;
}

long GenModelCoin::Init(string name)
{
	if(solverdata == NULL)
		solverdata = new CoinData();
	else
	{
		static_cast<CoinData*>(solverdata)->Delete();
		static_cast<CoinData*>(solverdata)->Reset();
	}

	CoinData* d = static_cast<CoinData*>(solverdata);

	d->model = new ClpSimplex();

	if(longParam.count("maxIter") > 0)
		d->model->setMaximumIterations(longParam["maxIter"]);
	else
		d->model->setMaximumIterations(99999999);

	if(longParam.count("maxTime") > 0)
		d->model->setMaximumSeconds(longParam["maxTime"]);
	else
		d->model->setMaximumSeconds(240*3600);

	if(dblParam.count("primalTol") > 0)
		d->model->setPrimalTolerance(dblParam["primalTol"]);
	else
		d->model->setPrimalTolerance(1e-9);

	if(dblParam.count("dualTol") > 0)
		d->model->setDualTolerance(dblParam["dualTol"]);
	else
		d->model->setDualTolerance(1e-9);

	if(boolParam.count("screenoff") > 0 && boolParam["screenoff"])
		d->model->setLogLevel(0);



	//d->modelPrimalDual->passInMessageHandler(mexprinter);

	return 0;
}

long GenModelCoin::AddSolverRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name)
{
	AddModelRow(ind, val, rhs, sense, name);
	AddCut(&ind[0], &val[0], int(ind.size()), rhs, sense, name.c_str());

	return 0;
}

long GenModelCoin::AddCut(int* cols, double* vals, int nz, double rhs, char sense, const char* name)
{
	CoinData* d = (CoinData*)solverdata;

	double lb = rhs;
	double ub = rhs;

	if(sense == 'L')
		lb = -COIN_DBL_MAX;
	else if(sense == 'G')
		ub = COIN_DBL_MAX;

	d->model->addRow(nz, cols, vals, lb, ub);
	d->nr++;

	return 0;
}

long GenModelCoin::AddSolverCol(vector<int>& ind, vector<double>& val, double obj, double lb, double ub, string name, char type)
{
	AddModelCol(ind, val, obj, lb, ub, name, type);
	// it's ok see: http://stackoverflow.com/questions/5241824/c-get-internal-pointer-to-a-vector
	AddCol(&ind[0], &val[0], int(ind.size()), obj, lb, ub, name.c_str(), type);

	return 0;
}


long GenModelCoin::AddCol(int* newi, double* newcol, int nz, double obj, double lb, double ub, const char* name, char type)
{
	CoinData* d = (CoinData*)solverdata;

	d->model->addColumn(nz, newi, newcol, lb, ub, obj);
	d->nc++;

	//delete[] start;

	return 0;
}

long GenModelCoin::Clean()
{
	if(solverdata != NULL)
		delete static_cast<CoinData*>(solverdata);

	return 0;
}

long CoinData::Reset()
{
	mat_beg = NULL;
	mat_r = NULL;
	mat_v = NULL;
	Q_beg = NULL;
	Q_r = NULL;
	Q_v = NULL;
	lrhs = NULL;
	urhs = NULL;
	ub = NULL;
	lb = NULL;
	type = NULL;
	obj = NULL;
	cname = NULL;
	rname = NULL;

	vnewi = NULL;
	vnewcol = NULL;

	return 0;
}

CoinData::CoinData()
{
	Reset();
}

CoinData::~CoinData()
{
	Delete();
}
long CoinData::Delete()
{
	if(mat_beg != NULL)
		delete[] mat_beg;
	if(mat_r != NULL)
		delete[] mat_r;
	if(mat_v != NULL)
		delete[] mat_v;
	if(Q_beg != NULL)
		delete[] Q_beg;
	if(Q_r != NULL)
		delete[] Q_r;
	if(Q_v != NULL)
		delete[] Q_v;
	if(lrhs != NULL)
		delete[] lrhs;
	if(obj != NULL)
		delete[] obj;
	if(urhs != NULL)
		delete[] urhs;
	if(ub != NULL)
		delete[] ub;
	if(lb != NULL)
		delete[] lb;
	if(type != NULL)
		delete[] type;

	if(cname != NULL)
	{
		for(long i = 0; i < nc; i++)
			delete[] cname[i];
	}
	delete[] cname;
	if(rname != NULL)
	{
		for(long i = 0; i < nr; i++)
			delete[] rname[i];
	}
	delete[] rname;

	if(vnewi != NULL)
		delete[] vnewi;
	if(vnewcol != NULL)
		delete[] vnewcol;

	return 0;
}
