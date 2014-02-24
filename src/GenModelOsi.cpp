/*
 * GenModelOsi.cpp
 *
 *  Created on: 2012-10-01
 *      Author: mbouchard
 */

#include "GenModelOsi.h"
#include "ProblemReader.h"
#include <limits>

using namespace std;

long GenModelOsi::Solve()
{
    OsiData* d = static_cast<OsiData*>(solverdata);
    if(boolParam.count("mip") > 0 && boolParam["mip"])
    {
        printf("B&b\n");
    	d->model->branchAndBound();
    }
    else if(boolParam.count("notroot") > 0 && boolParam["notroot"])
    	d->model->resolve();
    else
        d->model->initialSolve();
    
/*	CoinData* d = static_cast<CoinData*>(solverdata);

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
		d->model->primal();*/

	return 0;
}

long GenModelOsi::SetSol()
{
    OsiData* d = static_cast<OsiData*>(solverdata);
    
    vars.sol.clear();
	vars.sol.resize(vars.n,0);
	vars.rc.clear();
	vars.rc.resize(vars.n,0);
    
    const double* sol = d->model->getColSolution();
    const double* act = d->model->getRowActivity();
    
    if(boolParam.count("mip") == 0 || !boolParam["mip"])
    {
        const double* dual = d->model->getRowPrice();
        const double* rc = d->model->getReducedCost();
        for (uint i = 0; i < vars.n; i++)
            vars.rc[i] = rc[i];
        for(long i = 0; i < long(nr); i++)
            consts[i].dual = dual[i];
    }
    
    for (uint i = 0; i < vars.n; i++)
        vars.sol[i] = sol[i];
    for(long i = 0; i < long(nr); i++)
	{
        switch (consts[i].sense)
        {
            case 'L': consts[i].slack = consts[i].lrhs-act[i]; break;
            case 'G': consts[i].slack = act[i]-consts[i].lrhs; break;
            default: consts[i].slack = act[i]; break;
        }
	}
    
    solstat = d->model->isProvenOptimal();
    if (d->model->isProvenOptimal())
        solstat = 1;
    if (d->model->isProvenPrimalInfeasible())
        solstat = 2;
    if (d->model->isProvenDualInfeasible())
        solstat = 3;
    if (d->model->isIterationLimitReached())
        solstat = 4;
    
	objval = d->model->getObjValue();

	return 0;
}

long GenModelOsi::CreateModel()
{
    OsiData* d = static_cast<OsiData*>(solverdata);

    d->nc = nc;
	d->nr = nr;
    
	d->lrhs = new double[nr];
	d->urhs = new double[nr];
	d->ub = new double[nc];
	d->lb = new double[nc];
	d->obj = new double[nc];
	d->typei = new int[nc];
    d->typec = new int[nc];
	d->mat_r = new int[nz];
	d->mat_beg = new CoinBigIndex[nc+1];
	d->mat_v = new double[nz];
	d->cname = new char*[nc];
	d->rname = new char*[nr];
	memset(d->mat_beg, 0, (nc+1)*sizeof(CoinBigIndex));
    
    int numint = 0;
    int numcont = 0;
	
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
		//d->type[i] = vars.type[i];
        if(vars.type[i] == 'B')
        {
            if(d->lb[i] < 0.0)
                d->lb[i] = 0.0;
            if(d->ub[i] > 1.0)
                d->ub[i] = 1.0;
            d->typei[numint] = i;
            numint++;
        }
        else if(vars.type[i] == 'I')
        {
            d->typei[numint] = i;
            numint++;
        }
        else
        {
            d->typec[numcont] = i;
            numcont++;
        }
	}
    
	d->mat_beg[nc]=begcsum;
    
	d->model->loadProblem(nc,nr,d->mat_beg,d->mat_r,d->mat_v, d->lb,d->ub,d->obj,d->lrhs, d->urhs);
    
    d->model->setContinuous(d->typec, numcont);
    d->model->setInteger (d->typei, numint);
    
    if(boolParam.count("maximize") > 0 && boolParam["maximize"])
		d->model->setObjSense(-1.0);
    else
        d->model->setObjSense(1.0);
    
	return 0;
}

long GenModelOsi::CreateModel(string filename, int type, string dn)
{
    OsiData* d = static_cast<OsiData*>(solverdata);
    
    switch (type)
    {
        case 0: d->model->readMps(filename.c_str()); break;
        case 1: d->model->readLp(filename.c_str()); break;
        case 2: d->model->readGMPL(filename.c_str(), (dn == "" ? NULL : dn.c_str())); break;
        default: d->model->readMps(filename.c_str()); break;
    }

    //ReadFromObject(this, d->model);
    
    //ReadFromFile(static_cast<GenModel*>(this), filename, type);
    //SetNumbers();
    //CreateModel();
    
    return 0;
}

long GenModelOsi::Init(string name, int type)
{
	if(solverdata == NULL)
		solverdata = new OsiData();
	else
	{
		static_cast<OsiData*>(solverdata)->Delete();
		static_cast<OsiData*>(solverdata)->Reset();
	}

	OsiData* d = static_cast<OsiData*>(solverdata);

    d->solvertype = type;
    if(type==0)
        d->model = new OsiClpSolverInterface();
    else if(type==1)
        d->model = new OsiVolSolverInterface();
    else if(type==2)
        d->model = new OsiGlpkSolverInterface();
    else if(type==3)
        d->model = new OsiSpxSolverInterface();
    else if(type==4)
        d->model = new OsiCpxSolverInterface();
    else if(type==5)
        d->model = new OsiGrbSolverInterface();
    else 
        d->model = new OsiClpSolverInterface();

	if(longParam.count("maxIter") > 0)
        d->model->setIntParam( OsiMaxNumIteration, longParam["maxIter"]);
	else
        d->model->setIntParam( OsiMaxNumIteration, 99999999);
    
	if(dblParam.count("primalTol") > 0)
		d->model->setDblParam( OsiPrimalTolerance, dblParam["primalTol"]);
	else
		d->model->setDblParam( OsiPrimalTolerance, 1e-9);
	if(dblParam.count("dualTol") > 0)
		d->model->setDblParam( OsiDualTolerance, dblParam["dualTol"]);
	else
		d->model->setDblParam( OsiDualTolerance, 1e-9);

	//if(boolParam.count("screenoff") > 0 && boolParam["screenoff"])
	//	d->model->setLogLevel(0);
    
	return 0;
}

long GenModelOsi::AddSolverRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name)
{
	AddModelRow(ind, val, rhs, sense, name);
	AddCut(&ind[0], &val[0], int(ind.size()), rhs, sense, name.c_str());

	return 0;
}

long GenModelOsi::AddCut(int* cols, double* vals, int nz, double rhs, char sense, const char* name)
{
	OsiData* d = (OsiData*)solverdata;

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

long GenModelOsi::AddSolverCol(vector<int>& ind, vector<double>& val, double obj, double lb, double ub, string name, char type)
{
	AddModelCol(ind, val, obj, lb, ub, name, type);
	// it's ok see: http://stackoverflow.com/questions/5241824/c-get-internal-pointer-to-a-vector
	AddCol(&ind[0], &val[0], int(ind.size()), obj, lb, ub, name.c_str(), type);

	return 0;
}


long GenModelOsi::AddCol(int* newi, double* newcol, int nz, double obj, double lb, double ub, const char* name, char type)
{
	OsiData* d = (OsiData*)solverdata;
	d->model->addCol(nz, newi, newcol, lb, ub, obj);
	d->nc++;

	return 0;
}

long GenModelOsi::Clean()
{
	if(solverdata != NULL)
		delete static_cast<OsiData*>(solverdata);

	return 0;
}

long OsiData::Reset()
{
    model = NULL;
	mat_beg = NULL;
	mat_r = NULL;
	mat_v = NULL;
	lrhs = NULL;
	urhs = NULL;
	ub = NULL;
	lb = NULL;
	typei = NULL;
    typec = NULL;
	obj = NULL;
	cname = NULL;
	rname = NULL;

	return 0;
}

OsiData::OsiData()
{
	Reset();
}

OsiData::~OsiData()
{
	Delete();
}
long OsiData::Delete()
{
    //if(model != NULL)
	//	delete[] model;
	if(mat_beg != NULL)
		delete[] mat_beg;
	if(mat_r != NULL)
		delete[] mat_r;
	if(mat_v != NULL)
		delete[] mat_v;
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
	if(typei != NULL)
		delete[] typei;
    if(typec != NULL)
		delete[] typec;
    
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
	
	return 0;
}
