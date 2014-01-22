#include "GenModelCplex.h"
//#include "ProblemReader.h"
#include <limits>

GenModelCplex::~GenModelCplex()
{
    if (solverdata != NULL) delete static_cast<CplexData*>(solverdata);
}

GenModelCplex::GenModelCplex() {
    
}

long GenModelCplex::Solve()
{
	CplexData* d = static_cast<CplexData*>(solverdata);
	int status = 0;
	if(boolParam.count("qp") > 0 && boolParam["qp"])
		status = CPXqpopt(d->env, d->lp);
	else if(boolParam.count("mip") > 0 && boolParam["mip"])
		status = CPXmipopt(d->env, d->lp);
	else if(strParam.count("algo") > 0 && strParam["algo"] == "interior")
		status = CPXbaropt(d->env, d->lp);
	else if(strParam.count("algo") > 0 && strParam["algo"] == "dual")
		status = CPXdualopt(d->env, d->lp);
	else if(strParam.count("algo") > 0 && strParam["algo"] == "primal")
		status = CPXprimopt(d->env, d->lp);
	else if(strParam.count("algo") > 0 && strParam["algo"] == "concurrent")
	{
		//printf("choosing concurrent algo\n");
		CPXsetintparam (d->env, CPX_PARAM_LPMETHOD, CPX_ALG_CONCURRENT);
		status = CPXlpopt(d->env, d->lp);
	}
	else if(strParam.count("algo") > 0 && strParam["algo"] == "sifting")
	{
		CPXsetintparam (d->env, CPX_PARAM_LPMETHOD, CPX_ALG_SIFTING);
		status = CPXlpopt(d->env, d->lp);
	}
	else
		status = CPXlpopt(d->env, d->lp);

	return 0;
}

long GenModelCplex::SetSol()
{
	vars.sol.clear();
	vars.sol.resize(vars.n,0);
	vars.rc.clear();
	vars.rc.resize(vars.n,0);
	CplexData* d = (CplexData*)solverdata;
	int status = 0;

	if(d->x != NULL)
		delete[] d->x;
	if(d->dual != NULL)
		delete[] d->dual;
	if(d->rcost != NULL)
		delete[] d->rcost;
	if(d->slack != NULL)
		delete[] d->slack;

	d->x = new double[nc];
	d->dual = new double[nr];
	d->rcost = new double[nc];
	d->slack = new double[nr];

	int tempstat = CPXgetstat (d->env, d->lp);
	int tempfeas;
	int tempdualfeas;
	int temphassol;
	int currmeth = CPXgetmethod(d->env, d->lp);

	CPXsolninfo(d->env, d->lp, &currmeth, &temphassol, &tempfeas, &tempdualfeas);

	feasible = static_cast<bool>(tempfeas);
	dualfeasible = static_cast<bool>(tempdualfeas);
	hassolution= static_cast<bool>(temphassol);

	if(!hassolution)
		return 1;

	if(boolParam.count("mip") > 0 && boolParam["mip"])
		status = CPXsolution (d->env, d->lp, &solstat, &objval, d->x, NULL, NULL, NULL);
	else
		status = CPXsolution (d->env, d->lp, &solstat, &objval, d->x, d->dual, d->slack, d->rcost);

	solstat = tempstat;

	for(long i = 0; i < long(nc); i++)
	{
		vars.sol[i] = d->x[i];
		vars.rc[i] = d->rcost[i];
	}

	for(long i = 0; i < long(nr); i++)
	{
		consts[i].dual = d->dual[i];
		consts[i].slack = d->slack[i];
	}

	return 0;
}

long GenModelCplex::AddSolverRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name)
{
	AddModelRow(ind, val, rhs, sense, name);
	AddCut(&ind[0], &val[0], int(ind.size()), rhs, sense, name.c_str());

	return 0;
}

long GenModelCplex::AddCut(int* cols, double* vals, int nz, double rhs, char sense, const char* name)
{
	CplexData* d = (CplexData*)solverdata;
	int rmatbeg = 0;

	CPXaddrows(d->env, d->lp, 0, 1, nz, &rhs, &sense, &rmatbeg, cols, vals, NULL, (char**)(&name));
	d->nr++;

	return 0;
}

long GenModelCplex::AddSolverCol(vector<int>& ind, vector<double>& val, double obj, double lb, double ub, string name, char type)
{
	AddModelCol(ind, val, obj, lb, ub, name, type);
	AddCol(&ind[0], &val[0], int(ind.size()), obj, lb, ub, name.c_str(), type);

	return 0;
}

long GenModelCplex::AddCol(int* newi, double* newcol, int nz, double obj, double lb, double ub, const char* name, char type)
{
	CplexData* d = (CplexData*)solverdata;
	int cmatbeg = 0;

	double clb = lb;
	if(clb == numeric_limits<double>::infinity())
		clb = CPX_INFBOUND;
	else if(clb == -numeric_limits<double>::infinity())
		clb = -CPX_INFBOUND;

	double cub = ub;
		if(cub == numeric_limits<double>::infinity())
			cub = CPX_INFBOUND;
		else if(cub == -numeric_limits<double>::infinity())
			cub = -CPX_INFBOUND;

	CPXaddcols(d->env, d->lp, 1, nz, &obj, &cmatbeg, newi, newcol, &clb, &cub, (char**)(&name));
	if(type != 'C')
	{
		int cind = d->nc;
		CPXchgctype(d->env, d->lp, 1, &cind, &type);
	}
	d->nc++;

	return 0;
}

long GenModelCplex::CreateModel()
{
	CplexData* d = (CplexData*)solverdata;
	int status = 0;
	d->nc = nc;
	d->nc = nc;
	d->onc = nc;
	d->onr = nr;

	if(boolParam.count("maximize") > 0 && boolParam["maximize"])
		CPXchgobjsen (d->env, d->lp, CPX_MAX);
	else
		CPXchgobjsen (d->env, d->lp, CPX_MIN);
	d->lrhs = new double[nr];
	d->urhs = new double[nr];
	d->sense = new char[nr];
	d->ub = new double[nc];
	d->lb = new double[nc];
	d->obj = new double[nc];
	d->type = new char[nc];
	d->mat_r = new int[nz];
	d->mat_c = new int[nz];
	d->mat_v = new double[nz];
	d->cname = new char*[nc];
	d->rname = new char*[nr];


	nz=0;
	for(unsigned long i = 0; i < nr; i++)
	{
		d->rname[i] = new char[consts[i].name.length()+1];
		snprintf(d->rname[i], consts[i].name.length()+1, "%s", consts[i].name.c_str());
		//printf("%ld %s: ", i, consts[i].name.c_str());
		for(unsigned long j = 0; j < consts[i].nz; j++)
		{
			d->mat_r[nz] = i;
			d->mat_c[nz] = consts[i].cols[j];
			d->mat_v[nz] = consts[i].coefs[j];
			//if(i >= 198)
				//printf("(%ld,%ld(%s),%f) ", d->mat_r[nz], d->mat_c[nz], vars.name[d->mat_c[nz]].c_str(), d->mat_v[nz]);
			nz++;
		}

		if(consts[i].lrhs == numeric_limits<double>::infinity())
			d->lrhs[i] = CPX_INFBOUND;
		else if(consts[i].lrhs == -numeric_limits<double>::infinity())
			d->lrhs[i] = -CPX_INFBOUND;
		else
			d->lrhs[i] = consts[i].lrhs;
		if(consts[i].urhs == numeric_limits<double>::infinity())
			d->urhs[i] = CPX_INFBOUND;
		else if(consts[i].urhs == -numeric_limits<double>::infinity())
			d->urhs[i] = -CPX_INFBOUND;
		else
			d->urhs[i] = consts[i].urhs-consts[i].lrhs;
		d->sense[i] = consts[i].sense;
	//	printf("%ld/%ld -> %c\n", i, nr, d->sense[i]);
	}
	for(unsigned long i = 0; i < nc; i++)
	{
		d->cname[i] = new char[vars.name[i].length()+1];
		snprintf(d->cname[i], vars.name[i].length()+1, "%s", vars.name[i].c_str());
		d->obj[i] = vars.obj[i];
		if(vars.ub[i] == numeric_limits<double>::infinity())
			d->ub[i] = CPX_INFBOUND;
		else if(vars.ub[i] == -numeric_limits<double>::infinity())
			d->ub[i] = -CPX_INFBOUND;
		else
			d->ub[i] = vars.ub[i];
		if(vars.lb[i] == numeric_limits<double>::infinity())
			d->lb[i] = CPX_INFBOUND;
		else if(vars.lb[i] == -numeric_limits<double>::infinity())
			d->lb[i] = -CPX_INFBOUND;
		else
			d->lb[i] = vars.lb[i];
		d->type[i] = vars.type[i];

		//printf("%ld (%s) -> %f %f %f %c\n", i, vars.name[i].c_str(), d->obj[i], d->lb[i], d->ub[i], d->type[i]);
	}
	status = CPXnewrows (d->env, d->lp, nr, d->lrhs, d->sense, d->urhs, d->rname);
	if ( status )
	{
		char  errmsg[1024];
		fprintf (stderr, "Could not create new rows.\n");
		CPXgeterrorstring (d->env, status, errmsg);
		fprintf (stderr, "%s", errmsg);
		return 1;
	}
	//else
		//printf("Row added!\n");

	if(boolParam.count("mip") > 0 && boolParam["mip"])
		status = CPXnewcols (d->env, d->lp, nc, d->obj, d->lb, d->ub, d->type, d->cname);
	else
		status = CPXnewcols (d->env, d->lp, nc, d->obj, d->lb, d->ub, NULL, NULL);
	if ( status )
	{
		char  errmsg[1024];
		fprintf (stderr, "Could not create new cols.\n");
		CPXgeterrorstring (d->env, status, errmsg);
		fprintf (stderr, "%s", errmsg);
		return 1;
	}
	//status = CPXnewcols (env, lp, nc, obj, lb, ub, NULL, colname);
	if ( status )
		return 1;
	//else
		//printf("Col added!\n");
	status = CPXchgcoeflist (d->env, d->lp, nz, d->mat_r, d->mat_c, d->mat_v);
	if ( status )
		return 1;

	vector<long>::iterator iti;
	vector<long>::iterator itj = vars.qj.begin();
	vector<double>::iterator itv = vars.qobj.begin();

	vector<vector<pair<int,double> > > qptemp;
	qptemp.resize(nc);
	int* qpbeg = NULL;
	int* qpnum = NULL;
	int* qpind = NULL;
	double* qpv = NULL;
	int qpnz = 0;

	if(!vars.qi.empty())
	{
		qpbeg = new int[nc];
		qpnum = new int[nc];
	}
	for(iti = vars.qi.begin(); iti != vars.qi.end(); iti++, itj++, itv++)
	{
		boolParam["qp"] = true;
		//status = CPXchgqpcoef (d->env, d->lp, *iti, *itj, *itv);

		qptemp[*iti].push_back(pair<int, double>(*itj,*itv));
		qpnz++;
		if(*iti != *itj)
		{
			qptemp[*itj].push_back(pair<int, double>(*iti,*itv));
			qpnz++;
		}
	}
	if(!vars.qi.empty())
	{
		qpv = new double[qpnz];
		qpind = new int[qpnz];

		qpnz=0;
		for(int i = 0; i < int(nc); i++)
		{
			qpbeg[i] = qpnz;
			qpnum[i] = int(qptemp[i].size());
			for(int j = 0; j < int(qptemp[i].size()); j++)
			{
				qpind[qpnz] = qptemp[i][j].first;
				qpv[qpnz] = qptemp[i][j].second;
				qpnz++;
			}
		}
		status = CPXcopyquad(d->env, d->lp, qpbeg, qpnum, qpind, qpv);
		delete[] qpbeg;
		delete[] qpnum;
		delete[] qpind;
		delete[] qpv;
	}
	if ( status )
	{
		printf("QP problem!\n");
		return 1;
	}
	//else
		//printf("Coefs added!\n");

	return 0;
}

long GenModelCplex::CreateModel(string filename, int type, string dn)
{
    //ReadFromFile(static_cast<GenModel*>(this), filename, type);
	CreateModel();
    SetNumbers();
    CreateModel();
    
    return 0;
}

long GenModelCplex::ChangeBulkBounds(int count, int * ind, char * type, double * vals)
{
	CplexData* d = (CplexData*)solverdata;

	for(long i = 0; i < count; i++)
	{
		if (type[i] == 'L' || type[i] == 'B')
		{
			vars.lb[i] = vals[i];
		}
		if (type[i] == 'U' || type[i] == 'B')
		{
			vars.ub[i] = vals[i];
		}
	}

	CPXchgbds(d->env, d->lp, count, ind, type, vals);

	return 0;
}

long GenModelCplex::ChangeBulkObjectives(int count, int * ind, double * vals)
{
	CplexData* d = (CplexData*)solverdata;

	for(long i = 0; i < count; i++)
	{
		vars.obj[i] = vals[i];
	}

	CPXchgobj(d->env, d->lp, count, ind, vals);

	return 0;
}

long GenModelCplex::ChangeBulkNz(int count, int* rind, int* cind, double * vals)
{
    CplexData* d = (CplexData*)solverdata;
    
    for(long i = 0; i < count; i++)
    {
        bool found = false;
        for(long j = 0; j < int(consts[rind[i]].cols.size()); j++)
        {
            if(consts[rind[i]].cols[j] == cind[i])
            {
                consts[rind[i]].coefs[j] = vals[i];
                found = true;
                break;
            }
        }
        if(!found)
            consts[rind[i]].AddNz(cind[i], vals[i]);
    }
    
    CPXchgcoeflist(d->env, d->lp, count, rind, cind, vals);
    
    return 0;
}



long GenModelCplex::DeleteMipStarts()
{
	CplexData* d = (CplexData*)solverdata;
	int n = CPXgetnummipstarts(d->env, d->lp);
	if (n > 0)
		CPXdelmipstarts(d->env, d->lp, 0, n - 1);

	return 0;
}

double GenModelCplex::GetMIPRelativeGap()
{
	CplexData* d = (CplexData*)solverdata;
	double gap, bestobjval = 0;
	CPXgetbestobjval(d->env, d->lp, &bestobjval);
	if (bestobjval > 0)	// If the optimal solution is found by the presolve, the CPXgetbestobjval = 0, and the CPXgetmiprelgap ~ 1
		CPXgetmiprelgap(d->env, d->lp, &gap);
	
	return gap;
}

long GenModelCplex::SwitchToMip()
{
    vector<int> ind;
    vector<char> type;
    for(int i = 0; i < int(vars.type.size()); i++)
    {
        if(vars.type[i] == 'B' || vars.type[i] == 'I' || vars.type[i] == 'S' || vars.type[i] == 'N')
        {
            ind.push_back(i);
            type.push_back(vars.type[i]);
        }
    }
    CplexData* d = static_cast<CplexData*>(solverdata);
    CPXchgctype(d->env, d->lp, int(ind.size()), &(ind[0]), &(type[0]));
    boolParam["mip"] = true;
    
    return 0;
}

long GenModelCplex::SwitchToLp()
{
    vector<int> ind;
    vector<char> type;
    for(int i = 0; i < int(vars.type.size()); i++)
    {
        if(vars.type[i] == 'B' || vars.type[i] == 'I' || vars.type[i] == 'S' || vars.type[i] == 'N')
        {
            ind.push_back(i);
            type.push_back('C');
        }
    }
    CplexData* d = static_cast<CplexData*>(solverdata);
    CPXchgctype(d->env, d->lp, int(ind.size()), &(ind[0]), &(type[0]));
    boolParam["mip"] = false;
    
    return 0;
}

long GenModelCplex::Init(string name)
{
	if(solverdata == NULL)
		solverdata = new CplexData();
	else
	{
		static_cast<CplexData*>(solverdata)->Delete();
		static_cast<CplexData*>(solverdata)->Reset();
	}

	CplexData* d = static_cast<CplexData*>(solverdata);
	int status = 0;

	d->env = CPXopenCPLEX (&status);

	// If an error occurs
	if ( d->env == NULL )
	{
		char  errmsg[1024];
		fprintf (stderr, "Could not open CPLEX environment.\n");
		CPXgeterrorstring (d->env, status, errmsg);
		fprintf (stderr, "%s", errmsg);
		return 1;
	}

	hassolution = false;

	// Turn on output to the screen
	if(boolParam.count("screenoff") > 0 && boolParam["screenoff"])
		status = CPXsetintparam (d->env, CPX_PARAM_SCRIND, CPX_OFF);
	else
		status = CPXsetintparam (d->env, CPX_PARAM_SCRIND, CPX_ON);
	if ( status )
	{
		fprintf (stderr, "Failure to turn on screen indicator, error %d->\n", status);
		return 1;
	}

	// Set the time limit
	if(dblParam.count("timelimit") > 0)
		status = CPXsetdblparam (d->env, CPX_PARAM_TILIM, dblParam["timelimit"]);

	// Use cut callback
	if(boolParam.count("usecutcb") > 0 && boolParam["usecutcb"])
	{
		status = CPXsetintparam (d->env, CPX_PARAM_PRELINEAR, 0);
		status = CPXsetintparam (d->env, CPX_PARAM_MIPCBREDLP, 0);
	}

	// Cut pass
	if(longParam.count("cutpass"))
	{
		status = CPXsetintparam (d->env, CPX_PARAM_CUTPASS, longParam["cutpass"]);
	}

	// Pump level
	if(longParam.count("pumplevel"))
	{
		status = CPXsetintparam (d->env, CPX_PARAM_FPHEUR, longParam["pumplevel"]);
	}

	// MIP Emphasis
	if(longParam.count("mipemphasis"))
	{
		status = CPXsetintparam (d->env, CPX_PARAM_MIPEMPHASIS, longParam["mipemphasis"]);
	}

	// Probing level
	if(longParam.count("probinglevel"))
	{
		status = CPXsetintparam (d->env, CPX_PARAM_PROBE, longParam["probinglevel"]);
	}

	if ( status )
	{
		fprintf (stderr, "Failure to set cut callback parameters, error %d->\n", status);
		return 1;
	}

	if(dblParam.count("feastol"))
		status = CPXsetdblparam (d->env, CPX_PARAM_EPRHS, dblParam["feastol"]);
	if ( status )
	{
		fprintf (stderr, "Failure to change feasibility tolerance, error %d->\n", status);
		return 1;
	}
	if(dblParam.count("opttol"))
		status = CPXsetdblparam (d->env, CPX_PARAM_EPOPT, dblParam["opttol"]);
	if ( status )
	{
		fprintf (stderr, "Failure to change optimality tolerance, error %d->\n", status);
		return 1;
	}
	if(dblParam.count("marktol"))
		status = CPXsetdblparam (d->env, CPX_PARAM_EPMRK, dblParam["marktol"]);
	if ( status )
	{
		fprintf (stderr, "Failure to change Markowitz tolerance, error %d->\n", status);
		return 1;
	}
	
	if(longParam.count("threads"))
	{
		printf("Threads: %ld", longParam["threads"]);
		status = CPXsetintparam (d->env, CPX_PARAM_THREADS, longParam["threads"]);
		if ( status )
		{
			fprintf (stderr, "Failure to change the number of threads, error %d->\n", status);
			return 1;
		} 
	}

	// Turn off preprocessing
	if(boolParam.count("preprocoff") > 0 && boolParam["preprocoff"])
	{
		status = CPXsetintparam (d->env, CPX_PARAM_AGGFILL, 0);
		status = status && CPXsetintparam (d->env, CPX_PARAM_PREPASS, 0);
		status = status && CPXsetintparam (d->env, CPX_PARAM_AGGIND, CPX_OFF);
		status = status && CPXsetintparam (d->env, CPX_PARAM_DEPIND, 0);
		status = status && CPXsetintparam (d->env, CPX_PARAM_PRELINEAR, 0);
		status = status && CPXsetintparam (d->env, CPX_PARAM_PREDUAL, -1);
		status = status && CPXsetintparam (d->env, CPX_PARAM_REDUCE, 0);
		status = status && CPXsetintparam (d->env, CPX_PARAM_PREIND, CPX_OFF);
	}
	if ( status )
	{
		fprintf (stderr, "Failure to turn off preprocessing, error %d->\n", status);
		return 1;
	}

	// Turn on data checking
	if(boolParam.count("datacheckoff") > 0 && boolParam["datacheckoff"])
		status = CPXsetintparam (d->env, CPX_PARAM_DATACHECK, CPX_OFF);
	else
		status = CPXsetintparam (d->env, CPX_PARAM_DATACHECK, CPX_ON);

	if ( status )
	{
		fprintf (stderr, "Failure to turn on data checking, error %d->\n", status);
		return 1;
	}

	// Sets a relative tolerance on the gap between the best integer objective and the objective of the best node remaining (between 0.0 and 1.0)
	if(dblParam.count("epgap"))
	{
		//printf("setting epgap\n");
		status = CPXsetdblparam (d->env, CPX_PARAM_EPGAP, dblParam["epgap"]);
		if ( status )
		{
			fprintf (stderr, "Failure to set relative gap tolerance, error %d->\n", status);
			return 1;
		}
	}


	// Create the problem
	d->lp = CPXcreateprob (d->env, &status, name.c_str());
	if ( d->lp == NULL )
	{
		fprintf (stderr, "Failed to create LP.\n");
		return 1;
	}

	return 0;
}

long GenModelCplex::Clean()
{
	if(solverdata != NULL)
		delete static_cast<CplexData*>(solverdata);

	return 0;
}

long CplexData::Reset()
{
	mat_c = NULL;
	mat_r = NULL;
	mat_v = NULL;
	lrhs = NULL;
	urhs = NULL;
	sense = NULL;
	ub = NULL;
	lb = NULL;
	type = NULL;
	obj = NULL;
	x = NULL;
	dual = NULL;;
	rcost = NULL;
	slack = NULL;
	cname = NULL;
	rname = NULL;
	env = NULL;
	lp = NULL;

	return 0;
}

CplexData::CplexData()
{
	Reset();
}

CplexData::~CplexData()
{
	Delete();
}

long CplexData::ClearStructure()
{
	if(mat_c != NULL)
		delete[] mat_c;
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
	if(sense != NULL)
		delete[] sense;
	if(ub != NULL)
		delete[] ub;
	if(lb != NULL)
		delete[] lb;
	if(type != NULL)
		delete[] type;
	if(x != 0)
		delete[] x;
	if(dual != NULL)
		delete[] dual;
	if(rcost != NULL)
		delete[] rcost;
	if(slack != NULL)
		delete[] slack;
	if(cname != NULL)
	{
		for(long i = 0; i < onc; i++)
			delete[] cname[i];
		delete[] cname;
	}
	if(rname != NULL)
	{
		for(long i = 0; i < onr; i++)
			delete[] rname[i];
		delete[] rname;
	}
    Reset();
    
    return 0;
}

long CplexData::Delete()
{
	if(lp != NULL)
	{
		CPXfreeprob(env, &lp);
	}
	if(env != NULL)
	{
		CPXcloseCPLEX(&env);
	}

	ClearStructure();

	return 0;
}
