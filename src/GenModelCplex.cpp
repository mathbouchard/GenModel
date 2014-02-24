#include "GenModelCplex.h"
#ifdef OSI_MODULE
#include "ProblemReaderOsi.h"
#endif
#include <limits>

string getcplexerror(CPXENVptr env, int status)
{
    char  errmsg[4096];
    snprintf(errmsg, 4096, "CPLEX error %d : ", status);
    string prefix = string(errmsg);
    CPXgeterrorstring (env, status, errmsg);
    return prefix+string(errmsg);
}

GenModelCplex::~GenModelCplex()
{
    if (solverdata != NULL) delete static_cast<CplexData*>(solverdata);
}

GenModelCplex::GenModelCplex() {
    
}

long GenModelCplex::WriteProblemToLpFile(string filename)
{
    if(!bcreated)
        throw string("WriteProblemToLpFile() not available : Problem not created yet;");
    
    CplexData* d = static_cast<CplexData*>(solverdata);
    CPXwriteprob(d->env, d->lp, filename.c_str(), "LP");
    return 0;
}

long GenModelCplex::WriteSolutionToFile(string filename)
{
    if(!bcreated)
        throw string("WriteSolutionToFile() not available : Problem not created yet;");
    
    CplexData* d = static_cast<CplexData*>(solverdata);
    CPXsolwrite (d->env, d->lp, filename.c_str());
    return 0;
}

long GenModelCplex::Solve()
{
    if(!bcreated)
        throw string("Solve() not available : Problem not created yet");
	CplexData* d = static_cast<CplexData*>(solverdata);
	int status = 0;
    if(boolParam.count("mip") > 0 && boolParam["mip"])
        status = CPXmipopt(d->env, d->lp);
	else if(boolParam.count("qp") > 0 && boolParam["qp"])
		status = CPXqpopt(d->env, d->lp);
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
    if(!bcreated)
        throw string("SetSol() not available : Problem not created yet");
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

    
	status = CPXsolninfo(d->env, d->lp, &currmeth, &temphassol, &tempfeas, &tempdualfeas);
    if ( status )
        return ThrowError(getcplexerror(d->env, status)+string(". ")+string("Failure to set set solution (CPXsolninfo)"));

	feasible = static_cast<bool>(tempfeas);
	dualfeasible = static_cast<bool>(tempdualfeas);
	hassolution= static_cast<bool>(temphassol);

	if(!hassolution)
		return 0;

	if(boolParam.count("mip") > 0 && boolParam["mip"])
		status = CPXsolution (d->env, d->lp, &solstat, &objval, d->x, NULL, NULL, NULL);
	else
		status = CPXsolution (d->env, d->lp, &solstat, &objval, d->x, d->dual, d->slack, d->rcost);
    if ( status )
        return ThrowError(getcplexerror(d->env, status)+string(". ")+string("Failure to set set solution (CPXsolution)"));

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
    
    if(boolParam.count("print_version") > 0 && boolParam["print_version"])
        printf("*********** Genmodel version = %s ***********\n", version.c_str());

	return 0;
}

long GenModelCplex::AddSolverRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name)
{
    if(!bcreated)
        return ThrowError("AddSolverRow() not available : Problem not created yet");
	AddModelRow(ind, val, rhs, sense, name);
	AddCut(&ind[0], &val[0], int(ind.size()), rhs, sense, name.c_str());

	return 0;
}

long GenModelCplex::AddCut(int* cols, double* vals, int nz, double rhs, char sense, const char* name)
{
    if(!bcreated)
        return ThrowError("AddCut() not available : Problem not created yet");
	CplexData* d = (CplexData*)solverdata;
	int rmatbeg = 0;

	CPXaddrows(d->env, d->lp, 0, 1, nz, &rhs, &sense, &rmatbeg, cols, vals, NULL, (char**)(&name));
	d->nr++;

	return 0;
}

long GenModelCplex::AddSolverCol(vector<int>& ind, vector<double>& val, double obj, double lb, double ub, string name, char type)
{
    if(!bcreated)
        return ThrowError("AddSolverCol() not available : Problem not created yet");
	AddModelCol(ind, val, obj, lb, ub, name, type);
	AddCol(&ind[0], &val[0], int(ind.size()), obj, lb, ub, name.c_str(), type);

	return 0;
}

long GenModelCplex::AddCol(int* newi, double* newcol, int nz, double obj, double lb, double ub, const char* name, char type)
{
    if(!bcreated)
        return ThrowError("AddCol() not available : Problem not created yet");
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
    if(!binit)
        return ThrowError("CreateModel() not available : Problem not initialized yet");
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
        boolParam["qp"] = true;
		qpbeg = new int[nc];
		qpnum = new int[nc];
	}
    if(boolParam.count("qp_mat") == 0 || boolParam["qp_mat"])
    {
        for(iti = vars.qi.begin(); iti != vars.qi.end(); iti++, itj++, itv++)
        {
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
                    qpv[qpnz] = 2.0*qptemp[i][j].second;
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
	}
	//else
		//printf("Coefs added!\n");
    bcreated = true;

	return 0;
}

long AddQuadraticVector(vector<int>& ind, vector<double>& val)
{
    return 0;
}

long AddQuadraticMatrix(vector<int>& indi, vector<int>& indj, vector<double>& val)
{
    return 0;
}

long GenModelCplex::CreateModel(string filename, int type, string dn)
{
#ifdef OSI_MODULE
    ReadFromFile(static_cast<GenModel*>(this), filename, type);
    SetNumbers();
    CreateModel();
#else
    return ThrowError("Cannot use CreateModel(filenamem, type, dn) : Osi Module not present");
#endif
    return 0;
}

long GenModelCplex::ChangeBulkBounds(int count, int * ind, char * type, double * vals)
{
    if(!bcreated)
        return ThrowError("ChangeBulkBounds not available : Problem not created yet");
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
    if(!bcreated)
        return ThrowError("ChangeBulkObjectives() not available : Problem not created yet");
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
    if(!bcreated)
        return ThrowError("ChangeBulkNz() not available : Problem not created yet");
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
    if(!bcreated)
        return ThrowError("ChangeBulkNz() not available : Problem not created yet");
	CplexData* d = (CplexData*)solverdata;
	int n = CPXgetnummipstarts(d->env, d->lp);
	if (n > 0)
		CPXdelmipstarts(d->env, d->lp, 0, n - 1);

	return 0;
}

double GenModelCplex::GetMIPRelativeGap()
{
    if(!bcreated)
        return ThrowError("ChangeBulkNz() not available : Problem not created yet");
	CplexData* d = (CplexData*)solverdata;
	double gap = 0, bestobjval = 0;
	CPXgetbestobjval(d->env, d->lp, &bestobjval);
	if (bestobjval > 0)	// If the optimal solution is found by the presolve, the CPXgetbestobjval = 0, and the CPXgetmiprelgap ~ 1
		CPXgetmiprelgap(d->env, d->lp, &gap);
	
	return gap;
}

long GenModelCplex::SwitchToMip()
{
    if(!bcreated)
        return ThrowError("SwitchToMip() not available : Problem not created yet");
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
    if(!bcreated)
        return ThrowError("SwitchToLp() not available : Problem not created yet");
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
    
    //strParam.count("log_file")
    //dblParam.count("relative_mip_gap_tolerance")
    //dblParam.count("absolute_mip_gap_tolerance")
    //dblParam.count("time_limit")
    //dblParam.count("bounds_feasibility_tolerance")
    //dblParam.count("optimality_tolerance")
    //dblParam.count("markowitz_tolerance"))
    //longParam.count("threads")
    //longParam.count("cutpass")
    //longParam.count("pumplevel")
    //longParam.count("mipemphasis")
	//longParam.count("probinglevel")
    //longParam.count("max_iteration_limit");
	//boolParam.count("preprocoff") : turn on/off preprocessing
	//boolParam.count("datacheckoff")
    //boolParam.count("screen_output")
    //boolParam.count("usecutcb")
    
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
        return ThrowError(getcplexerror(d->env, status)+string(". ")+string("Could not open CPLEX environment"));
    
	
	hassolution = false;
    
    // Log file
    if(strParam.count("log_file") > 0)
    {
        d->cpxfileptr =  CPXfopen(strParam["log_file"].c_str(), "w");
        status = CPXsetlogfile(d->env, d->cpxfileptr);
        if ( status )
            return ThrowError(getcplexerror(d->env, status)+string(". ")+string("Failure to set the log file"));
    }
    
    // General settings
    boolParam["log_output_stdout"] = true;
    SetParam("log_output_stdout", CPX_PARAM_SCRIND, "bool", "Failure to turn on/off log output to stdout");
    SetParam("log_level", 0, "long", "Failure to set log level", false);
    SetParam("use_data_checking", CPX_PARAM_DATACHECK, "bool", "Failure to turn on/off data checking");
    SetParam("nb_threads", CPX_PARAM_THREADS, "long", "Failure to set the number of threads");
    if(boolParam.count("use_preprocessor") > 0 && !boolParam["use_preprocessor"])
	{
        SetDirectParam(CPX_PARAM_AGGFILL, long2param(0), "long", "Failure to use preprocessor (CPX_PARAM_AGGFILL)");
        SetDirectParam(CPX_PARAM_PREPASS, long2param(0), "long", "Failure to use preprocessor (CPX_PARAM_PREPASS)");
        SetDirectParam(CPX_PARAM_AGGIND, long2param(CPX_OFF), "long", "Failure to use preprocessor (CPX_PARAM_AGGIND)");
        SetDirectParam(CPX_PARAM_DEPIND, long2param(0), "long", "Failure to use preprocessor (CPX_PARAM_DEPIND)");
        SetDirectParam(CPX_PARAM_PRELINEAR, long2param(0), "long", "Failure to use preprocessor (CPX_PARAM_PRELINEAR)");
        SetDirectParam(CPX_PARAM_PREDUAL, long2param(-1), "long", "Failure to use preprocessor (CPX_PARAM_PREDUAL)");
        SetDirectParam(CPX_PARAM_REDUCE, long2param(0), "long", "Failure to use preprocessor (CPX_PARAM_REDUCE)");
        SetDirectParam(CPX_PARAM_PREIND, long2param(CPX_OFF), "long", "Failure to use preprocessor (CPX_PARAM_PREIND)");
	}
    
    // MIP settings
    SetParam("nb_cut_pass", CPX_PARAM_CUTPASS, "long", "Failure to set the number of cut pass");
    SetParam("feasibility_pump_level", CPX_PARAM_FPHEUR, "long", "Failure to set the feasibility pump level");
    SetParam("probing_level", CPX_PARAM_PROBE, "long", "Failure to set the probing level");
    SetParam("mip_emphasis", CPX_PARAM_MIPEMPHASIS, "long", "Failure to set the MIP emphasis");
    if(boolParam.count("use_cut_callback") > 0 && boolParam["use_cut_callback"])
    {
        SetDirectParam(CPX_PARAM_PRELINEAR, long2param(0), "long", "Failure to use cut callback (CPX_PARAM_PRELINEAR)");
        SetDirectParam(CPX_PARAM_MIPCBREDLP, long2param(0), "long", "Failure to use cut callback (CPX_PARAM_MIPCBREDLP)");
    }
    
    // Tolerance and limits
    SetParam("time_limit", CPX_PARAM_TILIM, "dbl", "Failure to set time limit");
    SetParam("max_iteration_limit", CPX_PARAM_ITLIM, "long", "Failure to set the maximal number of simplex iterations");
    SetParam("bounds_feasibility_tolerance", CPX_PARAM_EPRHS, "dbl", "Failure to set bounds feasibility tolerance");
    SetParam("optimality_tolerance", CPX_PARAM_EPOPT, "dbl", "Failure to set optimality tolerance");
    SetParam("markowitz_tolerance", CPX_PARAM_EPMRK, "dbl", "Failure to set Markowitz tolerance");
    SetParam("absolute_mip_gap_tolerance", CPX_PARAM_EPAGAP, "dbl", "Failure to set absolute gap tolerance");
    SetParam("relative_mip_gap_tolerance", CPX_PARAM_EPGAP, "dbl", "Failure to set relative gap tolerance");
    if(boolParam.count("maximize") > 0 && boolParam["maximize"])
        SetParam("lp_objective_limit", CPX_PARAM_OBJULIM, "dbl", "Failure to set lp objective limit");
    else
        SetParam("lp_objective_limit", CPX_PARAM_OBJLLIM, "dbl", "Failure to set lp objective limit");
    

	// Create the problem
	d->lp = CPXcreateprob (d->env, &status, name.c_str());
	if ( d->lp == NULL )
        return ThrowError(getcplexerror(d->env, status)+string(". ")+string("Failure to create Cplex optimization problem"));

    binit = true;
    
	return 0;
}

long GenModelCplex::SetDirectParam(int whichparam, genmodel_param value, string type, string message)
{
    int status = 0;
    if(type == "dbl")
        status = CPXsetdblparam (static_cast<CplexData*>(solverdata)->env, whichparam, value.dblval);
    else if(type == "long")
        status = CPXsetintparam (static_cast<CplexData*>(solverdata)->env, whichparam, value.longval);
    else if(type == "str")
        status = CPXsetstrparam (static_cast<CplexData*>(solverdata)->env, whichparam, value.strval);
    if ( status )
        return ThrowError(getcplexerror(static_cast<CplexData*>(solverdata)->env, status)+string(". ")+message);
    
    return 0;
}

long GenModelCplex::SetParam(string param, int whichparam, string type, string message, bool implemented)
{
    bool notimplmessage = boolParam.count("throw_on_unimplemeted_option") > 0 && boolParam["throw_on_unimplemeted_option"];
    
    if(type == "dbl")
    {
        if(dblParam.count(param) > 0 && implemented)
            SetDirectParam(whichparam, dbl2param(dblParam[param]), type, message);
        else if(notimplmessage && !implemented && dblParam.count(param) > 0)
            throw (string("Parameter ")+param+" not implemented in GenModelOsi");
    }
    else if(type == "long")
    {
        if(longParam.count(param) > 0 && implemented)
            SetDirectParam(whichparam, long2param(longParam[param]), type, message);
        else if(notimplmessage && !implemented && longParam.count(param) > 0)
            throw (string("Parameter ")+param+" not implemented in GenModelOsi");
    }
    else if(type == "str")
    {
        if(strParam.count(param) > 0 && implemented)
            SetDirectParam(whichparam, str2param(strParam[param]), type, message);
        else if(notimplmessage && !implemented && strParam.count(param) > 0)
            throw (string("Parameter ")+param+" not implemented in GenModelOsi");
    }
    else if(type == "bool")
    {
        if(boolParam.count(param) > 0 && implemented)
        {
            if(boolParam[param])
                SetDirectParam(whichparam, long2param(CPX_ON), "long", message);
            else
                SetDirectParam(whichparam, long2param(CPX_OFF), "long", message);
        }
        else if(notimplmessage && !implemented && boolParam.count(param) > 0)
            throw (string("Parameter ")+param+" not implemented in GenModelOsi");
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
    cpxfileptr = NULL;

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
    if(cpxfileptr != NULL)
    {
        CPXfclose(cpxfileptr);
    }

	ClearStructure();

	return 0;
}
