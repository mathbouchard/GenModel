#include "GenModelScip.h"
//#include "ProblemReader.h"
#include <limits>

GenModelScip::~GenModelScip()
{
    if (solverdata != NULL) delete static_cast<ScipData*>(solverdata);
}

GenModelScip::GenModelScip() {
    
}

long GenModelScip::WriteProblemToLpFile(string filename)
{
    ScipData* d = static_cast<ScipData*>(solverdata);
    SCIPwriteOrigProblem(d->scip, filename.c_str(), filename.substr(filename.find_last_of('.')).c_str(),false);

    return 0;
}

long GenModelScip::Solve()
{
	ScipData* d = static_cast<ScipData*>(solverdata);
	SCIP_CALL(SCIPsolve(d->scip));
    SCIP_CALL(SCIPfreeTransform(d->scip));

	return 0;
}

long GenModelScip::SetSol()
{
	vars.sol.clear();
	vars.sol.resize(vars.n,0);
	vars.rc.clear();
	vars.rc.resize(vars.n,0);
	ScipData* d = (ScipData*)solverdata;

    vector<double> x;
    vector<SCIP_VAR*> ic(nc);
    vector<SCIP_CONS*> ir(nr);
    
	/*int tempstat = CPXgetstat (d->env, d->lp);
	int tempfeas;
	int tempdualfeas;
	int temphassol;
	int currmeth = CPXgetmethod(d->env, d->lp);
     */
	//CPXsolninfo(d->env, d->lp, &currmeth, &temphassol, &tempfeas, &tempdualfeas);
    
    SCIPgetSolVals	(d->scip, SCIPgetBestSol(d->scip), nc, &(ic[0]), &(x[0]));

	//feasible = static_cast<bool>(tempfeas);
	//dualfeasible = static_cast<bool>(tempdualfeas);
	//hassolution= static_cast<bool>(temphassol);

	if(!hassolution)
		return 1;

	/*if(boolParam.count("mip") > 0 && boolParam["mip"])
		status = CPXsolution (d->env, d->lp, &solstat, &objval, d->x, NULL, NULL, NULL);
	else
		status = CPXsolution (d->env, d->lp, &solstat, &objval, d->x, d->dual, d->slack, d->rcost);

	solstat = tempstat;
    */
	for(long i = 0; i < long(nc); i++)
	{
		vars.sol[i] = x[i];
	//	vars.rc[i] = d->rcost[i];
	}

	/*for(long i = 0; i < long(nr); i++)
	{
		consts[i].dual = d->dual[i];
		consts[i].slack = d->slack[i];
	}*/

	return 0;
}

long GenModelScip::AddSolverRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name)
{
    /*
	AddModelRow(ind, val, rhs, sense, name);
	AddCut(&ind[0], &val[0], int(ind.size()), rhs, sense, name.c_str());
    */
	return 0;
}

long GenModelScip::AddCut(int* cols, double* vals, int nz, double rhs, char sense, const char* name)
{
    /*
	ScipData* d = (ScipData*)solverdata;
	int rmatbeg = 0;

	CPXaddrows(d->env, d->lp, 0, 1, nz, &rhs, &sense, &rmatbeg, cols, vals, NULL, (char**)(&name));
	d->nr++;
    */
	return 0;
}

long GenModelScip::AddSolverCol(vector<int>& ind, vector<double>& val, double obj, double lb, double ub, string name, char type)
{
    /*
	AddModelCol(ind, val, obj, lb, ub, name, type);
	AddCol(&ind[0], &val[0], int(ind.size()), obj, lb, ub, name.c_str(), type);
    */
	return 0;
}

long GenModelScip::AddCol(int* newi, double* newcol, int nz, double obj, double lb, double ub, const char* name, char type)
{
    /*
	ScipData* d = (ScipData*)solverdata;
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
    */
	return 0;
}

long GenModelScip::CreateModel()
{
	ScipData* d = (ScipData*)solverdata;
	int status = 0;
	
    SCIPsetObjsense(d->scip,SCIP_OBJSENSE_MAXIMIZE);

	if(boolParam.count("maximize") > 0 && boolParam["maximize"])
		SCIPsetObjsense(d->scip,SCIP_OBJSENSE_MAXIMIZE);
	else
		SCIPsetObjsense(d->scip,SCIP_OBJSENSE_MINIMIZE);
	double lb;
    double ub;
    double obj;
    int type;
	
    d->cons = new SCIP_CONS*[nr];
    d->vars = new SCIP_VAR*[nc];
    //d->cname = new char*[nc];
	//d->rname = new char*[nr];

	nz=0;
    
	for(unsigned long i = 0; i < nr; i++)
	{
        vector<SCIP_VAR*> cols;
		for(unsigned long j = 0; j < consts[i].nz; j++)
		{
            cols[i] = d->vars[i];
			nz++;
		}
		if(consts[i].lrhs == numeric_limits<double>::infinity())
			lb = SCIPinfinity(d->scip);
		else if(consts[i].lrhs == -numeric_limits<double>::infinity())
			lb = -SCIPinfinity(d->scip);
		else
			lb = consts[i].lrhs;
		if(consts[i].urhs == numeric_limits<double>::infinity())
			ub = SCIPinfinity(d->scip);
		else if(consts[i].urhs == -numeric_limits<double>::infinity())
			ub = -SCIPinfinity(d->scip);
		else
			ub = consts[i].urhs-consts[i].lrhs;
		if(consts[i].sense == 'L')
            lb = -SCIPinfinity(d->scip);
        
        SCIPcreateConsBasicLinear(d->scip, &(d->cons), consts[i].name.c_str(), consts[i].cols.size(),
                                  &(cols[0]), &(consts[i].coefs[0]), lb, ub);
	}
	for(unsigned long i = 0; i < nc; i++)
	{
		if(vars.ub[i] == numeric_limits<double>::infinity())
			ub = SCIPinfinity(d->scip);
		else if(vars.ub[i] == -numeric_limits<double>::infinity())
			ub = -SCIPinfinity(d->scip);
		else
			ub = vars.ub[i];
		if(vars.lb[i] == numeric_limits<double>::infinity())
			lb = SCIPinfinity(d->scip);
		else if(vars.lb[i] == -numeric_limits<double>::infinity())
			lb = -SCIPinfinity(d->scip);
		else
			lb = vars.lb[i];
		d->type[i] = vars.type[i];
        switch (vars.type[i])
        {
            case 'B':
                type = 0;
                break;
            case 'I':
                type = 1;
                break;
            case 'Z':
                type = 2;
                break;
            case 'C':
                type = 3;
                break;
            default:
                type = 3;
                break;
        }
        SCIP_CALL(SCIPcreateVarBasic(d->scip, &(d->vars[i]), vars.name[i].c_str(), lb, ub, obj, type));
        SCIP_CALL(SCIPaddVar(d->scip, d->vars[i]));
	}
	//if(boolParam.count("mip") > 0 && boolParam["mip"])
	
	vector<long>::iterator iti;
	vector<long>::iterator itj = vars.qj.begin();
	//vector<double>::iterator itv = vars.qobj.begin();

    vector<SCIP_VARS*> qpi;
    vector<SCIP_VARS*> qpj;
    //vector<SCIP_REAL> qpv;
    int qpnz=0;

	if(!vars.qi.empty())
	{
        boolParam["qp"] = true;
		qpbeg = new int[nc];
		qpnum = new int[nc];
	}
    if(boolParam["qp_mat"])
    {
        for(iti = vars.qi.begin(); iti != vars.qi.end(); iti++, itj++, itv++)
        {
            if(*iti < *itj)
            {
                qpi.push_back(*iti);
                qpj.push_back(*itj);
                //qpv.push_back(*itv);
                qpnz++;
                if(*iti != *itj)
                {
                    qpi.push_back(*itj);
                    qpj.push_back(*iti);
                    //qpv.push_back(*itv);
                    qpnz++;
                }
            }
        }
        if(!vars.qi.empty())
        {
            SCIP_CALL(SCIPcreateVarBasic(d->scip, &quad_var, "qvar", 0, SCIPinfinity(d->scip), 1.0, 3));
            SCIP_CALL(SCIPaddVar(d->scip, quad_var);
            nc++;
            
            // maximize => y <= x^2 : 0 <= x^2-y <= Inf
            // minimize => y >= x^2 : -Inf <= x^2-y <= 0
            
            double lincoefs = -1.0;
            double qlhs = (boolParam.count("maximize") > 0 && boolParam["maximize"] ? 0.0 : -SCIPinfinity(d->scip));
            double qrhs = (boolParam.count("maximize") > 0 && boolParam["maximize"] ? SCIPinfinity(d->scip) : 0.0);
            SCIPcreateConsBasicQuadratic(d->scip, quad_cons, "qobj", 1, &quad_var, &lincoefs, qpnz, &(qpi[0]),
                                         &(qpj[0]), &(vars.qobj) /*&(qpv[0])*/, qlhs, qrhs);
        }
	}
                      
	return 0;
}

long GenModelScip::CreateModel(string filename, int type, string dn)
{
    name = dn;
    Init(dn);
    ScipData* d = (ScipData*)solverdata;
    SCIPreadProb(d->scip, filename, filename.substr(filename.find_last_of('.'));
	CreateModel();
    SetNumbers();
    CreateModel();
    
    return 0;
}

long GenModelScip::ChangeBulkBounds(int count, int * ind, char * type, double * vals)
{
	/*
    ScipData* d = (ScipData*)solverdata;

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
    */

	return 0;
}

long GenModelScip::ChangeBulkObjectives(int count, int * ind, double * vals)
{
	/*
    ScipData* d = (ScipData*)solverdata;

	for(long i = 0; i < count; i++)
	{
		vars.obj[i] = vals[i];
	}

	CPXchgobj(d->env, d->lp, count, ind, vals);
    */

	return 0;
}

long GenModelScip::ChangeBulkNz(int count, int* rind, int* cind, double * vals)
{
    /*
    ScipData* d = (ScipData*)solverdata;
    
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
    */
    
    return 0;
}



long GenModelScip::DeleteMipStarts()
{
    /*
	ScipData* d = (ScipData*)solverdata;
	int n = CPXgetnummipstarts(d->env, d->lp);
	if (n > 0)
		CPXdelmipstarts(d->env, d->lp, 0, n - 1);
    */
	return 0;
}

double GenModelScip::GetMIPRelativeGap()
{
    /*
	ScipData* d = (ScipData*)solverdata;
	double gap, bestobjval = 0;
	CPXgetbestobjval(d->env, d->lp, &bestobjval);
	if (bestobjval > 0)	// If the optimal solution is found by the presolve, the CPXgetbestobjval = 0, and the CPXgetmiprelgap ~ 1
		CPXgetmiprelgap(d->env, d->lp, &gap);
	*/
	return gap;
}

long GenModelScip::SwitchToMip()
{
    /*
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
    ScipData* d = static_cast<ScipData*>(solverdata);
    CPXchgctype(d->env, d->lp, int(ind.size()), &(ind[0]), &(type[0]));
    boolParam["mip"] = true;
    */
    return 0;
}

long GenModelScip::SwitchToLp()
{
    /*
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
    ScipData* d = static_cast<ScipData*>(solverdata);
    CPXchgctype(d->env, d->lp, int(ind.size()), &(ind[0]), &(type[0]));
    boolParam["mip"] = false;
    */
    return 0;
}

long GenModelScip::Init(string name)
{
	if(solverdata == NULL)
		solverdata = new ScipData();
	else
	{
		static_cast<ScipData*>(solverdata)->Delete();
		static_cast<ScipData*>(solverdata)->Reset();
	}

	ScipData* d = static_cast<ScipData*>(solverdata);
	int status = 0;

	SCIP_CALL( SCIPcreate(&(d->scip)) );
    SCIP_CALL( SCIPincludeDefaultPlugins(d->scip) );

    // create empty problem
    SCIP_CALL( SCIPcreateProbBasic(d->scip, name.c_str()));
    
	return 0;
}

long GenModelScip::Clean()
{
	if(solverdata != NULL)
		delete static_cast<ScipData*>(solverdata);

	return 0;
}

long ScipData::Reset()
{
	vars = NULL;
    cons = NULL;

	return 0;
}

ScipData::ScipData()
{
	Reset();
}

ScipData::~ScipData()
{
	Delete();
}

long ScipData::ClearStructure()
{
    ScipData* d = static_cast<ScipData*>(solverdata);
    if(vars != NULL)
    {
        for(int i = 0; i < nc; i++)
            SCIP_CALL(SCIPreleaseVar(d->vars[i]));
    }
    if(cons != NULL)
    {
        for(int i = 0; i < nr; i++)
            SCIP_CALL(SCIPreleaseCons(d->vars[i]));
    }
	SCIP_CALL( SCIPfree(&(d->scip)) );
    Reset();
    
    return 0;
}

long ScipData::Delete()
{
	ClearStructure();

	return 0;
}
