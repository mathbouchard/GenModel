#include "GenModelScip.h"
#ifdef OSI_MODULE
#include "ProblemReaderOsi.h"
#endif
#include <limits>

GenModelScip::~GenModelScip()
{
    if (solverdata != NULL) delete static_cast<ScipData*>(solverdata);
}

GenModelScip::GenModelScip() {
    
}

long GenModelScip::WriteProblemToLpFile(string filename)
{
    try {
        ScipData* d = static_cast<ScipData*>(solverdata);
        FILE* file = fopen(filename.c_str(), "w");
        SCIPprintOrigProblem(d->scip, file, "lp", FALSE);
        fclose(file);
    } catch (string e) {
        printf("Error : %s\n", e.c_str());
    } catch (...) {
        printf("GenModel - Unexpected exception\n");
    }
    return 0;
}

long GenModelScip::WriteSolutionToFile(string filename)
{
    try {
        ScipData* d = static_cast<ScipData*>(solverdata);
        FILE* file = fopen(filename.c_str(), "w");
        SCIP_CALL(SCIPprintBestSol(d->scip, file, false));
        fclose(file);
    } catch (string e) {
        printf("Error : %s\n", e.c_str());
    } catch (...) {
        printf("GenModel - Unexpected exception\n");
    }
    return 0;
}


long GenModelScip::Solve()
{
	ScipData* d = static_cast<ScipData*>(solverdata);
	SCIP_CALL(SCIPsolve(d->scip));
    //SCIP_CALL(SCIPfreeTransform(d->scip));
	return 0;
}

long GenModelScip::SetSol()
{
	vars.sol.clear();
	vars.sol.resize(vars.n,0);
	vars.rc.clear();
	vars.rc.resize(vars.n,0);
	ScipData* d = (ScipData*)solverdata;
    
	//feasible = static_cast<bool>(tempfeas);
	//dualfeasible = static_cast<bool>(tempdualfeas);
    
    //printf("num cons = %d\n", SCIPgetNConss(d->scip));

    printf("SetSol - Init\n");
    
	hassolution =  SCIPgetNSols(d->scip) > 0;
    solstat = SCIPgetStatus(d->scip);

	if(!hassolution)
		return 1;
    
    printf("SCIP INVALID= %f %d\n",  SCIP_INVALID, vars.n);
    
    SCIP_SOL* sol = SCIPgetBestSol(d->scip);
    SCIP_CALL(SCIPgetSolVals(d->scip, sol, nc, d->vars, &(vars.sol[0])));
    
    //double quad_part = SCIPgetSolVal(d->scip, sol, d->vars[d->qvar_index]);
    
    objval =  SCIPgetSolOrigObj(d->scip, sol);
    
     SCIP_Bool transformed;
    
    for(long i = 0; i < long(nc); i++)
	{
        SCIP_VAR* temp; // = d->vars[i];
        if(transformed)
            SCIP_CALL(SCIPgetTransformedVar(d->scip, d->vars[i], &temp));
        else
            temp = d->vars[i];
        assert(temp != NULL);
        assert(transformed == SCIPvarIsTransformed(temp));
        if(boolParam.count("mip") == 0 || !boolParam["mip"])
            vars.rc[i] = SCIPgetVarRedcost(d->scip, temp);
	}
    
    printf("SetSol - Set columns\n");
    
    //SCIP_VAR* temp;
    //SCIP_CALL(SCIPgetTransformedVar(d->scip, d->vars[d->qvar_index], &temp));
    
	for(long i = 0; i < long(nr); i++)
	{
        SCIP_CONS* temp; // = d->cons[i];
        if(transformed)
            SCIP_CALL(SCIPgetTransformedCons(d->scip, d->cons[i], &temp));
        else
             temp = d->cons[i];
        assert(temp != NULL);
        assert(transformed == SCIPconsIsTransformed(temp));
        if(boolParam.count("mip") == 0 || !boolParam["mip"])
            consts[i].dual = SCIPgetDualsolLinear(d->scip, temp);
        consts[i].slack =  consts[i].lrhs-SCIPgetActivityLinear(d->scip,temp,sol);
        
        //consts[i].slack = SCIPgetRhsLinear(d->scip, d->cons[i]);
		//consts[i].slack = slack[i];
	}
    
    printf("SetSol - Set rows\n");
    
    //SCIPinfoMessage(d->scip, NULL, "\nSolution:\n");
    //SCIP_CALL(SCIPprintSol(d->scip, sol, NULL, false));
    printf("Variable set objval = %f\n", objval);
    
    if(boolParam.count("print_version") > 0 && boolParam["print_version"])
        printf("*********** Genmodel version = %s ***********\n", version.c_str());

    
    printf("SetSol - End\n");
    
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
	
    SCIP_CALL(SCIPsetObjsense(d->scip,SCIP_OBJSENSE_MAXIMIZE));

	if(boolParam.count("maximize") > 0 && boolParam["maximize"])
		SCIP_CALL(SCIPsetObjsense(d->scip,SCIP_OBJSENSE_MAXIMIZE));
	else
		SCIP_CALL(SCIPsetObjsense(d->scip,SCIP_OBJSENSE_MINIMIZE));
	double lb;
    double ub;
	
    d->qvar_index = -1;
    d->qcons_index = -1;
    d->nr = nr;
    d->nc = nc;
    if(boolParam.count("qp_mat") > 0 && boolParam["qp_mat"] && !vars.qi.empty())
    {
        ++(d->nc);
        d->qvar_index = nc;
        ++(d->nr);
        d->qcons_index = nr;
    }
    d->cons = new SCIP_CONS*[d->nr];
    d->vars = new SCIP_VAR*[d->nc];
    d->type = new SCIP_VARTYPE[d->nc];
    
	nz=0;
    
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
        switch (vars.type[i])
        {
            case 'B':
                d->type[i] = SCIP_VARTYPE_BINARY;
                break;
            case 'I':
                d->type[i] = SCIP_VARTYPE_INTEGER;
                break;
            case 'Z':
                d->type[i] = SCIP_VARTYPE_IMPLINT;
                break;
            case 'C':
                d->type[i] = SCIP_VARTYPE_CONTINUOUS;
                break;
            default:
                d->type[i] = SCIP_VARTYPE_CONTINUOUS;
                break;
        }
        SCIP_CALL(SCIPcreateVarBasic(d->scip, &(d->vars[i]), vars.name[i].c_str(), lb, ub, vars.obj[i], d->type[i]));
        SCIP_CALL(SCIPaddVar(d->scip, d->vars[i]));
	}
    
	for(unsigned long i = 0; i < nr; i++)
	{
        vector<SCIP_VAR*> cols(consts[i].nz);
		for(unsigned long j = 0; j < consts[i].nz; j++)
		{
            cols[j] = d->vars[consts[i].cols[j]];
			nz++;
		}
        if(consts[i].sense == 'R')
        {
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
        }
        else if(consts[i].sense == 'G')
        {
            ub = SCIPinfinity(d->scip);
            if(consts[i].lrhs == numeric_limits<double>::infinity())
                lb = SCIPinfinity(d->scip);
            else if(consts[i].lrhs == -numeric_limits<double>::infinity())
                lb = -SCIPinfinity(d->scip);
            else
                lb = consts[i].lrhs;
        }
        else if(consts[i].sense == 'E')
        {
            if(consts[i].lrhs == numeric_limits<double>::infinity())
                lb = SCIPinfinity(d->scip);
            else if(consts[i].lrhs == -numeric_limits<double>::infinity())
                lb = -SCIPinfinity(d->scip);
            else
                lb = consts[i].lrhs;
            ub = lb;
        }
        else
        {
            lb = -SCIPinfinity(d->scip);
            if(consts[i].lrhs == numeric_limits<double>::infinity())
                ub = SCIPinfinity(d->scip);
            else if(consts[i].lrhs == -numeric_limits<double>::infinity())
                ub = -SCIPinfinity(d->scip);
            else
                ub = consts[i].lrhs;
        }
        
        SCIP_CALL(SCIPcreateConsBasicLinear(d->scip, &d->cons[i], consts[i].name.c_str(), consts[i].nz, &(cols[0]), &(consts[i].coefs[0]), lb, ub));
        //SCIP_CALL(SCIPcreateConsLinear(d->scip, &d->cons[i], consts[i].name.c_str(), consts[i].nz, &(cols[0]), &(consts[i].coefs[0]), lb, ub,
        //                                    true,true,true,true,true,false,true,false,false,false));
        SCIP_CALL(SCIPaddCons(d->scip, d->cons[i]));
	}
	if(!vars.qi.empty())
        boolParam["qp"] = true;
    if(boolParam["qp_mat"])
    {
        vector<long>::iterator iti;
        vector<long>::iterator itj = vars.qj.begin();
        vector<double>::iterator itv = vars.qobj.begin();
        vector<SCIP_VAR*> qpi;
        vector<SCIP_VAR*> qpj;
        vector<double> qpv;
        int qpnz=0;
        
        for(iti = vars.qi.begin(); iti != vars.qi.end(); iti++, itj++, itv++)
        {
            if(*iti <= *itj)
            {
                qpi.push_back(d->vars[*iti]);
                qpj.push_back(d->vars[*itj]);
                qpv.push_back(*itv);
                qpnz++;
                printf("%d, %d = %f %d\n", *iti, *itj, *itv, qpnz);
                if(*iti != *itj)
                {
                    qpi.push_back(d->vars[*itj]);
                    qpj.push_back(d->vars[*iti]);
                    qpv.push_back(*itv);
                    qpnz++;
                    printf("%d, %d = %f %d\n", *iti, *itj, *itv, qpnz);
                }
            }
        }
        
        if(!vars.qi.empty())
        {
            SCIP_CALL(SCIPcreateVarBasic(d->scip, &(d->vars[d->qvar_index]), "quad_obj", -SCIPinfinity(d->scip), SCIPinfinity(d->scip), 1.0, SCIP_VARTYPE_CONTINUOUS));
            SCIP_CALL(SCIPaddVar(d->scip, d->vars[d->qvar_index]));
            // maximize => quad_obj = a_i x_i^2 : x^2-quad_obj = 0
            // minimize => quad_obj = a_i x_i^2 : x^2-quad_obj = 0
            
            double lincoefs = -1.0;
            //double qlhs = (boolParam.count("maximize") > 0 && boolParam["maximize"] ? 0.0 : -SCIPinfinity(d->scip));
            //double qrhs = (boolParam.count("maximize") > 0 && boolParam["maximize"] ? SCIPinfinity(d->scip) : 0.0);
            SCIP_CALL(SCIPcreateConsBasicQuadratic(d->scip, &(d->cons[d->qcons_index]), "qobj", 1, &(d->vars[d->qvar_index]), &lincoefs, qpnz, &(qpi[0]),
                                                   &(qpj[0]), &(qpv[0]), 0.0, 0.0));//qlhs, qrhs));
            SCIP_CALL(SCIPaddCons(d->scip, d->cons[d->qcons_index]));
        }
	}
    
	return 0;
}

long GenModelScip::CreateModel(string filename, int type, string dn)
{
#ifdef OSI_MODULE
    ReadFromFile(static_cast<GenModel*>(this), filename, type);
    SetNumbers();
    CreateModel();
#else
    throw string("Cannot use CreateModel(filenamem, type, dn) : Osi Module not present");
#endif
    return 0;
}

/*
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
 */

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
	
	return gap;*/
    return 0.0;
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
    try {
	if(solverdata == NULL)
		solverdata = new ScipData();
	else
	{
		static_cast<ScipData*>(solverdata)->Delete();
		static_cast<ScipData*>(solverdata)->Reset();
	}

	ScipData* d = static_cast<ScipData*>(solverdata);

        
	SCIP_CALL( SCIPcreate(&(d->scip)) );
    SCIP_CALL( SCIPincludeDefaultPlugins(d->scip) );
        
    // create empty problem
    SCIP_CALL( SCIPcreateProbBasic(d->scip, name.c_str()));
        
    SCIP_CALL( SCIPwriteParams(d->scip, "tmp/scip_before_params.set", TRUE, FALSE) );
    
    SetParam("log_file", "", "str", "Failure to set the log file", false);
    
    // General settings
    SetParam("log_output_stdout", "", "bool", "Failure to turn on/off log output to stdout", false);
    SetParam("log_level", "display/verblevel", "long", "Failure to set log level");
    SetParam("use_data_checking", "", "bool", "Failure to turn on/off data checking", false);
    SetParam("nb_threads", "lp/threads", "long", "Failure to set the number of threads");
    SetParam("use_preprocessor", "presolving/maxrounds", "bool", "Failure to use preprocessor");
        
    
    // MIP settings
    SetParam("nb_cut_pass", "", "long", "Failure to set the number of cut pass", false);
    SetParam("feasibility_pump_level", "", "long", "Failure to set the feasibility pump level", false);
    SetParam("probing_level", "", "long", "Failure to set the probing level", false);
    SetParam("mip_emphasis", "", "long", "Failure to set the MIP emphasis", false);
    SetParam("use_cut_callback", "", "bool", "Failure to use preprocessor", false);
    
    // Tolerance and limits
    SetParam("time_limit", "limits/time", "dbl", "Failure to set time limit");
    SetParam("max_iteration_limit", "lp/iterlim", "long", "Failure to set the maximal number of simplex iterations");
    SetParam("bounds_feasibility_tolerance", "numerics/feastol", "dbl", "Failure to set bounds feasibility tolerance");
    SetParam("bounds_feasibility_tolerance", "numerics/dualfeastol", "dbl", "Failure to set bounds feasibility tolerance");
    SetParam("optimality_tolerance", "", "dbl", "Failure to set optimality tolerance", false);
    SetParam("markowitz_tolerance", "", "dbl", "Failure to set Markowitz tolerance", false);
    SetParam("absolute_mip_gap_tolerance", "limits/absgap", "dbl", "Failure to set absolute gap tolerance");
    SetParam("relative_mip_gap_tolerance", "limits/gap", "dbl", "Failure to set relative gap tolerance");
    SetParam("lp_objective_limit", "", "dbl", "Failure to set lp objective limit", false);
    SetParam("lp_objective_limit", "", "dbl", "Failure to set lp objective limit", false);
        
    SCIP_CALL( SCIPwriteParams(d->scip, "tmp/scip_after_params.set", TRUE, FALSE) );
    
    //See http://scip.zib.de/doc/html/PARAMETERS.shtml
 
    } catch (string e) {
        printf("Error : %s\n", e.c_str());
    }
    
    
    
	return 0;
}

long GenModelScip::SetDirectParam(string whichparam, genmodel_param value, string type, string message)
{
    SCIP_RETCODE status = SCIP_OKAY;
    if(type == "dbl")
        status = SCIPsetRealParam((static_cast<ScipData*>(solverdata))->scip, whichparam.c_str(), value.dblval);
    else if(type == "long")
        SCIPsetIntParam((static_cast<ScipData*>(solverdata))->scip, whichparam.c_str(), value.longval);
    else if(type == "str")
        SCIPsetStringParam((static_cast<ScipData*>(solverdata))->scip, whichparam.c_str(), value.strval);
    else if(type == "char")
        SCIPsetCharParam((static_cast<ScipData*>(solverdata))->scip, whichparam.c_str(), (value.strval)[0]);
    if ( status != SCIP_OKAY )
        return ThrowError(message);
    return 0;
}

long GenModelScip::SetParam(string param, string whichparam, string type, string message, bool implemented)
{
    bool notimplmessage = boolParam.count("throw_on_unimplemeted_option") > 0 && boolParam["throw_on_unimplemeted_option"];
    
    if(type == "dbl")
    {
        if(dblParam.count(param) > 0 && implemented)
            SetDirectParam(whichparam, dbl2param(dblParam[param]), type, message);
        else if(notimplmessage && !implemented && dblParam.count(param) > 0)
            throw (string("Parameter ")+param+" not implemented in GenModelScip");
    }
    else if(type == "long")
    {
        if(longParam.count(param) > 0 && implemented)
            SetDirectParam(whichparam, long2param(longParam[param]), type, message);
        else if(notimplmessage && !implemented && longParam.count(param) > 0)
            throw (string("Parameter ")+param+" not implemented in GenModelScip");
    }
    else if(type == "str")
    {
        if(strParam.count(param) > 0 && implemented)
            SetDirectParam(whichparam, str2param(strParam[param]), type, message);
        else if(notimplmessage && !implemented && strParam.count(param) > 0)
            throw (string("Parameter ")+param+" not implemented in GenModelScip");
    }
    else if(type == "char")
    {
        if(strParam.count(param) > 0 && implemented)
            SetDirectParam(whichparam, str2param(strParam[param]), type, message);
        else if(notimplmessage && !implemented && strParam.count(param) > 0)
            throw (string("Parameter ")+param+" not implemented in GenModelScip");
    }
    else if(type == "bool")
    {
        if(boolParam.count(param) > 0 && implemented)
        {
            if(boolParam[param])
                SetDirectParam(whichparam, long2param(0), "long", message);
            else
                SetDirectParam(whichparam, long2param(-1), "long", message);
        }
        else if(notimplmessage && !implemented && boolParam.count(param) > 0)
            throw (string("Parameter ")+param+" not implemented in GenModelScip");
    }
    
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
    if(vars != NULL)
    {
        for(int i = 0; i < nc; i++)
            SCIP_CALL(SCIPreleaseVar(scip,&vars[i]));
    }
    if(cons != NULL)
    {
        for(int i = 0; i < nr; i++)
            SCIP_CALL(SCIPreleaseCons(scip,&cons[i]));
    }
    SCIP_CALL( SCIPfree(&scip) );
    Reset();
    
    return 0;
}

long ScipData::Delete()
{
	ClearStructure();

	return 0;
}
