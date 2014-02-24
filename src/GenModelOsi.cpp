/*
 * GenModelOsi.cpp
 *
 *  Created on: 2012-10-01
 *      Author: mbouchard
 */

#include "GenModelOsi.h"
#include "ProblemReaderOsi.h"
#include <limits>

using namespace std;

long GenModelOsi::WriteProblemToLpFile(string filename)
{
    if(!bcreated)
        throw string("WriteProblemToLpFile() not available : Problem not created yet;");
    
    OsiData* d = static_cast<OsiData*>(solverdata);
    d->model->writeLp(filename.c_str());
    return 0;
}

long GenModelOsi::WriteSolutionToFile(string filename)
{
    if(!bcreated)
        throw string("WriteSolutionToFile() not available : Problem not created yet;");
    
    FILE* f = fopen(filename.c_str(), "w");
    for(long i = 0; i < long(vars.n); i++)
	{
		if(fabs(vars.sol[i]) > 0.000001)
			fprintf(f, "%s :		%f\n", vars.name[i].c_str(), vars.sol[i]);
	}
    
    fclose(f);
    return 0;
}

long GenModelOsi::Solve()
{
    if(!bcreated)
        throw string("WriteSolutionToFile() not available : Problem not created yet;");
    OsiData* d = static_cast<OsiData*>(solverdata);
    d->model->setHintParam(OsiDoReducePrint,true,OsiHintTry);
    
    
    /*
    if(boolParam.count("mip") > 0 && boolParam["mip"])
    {
        printf("********* OSI Branch and bounds **************\n");
    	d->model->branchAndBound();
    }
    else if(boolParam.count("notroot") > 0 && boolParam["notroot"])
    {
        printf("********* OSI Resolve **************\n");
    	d->model->resolve();
    }
    else
    {
        printf("********* OSI Initial Solve **************\n");
        d->model->initialSolve();
    }
    */
    
    d->model->setHintParam(OsiDoReducePrint,true,OsiHintTry);

    // Set up some cut generators and defaults
    // Probing first as gets tight bounds on continuous

    CglProbing generator1;
    // can not use objective
    generator1.setUsingObjective(false);
    generator1.setMaxPass(3);
    generator1.setMaxProbe(100);
    generator1.setMaxLook(50);
    generator1.setRowCuts(3);
    
    CglGomory generator2;
    // try larger limit
    generator2.setLimit(30000);
    
    CglKnapsackCover generator3;
    
    CglOddHole generator4;
    generator4.setMinimumViolation(0.005);
    generator4.setMinimumViolationPer(0.00002);
    // try larger limit
    generator4.setMaximumEntries(200);
    
    CglClique generator5;
    generator5.setStarCliqueReport(false);
    generator5.setRowCliqueReport(false);
    
    CglMixedIntegerRounding mixedGen;
    CglFlowCover flowGen;
    // Add in generators (just at root)

    d->mipmodel = new CbcModel(*(d->model));
    // Add in generators
    // Only some generators work (and even then try without first)
    d->mipmodel->addCutGenerator(&generator1,1,"Probing");
    // Allow rounding heuristic

    d->mipmodel->addCutGenerator(&generator1,-99,"Probing",true,false,false,-1);
    d->mipmodel->addCutGenerator(&generator2,-99,"Gomory",true,false,false,-99);
    d->mipmodel->addCutGenerator(&generator3,-99,"Knapsack",true,false,false,-99);
    d->mipmodel->addCutGenerator(&generator4,-99,"OddHole",true,false,false,-99);
    d->mipmodel->addCutGenerator(&generator5,-99,"Clique",true,false,false,-99);
    d->mipmodel->addCutGenerator(&flowGen,-99,"FlowCover",true,false,false,-99);
    d->mipmodel->addCutGenerator(&mixedGen,-99,"MixedIntegerRounding",true,false,false,-100);
    
    CbcRounding heuristic1(*(d->mipmodel));
    
    // do not add yet as don't know how to deal with quadratic objective
    d->mipmodel->addHeuristic(&heuristic1);
    
    CbcHeuristicGreedyCover heuristic2(*(d->mipmodel));
    // Use original upper and perturb more
    heuristic2.setAlgorithm(11);
    d->mipmodel->addHeuristic(&heuristic2);

    // Redundant definition of default branching (as Default == User)
    //CbcBranchUserDecision branch;
    //d->mipmodel->setBranchingMethod(&branch);

    // Definition of node choice
    //CbcCompareUser compare;
    // breadth first
    //compare.setWeight(0.0);
    //d->mipmodel->setNodeComparison(compare);


    // Do initial solve to continuous
    d->mipmodel->initialSolve();
    
    // Could tune more
    d->mipmodel->setMinimumDrop(CoinMin(1.0, fabs(d->mipmodel->getMinimizationObjValue())*1.0e-3+1.0e-4));
    
    d->mipmodel->setMaximumCutPassesAtRoot(50);
    d->mipmodel->setMaximumCutPasses(100);
    
    // Switch off strong branching if wanted
    //d->mipmodel->setNumberStrong(5);
    
    d->mipmodel->solver()->setIntParam(OsiMaxNumIterationHotStart,10000);
    
    // If time is given then stop after that number of minutes
    if (dblParam.count("time_limit") > 0)
    {
        printf("Stopping after %f seconds\n", dblParam["time_limit"]);
        d->mipmodel->setDblParam(CbcModel::CbcMaximumSeconds, dblParam["time_limit"]);
    }
    // Switch off most output
    if (d->mipmodel->getNumCols()<3000)
    {
        d->mipmodel->messageHandler()->setLogLevel(1);
        //d->model->solver()->messageHandler()->setLogLevel(0);
    }
    else
    {
        d->mipmodel->messageHandler()->setLogLevel(2);
        d->mipmodel->solver()->messageHandler()->setLogLevel(1);
    }
    d->mipmodel->setPrintFrequency(50);

    double time1 = CoinCpuTime();
    
    // Do complete search
    d->mipmodel->branchAndBound();
    
    printf(" took %f seconds, %d nodes with objective %f, %s\n",
           CoinCpuTime()-time1, d->mipmodel->getNodeCount(), d->mipmodel->getObjValue(), (!d->mipmodel->status() ? " Finished" : " Not finished"));

    // Print more statistics
    printf("Cuts at root node changed objective from %f to %f\n", d->mipmodel->getContinuousObjective(), d->mipmodel->rootObjectiveAfterCuts());
    
    int numberGenerators = d->mipmodel->numberCutGenerators();
    for (int iGenerator=0;iGenerator<numberGenerators;iGenerator++)
    {
        CbcCutGenerator * generator = d->mipmodel->cutGenerator(iGenerator);
        printf("%s was tried %d times and created %d cuts of which %d were active after adding rounds of cuts\n", generator->cutGeneratorName(),
               generator->numberTimesEntered(),generator->numberCutsInTotal(), generator->numberCutsActive());
    }
    
    hassolution = d->mipmodel->getIntParam(CbcModel::CbcMaxNumSol) > 0;
    objval = d->mult*d->mipmodel->getObjValue();
    
    d->model = d->mipmodel->solver();
    
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
    if(!bcreated)
        throw string("WriteSolutionToFile() not available : Problem not created yet;");
    OsiData* d = static_cast<OsiData*>(solverdata);
    
    vars.sol.clear();
	vars.sol.resize(vars.n,0);
	vars.rc.clear();
	vars.rc.resize(vars.n,0);
    
    if(!hassolution)
		return 1;
    
    const double* sol;
    const double* act;
    
    sol = d->model->getColSolution();
    act = d->model->getRowActivity();
    
    if(boolParam.count("mip") == 0 || !boolParam["mip"])
    {
        const double* dual = d->model->getRowPrice();
        const double* rc = d->model->getReducedCost();
        for (unsigned int i = 0; i < vars.n; i++)
            vars.rc[i] = rc[i];
        for(long i = 0; i < long(nr); i++)
            consts[i].dual = dual[i];
    }
    
    for (unsigned int i = 0; i < vars.n; i++)
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
    
    if(boolParam.count("print_version") > 0 && boolParam["print_version"])
        printf("*********** Genmodel version = %s ***********\n", version.c_str());
    
    //objval = (dynamic_cast<OsiClpSolverInterface*>(d->model))->getModelPtr()->objectiveValue();
    
    delete[] sol;
    delete[] act;

	return 0;
}

long GenModelOsi::CreateModel()
{
    if(!binit)
        return ThrowError("CreateModel() not available : Problem not initialized yet");
    
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
    d->nq = 0;
	memset(d->mat_beg, 0, (nc+1)*sizeof(CoinBigIndex));
    
    d->mult = 1.0;
    if(boolParam.count("maximize") > 0 && boolParam["maximize"])
        d->mult = -1.0;
    
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
    long oldnz = nz;
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
		d->obj[i] = d->mult*vars.obj[i];
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
    
    /*if(boolParam.count("maximize") > 0 && boolParam["maximize"])
    {
		d->model->setObjSense(-1.0);
        d->model->writeMps("tmp/test", "mps", -1.0);
    }
    else
    {*/
        d->model->setObjSense(1.0);
    //    d->model->writeMps("tmp/test", "mps", 1.0);
    //}
    
    /*d->nq = 0;
    if(boolParam.count("qp_mat") > 0 && boolParam["qp_mat"] && !vars.qi.empty())
    {
        boolParam["qp"] = true;
        vector<long>::iterator iti;
        vector<long>::iterator itj = vars.qj.begin();
        vector<double>::iterator itv = vars.qobj.begin();
        vector<vector<pair<int,double> > > tvect;
        for(iti = vars.qi.begin(); iti != vars.qi.end(); iti++, itj++, itv++)
        {
            if(*iti )
            tvect[*iti]
            d->Q_beg[i]++;
            d->Q_r[d->nq] = *iti;
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
        
        printf("Adding quadratic obj\n");
        (dynamic_cast<GenModelOsiInterface*>(d->model))->setQuadraticObjective(d->nc,d->Q_beg,d->Q_r,d->Q_v);
        //(dynamic_cast<OsiClpSolverInterface*>(d->model))->getModelPtr()->loadQuadraticObjective(d->nc,d->Q_beg,d->Q_r,d->Q_v);
    }*/
    
    if(boolParam.count("qp_mat") == 0 || boolParam["qp_mat"])
    {
        vector<vector<pair<int,double> > > qptemp;
        qptemp.resize(nc);
        int* qpbeg = NULL;
        CoinBigIndex* qpind = NULL;
        double* qpv = NULL;
        int qpnz = 0;
        
        vector<long>::iterator iti;
        vector<long>::iterator itj = vars.qj.begin();
        vector<double>::iterator itv = vars.qobj.begin();
        
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
            boolParam["qp"] = true;
            qpbeg = new int[nc+1];
            
            qpv = new double[qpnz];
            qpind = new CoinBigIndex[qpnz];
            
            qpnz=0;
            for(int i = 0; i < int(nc); i++)
            {
                qpbeg[i] = qpnz;
                for(int j = 0; j < int(qptemp[i].size()); j++)
                {
                    qpind[qpnz] = qptemp[i][j].first;
                    qpv[qpnz] = d->mult*2*qptemp[i][j].second;
                    qpnz++;
                }
            }
            qpbeg[nc] = qpnz;
            printf("Adding quadratic obj\n");
            (dynamic_cast<GenModelOsiInterface*>(d->model))->setQuadraticObjective(nc,qpbeg,qpind,qpv);
            //(dynamic_cast<OsiClpSolverInterface*>(d->model))->getModelPtr()->loadQuadraticObjective(nc,qpbeg,qpind,qpv);
            d->nq = qpnz;
            
            delete[] qpbeg;
            delete[] qpind;
            delete[] qpv;
        }
	}
    nz = oldnz;
    
    bcreated = true;
    
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
    bcreated = true;
    
    return 0;
}

long GenModelOsi::Init(string name) //, int type)
{
    try {
    int type = 0;
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
    {
        if(boolParam.count("qp_mat") > 0 && boolParam["qp_mat"] && !vars.qi.empty())
            d->model = new GenModelOsiInterface();
        else
            d->model = new OsiClpSolverInterface();
    }
    /*else if(type==1)
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
        d->model = new OsiClpSolverInterface();*/

    d->model->setStrParam(OsiProbName, name);
    

    SetParam("log_file", 0, "str", "Failure to set the log file", false);
    
    // General settings
    SetParam("log_output_stdout", 0, "bool", "Failure to turn on/off log output to stdout", false);
    SetParam("log_level", 0, "long", "Failure to set log level", false);
    SetParam("use_data_checking", 0, "bool", "Failure to turn on/off data checking", false);
    SetParam("nb_threads", 0, "long", "Failure to set the number of threads", false);
    SetParam("use_preprocessor", 0, "bool", "Failure to use preprocessor", false);
    
    // MIP settings
    SetParam("nb_cut_pass", 0, "long", "Failure to set the number of cut pass", false);
    SetParam("feasibility_pump_level", 0, "long", "Failure to set the feasibility pump level", false);
    SetParam("probing_level", 0, "long", "Failure to set the probing level", false);
    SetParam("mip_emphasis", 0, "long", "Failure to set the MIP emphasis", false);
    SetParam("use_cut_callback", 0, "bool", "Failure to use preprocessor", false);
    
    // Tolerance and limits
    
    if(longParam.count("max_iteration_limit") ==  0) longParam["max_iteration_limit"] = 99999999;
        SetParam("max_iteration_limit", OsiMaxNumIteration, "long", "Failure to set the maximal number of simplex iterations");
    if(dblParam.count("bounds_feasibility_tolerance") ==  0) dblParam["bounds_feasibility_tolerance"] = 1e-9;
    SetParam("bounds_feasibility_tolerance", OsiPrimalTolerance, "dbl", "Failure to set bounds feasibility tolerance");
    SetParam("bounds_feasibility_tolerance", OsiDualTolerance, "dbl", "Failure to set bounds feasibility tolerance");
    SetParam("optimality_tolerance", 0, "dbl", "Failure to set optimality tolerance", false);
    SetParam("markowitz_tolerance", 0, "dbl", "Failure to set Markowitz tolerance", false);
    SetParam("absolute_mip_gap_tolerance", 0, "dbl", "Failure to set absolute gap tolerance", false);
    SetParam("relative_mip_gap_tolerance", 0, "dbl", "Failure to set relative gap tolerance", false);
    SetParam("lp_objective_limit", OsiPrimalObjectiveLimit, "dbl", "Failure to set lp objective limit");
    SetParam("lp_objective_limit", OsiDualObjectiveLimit, "dbl", "Failure to set lp objective limit");
    
    // int params
    //ClpMaxNumIterationHotStart : The maximum number of iterations Clp can execute in hotstart before terminating
        
    // double params
    // ClpDualObjectiveLimit : Set Dual objective limit. This is to be used as a termination criteria in methods where the dual objective monotonically changes (dual simplex)
    // ClpObjOffset : Objective function constant. This the value of the constant term in the objective function
    
    // OsiHintParam {
    //    OsiDoPresolveInInitial = 0, OsiDoDualInInitial, OsiDoPresolveInResolve, OsiDoDualInResolve,
    //    OsiDoScale, OsiDoCrash, OsiDoReducePrint, OsiDoInBranchAndCut,
    //    OsiLastHintParam
    //}
    // OsiHintStrength { OsiHintIgnore = 0, OsiHintTry, OsiHintDo, OsiForceDo }
    
    } catch (string e) {
        printf("Error : %s\n", e.c_str());
    }

    binit = true;
    
	return 0;
}

long GenModelOsi::SetDirectParam(int whichparam, genmodel_param value, string type, string message)
{
    bool status = true;
    if(type == "dbl")
        status = (static_cast<OsiData*>(solverdata))->model->setDblParam(OsiDblParam(whichparam), value.dblval);
    else if(type == "long")
        status = (static_cast<OsiData*>(solverdata))->model->setIntParam(OsiIntParam(whichparam), value.longval);
    else if(type == "str")
        status = (static_cast<OsiData*>(solverdata))->model->setStrParam(OsiStrParam(whichparam), value.strval);
    if ( !status )
        return ThrowError(message);
    return 0;
}

long GenModelOsi::SetParam(string param, int whichparam, string type, string message, bool implemented)
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
                SetDirectParam(whichparam, long2param(0), "long", message);
            else
                SetDirectParam(whichparam, long2param(1), "long", message);
        }
        else if(notimplmessage && !implemented && boolParam.count(param) > 0)
            throw (string("Parameter ")+param+" not implemented in GenModelOsi");
    }
    
    return 0;
}

long GenModelOsi::AddSolverRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name)
{
    if(!bcreated)
        throw string("WriteSolutionToFile() not available : Problem not created yet;");
    
	AddModelRow(ind, val, rhs, sense, name);
	AddCut(&ind[0], &val[0], int(ind.size()), rhs, sense, name.c_str());

	return 0;
}

long GenModelOsi::AddCut(int* cols, double* vals, int nz, double rhs, char sense, const char* name)
{
    if(!bcreated)
        throw string("WriteSolutionToFile() not available : Problem not created yet;");
    
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
    if(!bcreated)
        throw string("WriteSolutionToFile() not available : Problem not created yet;");
    
	AddModelCol(ind, val, obj, lb, ub, name, type);
	AddCol(&ind[0], &val[0], int(ind.size()), obj, lb, ub, name.c_str(), type);

	return 0;
}


long GenModelOsi::AddCol(int* newi, double* newcol, int nz, double obj, double lb, double ub, const char* name, char type)
{
    if(!bcreated)
        throw string("WriteSolutionToFile() not available : Problem not created yet;");
    
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
    mipmodel = NULL;
	mat_beg = NULL;
    /*Q_beg = NULL;
	Q_r = NULL;
	Q_v = NULL;*/
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
