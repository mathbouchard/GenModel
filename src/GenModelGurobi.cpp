#include "GenModelGurobi.h"
//#include "ProblemReader.h"
#include <limits>

long GenModelGurobi::Solve()
{
	GurobiData* d = static_cast<GurobiData*>(solverdata);

	try
	{
		d->model->optimize();
	}
	catch (GRBException e)
	{
		printf("Error code = %d\n", e.getErrorCode());
		printf("%s\n", e.getMessage().c_str());
	}
	catch (...)
	{
		printf("Exception during optimization\n");
	}


/*	if(boolParam.count("mip") > 0 && boolParam["mip"])
		glp_intopt(d->model, &(d->mip_param));
	else if(strParam.count("algo") > 0 && strParam["algo"] == "interior")
		glp_interior(d->model, &(d->interior_param));
	else if(strParam.count("algo") > 0 && strParam["algo"] == "dual")
	{
		d->simplex_param.meth = GLP_DUAL;
		glp_simplex(d->model, &(d->simplex_param));
	}
	else if(strParam.count("algo") > 0 && strParam["algo"] == "primal")
	{
		d->simplex_param.meth = GLP_PRIMAL;
		glp_simplex(d->model, &(d->simplex_param));
	}
	else
		glp_simplex(d->model, &(d->simplex_param));*/

	return 0;
}

long GenModelGurobi::SetSol()
{
	vars.sol.clear();
	vars.sol.resize(vars.n,0);
	vars.rc.clear();
	vars.rc.resize(vars.n,0);
	GurobiData* d = (GurobiData*)solverdata;

	try
	{
		solstat = d->model->get(GRB_IntAttr_Status);
		objval = d->model->get(GRB_DoubleAttr_ObjVal);

		GRBConstr*	row = d->model->getConstrs();

		for(unsigned long i = 0; i < (unsigned long)(d->nc); i++)
		{
			vars.sol[i] = d->model->getVar(i).get(GRB_DoubleAttr_X);
			//vars.sol[i] = d->grb_v[i].get(GRB_DoubleAttr_X);
			if(boolParam.count("mip") == 0 || !boolParam["mip"])
				vars.rc[i] = d->model->getVar(i).get(GRB_DoubleAttr_RC);
			else
				vars.rc[i] = 0.0;
			//vars.rc[i] = d->grb_v[i].get(GRB_DoubleAttr_RC);
		}
		for(unsigned long i = 0; i < (unsigned long)(d->nr); i++)
		{
			if(d->equiv[i] >= 0)
			{
				//consts[d->equiv[i]].dual = d->model->getConstr(i).get(GRB_DoubleAttr_Pi);
				if(boolParam.count("mip") == 0 || !boolParam["mip"])
					consts[d->equiv[i]-1].dual = row[i].get(GRB_DoubleAttr_Pi);
				else
					consts[d->equiv[i]-1].dual = 0.0;
				//consts[d->equiv[i]].slack = d->model->getConstr(i).get(GRB_DoubleAttr_Slack);
				consts[d->equiv[i]-1].slack = row[i].get(GRB_DoubleAttr_Slack);
			}
			else
			{
				consts[-d->equiv[i]-1].dual = 0.0;
				consts[-d->equiv[i]-1].slack = 0.0;
			}

			//consts[-d->equiv[i]].dual += d->model->getConstr(i).get(GRB_DoubleAttr_Pi);
			if(boolParam.count("mip") == 0 || !boolParam["mip"])
				consts[-d->equiv[i]-1].dual += row[i].get(GRB_DoubleAttr_Pi);
			//consts[-d->equiv[i]].slack += d->model->getConstr(i).get(GRB_DoubleAttr_Slack);
			consts[-d->equiv[i]-1].slack += row[i].get(GRB_DoubleAttr_Slack);
		}

		//if (solstat == GRB_OPTIMAL)
		//	solstat = 1;
	}
	catch (GRBException e)
	{
		printf("Error code = %d\n", e.getErrorCode());
		printf("%s\n", e.getMessage().c_str());
		exit(0);
	}
	catch (...)
	{
		printf("Exception during optimization\n");
	}

	return 0;
}

long GenModelGurobi::CreateModel()
{
	GurobiData* d = (GurobiData*)solverdata;

	try
	{
		if(boolParam.count("maximize") > 0 && boolParam["maximize"])
			d->model->set(GRB_IntAttr_ModelSense, -1);
		else
			d->model->set(GRB_IntAttr_ModelSense, 1);

		d->model->update();


		d->nc = nc;
		d->nr = nr;

		d->ub = new double[nc];
		d->lb = new double[nc];
		d->obj = new double[nc];
		d->type = new char[nc];
		d->cname = new string[nc];
		d->equiv = new long[2*nr];

		nz=0;

		for(unsigned long i = 0; i < nc; i++)
		{
			d->cname[i] = vars.name[i];
			d->obj[i] = vars.obj[i];
			d->ub[i] = vars.ub[i];
			d->lb[i] = vars.lb[i];
			if(vars.type[i] == 'I' && boolParam["mip"])
				d->type[i] = GRB_INTEGER;
			else if(vars.type[i] == 'B' && boolParam["mip"])
				d->type[i] = GRB_BINARY;
			else
				d->type[i] = GRB_CONTINUOUS;
		}
		d->grb_v = d->model->addVars(d->lb, d->ub, d->obj, d->type, d->cname, nc);
		d->model->update();

		long eqi = 0;
		for(long i = 0; i < long(nr); i++, eqi++)
		{
			GRBLinExpr lhs = 0;
			for(unsigned long j = 0; j < consts[i].nz; j++)
			{
				lhs += consts[i].coefs[j]*d->grb_v[consts[i].cols[j]];
				nz++;
			}
			if(consts[i].sense == 'L')
				d->model->addConstr(lhs, GRB_LESS_EQUAL, consts[i].lrhs, consts[i].name);
			else if(consts[i].sense == 'G')
				d->model->addConstr(lhs, GRB_GREATER_EQUAL, consts[i].lrhs, consts[i].name);
			else if(consts[i].sense == 'R')
			{
				d->model->addConstr(lhs, GRB_GREATER_EQUAL, consts[i].lrhs, consts[i].name+"_lb");
				d->model->addConstr(lhs, GRB_LESS_EQUAL, consts[i].urhs, consts[i].name+"_ub");
				d->nr++;
				d->equiv[eqi] = (i+1);
				eqi++;
			}
			else
				d->model->addConstr(lhs, GRB_EQUAL, consts[i].lrhs, consts[i].name);
			d->equiv[eqi] = -(i+1);
		}
		vector<long>::iterator iti;
		vector<long>::iterator itj = vars.qj.begin();
		vector<double>::iterator itv = vars.qobj.begin();

		GRBQuadExpr quad = 0;
		for(iti = vars.qi.begin(); iti != vars.qi.end(); iti++, itj++, itv++)
		{
			quad += 0.5*(*itv)*d->grb_v[*iti]*d->grb_v[*itj];
		}
		if(!vars.qi.empty())
		{
			printf("Quadratic obj !! \n");
			boolParam["qp"] = true;
			d->model->setObjective(quad);
			for(unsigned long i = 0; i < nc; i++)
				d->grb_v[i].set(GRB_DoubleAttr_Obj, vars.obj[i]);
		}
		d->model->update();

	}
	catch (GRBException e)
	{
		printf("Error code = %d\n", e.getErrorCode());
		printf("%s\n", e.getMessage().c_str());
	}
	catch (...)
	{
		printf("Exception during optimization\n");
	}

	return 0;
}

long GenModelGurobi::CreateModel(string filename, int type, string dn)
{
    //ReadFromFile(static_cast<GenModel*>(this), filename, type);
    SetNumbers();
    CreateModel();
    
    return 0;
}

long GenModelGurobi::AddSolverRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name)
{
	AddModelRow(ind, val, rhs, sense, name);
	AddCut(&ind[0], &val[0], int(ind.size()), rhs, sense, name.c_str());

	return 0;
}

long GenModelGurobi::AddCut(int* cols, double* vals, int nz, double rhs, char sense, const char* name)
{
	try
	{
	GurobiData* d = (GurobiData*)solverdata;

	char grbsens = GRB_EQUAL;
	if(sense == 'L')
		grbsens = GRB_LESS_EQUAL;
	else if(sense == 'G')
		grbsens = GRB_GREATER_EQUAL;

	GRBLinExpr lhs = 0;
	for(int i = 0; i < nz; i++)
		lhs += d->model->getVar(i)*vals[i];

	d->model->addConstr(lhs, grbsens, rhs, name);
	d->model->update();

	//GRBaddconstr (d->model,	nz, cols, vals, sense, rhs, name);
	d->nr++;
	}
	catch(GRBException e)
	{
		printf("GRBError: %s\n", e.getMessage().c_str());
		exit(0);
	}
	catch (...)
	{
		printf("Exception during optimization\n");
	}

	return 0;
}

long GenModelGurobi::AddSolverCol(vector<int>& ind, vector<double>& val, double obj, double lb, double ub, string name, char type)
{
	AddModelCol(ind, val, obj, lb, ub, name, type);
	AddCol(&ind[0], &val[0], int(ind.size()), obj, lb, ub, name.c_str(), type);

	return 0;
}


long GenModelGurobi::AddCol(int* newi, double* newcol, int nz, double obj, double lb, double ub, const char* name, char type)
{
	try
	{
	GurobiData* d = (GurobiData*)solverdata;

	char t = GRB_CONTINUOUS;
	if(type == 'I')
		t = GRB_INTEGER;
	else if(type == 'B')
		t = GRB_BINARY;

	GRBConstr* cons = new GRBConstr[nz];

	for(int i = 0; i < nz; i++)
		cons[i] = d->model->getConstr(newi[i]);


	//GRBColumn col;
	//col.addTerms(newcol, cons, nz);
	//d->model->addVar(lb, ub, obj, t, col, string(name));
	d->model->addVar(lb, ub, obj, t, nz, cons, newcol, string(name));
	d->model->update();

	d->nc++;

	delete[] cons;
	}
	catch(GRBException e)
	{
		printf("GRBError: %s\n", e.getMessage().c_str());
		exit(0);
	}
	catch (...)
	{
		printf("Exception during optimization\n");
	}

	return 0;
}

long GenModelGurobi::Init(string name)
{
	try
	{
	if(solverdata == NULL)
		solverdata = new GurobiData();
	else
	{
		static_cast<GurobiData*>(solverdata)->Delete();
		static_cast<GurobiData*>(solverdata)->Reset();
	}

	GurobiData* d = static_cast<GurobiData*>(solverdata);

	d->env = new GRBEnv();
	d->model = new GRBModel(*(d->env));
	d->model->set(GRB_StringAttr_ModelName, name);

	if(boolParam.count("screenoff") > 0 && boolParam["screenoff"])
		d->model->getEnv().set(GRB_IntParam_OutputFlag, 0);

	if(dblParam.count("timelimit") > 0)
		d->model->getEnv().set(GRB_DoubleParam_TimeLimit, dblParam["timelimit"]);

	if(dblParam.count("opttol"))
		d->model->getEnv().set(GRB_DoubleParam_OptimalityTol, dblParam["opttol"]);

	if(dblParam.count("epgap"))
		d->model->getEnv().set(GRB_DoubleParam_MIPGap, dblParam["epgap"]);

	}
	catch(GRBException e)
	{
		printf("GRBError: %s\n", e.getMessage().c_str());
		exit(0);
	}
	catch (...)
	{
		printf("Exception during optimization\n");
	}

	return 0;
}

long GenModelGurobi::Clean()
{
	if(solverdata != NULL)
		delete static_cast<GurobiData*>(solverdata);

	return 0;
}

long GurobiData::Reset()
{
	ub = NULL;
	lb = NULL;
	type = NULL;
	obj = NULL;
	cname = NULL;
	env = NULL;
	grb_v = NULL;
	equiv = NULL;

	return 0;
}

GurobiData::GurobiData()
{
	Reset();
}

GurobiData::~GurobiData()
{
	Delete();
}
long GurobiData::Delete()
{
	if(obj != NULL)
		delete[] obj;
	if(ub != NULL)
		delete[] ub;
	if(lb != NULL)
		delete[] lb;
	if(type != NULL)
		delete[] type;
	if(cname != NULL)
		delete[] cname;
	if(grb_v != NULL)
		delete[] grb_v;
	if(env != NULL)
		delete env;
	if(equiv != NULL)
		delete[] equiv;
	return 0;
}
