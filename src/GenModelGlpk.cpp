#include "GenModelGlpk.h"
#ifdef OSI_MODULE
#include "ProblemReaderOsi.h"
#endif
#include <limits>

long GenModelGlpk::Solve()
{
	if(boolParam.count("qp") > 0 && boolParam["qp"])
		return 1;
	GlpkData* d = static_cast<GlpkData*>(solverdata);

	if(boolParam.count("mip") > 0 && boolParam["mip"])
	{
		d->mip_param.presolve = GLP_ON;			// Set the presolve to GLP_ON which is disabled by default
		glp_intopt(d->model, &(d->mip_param));
	}
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
		glp_simplex(d->model, &(d->simplex_param));

	return 0;
}

long GenModelGlpk::SetSol()
{
	if(boolParam.count("qp") > 0 && boolParam["qp"])
		return 1;
	vars.sol.clear();
	vars.sol.resize(vars.n,0);
	vars.rc.clear();
	vars.rc.resize(vars.n,0);
	GlpkData* d = (GlpkData*)solverdata;

	if(boolParam.count("mip") > 0 && boolParam["mip"])
	{
		solstat = glp_mip_status(d->model);
		objval = glp_mip_obj_val(d->model);
		for(unsigned long i = 0; i < nc; i++)
			vars.sol[i] = glp_mip_col_val(d->model, i+1);
		for(unsigned long i = 0; i < nr; i++)
			consts[i].dual = glp_mip_row_val(d->model, i+1);
	}
	else if(strParam.count("algo") > 0 && strParam["algo"] == "interior")
	{
		solstat = glp_ipt_status(d->model);
		objval = glp_ipt_obj_val(d->model);
		for(unsigned long i = 0; i < nc; i++)
		{
			vars.sol[i] = glp_ipt_col_prim(d->model, i+1);
			vars.rc[i] = glp_ipt_col_dual(d->model, i+1);
		}
		for(unsigned long i = 0; i < nr; i++)
		{
			consts[i].dual = glp_ipt_row_prim(d->model, i+1);
			consts[i].slack = glp_ipt_row_dual(d->model, i+1);
		}
	}
	else
	{
		solstat = glp_get_status(d->model);
		objval = glp_get_obj_val(d->model);
		for(unsigned long i = 0; i < nc; i++)
		{
			vars.sol[i] = glp_get_col_prim(d->model, i+1);
			vars.rc[i] = glp_get_col_dual(d->model, i+1);
		}
		for(unsigned long i = 0; i < nr; i++)
		{
			//consts[i].dual = glp_get_row_prim(d->model, i+1);
			//consts[i].slack = glp_get_row_dual(d->model, i+1);
			consts[i].slack = consts[i].lrhs-glp_get_row_prim(d->model, i+1);
			consts[i].dual = glp_get_row_dual(d->model, i+1);
		}
	}
	//if (solstat == GLP_OPT)
	//	solstat = 1;
    if(boolParam.count("print_version") > 0 && boolParam["print_version"])
        printf("*********** Genmodel version = %s ***********\n", version.c_str());

	return 0;
}

long GenModelGlpk::CreateModel()
{
	if(!vars.qi.empty())
	{
		boolParam["qp"] = true;
		printf("No quadratic programming in glpk\n");
		return 1;
	}


	GlpkData* d = (GlpkData*)solverdata;

	d->nc = nc;
	d->nr = nr;

	glp_set_obj_name(d->model, "obj");

	if(boolParam.count("maximize") > 0 && boolParam["maximize"])
		glp_set_obj_dir(d->model, GLP_MAX);
	else
		glp_set_obj_dir(d->model, GLP_MIN);

	d->lrhs = new double[nr];
	d->urhs = new double[nr];
	d->sense = new char[nr];
	d->ub = new double[nc];
	d->lb = new double[nc];
	d->obj = new double[nc];
	d->type = new char[nc];

	d->mat_r = new int[nz+1];
	d->mat_c = new int[nz+1];
	d->mat_v = new double[nz+1];
	d->cname = new char*[nc];
	d->rname = new char*[nr];


	nz=0;
	for(unsigned long i = 0; i < nr; i++)
	{
		d->rname[i] = new char[consts[i].name.length()+1];
		snprintf(d->rname[i], consts[i].name.length()+1, "%s", consts[i].name.c_str());

		for(unsigned long j = 0; j < consts[i].nz; j++)
		{
			d->mat_r[nz+1] = i+1;
			d->mat_c[nz+1] = consts[i].cols[j]+1;
			d->mat_v[nz+1] = consts[i].coefs[j];
			nz++;
		}

		d->lrhs[i] = consts[i].lrhs;
		d->urhs[i] = consts[i].urhs;
		d->sense[i] = consts[i].sense;

	}
	for(unsigned long i = 0; i < nc; i++)
	{
		d->cname[i] = new char[vars.name[i].length()+1];
		snprintf(d->cname[i], vars.name[i].length()+1, "%s", vars.name[i].c_str());
		d->obj[i] = vars.obj[i];
		d->ub[i] = vars.ub[i];
		d->lb[i] = vars.lb[i];
		d->type[i] = vars.type[i];

	}

	glp_add_rows(d->model, nr);
	for(unsigned long i = 0; i < nr; i++)
	{
		glp_set_row_name(d->model, i+1, d->rname[i]);
		if(d->sense[i] == 'E')
			glp_set_row_bnds(d->model, i+1, GLP_FX, d->lrhs[i], d->urhs[i]);
		else if(d->sense[i] == 'G')
			glp_set_row_bnds(d->model, i+1, GLP_LO, d->lrhs[i], d->urhs[i]);
		else if(d->sense[i] == 'L')
			glp_set_row_bnds(d->model, i+1, GLP_UP, d->urhs[i], d->lrhs[i]);
		else if(d->sense[i] == 'R')
			glp_set_row_bnds(d->model, i+1, GLP_DB, d->lrhs[i], d->urhs[i]);
	}

	glp_add_cols(d->model, nc);
	for(unsigned long i = 0; i < nc; i++)
	{
		glp_set_col_name(d->model, i+1, d->cname[i]);
		if(d->lb[i] == -numeric_limits<double>::infinity() && d->ub[i] == numeric_limits<double>::infinity())
			glp_set_col_bnds(d->model, i+1, GLP_FR, d->lb[i], d->ub[i]);
		else if(d->lb[i] == -numeric_limits<double>::infinity())
			glp_set_col_bnds(d->model, i+1, GLP_UP, d->lb[i], d->ub[i]);
		else if(d->ub[i] == numeric_limits<double>::infinity())
			glp_set_col_bnds(d->model, i+1, GLP_LO, d->lb[i], d->ub[i]);
		else if(d->lb[i] == d->ub[i])
			glp_set_col_bnds(d->model, i+1, GLP_FX, d->lb[i], d->ub[i]);
		else
			glp_set_col_bnds(d->model, i+1, GLP_DB, d->lb[i], d->ub[i]);
		glp_set_obj_coef(d->model, i+1, d->obj[i]);

		if(d->type[i] == 'I')
			glp_set_col_kind(d->model, i+1, GLP_IV);
		else if(d->type[i] == 'B')
			glp_set_col_kind(d->model, i+1, GLP_BV);
		else
			glp_set_col_kind(d->model, i+1, GLP_CV);
	}

	glp_load_matrix(d->model, nz, d->mat_r,	d->mat_c, d->mat_v);

	return 0;
}

long GenModelGlpk::CreateModel(string filename, int type, string dn)
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

long GenModelGlpk::AddSolverRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name)
{
	AddModelRow(ind, val, rhs, sense, name);
	AddCut(&ind[0], &val[0], int(ind.size()), rhs, sense, name.c_str());

	return 0;
}

long GenModelGlpk::AddCut(int* cols, double* vals, int nz, double rhs, char sense, const char* name)
{
	GlpkData* d = (GlpkData*)solverdata;

	glp_add_rows(d->model, 1);
	for(unsigned long i = 0; i < nr; i++)
	{
		glp_set_row_name(d->model, d->nr+1, name);
		if(sense == 'E')
			glp_set_row_bnds(d->model, d->nr+1, GLP_FX, rhs, rhs);
		else if(d->sense[i] == 'G')
			glp_set_row_bnds(d->model, d->nr+1, GLP_LO, rhs, rhs);
		else if(d->sense[i] == 'L')
			glp_set_row_bnds(d->model, d->nr+1, GLP_UP, rhs, rhs);
	}

	for(long i = 0; i < nz; i++)
		cols[i]+=1;


	glp_set_mat_row(d->model, d->nr+1, nz, cols, vals);
	d->nr++;

	for(long i = 0; i < nz; i++)
		cols[i]-=1;

	return 0;
}

long GenModelGlpk::AddSolverCol(vector<int>& ind, vector<double>& val, double obj, double lb, double ub, string name, char type)
{
	AddModelCol(ind, val, obj, lb, ub, name, type);
	AddCol(&ind[0], &val[0], int(ind.size()), obj, lb, ub, name.c_str(), type);

	return 0;
}


long GenModelGlpk::AddCol(int* newi, double* newcol, int nz, double obj, double lb, double ub, const char* name, char type)
{
	GlpkData* d = (GlpkData*)solverdata;

	glp_add_cols(d->model, 1);

	glp_set_col_name(d->model, d->nc+1, name);
	if(lb == -numeric_limits<double>::infinity() && ub == numeric_limits<double>::infinity())
		glp_set_col_bnds(d->model, d->nc+1, GLP_FR, lb, ub);
	else if(lb == -numeric_limits<double>::infinity())
		glp_set_col_bnds(d->model, d->nc+1, GLP_UP, lb, ub);
	else if(ub == numeric_limits<double>::infinity())
		glp_set_col_bnds(d->model, d->nc+1, GLP_LO, lb, ub);
	else if(lb == ub)
		glp_set_col_bnds(d->model, d->nc+1, GLP_FX, lb, ub);
	else
		glp_set_col_bnds(d->model, d->nc+1, GLP_DB, lb, ub);

	glp_set_obj_coef(d->model, d->nc+1, obj);

	if(type == 'I')
		glp_set_col_kind(d->model, d->nc+1, GLP_IV);
	else if(type == 'B')
		glp_set_col_kind(d->model, d->nc+1, GLP_BV);
	else
		glp_set_col_kind(d->model, d->nc+1, GLP_CV);

	int* vnewi = new int[nz+1];
	double* vnewcol = new double[nz+1];
	vnewi[0] = 0;
	vnewcol[0] = 0.0;

	for(long i = 0; i < nz; i++)
	{
		vnewi[i+1] = newi[i]+1;
		vnewcol[i+1] = newcol[i];
	}

	glp_set_mat_col(d->model, d->nc+1, nz, vnewi, vnewcol);
	d->nc++;

	delete[] vnewi;
	delete[] vnewcol;

	return 0;
}

long GenModelGlpk::Init(string name)
{
	if(solverdata == NULL)
		solverdata = new GlpkData();
	else
	{
		static_cast<GlpkData*>(solverdata)->Delete();
		static_cast<GlpkData*>(solverdata)->Reset();
	}

	GlpkData* d = static_cast<GlpkData*>(solverdata);

	d->model = glp_create_prob();
	glp_set_prob_name(d->model, name.c_str());
	glp_init_iptcp(&(d->interior_param));
	glp_init_smcp(&(d->simplex_param));
	glp_init_iocp(&(d->mip_param));

	if(boolParam.count("log_output_stdout") > 0 && !boolParam["log_output_stdout"])
	{
		d->simplex_param.msg_lev = GLP_MSG_OFF;
		d->mip_param.msg_lev = GLP_MSG_OFF;
	}
	if(dblParam.count("time_limit") > 0)
	{
		d->simplex_param.tm_lim = dblParam["time_limit"] * 1000;
		d->mip_param.tm_lim = dblParam["time_limit"] * 1000;
	}
    if(dblParam.count("relative_mip_gap_tolerance") > 0)
	{
		//d->simplex_param.tm_lim = dblParam["relative_mip_gap_tolerance"];
		d->mip_param.mip_gap = dblParam["relative_mip_gap_tolerance"];
	}
    
    
	return 0;
}

long GenModelGlpk::Clean()
{
	if(solverdata != NULL)
		delete static_cast<GlpkData*>(solverdata);

	return 0;
}

long GlpkData::Reset()
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
	cname = NULL;
	rname = NULL;
	model = NULL;

	return 0;
}

GlpkData::GlpkData()
{
	Reset();
}

GlpkData::~GlpkData()
{
	Delete();
}
long GlpkData::Delete()
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
	if(model != NULL)
		glp_delete_prob(model);

	delete[] rname;

	return 0;
}
