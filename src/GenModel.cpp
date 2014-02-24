#include "GenModel.h"
//#include "SqlCaller.h"
#include <math.h>
#include <time.h>

GenModel::GenModel()
{
    version = "genmodelgrb-0.0.7 build 0001";
    hassolution = false;
    bcreated = false;
    binit = false;
	nc=0;
	nr=0;
	solverdata = NULL;
}


double GenModel::FindConstraintMaxLhs(long row)
{
    double total = 0.0;
    for(int i = 0; i < consts[row].nz; i++)
        total += (consts[row].coefs[i] >= 0 ? vars.ub[consts[row].cols[i]] : vars.lb[consts[row].cols[i]])*consts[row].coefs[i];
    return total;
}

double GenModel::FindConstraintMinLhs(long row)
{
    double total = 0.0;
    for(int i = 0; i < consts[row].nz; i++)
        total += (consts[row].coefs[i] >= 0 ? vars.lb[consts[row].cols[i]] : vars.ub[consts[row].cols[i]])*consts[row].coefs[i];
    return total;
}

long GenModel::MakeConstraintFeasible(long row)
{
    if(consts[row].sense == 'L')
    {
        double min = FindConstraintMinLhs(row);
        if(min > consts[row].lrhs)
            consts[row].lrhs = min;
    }
    else if(consts[row].sense == 'G')
    {
        double max = FindConstraintMaxLhs(row);
        if(max < consts[row].lrhs)
            consts[row].lrhs = max;
    }
    else if(consts[row].sense == 'R')
    {
        double min = FindConstraintMinLhs(row);
        double max = FindConstraintMaxLhs(row);
        if(max < consts[row].lrhs)
            consts[row].lrhs = max;
        else if(min > consts[row].urhs)
            consts[row].urhs = min;
    }
    return 0;
}

long GenModel::SetLongParam(string param, long val)
{
	longParam[param] = val;
	return 0;
}

long GenModel::SetDblParam(string param, double val)
{
	dblParam[param] = val;
	return 0;
}

long GenModel::SetBoolParam(string param, bool val)
{
	boolParam[param] = val;
	return 0;
}

long GenModel::SetStrParam(string param, string val)
{
	strParam[param] = val;
	return 0;
}

long GenModel::ThrowError(string error)
{
    printf("%s\n", error.c_str());
    throw error;
}

long GenModel::PrintObjVal()
{
	printf("obj :		%f\n", objval);

	return 0;
}

long GenModel::PrintSol()
{
	//printf("obj :		%f\n", objval);
	for(long i = 0; i < long(vars.n); i++)
	{
		printf("%s :		%f\n", vars.name[i].c_str(), vars.sol[i]);
	}
	return 0;
}

long GenModel::PrintSol(string v)
{
	//printf("obj :		%f\n", objval);
	if(vars.offset.count(v) == 0)
		return 0;
	map<string, long>::iterator it = vars.offset.find(v);
	long deb = it->second;
	long fin = vars.n;
	for(it = vars.offset.begin(); it != vars.offset.end(); it++)
	{
		if(it->second > deb && it->second < fin)
			fin = it->second;
	}
	for(long i = deb; i < fin; i++)
	{
		printf("%s :		%f\n", vars.name[i].c_str(), vars.sol[i]);
	}
	return 0;
}

long GenModel::PrintSolNz()
{
	//printf("obj :		%f\n", objval);
	for(long i = 0; i < long(vars.n); i++)
	{
		if(fabs(vars.sol[i]) > 0.000001)
			printf("%s :		%f\n", vars.name[i].c_str(), vars.sol[i]);
	}
	return 0;
}

long GenModel::PrintSolNz(string v)
{
	//printf("obj :		%f\n", objval);
	if(vars.offset.count(v) == 0)
		return 0;
	map<string, long>::iterator it = vars.offset.find(v);
	long deb = it->second;
	long fin = vars.n;
	for(it = vars.offset.begin(); it != vars.offset.end(); it++)
	{
		if(it->second > deb && it->second < fin)
			fin = it->second;
	}
	for(long i = deb; i < fin; i++)
	{
		if(fabs(vars.sol[i]) > 0.000001)
			printf("%s :		%f\n", vars.name[i].c_str(), vars.sol[i]);
	}
	return 0;
}

long ModConsts::AddNz(long c, double v)
{
	cols.push_back(c);
	coefs.push_back(v);
	nz++;
	return 0;
}

long GenModel::AddConst(string cname)
{
	ci[cname] = nc;
	consts.push_back(ModConsts());
	consts.back().name = cname;
	nc++;
	return 0;
}

long GenModel::AddConst(string cname, double rhs, char sense)
{
	AddConst(cname);
	consts.back().lrhs = rhs;
	consts.back().sense = sense;
	return 0;
}

long GenModel::AddNz(long row, long col, double val)
{
	consts[row].AddNz(col, val);

	return 0;
}
long GenModel::AddNzToLast(long col, double val)
{
	consts.back().AddNz(col, val);
	return 0;
}

long GenModel::SetNumbers()
{
	nr = consts.size();
	nc = vars.n; //vars.obj.size();
	nz = 0;
	for(long i = 0; i < long(nr); i++)
	{
		consts[i].nz = consts[i].cols.size();
		nz+=consts[i].nz;
	}
	return 0;
}

long GenModel::AddVar(string nn, double o, double l, double u, char t)
{
	vars.AddVar(nn, o, l, u, t);
	return 0;
}

long GenModel::AddVars(string nn, long size, double o, double l, double u, char t)
{
	vars.AddVars(nn, size, o, l, u, t);
	return 0;
}

long GenModel::SetQpCoef(long i, long j, double val)
{
	vars.SetQpCoef(i, j, val);
	return 0;
}

long GenModel::AddModelCol(vector<int>& ind, vector<double>& val, double obj, double lb, double ub, string name, char type)
{
	AddVar(name, obj, lb, ub, type);
	for(long i = 0; i < long(ind.size()); i++)
		AddNz(ind[i], vars.n-1, val[i]);

	return 0;
}

long GenModel::AddModelRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name)
{
	AddConst(name, rhs, sense);
	for(long i = 0; i < long(ind.size()); i++)
		AddNzToLast(ind[i], val[i]);

	return 0;
}

long GenModel::ChangeBulkBounds(int count, int * ind, char * type, double * vals)
{
	return 0;
}

long GenModel::ChangeBulkObjectives(int count, int * ind, double * vals)
{
	return 0;
}

long GenModel::WriteProblemToLpFile(string filename)
{
    throw string("WriteProblemToLpFile() Not implemented");
	return 0;
}

long GenModel::WriteSolutionToFile(string filename)
{
    throw string("WriteSolutionToFile() Not implemented");
	return 0;
}

long GenModel::DeleteMipStarts()
{
    throw string("DeleteMipStarts() Not implemented");
	return 0;
}

double GenModel::GetMIPRelativeGap()
{
	throw string("DeleteMipStarts() Not implemented");
    return 0.0;
}

long ModVars::AddVar(string nn, double o, double l, double u, char t)
{
	//offset[nn] = n;
	name.push_back(nn);
	obj.push_back(o);
	ub.push_back(u);
	lb.push_back(l);
	type.push_back(t);
	n++;

	return 0;
}

long ModVars::AddVars(string nn, long size, double o, double l, double u, char t)
{
	offset[nn] = n;
	char c[4096];
	for(long i = 0; i < size; i++)
	{
		snprintf(c, 4096, "%s_%ld", nn.c_str(), i);
		AddVar(string(c), o, l, u, t);
	}

	return 0;
}

long ModVars::SetQpCoef(long i, long j, double val)
{
	qi.push_back(i);
	qj.push_back(j);
	qobj.push_back(val);

	return 0;
}

long ModVars::Print()
{
	for(long i = 0; i < long(n); i++)
	{
		printf("%ld: %s obj=%f [%c]\n", i, name[i].c_str(), obj[i], type[i]);
	}

	return 0;
}

template <typename T>
void _freeall(T& t)
{
    T tmp;
    t.swap(tmp);
}

long GenModel::ClearStructure()
{
    consts.clear();
    _freeall(consts);
    ci.clear();
    _freeall(ci);
    vars.name.clear();
    _freeall(vars.name);
	vars.obj.clear();
    _freeall(vars.obj);
	vars.type.clear();
    _freeall(vars.type);
	vars.offset.clear();
    _freeall(vars.offset);
	vars.ub.clear();
    _freeall(vars.ub);
	vars.lb.clear();
	_freeall(vars.lb);
    vars.sol.clear();
    _freeall(vars.sol);
	vars.rc.clear();
    _freeall(vars.rc);
	vars.qobj.clear();
    _freeall(vars.qobj);
	vars.qi.clear();
    _freeall(vars.qi);
	vars.qj.clear();
    _freeall(vars.qj);
    
    return 0;
}


genmodel_param dbl2param(double val)
{
    genmodel_param ret;
    ret.dblval = val;
    return ret;
}

genmodel_param long2param(long val)
{
    genmodel_param ret;
    ret.longval = val;
    return ret;
}

genmodel_param str2param(string val)
{
    genmodel_param ret;
    ret.strval = val.c_str();
    return ret;
}
