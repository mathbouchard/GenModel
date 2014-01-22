// GenModelDLL.cpp : Defines the initialization routines for the DLL.
//

#include "GenModelDLL.h"

#ifdef CPLEX_SOLVER
	#include "GenModelCplex.h"
#endif

#ifdef GUROBI_SOLVER
	#include "GenModelGurobi.h"
#endif

#ifdef HG_SOLVER
	#include "GenModelHG.h"
#endif

#ifdef GLPK_SOLVER
	#include "GenModelGlpk.h"
#endif

#include <vector>
#include <memory.h>

using namespace std;

vector<GenModel*> gmmap;

GENMODEL_EXTERN_C_BEGIN

GenModelDLL_API long STDCALL _AddConst(char* cname, double rhs, char sense, long token)
{
	string cnameAsString(cname);
	return gmmap[token]->AddConst(cnameAsString, rhs, sense);
}

GenModelDLL_API bool STDCALL _AddConstBulk(char* cname, double* rhs, long length, char sense, long token)
{
	string cnameAsString(cname);
	for (int i = 0; i < length; i++)
	{
		gmmap[token]->AddConst(cnameAsString, rhs[i], sense);
	}
	return true;
}

GenModelDLL_API long STDCALL _AddVar(char* nn, double o, double l, double u, char t, long token)
{
	string nnAsString(nn);
	return gmmap[token]->AddVar(nnAsString, o, l, u, t);
}

GenModelDLL_API bool STDCALL _AddVarBulk(char* nn, double* o, long length, double l, double u, char t, long token)
{
	string nnAsString(nn);
	for (int i = 0; i < length; i++)
	{
		gmmap[token]->AddVar(nnAsString, o[i], l, u, t);
	}
	return true;
}

GenModelDLL_API long STDCALL _AddNz(long row, long col, double val, long token)
{
	return gmmap[token]->AddNz(row, col, val);
}

GenModelDLL_API long STDCALL _AddNzToLast(long col, double val, long token)
{
	return gmmap[token]->AddNzToLast(col, val);
}

GenModelDLL_API long STDCALL _AddNzBulk(long* rows, long* cols, double* values, long valuesLength, long rowCount, long colCount, long iterations, long token)
{
	//printf ("AddNzBulk %ld,%ld,%ld,%ld\n",valuesLength, rowCount, colCount, iterations);

	int rowIndex = 0;
	int colIndex = 0;
	for(int iter = 1; iter <= iterations; iter++)
	{
		//printf ("Iter %ld\n",iter);	
		for(int i = 0; i < valuesLength; i++)
		{
			//printf ("i ri ci %ld,%ld,%ld\n",i,rows[i] + rowIndex, cols[i] + colIndex);
			gmmap[token]->AddNz(rows[i] + rowIndex, cols[i] + colIndex, values[i]);
		}
		rowIndex += rowCount;
		colIndex += colCount;
		//printf ("ri ci %ld,%ld\n", rowIndex, colIndex);	
	}
	return true;
}

GenModelDLL_API long STDCALL _SetNumbers(long token)
{
	return gmmap[token]->SetNumbers();
}

GenModelDLL_API long STDCALL _SetLongParam(char* param, long val, long token)
{
	string paramAsString(param);
	return gmmap[token]->SetLongParam(paramAsString, val);
}

GenModelDLL_API long STDCALL _SetDblParam(char* param, double val, long token)
{
	string paramAsString(param);
	return gmmap[token]->SetDblParam(paramAsString, val);
}

GenModelDLL_API long STDCALL _SetBoolParam(char* param, bool val, long token)
{
	string paramAsString(param);
	return gmmap[token]->SetBoolParam(paramAsString, val);
}

GenModelDLL_API long STDCALL _SetStrParam(char* param, char* val, long token)
{
	string paramAsString(param);
	string valueAsString(val);
	return gmmap[token]->SetStrParam(paramAsString, valueAsString);
}

GenModelDLL_API long STDCALL _CreateNewModel(char model)
{
	long ret = gmmap.size();
	if (model == 'C') // Cplex
	{
		CPLEX_CREATE
	}
	else if (model == 'G') // Gurobi
	{
		GUROBI_CREATE
	}
	else if (model == 'H') // Hyper-graph
	{
		HG_CREATE
	}
	else if (model == 'K') // Glpk
	{
		GLPK_CREATE
	}
	else
	{
		throw new exception();
	}
	return ret;
}

GenModelDLL_API long STDCALL _CopyOrder(long token, int count, int* ind, int* weight)
{
#ifdef CPLEX_SOLVER
	CplexData* d = static_cast<CplexData*>(gmmap[token]->solverdata);
	CPXcopyorder(d->env, d->lp, count, ind, weight, NULL);
#else
	throw new exception();
#endif
	return 0;
}

GenModelDLL_API long STDCALL _DeleteModel(long token)
{
	delete gmmap[token];
	return 0;
}

GenModelDLL_API long STDCALL _CreateModel(long token)
{
	gmmap[token]->SetNumbers();
	gmmap[token]->Init("GenModel");
	gmmap[token]->CreateModel();
	return 0;
}

GenModelDLL_API long STDCALL _SolveModel(long token)
{
	gmmap[token]->Solve();
	gmmap[token]->SetSol();
	return long(gmmap[token]->solstat);
}

GenModelDLL_API bool STDCALL _GetSolVars(double* values, long length, long token)
{
	if (length != gmmap[token]->nc)
	{
		return false;
	}

	for(int i = 0; i < gmmap[token]->nc; i++)
	{
		values[i] = gmmap[token]->vars.sol[i];	// Solution des variables
	}

	return true;
}

GenModelDLL_API bool __stdcall _HasSolution(long token)
{
	//printf("before return %d\n", gmmap[token]->hassolution);
	return gmmap[token]->hassolution;
}

GenModelDLL_API bool STDCALL _GetDualPrices(double* values, long length, long token)
{
	if (length != gmmap[token]->nr)
	{
		return false;
	}

	for(int i = 0; i < gmmap[token]->nr; i++)
	{
		values[i] = gmmap[token]->consts[i].dual;	// Dual prices des contraintes
	}

	return true;
}

GenModelDLL_API bool STDCALL _GetReducedCosts(double* values, long length, long token)
{
	if (length != gmmap[token]->nc)
	{
		return false;
	}

	for(int i = 0; i < gmmap[token]->nc; i++)
	{
		values[i] = gmmap[token]->vars.rc[i];	// Reduced cost for variables
	}

	return true;
}

GenModelDLL_API bool STDCALL _GetRowValues(double* values, long length, long rowIndex, long token)
{
	if (length != gmmap[token]->nc || rowIndex >= gmmap[token]->nr)
	{
		return false;
	}

	memset(values, 0, sizeof(double)*length);
	for(int i = 0; i < int(gmmap[token]->consts[rowIndex].cols.size()); i++)
		values[gmmap[token]->consts[rowIndex].cols[i]] = gmmap[token]->consts[rowIndex].coefs[i];	
	
	return true;
}

GenModelDLL_API bool STDCALL _GetObjCoef(double* values, long length, long token)
{
	if (length != gmmap[token]->nc)
	{
		return false;
	}

	for(int i = 0; i < gmmap[token]->nc; i++)
	{
		values[i] = gmmap[token]->vars.obj[i];	
	}

	return true;
}

GenModelDLL_API bool STDCALL _GetBounds(double* lb, double* ub, long length, long token)
{
	if (length != gmmap[token]->nc)
	{
		return false;
	}

	for(int i = 0; i < gmmap[token]->nc; i++)
	{
		lb[i] = gmmap[token]->vars.lb[i];
		ub[i] = gmmap[token]->vars.ub[i];
	}

	return true;
}

GenModelDLL_API double STDCALL _GetRHS(long row, long token)
{
	return gmmap[token]->consts[row].lrhs;
}

GenModelDLL_API bool STDCALL _SetRHS(long row, double val, long token)
{
	gmmap[token]->consts[row].lrhs = val;
	return true;
}

GenModelDLL_API char STDCALL _GetSense(long row, long token)
{
	return gmmap[token]->consts[row].sense;
}

GenModelDLL_API bool STDCALL _SetSense(long row, char sense, long token)
{
	gmmap[token]->consts[row].sense = sense;
	return true;
}

GenModelDLL_API double STDCALL _GetObjVal(long token)
{
	return gmmap[token]->objval;
}

GenModelDLL_API long STDCALL _ChangeBulkBounds(int count, int* indices, char* types, double* values, long token)
{
	return gmmap[token]->ChangeBulkBounds(count, indices, types, values);
}

GenModelDLL_API long STDCALL _ChangeBulkObjectives(int count, int* indices, double* values, long token)
{
	return gmmap[token]->ChangeBulkObjectives(count, indices, values);
}

GenModelDLL_API long STDCALL _DeleteMipStarts(long token)
{
	return gmmap[token]->DeleteMipStarts();
}

GenModelDLL_API double STDCALL _GetMIPRelativeGap(long token)
{
	return gmmap[token]->GetMIPRelativeGap();
}

GENMODEL_EXTERN_C_END
