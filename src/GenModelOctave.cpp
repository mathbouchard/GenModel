// GenModelDLL.cpp : Defines the initialization routines for the DLL.
//

#include <octave/oct.h>
#include "GenModelOctave.h"

//#ifdef CPLEX_SOLVER
	#include "GenModelCplex.h"
//#endif

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

static bool checkargs(const octave_value_list& args);

map<string,GenModel*> gmmap;

DEFUN_DLD(gensolve, args, , "Return the solution of Ax=b")
{
    printf("\toct : Checking args... ");
    char tmp[4096];
    if (checkargs(args))
        return octave_value_list();
    printf("done.\n");

    printf("\toct : Initiating solver... ");
    gmmap["cplex"] = new GenModelCplex();
    GenModel* gmp = gmmap["cplex"];
    gmp->Init("OctaveModel");
    printf("done.\n");

    printf("\toct : Preparing data... \n");
    int nr = args(0).rows();
    int nc = args(0).columns();
    SparseMatrix sm = args(0).sparse_matrix_value();
    int nz = sm.nnz();

    printf("\t\toct : row = %d, col = %d, nz = %d\n", nr, nc, nz);

    RowVector c_obj = (args(1).rows() == args(0).columns() ? args(1).matrix_value().column(0).transpose() : args(1).matrix_value().row(0));
    RowVector c_lb = (args(2).rows() == args(0).columns() ? args(2).matrix_value().column(0).transpose() : args(2).matrix_value().row(0));
    RowVector c_ub = (args(3).rows() == args(0).columns() ? args(3).matrix_value().column(0).transpose() : args(3).matrix_value().row(0));
    string c_type = (args(4).rows() == args(0).columns() ? args(4).char_matrix_value().transpose().row_as_string(0) : args(4).char_matrix_value().row_as_string(0));
    string r_type = (args(5).columns() == args(0).rows() ? args(5).char_matrix_value().row_as_string(0) : args(5).char_matrix_value().transpose().row_as_string(0));
    ColumnVector r_rhs = (args(6).columns() == args(0).rows() ? args(6).matrix_value().row(0).transpose() : args(6).matrix_value().column(0));
    bool max = args(7).int_value();
    bool mip = false;

    gmp->SetBoolParam("maximize", max);

    printf("\toct : done.\n");
    printf("\toct : Creating problem from octave data... ");
    for(int i = 0; i < nc; i++)
    {
	snprintf(tmp, 4096, "Vars_%d", i);
	//printf("\t\toct : Adding col %s %f %f %f %c\n", tmp, c_obj(i), c_lb(i), c_ub(i), c_type[i]); 
	gmp->AddVar(tmp, c_obj(i), c_lb(i), c_ub(i), c_type[i]);
	if(c_type[i] == 'I' || c_type[i] == 'B')
	    mip = true;
    }
    for(int i = 0; i < nr; i++)
    {
        snprintf(tmp, 4096, "Consts_%d", i);
        //printf("\t\toct : Adding consts %s %f %c\n", tmp, r_rhs(i), r_type[i]); 
        gmp->AddConst(tmp, r_rhs(i), r_type[i]);
    }

    int tot=0;
    for(int j = 0; j < nc; j++)
    {
	for (int i = sm.cidx(j); i < sm.cidx(j+1); i++)
        {
            //printf("\t\toct : %d %d(%d) => %f\n", sm.ridx(i), j, sm.cidx(i), sm.data(i));
            gmp->AddNz(sm.ridx(i), j, sm.data(i));
	}
    }
    printf("done.\n");

    printf("\toct : Setting numbers... ");
    gmp->SetNumbers();
    printf("done.\n");
    printf("\toct : Creating model... ");
    gmp->CreateModel();
    printf("done.\n");
    printf("\toct : Solving... \n");
    gmp->Solve();
    printf("done.\n");
    printf("\toct : Setting solution... ");
    gmp->SetSol();
    printf("done.\n");

    intNDArray<octave_int32> c_stat(dim_vector(1,nc));
    intNDArray<octave_int32> r_stat(dim_vector(nr,1));
    int* c_stat_vec = static_cast<int*>(c_stat.mex_get_data());
    int* r_stat_vec = static_cast<int*>(r_stat.mex_get_data());

    CPXgetbase(((CplexData*)(gmp->solverdata))->env, ((CplexData*)(gmp->solverdata))->lp, c_stat_vec, r_stat_vec);

    double* binvci = new double[nr];
    int invnz = 0;
    for(int i = 0; i < nr; i++)
    {
        CPXbinvcol(((CplexData*)(gmp->solverdata))->env, ((CplexData*)(gmp->solverdata))->lp, i, binvci);
	for(int j = 0; j < nr; j++)
	    if(binvci[j] != 0.0)
		++invnz;
    }
    SparseMatrix binv(nr,nr,invnz);
    for(int i = 0; i < nr; i++)
    {
        CPXbinvcol(((CplexData*)(gmp->solverdata))->env, ((CplexData*)(gmp->solverdata))->lp, i, binvci);
        for(int j = 0; j < nr; j++)
        {
            if(binvci[j] != 0.0)
            {
		binv(j,i) = binvci[j];
            }
        }
    }
    delete[] binvci;

    RowVector sol(c_obj.columns());
    RowVector rc(c_obj.columns());
    for(int i = 0; i < nc; i++)
    {
	//printf("sol %d = %f\n", i, gmp->vars.sol[i]);
	sol(i) = gmp->vars.sol[i];
    }
    for(int i = 0; i < nc && !mip; i++)
        rc(i) = gmp->vars.rc[i];
    ColumnVector dual(r_rhs);
    ColumnVector slack(r_rhs);
    for(int i = 0; i < nr && !mip; i++)
    {
        dual(i) = gmp->consts[i].dual;
	slack(i) = gmp->consts[i].slack;
    }

    octave_value_list ret = octave_value_list();
    ret.append(octave_value(sol));
    ret.append(octave_value(rc));
    ret.append(octave_value(dual));
    ret.append(octave_value(slack));
    ret.append(octave_value(c_stat));
    ret.append(octave_value(r_stat));
    ret.append(octave_value(binv));

    return ret;
}

static bool checkargs(const octave_value_list& args)
{
    if(args.length() < 8)
    {
        error("gensolve: expecting eight arguments : \n"
		"\tA : (n,m) sparse matrix\n"
		"\tc_obj : (1,m) or (m,1) vector\n"  
		"c_lb : (1,m) or (m,1) vector\n"
		"c_ub : (1,m) or (m,1) vector\n"
		"c_type : (1,m) or (m,1) vector\n"
		"r_type (r_lb) : (n,1) or (1,n) vector\n"
		"r_rhs (r_ub) : (n,1) or (1,n) vector\n"
		"minormax : 0 = minimize, 1 = maximize\n");
        return true;
    }
    if (!args(0).is_sparse_type())
    {
        error("gensolve: expecting matrix A (arg 1) to be a rsparsematrix");
        return true;
    }

    if (args(1).rows() != args(0).columns() && args(1).columns() != args(0).columns())
    {
        error("gensolve: expecting matrix A (arg 1) to have the same number of columns as c_obj (arg 2)");
        return true;
    }

    if (!args(1).is_real_matrix())
    {
        error("gensolve: excpecting c_obj (arg 2) must be a real svector");
        return true;
    }

    if (args(2).rows() != args(0).columns() && args(2).columns() != args(0).columns())
    {
        error("gensolve: expecting matrix A (arg 1) to have the same number of columns as c_lb (arg 3)");
        return true;
    }

    if (!args(2).is_real_matrix())
    {
        error("gensolve: excpecting c_lb (arg 3) must be a real vector");
        return true;
    }

    if (args(3).rows() != args(0).columns() && args(3).columns() != args(0).columns())
    {
        error("gensolve: expecting matrix A (arg 1) to have the same number of columns as c_ub (arg 4)");
        return true;
    }

    if (!args(3).is_real_matrix())
    {
        error("gensolve: excpecting c_ub (arg 4) must be a real vector");
        return true;
    }

    if (args(4).rows() != args(0).columns() && args(4).columns() != args(0).columns())
    {
        error("gensolve: expecting matrix A (arg 1) to have the same number of columns as c_type (arg 5)");
        return true;
    }

    if (args(5).rows() != args(0).rows() && args(5).columns() != args(0).rows())
    {
        error("gensolve: expecting matrix A (arg 1) to have the same number of crow as r_type (or r_lb) (arg 6)");
        return true;
    }

    if (args(6).rows() != args(0).rows() && args(6).columns() != args(0).rows())
    {
        error("gensolve: expecting matrix A (arg 1) to have the same number of rows as r_rhs (or r_ub) (arg 7)");
        return true;
    }

    if (args(7).scalar_value() < 0.0 || args(7).scalar_value() > 1.0) 
    {
        error("matpow: objective diretion (arg 8) must be 0 (min) or 1 (max)");
        return true;
    }

    return false;
}

/*
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
*/
