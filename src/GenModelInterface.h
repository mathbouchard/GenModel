// GenModelInterface.h : main header file for the GenModelDLL DLL
//

#pragma once

#include "GenModel.h"

#ifdef WIN64
    #define GenModelDLL_API __declspec(dllexport)
    #define STDCALL __stdcall
    #define GENMODEL_EXTERN_C_BEGIN 
    #define GENMODEL_EXTERN_C_END 
#else
    #define GenModelDLL_API
    #define STDCALL
    #define GENMODEL_EXTERN_C_BEGIN extern "C" {
    #define GENMODEL_EXTERN_C_END }
#endif

#ifdef CPLEX_MODULE
	#define CPLEX_CREATE gmmap.push_back(new GenModelCplex())
    #define CPLEX_EXIST true
#else
	#define CPLEX_CREATE throw string("Cplex solver not available")
    #define CPLEX_EXIST false
#endif

#ifdef GUROBI_MODULE
	#define GUROBI_CREATE gmmap.push_back(new GenModelGurobi())
    #define GUROBI_EXIST true
#else
	#define GUROBI_CREATE throw string("Gurobi solver not available")
    #define GUROBI_EXIST false
#endif

#ifdef HG_MODULE
	#define HG_CREATE gmmap.push_back(new GenModelHG())
    #define HG_EXIST true
#else
	#define HG_CREATE throw string("Hypergraph solver not available")
    #define HG_EXIST false
#endif

#ifdef GLPK_MODULE
	#define GLPK_CREATE gmmap.push_back(new GenModelGlpk())
    #define GLPK_EXIST true
#else
	#define GLPK_CREATE throw string("Glpk solver not available")
    #define GLPK_EXIST false
#endif

#ifdef SCIP_MODULE
    #define SCIP_CREATE gmmap.push_back(new GenModelScip())
    #define SCIP_EXIST true
#else
    #define SCIP_CREATE throw string("Scip solver not available")
    #define SCIP_EXIST false
#endif

#ifdef OSI_MODULE
    #define OSI_CREATE gmmap.push_back(new GenModelOsi())
    #define OSI_EXIST true
#else
    #define OSI_CREATE throw string("Coin solver not available")
    #define OSI_EXIST false
#endif

double FindConstraintMaxLhs(long row, long token);
double FindConstraintMinLhs(long row, long token);
long MakeConstraintFeasible(long row, long token);
long WriteProblemToLpFile(char* filename, long token);
long WriteSolutionToFile(char* filename, long token);
long AddConst(char* cname, double rhs, char sense, long token);
bool AddConstBulk(char* cname, double* rhs, long length, char sense, long token);
long AddVar(char* nn, double o, double l, double u, char t, long token);
bool AddVarBulk(char* nn, double* o, long length, double l, double u, char t, long token);
long AddNz(long row, long col, double val, long token);
long AddNzToLast(long col, double val, long token);
long AddNzBulk(long* rows, long* cols, double* values, long valuesLength, long rowCount, long colCount, long iterations, long token);
long SetQpCoef(long i, long j, double val, long token);
long SetNumbers(long token);
long SetLongParam(char* param, long val, long token);
long SetDblParam(char* param, double val, long token);
long SetBoolParam(char* param, bool val, long token);
long SetStrParam(char* param, char* val, long token);
long CreateNewModel(char type, char* name);
bool IsSolverAvailable(char type);
long CopyOrder(long token, int count, int* ind, int* weight);
long DeleteModel(long token);
long CreateModel(long token);
long SolveModel(long token);
bool GetSolVars(double* values, long length, long token);
bool HasSolution(long token);
bool GetDualPrices(double* values, long length, long token);
bool GetSlacks(double* values, long length, long token);
bool GetReducedCosts(double* values, long length, long token);
bool GetRowValues(double* values, long length, long rowIndex, long token);
bool GetObjCoef(double* values, long length, long token);
bool GetBounds(double* lb, double* ub, long length, long token);
double GetLowerBound(long col, long token);
double GetUpperBound(long col, long token);
bool SetLowerBound(long col, double val, long token);
bool SetUpperBound(long col, double val, long token);
double GetRHS(long row, long token);
bool SetRHS(long row, double val, long token);
char GetSense(long row, long token);
bool SetSense(long row, char sense, long token);
double GetObjVal(long token);
long ChangeBulkBounds(int count, int* indices, char* types, double* values, long token);
long ChangeBulkObjectives(int count, int* indices, double* values, long token);
long DeleteMipStarts(long token);
double GetMIPRelativeGap(long token);

template<class T> class InterfaceVector {
public:
    InterfaceVector();
    InterfaceVector(int size);
    ~InterfaceVector();
    void SetSize(int size);
    void Delete();
    void Set(int index, T val);
    T Get(int index);
    T* Ptr();
    
private:
    int size;
    T * val;
};