// GenModelDLL.h : main header file for the GenModelDLL DLL
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

#ifdef CPLEX_SOLVER
	#define CPLEX_CREATE gmmap.push_back(new GenModelCplex());
#else
	#define CPLEX_CREATE throw new exception();
#endif

#ifdef GUROBI_SOLVER
	#define GUROBI_CREATE gmmap.push_back(new GenModelGurobi());
#else
	#define GUROBI_CREATE throw new exception();
#endif

#ifdef HG_SOLVER
	#define HG_CREATE gmmap.push_back(new GenModelHG());
#else
	#define HG_CREATE throw new exception();
#endif

#ifdef GLPK_SOLVER
	#define GLPK_CREATE gmmap.push_back(new GenModelGlpk());
#else
	#define GLPK_CREATE throw new exception();
#endif


