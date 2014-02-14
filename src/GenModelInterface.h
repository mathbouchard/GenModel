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
    #define SCIP_CREATE gmmap.push_back(new GenModelGlpk())
    #define SCIP_EXIST true
#else
    #define SCIP_CREATE throw string("Scip solver not available")
    #define SCIP_EXIST false
#endif

#ifdef COIN_MODULE
    #define COIN_CREATE gmmap.push_back(new GenModelCoin())
    #define COIN_EXIST true
#else
    #define COIN_CREATE throw string("Scip solver not available")
    #define COIN_EXIST false
#endif

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