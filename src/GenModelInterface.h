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
	#define CPLEX_CREATE gmmap.push_back(new GenModelCplex());
#else
	#define CPLEX_CREATE throw new exception();
#endif

#ifdef GUROBI_MODULEZ
	#define GUROBI_CREATE gmmap.push_back(new GenModelGurobi());
#else
	#define GUROBI_CREATE throw new exception();
#endif

#ifdef HG_MODULEZ
	#define HG_CREATE gmmap.push_back(new GenModelHG());
#else
	#define HG_CREATE throw new exception();
#endif

#ifdef GLPK_MODULEZ
	#define GLPK_CREATE gmmap.push_back(new GenModelGlpk());
#else
	#define GLPK_CREATE throw new exception();
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