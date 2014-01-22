// GenModelDLL.h : main header file for the GenModelDLL DLL
//

#include "GenModel.h"

#ifdef CPLEX_SOLVER
	#define CPLEX_CREATE gmmap["cplex"] = new GenModelCplex();
#else
	#define CPLEX_CREATE throw new exception();
#endif

#ifdef GUROBI_SOLVER
	#define GUROBI_CREATE gmmap["gurobi"] = new GenModelGurobi())
#else
	#define GUROBI_CREATE throw new exception();
#endif

#ifdef HG_SOLVER
    #define HG_CREATE gmmap["hg"] = new GenModelHG();
#else
	#define HG_CREATE throw new exception();
#endif

#ifdef GLPK_SOLVER
	#define GLPK_CREATE gmmap["glpk"] = new GenModelGlpk();
#else
	#define GLPK_CREATE throw new exception();
#endif


