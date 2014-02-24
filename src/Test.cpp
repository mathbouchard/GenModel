/***************************************************************************
 *  Test.cpp
 *  An application to test GenModel
 *
 *  January 8 11:32 2014
 *  Copyright  2014  Mathieu Bouchard
 *  mathbouchard@gmail.com
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "GenModelInterface.h"

using namespace std;

int main(int argc, char** argv)
{
    char solver = 'C';
    bool quadratic = false;
    bool mip = true;
    
    printf("Solving problem from Integer Programming book p.130\n\n");
    printf("z = max 4 x1 - 1 x2\n");
    printf("        7 x1 - 2 x2 <= 14\n");
    printf("               1 x2 <= 3\n");
    printf("        2 x1 - 2 x2 <= 3\n\n");
    printf("The problem has the optimal solution x=(2,1/2) with objective value 7.5\n\n");
    
    if(argc > 1)
    {
        printf("Using solver option %c\n", argv[1][0]);
        solver = argv[1][0];
    }
    if(argc > 2)
    {
        quadratic = argv[2][0] != '0';
        if(quadratic)
            printf("Solving quadratic problem\n");
    }
    if(argc > 3)
    {
        mip = argv[3][0] != '0';
        if(mip)
            printf("Solving mip problem\n");
    }
    
    if(IsSolverAvailable('C')) printf("Cplex available\n");
    if(IsSolverAvailable('G')) printf("Gurobi available\n");
    if(IsSolverAvailable('H')) printf("Hypergraph available\n");
    if(IsSolverAvailable('K')) printf("Glpk available\n");
    if(IsSolverAvailable('S')) printf("Scip available\n");
    if(IsSolverAvailable('O')) printf("Coin-OR available\n");
    
    char version[] = "print_version";
    char epgap[] = "relative_mip_gap_tolerance";
    char timelimit[] = "time_limit";
    char probname[] = "Genmodel Test";
    char problem_type[] = "mip";
    char problem_objective[] = "maximize";
    char log_param[] = "log_file";
    char log_file[] = "tmp/test.log";
    char lp_file[] = "tmp/test.lp";
    char sol_file[] = "tmp/test.sol";
    char qp_string[] = "qp_mat";
    char x1[] = "x1";
    char x2[] = "x2";
    char c1[] = "c1";
    char c2[] = "c2";
    char c3[] = "c3";
    
    printf("Initializing problem... ");
    long model_id = CreateNewModel(solver, probname);
    printf("done.\n");
    printf("Setting problem... ");
    if(mip)
        SetBoolParam(problem_type, true, model_id);
    if(quadratic)
        SetBoolParam(qp_string, true, model_id);
    SetBoolParam(version, true, model_id);
    SetBoolParam(problem_objective, true, model_id);
    SetDblParam(timelimit, 120, model_id);
    SetDblParam(epgap, 0.01, model_id);
    SetStrParam(log_param,log_file, model_id);
    AddVar(x1, 4.0, 0.0, numeric_limits<double>::infinity(), (mip ? 'I' : 'C'), model_id);
    if(quadratic)
        SetQpCoef(0,0,-2.5,model_id);
    AddVar(x2, -1.0, 0.0, numeric_limits<double>::infinity(), 'C', model_id);
    
    AddConst(c1, 14.0, 'L', model_id);
    AddNzToLast(0, 7.0, model_id);
    AddNzToLast(1, -2.0, model_id);
    
    AddConst(c2, 3.0, 'L', model_id);
    AddNzToLast(1, 1.0, model_id);
    
    AddConst(c3, 3.0, 'L', model_id);
    AddNzToLast(0, 2.0, model_id);
    AddNzToLast(1, -2.0, model_id);
    
    //printf("Min => #{Genmodel.FindConstraintMinLhs(0, model_id)}\n");
    //printf("Max => #{Genmodel.FindConstraintMaxLhs(0, model_id)}\n");
    //MakeConstraintFeasible(0, model_id);
    //MakeConstraintFeasible(2, model_id);
    printf("done.\n");
    printf("Building problem... ");
    CreateModel(model_id);
    printf("done.\n");
    printf("Solving problem... ");
    long solution_status=SolveModel(model_id);
    printf("done.\n");
    bool has_solution = HasSolution(model_id);
    printf("Retrieving objective value... ");
    double objective_value = GetObjVal(model_id);
    printf("done.\n");
    printf("solution_status=%ld\nhas_solution=%d\nobjective_value=%f\n", solution_status, has_solution, objective_value);
    
    double solutions[2];
    double rcs[2];
    double duals[3];
    double slacks[3];
    
    try {
        printf("Writing to lp file... ");
        WriteProblemToLpFile(lp_file, model_id);
        printf("done.\n");
        printf("Writing to solution file... ");
        WriteSolutionToFile(sol_file, model_id);
        printf("done.\n");
    } catch (string e) {
        printf("%s\n",e.c_str());
    }
    
    
    printf("Retrieving solution variables... ");
    GetSolVars(solutions,2,model_id);
    GetReducedCosts(rcs, 2, model_id);
    printf("done.\n");
    for(int i = 0; i < 2; i++)
        printf("x_%d -> %f [rc = %f]\n", i, solutions[i], rcs[i]);
    
    printf("Retrieving solution dual prices... ");
    GetDualPrices(duals,3,model_id);
    GetSlacks(slacks,3,model_id);
    printf("done.\n");
    for(int i = 0; i < 3; i++)
        printf("c_%d -> %f [slack = %f]\n", i, duals[i], slacks[i]);
    
    return 0;
}
