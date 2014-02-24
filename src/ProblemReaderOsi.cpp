/*
 * ProblemReaderOsi.cpp
 *
 *  Created on: 2012-10-02
 *      Author: mbouchard
 */

#include "ProblemReaderOsi.h"

int ReadFromFile(GenModel* pgm, string fn, int type, string dn)
{
    OsiSolverInterface* pm = new OsiClpSolverInterface();
    //OsiSolverInterface* pm = new OsiGlpkSolverInterface();

    switch (type)
    {
        case 0: pm->readMps(fn.c_str()); break;
        case 1: pm->readLp(fn.c_str()); break;
        case 2: pm->readGMPL(fn.c_str(), (dn == "" ? NULL : dn.c_str())); break;
        default: pm->readMps(fn.c_str()); break;
    }
    
    return ReadFromObject(pgm, pm);
}

int ReadFromObject(GenModel* pgm, OsiSolverInterface* pm)
{
    char tmp[64];
    int numr = pm->getNumRows();
    int numc = pm->getNumCols();
    pgm->nz = 0;
    
    const double* obj = pm->getObjCoefficients();
    const double* lb = pm->getColLower();
    const double* ub = pm->getColUpper();
    const char* type = pm->getColType();
    for(int i = 0; i < numc; i++)
    {
        char gmtype;
        snprintf(tmp,64,"C%d", i);
        switch (type[i])
        {
            case 0: gmtype = 'C'; break;
            case 1: gmtype = 'B'; break;
            case 2: gmtype = 'I'; break;
            default: gmtype = 'C'; break;
        }
        pgm->AddVar(string(tmp), obj[i],
                    (lb[i]==COIN_DBL_MAX?numeric_limits<double>::infinity():(lb[i]==-COIN_DBL_MAX?-numeric_limits<double>::infinity():lb[i])),
                    (ub[i]==COIN_DBL_MAX?numeric_limits<double>::infinity():(ub[i]==-COIN_DBL_MAX?-numeric_limits<double>::infinity():ub[i])), gmtype);
    }
    const char* sense = pm->getRowSense();
    const double* lrhs = pm->getRowLower();
    const double* urhs = pm->getRowUpper();
    
    const CoinPackedMatrix* pcpm = pm->getMatrixByRow();
    const double* elem = pcpm->getElements();
    const int* ind = pcpm->getIndices();
    
    for(int i = 0; i < numr; i++)
    {
        snprintf(tmp,64,"R%d", i);
        if(sense[i] == 'L')
            pgm->AddConst(string(tmp),
                (urhs[i]==COIN_DBL_MAX?numeric_limits<double>::infinity():(urhs[i]==-COIN_DBL_MAX?-numeric_limits<double>::infinity():urhs[i])), sense[i]);
        else
            pgm->AddConst(string(tmp),
                (lrhs[i]==COIN_DBL_MAX?numeric_limits<double>::infinity():(lrhs[i]==-COIN_DBL_MAX?-numeric_limits<double>::infinity():lrhs[i])), sense[i]);
        pgm->consts.back().urhs = (urhs[i]==COIN_DBL_MAX?numeric_limits<double>::infinity():(urhs[i]==-COIN_DBL_MAX?-numeric_limits<double>::infinity():urhs[i]));
        //if(fabs(lrhs[i]) > 1e-4 || fabs(urhs[i]) > 1e-4) //What is that ???
        for(int j = pcpm->getVectorFirst(i); j < pcpm->getVectorLast(i); j++)
        {
            pgm->AddNzToLast(ind[j], elem[j]);
        }
    }
    
    return 0;
}