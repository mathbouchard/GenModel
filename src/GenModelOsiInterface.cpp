#include <cassert>

#include "CoinHelperFunctions.hpp"
//#include "CoinIndexedVector.hpp"
#include "ClpQuadraticObjective.hpp"
#include "ClpLinearObjective.hpp"
#include "ClpObjective.hpp"
//#include "ClpSimplex.hpp"
#include "GenModelOsiInterface.h"

long GenModelOsiInterface::setQuadraticObjective(long nc, CoinBigIndex* Q_beg, int* Q_r, double* Q_v)
{
    double offset;
    double csol[2];
    csol[0] = 2.2;
    csol[1] = 0.7;
    
    printf("Set Quadratic\n");
    //ClpObjective * saveObjective  = modelPtr_->objectiveAsObject();
    //modelPtr_->loadQuadraticObjective(nc, Q_beg, Q_r, Q_v);
    
    //double* gradz = modelPtr_->objective(csol, offset, true);
    //for(int i = 0; i < nc; i++)
    //    printf("gradz = %f (%d)\n", gradz[i], Q_beg[i]);
    
    if(quadraticObjective_ != NULL)
       delete quadraticObjective_;

    quadraticObjective_ = new ClpQuadraticObjective(modelPtr_->objective(),nc,Q_beg,Q_r,Q_v);
    double* buf = modelPtr_->objective();
    for(int i = 0; i < nc; i++)
        printf("obj = %f\n", buf[i]);
    
    double* grad = ((ClpQuadraticObjective*)quadraticObjective_)->gradient(modelPtr_, csol, offset, true);
    for(int i = 0; i < nc; i++)
        printf("grad = %f\n", grad[i]);
    
    ClpLinearObjective * linearObjective = new ClpLinearObjective(NULL,nc);
    modelPtr_->setObjectivePointer(linearObjective);
    
    return 0;
}

void GenModelOsiInterface::initialSolve()
{
    // save cutoff
    double cutoff = modelPtr_->dualObjectiveLimit();
    modelPtr_->setDualObjectiveLimit(1.0e50);
    modelPtr_->scaling(0);
    //modelPtr_->setLogLevel(0);
    // solve with no objective to get feasible solution
    setBasis(basis_,modelPtr_);
    printf("before dual\n");
    modelPtr_->dual();
    printf("after dual\n");
    basis_ = getBasis(modelPtr_);
    modelPtr_->setDualObjectiveLimit(cutoff);
    if (modelPtr_->problemStatus())
    {
        assert (modelPtr_->problemStatus()==1);
        return;
    }
    ClpObjective * saveObjective  = modelPtr_->objectiveAsObject();
    modelPtr_->setObjectivePointer(quadraticObjective_);
    //modelPtr_->setLogLevel(1);
    // Could load up any data into a solver
    printf("before primal\n");
    modelPtr_->primal();
    printf("after primal\n");
    modelPtr_->setDualObjectiveLimit(cutoff);
    if (modelPtr_->objectiveValue()>cutoff)
        modelPtr_->setProblemStatus(1);
    // zero reduced costs
    // Should not have to as convex
    //CoinZeroN(modelPtr_->dualRowSolution(),modelPtr_->numberRows());
    //CoinZeroN(modelPtr_->dualColumnSolution(),modelPtr_->numberColumns());
    modelPtr_->setObjectivePointer(saveObjective);
}

void GenModelOsiInterface::resolve()
{
    initialSolve();
}


GenModelOsiInterface::GenModelOsiInterface () : OsiClpSolverInterface()
{
    quadraticObjective_ = NULL;
}


OsiSolverInterface* GenModelOsiInterface::clone(bool CopyData) const
{
    if (CopyData)
        return new GenModelOsiInterface(*this);
    else
        throw string("warning GenModelOsiInterface clone with copyData false\n");
}


GenModelOsiInterface::GenModelOsiInterface (const GenModelOsiInterface & rhs) : OsiClpSolverInterface(rhs)
{
    if (rhs.quadraticObjective_)
        quadraticObjective_=rhs.quadraticObjective_->clone();
    else
        quadraticObjective_=NULL;
}

GenModelOsiInterface::~GenModelOsiInterface()
{
    delete quadraticObjective_;
}

GenModelOsiInterface & GenModelOsiInterface::operator=(const GenModelOsiInterface& rhs)
{
    if (this != &rhs)
    {
        OsiClpSolverInterface::operator=(rhs);
        if (rhs.quadraticObjective_)
            quadraticObjective_=rhs.quadraticObjective_->clone();
        else
            quadraticObjective_=NULL;
    }
    return *this;
}

// Real initializer
void GenModelOsiInterface::initialize ()
{
    printf("!!!! @@@@@@@@@ Initialize @@@@@@@@@ !!!!!!!\n");
    // Save true objective and create a fake one
    delete quadraticObjective_;
    quadraticObjective_ = modelPtr_->objectiveAsObject();
    ClpLinearObjective * linearObjective = new ClpLinearObjective(NULL,modelPtr_->numberColumns());
    modelPtr_->setObjectivePointer(linearObjective);
}

// Get objective function value (can't use default)
double GenModelOsiInterface::getObjValue() const
{
    // first try easy way
    return modelPtr_->objectiveValue();
}