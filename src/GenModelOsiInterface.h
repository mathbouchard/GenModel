/***************************************************************************
 *  GenModelOsiInterface.h
 *  Replacememnt for the default OsiClpSolverInterface
 *
 *  February 14 9:32 2014
 *  Copyright  2014  Mathieu Bouchard
 *  mathbouchard@gmail.com
 ****************************************************************************/

#ifndef GENMODELOSIINTERFACE_H_
#define GENMODELOSIINTERFACE_H_

#ifndef CBC_THREAD
    #define CBC_THREAD
#endif

#include "OsiClpSolverInterface.hpp"
#include <string>

using namespace std;

// This is to allow the user to replace initialSolve and resolve

class GenModelOsiInterface : public OsiClpSolverInterface {

public:
    
    GenModelOsiInterface ();
    GenModelOsiInterface (const GenModelOsiInterface &);
    virtual ~GenModelOsiInterface ();
    virtual OsiSolverInterface * clone(bool CopyData=true) const;
    GenModelOsiInterface & operator=(const GenModelOsiInterface& rhs);
    
    // Solve initial LP relaxation
    virtual void initialSolve();
    /// Resolve an LP relaxation after problem modification
    virtual void resolve();
    // Setup fake objective or somehow get nonlinear  info
    void initialize();
    virtual double getObjValue() const;
    long setQuadraticObjective(long nc, CoinBigIndex* Q_beg, int* Q_r, double* Q_v);
    
private:
    ClpObjective * quadraticObjective_;
};

#endif //GENMODELOSIINTERFACE_H_
