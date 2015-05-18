/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef ENTROPYFLUXAUX_H
#define ENTROPYFLUXAUX_H

#include "AuxKernel.h"

class EntropyFluxAux;

template<>
InputParameters validParams<EntropyFluxAux>();

class EntropyFluxAux : public AuxKernel
{
public:

  EntropyFluxAux(const std::string & name, InputParameters parameters);

protected:
    virtual Real computeValue();
    
    // Implicit:
    //const bool & _isImplicit;
    
    // Aux variables:
    VariableValue & _curr_momentum;
    VariableValue & _h;
    VariableValue & _q_x;
    VariableValue & _q_y;
    VariableValue & _bathymetry;
    
    // Parameters:
    Real _gravity;
};

#endif // ENTROPYFLUXAUX_H
