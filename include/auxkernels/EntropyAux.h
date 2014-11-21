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

#ifndef ENTROPYAUX_H
#define ENTROPYAUX_H

#include "AuxKernel.h"

class EntropyAux;

template<>
InputParameters validParams<EntropyAux>();

class EntropyAux : public AuxKernel
{
public:

  EntropyAux(const std::string & name, InputParameters parameters);

protected:
    virtual Real computeValue();
    
    // Implicit:
    const bool & _isImplicit;
    
    // Aux variables:
    VariableValue & _h;
    VariableValue & _q_x;
    VariableValue & _q_y;
    VariableValue & _q_z;
    
    // Parameters:
    RealVectorValue _gravity;
};

#endif // EntropyAux_H
