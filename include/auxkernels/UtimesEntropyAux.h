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

#ifndef UTIMESENTROPYAUX_H
#define UTIMESENTROPYAUX_H

#include "AuxKernel.h"

class UtimesEntropyAux;

template<>
InputParameters validParams<UtimesEntropyAux>();

class UtimesEntropyAux : public AuxKernel
{
public:

  UtimesEntropyAux(const std::string & name, InputParameters parameters);

protected:
    virtual Real computeValue();
    
    // Implicit:
    const bool & _isImplicit;
    
    // Coupled variable:
    VariableValue & _u;
    VariableValue & _s;
};

#endif //UTIMESENTROPYAUX_H
