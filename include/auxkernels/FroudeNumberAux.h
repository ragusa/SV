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

#ifndef FROUDENUMBERAUX_H
#define FROUDENUMBERAUX_H

#include "AuxKernel.h"
#include "HydrostaticPressure.h" 

class FroudeNumberAux;

template<>
InputParameters validParams<FroudeNumberAux>();

class FroudeNumberAux : public AuxKernel
{
public:

  FroudeNumberAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  // Coupled variables
  VariableValue & _h;
  VariableValue & _q_x;
  VariableValue & _q_y;

  // Equation of state:
  const HydrostaticPressure & _eos;

  // Parameters:
  Real _gravity;

};

#endif // FROUDENUMBERAUX_H
