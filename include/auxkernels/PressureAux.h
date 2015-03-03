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

#ifndef PRESSUREAUX_H
#define PRESSUREAUX_H

#include "AuxKernel.h"
#include "HydroStaticPressure.h"

class PressureAux;

template<>
InputParameters validParams<PressureAux>();

class PressureAux : public AuxKernel
{
public:

  PressureAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  // Coupled variables
  VariableValue & _h;
  VariableValue & _q_x;
  VariableValue & _q_y;

  // Equation of state
  const HydroStaticPressure & _eos;
};

#endif //PRESSUREAUX_H
