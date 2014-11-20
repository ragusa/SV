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

#include "SV_TimeDerivative.h"

template<>
InputParameters validParams<SV_TimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  return params;
}

SV_TimeDerivative::SV_TimeDerivative(const std::string & name,
                                             InputParameters parameters) :
    TimeDerivative(name,parameters)
{}

Real SV_TimeDerivative::computeQpResidual()
{
    return TimeDerivative::computeQpResidual();
}

Real SV_TimeDerivative::computeQpJacobian()
{
    return TimeDerivative::computeQpJacobian();
}
