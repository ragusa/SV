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
/**
This function computes the velocity components from the water height and the momentum components.
**/
#include "VelocityAux.h"

template<>
InputParameters validParams<VelocityAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("h", "h");
  params.addRequiredCoupledVar("q", "q_{x,y or z}");
  return params;
}

VelocityAux::VelocityAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    _h(coupledValue("h")),
    _q(coupledValue("q"))
{}

Real
VelocityAux::computeValue()
{
  return _q[_qp] / _h[_qp];
}
