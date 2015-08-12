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
This function computes the velocity components from the water height 
and the momentum components.
**/
#include "VelocityAux.h"

template<>
InputParameters validParams<VelocityAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("h", "h");
  params.addRequiredCoupledVar("q_x", "q_{x}");
  params.addCoupledVar("q_y", "q_{y}");
  return params;
}

VelocityAux::VelocityAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    _h(coupledValue("h")),
    _q_x(coupledValue("q_x")),
    _q_y(_mesh.dimension() == 2 ? coupledValue("q_y") : _zero)
{}

Real VelocityAux::computeValue()
{
  if(_mesh.dimension() == 2)
  {
	return std::sqrt(_q_x[_qp]*_q_x[_qp]+_q_y[_qp]*_q_y[_qp]) / _h[_qp];
  }
  else
  {
    return _q_x[_qp] / _h[_qp];
  }
}
