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
This function computes the entropy used in the definition of 
the entropy viscosity coefficient.
**/
#include "EntropyAux.h"

template<>
InputParameters validParams<EntropyAux>()
{
  InputParameters params = validParams<AuxKernel>();
  // params.addParam<bool>("isImplicit", true, "implicit or explicit schemes.");
  // Coupled variables
  params.addRequiredCoupledVar("h", "height of the fluid");
  params.addRequiredCoupledVar("q_x", "x-component of the momentum");
  params.addCoupledVar("q_y", "y-component of the momentum");
  // Gravity
  params.addRequiredParam<Real>("gravity", "gravity magnitude");

  return params;
}

EntropyAux::EntropyAux(const InputParameters & parameters) :
                       AuxKernel(parameters),
  _h(coupledValue("h")),
  _q_x(coupledValue("q_x")),
  _q_y(_mesh.dimension() == 2 ? coupledValue("q_y") : _zero),
  // Gravity:
  _gravity(getParam<Real>("gravity"))
{
}

Real EntropyAux::computeValue()
{
  RealVectorValue _vector_q( _q_x[_qp], _q_y[_qp], 0. );
  RealVectorValue _vector_vel = _vector_q / _h[_qp];

  // size_sq returns the square of the vector, hence u^2+v^2
  Real entr = 0.5*_h[_qp]* ( _gravity*_h[_qp] + _vector_vel.size_sq() );
  return entr;
}

