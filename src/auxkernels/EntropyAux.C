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
  params.addRequiredCoupledVar("q_x", "x component of the momentum");
  params.addCoupledVar("q_y", "y component of the momentum");
  // Gravity
  params.addParam<Real>("gravity", 9.81, "gravity magnitude");

  return params;
}

EntropyAux::EntropyAux(const std::string & name, InputParameters parameters) :
                       AuxKernel(name, parameters),
  // Implicit
  //_isImplicit(getParam<bool>("isImplicit")),
  // Coupled variable:
  //_h(_isImplicit ? coupledValue("h") : coupledValueOld("h")),  // jcr: et les autres variables?
  _h(coupledValue("h")),
  _q_x(coupledValue("q_x")),
  _q_y(_mesh.dimension() == 2 ? coupledValue("q_y") : _zero),
  //_q_y_old(_mesh.dimension() == 2 ? coupledValueOld("q_y") : _zero),
  // Gravity:
  _gravity(getParam<Real>("gravity"))
{
}

Real EntropyAux::computeValue()
{
  //Real q_y = (_isImplicit ? _q_y : _q_y_old);
  RealVectorValue _vector_q( _q_x[_qp], _q_y[_qp], 0. );
  RealVectorValue _vector_vel = _vector_q / _h[_qp];

  Real entr = 0.5*_h[_qp]* ( _gravity*_h[_qp] + _vector_vel.size_sq() );
  return entr;
}

