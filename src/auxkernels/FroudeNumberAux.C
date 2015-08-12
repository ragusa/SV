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
This function computes the density of the fluid.
**/
#include "FroudeNumberAux.h"

template<>
InputParameters validParams<FroudeNumberAux>()
{
  InputParameters params = validParams<AuxKernel>();
  
  // Coupled variables
  params.addRequiredCoupledVar("h", "water height");
  params.addRequiredCoupledVar("q_x", "x-component of h*vec{u}");
  params.addCoupledVar("q_y", "y-component of h*vec{u}");
  // Eos 
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  // Gravity
  params.addRequiredParam<Real>("gravity", "Gravity magnitude");
  
  return params;
}

FroudeNumberAux::FroudeNumberAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  // Coupled variables
  _h(coupledValue("h")),
  _q_x(coupledValue("q_x")),
  _q_y(_mesh.dimension() == 2 ? coupledValue("q_y") : _zero),
  // Equation of state:
  _eos(getUserObject<HydrostaticPressure>("eos")),
  // gravity
  _gravity(getParam<Real>("gravity"))
{}

Real
FroudeNumberAux::computeValue()
{
  // Compute the momentum vector
  RealVectorValue _vector_q( _q_x[_qp], _q_y[_qp], 0. );

  // Return the value of the Froude number: Fr = u/c
  return  _vector_q.size()/( _h[_qp] * std::sqrt(_eos.c2(_h[_qp], _vector_q)) );
}
