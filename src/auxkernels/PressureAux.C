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
#include "PressureAux.h"

template<>
InputParameters validParams<PressureAux>()
{
  InputParameters params = validParams<AuxKernel>();
  
  // Coupled variables
  params.addRequiredCoupledVar("h", "water height");
  params.addRequiredCoupledVar("q_x", "x component of momentum");
  params.addCoupledVar("q_y", "y component of momentum");
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");

  return params;
}

PressureAux::PressureAux(const std::string & name, InputParameters parameters) :
                         AuxKernel(name, parameters),
    // Coupled variables
    _h(coupledValue("h")),
    _q_x(coupledValue("q_x")),
    _q_y(_mesh.dimension() == 2 ? coupledValue("q_y") : _zero),
    // Equation of state:
    _eos(getUserObject<HydrostaticPressure>("eos"))
{}

Real
PressureAux::computeValue()
{
  RealVectorValue _vector_q(_q_x[_qp], _q_y[_qp], 0.);
  return _eos.pressure(_h[_qp], _vector_q);
}
