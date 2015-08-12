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
This function computes the entropy used in the definition of the entropy viscosity coefficient.
**/
#include "EntropyFluxAux.h"

template<>
InputParameters validParams<EntropyFluxAux>()
{
  InputParameters params = validParams<AuxKernel>();
  // params.addParam<bool>("isImplicit", true, "implicit or explicit schemes.");
  // Coupled aux variables
  params.addRequiredCoupledVar("momentum_component", "current momentum component");
  params.addRequiredCoupledVar("h", "height of the fluid");
  params.addRequiredCoupledVar("q_x", "x-component of the momentum");
  params.addCoupledVar("q_y", "y-component of the momentum");
  params.addCoupledVar("B", "bathymetry data");
  // Gravity
  params.addRequiredParam<Real>("gravity", "gravity magnitude");
  
  return params;
}

EntropyFluxAux::EntropyFluxAux(const InputParameters & parameters) :
                               AuxKernel(parameters),
  // Coupled variables:
  _momentum_cmp(coupledValue("momentum_component")),
  _h(coupledValue("h")),
  _q_x(coupledValue("q_x")),
  _q_y(_mesh.dimension() == 2 ? coupledValue("q_y") : _zero),
  _bathymetry(isCoupled("B") ? coupledValue("B") : _zero),
  _gravity(getParam<Real>("gravity"))
{
}

Real EntropyFluxAux::computeValue()
{
  // Compute the momentum vector q:
  RealVectorValue _vector_q(_q_x[_qp], _q_y[_qp], 0.);
  // Compute the velocity vector
  RealVectorValue _vector_vel = _vector_q / _h[_qp];
  // component of the entropy flux:
  // momentum_component * ( h+B + 0.5*(velocity)^2 )
  return _momentum_cmp[_qp]*( _gravity*(_h[_qp]+_bathymetry[_qp]) + 0.5*_vector_vel.size_sq() );

}
