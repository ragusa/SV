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
    params.addParam<bool>("isImplicit", true, "implicit or explicit schemes.");
    // Coupled aux variables
    params.addRequiredCoupledVar("h", "height of the fluid");
    params.addRequiredCoupledVar("q_x", "x component of the momentum");
    params.addCoupledVar("q_y", "y component of the momentum");
    params.addCoupledVar("B", "bathymetry data");
    // Gravity
    params.addParam<Real>("gravity", 9.81, "gravity magnitude");
  return params;
}

EntropyFluxAux::EntropyFluxAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Implicit
    _isImplicit(getParam<bool>("isImplicit")),
    // Coupled variable:
    _h(_isImplicit ? coupledValue("h") : coupledValueOld("h")),
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
  // Compute ||q||^2/h^2/h^2
  Real norm_of_q_squared_divdided_h2 = _vector_q.size_sq() / (std::pow(_h[_qp],2));
  // Compute \vec{q}* (g(h+B) + 0.5*||q||^2/h^2/h^2 )
  return _vector_q*( _gravity*(_h[_qp] +_bathymetry[_qp]) + 0.5*norm_of_q_squared_divdided_h2 );
}
