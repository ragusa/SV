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
    params.addParam<Real>("gravity", 9.81, "value of the gravity");
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
    _gravity(getParam<Real>("gravity")),
{
}

Real EntropyFluxAux::computeValue()
{
    Real eflx = _gravity(0)*std::pow(_h[_qp],3) + 
                ( std::pow(_q_x[_qp],3) + std::pow(_q_y[_qp],3) ) / ( std::pow(_h[_qp],2) * 2.);
    return eflx;
}
