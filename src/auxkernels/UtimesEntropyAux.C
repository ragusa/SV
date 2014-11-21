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
#include "UtimesEntropyAux.h"

template<>
InputParameters validParams<UtimesEntropyAux>()
{
  InputParameters params = validParams<AuxKernel>();
    params.addParam<bool>("isImplicit", true, "implicit or explicit schemes.");
    params.addRequiredCoupledVar("u", "variable it is solved for.");
    params.addRequiredCoupledVar("s", "Entropy function.");
  return params;
}

UtimesEntropyAux::UtimesEntropyAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Implcit
    _isImplicit(getParam<bool>("isImplicit")),
    // Coupled variable:
    _u(_isImplicit ? coupledValue("u") : coupledValueOld("u")),
    _s(_isImplicit ? coupledValue("s") : coupledValueOld("s"))
{}

Real
UtimesEntropyAux::computeValue()
{
    return _u[_qp]*_u[_qp]*_u[_qp]/3.;
}
