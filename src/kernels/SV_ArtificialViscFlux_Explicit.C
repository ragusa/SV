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

#include "SV_ArtificialViscFlux_Explicit.h"
/**
This function computes the dissipative terms for all of the equations. It is 
dimension-agnostic and equation-agnostic.
 */
template<>
InputParameters validParams<SV_ArtificialViscFlux_Expl>()
{
  InputParameters params = validParams<Kernel>();
  // Equation name:
  params.addParam<std::string>("equation_name", "INVALID", "Name of the equation.");
  
  return params;
}

SV_ArtificialViscFlux_Expl::SV_ArtificialViscFlux_Expl(const InputParameters parameters) :
  Kernel(parameters),
  // Declare equation types
  _equ_type("CONTINUITY X_MOMENTUM Y_MOMENTUM INVALID", getParam<std::string>("equation_name")),
  // Material property: viscosity coefficient.
  _kappa(getMaterialPropertyOld<Real>("kappa"))
{
}

Real SV_ArtificialViscFlux_Expl::computeQpResidual()
{
  // Determine if cell is on boundary or not:
  if (_mesh.isBoundaryNode(_current_elem->node(_i))==true) 
  {
    return 0.;
  }
  else 
  { 
    switch (_equ_type)  // jcr all the same, single line then?
    {
    case CONTINUITY:
      return _kappa[_qp]*_grad_u[_qp]*_grad_test[_i][_qp];
      break;
    case X_MOMENTUM:
      return _kappa[_qp]*_grad_u[_qp]*_grad_test[_i][_qp];
      break;
    case Y_MOMENTUM:
      return _kappa[_qp]*_grad_u[_qp]*_grad_test[_i][_qp];
      break;
    default:
      mooseError("INVALID equation name.");
    }
  }
}

Real SV_ArtificialViscFlux_Expl::computeQpJacobian()
{
  return 0.;
}

Real SV_ArtificialViscFlux_Expl::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}
