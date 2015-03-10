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

#include "SV_ArtificialViscFlux.h"
/**
This function computes the dissipative terms for all of the equations. It is dimension agnostic.
 */
template<>
InputParameters validParams<SV_ArtificialViscFlux>()
{
  InputParameters params = validParams<Kernel>();
    // Equation name:
    params.addParam<std::string>("equation_name", "INVALID", "Name of the equation.");
    // Coupled aux variables
//    params.addRequiredCoupledVar("h", "height of the fluid");
//    params.addRequiredCoupledVar("velocity_x", "x component of the velocity");
//    params.addCoupledVar("velocity_y", "y component of the velocity");
  return params;
}

SV_ArtificialViscFlux::SV_ArtificialViscFlux(const std::string & name,
                                     InputParameters parameters) :
  Kernel(name, parameters),
    // Declare equation types
    _equ_type("CONTINUITY, X_MOMENTUM, Y_MOMENTUM, INVALID", getParam<std::string>("equation_name")),
    // Coupled auxiliary variables
//    _grad_h(coupledGradient("h")),
//    _vel_x(coupledValue("velocity_x")),
//    _vel_y(_mesh.dimension()>=2 ? coupledValue("velocity_y") : _zero),
//    _grad_vel_x(coupledGradient("velocity_x")),
//    _grad_vel_y(_mesh.dimension()>=2 ? coupledGradient("velocity_y") : _grad_zero),
    // Material property: viscosity coefficient.
    _kappa(getMaterialProperty<Real>("kappa"))
{
}

Real SV_ArtificialViscFlux::computeQpResidual()
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

Real SV_ArtificialViscFlux::computeQpJacobian()
{
    // We assumed that the all of the above regularization can be approximated by the parabolic regularization:
    return _kappa[_qp]*_grad_phi[_j][_qp]*_grad_test[_i][_qp];
}

Real SV_ArtificialViscFlux::computeQpOffDiagJacobian( unsigned int _jvar)
{
    // With above assumption, there is no contribution from the dissipative terms to the off diagonal terms of the jacobian matrix.
    return 0.*_jvar;
}
