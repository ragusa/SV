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

#include "SV_ArtificialVisc.h"
/**
This function computes the dissipative terms for all of the equations. It is dimension agnostic.
 */
template<>
InputParameters validParams<SV_ArtificialVisc>()
{
  InputParameters params = validParams<Kernel>();
    // Equation and diffusion names:
    params.addParam<std::string>("equation_name", "INVALID", "Name of the equation.");
    params.addParam<std::string>("diffusion_name", "PARABOLIC", "Name of the diffusion.");
    // Coupled aux variables
    params.addRequiredCoupledVar("h", "height of the fluid");
    params.addRequiredCoupledVar("velocity_x", "x component of the velocity");
    params.addCoupledVar("velocity_y", "y component of the velocity");
  return params;
}

SV_ArtificialVisc::SV_ArtificialVisc(const std::string & name,
                                     InputParameters parameters) :
  Kernel(name, parameters),
    // Declare equation types
    _equ_name(getParam<std::string>("equation_name")),
    _diff_name(getParam<std::string>("diffusion_name")),
    _equ_type("CONTINUITY, XMOMENTUM, YMOMENTUM, INVALID", _equ_name),
    _diff_type("ENTROPY, PARABOLIC, NONE, INVALID",_diff_name),
    // Coupled auxiliary variables
    //_h(coupledValue("h")),
    _grad_h(coupledGradient("h")),
    _vel_x(coupledValue("velocity_x")),
    _vel_y(_mesh.dimension()>=2 ? coupledValue("velocity_y") : _zero),
    _grad_vel_x(coupledGradient("velocity_x")),
    _grad_vel_y(_mesh.dimension()>=2 ? coupledGradient("velocity_y") : _grad_zero),
    // Material property: viscosity coefficient.
    _mu(getMaterialProperty<Real>("mu")),
    _kappa(getMaterialProperty<Real>("kappa"))
{
//    _equ_type = _equ_name;
//    _diff_type = _diff_name;
}

Real SV_ArtificialVisc::computeQpResidual()
{
    // Determine if cell is on boundary or not and then compute a unit vector 'l=grad(norm(vel))/norm(grad(norm(vel)))':
    Real isonbnd = 1.;
    if (_mesh.isBoundaryNode(_current_elem->node(_i))==true) {
        isonbnd = 0.; //jcr note: not used???
    }

    // If statement on diffusion type:
    if (_diff_type == 1) {
        switch (_equ_type) {
            case CONTINUITY:
                return _kappa[_qp]*_grad_h[_qp]*_grad_test[_i][_qp];
                break;
            case XMOMENTUM:
                return _kappa[_qp]*(_h[_qp]*_grad_vel_x[_qp]+_vel_x[_qp]*_grad_h[_qp])*_grad_test[_i][_qp];
                break;
            case YMOMENTUM:
                return _kappa[_qp]*(_h[_qp]*_grad_vel_y[_qp]+_vel_y[_qp]*_grad_h[_qp])*_grad_test[_i][_qp];
                break;
            default:
                mooseError("INVALID equation name.");
        }
    }
    else if (_diff_type == 0) {
        // Compute 0.5*rho*mu*grad(vel)_symmetric: (get a symmetric tensor)
        TensorValue<Real> grad_vel_tensor(_grad_vel_x[_qp], _grad_vel_y[_qp], _grad_vel_z[_qp]);
        grad_vel_tensor = ( grad_vel_tensor + grad_vel_tensor.transpose() );
        grad_vel_tensor *= 0.5 * _h[_qp] * _mu[_qp];
        
        // Compute f = kappa * grad(rho):
        RealVectorValue f(_kappa[_qp]*_grad_h[_qp](0), _kappa[_qp]*_grad_h[_qp](1), _kappa[_qp]*_grad_h[_qp](2));
        
        // return the dissipative terms:
            switch (_equ_type) {
                case CONTINUITY: // div(kappa grad(rho))
                    return isonbnd* f * _grad_test[_i][_qp];
                    break;
                case XMOMENTUM:
                    return isonbnd*( _vel_x[_qp]*f + grad_vel_tensor.row(0) ) * _grad_test[_i][_qp];
                    break;
                case YMOMENTUM:
                    return isonbnd*( _vel_y[_qp]*f + grad_vel_tensor.row(1) ) * _grad_test[_i][_qp];
                    break;
                default:
                    mooseError("INVALID equation name.");
            }
    }
    else {
        mooseError("INVALID dissipation terms.");
    }
}

Real SV_ArtificialVisc::computeQpJacobian()
{
    // We assumed that the all of the above regularization can be approximated by the parabolic regularization:
    Real _mu_jac = std::max(_mu[_qp], _kappa[_qp]);
    return _mu_jac*_grad_phi[_j][_qp]*_grad_test[_i][_qp];
}

Real SV_ArtificialVisc::computeQpOffDiagJacobian( unsigned int _jvar)
{
    // With above assumption, there is no contribution from the dissipative terms to the off diagonal terms of the jacobian matrix.
    return 0.*_jvar;
}
