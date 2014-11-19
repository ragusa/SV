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

#include "SV_Momentum.h"
/**
This function computes the x, y and z momentum equationS. It is dimension agnostic. 
 */
template<>
InputParameters validParams<EelMomentum>()
{
  InputParameters params = validParams<Kernel>();
    params.addRequiredCoupledVar("q_x", "x component of momentum");
    params.addCoupledVar("q_y", "y component of momentum");
    params.addCoupledVar("q_z", "z component of momentum");
    params.addCoupledVar("h"  , "h: water height");
//    params.addRequiredCoupledVar("pressure", "pressure");
//    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    params.addParam<int>("component", 0, "component of the momentum equation to compute (0,1,2)->(x,y,z)");
    params.addParam<RealVectorValue>("gravity", (0., 0., 0.), "Gravity vector.");
  return params;
}

SV_Momentum::SV_Momentum(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Coupled auxilary variables
    _q_x(coupledValue("q_x")),
    _q_y(_mesh.dimension()>=2 ? coupledValue("q_y") : _zero),
    _q_z(_mesh.dimension()==3 ? coupledValue("q_z") : _zero),
    //_pressure(coupledValue("pressure")),
    // Equation of state:
    //_eos(getUserObject<EquationOfState>("eos")),
    _h(coupledValue("h")),
    // Parameters:
    _component(getParam<int>("component")),
    _gravity(getParam<RealVectorValue>("gravity")),
    // Parameters for jacobian:
    _rhoA_nb(coupled("rhoA")),
    _q_x_nb(coupled("q_x")),
    _q_y_nb(isCoupled("q_y") ? coupled("q_y") : -1),
    _q_z_nb(isCoupled("q_z") ? coupled("q_z") : -1),


{
    if ( _component > 2 )
        mooseError("ERROR: the integer variable 'component' can only take three values: 0, 1 and 2 that correspond to x, y and z momentum components, respectively.");
}

Real SV_Momentum::computeQpResidual()
{
  // Convection term: _u = rho*vel*vel*A
    RealVectorValue _vector_vel( _q_x[_qp]/_rhoA[_qp], _q_y[_qp]/_rhoA[_qp], _q_z[_qp]/_rhoA[_qp] );
    RealVectorValue _advection = _u[_qp] * _vector_vel;
    
  // Pressure term:
    Real _press = _pressure[_qp]*_area[_qp];
    
  // Source term: P*dA/dx
    Real _PdA = _pressure[_qp]*_grad_area[_qp](_component);
    
  // Wall friction term:
    Real _wall_friction = 0.5 * _friction * _rhoA[_qp] * _vector_vel.size() * _vector_vel(_component) / _Dh;
    
  // Gravity force:
    Real _gravity_force = _gravity(_component) * _rhoA[_qp];
    
  // Return the kernel value:
    return -( _advection*_grad_test[_i][_qp] + _press*_grad_test[_i][_qp](_component) + (_PdA - _wall_friction - _gravity_force)*_test[_i][_qp] );
}

Real SV_Momentum::computeQpJacobian()
{
    // Compute the momentum vector q:
    RealVectorValue _q_vec(_q_x[_qp], _q_y[_qp], _q_z[_qp]);
    
    // Compute the velocity vector:
    RealVectorValue _vector_vel( _q_x[_qp]/_rhoA[_qp], _q_y[_qp]/_rhoA[_qp], _q_z[_qp]/_rhoA[_qp] );
    _vector_vel(_component) *= 2.;
    
    // Compute the derivative of \partial_(x,y,z) (AP) - P \partial_(x,y,z) A:
    Real _press_term = _eos.dAp_dq(_rhoA[_qp], _u[_qp], _rhoEA[_qp])*(_grad_area[_qp](_component)/_area[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](_component));
    
    // Return the value of the jacobian:
    return -_phi[_j][_qp] * ( _vector_vel * _grad_test[_i][_qp] + _press_term );
}

Real SV_Momentum::computeQpOffDiagJacobian( unsigned int _jvar)
{
    // Compute q_vec:
    RealVectorValue _q_vec(_q_x[_qp], _q_y[_qp], _q_z[_qp]);
    
    // Compute the velocity vector:
    RealVectorValue _vector_vel( _q_x[_qp]/_rhoA[_qp], _q_y[_qp]/_rhoA[_qp], _q_z[_qp]/_rhoA[_qp] );
    
    // density (rho*A):
    if (_jvar == _rhoA_nb) {
        Real _press_term = _eos.dAp_drhoA(_rhoA[_qp], _q_vec.size(), _rhoEA[_qp])*(_grad_area[_qp](_component)/_area[_qp]*_test[_i][_qp]+_grad_test[_i][_qp](_component));
        return _phi[_j][_qp] * (_u[_qp] / _rhoA[_qp] * _vector_vel * _grad_test[_i][_qp] - _press_term );
    }
    
    // x-momentum component:
    else if (_jvar == _q_x_nb ) {
        Real _press_term = _eos.dAp_dq(_rhoA[_qp], _q_x[_qp], _rhoEA[_qp])*(_grad_area[_qp](_component)/_area[_qp]*_test[_i][_qp]+_grad_test[_i][_qp](_component));
        return -_phi[_j][_qp] * ( _grad_test[_i][_qp](0)*_u[_qp]/_rhoA[_qp] + _press_term );
    }
    // y-momentum component:
    else if (_jvar == _q_y_nb && _mesh.dimension()>=2 ) {
        Real _press_term = _eos.dAp_dq(_rhoA[_qp], _q_y[_qp], _rhoEA[_qp])*(_grad_area[_qp](_component)/_area[_qp]*_test[_i][_qp]+_grad_test[_i][_qp](_component));
        return -_phi[_j][_qp] * ( _grad_test[_i][_qp](1)*_u[_qp]/_rhoA[_qp] + _press_term );
    }
    // z-momentum component:
    else if (_jvar == _q_z_nb && _mesh.dimension()==3 ) {
        Real _press_term = _eos.dAp_dq(_rhoA[_qp], _q_z[_qp], _rhoEA[_qp])*(_grad_area[_qp](_component)/_area[_qp]*_test[_i][_qp]+_grad_test[_i][_qp](_component));
        return -_phi[_j][_qp] * ( _grad_test[_i][_qp](2)*_u[_qp]/_rhoA[_qp] + _press_term );
    }
    
    // energy (rho*E*A):
    else if (_jvar == _rhoEA_nb) {
        return -_phi[_j][_qp] * _eos.dAp_drhoEA(_rhoA[_qp], _q_vec.size(), _rhoEA[_qp]) * (_grad_area[_qp](_component)/_area[_qp]*_test[_i][_qp]+_grad_test[_i][_qp](_component));
    }
    else
        return 0.;
}
