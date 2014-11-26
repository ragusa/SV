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
This function computes the x, y and z momentum equations for the SV system. 
It is dimension-agnostic. 
*/
 
template<>
InputParameters validParams<EelMomentum>()
{
  InputParameters params = validParams<Kernel>();
    params.addRequiredCoupledVar("q_x", "x component of momentum");
    params.addCoupledVar("q_y", "y component of momentum");
    params.addCoupledVar("h"  , "h: water height");
//    params.addRequiredCoupledVar("pressure", "pressure");
//    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    params.addParam<int>("component", 0, "component of the momentum equation to compute (0,1)->(x,y)");
    params.addParam<RealVectorValue>("gravity", (0., 0., 0.), "Gravity vector.");
  return params;
}

SV_Momentum::SV_Momentum(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Coupled auxiliary variables
    _q_x(coupledValue("q_x")),
    _q_y(_mesh.dimension()>=2 ? coupledValue("q_y") : _zero),
    //_pressure(coupledValue("pressure")),
    // Equation of state:
    //_eos(getUserObject<EquationOfState>("eos")),
    _h(coupledValue("h")),
    // Parameters:
    _component(getParam<int>("component")),
    _gravity(  getParam<RealVectorValue>("gravity")),
    // Parameters for jacobian:
    _h_nb(coupled("h")),
    _q_x_nb(coupled("q_x")),
    _q_y_nb(isCoupled("q_y") ? coupled("q_y") : -1),


{
    if ( _component > 1 )
        mooseError("ERROR: the integer variable 'component' can only take values: 0 or 1 for the shallow water system!");
}

Real SV_Momentum::computeQpResidual()
{
  // vector q
  RealVectorValue _vector_q( _q_x[_qp], _q_y[_qp], 0. );
  // _u/h is one of the momentum components {qx, qy, or qz} divided by h
  if( _h[_qp] < 0. )
    mooseError("h < 0");

  // the "-" comes from the integration by parts
  RealVectorValue _advection = -_u[_qp]/_h[_qp] * _vector_q ;
        
  // Hydrostatic pressure term:
  // the "-" comes from the integration by parts
  // jcr_note: not sure about the gravity here. it may depend on the orientation
  Real _P = -0.5 * _gravity(_component) * std::pow(_h[_qp],2);
  //  Real _PdA = _pressure[_qp]*_grad_area[_qp](_component);

  // Source term: gravity * h * grad(B)
  // jcr_note: we could also integrate this term by parts. 
  Real _source_term = _gravity(_component) * _h[_qp] * grad_bathymetry[_qp](_component);
     
  // Return the kernel value (convention: lhs of the = sign):
  return ( _advection*_grad_test[_i][_qp] + _P*_grad_test[_i][_qp](_component) + _source_term*test[_i][_qp] ); 
}

Real SV_Momentum::computeQpJacobian()
{
  // Compute the momentum vector q:
  RealVectorValue _vector_q(_q_x[_qp], _q_y[_qp], 0.);

  // Compute the velocity vector:
  RealVectorValue _vector_vel = _vector_q / h[_qp];
  _vector_vel(_component) *= 2.;
    
  // Compute the derivative of hydrostatic pressure
  //Real _press_term = _eos.dAp_dq(_rhoA[_qp], _u[_qp], _rhoEA[_qp])*(...
  
  // Return the value of the jacobian:
  return -_phi[_j][_qp] * ( _vector_vel * _grad_test[_i][_qp] );
}

Real SV_Momentum::computeQpOffDiagJacobian( unsigned int _jvar)
{
  // Compute the momentum vector q:
  RealVectorValue _vector_q(_q_x[_qp], _q_y[_qp], 0.);

  // Compute the velocity vector:
  RealVectorValue _vector_vel = _vector_q / h[_qp];

  // density h:
  if (_jvar == _h_nb) {
    Real _Psrc_term = _gravity(_component) * (_h[_qp] - grad_bathymetry[_qp](_component));
    return _phi[_j][_qp] * (_u[_qp]/_h[_qp] * _vector_vel * _grad_test[_i][_qp] + _Psrc_term );
  }
    
  // x-momentum component:
  else if (_jvar == _q_x_nb ) {
   return -_phi[_j][_qp] * ( _grad_test[_i][_qp](0)*_u[_qp]/_h[_qp] );
  }
  // y-momentum component:
  else if (_jvar == _q_y_nb && _mesh.dimension()>=2 ) {
   return -_phi[_j][_qp] * ( _grad_test[_i][_qp](1)*_u[_qp]/_h[_qp] );
  }
  //
  else
    return 0.;
}
