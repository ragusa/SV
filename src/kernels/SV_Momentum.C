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
This function computes the x and y momentum equations for the SV system. 
It is dimension-agnostic. 
**/
 
template<>
InputParameters validParams<SV_Momentum>()
{
  InputParameters params = validParams<Kernel>();
  // Coupled variables
  params.addRequiredCoupledVar("h"  , "h: water height");
  params.addRequiredCoupledVar("q_x", "x component of momentum");
  params.addCoupledVar(        "q_y", "y component of momentum");
  // Coupled aux variables
  params.addCoupledVar("B", "bathymetry data");// jcr how about addRequiredCoupledVar instead
  // Constants and parameters
  params.addRequiredParam<Real>("gravity", "Gravity magnitude");
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  params.addRequiredParam<unsigned int>("component", "component of the momentum equation to compute (0,1)->(x,y)");
    
  return params;
}

SV_Momentum::SV_Momentum(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
  // Coupled variables
  _h(coupledValue("h")),
  _q_x(coupledValue("q_x")),
	// the test with is_implicit is not needed. just use implicit=T/F in the kernel block itself!
	// _q_y( _mesh.dimension() == 2 ? ( _is_implicit ? coupledValue("q_y") : coupledValueOld("q_y") ) : _zero),
  _q_y(_mesh.dimension() == 2 ? coupledValue("q_y") : _zero),
  // Coupled auxiliary variables
  _grad_bathymetry(isCoupled("B") ? coupledGradient("B") : _grad_zero),
  //_pressure(coupledValue("pressure")),
  // Equation of state:
  _eos(getUserObject<HydrostaticPressure>("eos")),
  // Parameters:
  _component(getParam<unsigned int>("component")),
  _gravity(getParam<Real>("gravity")),
  // Parameters for jacobian:
  _h_ivar(coupled("h")),
  _q_x_ivar(coupled("q_x")),
  _q_y_ivar(_mesh.dimension() == 2 ? coupled("q_y") : 0)

{
    if ( _component > 1 )
        mooseError("ERROR in "<<this->name()<<": the integer variable 'component' can only take values: 0 or 1 for the shallow water system!");
}

Real SV_Momentum::computeQpResidual()
{
  if ( !_is_implicit )
  	return 0.;  // nothing to do for explicit time integration
  	
  // vector q
  RealVectorValue _vector_q( _q_x[_qp], _q_y[_qp], 0. );
  if( _h[_qp] < 0. ){
    Moose::out << "NEGATIVE: " << _q_point[0] << " " << _h[_qp] << std::endl;
    //mooseError("ERROR in "<<this->name()<<":h < 0");
	}

	
  // advection term q^2/h, integrated by parts:
  // component 0: -\int (qxqx/h dbdx + qxqy/h dbdy) = -\int qx/h \vec{q} \vec{grad}b
  // component 1: -\int (qxqy/h dbdx + qyqy/h dbdy) = -\int qy/h \vec{q} \vec{grad}b
  
  // caveat: _u: is the current Moose variable. Here, it is q_x or q_y, depending on the component  
  // _u/h is one of the momentum components {qx or qy} divided by h
  RealVectorValue _advection = _u[_qp]/_h[_qp] * _vector_q ;
        
  // Hydrostatic pressure term, integrate by parts:
  // component 0: -\int P dbdx
  // component 1: -\int P dbdy
  // old: Real _P = -0.5 * _gravity(_component) * std::pow(_h[_qp],2);
  Real _P = _eos.pressure(_h[_qp], _vector_q);

  // Source term (placed on the LHS): gravity * h * grad(B)
  Real _source_term = _gravity * _h[_qp] * _grad_bathymetry[_qp](_component);
     
  // Return the kernel value (convention: LHS of the = sign):
  // the "-" in front of advection and pressure comes from the integration by parts
  return ( -_advection*_grad_test[_i][_qp] - _P*_grad_test[_i][_qp](_component) + _source_term*_test[_i][_qp] ); 
}

Real SV_Momentum::computeQpJacobian()
{
  if ( !_is_implicit )
  	return 0.;  // nothing to do for explicit time integration

  // derivative of component 0 with respect to component 0:
  //  -\int (2qx/h dbdx + qy/h dbdy) 
  // derivative of component 1 with respect to component 1:
  //  -\int (qx/h dbdx + 2qy/h dbdy) 

  // Compute the momentum vector q:
  RealVectorValue _vector_q(_q_x[_qp], _q_y[_qp], 0.);

  // Compute the velocity vector:
  RealVectorValue _vector_vel = _vector_q / _h[_qp];
  _vector_vel(_component) *= 2.;
  
  // Compute the derivative of the advection term
  Real  d_advection_du = _phi[_j][_qp] * ( _vector_vel * _grad_test[_i][_qp] );
  // Compute the derivative of hydrostatic pressure
  Real d_pressure_du = _phi[_j][_qp]*_grad_test[_i][_qp](_component);
  if (_component == 0)
    d_pressure_du *= _eos.dp_dqx(_h[_qp], _vector_q);
  else if (_component == 1)
    d_pressure_du *= _eos.dp_dqy(_h[_qp], _vector_q);
  else 
    mooseError("ERROR in "<<this->name()<<": the integer variable 'component' can only take values 0 or 1 for the shallow water system!");

  // Return the value of the jacobian:
  // the "-" signs come from IBP
  return -d_advection_du -d_pressure_du;
}

Real SV_Momentum::computeQpOffDiagJacobian(unsigned int _jvar)
{
  if ( !_is_implicit )
  	return 0.;  // nothing to do for explicit time integration

  // Compute the momentum vector q:
  RealVectorValue _vector_q(_q_x[_qp], _q_y[_qp], 0.);

  // Compute the velocity vector:
  RealVectorValue _vector_vel = _vector_q / _h[_qp];

  // density h:
  // 
  if (_jvar == _h_ivar) 
  {
	// advection off-diag derivative (caveat in comment: qvec/h is later written as vvec)
	// momentum-0: -d(qx/h qvec.gradb)/dh = -(-qx/h) vvec.gradb
	// momentum-1: -d(qy/h qvec.gradb)/dh = -(-qy/h) vvec.gradb
	Real d_advection_dh = -_phi[_j][_qp] * _u[_qp]/_h[_qp] * (_vector_vel * _grad_test[_i][_qp]);
	// pressure off-diag derivative: -d(0.5gh^2 gradb)/dh = -gh gradb
	Real d_pressure_dh = _phi[_j][_qp]*_eos.dp_dh(_h[_qp], _vector_q)*_grad_test[_i][_qp](_component);
	// bathymetry off-diag derivative: d(gh gradB b)/dh = g gradB b
    Real d_bathy_dh  = _phi[_j][_qp] * _gravity *_grad_bathymetry[_qp](_component) *_test[_i][_qp];
    return -d_advection_dh -d_pressure_dh +d_bathy_dh;
    // _phi[_j][_qp] * (_u[_qp]/_h[_qp] * _vector_vel * _grad_test[_i][_qp] + _Psrc_term );
  }
    
  // derivative wrt x-momentum component (thus, this is the q_y equation)
  else if (_jvar == _q_x_ivar ) 
  {
    return -_phi[_j][_qp] * ( _grad_test[_i][_qp](0)*_u[_qp]/_h[_qp] );
  }
  // y-momentum component:
  else if (_jvar == _q_y_ivar ) 
  {
    return -_phi[_j][_qp] * ( _grad_test[_i][_qp](1)*_u[_qp]/_h[_qp] );
  }
  //
  else
    return 0.;
}
