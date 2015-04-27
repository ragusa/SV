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

#include "SV_Momentum_Explicit.h"

/**
This function computes the x and y momentum equations for the SV system. 
It is dimension-agnostic. 
**/
 
template<>
InputParameters validParams<SV_Momentum_Expl>()
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

SV_Momentum_Expl::SV_Momentum_Expl(const std::string & name,
                                   InputParameters parameters) :
  Kernel(name, parameters),
  // Coupled variables
  _h(coupledValueOld("h")),
  _q_x(coupledValueOld("q_x")),
  _q_y(_mesh.dimension() == 2 ? coupledValueOld("q_y") : _zero),
  // Coupled auxiliary variables
  _grad_bathymetry(isCoupled("B") ? coupledGradient("B") : _grad_zero),
  //_pressure(coupledValue("pressure")),
  // Equation of state:
  _eos(getUserObject<HydrostaticPressure>("eos")),
  // Parameters:
  _component(getParam<unsigned int>("component")),
  _gravity(getParam<Real>("gravity"))
{
  if ( _component > 1 )
    mooseError("ERROR: the integer variable 'component' can only take values: 0 or 1 for the shallow water system!");
}

Real SV_Momentum_Expl::computeQpResidual()
{
  // vector q
  RealVectorValue _vector_q( _q_x[_qp], _q_y[_qp], 0. );
  if( _h[_qp] < 0. )
    mooseError("h < 0");

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

Real SV_Momentum_Expl::computeQpJacobian()
{
  return 0.;
}

Real SV_Momentum_Expl::computeQpOffDiagJacobian(unsigned int _jvar)
{
  return 0.;
}
