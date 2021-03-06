/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                               */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "SV_Continuity.h"

/**
This Kernel computes the convection flux of the continuity equation :
div ( q ).
**/

template<>
InputParameters validParams<SV_Continuity>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("q_x", "x component of q");
  params.addCoupledVar("q_y", "y component of q");
  
  return params;
}

SV_Continuity::SV_Continuity(const InputParameters & parameters) :
  Kernel(parameters),
    // Coupled variables
    _q_x(coupledValue("q_x")),
    _q_y(_mesh.dimension() == 2 ? coupledValue("q_y") : _zero ),
    // Parameters for jacobian
    _q_x_ivar(coupled("q_x")),
    _q_y_ivar(_mesh.dimension() == 2  ? coupled("q_y") :0)
{}

Real SV_Continuity::computeQpResidual()
{
    // Compute inviscid flux for the continuity equation:
    RealVectorValue _q_vector(_q_x[_qp], _q_y[_qp], 0.);
    
    // Returns: - \vec{q} \cdot grad_test[i] (- sign comes from int.by.parts):
/*    if(_q_point[_qp](0) >-1. & _q_point[_qp](0) <0.)
    Moose::out << "qp="<<_q_point[_qp](0) << " Continuity:" << "\t h  ="  <<_u[_qp] 
                                                            << "\t qx =" << _q_x[_qp] 
                                                            << "\t RES="<<
                                -_q_vector * _grad_test[_i][_qp] << std::endl; */
    return -_q_vector * _grad_test[_i][_qp];
}

Real SV_Continuity::computeQpJacobian()
{
  // here we deal with the continuity row of the Jacobian matrix
  return 0.;
}

Real SV_Continuity::computeQpOffDiagJacobian( unsigned int _jvar)
{
    // here we deal with the continuity row of the Jacobian matrix (off-diag terms)
    // off-diag term is momentum in x
    if (_jvar == _q_x_ivar)
        return -_phi[_j][_qp]*_grad_test[_i][_qp](0);
    // off-diag term is momentum in y
    else if (_jvar == _q_y_ivar)
        return -_phi[_j][_qp]*_grad_test[_i][_qp](1);
    else
        return 0.;
}
