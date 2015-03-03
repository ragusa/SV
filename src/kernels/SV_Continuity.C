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

SV_Continuity::SV_Continuity(const std::string & name,
                             InputParameters parameters) :
  Kernel(name, parameters),
    // Coupled aux variables jcr: why aux?
    _q_x(coupledValue("q_x")),
    _q_y(_mesh.dimension() == 2 ? coupledValue("q_y") : _zero ),
    // Parameters for jacobian
    _q_x_ivar(coupled("q_x")),
    // jcr: another way
    // _q_y_ivar(isCoupled("q_y") ? coupled("q_y") : 0)
    _q_y_ivar(_mesh.dimension() == 2  ? coupled("q_y") :0)
{}

Real SV_Continuity::computeQpResidual()
{
    // Compute inviscid flux for the continuity equation:
    RealVectorValue _q_vector(_q_x[_qp], _q_y[_qp], 0.);
    
    // Returns: - \vec{q} \cdot grad_test[i] (- sign comes from int.by.parts):
    return -_q_vector * _grad_test[_i][_qp];
}

Real SV_Continuity::computeQpJacobian()
{
  return 0.;
}

Real SV_Continuity::computeQpOffDiagJacobian( unsigned int _jvar)
{
    if (_jvar == _q_x_ivar)
        return -_phi[_j][_qp]*_grad_test[_i][_qp](0);
    else if (_jvar == _q_y_ivar)
        return -_phi[_j][_qp]*_grad_test[_i][_qp](1);
    else
        return 0.;
}
