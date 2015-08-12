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

#include "SV_Continuity_Explicit.h"

/**
This Kernel computes the convection flux of the continuity equation :
div ( q ).
**/

template<>
InputParameters validParams<SV_Continuity_Expl>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("q_x", "x component of q");
  params.addCoupledVar("q_y", "y component of q");
  
  return params;
}

SV_Continuity_Expl::SV_Continuity_Expl(const InputParameters parameters) :
  Kernel(parameters),
  // Coupled variables
  _q_x(coupledValueOld("q_x")),
  _q_y(_mesh.dimension() == 2 ? coupledValueOld("q_y") : _zero )
{}

Real SV_Continuity_Expl::computeQpResidual()
{
  // Compute inviscid flux for the continuity equation:
  RealVectorValue _q_vector(_q_x[_qp], _q_y[_qp], 0.);
  
  // Returns: - \vec{q} \cdot grad_test[i] (- sign comes from int.by.parts):
  return -_q_vector * _grad_test[_i][_qp];
}

Real SV_Continuity_Expl::computeQpJacobian()
{
  // it's explicit (and even for implicit, this diagonal block is 0)
  return 0.;
}

Real SV_Continuity_Expl::computeQpOffDiagJacobian( unsigned int _jvar)
{
  // it's explicit
  return 0.;
}
