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
/**
This function aims at computing the bathymetry data at the nodes and its gradient. This auxiliary variable is coupled
to rho_bar, m_bar and E_bar defined as the product of the usual density, momentum and energy, and the cross section
A computed by the function AreaFunction.
**/
#include "Bathymetry_Aux.h"

template<>
InputParameters validParams<Bathymetry_Aux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<FunctionName>("area", "function to compute the cross section");
  return params;
}

Bathymetry_Aux::Bathymetry_Aux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    _bathymetry(getFunction("bathymetry"))
{}

/**
 * Auxiliary Kernels override computeValue() instead of computeQpResidual().  Aux Variables
 * are calculated either one per element or one per node depending on whether we declare
 * them as "Elemental (Constant Monomial)" or "Nodal (First Lagrange)".  No changes to the
 * source are necessary to switch from one type or the other.
 */
Real Bathymetry_Aux::computeValue()
{
  Point _node_coord = (*_current_node);
  return _bathymetry.value(0.0,_node_coord); 
}
