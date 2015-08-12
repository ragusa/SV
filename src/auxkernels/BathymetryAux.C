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
This function aims at computing the bathymetry data.
**/
#include "BathymetryAux.h"  // jcr eliminate and use FunctionAux in input file

template<>
InputParameters validParams<BathymetryAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<FunctionName>("Bathymetry", "function to compute the bathymetry");
  return params;
}

BathymetryAux::BathymetryAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _bathymetry(getFunction("bathymetry"))
{}

/**
 * Auxiliary Kernels override computeValue() instead of computeQpResidual().  Aux Variables
 * are calculated either one per element or one per node depending on whether we declare
 * them as "Elemental (Constant Monomial)" or "Nodal (First Lagrange)".  No changes to the
 * source are necessary to switch from one type or the other.
 */
Real BathymetryAux::computeValue()
{
  Point _node_coord = (*_current_node);
  return _bathymetry.value(0.0,_node_coord); 
}
