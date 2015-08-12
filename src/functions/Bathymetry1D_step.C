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

#include "Bathymetry1D_step.h"

template<>
InputParameters validParams<Bathymetry1D_step>()
{
  InputParameters params = validParams<Function>();
  params.addParam<Real>("left"    , 0., "Value to the left of the membrane" );
  params.addParam<Real>("right"   , 1., "Value to the right of the membrane");
  params.addParam<Real>("membrane",0.5, "Location of membrane");
  return params;
}

Bathymetry1D_step::Bathymetry1D_step(const InputParameters & parameters) :
  Function(parameters),
  _left(getParam<Real>("left")),
  _right(getParam<Real>("right")),
  _membrane(getParam<Real>("membrane"))
{}

Real
Bathymetry1D_step::value(Real /*t*/, const Point & p)
{
  if ( p(0)<_membrane )
    return _left;
  else 
    return _right;
}

RealVectorValue 
Bathymetry1D_step::gradient(Real /*t*/, const Point & p)
{
  if ( p(0)<_membrane )
    return 0.;
  else 
    return 0.;
}
