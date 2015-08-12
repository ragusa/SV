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

#include "RectangleFunction1D.h"

template<> 
InputParameters validParams<RectangleFunction1D>()
{
  InputParameters params = validParams<Function>();
  params.addParam<Real>("xmin", 0.4, "Beginning of the rectangle" );
  params.addParam<Real>("xmax", 0.6, "End of the rectangle" );
  params.addParam<Real>("value_before", 0., "Value of the function before the rectangle" );
  params.addParam<Real>("value_rect"  , 1., "Value of the function in the rectangle" );
  params.addParam<Real>("value_after" , 2., "Value of the function after the rectangle" );
  return params;
}

RectangleFunction1D::RectangleFunction1D(const InputParameters & parameters) :
  Function(parameters),
    _xmin(getParam<Real>("xmin")),
    _xmax(getParam<Real>("xmax")),
    _before_rect(getParam<Real>("value_before")),
    _in_rect(getParam<Real>("value_rect")),
    _after_rect(getParam<Real>("value_after"))
{
}

Real
RectangleFunction1D::value(Real /*t*/, const Point & p)
{
  if ( p(0)<_xmin )
    return _before_rect;
  else if ( p(0)>_xmax )
    return _after_rect;
  else
    return _in_rect;
}

RealVectorValue 
RectangleFunction1D::gradient(Real /*t*/, const Point & p)
{
  return 0.;
}
