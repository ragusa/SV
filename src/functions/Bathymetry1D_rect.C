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

#include "Bathymetry1D_rect.h"

template<>
InputParameters validParams<Bathymetry1D_rect>()
{
  InputParameters params = validParams<Function>();
  params.addParam<Real>("domain_length", 1500., "Length of the domain" );
  params.addParam<Real>("step_height", 8., "Height of the step" );
  return params;
}

Bathymetry1D_rect::Bathymetry1D_rect(const InputParameters & parameters) :
  Function(parameters),
  _length(getParam<Real>("domain_length")),
  _height(getParam<Real>("step_height"))
{
  if (_height <= 0.)
    mooseError("'" << this->name() << "' : height of step cannot be <=0.");
}

Real
Bathymetry1D_rect::value(Real /*t*/, const Point & p)
{
  if ( std::abs(p(0)-_length/2) < _length/_height )
    return _height;
  else 
    return 0.;
}

RealVectorValue 
Bathymetry1D_rect::gradient(Real /*t*/, const Point & p)
{
  return 0.;
}
