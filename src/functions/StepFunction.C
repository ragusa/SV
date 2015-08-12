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

#include "StepFunction.h"
// #include "MooseMesh.h"
//#include "MooseMesh.h"
//#include "FEProblem.h"
//#include "Coupleable.h"

template<>
InputParameters validParams<StepFunction>()
{
  InputParameters params = validParams<Function>();
  params.addParam<Real>("value_before_step"    , 0., "Value before the step" );
  params.addParam<Real>("value_after_step"     , 1., "Value after the step");
  params.addParam<Real>("x_step",0.5, "x-location of step");
  params.addParam<Real>("y_step",0.5, "y-location of step");
  params.addParam<Real>("z_step",0.5, "z-location of step");
  return params;
}

StepFunction::StepFunction(const InputParameters & parameters) :
  Function(parameters),
//  _mesh(_sc_fe_problem.mesh()),
//  _mesh(*parameters.get<MooseMesh *>("mesh")),
  _dim(_sc_fe_problem.mesh().dimension()),
  _value_before_step(getParam<Real>("value_before_step")),
  _value_after_step(getParam<Real>("value_after_step")),
  _x0(getParam<Real>("x_step")),
  _y0(getParam<Real>("y_step")),
  _z0(getParam<Real>("z_step"))
{}

Real
StepFunction::value(Real /*t*/, const Point & p)
{
//  _mesh=_sc_fe_problem.mesh();
//  MooseMesh & _mesh=_sc_fe_problem.mesh();
//  switch( _mesh.dimension() ) 
  switch( _dim ) 
  {
  case 1:
    if ( p(0)<_x0 )
      return _value_before_step;
    else 
      return _value_after_step;
    break;
  case 2:
    if ( p(0)<_x0 && p(1)<_y0 )
      return _value_before_step;
    else 
      return _value_after_step;
    break;
  case 3:
    if ( p(0)<_x0 && p(1)<_y0 && p(2)<_z0 )
      return _value_before_step;
    else 
      return _value_after_step;
    break;
  default:
    mooseError("ERROR in "<<this->name()<<": wrong mesh dimension.");
    
  } 
}

RealVectorValue 
StepFunction::gradient(Real /*t*/, const Point & p)
{
    return RealVectorValue(0, 0, 0);
}
