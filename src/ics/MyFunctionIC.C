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

#include "MyFunctionIC.h"
#include "Function.h"

template<>
InputParameters validParams<MyFunctionIC>()
{
  InputParameters params = validParams<FunctionIC>();
  params.addRequiredParam<FunctionName>("function", "The function needed to define the initial condition.");
  params.addRequiredParam<Real>("real_input", "Real value needed to define the IC function");
  
  return params;
}

MyFunctionIC::MyFunctionIC(const InputParameters & parameters) :
    FunctionIC(parameters),
    _func(getFunction("function")),
    _input_real_value(getParam<Real>("real_input"))
{
}

Real
MyFunctionIC::value(const Point & p)
{
  return _input_real_value - _func.value(_t, p);
}

RealGradient
MyFunctionIC::gradient(const Point & p)
{
  return -1.0*_func.gradient(_t, p);
}
