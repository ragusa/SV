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

#include "MyFunctionIC_f1_minus_f2.h"
#include "Function.h"

template<>
InputParameters validParams<MyFunctionIC_f1_minus_f2>()
{
  InputParameters params = validParams<FunctionIC>();
  params.addRequiredParam<FunctionName>("function1", "The 1st function needed to define the initial condition.");
  params.addRequiredParam<FunctionName>("function2", "The 2nd function needed to define the initial condition.");
  
  return params;
}

MyFunctionIC_f1_minus_f2::MyFunctionIC_f1_minus_f2(const InputParameters & parameters) :
    FunctionIC(parameters),
    _func1(getFunction("function1")),
    _func2(getFunction("function2"))
{
}

Real
MyFunctionIC_f1_minus_f2::value(const Point & p)
{
  return _func1.value(_t, p) - _func2.value(_t, p);
}

RealGradient
MyFunctionIC_f1_minus_f2::gradient(const Point & p)
{
  return _func1.gradient(_t, p) * _func2.gradient(_t, p);
}
