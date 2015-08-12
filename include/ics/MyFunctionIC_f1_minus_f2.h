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

#ifndef MYFUNCTIONIC_H
#define MYFUNCTIONIC_H

#include "FunctionIC.h"

//Forward Declarations
class MyFunctionIC_f1_minus_f2;
class Function;

template<>
InputParameters validParams<FunctionIC>();

/**
 * Defines a boundary condition that forces the value to be a user specified
 * function at the boundary.
 */
class MyFunctionIC_f1_minus_f2 : public FunctionIC
{
public:
  MyFunctionIC_f1_minus_f2(const InputParameters & parameters);

protected:
  /**
   * Evaluate the function at the current quadrature point and timestep.
   */
  Real f();

  /**
   * The value of the variable at a point.
   */
  virtual Real value(const Point &p);

  /**
   * The value of the gradient at a point.
   */
  virtual RealGradient gradient(const Point &p);

  Function & _func1;
  Function & _func2;


};

#endif // MYFUNCTIONIC_H
