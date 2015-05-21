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

#ifndef RECTANGLE_FUNCTION_1D_H
#define RECTANGLE_FUNCTION_1D_H

#include "Function.h"
#include "FunctionInterface.h" // jcr what is this for?

class RectangleFunction1D;

template<>
InputParameters validParams<RectangleFunction1D>();

class RectangleFunction1D : public Function
{
public:
  RectangleFunction1D(const std::string & name, InputParameters parameters);

  virtual Real value(Real t, const Point & p);

  virtual RealVectorValue gradient(Real t, const Point & p);

protected:
    
  // limits of the rectangle
  Real _xmin;
  Real _xmax;
  // values of the function
  Real _before_rect;
  Real _in_rect;
  Real _after_rect;
};

#endif //RECTANGLE_FUNCTION_1D_H
