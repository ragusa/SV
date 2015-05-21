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

#ifndef PARTIAL_COSINE_FUNCTION_1D_H
#define PARTIAL_COSINE_FUNCTION_1D_H

#include "Function.h"
#include "FunctionInterface.h" // jcr what is this for?

class PartialCosineFunction1D;

template<>
InputParameters validParams<PartialCosineFunction1D>();

class PartialCosineFunction1D : public Function
{
public:
  PartialCosineFunction1D(const std::string & name, InputParameters parameters);

  virtual Real value(Real t, const Point & p);

  virtual RealVectorValue gradient(Real t, const Point & p);

protected:
    
  // limits of the cosine hump
  Real _xmin;
  Real _xmax;

};

#endif //PARTIAL_COSINE_FUNCTION_1D_H
