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

#ifndef BATHYMETRY1D_STEP_H
#define BATHYMETRY1D_STEP_H

#include "Function.h"
#include "FunctionInterface.h" // jcr what is this for?

class Bathymetry1D_step;

template<>
InputParameters validParams<Bathymetry1D_step>();

class Bathymetry1D_step : public Function
{
public:
  Bathymetry1D_step(const std::string & name, InputParameters parameters);

  virtual Real value(Real t, const Point & p);

  virtual RealVectorValue gradient(Real t, const Point & p);

protected:
    
    // left/right value of a step in bathymetry
    Real _left;
    Real _right;
    Real _membrane; // location of the membrane

};

#endif //BATHYMETRY1D_STEP_H
