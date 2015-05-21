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

#ifndef BATHYMETRY1D_RECT_H
#define BATHYMETRY1D_RECT_H

#include "Function.h"
// #include "FunctionInterface.h" // jcr what is this for?

class Bathymetry1D_rect;

template<>
InputParameters validParams<Bathymetry1D_rect>();

class Bathymetry1D_rect : public Function
{
public:
  Bathymetry1D_rect(const std::string & name, InputParameters parameters);

  virtual Real value(Real t, const Point & p);

  virtual RealVectorValue gradient(Real t, const Point & p);

protected:
    
    // domain length
    Real _length;
    Real _height;

};

#endif //BATHYMETRY1D_RECT_H
