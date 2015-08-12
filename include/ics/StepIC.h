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

#ifndef STEPIC_H
#define STEPIC_H

// MOOSE Includes
#include "InitialCondition.h"

// Forward Declarations
class StepIC;

template<>
InputParameters validParams<StepIC>();

class StepIC : public InitialCondition
{
public:

  StepIC(const InputParameters & parameters);

  virtual Real value(const Point & p);

private:
  // Initial conditions for left and right values:
  Real _h_left;
  Real _h_right;
  RealVectorValue _U_left;
  RealVectorValue _U_right;

  // Position of the membrane:
  Real _x_membrane;
  Real _y_membrane;
  Real _radius;
};

#endif //STEPIC_H
