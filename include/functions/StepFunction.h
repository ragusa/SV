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

#ifndef STEPFUNCTION_H
#define STEPFUNCTION_H

#include "Function.h"
#include "Coupleable.h"
#include "FEProblem.h"
#include "MooseMesh.h"

// #include "FunctionInterface.h" // jcr what is this for?

class StepFunction;

template<>
InputParameters validParams<StepFunction>();

class StepFunction : public Function
{
public:
  StepFunction(const std::string & name, InputParameters parameters);

  virtual Real value(Real t, const Point & p);

  virtual RealVectorValue gradient(Real t, const Point & p);

protected:
    
//  MooseMesh & _mesh=_sc_fe_problem.mesh();
//  MooseMesh & _mesh;

  // values for a step function
  Real _value_before_step;
  Real _value_after_step;
  Real _x0; // location of the step
  Real _y0; // location of the step
  Real _z0; // location of the step

};

#endif //STEPFUNCTION_H
