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

#ifndef SV_TIMEDERIVATIVE
#define SV_TIMEDERIVATIVE

#include "TimeDerivative.h"
#include "Function.h"

class SVTimeDerivative;

template<>
InputParameters validParams<SVTimeDerivative>();

class SVTimeDerivative : public TimeDerivative
{
public:

  SVTimeDerivative(const std::string & name,
                        InputParameters parameters);

protected:
    virtual Real computeQpResidual();

    virtual Real computeQpJacobian();

};

#endif // SV_TIMEDERIVATIVE
