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

#ifndef SV_ARTIFICIALVISCFLUX_H
#define SV_ARTIFICIALVISCFLUX_H

#include "Kernel.h"

// Forward Declarations
class SV_ArtificialViscFlux;

template<>
InputParameters validParams<SV_ArtificialViscFlux>();

class SV_ArtificialViscFlux : public Kernel
{
public:

  SV_ArtificialViscFlux(const InputParameters & parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int _jvar);
    
private:
  // Equations types
  enum EquationType
  {
    CONTINUITY = 0,
    X_MOMENTUM = 1,
    Y_MOMENTUM = 2
  };
  // Equations type
  MooseEnum _equ_type;

  // Material property: viscosity coefficient.
  const MaterialProperty<Real> & _kappa;
};

#endif // SV_ARTIFICIALVISCFLUX_H
