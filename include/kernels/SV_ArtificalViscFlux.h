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

  SV_ArtificialViscFlux(const std::string & name,
                   InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int _jvar);
    
private:
    // Equations types
    enum EquationType
    {
        CONTINUITY = 0,
        XMOMENTUM = 1,
        YMOMENTUM = 2
    };
    // Diffusion name
    std::string _equ_name;
    // Diffusion type
    MooseEnum _equ_type;

    // Material property: viscosity coefficient.
    MaterialProperty<Real> & _kappa;
};

#endif // SV_ARTIFICIALVISCFLUX_H
