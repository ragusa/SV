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

#ifndef BATHYMETRY_AUX_H
#define BATHYMETRY_AUX_H

#include "AuxKernel.h"
#include "Bathymetry1D_step.h" //jcr note" ?

//Forward Declarations
class BathymetryAux;

template<>
InputParameters validParams<BathymetryAux>();

/**
 * Coupled auxiliary value
 */
class BathymetryAux : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  BathymetryAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  Function & _bathymetry;
};

#endif // BATHYMETRY_AUX_H
