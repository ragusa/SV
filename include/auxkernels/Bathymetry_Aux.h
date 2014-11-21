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
#include "AreaFunction.h" //jcr note" ?

//Forward Declarations
class Bathymetry_Aux;

template<>
InputParameters validParams<Bathymetry_Aux>();

/**
 * Coupled auxiliary value
 */
class Bathymetry_Aux : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  Bathymetry_Aux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();
  Function & _bathy;
};

#endif // BATHYMETRY_AUX_H
