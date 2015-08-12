#ifndef SV_MOMENTUM_EXPLICIT_H
#define SV_MOMENTUM_EXPLICIT_H

#include "Kernel.h"
// We may want to use the hydrostatic expression and keep P here
#include "HydrostaticPressure.h" 

// Forward Declarations
class SV_Momentum_Expl;

template<>
InputParameters validParams<SV_Momentum_Expl>();

class SV_Momentum_Expl : public Kernel
{
public:

  SV_Momentum_Expl(const InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:
  //  variables:
  VariableValue & _h;
  VariableValue & _q_x;
  VariableValue & _q_y;

  VariableGradient & _grad_bathymetry;

  // Equation of state:
  const HydrostaticPressure & _eos;

  // Parameters:
  unsigned int _component;
  Real _gravity;
};

#endif // SV_MOMENTUM_EXPLICIT_H
