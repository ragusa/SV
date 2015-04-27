#ifndef SV_MOMENTUM_H
#define SV_MOMENTUM_H

#include "Kernel.h"
// We may want to use the hydrostatic expression and keep P here
#include "HydrostaticPressure.h" 

// Forward Declarations
class SV_Momentum;

template<>
InputParameters validParams<SV_Momentum>();

class SV_Momentum : public Kernel
{
public:

  SV_Momentum(const std::string & name,
              InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private: // jcr protected?, no
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

  // Parameters for jacobian:
  unsigned int _h_ivar;
  unsigned int _q_x_ivar;
  unsigned int _q_y_ivar;
};

#endif // SV_MOMENTUM_H
