#ifndef SVSETWATERVELOCITY_H
#define SVSETWATERVELOCITY_H

#include "IntegratedBC.h"
#include "Function.h"

// Forward Declarations
class SVSetWaterVelocity;
class HydrostaticPressure;

template<>
InputParameters validParams<SVSetWaterVelocity>();

class SVSetWaterVelocity : public IntegratedBC
{
public:
  SVSetWaterVelocity(const InputParameters & parameters);
  virtual ~SVSetWaterVelocity(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Equation type
  enum EquationType
  {
    CONTINUITY = 0,
    X_MOMENTUM = 1,
    Y_MOMENTUM = 2
  };
  MooseEnum _equ_type;

  // Coupled variables
  VariableValue & _h;

  // Constants and parameters
  Real _u_bc;

  // Equation of state
  const HydrostaticPressure & _eos;

  // Integers for jacobian terms
  unsigned _h_var;
  unsigned _q_x_var;
};

#endif // SVSETWATERVELOCITY_H

