#ifndef SVSETWATERVELOCITYINLETBC_H
#define SVSETWATERVELOCITYINLETBC_H

#include "IntegratedBC.h"
#include "Function.h"

// Forward Declarations
class SVSetWaterVelocityInletBC;
class HydrostaticPressure;

template<>
InputParameters validParams<SVSetWaterVelocityInletBC>();

class SVSetWaterVelocityInletBC : public IntegratedBC
{
public:
  SVSetWaterVelocityInletBC(const InputParameters & parameters);
  virtual ~SVSetWaterVelocityInletBC(){}

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
  Real _h_bc;

  // Equation of state
  const HydrostaticPressure & _eos;

  // Integers for jacobian terms
  unsigned _h_var;
  unsigned _q_x_var;

  // boolean
  bool _is_h_bc_specified;
};

#endif // SVSETWATERVELOCITYINLETBC_H

