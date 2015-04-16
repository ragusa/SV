#ifndef SOLIDWALLBC_H
#define SOLIDWALLBC_H

#include "IntegratedBC.h"

// Forward Declarations
class SolidWallBC;
class HydrostaticPressure;

template<>
InputParameters validParams<SolidWallBC>();

class SolidWallBC : public IntegratedBC
{
public:
  SolidWallBC(const std::string & name, InputParameters parameters);
  virtual ~SolidWallBC(){}

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

  // Equation of state
  const HydrostaticPressure & _eos;

  // Integers for jacobian terms
  unsigned _h_var;
  unsigned _q_x_var;
  unsigned _q_y_var;
};

#endif // SOLIDWALLBC_H

