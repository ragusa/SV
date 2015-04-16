#ifndef SVSETWATERHEIGHTINLETBC_H
#define SVSETWATERHEIGHTINLETBC_H

#include "IntegratedBC.h"
#include "Function.h"

// Forward Declarations
class SVSetWaterHeightInletBC;
class HydrostaticPressure;

template<>
InputParameters validParams<SVSetWaterHeightInletBC>();


/**
**/
class SVSetWaterHeightInletBC : public IntegratedBC
{
public:
  SVSetWaterHeightInletBC(const std::string & name, InputParameters parameters);
  virtual ~SVSetWaterHeightInletBC(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Equation type
  enum EquationType
  enum EquationType
  {
    CONTINUITY = 0,
    X_MOMENTUM = 1,
    Y_MOMENTUM = 2
  };
  MooseEnum _equ_type;

  // Coupled variables
  VariableValue & _q_x;

  // Constants and parameters
  Real _h_bc;
  Real _u_bc;

  // Equation of state
  const HydrostaticPressure & _eos;

  // Integers for jacobian terms
  unsigned _q_x_var;

  // boolean
  bool _u_bc_specified;
};

#endif // SVSETWATERHEIGHTINLETBC_H

