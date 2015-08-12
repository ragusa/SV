#ifndef SVSETWATERHEIGHTOUTLETBC_H
#define SVSETWATERHEIGHTOUTLETBC_H

#include "IntegratedBC.h"
#include "Function.h"

// Forward Declarations
class SVSetWaterHeightOutletBC;
class HydrostaticPressure;

template<>
InputParameters validParams<SVSetWaterHeightOutletBC>();


/**
**/
class SVSetWaterHeightOutletBC : public IntegratedBC
{
public:
  SVSetWaterHeightOutletBC(const InputParameters & parameters);
  virtual ~SVSetWaterHeightOutletBC(){}

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
  VariableValue & _q_x;

  // Constants and parameters
  Real _h_bc;

  // Equation of state
  const HydrostaticPressure & _eos;

  // Integers for jacobian terms
  unsigned _h_var;  
  unsigned _q_x_var;
};

#endif // SVSETWATERHEIGHTOUTLETBC_H

