#ifndef SV_MOMENTUM_H
#define SV_MOMENTUM_H

#include "Kernel.h"
// We may want to use the hydrostatic expression and keep P here
// #include "EquationOfState.h". 

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

private:
    // Aux variables:
    VariableValue & _q_x;
    VariableValue & _q_y;
    // VariableValue & _pressure;
    VariableValue & _h;

    
    // Equation of state:
    // const EquationOfState & _eos;
    
    // Parameters:
    int _component;
    RealVectorValue _gravity;
    
    // Parameters for jacobian:
    unsigned int _h_nb;
    unsigned int _q_x_nb;
    unsigned int _q_y_nb;
};

#endif // SV_MOMENTUM_H
