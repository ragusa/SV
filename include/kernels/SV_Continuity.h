#ifndef SV_CONTINUITY_H
#define SV_CONTINUITY_H

#include "Kernel.h"

class SV_Continuity;

template<>
InputParameters validParams<SV_Continuity>();
class SV_Continuity : public Kernel
{
public:

  SV_Continuity(const std::string & name,
                      InputParameters parameters);

protected:
 
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:
    // Coupled variables:
    VariableValue & _q_x;
    VariableValue & _q_y;
    
    // Parameters for jacobian (off-diagonal blocks):
    unsigned int _q_x_ivar;
    unsigned int _q_y_ivar;
};

#endif // SV_CONTINUITY_H
