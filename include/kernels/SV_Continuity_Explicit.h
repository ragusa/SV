#ifndef SV_CONTINUITY_EXPLICIT_H
#define SV_CONTINUITY_EXPLICIT_H

#include "Kernel.h"

class SV_Continuity_Expl;

template<>
InputParameters validParams<SV_Continuity_Expl>();
class SV_Continuity_Expl : public Kernel
{
public:

  SV_Continuity_Expl(const InputParameters parameters);

protected:
 
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:
  // Coupled variables:
  VariableValue & _q_x;
  VariableValue & _q_y;

};

#endif // SV_CONTINUITY_EXPLICIT_H
