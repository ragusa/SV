#ifndef JUMPINTERFACE_H
#define JUMPINTERFACE_H

#include "InternalSideUserObject.h"

// Forward Declarations
class JumpInterface;

template<>
InputParameters validParams<JumpInterface>();

class JumpInterface : public InternalSideUserObject
{
public:
  // Constructor
  JumpInterface(const std::string & name, InputParameters parameters);

  // Destructor
  virtual ~JumpInterface();

  virtual void initialize();
  virtual void execute();
  virtual void destroy();
  virtual void finalize();
  virtual void threadJoin(const UserObject & uo);

  Real getValue() const { return _value; }  

protected:
  // Auxiliary system variable:
  AuxiliarySystem & _aux;

  // Gradient value:
  VariableGradient & _grad_F_x;
  VariableGradient & _grad_F_x_neighbor;
  VariableGradient & _grad_F_y;
  VariableGradient & _grad_F_y_neighbor;

  // Name of the variable storing the jump:
  std::string _var_name;

  // Temporary variable:
  Real _value;
};

#endif // JUMPINTERFACE_H

