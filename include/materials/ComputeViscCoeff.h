#ifndef COMPUTEVISCCOEFF_H
#define COMPUTEVISCCOEFF_H

#include "Material.h"
#include "MaterialProperty.h"
// #include "EquationOfState.h"

//Forward Declarations
class ComputeViscCoeff;

template<>
InputParameters validParams<ComputeViscCoeff>();

class ComputeViscCoeff : public Material
{
public:
  ComputeViscCoeff(const std::string & name, InputParameters parameters);

protected:
  virtual void computeQpProperties();

private:
    // Viscosity types
    enum ViscosityType
    {
        NONE = 0,
        FIRST_ORDER = 1,
        ENTROPY = 2
    };
    std::string _visc_name;
    MooseEnum _visc_type;
    
    
    // Boolean for jump
    bool _isJumpOn;
    bool _isShock;
    
    // Coupled aux variables: velocity
    VariableValue & _vel_x;
    VariableValue & _vel_y;
    VariableValue & _vel_z;
    VariableGradient & _grad_vel_x;
    VariableGradient & _grad_vel_y;
    VariableGradient & _grad_vel_z;
    
    // Coupled aux variables: entropy pair
    VariableValue    & _s;
    VariableValue    & _s_dot;
    VariableValue    & _s_old;
    VariableValue    & _s_older;
    VariableGradient & _grad_psi;
    VariableGradient & _grad_psi_old;

    
    // Coupled aux variable: norm of velocity
    VariableValue & _norm_vel;
    VariableGradient & _grad_norm_vel;
    
    // Jump of entropy:
    VariableValue & _jump_grad_s;
        
    // Material properties
    MaterialProperty<Real> & _mu;
    MaterialProperty<Real> & _mu_max;
    MaterialProperty<Real> & _kappa;
    MaterialProperty<Real> & _kappa_max;
    MaterialProperty<RealVectorValue> & _l;
    
    
    // Multiplicative coefficient for viscosity:
    double _Ce;
    double _Cjump;
    double _Cmax;
    
    // UserObject: equation of state
    // const EquationOfState & _eos;
    
    // Name of the postprocessor:
    std::string _entropy_pps_name;

};

#endif // ComputeViscCoeff_H
