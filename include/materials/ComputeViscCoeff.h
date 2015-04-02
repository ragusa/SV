#ifndef COMPUTEVISCCOEFF_H
#define COMPUTEVISCCOEFF_H

#include "Material.h"
#include "MaterialProperty.h"  // jcr_new why no longer
#include "HydrostaticPressure.h"

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
        NONE        = 0,
        FIRST_ORDER = 1,
        ENTROPY     = 2
    };
    MooseEnum _visc_type;
      
    // Coupled variables:
    VariableValue & _h;
    VariableValue & _q_x;
    VariableValue & _q_y;
    
    // Coupled aux variables: entropy pair
    VariableValue    & _E;
    VariableValue    & _E_old;
    VariableValue    & _E_older;
    //VariableValue    & _E_dot;
    VariableGradient & _grad_F;
    //VariableGradient & _grad_F_old;
    VariableGradient & _grad_G;
    //VariableGradient & _grad_G_old;

        // Coupled aux variable: bathymetry
    VariableValue & _bathymetry;
    Real _gravity;
    
    // Jump of entropy:
    //VariableValue & _jump_grad_s;
        
    // Material properties
    MaterialProperty<Real> & _kappa;
    MaterialProperty<Real> & _kappa_max;
    MaterialProperty<Real> & _residual;    
    
    // Boolean for jump
    // bool _isJumpOn;
    // Boolean for first order
    //bool _is_first_order;    
    
    // Multiplicative coefficient for viscosity:
    double _Ce;
    double _Cjump;
    double _Cmax;
    
    // UserObject: equation of state
    const HydrostaticPressure & _eos;
    
    // Name of the post-processor:
    //std::string _entropy_pps_name;

};

#endif // ComputeViscCoeff_H
