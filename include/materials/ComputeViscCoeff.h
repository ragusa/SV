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
  ComputeViscCoeff(const InputParameters & parameters);

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
    
    // Coupled aux variables: entropy
    VariableValue    & _E;
    VariableValue    & _E_old;
    VariableValue    & _E_older;
    //VariableValue    & _E_dot;
    
    // Coupled aux variables: entropy flux
    VariableGradient & _grad_F;
    //VariableGradient & _grad_F_old;
    VariableGradient & _grad_G;
    //VariableGradient & _grad_G_old;

    // Coupled aux variable: bathymetry
    VariableValue & _bathymetry;

    // Coupled aux variables: jumps
    VariableValue & _jump;
    
    // Multiplicative coefficient for viscosity:
    double _Ce;
    double _Cjump;
    double _Cmax;
    
    // UserObject: equation of state
    const HydrostaticPressure & _eos;
    Real _gravity;
    
    // Material properties
    MaterialProperty<Real> & _kappa;
    MaterialProperty<Real> & _kappa_max;
    MaterialProperty<Real> & _residual;    
    
    // Name of the post-processor:
    //std::string _entropy_pps_name;

};

#endif // ComputeViscCoeff_H
