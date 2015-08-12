#include "HydrostaticPressure.h"

template<>
InputParameters validParams<HydrostaticPressure>()
{
  InputParameters params = validParams<UserObject>();
  params.addRequiredParam<Real>("gravity","gravity magnitude");

  return params;
}

HydrostaticPressure::HydrostaticPressure(const InputParameters & parameters) :
    GeneralUserObject(parameters),
    _gravity(getParam<Real>("gravity"))
{}

Real HydrostaticPressure::pressure(Real h, RealVectorValue vector_q) const
{
  return 0.5*_gravity*h*h;
}

// The derivative of pressure wrt h
Real HydrostaticPressure::dp_dh(Real h, RealVectorValue vector_q) const
{
  return _gravity*h;
}

// The derivative of pressure wrt q_x
Real HydrostaticPressure::dp_dqx(Real h, RealVectorValue vector_q) const
{
  return 0.;
}

// The derivative of pressure wrt q_y
Real HydrostaticPressure::dp_dqy(Real h, RealVectorValue vector_q) const
{
  return 0.;
}

// Sound speed squared
Real
HydrostaticPressure::c2(Real h, RealVectorValue vector_q) const
{
  return _gravity*h;
}