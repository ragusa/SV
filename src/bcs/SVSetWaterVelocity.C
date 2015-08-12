#include "SVSetWaterVelocity.h"
#include "HydrostaticPressure.h"

/** Set the fluid velocity at the boundary **/

template<>
InputParameters validParams<SVSetWaterVelocity>()
{
  InputParameters params = validParams<IntegratedBC>();

  // Equation name:
  params.addParam<std::string>("equ_name", "invalid", "Name of the equation.");
  // Coupled variables
  params.addRequiredCoupledVar("h", "water height");
  params.addRequiredCoupledVar("q_x", "x-mom of h*vec{u}");  
  // Constants and parameters
  params.addRequiredParam<Real>("u_bc", "boundary value of the value");
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");

  return params;
}


SVSetWaterVelocity::SVSetWaterVelocity(const InputParameters & parameters) :
    IntegratedBC(parameters),
    // Equation name
    _equ_type("CONTINUITY X_MOMENTUM INVALID", getParam<std::string>("equ_name")),
    // Coupled variables
    _h(coupledValue("h")),
    // Constants and parameters
    _u_bc(getParam<Real>("u_bc")),
    // Equation of state:
    _eos(getUserObject<HydrostaticPressure>("eos")),
    // Integer for jacobian terms
    _h_var(coupled("h")),
    _q_x_var(coupled("q_x"))
{
  if (_mesh.dimension() > 1)
    mooseError("'" << this->name() << "' can only be used with 1-D mesh since it is designed for the Saint-Venant equations.");
}

Real
SVSetWaterVelocity::computeQpResidual()
{
  RealVectorValue q_bc(_h[_qp]*_u_bc, 0., 0.);
  Real p = _eos.pressure(_h[_qp], q_bc);
  switch (_equ_type)
  {
    case CONTINUITY:
      return _h[_qp]*_u_bc*_normals[_qp](0)*_test[_i][_qp];
      break;
    case X_MOMENTUM:
      return (_u_bc*_u_bc*_h[_qp]+p)*_normals[_qp](0)*_test[_i][_qp];
      break;
    default:
      mooseError("'" << this->name() << "' Invalid equation name.");
  }
}

Real
SVSetWaterVelocity::computeQpJacobian()
{
  return 0.; // jcr really ?
}

Real
SVSetWaterVelocity::computeQpOffDiagJacobian(unsigned jvar)
{
  return 0.; // jcr really ?
}

