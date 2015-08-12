#include "SVSetWaterVelocityInletBC.h"
#include "HydrostaticPressure.h"

/** Set the fluid velocity at the boundary **/

template<>
InputParameters validParams<SVSetWaterVelocityInletBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  // Equation name:
  params.addParam<std::string>("equ_name", "invalid", "Name of the equation.");
  // Coupled variables
  params.addRequiredCoupledVar("h", "water height");
  params.addRequiredCoupledVar("q_x", "x-mom of h*vec{u}");  
  // Constants and parameters
  params.addRequiredParam<Real>("u_bc", "boundary value of the value");
  params.addParam<Real>("h_bc", "boundary value of the water height (only used if fluid becomes supersonic at the inlet.");  
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");

  return params;
}


SVSetWaterVelocityInletBC::SVSetWaterVelocityInletBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    // Equation name
    _equ_type("CONTINUITY X_MOMENTUM INVALID", getParam<std::string>("equ_name")),
    // Coupled variables
    _h(coupledValue("h")),
    // Constants and parameters
    _u_bc(getParam<Real>("u_bc")),
    _h_bc(isCoupled("h_bc") ? getParam<Real>("h_bc") : 0.),
    // Equation of state:
    _eos(getUserObject<HydrostaticPressure>("eos")),
    // Integer for jacobian terms
    _h_var(coupled("h")),
    _q_x_var(coupled("q_x"))
{
  // Assert mesh dimension
  if (_mesh.dimension() > 1)
    mooseError("'" << this->name() << "' can only be used with 1-D mesh since it is designed for the Saint-Venant equations.");
  // Determine whether or not u_bc is specified in the input file
  _is_h_bc_specified = isCoupled("h_bc") ? true : false;
}

Real
SVSetWaterVelocityInletBC::computeQpResidual()
{
  // Check that the bc is an inlet bc
  if (_u_bc*_normals[_qp](0)>0)
    mooseError("'" << this->name() << "' is not/no longer an inlet bc: 'vec{u} dot vec{normal}' is greater than zero");

  // Current bc values of the momentum, sound speed and pressure
  RealVectorValue q_bc(_h[_qp]*_u_bc, 0., 0.);
  Real c2 = _eos.c2(_h[_qp], q_bc);
  Real Mach = std::fabs(_u_bc)/std::sqrt(c2);
  Real p_bc = _eos.pressure(_h[_qp], q_bc);
  Real h_bc = _h[_qp];

  // If the fluid is supercritical u_bc is used to evaluate q_bc at the boundary
  if (Mach>1)
  {
    if (!_is_h_bc_specified)
      mooseError("'" << this->name() << "': the fluid becomes supercritical but you did not specify an inlet water height value in the input file.");
    q_bc(0) = _h_bc*_u_bc;
    p_bc = _eos.pressure(_h_bc, q_bc);
    h_bc = _h_bc;
  }
  else
  mooseError("'" << this->name() << "': the fluid is SUBcritical. You should impose h, not velocity. This is the wrong BC type for this problem");

  // Return flux
  switch (_equ_type)
  {
    case CONTINUITY:
  // \int test div(q) = -\int grad(test) q + \int_bd q dot n test
      return h_bc*_u_bc*_normals[_qp](0)*_test[_i][_qp];
      break;
    case X_MOMENTUM:
      return ( _u_bc*_u_bc*h_bc + p_bc )*_normals[_qp](0)*_test[_i][_qp];
      break;
    default:
      mooseError("'" << this->name() << "' Invalid equation name.");
  }
}

Real
SVSetWaterVelocityInletBC::computeQpJacobian()
{
  // for Mach>1, both h_bc and u_bc are specified, no jac terms
  // for Mach<1, this is the wrong BC type
  return 0.;
}

Real
SVSetWaterVelocityInletBC::computeQpOffDiagJacobian(unsigned jvar)
{
  // for Mach>1, both h_bc and u_bc are specified, no jac terms
  // for Mach<1, this is the wrong BC type
  return 0.;
}

