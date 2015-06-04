#include "SVSetWaterHeightInletBC.h"
#include "HydrostaticPressure.h"

/** Set the water height at the inlet boundary. This function can only be used with 1-D mesh. 
 This bc function can handle both subsonic and supersonic boundaries. **/

template<>
InputParameters validParams<SVSetWaterHeightInletBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  // Equation name:
  params.addParam<std::string>("equ_name", "invalid", "Name of the equation.");
  // Coupled variables
  params.addRequiredCoupledVar("q_x", "x component of h*vec{u}");
  // Constants and parameters
  params.addRequiredParam<Real>("h_bc", "boundary value of the density/water high h");
  params.addParam<Real>("u_bc", "boundary value of the fluid velocity (only used if fluid becomes supersonic at the inlet.");
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");

  return params;
}


SVSetWaterHeightInletBC::SVSetWaterHeightInletBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Equation name
    _equ_type("CONTINUITY X_MOMENTUM INVALID", getParam<std::string>("equ_name")),
    // Coupled variables
    _q_x(coupledValue("q_x")),
    // Constants and parameters
    _h_bc(getParam<Real>("h_bc")),
    _u_bc(isCoupled("u_bc") ? getParam<Real>("u_bc") : 0.),
    // Equation of state:
    _eos(getUserObject<HydrostaticPressure>("eos")),
    // Integer for jacobian terms
    _q_x_var(coupled("q_x"))
{
  // Assert mesh dimension
  if (_mesh.dimension() >= 2)
    mooseError("'" << this->name() << "' can only be used with 1-D mesh since it is designed for the Saint-Venant equations.");
  // Determine whether or not u_bc is specified in the input file
  _is_u_bc_specified = isCoupled("u_bc") ? true : false;
}

Real
SVSetWaterHeightInletBC::computeQpResidual()
{
  // Check that the bc is an inlet bc
  Real vel = _q_x[_qp]/_h_bc;
  if (vel*_normals[_qp](0)>0)
    mooseError("'" << this->name() << "' is not/no longer an inlet bc: 'vec{u} dot vec{normal}' is greater than zero");

  // Current bc values of the momentum, sound speed and pressure
  RealVectorValue q_bc(_q_x[_qp], 0., 0.);
  Real c2 = _eos.c2(_h_bc, q_bc);
  Real Mach = std::fabs(vel)/std::sqrt(c2);

  // If the fluid is supercritical u_bc is used to evaluate q_bc at the boundary
  if (Mach>1)
  {
    if (!_is_u_bc_specified)
      mooseError("'" << this->name() << "': the fluid becomes supersonic but you did not specify an inlet fluid velocity value in the input file.");
    q_bc(0) = _h_bc*_u_bc;
  }
  Real p = _eos.pressure(_h_bc, q_bc);

  // Return the value of the bc flux for each equation
  switch (_equ_type)
  {
    case CONTINUITY:
  // \int test div(q) = -\int grad(test) q + \int_bd q dot n test
      return q_bc*_normals[_qp]*_test[_i][_qp];
      break;
    case X_MOMENTUM:
      return ( _u[_qp]*_u[_qp]/_h_bc + p )*_normals[_qp](0)*_test[_i][_qp];
      break;
    default:
      mooseError("'" << this->name() << "' Invalid equation name.");
  }
}

Real
SVSetWaterHeightInletBC::computeQpJacobian()
{
  // recall that h_bc is imposed in this class
  
  // Precompute some values
  Real vel = _q_x[_qp]/_h_bc;
  RealVectorValue q_bc(_q_x[_qp], 0., 0.);
  Real c2 = _eos.c2(_h_bc, q_bc);
  Real Mach = std::fabs(vel)/std::sqrt(c2);
  Real dpdqx;

  // Return diagonal jacobian terms
  if (Mach<1)
  {
    switch (_equ_type)
    {
      case X_MOMENTUM:
        dpdqx = _eos.dp_dqx(_h_bc, q_bc);
        // recall that _u is the moose var, here _u is q_x
        return _phi[_j][_qp]*(2.*_u[_qp]/_h_bc+dpdqx)*_normals[_qp](0)*_test[_i][_qp];
        break;
      default:
        return 0.;
    }
  }
  else
    // if Mach>1, both h_bc and q_bc=h_bc*u_bc are imposed
    return 0.;
}

Real
SVSetWaterHeightInletBC::computeQpOffDiagJacobian(unsigned jvar)
{
  // recall that h_bc is imposed in this class

  // Precompute some values
  Real vel = _q_x[_qp]/_h_bc;
  RealVectorValue q_bc(_q_x[_qp], 0., 0.);
  Real c2 = _eos.c2(_h_bc, q_bc);
  Real Mach = std::fabs(vel)/std::sqrt(c2);

  // Non-zero jacobian term only if subsonic fluid
  if (jvar == _q_x_var && Mach < 1)
  {
    switch (_equ_type)
    {
      case CONTINUITY:
      // residual: \int_bc q n test --> off-diag jac below
        return _phi[_j][_qp]*_normals[_qp](0)*_test[_i][_qp];
        break;
      default:
    // since h_bc is imposed in this class, there is no off-diag jac for the momentum eq in 1d
        return 0.;
        break;
    }
  }
  else
    return 0.;
}

