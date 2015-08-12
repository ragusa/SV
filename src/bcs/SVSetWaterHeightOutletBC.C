#include "SVSetWaterHeightOutletBC.h"
#include "HydrostaticPressure.h"

/** Set the water height at the outlet boundary. This function can only be used with 1-D mesh. 
 This bc function can handle both subsonic and supersonic boundaries. **/

template<>
InputParameters validParams<SVSetWaterHeightOutletBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  // Equation name:
  params.addParam<std::string>("equ_name", "invalid", "Name of the equation.");
  // Coupled variables
  params.addRequiredCoupledVar("h", "water height");  
  params.addRequiredCoupledVar("q_x", "x component of h*vec{u}");
  // Constants and parameters
  params.addRequiredParam<Real>("h_bc", "boundary value of the density/water high h");
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");

  return params;
}

SVSetWaterHeightOutletBC::SVSetWaterHeightOutletBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    // Equation name
    _equ_type("CONTINUITY X_MOMENTUM INVALID", getParam<std::string>("equ_name")),
    // Coupled variables
    _h(coupledValue("h")),
    _q_x(coupledValue("q_x")),
    // Constants and parameters
    _h_bc(getParam<Real>("h_bc")),
    // Equation of state:
    _eos(getUserObject<HydrostaticPressure>("eos")),
    // Integer for jacobian terms
    _h_var(coupled("h")),
    _q_x_var(coupled("q_x"))
{
  // Assert mesh dimension
  if (_mesh.dimension() >= 2)
    mooseError("'" << this->name() << "' can only be used with 1-D mesh since it is designed for the Saint-Venant equations.");
}

Real
SVSetWaterHeightOutletBC::computeQpResidual()
{
  // Check that the bc is an outlet bc
  Real vel = _q_x[_qp]/_h[_qp];
  if (vel*_normals[_qp](0)<0)
    mooseError("'" << this->name() << "' is not/no longer an outlet bc: 'vec{u} dot vec{normal}' is negative");

  // Current bc values of the momentum, sound speed and pressure
  RealVectorValue q_bc(_q_x[_qp], 0., 0.);
  Real c2 = _eos.c2(_h[_qp], q_bc);
  Real Mach = std::fabs(vel)/std::sqrt(c2);

  // Compute bc pressure and water height
  Real p_bc = _eos.pressure(_h_bc, q_bc);
  Real h_bc = _h_bc;
  // use computed value, not imposed value, if Mach>1
  // because we need to impose nothing if Mach>1  (both characteristics leave)
  if (Mach>1)
  {
    p_bc = _eos.pressure(_h[_qp], q_bc);
    h_bc = _h[_qp];
  }

  // Return the value of the bc flux for each equation
  switch (_equ_type)
  {
    case CONTINUITY:
      return q_bc*_normals[_qp]*_test[_i][_qp];
      break;
    case X_MOMENTUM:
      return (_u[_qp]*_u[_qp]/h_bc+p_bc)*_normals[_qp](0)*_test[_i][_qp];
      break;
    default:
      mooseError("'" << this->name() << "' Invalid equation name.");
  }
}

Real
SVSetWaterHeightOutletBC::computeQpJacobian()
{
  // Current bc values of the momentum, sound speed and pressure
  RealVectorValue q_bc(_q_x[_qp], 0., 0.);
  Real vel = _q_x[_qp]/_h[_qp];  
  Real c2 = _eos.c2(_h[_qp], q_bc);
  Real Mach = std::fabs(vel)/std::sqrt(c2);
  Real h_bc = _h_bc;
  if (Mach>1)
    h_bc = _h[_qp];

  // Return the value of the bc flux for each equation
  switch (_equ_type)
  {
  case X_MOMENTUM:
      // recall that _u is the moose var, here _u is q_x
      return _phi[_j][_qp]*(2.*_u[_qp]/h_bc+_eos.dp_dqx(h_bc, q_bc))*_normals[_qp](0)*_test[_i][_qp];
      break;
    default:
      return 0.;
  }
}

Real
SVSetWaterHeightOutletBC::computeQpOffDiagJacobian(unsigned jvar)
{
  // Current bc values of the momentum, sound speed and pressure
  RealVectorValue q_bc(_q_x[_qp], 0., 0.);
  Real vel = _q_x[_qp]/_h[_qp];
  Real c2 = _eos.c2(_h[_qp], q_bc);
  Real Mach = std::fabs(vel)/std::sqrt(c2);
  Real h_bc = _h_bc;
  if (Mach>1)
    h_bc = _h[_qp];

  if (jvar == _h_var && Mach>1) //jcr what about Mach<1?
  {
    switch (_equ_type)
    {
      case X_MOMENTUM:
      // set Mach>1, h is the code var and is not set by the BC
    // we need to take the derivative of x-mom wrt to h. below, _u )moose var) is q_x
        return _phi[_j][_qp]*(-_u[_qp]*_u[_qp]/(_h[_qp]*_h[_qp])+_eos.dp_dh(h_bc, q_bc))*_normals[_qp](0)*_test[_i][_qp];
        break;
      default:
        return 0.;
        break;
    }
  }
  else if (jvar == _q_x_var)
  {
    switch (_equ_type)
    {
      case CONTINUITY:
        return _phi[_j][_qp]*_normals[_qp](0)*_test[_i][_qp];
        break;
      default:
        return 0.;
        break;
    }
  }
  else
    return 0.;
}

