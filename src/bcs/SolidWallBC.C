#include "SolidWallBC.h"
#include "HydrostaticPressure.h"

/** Wall boundary condition: the normal velocity is assumed to be zero. **/

template<>
InputParameters validParams<SolidWallBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  // Equation name:
  params.addParam<std::string>("equ_name", "invalid", "Name of the equation.");
  // Coupled variables
  params.addRequiredCoupledVar("h", "water height");
  params.addRequiredCoupledVar("q_x", "x-mom of h*vec{u}");
  params.addCoupledVar("q_y", "y-mom of h*vec{u}");   
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");

  return params;
}


SolidWallBC::SolidWallBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Equation name
    _equ_type("CONTINUITY X_MOMENTUM Y_MOMENTUM INVALID", getParam<std::string>("equ_name")),
    // Coupled variables
    _h(coupledValue("h")),
    // Equation of state:
    _eos(getUserObject<HydrostaticPressure>("eos")),
    // Integer for jacobian terms
    _h_var(coupled("h")),
    _q_x_var(coupled("q_x")),
    _q_y_var(_mesh.dimension() == 2 ? coupled("hv") : 0)
{}

Real
SolidWallBC::computeQpResidual()
{
  RealVectorValue hU(0., 0., 0.);
  Real p = _eos.pressure(_h[_qp], hU);
  switch (_equ_type)
  {
    case CONTINUITY:
      return 0.;
      break;
    case X_MOMENTUM:
      return p*_normals[_qp](0)*_test[_i][_qp];
    case Y_MOMENTUM:
      return p*_normals[_qp](1)*_test[_i][_qp];
      break;
    default:
      mooseError("'" << this->name() << "' Invalid equation name.");
  }
}

Real
SolidWallBC::computeQpJacobian()
{
  RealVectorValue q_bc(0., 0., 0.);
  switch (_equ_type)
  {
    case CONTINUITY:
      return 0.;
      break;
    case X_MOMENTUM:
      return _eos.dp_dqx(_h[_qp], q_bc)*_normals[_qp](0)*_test[_i][_qp];
      break;
    case Y_MOMENTUM:
      return _eos.dp_dqy(_h[_qp], q_bc)*_normals[_qp](1)*_test[_i][_qp];
      break;
    default:
      mooseError("'" << this->name() << "' Invalid equation name.");
  }
}

Real
SolidWallBC::computeQpOffDiagJacobian(unsigned jvar)
{
  RealVectorValue q_bc(0., 0., 0.);  
  if (jvar == _q_x_var)
  {
    switch (_equ_type)
    {
      case Y_MOMENTUM:
        return _eos.dp_dqx(_h[_qp], q_bc)*_normals[_qp](1)*_test[_i][_qp];
        break;
      default:
        return 0.;
        break;
    }
  }
  else if (jvar == _q_y_var)
  {
    switch (_equ_type)
    {
      case X_MOMENTUM:
        return _eos.dp_dqy(_h[_qp], q_bc)*_normals[_qp](0)*_test[_i][_qp];
        break;
      default:
        return 0.;
        break;
    }
  }
  else
    return 0.;
}
