#include "SolidWallBC.h"
#include "EquationOfState.h"

/** Wall boundary condition: the normal velocity is assumed to be zero. **/

template<>
InputParameters validParams<SolidWallBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  // Equation name:
  params.addParam<std::string>("equ_name", "invalid", "Name of the equation.");
  // Coupled variables
  params.addRequiredCoupledVar("h", "water height");
  params.addRequiredCoupledVar("hu", "x-mom of h*vec{u}");
  params.addCoupledVar("hv", "y-mom of h*vec{u}");   
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");

  return params;
}


SolidWallBC::SolidWallBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Equation name
    _equ_type("continuity x_mom y_mom invalid", getParam<std::string>("equ_name")),
    // Coupled variables
    _h(coupledValue("h")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos")),
    // Integer for jacobian terms
    _h_var(coupled("h")),
    _hu_var(coupled("hu")),
    _hv_var(_mesh.dimension() == 2 ? coupled("hv") : 0)
{}

Real
SolidWallBC::computeQpResidual()
{
  RealVectorValue hU(0., 0., 0.);
  Real p = _eos.pressure(_h[_qp], hU);
  switch (_equ_type)
  {
    case continuity:
      return 0.;
      break;
    case x_mom:
      return p*_normals[_qp](0)*_test[_i][_qp];
    case y_mom:
      return p*_normals[_qp](1)*_test[_i][_qp];
      break;
    default:
      mooseError("'" << this->name() << "' Invalid equation name.");
  }
}

Real
SolidWallBC::computeQpJacobian()
{
  RealVectorValue hU(0., 0., 0.);
  switch (_equ_type)
  {
    case continuity:
      return 0.;
      break;
    case x_mom:
      return _eos.dp_dhu(_h[_qp], hU)*_normals[_qp](0)*_test[_i][_qp];
      break;
    case y_mom:
      return _eos.dp_dhv(_h[_qp], hU)*_normals[_qp](1)*_test[_i][_qp];
      break;
    default:
      mooseError("'" << this->name() << "' Invalid equation name.");
  }
}

Real
SolidWallBC::computeQpOffDiagJacobian(unsigned jvar)
{
  RealVectorValue hU(0., 0., 0.);  
  if (jvar == _hu_var)
  {
    switch (_equ_type)
    {
      case y_mom:
        return _eos.dp_dhu(_h[_qp], hU)*_normals[_qp](1)*_test[_i][_qp];
        break;
      default:
        return 0.;
        break;
    }
  }
  else if (jvar == _hv_var)
  {
    switch (_equ_type)
    {
      case x_mom:
        return _eos.dp_dhv(_h[_qp], hU)*_normals[_qp](0)*_test[_i][_qp];
        break;
      default:
        return 0.;
        break;
    }
  }
  else
    return 0.;
}