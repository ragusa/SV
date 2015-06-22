#include "ComputeViscCoeff.h"

template<>
InputParameters validParams<ComputeViscCoeff>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredParam<std::string>("viscosity_name", "Name of the viscosity definition to use.");
  // Coupled variables
  params.addRequiredCoupledVar("h"  , "h: water height");
  params.addRequiredCoupledVar("q_x", "x-component of momentum");
  params.addCoupledVar("q_y", "y-component of momentum");
  // Coupled aux variables
  params.addRequiredCoupledVar("entropy", "entropy function");
  params.addRequiredCoupledVar("F", "x-component of the entropy flux ");
  params.addCoupledVar("G", "y-component of the entropy flux ");
  params.addCoupledVar("B", "bathymetry data");  
  params.addRequiredParam<Real>("gravity", "gravity magnitude");
  // constant parameters:
  //params.addParam<bool>("is_first_order", false, "if true, use the first-order viscosity coefficient");
  params.addParam<double>("Ce"   , 1.0, "Coefficient for entropy viscosity");
  params.addParam<double>("Cjump", 1.0, "Coefficient for jumps");
  params.addParam<double>("Cmax" , 0.5, "Coefficient for first-order viscosity");
  // Userobject:
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  // PPS names:
  //params.addParam<std::string>("press_PPS_name", "name of the pps computing pressure");
    
  return params;
}

ComputeViscCoeff::ComputeViscCoeff(const std::string & name, 
                                   InputParameters parameters) :
  Material(name, parameters),
  // Declare viscosity types
  _visc_type("NONE FIRST_ORDER ENTROPY INVALID", getParam<std::string>("viscosity_name")),   
  // Declare variables
  _h(coupledValue("h")),
  _q_x(coupledValue("q_x")),
  _q_y(_mesh.dimension() == 2 ? coupledValue("q_y") : _zero),
  // entropy:
  _E(coupledValue("entropy")),
  _E_old(coupledValueOld("entropy")),
  _E_older(coupledValueOlder("entropy")),
  // entropy flux:
  _grad_F(coupledGradient("F")),
  _grad_G(_mesh.dimension() == 2 ? coupledGradient("G") : _grad_zero),
  // bathymetry:
  _bathymetry(isCoupled("B") ? coupledValue("B") : _zero),
  _gravity(getParam<Real>("gravity")),
  // Jump of entropy gradients:
  //_jump_grad_entropy(isCoupled("jump_grad_entropy") ? coupledValue("jump_grad_entropy") : _zero),
  // Declare material properties
  _kappa(declareProperty<Real>("kappa")),
  _kappa_max(declareProperty<Real>("kappa_max")),
  _residual(declareProperty<Real>("residual")), // jcr: why declare property for residual?, for output
  // Get constant parameters
  _Ce(getParam<double>("Ce")),
  _Cjump(getParam<double>("Cjump")),
  _Cmax(getParam<double>("Cmax")),
  // UserObject:
  _eos(getUserObject<HydrostaticPressure>("eos")) //,
  // PPS name:
  //_entropy_pps_name(getParam<std::string>("PPS_name"))
{
//  if (_Ce <= 0. || _Ce > 2.)
//    mooseError("ERROR in "<<this->name()<<": the coefficient Ce has to be positive and should not be larger than 2.");
}

void
ComputeViscCoeff::computeQpProperties()
{
  // Determine h (length used in definition of first and second order viscosities):
  Real _h_min = _current_elem->hmin();// /_qrule->get_order(); jcr why called min, not cell?
    
  // vector q
  RealVectorValue _vector_q( _q_x[_qp], _q_y[_qp], 0. );

  // Compute first order viscosity:
  Real c = std::sqrt(_eos.c2(_h[_qp], _vector_q));
  _kappa_max[_qp] = _Cmax*_h_min*(_vector_q.size()/_h[_qp] + c);
    
  

  // Epsilon value normalization of unit vectors:
  Real _eps = std::sqrt(std::numeric_limits<Real>::min());
    
  // Compute Mach number and velocity variable to use in the normalization parameter:
  // Real entropy_pps = std::max(getPostprocessorValueByName(_entropy_pps_name), eps);
        
  // Initialize some vector, values, ... for entropy viscosity method:
  RealVectorValue _vector_vel = _vector_q / _h[_qp];

  Real jump=0.;
  
  // Switch statement over viscosity type:
  switch (_visc_type) {
  case NONE:             
    _kappa[_qp] = 0.;
    break;
  case FIRST_ORDER:
    _kappa[_qp] = _kappa_max[_qp];
    break;
  case ENTROPY:
    // Compute the viscosity coefficients:
    if (_t_step == 1) {
      _kappa[_qp] = _kappa_max[_qp];
    }
    else {
      // Weights for BDF2
      Real w0 = _t_step > 2 ? (2.*_dt+_dt_old)/(_dt*(_dt+_dt_old)) :  1. / _dt;
      Real w1 = _t_step > 2 ? -(_dt+_dt_old)/(_dt*_dt_old)         : -1. / _dt;
      Real w2 = _t_step > 2 ? _dt/(_dt_old*(_dt+_dt_old))          :  0.      ;

      // Entropy residual
      Real residual = w0*_E[_qp] + w1*_E_old[_qp] + w2*_E_older[_qp];
      residual += _grad_F[_qp](0)+_grad_G[_qp](1);
      // store at qp
      _residual[_qp] = std::fabs(residual);
    
      // Compute kappa_e:
      /*if (_isJumpOn)
          jump = _Cjump*_norm_vel[_qp]*std::max( _jump_grad_entropy[_qp], c*c*_jump_grad_dens[_qp] );
      else
          jump = _Cjump*_norm_vel[_qp]*std::max( _grad_press[_qp].size(), c*c*_grad_rho[_qp].size() );*/

      // Froude number (use from Marco while I figure out |s-save|)
      Real Froude = _vector_q.size()/_h[_qp]/std::sqrt(_gravity*(_h[_qp]+_eps));
      Real _norm = _gravity*(std::fabs(_h[_qp])+_bathymetry[_qp]+_eps);
      Real kappa_e = _Ce*_h_min*_h_min*(std::fabs(residual) + jump) / _norm;
      mooseAssert(kappa_e<0.,"Entropy viscosity <0 in "<<this->name());

      //jump = _Cjump*_norm_vel[_qp]*std::max( _grad_press[_qp].size(), c*c*_grad_rho[_qp].size() );
                   
      // Compute kappa:
      _kappa[_qp] = std::min( _kappa_max[_qp], kappa_e);
    }
    break;
  default:
    mooseError("Error in "<<this->name()<<" The viscosity type entered in the input file is not implemented.");
    break;
  }
}
