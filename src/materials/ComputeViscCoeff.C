#include "ComputeViscCoeff.h"

template<>
InputParameters validParams<ComputeViscCoeff>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredParam<std::string>("viscosity_name", "Name of the viscosity definition to use.");
  // Coupled variables
  params.addRequiredCoupledVar("h"  , "h: water height");
  params.addRequiredCoupledVar("q_x", "x-component of momentum");
  params.addCoupledVar(        "q_y", "y-component of momentum");
  // Coupled aux variables: entropy pair
  params.addRequiredCoupledVar("entropy", "entropy function");
  params.addRequiredCoupledVar("F", "x-component of the entropy flux ");
  params.addCoupledVar(        "G", "y-component of the entropy flux ");
  // Coupled aux variables: bathymetry
  params.addCoupledVar(        "B", "bathymetry data");  
  // Coupled aux variables: jump
  params.addCoupledVar("jump_entropy_flux", "Jump of the of the entropy flux");
  // constant parameters:
  params.addParam<double>("Ce"   , 1.0, "Coefficient for entropy viscosity");
  params.addParam<double>("Cjump", 1.0, "Coefficient for jumps");
  params.addParam<double>("Cmax" , 0.5, "Coefficient for first-order viscosity");
  // Userobject:
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  params.addRequiredParam<Real>("gravity", "gravity magnitude");
  // PPS names:
  //params.addParam<std::string>("press_PPS_name", "name of the pps computing pressure");
    
  return params;
}

ComputeViscCoeff::ComputeViscCoeff(const InputParameters & parameters) :
  Material(parameters),
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
  // Jump 
  _jump(isCoupled("jump_entropy_flux") ? coupledValue("jump_entropy_flux") : _zero),
  // Get constant parameters
  _Ce(   getParam<double>("Ce")),
  _Cjump(getParam<double>("Cjump")),
  _Cmax( getParam<double>("Cmax")),
  // UserObject:
  _eos(getUserObject<HydrostaticPressure>("eos")),
  // gravity
  _gravity(getParam<Real>("gravity")),
  // Declare material properties
  _kappa(    declareProperty<Real>("kappa")),
  _kappa_max(declareProperty<Real>("kappa_max")),
  _residual( declareProperty<Real>("residual")) // jcr: why declare property for residual?, for output
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
  // Real _h_min = _current_elem->hmin();// /_qrule->get_order(); jcr why called min, not cell?
  // Determine h (characteristic length of the current cel):
  Real h_cell = std::pow(_current_elem->volume(), 1./_mesh.dimension());
      
  // vector q
  RealVectorValue _vector_q( _q_x[_qp], _q_y[_qp], 0. );

  // Compute sound speed:
  Real c = std::sqrt(_eos.c2(_h[_qp], _vector_q));
  
  // Compute first-order viscosity: Cmax ( |u| + c ) = Cmax ( |q|/h + c )
  _kappa_max[_qp] = _Cmax*h_cell*(_vector_q.size()/_h[_qp] + c);
    
  // Epsilon value normalization of unit vectors:
  Real _eps = std::sqrt(std::numeric_limits<Real>::min());
    
  // Compute Mach number and velocity variable to use in the normalization parameter:
  // Real entropy_pps = std::max(getPostprocessorValueByName(_entropy_pps_name), eps);
        
  // Initialize some vector, values, ... for entropy viscosity method:
  RealVectorValue _vector_vel = _vector_q / _h[_qp];
  
  // Switch statement over viscosity type:
  switch (_visc_type) {
  case NONE:             
    _kappa[_qp] = 3.e-2;
    break;
  case FIRST_ORDER:
    _kappa[_qp] = _kappa_max[_qp];
    break;
  case ENTROPY:
    // Compute the viscosity coefficients:
    // _t_step=1 is the first time step, see Transient.C
    if (_t_step == 1) {
      _kappa[_qp] = _kappa_max[_qp];
/*      Moose::out << " Compute_Viscosity:qp=" << _q_point[_qp](0) 
                 << ",\t FOV= " <<  _kappa_max[_qp]
                 << std::endl;
*/
    }
    else {
      // Weights for BDF2      
      /*
      Real sum = _dt + _dt_old;
      w0 =  1./_dt + 1./ sum;
      w1 = -sum / (_dt*_dt_old);
      w2 = _dt / (_dt_old * sum);
      */
      Real w0 = _t_step > 1 ? (2.*_dt+_dt_old)/(_dt*(_dt+_dt_old)) :  1. / _dt;
      Real w1 = _t_step > 1 ? -(_dt+_dt_old)/(_dt*_dt_old)         : -1. / _dt;
      Real w2 = _t_step > 1 ? _dt/(_dt_old*(_dt+_dt_old))          :  0.      ;

      // Entropy residual (in absolute value)
      Real residual = w0*_E[_qp] + w1*_E_old[_qp] + w2*_E_older[_qp];
      residual += _grad_F[_qp](0)+_grad_G[_qp](1);
      residual = std::fabs(residual);
      // store at qp
      _residual[_qp] = residual;
    
      // Froude number (use from Marco while I figure out |s-save|)
      Real Froude = _vector_vel.size() / c;

      // normalization is c^2
      Real _norm = std::fabs( _gravity * ( std::fabs(_h[_qp]) + _bathymetry[_qp]+ 1.e-06 ) ); // little effect on result
//      Real _norm = std::fabs( _gravity * ( std::fabs(_h[_qp]) + _bathymetry[_qp]+ _eps ) );

      // Entropy viscosity
      Real kappa_e = h_cell*h_cell* std::max( _Ce*_residual[_qp], _Cjump*_jump[_qp] ) / _norm;
//      Real kappa_e = h_cell*h_cell* ( _Ce*_residual[_qp] + _Cjump*_jump[_qp] ) / _norm;
      mooseAssert(kappa_e<0.,"Entropy viscosity <0 in "<<this->name());

      // Compute kappa:
      _kappa[_qp] = std::min( _kappa_max[_qp], kappa_e);

/*      Moose::out << " Compute_Viscosity:qp=" << _q_point[_qp](0) 
                 << ",\t FOV= "           <<  _kappa_max[_qp]
                 << ",\t entr_res[qp]= "  <<  _residual[_qp]
                 << ",\t jump[_qp]= "     <<  _jump[_qp]
                 << ",\t norm= "          <<  _norm
                 << ",\t kappa_e= "       << kappa_e
                 << ",\t kappa = "        << _kappa[_qp] 
                 << std::endl;
*/
    }
    break;
  default:
    mooseError("Error in "<<this->name()<<" The viscosity type entered in the input file is not implemented.");
    break;
  }
}
