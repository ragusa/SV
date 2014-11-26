#include "ComputeViscCoeff.h"

template<>
InputParameters validParams<ComputeViscCoeff>()
{
  InputParameters params = validParams<Material>();
    params.addParam<std::string>("viscosity_name", "FIRST_ORDER", "Name of the viscosity definition to use: set to first order by default.");
    params.addParam<bool>("isJumpOn", true, "Is jump on?.");
    params.addParam<bool>("isShock", false, "Is a low Mach shock?.");
    params.addRequiredCoupledVar("velocity_x", "x component of the velocity");
    params.addCoupledVar("velocity_y", "y component of the velocity");
    params.addCoupledVar("jump_grad_entropy", "jump of entropy gradient");
    // constant parameters:
    params.addParam<double>("Ce"   , 1.0, "Coefficient for entropy viscosity");
    params.addParam<double>("Cjump", 1.0, "Coefficient for jumps");
    params.addParam<double>("Cmax" , 0.5, "Coefficient for first-order viscosity");
    // Userobject:
    // params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    // PPS names:
    params.addParam<std::string>("press_PPS_name", "name of the pps computing pressure");
    return params;
}

ComputeViscCoeff::ComputeViscCoeff(const std::string & name, InputParameters parameters) :
    Material(name, parameters),
    // Declare viscosity types
    _visc_name(getParam<std::string>("viscosity_name")),
    _visc_type("NONE, FIRST_ORDER, ENTROPY, INVALID", _visc_name), // jcr so what about that enum?
    // Booleans
    _isJumpOn(getParam<bool>("isJumpOn")), // jcr note: purpose?
    _isShock(getParam<bool>("isShock")),
    // Declare aux variables: velocity
    _vel_x(coupledValue("velocity_x")),
    _vel_y(_mesh.dimension()>=2 ? coupledValue("velocity_y") : _zero),
    _grad_vel_x(coupledGradient("velocity_x")),
    _grad_vel_y(_mesh.dimension()>=2 ? coupledGradient("velocity_y") : _grad_zero),
    // entropy:
    _entropy(coupledValue("entropy")),
    _entropy_old(coupledValueOld("entropy")),
    _entropy_older(coupledValueOlder("entropy")),
    _grad_entropy(coupledGradient("entropy")),
    // Jump of entropy gradients:
    _jump_grad_entropy(isCoupled("jump_grad_entropy") ? coupledValue("jump_grad_entropy") : _zero),
    //jcr note: _area(coupledValue("area")),
    //_grad_area(isCoupled("area") ? coupledGradient("area") : _grad_zero),
    // Declare material properties
    _mu(declareProperty<Real>("mu")),
    _mu_max(declareProperty<Real>("mu_max")),
    _kappa(declareProperty<Real>("kappa")),
    _kappa_max(declareProperty<Real>("kappa_max")),
//    _residual(declareProperty<Real>("residual")),
    // Get parameter Ce
    _Ce(getParam<double>("Ce")),
    _Cjump(getParam<double>("Cjump")),
    _Cmax(getParam<double>("Cmax")),
    // UserObject:
    //_eos(getUserObject<EquationOfState>("eos")),
    // PPS name:
    _entropy_pps_name(getParam<std::string>("PPS_name"))
{
    if (_Ce < 0.) || (_Ce > 2.)
        mooseError("The coefficient Ce has to be positive and cannot be larger than 2.");
}

void
ComputeViscCoeff::computeQpProperties()
{
    // Determine h (length used in definition of first and second order viscosities):
    Real _h_min = _current_elem->hmin();// /_qrule->get_order();
    
    // Compute first order viscosity:
    Real c = std::sqrt(_gravity(0)*_h[_qp]));
    _mu_max[_qp]    = _Cmax*_h_min*_norm_vel[_qp];
    _kappa_max[_qp] = _Cmax*_h_min*(_norm_vel[_qp] + c);
    
    // Epsilon value normalization of unit vectors:
    Real eps = std::sqrt(std::numeric_limits<Real>::min());
    
    // Compute Mach number and velocity variable to use in the normalization parameter:
    Real entropy_pps = std::max(getPostprocessorValueByName(_entropy_pps_name), eps);
        
    // Initialize some vector, values, ... for entropy viscosity method:
    RealVectorValue vel(_vel_x[_qp], _vel_y[_qp], _vel_z[_qp]);
    
    Real norm = 0.;
    Real jump = 0.; Real residual = 0.;
    Real kappa_e = 0.; Real mu_e = 0.;
    Real weight0 = 0.; Real weight1 = 0.; Real weight2 = 0.;
//    bool isentropic = false;
    
    // Switch statement over viscosity type:
    switch (_visc_type) {
        case NONE:             
            _mu[_qp]    = 0.;
            _kappa[_qp] = 0.;
            break;
        case FIRST_ORDER:
            _mu[_qp]    = _mu_max[_qp];
            _kappa[_qp] = _kappa_max[_qp];
            break;
        case ENTROPY:
            // Compute the viscosity coefficients:
            if (_t_step == -1) {
                _mu[_qp]    = _kappa_max[_qp];
                _kappa[_qp] = _kappa_max[_qp];
            }
            else {
                // Compute the weight for BDF2
                weight0 = (2.*_dt+_dt_old)/(_dt*(_dt+_dt_old));
                weight1 = -(_dt+_dt_old)/(_dt*_dt_old);
                weight2 = _dt/(_dt_old*(_dt+_dt_old));
                
                // Compute the characteristic equation u:
                residual = 0.;
                residual = vel*_grad_press[_qp];
                residual += (weight0*_pressure[_qp]+weight1*_pressure_old[_qp]+weight2*_pressure_older[_qp]);
                residual -= c*c*vel*_grad_rho[_qp];
                residual -= c*c*(weight0*_rho[_qp]+weight1*_rho_old[_qp]+weight2*_rho_older[_qp]);
                residual *= _Ce;
                
                // Compute kappa_e:
                if (_isJumpOn)
                    jump = _Cjump*_norm_vel[_qp]*std::max( _jump_grad_entropy[_qp], c*c*_jump_grad_dens[_qp] );
                else
                    jump = _Cjump*_norm_vel[_qp]*std::max( _grad_press[_qp].size(), c*c*_grad_rho[_qp].size() );
//                norm = 0.5*(std::fabs(1.-Mach)*_rho[_qp]*c*c + Mach*_rho[_qp]*std::min(_norm_vel[_qp]*_norm_vel[_qp], c*c));
                norm = 0.5 * _rho[_qp] * c * c;
                kappa_e = _h_min*_h_min*(std::fabs(residual) + jump) / norm;
                kappa_e += _h_min*_h_min*std::fabs(vel*_grad_area[_qp])/_area[_qp];

                // Compute mu_e:
                if (_isJumpOn)
                    jump = _Cjump*_norm_vel[_qp]*std::max( _jump_grad_press[_qp], c*c*_jump_grad_dens[_qp] );
                else
                    jump = _Cjump*_norm_vel[_qp]*std::max( _grad_press[_qp].size(), c*c*_grad_rho[_qp].size() );
                
                if (_isShock)
                    norm = 0.5 * std::max(_rho[_qp]*std::min(_norm_vel[_qp]*_norm_vel[_qp], c*c), (1.-Mach)*rhov2_pps );
                
                mu_e = _h_min*_h_min*(std::fabs(residual) + jump) / norm;
                mu_e += _h_min*_h_min*std::fabs(vel*_grad_area[_qp])/_area[_qp];

                // Compute mu and kappa:
                _mu[_qp]    = std::min( _kappa_max[_qp], mu_e);
                _kappa[_qp] = std::min( _kappa_max[_qp], kappa_e);
            }
            break;
        default:
            mooseError("The viscosity type entered in the input file is not implemented.");
            break;
    }
}
