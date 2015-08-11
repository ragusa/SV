#include "JumpInterface.h"
#include "MooseError.h"

template<>
InputParameters validParams<JumpInterface>()
{
  InputParameters params = validParams<InternalSideUserObject>();

  // Components of the entropy flux
  params.addRequiredCoupledVar("entropy_flux_x", "The x-component of the entropy flux.");
  params.addCoupledVar(        "entropy_flux_y", "The y-component of the entropy flux.");
  // Variable storing the jump computed in this function
  params.addRequiredParam<std::string>("var_name_jump", "the name of the variable that will store the jump variable.");

  return params;
}

JumpInterface::JumpInterface(const std::string & name, InputParameters parameters) :
    InternalSideUserObject(name, parameters),
    // Auxiliary system
    _aux(_fe_problem.getAuxiliarySystem()),
    // Components of the entropy flux
    _grad_F_x(coupledGradient("entropy_flux_x")),
    _grad_F_x_neighbor(coupledNeighborGradient("entropy_flux_x")),
    _grad_F_y(_mesh.dimension() == 2 ? coupledGradient("entropy_flux_y") : _grad_zero),
    _grad_F_y_neighbor(_mesh.dimension() == 2 ? coupledNeighborGradient("entropy_flux_y") : _grad_zero),
    // Variable storing the jump computed in this function
    _var_name(getParam<std::string>("var_name_jump")),
    _value(0.)
{}

JumpInterface::~JumpInterface()
{
}

void
JumpInterface::initialize()
{
  NumericVector<Number> & sln = _aux.solution();
  _aux.system().zero_variable(sln, _aux.getVariable(_tid, _var_name).number());
}

void
JumpInterface::execute()
{
  //  Get degree of freedom (dof)
  dof_id_type dof_nb_aux = _current_elem->n_dofs(_aux.number(), _fe_problem.getVariable(_tid, _var_name).number());
  dof_id_type dof_nb = 0.;
  dof_id_type dof_nb_neighbor = 0.;

  // Compute the jump if the variable is defined at the node
  if (dof_nb_aux != 0) {
    _value = 0.;
    NumericVector<Number> & sln = _aux.solution();

    for (unsigned int qp = 0; qp < _q_point.size(); ++qp)
    {
      RealVectorValue grad_F(         _grad_F_x[qp](0)         , _grad_F_y[qp](1)         , 0.);
      RealVectorValue grad_F_neighbor(_grad_F_x_neighbor[qp](0), _grad_F_y_neighbor[qp](1), 0.);
//    std::cout << "print grad_F = " << grad_F << std::endl;
//    std::cout << "print grad_F_neighbor = " << grad_F_neighbor << std::endl;
      Real _value_temp = std::fabs(grad_F*_normals[qp] - grad_F_neighbor*_normals[qp]);
      _value = std::max(_value_temp, _value);
    }

    dof_nb          = _current_elem->dof_number(_aux.number() , _fe_problem.getVariable(_tid, _var_name).number(), 0);
    dof_nb_neighbor = _neighbor_elem->dof_number(_aux.number(), _fe_problem.getVariable(_tid, _var_name).number(), 0);

    sln.add(dof_nb,          _value*0.5);
    sln.add(dof_nb_neighbor, _value*0.5);
//    std::cout << "print value = " << _value << std::endl;
  }
}

void
JumpInterface::destroy()
{
}

void
JumpInterface::finalize()
{
  _aux.solution().close();

}

void
JumpInterface::threadJoin(const UserObject & uo)
{
}