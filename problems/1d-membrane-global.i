########################
### global parameters
########################
# for simplicity, gravity is put as a global parameter!
# Caveat: do not forget to change the status of implicit for your time discretization
#   note: we wish to have the boundary conditions solved implicitly, 
#         regardless of the time integration
[GlobalParams]
  gravity = 9.81
  implicit=true
[]

########################
### mesh 
########################
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
  xmin = 0.
  xmax = 2000.
[]

########################
### mesh 
########################
[Functions]
  [./ic_func_height]
    axis = 0
    type = PiecewiseLinear
    x = '0  1000  1000.01 2000'
    y = '10 10    0.5     0.5'
  [../]
[]

########################
### user objects
########################
# equation of state
[UserObjects]
  [./hydro]
    type = HydrostaticPressure
  [../]
[]

########################
### variables
########################
# h: water height
# q_d: d-component of momentum
# the initial conditions is typically what needs to be changed for different testcases
[Variables]
  [./h]
    family = LAGRANGE
    order = FIRST
    [./InitialCondition]
      type = FunctionIC
      function = ic_func_height
    [../]
  [../]

  [./q_x]
    family = LAGRANGE
    order = FIRST  
    [./InitialCondition]
      type = ConstantIC
      value = 0.
    [../]
  [../]
[]

########################
### FEM kernels 
########################
[Kernels]
  [./Continuity_Time]
    type = TimeDerivative
    variable = h
  [../]

  [./Continuity_InviscidFlx]
    type = SV_Continuity
    variable = h
    q_x = q_x
  [../]

  [./Continuity_ViscousFlx]
    type = SV_ArtificialViscFlux
    variable = h
    equation_name = CONTINUITY
  [../]

  [./Momentum_Time]
    type = TimeDerivative
    variable = q_x
  [../]

  [./Momentum_InviscidFlx]
    type = SV_Momentum
    variable = q_x
    h = h
    q_x = q_x
    component = 0
    eos = hydro
  [../]
  
  [./Momentum_ViscousFlx]
    type = SV_ArtificialViscFlux
    variable = q_x
    equation_name = X_MOMENTUM  # why equation_name and component? better to only have 1
  [../]
[]

########################
### auxiliary variables, as needed
########################
### velocity: for output
### entropy and entropy flux for the EVM
### kappa: viscosity coefficient
[AuxVariables]
  [./vel_aux]
    family = LAGRANGE
    order = FIRST
  [../]

  [./entropy_aux]
    family = LAGRANGE
    order = FIRST
  [../]

  [./F_aux]
    family = LAGRANGE
    order = FIRST
  [../]

  [./kappa_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]
  
  [./kappa_max_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

########################
### auxiliary kernels, as needed
########################
[AuxKernels]
  [./vel_ak]
    type = VelocityAux
    variable = vel_aux
    h = h
    q_x = q_x
  [../]

 [./entropy_ak]
    type = EntropyAux
    variable = entropy_aux
    h = h
    q_x = q_x
  [../]

  [./F_ak]
    type = EntropyFluxAux
    variable = F_aux
    momentum_component = q_x
    h = h
    q_x = q_x
  [../]

  [./kappa_ak]
    type = MaterialRealAux
    variable = kappa_aux
    property = kappa
  [../]

  [./kappa_max_ak]
    type = MaterialRealAux
    variable = kappa_max_aux
    property = kappa_max
  [../]
[]

########################
### materials
########################
[Materials]
  [./ViscosityCoeff]
    type = ComputeViscCoeff
    block = 0
    h = h
    q_x = q_x
    entropy = entropy_aux
    F = F_aux
    eos = hydro
  [../]
[]

########################
### boundary conditions
########################
[BCs]
  [./left_h]
    type = DirichletBC
    variable = h
    boundary = left
    value = 10.
  [../]

  [./right_h]
    type = DirichletBC
    variable = h
    boundary = right
    value = 0.5
  [../]

  [./left_q_x]
    type = DirichletBC
    variable = q_x
    boundary = left
    value = 0.
  [../]

  [./right_q_x]
    type = DirichletBC
    variable = q_x
    boundary = right
    value = 0.
  [../]
[]

########################
### preconditioner
########################
[Preconditioning]
  [./FDP]
    type = FDP
    full = true
    solve_type = 'PJFNK'
  [../]
[]

########################
### run options
########################
[Executioner]
  type = Transient
  scheme = bdf2

  dt = 1.e-0

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-6
  nl_max_its = 10

#  end_time = 50
  num_steps = 1000
  [./Quadrature]
    type = TRAP
  [../]

[]

########################
### output
########################
[Outputs]
  file_base = memb_test
  output_initial = true
  exodus = true
  print_linear_residuals = false
  print_perf_log = true
[]

########################
### debugging
########################
#[Debug]
#  show_var_residual = 'h q_x'
#  show_var_residual_norms = true
#[]