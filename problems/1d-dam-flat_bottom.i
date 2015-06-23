### This is example #3 from "Adaptive artificial viscosity for the Saint-Venant equations",
###     by Chen, Kurganov, Lei, Liu. 

########################
### global parameters
########################
# for simplicity, gravity is put as a global parameter!
# Caveat: do not forget to change the status of implicit for your time discretization
#   note: we wish to have the boundary conditions solved implicitly, 
#         regardless of the time integration
[GlobalParams]
  gravity = 1.0
  implicit=true
[]

########################
### mesh 
########################
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 400
  xmin = -5.
  xmax =  5.
[]

########################
### mesh 
########################
[Functions]
  [./ic_func_height]
    type = StepFunction
    x0                 = 0.
    value_before_step  = 3.
    value_after_step   = 1.
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
  [./jump]
    type = JumpInterface
    entropy_flux_x = F_aux
    var_name_jump = jump_aux
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

  [./residual_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_aux]
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

  [./residual_ak]
    type = MaterialRealAux
    variable = residual_aux
    property = residual
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
    viscosity_name = ENTROPY
    Ce = 2
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
    value = 3.
  [../]

  [./right_h]
    type = DirichletBC
    variable = h
    boundary = right
    value = 1.0
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
    petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
    petsc_options_value = '1.e-10       ds             ds'
    line_search = 'default'
  [../]
[]

########################
### run options
########################
[Executioner]
  type = Transient
  scheme = bdf2

  dt = 5.e-3

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-6
  nl_max_its = 10

  end_time = 2
#  num_steps = 2
  [./Quadrature]
    type = TRAP
  [../]

[]

########################
### output
########################
[Outputs]
  file_base = test1
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