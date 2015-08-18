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
  implicit=false
[]

########################
### mesh 
########################
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx =  1000
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
  [./ic_func]
    axis = 0
    type = PiecewiseLinear
    x = '-5.  0.   0.00001 5.'
    y = ' 3.  3.   1.      1.'
  [../]
#   [./ic_func]
#     axis = 0
#     type = PiecewiseLinear
#     x = '-5.  5.'
#     y = ' 3.  1.'
#   [../]`
[]

########################
### user objects
########################
# equation of state
[UserObjects]
  [./hydro]
    type = HydrostaticPressure
  [../]
#   [./jump]
#     type = JumpInterface
#     entropy_flux_x = F_aux
#     var_name_jump = jump_aux
#   [../]
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
      function = ic_func
#      function = ic_func_height
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
    implicit=true
  [../]

  [./Continuity_InviscidFlx]
    type = SV_Continuity
    variable = h
    q_x = h
  [../]

  [./Continuity_ViscousFlx]
    type = SV_ArtificialViscFlux
    variable = h
    equation_name = CONTINUITY
  [../]

  [./Momentum_Time]
    type = TimeDerivative
    variable = q_x
    implicit=true
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
 [./entropy_ak]
    type = ConstantAux
    variable = entropy_aux
    value =2.
#     h = h
#     q_x = q_x
  [../]

  [./F_ak]
    type = ConstantAux
    variable = F_aux
    value = 6.
#     momentum_component = q_x
#     h = h
#     q_x = q_x
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
#    jump_entropy_flux = jump_aux
    jump_entropy_flux = entropy_aux
    eos = hydro
#    viscosity_name = ENTROPY
#    viscosity_name = FIRST_ORDER
    viscosity_name = NONE
    Ce = 5.
    Cjump = 5.
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
    value = 1.
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
### postprocessor
########################
[Postprocessors]
  [./dt]
    type = TimeStepCFL
    h = h
    q_x = q_x
    eos = hydro
    cfl = 0.1
#     outputs = none
    execute_on = timestep_end
  [../]
[]

########################
### preconditioner
########################
[Preconditioning]
  active = 'fdp'
  [./dbg]
    type = SMP
    full = true
  [../]
    [./fdp]
    type = FDP
    full = true
    solve_type = 'PJFNK'
    petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
    petsc_options_value = '1.e-9                   ds             ds'
    line_search = 'default'
  [../]
  [./smp]
    type = SMP
    full = true
    solve_type = NEWTON               # Use "regular" Newton's method (Jacobian MUST be correct)
    petsc_options_iname = '-pc_type'  # LU is essentially a "direct solve"
    petsc_options_value = 'lu'
    line_search = 'none'              # Don't allow line search to cut back the steps
  [../]
[]

########################
### run options
########################
[Executioner]
  type = Transient
  scheme = explicit-euler
  solve_type = 'LINEAR'
  l_tol =1.e-10
  
  dt =1.e-4
  
#   [./TimeStepper]
#   type = PostprocessorDT
#   postprocessor = dt
#   dt = 1.e-6
# #    type = FunctionDT
# #    time_t = '0 50'
# #    time_dt= '1e-1 1e-1'
#   [../]

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-6
#   nl_max_its = 10
  
#  end_time = 2
 num_steps = 5000

 [./Quadrature]
    type = GAUSS
    order = SECOND
  [../]

[]

########################
### output
########################
[Outputs]
  [./out]
    type = CSV
    file_base = cfl
  [../]
  file_base = explicit
  output_initial = true
  exodus = true
  print_linear_residuals = false
  print_perf_log = true
[]

########################
### debugging
########################
[Debug]
  show_var_residual = 'h q_x'
  show_var_residual_norms = true
[]