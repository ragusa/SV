[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
[]

[Functions]
  [./ic_func]
    axis = 0
    type = PiecewiseLinear
    x = '0  0.5  0.5001  1'
    y = '1  1    0.5     0.5'
  [../]
[]

[UserObjects]
  [./hydro]
    type = HydrostaticPressure
    gravity = 9.81
  [../]
[]

[Variables]
  [./h]
    family = LAGRANGE
    order = FIRST
    [./InitialCondition]
      type = FunctionIC
      function = ic_func
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

[Kernels]
  [./TimeContinuity]
    type = TimeDerivative
    variable = h
  [../]

  [./InviscidFlxContinuity]
    type = SV_Continuity
    variable = h
    q_x = q_x
  [../]

  [./ViscousFlxContinuity]
    type = SV_ArtificalViscFlux
    variable = h
    equ_name = CONTINUITY
  [../]

  [./TimeMomentum]
    type = TimeDerivative
    variable = q_x
  [../]

  [./InviscidFlxMomentum]
    type = SV_Momentum
    variable = q_x
    h = h
    q_x = q_x
    gravity = 9.81
    component = 0
    eos = hydro
  [../]
  
  [./ViscousFlxMomentum]
    type = SV_ArtificalViscFlux
    variable = q_x
    equ_name = XMOMENTUM
  [../]
[]

[AuxVariables]
  [./vel_aux]
    family = LAGRANGE
    order = FIRST
  [../]

  [./kappa_max_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[AuxKernels]
  [./vel_ak]
    type = VelocityAux
    variable = vel_aux
    h = h
    q_x = q_x
  [../]

#  [./kappa_max_ak]
#    type = MaterialRealAux
#    variable = kappa_max_aux
#    property = kappa_max
#  [../]
[]

[Materials]
  [./ViscosityCoeff]
    type = ComputeViscCoeff
    block = 0
    h = h
    q_x = q_x
    eos = hydro
  [../]
[]

[BCs]
  [./left_h]
    type = DirichletBC
    variable = h
    boundary = left
    value = 1.
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

[Preconditioning]
  [./FDP]
    type = FDP
    full = true
    solve_type = 'PJFNK'
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2

  dt = 1.e-4

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-6
  nl_max_its = 10

  end_time = 0.1
#  num_steps = 10

[]

[Outputs]
  output_initial = true
  exodus = true
  print_linear_residuals = false
  print_perf_log = true
[]

#[Debug]
#  show_var_residual = 'h q_x'
#  show_var_residual_norms = true
#[]