[Mesh]
  dim = 2
#  file=partial-dam-break-small.e
  file=partial-dam-break-line.e
  uniform_refine = 2
[]

[Functions]
  [./topology]
    type = ConstantFunction
    value = 0.
  [../]
  
  [./ic_func]
    axis = 0
    type = PiecewiseLinear
    x = '-700  -0.5   0.5  700'
    y = '10    10   9.5    9.5'
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
    order = first
    [./InitialCondition]
      type = FunctionIC
      function = ic_func
    [../]
  [../]

  [./q_x]
    family = LAGRANGE
    order = first  
    [./InitialCondition]
    type = ConstantIC
    value = 0.
    [../]
  [../]
  
  [./q_y]
    family = LAGRANGE
    order = first  
    [./InitialCondition]
    type = ConstantIC
    value = 0.
    [../]
  [../]
[]

[Kernels]
  [./Time_Continuity]
    type = TimeDerivative
    variable = h
  [../]

  [./InviscidFlx_Continuity]
    type = SV_Continuity
    variable = h
    q_x = q_x
    q_y = q_y
  [../]
  
  [./ViscousFlx_Continuity]
    type = SV_ArtificialViscFlux
    variable = h
    equation_name = CONTINUITY
  [../]

  # x-momentum equation
  [./Time_X-Momentum]
    type = TimeDerivative
    variable = q_x
  [../]

  [./InviscidFlx_X-Momentum]
    type = SV_Momentum
    variable = q_x
    h = h
    q_x = q_x
    q_y = q_y
    gravity = 9.81
    component = 0
    eos = hydro
  [../]

  [./ViscousFlx_X-Momentum]
    type = SV_ArtificialViscFlux
    variable = q_x
    equation_name = X_MOMENTUM
  [../]

  # y-momentum equation  
  [./Time_Y-Momentum]
    type = TimeDerivative
    variable = q_y
  [../]
  
  [./InviscidFlx_Y-Momentum]
    type = SV_Momentum
    variable = q_y
    h = h
    q_x = q_x
    q_y = q_y
    gravity = 9.81
    component = 1
    eos = hydro
  [../]

  [./ViscousFlx_Y-Momentum]
    type = SV_ArtificialViscFlux
    variable = q_y
    equation_name = Y_MOMENTUM
  [../]

[]

[AuxVariables]
  [./entropy_aux]
    family = LAGRANGE
    order = FIRST
  [../]

  [./F_aux]
    family = LAGRANGE
    order = FIRST
  [../]

  [./G_aux]
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

[AuxKernels]
  [./entropy_ak]
    type = EntropyAux
    variable = entropy_aux
    h = h
    q_x = q_x
    q_y = q_y
  [../]

  [./F_ak]
    type = EntropyFluxAux
    variable = F_aux
    current_momentum = q_x
    h = h
    q_x = q_x
    q_y = q_y
  [../]

  [./G_ak]
    type = EntropyFluxAux
    variable = G_aux
    current_momentum = q_y
    h = h
    q_x = q_x
    q_y = q_y
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

[Materials]
  [./ViscosityCoeff]
    type = ComputeViscCoeff
    block = 1
    h = h
    q_x = q_x
    q_y = q_y
    entropy = entropy_aux
    F = F_aux
    G = G_aux    
    eos = hydro
  [../]
[]

[BCs]
  [./bc_h]
    type = SolidWallBC
    variable = h
    boundary = '1 2 3'
    equ_name = CONTINUITY
    h = h
    q_x = q_x
    q_y = q_y
    eos = hydro
  [../]

  [./bc_q_x]
    type = SolidWallBC
    variable = q_x
    boundary = '1 2 3'
    equ_name = X_MOMENTUM
    h = h
    q_x = q_x
    q_y = q_y
    eos = hydro
  [../]
  
  [./bc_q_y]
    type = SolidWallBC
    variable = q_y
    boundary = '1 2 3'
    equ_name = Y_MOMENTUM
    h = h
    q_x = q_x
    q_y = q_y
    eos = hydro
  [../]
[]

[Preconditioning]
  [./FDP]
    type = SMP # FDP
    full = true
    solve_type = 'PJFNK'
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  
  dt = 1.e-2
  
  [./TimeStepper]
    type = FunctionDT
    time_t =  '0  0.1.'
    time_dt = '5.e-3  5.e-2'
  [../]

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-8
  nl_max_its = 30

  end_time = 10.
#  num_steps = 10

#  [./Adaptivity]
#    refine_fraction = 0.5
#    coarsen_fraction = 0.05
#    max_h_level = 4
#    error_estimator = PatchRecoveryErrorEstimator
#  [../]

  [./Quadrature]
    type = GAUSS
    order = SECOND
  [../]
[]

[Outputs]
  output_initial = true
  exodus = true
  print_linear_residuals = false
  print_perf_log = true
[]
