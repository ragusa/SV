[GlobalParams]
  gravity = 1.
[]

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  xmin = 0.
  xmax = 1.
[]

[Functions]
  [./ic_func_height]
    axis = 0
    type = MyFunctionIC_f1_minus_f2
    function1 = total_height
    function2 = bathymetry
  [../]

  [./bathymetry]
    type = Bathymetry1D_partial_cosine
    xmin  = 0.4
    xmax  = 0.6
  [../]

    [./total_height]
    type = ConstantFunction
    xmin  = 0.1
    xmax  = 0.2
    value = 1.
  [../]
[]

[UserObjects]
  [./hydro]
    type = HydrostaticPressure
  [../]
[]

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
      value = 0.3
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
    type = SV_ArtificialViscFlux
    variable = h
    equation_name = CONTINUITY
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
    B = bathy_aux
    component = 0
    eos = hydro
  [../]
  
  [./ViscousFlxMomentum]
    type = SV_ArtificialViscFlux
    variable = q_x
    equation_name = X_MOMENTUM
  [../]
[]

[AuxVariables]
  [./bathy_aux]
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

[AuxKernels]
  [./bathy_ak]
    type = FunctionAux
    variable = bathy_aux
    function = bathymetry
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
    current_momentum = q_x
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

[BCs]
  [./left_h]
    type = DirichletBC
    variable = h
    boundary = left
    value = 20.
  [../]

  [./right_h]
    type = DirichletBC
    variable = h
    boundary = right
    value = 15.
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

  dt = 1.e-2

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-6
  nl_max_its = 10

#  end_time = 50
  num_steps = 1000
  [./Quadrature]
    type = TRAP
  [../]

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