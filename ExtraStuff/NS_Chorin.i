[Mesh]
  second_order = true
  [fmg]
    type = FileMeshGenerator
    file = exodus.exo
  []
[]

[Variables]
  # x-velocity
  [./u]
    order = SECOND
    family = LAGRANGE

    [./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]

  # y-velocity
  [./v]
    order = SECOND
    family = LAGRANGE

    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]

  # x-star velocity
  [./u_star]
    order = SECOND
    family = LAGRANGE

    [./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]

  # y-star velocity
  [./v_star]
    order = SECOND
    family = LAGRANGE

    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]

  # Pressure
  [./p]
    order = FIRST
    family = LAGRANGE

    [./InitialCondition]
      type = ConstantIC
      value = 0
    [../]
  [../]
[]

[Kernels]
  [./x_chorin_predictor]
    type = INSChorinPredictor
    variable = u_star
    u = u
    v = v
    u_star = u_star
    v_star = v_star
    component = 0
    predictor_type = 'STAR'
  [../]

  [./y_chorin_predictor]
    type = INSChorinPredictor
    variable = v_star
    u = u
    v = v
    u_star = u_star
    v_star = v_star
    component = 1
    predictor_type = 'STAR'
  [../]

  [./x_chorin_corrector]
    type = INSChorinCorrector
    variable = u
    u_star = u_star
    v_star = v_star
    p = p
    component = 0
  [../]

  [./y_chorin_corrector]
    type = INSChorinCorrector
    variable = v
    u_star = u_star
    v_star = v_star
    p = p
    component = 1
  [../]

  [./chorin_pressure_poisson]
    type = INSChorinPressurePoisson
    variable = p
    u_star = u_star
    v_star = v_star
  [../]
[]

[BCs]
  [./u_no_slip]
    type = DirichletBC
    variable = u
    boundary = 'WALL AIRFOIL'
    value = 0.0
  [../]

  [./v_no_slip]
    type = DirichletBC
    variable = v
    boundary = 'WALL AIRFOIL'
    value = 0.0
  [../]

  [./u_star_no_slip]
    type = DirichletBC
    variable = u_star
    boundary = 'WALL AIRFOIL'
    value = 0.0
  [../]

  [./v_star_no_slip]
    type = DirichletBC
    variable = v_star
    boundary = 'WALL AIRFOIL'
    value = 0.0
  [../]

  [./u_inlet]
    type = DirichletBC
    variable = u
    boundary = 'INLET'
    value = 1.0
  []

  [./v_inlet]
    type = DirichletBC
    variable = v
    boundary = 'INLET'
    value = 0.0
  []

  [./u_star_inlet]
    type = DirichletBC
    variable = u_star
    boundary = 'INLET'
    value = 1.0
  []

  [./v_star_inlet]
    type = DirichletBC
    variable = v_star
    boundary = 'INLET'
    value = 0.0
  []

  [./pressure_outlet]
    type = DirichletBC
    variable = p
    boundary = 'OUTLET'
    value = 0
  [../]
[]

[Materials]
  [./const]
    type = GenericConstantMaterial
    block = 'FLUID'
    prop_names = 'rho mu'
    prop_values = '1  0.002'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    solve_type = 'NEWTON'
  [../]
[]

[Executioner]
  type = Transient
  # scheme = bdf2
  num_steps = 4
  dt = .2
  dtmin = .2
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  #petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_levels'
  #petsc_options_value = 'asm      2               ilu          4'
  #line_search = 'none'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  nl_max_its = 6
  l_tol = 1e-6
  l_max_its = 500
[]

[Outputs]
  file_base = Chorin
  exodus = true
[]
