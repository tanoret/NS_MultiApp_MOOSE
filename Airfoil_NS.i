[Mesh]
  second_order = true
  [fmg]
    type = FileMeshGenerator
    file = Mesh3.exo
  []
[]

[Variables]
  [./vel_x]
    order = SECOND
    family = LAGRANGE
  [../]

  [./vel_y]
    order = SECOND
    family = LAGRANGE
  [../]

  [./pressure]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[ICs]
  [./x_velocity]
    type = ConstantIC
    variable = vel_x
    value = 1.0
  [../]
  [./y_velocity]
    type = ConstantIC
    variable = vel_y
    value = 0.0
  [../]
  [./pressure]
    type = ConstantIC
    variable = pressure
    value = 0
  [../]
[]

[Kernels]
  [./mass]
    type = INSMass
    variable = pressure
    u = vel_x
    v = vel_y
    pressure = pressure
  [../]

  [./x_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_x
  [../]
  [./y_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_y
  [../]

  [./x_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_x
    u = vel_x
    v = vel_y
    pressure = pressure
    component = 0
  [../]
  [./y_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_y
    u = vel_x
    v = vel_y
    pressure = pressure
    component = 1
  [../]
[]

[BCs]
  [./x_no_slip]
    type = DirichletBC
    variable = vel_x
    boundary = 'WALL AIRFOIL'
    value = 0.0
  [../]
  [./y_no_slip]
    type = DirichletBC
    variable = vel_y
    boundary = 'WALL AIRFOIL'
    value = 0.0
  [../]

  [./x_inlet]
    type = DirichletBC
    variable = vel_x
    boundary = 'INLET'
    value = 1.0
  [../]
  [./y_inlet]
    type = DirichletBC
    variable = vel_y
    boundary = 'INLET'
    value = 0.0
  [../]

  [./pressure_outlet]
    type = DirichletBC
    variable = pressure
    boundary = 'OUTLET'
    value = 0
  [../]
[]

[Materials]
  [./const]
    type = GenericConstantMaterial
    #block = 'FLUID'
    block = 'SOLID'
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
  num_steps = 100
  dt = .1
  dtmin = .1

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  #petsc_options_iname = ' -pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart -sub_pc_factor_shift_type'
  #petsc_options_value = ' asm lu 2 200 NONZERO'



  #petsc_options_iname = '-ksp_type -ksp_gmres_restart -pc_type -sub_pc_type -snes_max_it -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol '
  #petsc_options_value = 'gmres 30 asm lu 100 NONZERO 2 1E-14 1E-12'

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
  file_base = NACA_airfoil
  exodus = true
[]
