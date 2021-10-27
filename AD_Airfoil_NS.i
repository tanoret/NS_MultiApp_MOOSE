[Mesh]
  second_order = true
  [fmg]
    type = FileMeshGenerator
    file = Mesh4.exo
  []
[]

[AuxVariables]
  [vel_x]
    order = SECOND
  []
  [vel_y]
    order = SECOND
  []
[]

[AuxKernels]
  [vel_x]
    type = VectorVariableComponentAux
    variable = vel_x
    vector_variable = velocity
    component = 'x'
  []
  [vel_y]
    type = VectorVariableComponentAux
    variable = vel_y
    vector_variable = velocity
    component = 'y'
  []
[]

[Variables]
  [./velocity]
    order = SECOND
    family = LAGRANGE_VEC
  [../]

  [./p]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[ICs]
  [./velocity]
    type = VectorConstantIC
    variable = velocity
    x_value = 1
    y_value = 0
  [../]
  [./p]
    type = ConstantIC
    variable = p
    value = 0
  [../]
[]

[Kernels]
  [./mass]
    type = INSADMass   #INSMass
    variable = p
  [../]

  [./momentum_time]
    type = INSADMomentumTimeDerivative #INSMomentumTimeDerivative
    variable = velocity
  [../]

  [./momentum_convection]
    type = INSADMomentumAdvection
    variable = velocity
  [../]

  [./momentum_viscous]
    type = INSADMomentumViscous #INSMomentumLaplaceForm
    variable = velocity
    viscous_form = 'traction'
  [../]

  [./momentum_pressure]
    type = INSADMomentumPressure
    variable = velocity
    p = p
    integrate_p_by_parts = false
  [../]
[]

[BCs]
  [./no_slip]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = 'WALL AIRFOIL'
  [../]

  [./velocity_inlet]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = 'INLET'
    function_x = 1.0
    function_y = 0.0
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
    type = ADGenericConstantMaterial
    block = 'SOLID'
    #block = 'FLUID'
    prop_names = 'rho mu'
    prop_values = '1  0.002'
  [../]
  [ins_mat]
    type = INSADMaterial
    velocity = velocity
    pressure = p
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
  num_steps = 50
  dt = .1
  dtmin = .1

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  #petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -snes_max_it -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol '
  #petsc_options_value = 'gmres asm lu 100 NONZERO 2 1E-14 1E-12'

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
