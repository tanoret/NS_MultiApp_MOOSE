[Mesh]
  second_order = true
  [fmg]
    type = FileMeshGenerator
    file = exodus.exo
  []
[]
#[Mesh]
#  file = NACA_airfoil_Pred.e
#[]

[Variables]
  [./u]
    order = SECOND
    family = LAGRANGE_VEC
  [../]
[]

[AuxVariables]
  [./p]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[ICs]
  [./velocity]
    type = VectorConstantIC
    variable = u
    x_value = 1
    y_value = 0
  [../]
[]

[Kernels]
  [./predictor_timederivative]
    type = INSADMomentumTimeDerivative
    variable = u
  [../]

  [./predictor_advection]
    type = INSADMomentumAdvection
    variable = u
  [../]

  [./predictor_viscous]
    type = INSADMomentumViscous
    variable = u
  [../]
[]

[BCs]
  [./no_slip]
    type = VectorFunctionDirichletBC
    variable = u
    boundary = 'WALL AIRFOIL'
  [../]

  [./velocity_inlet]
    type = VectorFunctionDirichletBC
    variable = u
    boundary = 'INLET'
    function_x = 1.0
    function_y = 0.0
  []
[]

[Materials]
  [./const]
    type = ADGenericConstantMaterial
    block = 'FLUID'
    prop_names = 'rho mu'
    prop_values = '1  0.0004'
  [../]
  [ins_mat]
    type = INSADMaterial
    velocity = u
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
  num_steps = 50
  dt = .1
  dtmin = .1
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  #petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_levels'
  #petsc_options_value = 'asm      2               ilu          4'
  #line_search = 'none'
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-8
  nl_max_its = 40
  l_tol = 1e-6
  l_max_its = 500
[]

[Outputs]
  file_base = NACA_airfoil_Pred
  exodus = true
  checkpoint = true
[]
