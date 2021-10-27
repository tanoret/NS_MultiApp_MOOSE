[Mesh]
  second_order = true
  [fmg]
    type = FileMeshGenerator
    file = exodus.exo
  []
[]
#[Mesh]
#  file = NACA_airfoil_Corr.e
#[]

[Variables]
  [./u]
    order = SECOND
    family = LAGRANGE_VEC
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

[AuxVariables]
  [./p]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./x_corrector]
    type = ADNavStokesCorrector
    variable = u
    u = u
    p = p
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
  nl_max_its = 6
  l_tol = 1e-6
  l_max_its = 500
[]

[Outputs]
  file_base = NACA_airfoil_Corr
  exodus = true
  checkpoint = true
[]
