[Mesh]
  second_order = true
  [fmg]
    type = FileMeshGenerator
    file = exodus.exo
  []
[]
#[Mesh]
#  file = NACA_airfoil_PP.e
#[]

[Variables]
  # Pressure
  [./p]
    order = FIRST
    family = LAGRANGE
    #initial_from_file_var = p
    #initial_from_file_timestep = LATEST

    [./InitialCondition]
      type = ConstantIC
      value = 0
    [../]
  [../]
[]

[AuxVariables]
  [./u]
    order = SECOND
    family = LAGRANGE_VEC
    #initial_from_file_var = u
    #initial_from_file_timestep = LATEST
  [../]
[../]

[Kernels]
  [./PressurePoisson]
    type = ADNavStokesPressurePoisson
    variable = p
    u = u
  [../]
[]

[BCs]
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
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  nl_max_its = 6
  l_tol = 1e-6
  l_max_its = 500
[]

[Outputs]
  file_base = NACA_airfoil_PP
  exodus = true
  checkpoint = true
[]

[MultiApps]
  [./sub_predictor]
    type = TransientMultiApp
    input_files = NS_Predictor_sub_vector.i
    execute_on = TIMESTEP_BEGIN
  [../]
  [./sub_corrector]
    type = TransientMultiApp
    input_files = NS_Corrector_sub_vector.i
    execute_on = TIMESTEP_END
  [../]
[]

[Transfers]
  [./u_from_sub_predictor]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_predictor
    source_variable = u
    variable = u
  [../]

  [./u_from_sub_corrector]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_corrector
    source_variable = u
    variable = u
  [../]

  [./u_to_sub_predictor]
    type = MultiAppCopyTransfer_old
    direction = to_multiapp
    multi_app = sub_predictor
    source_variable = u
    variable = u
  [../]

  #[./p_to_sub_predictor]
  #  type = MultiAppCopyTransfer
  #  direction = to_multiapp
  #  multi_app = sub_predictor
  #  source_variable = p
  #  variable = p
  #[../]

  [./p_to_sub_corrector]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_corrector
    source_variable = p
    variable = p
  [../]

  [./u_to_sub_corrector]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_corrector
    source_variable = u
    variable = u
  [../]

[]
