[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 10
  xmax = 0.304
  ymax = 0.0257
  second_order = true
[]

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
  # x-star velocity
  [./u_star]
    order = SECOND
    family = LAGRANGE
    #initial_from_file_var = u_star
    #initial_from_file_timestep = LATEST

    [./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]

  # y-star velocity
  [./v_star]
    order = SECOND
    family = LAGRANGE
    #initial_from_file_var = v_star
    #initial_from_file_timestep = LATEST

    [./InitialCondition]
      type = ConstantIC
      value = 0.05
    [../]
  [../]

  [./u]
    order = SECOND
    family = L2_LAGRANGE
    #initial_from_file_var = u
    #initial_from_file_timestep = LATEST

    [./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  [./v]
    order = SECOND
    family = L2_LAGRANGE
    #initial_from_file_var = v
    #initial_from_file_timestep = LATEST

    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]
  [./p_current]
    order = FIRST
    family = LAGRANGE
    #initial_from_file_var = p
    #initial_from_file_timestep = LATEST
  [../]

  [./p_old]
    order = FIRST
    family = LAGRANGE
    #initial_from_file_var = p
    #initial_from_file_timestep = LATEST
  [../]
[]

[Kernels]
  [./PressurePoisson]
    type = NavStokesPressurePoisson_p
    variable = p
    u_star = u_star
    v_star = v_star
    p = p
  [../]
[]

[AuxKernels]
  [corrector_x]
    type = Corrector
    variable = u
    execute_on = timestep_end
    p = p
    u_star = u_star
    v_star = v_star
    component = 0
  []
  [corrector_y]
    type = Corrector
    variable = v
    execute_on = timestep_end
    p = p
    u_star = u_star
    v_star = v_star
    component = 1
  []
[]

[BCs]
  [./pressure_outlet]
    type = DirichletBC
    variable = p
    boundary = 'right'
    value = 0
  [../]
[]

[Materials]
  [./const]
    type = GenericConstantMaterial
    #block = 'FLUID'
    #block = 'SOLID'
    prop_names = 'rho mu'
    prop_values = '1  0.002'  #Re 100 : 0.002; Re 500 : 0.0004
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    #solve_type = 'NEWTON'
    solve_type = 'LINEAR'
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 100
  dt = .1
  dtmin = .1

  petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter' #USED FOR RE 100
  petsc_options_value = 'hypre boomeramg 6'

  #petsc_options_iname = '-pc_type'
  #petsc_options_value = 'lu'

  #petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_levels'
  #petsc_options_value = 'asm      2               ilu          4'

  #line_search = 'none'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  nl_max_its = 6
  l_tol = 1e-6
  l_max_its = 500
  picard_max_its = 1
[]

[MultiApps]
  [./sub_predictor]
    type = TransientMultiApp
    input_files = NS_Predictor.i
    execute_on = TIMESTEP_BEGIN
  [../]
[]

[Transfers]
  [./u_star_from_sub_predictor]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_predictor
    source_variable = u_star
    variable = u_star
  [../]

  [./v_star_from_sub_predictor]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_predictor
    source_variable = v_star
    variable = v_star
  [../]

  [./u_to_sub_predictor]
    type = MultiAppCopyTransfer
    #type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_predictor
    source_variable = u
    variable = u
  [../]

  [./v_to_sub_predictor]
    type = MultiAppCopyTransfer
    #type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_predictor
    source_variable = v
    variable = v
  [../]

  [./p_to_sub_predictor]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_predictor
    source_variable = p
    variable = p
  [../]
[]

[Outputs]
  file_base = channel
  exodus = true
  checkpoint = true
[]
