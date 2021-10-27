[Mesh]
  second_order = true
  [fmg]
    type = FileMeshGenerator
    file = BFSMeshes/exodus.exo   #Mesh is in mm! REMEMBER
  []
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
      value = 100.0
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

  [./w_star]
    order = SECOND
    family = LAGRANGE
    #initial_from_file_var = w_star
    #initial_from_file_timestep = LATEST

    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]

  [./u]
    order = SECOND
    family = LAGRANGE
    #initial_from_file_var = u
    #initial_from_file_timestep = LATEST

    [./InitialCondition]
      type = ConstantIC
      value = 100.0
    [../]
  [../]
  [./v]
    order = SECOND
    family = LAGRANGE
    #initial_from_file_var = v
    #initial_from_file_timestep = LATEST

    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]
  [./w]
    order = SECOND
    family = LAGRANGE
    #initial_from_file_var = w
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
    w_star = w_star
    p = p
  [../]
[]

[AuxKernels]
  [./normalization_auxkernel]
    type = NormalizationAuxOld
    variable = p_old
    source_variable = p
    normal_factor = 1.0
    execute_on = timestep_end
    # Note: 'normalization' or 'shift' are provided as CLI args
  [../]

  [./normalization_auxkernel2]
    type = NormalizationAux
    variable = p_current
    source_variable = p
    normal_factor = 1.0
    execute_on = timestep_end
    # Note: 'normalization' or 'shift' are provided as CLI args
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
    type = GenericConstantMaterial
    block = 'FLUID'
    prop_names = 'rho mu'
    prop_values = '1  14.8'
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
  num_steps = 1
  dt = .001
  dtmin = .001

  #petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  #petsc_options_value = 'hypre boomeramg 6'

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
  picard_max_its = 1
[]

#[Outputs]
#  file_base = NACA_airfoil_PP
#  exodus = true
#  checkpoint = true
#[]

[MultiApps]
  [./sub_predictor]
    type = TransientMultiApp
    input_files = NS_Predictor_BFS.i
    execute_on = TIMESTEP_BEGIN
  [../]
  [./sub_corrector]
    type = TransientMultiApp
    input_files = NS_Corrector_BFS.i
    execute_on = TIMESTEP_END
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

  [./w_star_from_sub_predictor]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_predictor
    source_variable = w_star
    variable = w_star
  [../]

  [./u_from_sub_corrector]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_corrector
    source_variable = u
    variable = u
  [../]

  [./v_from_sub_corrector]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_corrector
    source_variable = v
    variable = v
  [../]

  [./w_from_sub_corrector]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_corrector
    source_variable = w
    variable = w
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

  [./w_to_sub_predictor]
    type = MultiAppCopyTransfer
    #type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_predictor
    source_variable = w
    variable = w
  [../]

  [./p_to_sub_predictor]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_predictor
    source_variable = p
    variable = p
  [../]

  [./p_to_sub_corrector]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_corrector
    source_variable = p
    variable = p
  [../]

  [./p_old_to_sub_corrector]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_corrector
    source_variable = p_old
    variable = p_old
  [../]

  [./u_star_to_sub_corrector]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_corrector
    source_variable = u_star
    variable = u_star
  [../]

  [./v_star_to_sub_corrector]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_corrector
    source_variable = v_star
    variable = v_star
  [../]

  [./w_star_to_sub_corrector]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_corrector
    source_variable = w_star
    variable = w_star
  [../]

[]
