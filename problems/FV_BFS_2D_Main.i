[Mesh]
  file = BFS_2D_flat.e
  uniform_refine = 1
[]

[Problem]
  fv_bcs_integrity_check = true
[]

[Variables]
  [pressure_p]
    type = INSFVPressureVariable
  []
[]

[AuxVariables]
  [u_star]
    type = INSFVVelocityVariable
  []
  [v_star]
    type = INSFVVelocityVariable
  []
  [Ainv_x]
    type = MooseVariableFVReal
  []
  [Hu_x]
    type = MooseVariableFVReal
  []
  [Hhat_x]
    type = MooseVariableFVReal
  []
  [RHS_x]
    type = MooseVariableFVReal
  []
  [Ainv_y]
    type = MooseVariableFVReal
  []
  [Hu_y]
    type = MooseVariableFVReal
  []
  [Hhat_y]
    type = MooseVariableFVReal
  []
  [RHS_y]
    type = MooseVariableFVReal
  []
  [pressure_old]
    type = INSFVPressureVariable
  []
  [pressure_relaxed]
    type = INSFVPressureVariable
  []
[]

[FVKernels]
  [pressure_poisson_predictor]
    type = FVNavStokesPressurePredictor_p
    variable = pressure_p
    Ainv_x = Ainv_x
    Ainv_y = Ainv_y
    Hu_x = Hhat_x
    Hu_y = Hhat_y
    boundaries_to_force = 'Outlet'
  []
  [divergence]
    type = FVFunctorDivergence
    sign = -1
    x_functor = Hhat_x
    y_functor = Hhat_y
    variable = pressure_p
  []
[]

[AuxKernels]
  [Hhat_x]
    type = FVHhat
    variable = Hhat_x
    execute_on = timestep_begin
    pressure = pressure_old
    Ainv = Ainv_x
    Hu = Hu_x
    rhs = RHS_x
    momentum_component = 'x'
  []
  [Hhat_y]
    type = FVHhat
    variable = Hhat_y
    execute_on = timestep_begin
    pressure = pressure_p
    Ainv = Ainv_y
    Hu = Hu_y
    rhs = RHS_y
    momentum_component = 'y'
  []
  [pressure_relaxation]
    type = pressureRelaxation
    variable = pressure_relaxed
    execute_on = timestep_end
    pressure = pressure_p
    pressure_old = pressure_old
    pressure_relaxation = 0.1
  []
[]


[FVBCs]
  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = 'Outlet'
    variable = pressure_p
    function = 0
  []
  [outlet_p_zerograd]
    type = INSFVSymmetryPressureBC
    boundary = 'Inlet Wall'
    variable = pressure_p
  []
[]

# [Executioner]
#   type = Transient
#   solve_type = 'NEWTON'
#   petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
#   petsc_options_value = 'asm      200                lu           NONZERO'
#   line_search = 'none'
#   nl_rel_tol = 1e-3
#   num_steps = 1
# []

# [Preconditioning]
#   [./SMP]
#     type = SMP
#     full = true
#     #solve_type = 'NEWTON'
#     solve_type = 'LINEAR'
#   [../]
# []

[Executioner]
  type = Transient
  num_steps = 200
  dt = 0.1
  dtmin = 0.1
  solve_type = 'LINEAR'
  # petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  # petsc_options_value = 'asm      200                lu           NONZERO'
  # petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -snes_max_it -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol '
  # petsc_options_value = 'gmres asm lu 100 NONZERO 2 1E-14 1E-12'
  petsc_options_iname = '-ksp_type -pc_type -pc_sub_type -sub_pc_factor_levels'
  petsc_options_value = 'gmres asm ilu 4'
  nl_rel_tol = 1e-2
  nl_abs_tol = 1e-6
  nl_max_its = 40
  l_tol = 1e-2
  l_max_its = 100
[]

[MultiApps]
  [sub_predictor]
    #type = TransientMultiApp
    type = FullSolveMultiApp
    input_files = FV_BFS_2D_Predictor.i
    execute_on = TIMESTEP_BEGIN
    sub_cycling = false
  []
[]

[Transfers]
  [u_star_from_sub_predictor]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_predictor
    source_variable = u
    variable = u_star
  []

  [v_star_from_sub_predictor]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_predictor
    source_variable = v
    variable = v_star
  []

  [Ainv_x_from_sub_predictor]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_predictor
    source_variable = Ainv_x
    variable = Ainv_x
  []

  [Ainv_y_from_sub_predictor]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_predictor
    source_variable = Ainv_y
    variable = Ainv_y
  []

  [Hhat_x_from_sub_predictor]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_predictor
    source_variable = Hu_x
    variable = Hu_x
  []

  [Hhat_y_from_sub_predictor]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_predictor
    source_variable = Hu_y
    variable = Hu_y
  []

  [RHS_x_from_sub_predictor]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_predictor
    source_variable = RHS_x
    variable = RHS_x
  []

  [RHS_y_from_sub_predictor]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_predictor
    source_variable = RHS_y
    variable = RHS_y
  []

  [p_old_from_sub_predictor]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_predictor
    source_variable = pressure_mom
    variable = pressure_old
  []

  [p_to_sub_predictor]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_predictor
    source_variable = pressure_relaxed
    variable = pressure_mom
  []

  [Hhat_x_to_sub_predictor]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_predictor
    source_variable = Hhat_x
    variable = Hu_x
  []

  [Hhat_y_to_sub_predictor]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_predictor
    source_variable = Hhat_y
    variable = Hu_y
  []

[]

[Outputs]
  file_base = BFS_2D
  exodus = true
  checkpoint = true
  # perf_graph = true # prints a performance report to the terminal
[]
