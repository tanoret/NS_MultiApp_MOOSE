[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 3
    ny = 3
    xmin = 0.0
    xmax = 3.0
    ymin = 0.0
    ymax = 3.0
  []
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
  # [u_adv]
  #   type = INSFVVelocityVariable
  # []
  # [v_adv]
  #   type = INSFVVelocityVariable
  # []
[]

[FVKernels]
  [pressure_poisson_predictor]
    type = FVNavStokesPressurePredictor_p
    variable = pressure_p
    Ainv_x = Ainv_x
    Ainv_y = Ainv_y
    Hu_x = Hhat_x
    Hu_y = Hhat_y
    boundaries_to_force = 'right'
  []
  [divergence]
    type = FVFunctorDivergence
    sign = -1
    x_functor = Hhat_x
    y_functor = Hhat_y
    variable = pressure_p
  []
  # [body_v]
  #   type = FVBodyForce
  #   variable = pressure_p
  #   function = '-1.0 * (max(1.0 - x, 0) * max(0.75 - y, 0.0) * max(y - 0.25, 0.0))'
  # []
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
    pressure_relaxation = 0.01
  []
  # [corrector_x]
  #   type = FVCorrector
  #   variable = u_adv
  #   execute_on = timestep_end
  #   pressure = pressure_relaxed
  #   pressure_old = pressure_old
  #   Ainv = Ainv_x
  #   Hhat = Hhat_x
  #   momentum_component = 'x'
  #   pressure_relaxation = 1.0
  # []
  # [corrector_y]
  #   type = FVCorrector
  #   variable = v_adv
  #   execute_on = timestep_end
  #   pressure = pressure_relaxed
  #   pressure_old = pressure_old
  #   Ainv = Ainv_y
  #   Hhat = Hhat_y
  #   momentum_component = 'y'
  #   pressure_relaxation = 1.0
  # []
[]

[FVBCs]
  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = 'right'
    variable = pressure_p
    function = 0
  []
  [outlet_p_zerograd]
    type = INSFVSymmetryPressureBC
    boundary = 'left top bottom'
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
  num_steps = 500
  dt = 0.1
  dtmin = 0.1
  solve_type = 'LINEAR'
  # petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  # petsc_options_value = 'asm      200                lu           NONZERO'
  # petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -snes_max_it -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol '
  # petsc_options_value = 'gmres asm lu 100 NONZERO 2 1E-14 1E-12'
  petsc_options_iname = '-ksp_type -pc_type -pc_sub_type -sub_pc_factor_levels'
  petsc_options_value = 'gmres asm ilu 4'
  momentum_predictor_bool = false
  verbose_print = true
[]

[MultiApps]
  [sub_predictor]
    #type = TransientMultiApp
    type = FullSolveMultiApp
    input_files = FV_Channel_Momentum_Predictor.i
    execute_on = TIMESTEP_BEGIN
    sub_cycling = false
  []
  # [sub_corrector]
  #   # type = TransientMultiApp
  #   type = FullSolveMultiApp
  #   input_files = FV_Channel_Momentum_Corrector.i
  #   execute_on = TIMESTEP_END
  #   sub_cycling = false
  # []
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

  # [u_to_sub_predictor]
  #   type = MultiAppCopyTransfer
  #   #type = MultiAppCopyTransfer
  #   direction = to_multiapp
  #   multi_app = sub_predictor
  #   source_variable = u_adv
  #   variable = u_adv
  # []

  # [v_to_sub_predictor]
  #   type = MultiAppCopyTransfer
  #   #type = MultiAppCopyTransfer
  #   direction = to_multiapp
  #   multi_app = sub_predictor
  #   source_variable = v_adv
  #   variable = v_adv
  # []

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

  # [Ainv_x_to_sub_corrector]
  #   type = MultiAppCopyTransfer
  #   direction = to_multiapp
  #   multi_app = sub_corrector
  #   source_variable = Ainv_x
  #   variable = Ainv_x
  # []
  #
  # [Ainv_y_to_sub_corrector]
  #   type = MultiAppCopyTransfer
  #   direction = to_multiapp
  #   multi_app = sub_corrector
  #   source_variable = Ainv_y
  #   variable = Ainv_y
  # []
  #
  # [Hhat_x_to_sub_corrector]
  #   type = MultiAppCopyTransfer
  #   direction = to_multiapp
  #   multi_app = sub_corrector
  #   source_variable = Hhat_x
  #   variable = Hhat_x
  # []
  #
  # [Hhat_y_to_sub_corrector]
  #   type = MultiAppCopyTransfer
  #   direction = to_multiapp
  #   multi_app = sub_corrector
  #   source_variable = Hhat_y
  #   variable = Hhat_y
  # []
  #
  # [p_to_sub_corrector]
  #   type = MultiAppCopyTransfer
  #   direction = to_multiapp
  #   multi_app = sub_corrector
  #   source_variable = pressure_relaxed
  #   variable = pressure_corr
  # []
  #
  # [p_old_to_sub_corrector]
  #   type = MultiAppCopyTransfer
  #   direction = to_multiapp
  #   multi_app = sub_corrector
  #   source_variable = pressure_old
  #   variable = pressure_old
  # []
  #
  # [u_adv_old_to_sub_corrector]
  #   type = MultiAppCopyTransfer
  #   direction = to_multiapp
  #   multi_app = sub_corrector
  #   source_variable = u_adv
  #   variable = u_adv_old
  # []
  #
  # [v_adv_old_to_sub_corrector]
  #   type = MultiAppCopyTransfer
  #   direction = to_multiapp
  #   multi_app = sub_corrector
  #   source_variable = v_adv
  #   variable = v_adv_old
  # []
  #
  # [u_adv_from_sub_corrector]
  #   type = MultiAppCopyTransfer
  #   direction = from_multiapp
  #   multi_app = sub_corrector
  #   source_variable = u_adv_corr
  #   variable = u_adv
  # []
  #
  # [v_adv_from_sub_corrector]
  #   type = MultiAppCopyTransfer
  #   direction = from_multiapp
  #   multi_app = sub_corrector
  #   source_variable = v_adv_corr
  #   variable = v_adv
  # []
[]

[Outputs]
  file_base = channel
  exodus = true
  checkpoint = true
  # perf_graph = true # prints a performance report to the terminal
[]
