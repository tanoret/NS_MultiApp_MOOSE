mu=1.1
rho=1.0
U=1.0
advected_interp_method='average' #'average'
#velocity_interp_method='rc'
velocity_interp_method='average'

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
  [u]
    type = INSFVVelocityVariable
    initial_condition = 0
  []
  [v]
    type = INSFVVelocityVariable
    initial_condition = 0
  []
[]

[AuxVariables]
  [u_adv]
    type = INSFVVelocityVariable
    initial_condition = 0
  []
  [v_adv]
    type = INSFVVelocityVariable
    initial_condition = 0
  []
  [pressure_mom]
    type = INSFVPressureVariable
  []
  [Ainv_x]
    type = MooseVariableFVReal
    initial_condition = 0
  []
  [Ainv_x_old]
    type = MooseVariableFVReal
    initial_condition = 0
  []
  [Hu_x]
    type = MooseVariableFVReal
  []
  [RHS_x]
    type = MooseVariableFVReal
  []
  [Ainv_y]
    type = MooseVariableFVReal
    initial_condition = 0
  []
  [Ainv_y_old]
    type = MooseVariableFVReal
    initial_condition = 0
  []
  [Hu_y]
    type = MooseVariableFVReal
  []
  [RHS_y]
    type = MooseVariableFVReal
  []
[]

[FVKernels]
  [u_adv_diff_residual]
    type = FVNavStokesPredictor_p
    variable = u
    advected_quantity = 'rhou'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    pressure = pressure_mom
    u = u_adv
    v = v_adv
    mu = ${mu}
    rho = ${rho}
    momentum_component = 'x'
    velocity_relaxation = 1.0
  []

  [u_time_derivative_and_relax]
    type = FVNavierStokesTimeRelax
    variable = u
    add_time_derivative = true
    velocity_relaxation = 0.5
    Ainv = Ainv_x
    Ainv_old = Ainv_x_old
  []

  [v_adv_diff_residual]
    type = FVNavStokesPredictor_p
    variable = v
    advected_quantity = 'rhov'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    pressure = pressure_mom
    u = u_adv
    v = v_adv
    mu = ${mu}
    rho = ${rho}
    momentum_component = 'y'
    velocity_relaxation = 1.0
  []

  [v_time_derivative_and_relax]
    type = FVNavierStokesTimeRelax
    variable = v
    add_time_derivative = true
    velocity_relaxation = 0.5
    Ainv = Ainv_y
    Ainv_old = Ainv_y_old
  []

  # [u_advection]
  #   type = INSFVMomentumAdvection
  #   variable = u
  #   advected_quantity = 'rhou'
  #   vel = 'velocity'
  #   advected_interp_method = ${advected_interp_method}
  #   velocity_interp_method = ${velocity_interp_method}
  #   pressure = pressure_mom
  #   u = u
  #   v = v
  #   mu = ${mu}
  #   rho = ${rho}
  # []
  # [u_viscosity]
  #   type = FVDiffusion
  #   variable = u
  #   coeff = ${mu}
  # []
  # [u_pressure]
  #   type = INSFVMomentumPressure
  #   variable = u
  #   momentum_component = 'x'
  #   pressure = pressure_mom
  # []
  #
  # [v_advection]
  #   type = INSFVMomentumAdvection
  #   variable = v
  #   advected_quantity = 'rhov'
  #   vel = 'velocity'
  #   advected_interp_method = ${advected_interp_method}
  #   velocity_interp_method = ${velocity_interp_method}
  #   pressure = pressure_mom
  #   u = u
  #   v = v
  #   mu = ${mu}
  #   rho = ${rho}
  # []
  # [v_viscosity]
  #   type = FVDiffusion
  #   variable = v
  #   coeff = ${mu}
  # []
  # [v_pressure]
  #   type = INSFVMomentumPressure
  #   variable = v
  #   momentum_component = 'y'
  #   pressure = pressure_mom
  # []

[]

[AuxKernels]
  [corrector_x]
    type = FVCorrectorPredictor
    variable = u_adv
    execute_on = timestep_begin
    pressure = pressure_mom
    Ainv = Ainv_x
    Hhat = Hu_x
    momentum_component = 'x'
    pressure_relaxation = 0.3
  []
  [corrector_y]
    type = FVCorrectorPredictor
    variable = v_adv
    execute_on = timestep_begin
    pressure = pressure_mom
    Ainv = Ainv_y
    Hhat = Hu_y
    momentum_component = 'y'
    pressure_relaxation = 0.3
  []
  [copy_Ainv_x]
    type = FVCopyKernel
    variable = Ainv_x_old
    execute_on = timestep_end
    FVVar = Ainv_x
  []
  [copy_Ainv_y]
    type = FVCopyKernel
    variable = Ainv_y_old
    execute_on = timestep_end
    FVVar = Ainv_y
  []
[]

[FVBCs]
  [inlet-u]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = u
    function = ${U}
  []
  [inlet-v]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = v
    function = '0'
  []
  [walls-u]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom'
    variable = u
    function = 0
  []
  [walls-v]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom'
    variable = v
    function = 0
  []
  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = 'right'
    variable = pressure_mom
    function = 0
  []
  # [outlet_u]
  #   type = INSFVMomentumAdvectionOutflowBC
  #   variable = u
  #   advected_quantity = 'rhou'
  #   vel = 'velocity'
  #   advected_interp_method = ${advected_interp_method}
  #   u = u_adv
  #   v = v_adv
  #   boundary = 'right'
  # []
  # [outlet_v]
  #   type = INSFVMomentumAdvectionOutflowBC
  #   variable = v
  #   advected_quantity = 'rhov'
  #   vel = 'velocity'
  #   advected_interp_method = ${advected_interp_method}
  #   u = u_adv
  #   v = v_adv
  #   boundary = 'right'
  # []
[]

[Materials]
  [ins_fv]
    type = INSFVMaterial
    u = 'u_adv'
    v = 'v_adv'
    pressure = 'pressure_mom'
    rho = ${rho}
  []
[]

[Executioner]
  type = CustomTransient
  num_steps = 1
  dt = 0.1
  dtmin = 0.1

  solve_type = 'LINEAR'

  #petsc_options_iname = '-pc_type'
  #petsc_options_value = 'lu'

  #petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type -sub_pc_factor_levels'
  #petsc_options_value = '300                bjacobi  ilu          4'

  # petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -snes_max_it -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol '
  # petsc_options_value = 'gmres asm lu 100 NONZERO 2 1E-14 1E-12'

  petsc_options_iname = '-ksp_type -pc_type -pc_sub_type -sub_pc_factor_levels'
  petsc_options_value = 'gmres asm ilu 4'

  #petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_levels'   #Option 2
  #petsc_options_value = 'asm      2               ilu          4'
  #line_search = 'none'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  nl_max_its = 40
  l_tol = 1e-8
  l_max_its = 500

  momentum_predictor_bool = true
  verbose_print = false
  velocity_relaxation = 1.0
[]

[Outputs]
  file_base = channel_sub
  exodus = true
  # perf_graph = true # prints a performance report to the terminal
[]
