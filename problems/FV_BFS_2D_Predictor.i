mu=0.001
rho=1.0
U=1.0
advected_interp_method='average'
#velocity_interp_method='rc'
velocity_interp_method='average'

[Mesh]
  file = BFS_2D_flat.e
  # uniform_refine = 1
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
  []
  [Hu_x]
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
  []

  [u_time_derivative_and_relax]
    type = FVNavierStokesTimeRelax
    variable = u
    add_time_derivative = false
    velocity_relaxation = 0.7
    Ainv = Ainv_x
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
  []

  [v_time_derivative_and_relax]
    type = FVNavierStokesTimeRelax
    variable = v
    add_time_derivative = false
    velocity_relaxation = 0.7
    Ainv = Ainv_y
  []

[]

[FVBCs]
  [inlet-u]
    type = INSFVInletVelocityBC
    boundary = 'Inlet'
    variable = u
    function = ${U}
  []
  [inlet-v]
    type = INSFVInletVelocityBC
    boundary = 'Inlet'
    variable = v
    function = '0'
  []
  [walls-u]
    type = INSFVNoSlipWallBC
    boundary = 'Wall'
    variable = u
    function = 0
  []
  [walls-v]
    type = INSFVNoSlipWallBC
    boundary = 'Wall'
    variable = v
    function = 0
  []
  # [outlet_p]
  #   type = INSFVOutletPressureBC
  #   boundary = 'Outlet'
  #   variable = pressure_mom
  #   function = 0
  # []
  [outlet_u]
    type = INSFVMomentumAdvectionOutflowBC
    variable = u
    advected_quantity = 'rhou'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    u = u_adv
    v = v_adv
    boundary = 'Outlet'
  []
  [outlet_v]
    type = INSFVMomentumAdvectionOutflowBC
    variable = v
    advected_quantity = 'rhov'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    u = u_adv
    v = v_adv
    boundary = 'Outlet'
  []
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

  petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -snes_max_it -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol '
  petsc_options_value = 'gmres asm lu 100 NONZERO 2 1E-14 1E-12'

  # petsc_options_iname = '-ksp_type -pc_type -pc_sub_type -sub_pc_factor_levels'
  # petsc_options_value = 'gmres asm ilu 4'

  #petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_levels'   #Option 2
  #petsc_options_value = 'asm      2               ilu          4'
  #line_search = 'none'

  nl_rel_tol = 1e-2
  nl_abs_tol = 1e-6
  nl_max_its = 40
  l_tol = 1e-2
  l_max_its = 100

  momentum_predictor_bool = true
  verbose_print = false
[]

[Outputs]
  file_base = BFS_2D_sub
  exodus = true
  # perf_graph = true # prints a performance report to the terminal
[]
