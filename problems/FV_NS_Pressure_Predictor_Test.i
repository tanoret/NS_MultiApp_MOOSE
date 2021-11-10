mu=1.1
rho=1.1
advected_interp_method='average'
velocity_interp_method='rc'


[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 2
    ymin = 0
    ymax = 2
    nx = 2
    ny = 2
  []
[]

[Problem]
  fv_bcs_integrity_check = true
[]

[Variables]
  [pressure]
    type = INSFVPressureVariable
    initial_condition = 1
  []
[]

[AuxVariables]
  [u_star]
    type = INSFVVelocityVariable
    initial_condition = 1
  []
  [v_star]
    type = INSFVVelocityVariable
    initial_condition = 1
  []
  [pressure_old]
    type = INSFVPressureVariable
    initial_condition = 1
  []
  [Ainv_x]
    type = MooseVariableFVReal
    initial_condition = 1
  []
  [Hhat_x]
    type = MooseVariableFVReal
    initial_condition = 1
  []
  [RHS_x]
    type = MooseVariableFVReal
    initial_condition = 1
  []
  [Ainv_y]
    type = MooseVariableFVReal
    initial_condition = 1
  []
  [Hhat_y]
    type = MooseVariableFVReal
    initial_condition = 1
  []
  [RHS_y]
    type = MooseVariableFVReal
    initial_condition = 1
  []
[]

[FVKernels]
  [pressure_poisson]
    type = FVNavStokesPressurePredictor_p
    variable = pressure
    #pressure = pressure
    pressure_old = pressure_old
    u_star = u_star
    v_star = v_star
    Ainv_x = Ainv_x
    Ainv_y = Ainv_y
    Hu_x = Hhat_x
    Hu_y = Hhat_y
    rhs_x = RHS_x
    rhs_y = RHS_y
  []
[]

[AuxKernels]
  [corrector_x]
    type = FVCorrector
    variable = u_star
    execute_on = timestep_end
    pressure = pressure
    pressure_old = pressure_old
    vel_star = u_star
    Ainv = Ainv_x
    Hu = Hhat_x
    rhs = RHS_x
    momentum_component = 'x'
  []
  [corrector_y]
    type = FVCorrector
    variable = v_star
    execute_on = timestep_end
    pressure = pressure
    pressure_old = pressure_old
    vel_star = v_star
    Ainv = Ainv_y
    Hu = Hhat_y
    rhs = RHS_y
    momentum_component = 'y'
  []
[]

[FVBCs]
  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = 'right'
    variable = pressure
    function = '0'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      200                lu           NONZERO'
  line_search = 'none'
  nl_rel_tol = 1e-3
  num_steps = 1
[]

[Outputs]
  exodus = true
[]
