mu=1.0
rho=1.1
U=1.0
advected_interp_method='average'
#velocity_interp_method='rc'
velocity_interp_method='average'

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 200
    ny = 50
    xmin = 0.0
    xmax = 5.0
    ymin = 0.0
    ymax = 1.0
  []
[]

[Problem]
  fv_bcs_integrity_check = true
[]

[Variables]
  [u_adv_corr]
    type = INSFVVelocityVariable
  []
  [v_adv_corr]
    type = INSFVVelocityVariable
  []
[]

[AuxVariables]
  [pressure_corr]
    type = INSFVPressureVariable
  []
  [pressure_old]
    type = INSFVPressureVariable
  []
  [u_adv_old]
    type = INSFVVelocityVariable
  []
  [v_adv_old]
    type = INSFVVelocityVariable
  []
  [Ainv_x]
    type = MooseVariableFVReal
  []
  [Hhat_x]
    type = MooseVariableFVReal
  []
  [Ainv_y]
    type = MooseVariableFVReal
  []
  [Hhat_y]
    type = MooseVariableFVReal
  []
[]

[FVKernels]
  [corrector_x]
    type = FVCorrectorMultiApp
    variable = u_adv_corr
    pressure = pressure_corr
    pressure_old = pressure_old
    Ainv = Ainv_x
    Hhat = Hhat_x
    vel_adv_old = u_adv_old
    momentum_component = 'x'
    advection_relaxation = 1.0 # This is the gradient relaxation - try to keep in 1.0
  []
  [corrector_y]
      type = FVCorrectorMultiApp
      variable = v_adv_corr
      pressure = pressure_corr
      pressure_old = pressure_old
      Ainv = Ainv_y
      Hhat = Hhat_y
      vel_adv_old = v_adv_old
      momentum_component = 'y'
      advection_relaxation = 1.0 # This is the gradient relaxation - try to keep in 1.0
    []
  []

[FVBCs]
  [inlet-u]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = u_adv_corr
    function = ${U}
  []
  [inlet-v]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = v_adv_corr
    function = '0'
  []
  [walls-u]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom'
    variable = u_adv_corr
    function = 0
  []
  [walls-v]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom'
    variable = v_adv_corr
    function = 0
  []
  # [outlet_u]
  #   type = INSFVMomentumAdvectionOutflowBC
  #   variable = u_adv_corr
  #   advected_quantity = 'rhou'
  #   vel = 'velocity'
  #   advected_interp_method = ${advected_interp_method}
  #   u = u_adv_corr
  #   v = v_adv_corr
  #   boundary = 'right'
  # []
  # [outlet_v]
  #   type = INSFVMomentumAdvectionOutflowBC
  #   variable = v_adv_corr
  #   advected_quantity = 'rhov'
  #   vel = 'velocity'
  #   advected_interp_method = ${advected_interp_method}
  #   u = u_adv_corr
  #   v = v_adv_corr
  #   boundary = 'right'
  # []
[]

[Materials]
  [ins_fv]
    type = INSFVMaterial
    u = 'u_adv_corr'
    v = 'v_adv_corr'
    pressure = 'pressure_corr'
    rho = ${rho}
  []
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    solve_type = 'NEWTON'
  [../]
[]

[Executioner]
  type = Steady
  num_steps = 5
  # dt = .1
  # dtmin = .1

  #petsc_options_iname = '-pc_type'
  #petsc_options_value = 'lu'

  #line_search = 'none'
  nl_rel_tol = 1e-2
  l_max_its = 500
[]

# [Executioner]
#   type = Transient
#   solve_type = 'NEWTON'
#   # petsc_options = '-pc_svd_monitor'
#   # petsc_options_iname = '-pc_type'
#   # petsc_options_value = 'svd'
# []

[Outputs]
  file_base = channel_corr
  exodus = true
  checkpoint = true
[]
