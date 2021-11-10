mu=1.1
rho=1.1
advected_interp_method='average'
velocity_interp_method='rc'

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 5
    ymin = -1
    ymax = 1
    nx = 50
    ny = 10
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
    initial_condition = 0
  []
  [v_star]
    type = INSFVVelocityVariable
    initial_condition = 0
  []
  [p_old]
    type = INSFVPressureVariable
    initial_condition = 0
  []
[]

[FVKernels]
  [pressure_correction]
    type = FVNavStokesPressurePoisson_p
    variable = pressure
    u_star = u_star
    v_star = v_star
    pressure_old = p_old
    mu = ${mu}
    rho = ${rho}
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

# [Materials]
#   [ins_fv]
#     type = INSFVMaterial
#     u = 'u_star'
#     v = 'v_star'
#     pressure = 'pressure'
#     rho = ${rho}
#   []
# []

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      200                lu           NONZERO'
  line_search = 'none'
  nl_rel_tol = 1e-12
[]

[Outputs]
  exodus = true
[]
