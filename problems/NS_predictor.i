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
[]

[AuxVariables]
  # x-velocity
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

  # y-velocity
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

  [./u_old]
      order = SECOND
      family = L2_LAGRANGE
      #initial_from_file_var = u
      #initial_from_file_timestep = LATEST
  [../]

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

  [./p_old]
    order = FIRST
    family = LAGRANGE
    #initial_from_file_var = p
    #initial_from_file_timestep = LATEST
  [../]
[]

[Kernels]
  [./x_predictor]
    type = NavStokesPredictor_p
    variable = u_star
    u = u
    v = v
    u_star = u_star
    v_star = v_star
    p = p
    component = 0
  [../]

  [./y_predictor]
    type = NavStokesPredictor_p
    variable = v_star
    u = u
    v = v
    u_star = u_star
    v_star = v_star
    p = p
    component = 1
  [../]
[]

[AuxKernels]
  [./normalization_auxkernel]
    type = NormalizationAux
    variable = u_old
    source_variable = u
    normal_factor = 1.0
    execute_on = timestep_end
    # Note: 'normalization' or 'shift' are provided as CLI args
  [../]

  [./normalization_auxkernel2]
    type = NormalizationAux
    variable = p_old
    source_variable = p
    normal_factor = 1.0
    execute_on = timestep_end
    # Note: 'normalization' or 'shift' are provided as CLI args
  [../]
[]

[BCs]
  [./x_no_slip]
    type = DirichletBC
    variable = u_star
    boundary = 'top bottom'
    value = 0.0
  [../]

  [./y_no_slip]
    type = DirichletBC
    variable = v_star
    boundary = 'top bottom'
    value = 0.0
  [../]

  [./velocity_inlet_x]
    type = DirichletBC
    variable = u_star
    boundary = 'left'
    value = 1.0
  []

  [./velocity_inlet_y]
    type = DirichletBC
    variable = v_star
    boundary = 'left'
    value = 0.0
  []
[]

[Materials]
  [./const]
    type = GenericConstantMaterial
    #block = 'FLUID'
    #block = 'SOLID'
    prop_names = 'rho mu'
    prop_values = '1  0.002'
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
  #num_steps = 10
  #dt = .06
  #dtmin =

  #petsc_options_iname = '-pc_type'
  #petsc_options_value = 'lu'

  #petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type -sub_pc_factor_levels'
  #petsc_options_value = '300                bjacobi  ilu          4'

  petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -snes_max_it -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol '
  petsc_options_value = 'gmres asm lu 100 NONZERO 2 1E-14 1E-12'

  #petsc_options_iname = '-ksp_type -pc_type -pc_sub_type -sub_pc_factor_levels'
  #petsc_options_value = 'gmres asm ilu 4'

  #petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_levels'   #Option 2
  #petsc_options_value = 'asm      2               ilu          4'
  #line_search = 'none'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  nl_max_its = 40
  l_tol = 1e-8
  l_max_its = 500
[]
