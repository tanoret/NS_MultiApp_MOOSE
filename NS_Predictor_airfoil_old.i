[Mesh]
  second_order = true
  [fmg]
    type = FileMeshGenerator
    file = AirfoilMeshes/Re100MediumMesh.exo
  []
[]
#[Mesh]
#  file = NACA_airfoil_Pred.e
#[]

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
      family = LAGRANGE
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
    family = LAGRANGE
    #initial_from_file_var = v
    #initial_from_file_timestep = LATEST

    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]

  [./u_old]
      order = SECOND
      family = LAGRANGE
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
    boundary = 'WALL AIRFOIL'
    value = 0.0
  [../]

  [./y_no_slip]
    type = DirichletBC
    variable = v_star
    boundary = 'WALL AIRFOIL'
    value = 0.0
  [../]

  [./velocity_inlet_x]
    type = DirichletBC
    variable = u_star
    boundary = 'INLET'
    value = 1.0
  []

  [./velocity_inlet_y]
    type = DirichletBC
    variable = v_star
    boundary = 'INLET'
    value = 0.0
  []
[]

[Materials]
  [./const]
    type = GenericConstantMaterial
    block = 'FLUID'
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
  #num_steps = 50
  #dt = .1
  #dtmin = .1

  #petsc_options_iname = '-pc_type'
  #petsc_options_value = 'lu'

  #petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type -sub_pc_factor_levels'
  #petsc_options_value = '300                bjacobi  ilu          4'

  petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -snes_max_it -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol '
  petsc_options_value = 'gmres asm lu 100 NONZERO 2 1E-14 1E-12'

  #petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_levels'   #Option 2
  #petsc_options_value = 'asm      2               ilu          4'
  #line_search = 'none'

  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-8
  nl_max_its = 40
  l_tol = 1e-6
  l_max_its = 500
[]

#[Outputs]
#  file_base = NACA_airfoil_Pred
#  exodus = true
#  checkpoint = true
#[]
