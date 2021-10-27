[Mesh]
  second_order = true
  [fmg]
    type = FileMeshGenerator
    file = BFSMeshes/exodus.exo   #Mesh is in mm! REMEMBER
  []
[]

[Variables]
  [./velocity]
    order = SECOND
    family = LAGRANGE_VEC
  [../]

  [./p]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[ICs]
  [./velocity]
    type = VectorConstantIC
    variable = velocity
    x_value = 100   #mm/s, max velocity is ~800mm/s
    y_value = 0
    z_value = 0
  [../]
  [./p]
    type = ConstantIC
    variable = p
    value = 0
  [../]
[]

[Kernels]
  [./mass]
    type = INSADMass
    variable = p
  [../]

  [./momentum_time]
    type = INSADMomentumTimeDerivative
    variable = velocity
  [../]

  [./momentum_convection]
    type = INSADMomentumAdvection
    variable = velocity
  [../]

  [./momentum_viscous]
    type = INSADMomentumViscous
    variable = velocity
    viscous_form = 'traction'
  [../]

  [./momentum_pressure]
    type = INSADMomentumPressure
    variable = velocity
    p = p
    integrate_p_by_parts = false
  [../]
[]

[BCs]
  [./no_slip]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = 'WALL'
  [../]

  [./velocity_inlet]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = 'INLET'
    function_x = 'Inlet_Function'
    function_y = 0.0
    function_z = 0.0
  []

  [./pressure_outlet]
    type = DirichletBC
    variable = p
    boundary = 'OUTLET'
    value = 0
  [../]
[]

[Materials]
  [./const]
    type = ADGenericConstantMaterial
    block = 'FLUID'
    prop_names = 'rho mu'
    prop_values = '1  14.8'
    #use_legacy_material_output = false
  [../]
  [ins_mat]
    type = INSADMaterial
    velocity = velocity
    pressure = p
  [../]
[]

[Functions]
  [./Inlet_Function]
    type = ParsedFunction
    value = '872.75*( (1 - (cosh(pi*z/5.2))/(cosh(pi*90/5.2)))*(cos(pi*(y-7.5)/5.2)) -
                      (1 - (cosh(3*pi*z/5.2))/(cosh(3*pi*90/5.2)))*(cos(3*pi*(y-7.5)/5.2))/(3^3) +
                      (1 - (cosh(5*pi*z/5.2))/(cosh(5*pi*90/5.2)))*(cos(5*pi*(y-7.5)/5.2))/(5^3) )'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    solve_type = 'NEWTON'
  [../]
[]

[Executioner]
  type = Transient
  # scheme = bdf2
  num_steps = 1
  dt = .001
  dtmin = .001

  #petsc_options_iname = '-pc_type'
  #petsc_options_value = 'lu'

  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm         lu                    NONZERO'

  #petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -snes_max_it -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol '
  #petsc_options_value = 'gmres asm lu 100 NONZERO 2 1E-14 1E-12'

  #petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_levels'
  #petsc_options_value = 'asm      2               ilu          4'

  #line_search = 'none'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  nl_max_its = 6
  l_tol = 1e-6
  l_max_its = 500
[]

[Outputs]
  file_base = BFS_results
  exodus = true
[]
