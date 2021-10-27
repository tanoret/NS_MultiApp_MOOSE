[Mesh]
  second_order = true
  [fmg]
    type = FileMeshGenerator
    file = NACA_airfoil.exo
  []
[]

[AuxVariables]
  [vel_x]
    order = SECOND
  []
  [vel_y]
    order = SECOND
  []
[]

[AuxKernels]
  [vel_x]
    type = VectorVariableComponentAux
    variable = vel_x
    vector_variable = velocity
    component = 'x'
  []
  [vel_y]
    type = VectorVariableComponentAux
    variable = vel_y
    vector_variable = velocity
    component = 'y'
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
  [../]

  [./momentum_pressure]
    type = INSADMomentumPressure
    variable = velocity
    p = p
    integrate_p_by_parts = true
  [../]
[]

[BCs]
  [./no_slip]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = 'WALL AIRFOIL'
  [../]

  [./velocity_inlet]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = 'INLET'
    function_x = 1.0
    function_y = 0.0
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
    prop_values = '1  0.0004'
  [../]
  [ins_mat]
    type = INSADMaterial
    velocity = velocity
    pressure = p
  [../]
[]

#[Preconditioning]
#  active = FSP
#  [./FSP]
#    type = FSP
#    # It is the starting point of splitting
#    topsplit = 'up' # 'up' should match the following block name
#    [./up]
#      splitting = 'u p' # 'u' and 'p' are the names of subsolvers
#      splitting_type  = schur
#      petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition'
#      petsc_options_value = 'full selfp'
#    [../]
#    [./u]
#      vars = 'velocity'
#      petsc_options_iname = '-pc_type -ksp_type -ksp_rtol'
#      petsc_options_value = '     hypre gmres 1e-4'
#    [../]
#    [./p]
#      vars = 'p'
#      petsc_options = '-ksp_monitor'
#      petsc_options_iname = '-pc_type -ksp_type -ksp_rtol'
#      petsc_options_value = '   jacobi    gmres     1e-4'
#    [../]
#  [../]
#[]

[Preconditioning]
  active = FSP
  [./FSP]
    type = FSP
    # It is the starting point of splitting
    topsplit = 'up' # 'up' should match the following block name
    [./up]
      splitting = 'u p' # 'u' and 'p' are the names of subsolvers
      splitting_type  = schur
      # Splitting type is set as schur, because the pressure part of Stokes-like systems
      # is not diagonally dominant. CAN NOT use additive, multiplicative and etc.
      # Original system:
      # | A B | | u | = | f_u |
      # | C 0 | | p |   | f_v |
      # is factorized into
      # |I        0 | | A    0|  | I  A^{-1}B | | u | = | f_u |
      # |CA^{-1}  I | | 0   -S|  | 0    I     | | p |   | f_v |
      # S = CA^{-1}B
      # The preconditioning is accomplished via the following steps
      # (1) p^{(0)} = f_v - CA^{-1}f_u,
      # (2) p = (-S)^{-1} p^{(0)}
      # (3) u = A^{-1}(f_u-Bp)
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition'
      petsc_options_value = 'full selfp'
      # Factorization type here is full, which means we approximate the original system
      # exactly. There are three other options:
      # diag:
      # | A    0|
      # | 0   -S|
      # lower:
      # |I        0  |
      # |CA^{-1}  -S |
      # upper:
      # | I  A^{-1}B |
      # | 0    -S    |
      # The preconditioning matrix is set as selfp, which means we explicitly form a
      # matrix \hat{S} = C(diag(A))^{-1}B. We do not compute the inverse of A, but instead, we compute
      # the inverse of diag(A).
    [../]
    [./u]
      vars = 'velocity'
      # PETSc options for this subsolver
      # A prefix will be applied, so just put the options for this subsolver only
      petsc_options_iname = '-pc_type -ksp_type -ksp_rtol'
      petsc_options_value = '     hypre gmres 1e-4'
      # Specify options to solve A^{-1} in the steps (1), (2) and (3).
      # Solvers for A^{-1} could be different in different steps. We could
      # choose in the following pressure block.
    [../]
    [./p]
      vars = 'p'
      # PETSc options for this subsolver in the step (2)
      petsc_options = '-ksp_monitor'
      petsc_options_iname = '-pc_type -ksp_type -ksp_rtol'
      petsc_options_value = '   lu    gmres     1e-4'
      # Use -inner_ksp_type and -inner_pc_type to override A^{-1} in the step (2)
      # Use -lower_ksp_type and -lower_pc_type to override A^{-1} in the step (1)
    [../]
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  num_steps = 20
  dt = .1
  dtmin = .1
  line_search = 'basic'
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-13
  nl_max_its = 6
  l_tol = 1e-6
  l_max_its = 500
[]


[Outputs]
  file_base = NACA_airfoil
  exodus = true
[]
