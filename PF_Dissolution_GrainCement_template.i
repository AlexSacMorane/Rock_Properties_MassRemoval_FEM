[Mesh]
  type = GeneratedMesh
  dim = 2
  nx =
  ny =
  nz = 0
  xmin =
  xmax =
  ymin =
  ymax =
  elem_type = QUAD4
[]

[Variables]
  active = 'grain cement'

  [./grain]
    order = FIRST
    family = LAGRANGE
    outputs = exodus
    [./InitialCondition]
      type = FunctionIC
      function = grain_txt
    [../]
  [../]
  [./cement]
    order = FIRST
    family = LAGRANGE
    outputs = exodus
    [./InitialCondition]
      type = FunctionIC
      function = cement_txt
    [../]
  [../]
[]

[Kernels]
  # Order parameter grain
  [./dgraindt]
    type = TimeDerivative
    variable = grain
  [../]
  [./ACBulk_grain]
    type = AllenCahn
    variable = grain
    mob_name = L_grain
    f_name = g_grain
  [../]
  [./ACInterface_grain]
    type = ACInterface
    variable = grain
    mob_name = L_grain
    kappa_name = kappa_grain
  [../]

  # Order parameter cement
  [./dcementdt]
    type = TimeDerivative
    variable = cement
  [../]
  [./ACBulk_cement]
    type = AllenCahn
    variable = cement
    mob_name = L_cement
    f_name = g_cement
  [../]
  [./ACInterface_cement]
    type = ACInterface
    variable = cement
    mob_name = L_cement
    kappa_name = kappa_cement
  [../]
[]

[BCs]
  [./Periodic]
    [./per_grain_xy]
      variable = grain
      auto_direction = 'x y'
    [../]
    [./per_cement_xy]
      variable = cement
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  [./consts]
    # L_cement or L_grain can be changed to play on the influence of the dissolution kinetics
    type = GenericConstantMaterial
    prop_names  = 'L_cement kappa_cement L_grain kappa_grain'
    prop_values =
  [../]

  [./ed_map]
    type = GenericFunctionMaterial
    block = 0
    prop_names = ed
    prop_values = ed_txt
    outputs = exodus
  [../]
  [./free_energy_grain]
    type = DerivativeParsedMaterial
    block = 0
    property_name = g_grain
    coupled_variables = 'grain'
    material_property_names = 'ed'
    constant_names = 'W'
    constant_expressions =
    expression = 'W*(grain^2)*((1-grain)^2) + ed*(3*grain^2-2*grain^3)'
    enable_jit = true
    derivative_order = 1
    #outputs = exodus
  [../]
  [./free_energy_cement]
    type = DerivativeParsedMaterial
    block = 0
    property_name = g_cement
    coupled_variables = 'cement'
    material_property_names = 'ed'
    constant_names = 'W'
    constant_expressions =
    expression = 'W*(cement^2)*((1-cement)^2) + ed*(3*cement^2-2*cement^3)'
    enable_jit = true
    derivative_order = 2
    #outputs = exodus
  [../]
[]

[Functions]
  [grain_txt]
    type = PiecewiseMultilinear
    data_file = data/grain.txt
  []
  [cement_txt]
    type = PiecewiseMultilinear
    data_file = data/cement.txt
  []
  [ed_txt]
    type = PiecewiseMultilinear
    data_file = data/ed.txt
  []
[]

[Preconditioning]
  # This preconditioner makes sure the Jacobian Matrix is fully populated. Our
  # kernels compute all Jacobian matrix entries.
  # This allows us to use the Newton solver below.
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  # Automatic differentiation provides a _full_ Jacobian in this example
  # so we can safely use NEWTON for a fast solve
  solve_type = 'NEWTON'

  l_max_its = 20
  l_tol = 
  l_abs_tol = 

  nl_max_its = 10
  nl_rel_tol = 
  nl_abs_tol = 

  start_time = 0.0
  num_steps = 
  end_time = 

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 
  [../]
[]

[Postprocessors]
  [grain_pp]
    type = ElementAverageValue
    variable = grain
  []
  [cement_pp]
    type = ElementAverageValue
    variable = cement
  []
[]

[Outputs]
  execute_on = 'initial timestep_end'
  exodus = true
  [./other]
    type = VTK
    execute_on = 'TIMESTEP_END'
  [../]
  [console]
    type = Console
    execute_on = 'nonlinear'
    all_variable_norms = true
  []
  [./csv]
    type = CSV
    execute_on = 'TIMESTEP_END'
    show = 'grain_pp cement_pp'
  [../]
[]
