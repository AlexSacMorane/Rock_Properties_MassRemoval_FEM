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
  active = 'matter'

  [./matter]
    order = FIRST
    family = LAGRANGE
    outputs = exodus
    [./InitialCondition]
      type = FunctionIC
      function = matter_txt
    [../]
  [../]
[]

[Kernels]
  # Order parameter matter
  [./dmatterdt]
    type = TimeDerivative
    variable = matter
  [../]
  [./ACBulk_matter]
    type = AllenCahn
    variable = matter
    mob_name = L_matter
    f_name = g_matter
  [../]
  [./ACInterface_matter]
    type = ACInterface
    variable = matter
    mob_name = L_matter
    kappa_name = kappa_matter
  [../]
[]

[BCs]
  [./Periodic]
    [./per_matter_xy]
      variable = matter
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  [./consts]
    type = GenericConstantMaterial
    prop_names  = 'L_matter kappa_matter'
    prop_values =
  [../]

  [./ed_map]
    type = GenericFunctionMaterial
    block = 0
    prop_names = ed
    prop_values = ed_txt
    outputs = exodus
  [../]
  [./free_energy_matter]
    type = DerivativeParsedMaterial
    block = 0
    property_name = g_matter
    coupled_variables = 'matter'
    material_property_names = 'ed'
    constant_names = 'W'
    constant_expressions =
    expression = 'W*(matter^2)*((1-matter)^2) + ed*(3*matter^2-2*matter^3)'
    enable_jit = true
    derivative_order = 1
    #outputs = exodus
  [../]
[]

[Functions]
  [matter_txt]
    type = PiecewiseMultilinear
    data_file = data/matter.txt
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
  [matter_pp]
    type = ElementAverageValue
    variable = matter
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
    show = 'matter_pp'
  [../]
[]
