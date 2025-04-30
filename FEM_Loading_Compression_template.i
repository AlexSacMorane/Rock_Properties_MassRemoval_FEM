[Mesh]
  [generated]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 
    ny = 
    xmin = 
    xmax = 
    ymin = 
    ymax = 
  []
  # assign three subdomains
  [sub_domains]
    input = generated
    type = ImageSubdomainGenerator
    file = data/microstructure.png
    component = 0
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Modules/TensorMechanics/Master]
  [all]
    add_variables = true
    generate_output = 'stress_xx stress_xy stress_yx stress_yy'
  []
[]

[BCs]
  [bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = bottom
    value = 0
  []
  [bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [top_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = 
  []
  [./Periodic]
    [./per_x_x]
      variable = disp_x
      auto_direction = 'x'
    [../]
    [./per_y_x]
      variable = disp_y
      auto_direction = 'x'
    [../]
  [../]
[]

[Materials]
  # pore
  [./pore_elastic]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus =
    poissons_ratio =
    block = 0
  [../]
  [./pore_stress_elastic]
    type = ComputeLinearElasticStress
    block = 0
  [../]
  # grain
  [./grain_elastic]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus =
    poissons_ratio =
    block = 1
  [../]
  [./grain_stress_elastic]
    type = ComputeLinearElasticStress
    block = 1
  [../]
  # cement
  [./cement_elastic]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 
    poissons_ratio = 
    block = 2
  [../]
  [./cement_stress_elastic]
    type = ComputeLinearElasticStress
    block = 2
  [../]
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient

  l_max_its  = 20
  l_tol      = 
  nl_max_its = 30
  nl_rel_tol = 
  nl_abs_tol = 
  
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  
  end_time = 1
  dt = 
[]

[Postprocessors]
  # grain
  [stress_xx_gra_pp]
    type = ElementAverageValue
    variable = stress_xx
    block = 1
  []
  [stress_xy_gra_pp]
    type = ElementAverageValue
    variable = stress_xy
    block = 1
  []
  [stress_yy_gra_pp]
    type = ElementAverageValue
    variable = stress_yy
    block = 1
  []
  # cement
  [stress_xx_cem_pp]
    type = ElementAverageValue
    variable = stress_xx
    block = 2
  []
  [stress_xy_cem_pp]
    type = ElementAverageValue
    variable = stress_xy
    block = 2
  []
  [stress_yy_cem_pp]
    type = ElementAverageValue
    variable = stress_yy
    block = 2
  []
[]

[Outputs]
  exodus = true
  [console]
    type = Console
    execute_on = 'nonlinear'
  []
  [./csv]
    type = CSV
    show = 'stress_xx_gra_pp stress_xy_gra_pp stress_yy_gra_pp stress_xx_cem_pp stress_xy_cem_pp stress_yy_cem_pp'
  [../]
[]