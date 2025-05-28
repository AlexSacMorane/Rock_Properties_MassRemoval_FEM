[Mesh]
  [generated]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 
    ny = 
    nz = 
    xmin = 
    xmax = 
    ymin = 
    ymax = 
    zmin = 
    zmax = 
  []
  # assign three subdomains
  [sub_domains]
    input = generated
    type = ImageSubdomainGenerator
    file_base = data/microstructure
    file_range = 
    file_suffix = 'png'
    component = 0
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Modules/TensorMechanics/Master]
  [all]
    add_variables = true
    generate_output = 'stress_xx stress_xy stress_xz stress_yx stress_yy stress_yz stress_zx stress_zy stress_zz'
  []
[]

[BCs]
  [bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = bottom
    value = 0
  []
  [top_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = top
    function = 
  []
  [./Periodic]
    [./per_x_x]
      variable = disp_x
      auto_direction = 'x'
    [../]
    [./per_x_y]
      variable = disp_x
      auto_direction = 'y'
    [../]
    [./per_y_x]
      variable = disp_y
      auto_direction = 'x'
    [../]
    [./per_y_y]
      variable = disp_y
      auto_direction = 'y'
    [../]
    [./per_z_x]
      variable = disp_z
      auto_direction = 'x'
    [../]
    [./per_z_y]
      variable = disp_z
      auto_direction = 'y'
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
  
  end_time = 0.05 # 1
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
  [stress_xz_gra_pp]
    type = ElementAverageValue
    variable = stress_xz
    block = 1
  []
  [stress_yy_gra_pp]
    type = ElementAverageValue
    variable = stress_yy
    block = 1
  []
  [stress_yz_gra_pp]
    type = ElementAverageValue
    variable = stress_yz
    block = 1
  []
  [stress_zz_gra_pp]
    type = ElementAverageValue
    variable = stress_zz
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
  [stress_xz_cem_pp]
    type = ElementAverageValue
    variable = stress_xz
    block = 2
  []
  [stress_yy_cem_pp]
    type = ElementAverageValue
    variable = stress_yy
    block = 2
  []
  [stress_yz_cem_pp]
    type = ElementAverageValue
    variable = stress_yz
    block = 2
  []
  [stress_zz_cem_pp]
    type = ElementAverageValue
    variable = stress_zz
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
    show = 'stress_xx_gra_pp stress_xy_gra_pp stress_xz_gra_pp stress_yy_gra_pp stress_yz_gra_pp stress_zz_gra_pp stress_xx_cem_pp stress_xy_cem_pp stress_xz_cem_pp stress_yy_cem_pp stress_yz_cem_pp stress_zz_cem_pp'
  [../]
[]