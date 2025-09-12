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
    generate_output = 'stress_xx stress_xy stress_xz stress_yx stress_yy stress_yz stress_zx stress_zy stress_zz
                       strain_xx strain_xy strain_xz strain_yx strain_yy strain_yz strain_zx strain_zy strain_zz'
  []
[]

[BCs]
  [left_x]  
    type = DirichletBC
    variable = disp_x
    boundary = left 
    value = 0
  []
  [right_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right 
    function = 
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
  [back_z]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
  [front_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = 
  []
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
  
  end_time = 
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
  [strain_xx_gra_pp]
    type = ElementAverageValue
    variable = strain_xx
    block = 1
  []
  [strain_xy_gra_pp]
    type = ElementAverageValue
    variable = strain_xy
    block = 1
  []
  [strain_xz_gra_pp]
    type = ElementAverageValue
    variable = strain_xz
    block = 1
  []
  [strain_yy_gra_pp]
    type = ElementAverageValue
    variable = strain_yy
    block = 1
  []
  [strain_yz_gra_pp]
    type = ElementAverageValue
    variable = strain_yz
    block = 1
  []
  [strain_zz_gra_pp]
    type = ElementAverageValue
    variable = strain_zz
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
  [strain_xx_cem_pp]
    type = ElementAverageValue
    variable = strain_xx
    block = 2
  []
  [strain_xy_cem_pp]
    type = ElementAverageValue
    variable = strain_xy
    block = 2
  []
  [strain_xz_cem_pp]
    type = ElementAverageValue
    variable = strain_xz
    block = 2
  []
  [strain_yy_cem_pp]
    type = ElementAverageValue
    variable = strain_yy
    block = 2
  []
  [strain_yz_cem_pp]
    type = ElementAverageValue
    variable = strain_yz
    block = 2
  []
  [strain_zz_cem_pp]
    type = ElementAverageValue
    variable = strain_zz
    block = 2
  []
  [disp_x_pp]
    type = NormalBoundaryDisplacement
    boundary = right
    displacements = 'disp_x disp_y disp_z'
    value_type = average
  []
  [disp_y_pp]
    type = NormalBoundaryDisplacement
    boundary = top
    displacements = 'disp_x disp_y disp_z'
    value_type = average
  []
  [disp_z_pp]
    type = NormalBoundaryDisplacement
    boundary = front
    displacements = 'disp_x disp_y disp_z'
    value_type = average
  []
[]

[Outputs]
  exodus = true
  [console]
    type = Console
    execute_on = 'nonlinear'
    max_rows = 3
    show = 'disp_x_pp disp_y_pp disp_z_pp 
            stress_xx_gra_pp stress_xx_cem_pp stress_yy_gra_pp stress_yy_cem_pp stress_zz_gra_pp stress_zz_cem_pp'
  []
  [./csv]
    type = CSV
    show = 'disp_x_pp disp_y_pp disp_z_pp
            stress_xx_gra_pp stress_xy_gra_pp stress_xz_gra_pp stress_yy_gra_pp stress_yz_gra_pp stress_zz_gra_pp stress_xx_cem_pp stress_xy_cem_pp stress_xz_cem_pp stress_yy_cem_pp stress_yz_cem_pp stress_zz_cem_pp
            strain_xx_gra_pp strain_xy_gra_pp strain_xz_gra_pp strain_yy_gra_pp strain_yz_gra_pp strain_zz_gra_pp strain_xx_cem_pp strain_xy_cem_pp strain_xz_cem_pp strain_yy_cem_pp strain_yz_cem_pp strain_zz_cem_pp'
  [../]
[]