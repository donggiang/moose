[GlobalParams]
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = false
[]

[Mesh]
   #file = edge_crack.e
  file =  square3_3.e
[]




[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
 [./saved_x]
  [../]
  [./saved_y]
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./TensorMechanics]
  use_displaced_mesh = false
    displacements = 'disp_x disp_y'
        add_variables = true
    incremental = true
    generate_output = 'stress_xx stress_yy'
    use_automatic_differentiation = true
  [../]
[]


[Functions]
  [./rampConstant]
    type = PiecewiseLinear
    x = '0. 1e-12 1.0'
    y = '1.0 1.0 1.0'
    scale_factor = -1e2
  [../]
 [dt_func]
    type = PiecewiseLinear
    x = '0.    1e-12 1e20'
    y = '0.1  0.1 0.1'
  [../]
[]

[BCs]
  [./fix]
    type = DirichletBC
    boundary = 'top bottom fix'
    variable = disp_x
    value = 0.0
  [../]
  [./fixy]
    type = DirichletBC
    boundary = 'fix'
    variable = disp_y
    value = 0.0
  [../]
  [./Pressure]
    [./Side1]
      boundary = 'top'
      function = rampConstant
    [../]
    [./Side2]
      boundary = 'bottom'
      function = rampConstant
    [../]
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
  youngs_modulus = 200000
  poissons_ratio = 0.3
  [../]
[]

[Materials]
  [./strain]
    type = ADComputeGreenLagrangeStrain
  #  displacements = 'disp_x disp_y'
  #crack_front_definition = crackFrontDefinition

  #  enrichment_displacements = 'enrich1_x enrich2_x enrich3_x enrich4_x enrich1_y enrich2_y enrich3_y enrich4_y'
  [../]
  [./stress]
    type = ADComputeLinearElasticStress
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'Newton'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  line_search = 'none'
  [./Quadrature]
    type = GAUSS
    order =SECOND
  [../]

#  [./Predictor]
#    type = SimplePredictor
#    scale = 0.0
#  [../]

  l_max_its = 50
  nl_max_its = 500

  nl_abs_tol = 1e-0
  nl_rel_tol = 1e-1
  #l_tol=1e-3

# time control
  start_time = 0.0
  dt = 0.1
  end_time =0.3

  max_xfem_update = 1
[]


[Outputs]
  file_base = CCG_en_cr_cl_h0p16_dt2_rd0p59_ri2_ro6p5_0p18em3_e1p35_c6em24_n6p8_v2
  exodus = true
  csv = true
  execute_on = timestep_end
  [./console]
    type = Console
    output_linear = true
  [../]
[]
