[GlobalParams]
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = false
[]
[XFEM]
  geometric_cut_userobjects = 'cut_mesh2'
  #qrule = volfrac
  qrule = moment_fitting
  output_cut_plane = true
  use_crack_tip_enrichment =true
  use_AD=true
  crack_front_definition = crackFrontDefinition
  enrichment_displacements = 'enrich1_x enrich2_x enrich3_x enrich4_x enrich1_y enrich2_y enrich3_y enrich4_y'
  displacements = 'disp_x disp_y'
  cut_off_boundary = all
  cut_off_radius =0.3726 #0.10984 #0.3004 #
[]



[Mesh]
   #file = edge_crack.e
  file =  square3_3.e
[]
[DomainIntegral]
  integrals = 'JIntegral InteractionIntegralKI InteractionIntegralKII'
  displacements = 'disp_x disp_y'
  crack_front_points_provider = cut_mesh2
  2d = true
  number_points_from_provider = 1
  crack_direction_method = CurvedCrackFront
  radius_inner = '0.3'
  radius_outer = '0.5'
  youngs_modulus = 200000
  poissons_ratio = 0.3
  output_q = true
  output_vpp = false
  incremental = true
  used_by_xfem_to_grow_crack = false
  use_automatic_differentiation = true
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
[AuxKernels]
  [./stress_xx]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
  [./stress_yy]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./stress_xy]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
    execute_on = timestep_end
  [../]
  [./vonmises]
    type = ADRankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    scalar_type = vonmisesStress
    execute_on = timestep_end
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
    type =ADComputeCrackTipEnrichmentIncrementalStrain #ADComputeGreenLagrangeStrain
    displacements = 'disp_x disp_y'
    crack_front_definition = crackFrontDefinition
    enrichment_displacements = 'enrich1_x enrich2_x enrich3_x enrich4_x enrich1_y enrich2_y enrich3_y enrich4_y'
  [../]
    [./radial_return_stress]
      type = ADComputeMultipleInelasticStress
      inelastic_models = 'powerlawcrp'
    [../]
    [./powerlawcrp]
      type = ADPowerLawCreepStressUpdate
      coefficient =2e-23 #6e-24#2e-23 #
      n_exponent =6.3#6.25## 5.4
      m_exponent = 0.0
      activation_energy = 0.0
     #relative_tolerance=1e-6
     # absolute_tolerance=1e-4
    [../]
#  [./stress]
#    type = ComputeFiniteStrainElasticStress
#  [../]
[]
[UserObjects]
  [cut_mesh2]
    type = MeshCut2DFractureUserObject
    mesh_file = line0p45.e
    growth_increment =0.15
    ki_vectorpostprocessor = "II_KI_1"
    k_critical =0
  []
[]

[Executioner]
  type = Transient

  solve_type = 'Newton'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  line_search = 'none'
  [./Quadrature]
    type = GAUSS
    order =SIXTH
  [../]

#  [./Predictor]
#    type = SimplePredictor
#    scale = 0.0
#  [../]

  l_max_its = 50
  nl_max_its = 500

  nl_abs_tol = 1e-4
  nl_rel_tol = 1e-7
  #l_tol=1e-3

# time control
  start_time = 0.0
  dt = 0.1
  end_time =0.5

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
