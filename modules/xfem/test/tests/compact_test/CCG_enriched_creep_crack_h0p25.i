[Problem]

  type = ReferenceResidualProblem
  reference_vector = 'ref'
  extra_tag_vectors = 'ref'

[]
[GlobalParams]
  displacements = 'disp_x disp_y'
  volumetric_locking_correction =false
  absolute_value_vector_tags = ref

[]

[XFEM]
###  geometric_cut_userobjects = 'cut_mesh2'
  #qrule = volfrac
  qrule = moment_fitting
  output_cut_plane = true
  use_AD=false
  use_crack_tip_enrichment =true
  crack_front_definition = crackFrontDefinition
  enrichment_displacements = 'enrich1_x enrich2_x enrich3_x enrich4_x enrich1_y enrich2_y enrich3_y enrich4_y'
  displacements = 'disp_x disp_y'
  cut_off_boundary = all
  cut_off_radius =0.2693  ####0.527 #0.1752 #0.2828 #0.10984 # 0.175219
[]


####1/11 for h=0.1, rd=0.10975



[Mesh]
  file = ct_msh_h0p25_ref_v1.e #compact_test2d_h05_4.e
[]

[DomainIntegral]
  integrals = 'CIntegral InteractionIntegralKI InteractionIntegralKII'
  displacements = 'disp_x disp_y'
  crack_front_points_provider = cut_mesh2
  2d = true
  number_points_from_provider = 1
  crack_direction_method = CurvedCrackFront
  radius_inner = '2.5' #'1.4'
  radius_outer = '6.5' #'4.0'
  youngs_modulus = 120000.0
  poissons_ratio = 0.3
 inelastic_models = 'powerlawcrp'
  output_q = true
  output_vpp = false
  incremental = true
  used_by_xfem_to_grow_crack = false
  use_automatic_differentiation = false
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
   #use_displaced_mesh = false
    displacements = 'disp_x disp_y'
    add_variables = true
    incremental = true
    generate_output = 'stress_xx stress_yy'
    use_automatic_differentiation =false
  []
[]

[AuxKernels]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
    execute_on = timestep_end
  [../]
  [./vonmises]
    type = RankTwoScalarAux
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
    y = '0. 1.0 1.0'
    scale_factor = -250.0
  [../]
 [dt_func]
    type = PiecewiseLinear
    x = '0.    1e20'
    y = '1.6 1.6'
  [../]
[]


[DiracKernels]
  [point1]
    type = ConstantPointSource
    variable = disp_y
    point = '5.0 7.1'
    value = 2050.0 # P = 10
  []
  [point2]
    type = ConstantPointSource
    variable = disp_y
    point = '5.0 -7.1'
    value = -2050.0 # P = 10
  []
[]

[BCs]
  [ydisp2]
    type = DirichletBC
    variable = disp_y
    boundary = 'middle_right'
    value = 0.0
  []
  [ydisp3]
    type = DirichletBC
    variable = disp_y
    boundary = 'middle_right'
    value = 0.0
  []
  [xdisp2]
    type = DirichletBC
    variable = disp_x
    boundary = 'top_load'
    value = 0.0
  []
  [xdisp3]
    type = DirichletBC
    variable = disp_x
    boundary = 'bottom_load'
    value = 0.0
  []
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
  youngs_modulus = 120000
  poissons_ratio = 0.3
  [../]
[]

[Materials]
  [./strain]
    type =ComputeCrackTipEnrichmentIncrementalStrain #ADComputeGreenLagrangeStrain
    displacements = 'disp_x disp_y'
    crack_front_definition = crackFrontDefinition
    enrichment_displacements = 'enrich1_x enrich2_x enrich3_x enrich4_x enrich1_y enrich2_y enrich3_y enrich4_y'
  [../]
  #[./stress]
  #  type = ADComputeFiniteStrainElasticStress
  #[../]
  [./radial_return_stress]
    type = ComputeMultipleInelasticStress
    inelastic_models ='powerlawcrp'
  [../]
  [./powerlawcrp]
    type = PowerLawCreepStressUpdate
    coefficient =1.0e-23#2e-23 #
    n_exponent =7.1#6.25## 5.4
    m_exponent = 0.0
    activation_energy = 0.0
   #relative_tolerance=1e-6
   # absolute_tolerance=1e-4
  [../]
[]

[UserObjects]
  [cut_mesh2]
    type = MeshCut2DCCGUserObject
    mesh_file = initialcrack14p86.e
    growth_increment = 0.000
    ki_vectorpostprocessor = "II_KI_1"
    kii_vectorpostprocessor = "II_KII_1"
    c_vectorpostprocessor="C_1"
    paris_coeff =3e-4 #1.5e-3#
    paris_exponent =0.9#### 1.03#1.25
  []
[]

[Executioner]
  type = Transient

  solve_type = 'Newton'


  #petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type '
  petsc_options_value = 'lu '
  [./Quadrature]
    type = GAUSS
    order = SIXTH
  [../]


  l_max_its = 50
  nl_max_its = 300

  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-6
 #line_search=none
# time control
  start_time = 0.0
 # dt = 1.6
  end_time =500
  automatic_scaling = true
  [TimeStepper]
    type = FunctionDT
    function = '0.8'
    growth_factor = 2
    cutback_factor_at_failure = 0.5
  []
  max_xfem_update = 1
[]



[Outputs] ###dt3_0p225.

  file_base = CCG_en_cr_cl_h0p25_d4_m5em4_n0p9
  exodus = true
  csv = true
  execute_on = timestep_end
  [./console]
    type = Console
    output_linear = true
  [../]
[]
