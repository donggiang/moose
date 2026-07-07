nl = 2
h=1
ri = 0
ro = 4
[GlobalParams]
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = false
[]

[XFEM]
  geometric_cut_userobjects = 'cut_mesh2'
  qrule = vol_frac #moment_fitting #
  output_cut_plane = true
  use_crack_tip_enrichment = true
  crack_front_definition = crackFrontDefinition
  enrichment_displacements = 'enrich1_x enrich2_x enrich3_x enrich4_x enrich1_y enrich2_y enrich3_y enrich4_y'
  cut_off_boundary = enriched_interface
  use_AD=true
  block=2
[]

[Mesh]
  [cutter_mesh]
    type = PolyLineMeshGenerator
    points = '-40.01 0 0
              -17.5 0 0'
    loop = false
    num_edges_between_points = 1
    save_with_name = cut_mesh
  []
  [gen]
    type = FileMeshGenerator
    file = msh_ct_h2.e
  []
  [main_domain]
    type = SubdomainBoundingBoxGenerator
    input = gen
    block_id = 1
    bottom_left = '-41 -20 0'
    top_right = '0.1 20 0'
  []
  [enriched_domain]
    type = SubdomainBoundingBoxGenerator
    input = main_domain
    block_id = 2
    bottom_left = '-29.1 -3 0'
    top_right = '-10 3 0'
  []
  [nset_enriched]
    type = SideSetsAroundSubdomainGenerator
    input = enriched_domain
    block = '2'
    new_boundary = 'enriched_interface'
  []
  [upper_node]
    type = ExtraNodesetGenerator
    coord = '-32 8.8 0'
    input = nset_enriched
    new_boundary = 'upper_node'
    use_closest_node = true
  []
  [lower_node]
    type = ExtraNodesetGenerator
    coord = '-32 -8.8 0'
    input = upper_node
    new_boundary = 'lower_node'
    use_closest_node = true
  []
  final_generator = 'lower_node'
[]

[DomainIntegral]
  integrals = 'CIntegral InteractionIntegralKI InteractionIntegralKII'
  crack_front_points_provider = cut_mesh2
  2d = true
  number_points_from_provider = 1
  crack_direction_method = CurvedCrackFront
  radius_inner = '${fparse h*ri}'
  radius_outer = '${fparse h*ro}'
  youngs_modulus = 132000.0
  poissons_ratio = 0.3
  inelastic_models = 'powerlawcrp'
  output_q = true
  output_vpp = false
  incremental = true
  use_automatic_differentiation =true
[]

[UserObjects]
  [cut_mesh2]
    type = MeshCut2DCCGUserObject
    mesh_generator_name = 'cut_mesh'
    growth_increment = 0.000
    ki_vectorpostprocessor = "II_KI_1"
    kii_vectorpostprocessor = "II_KII_1"
    c_vectorpostprocessor="C_1"
    paris_coeff =0.032#1.5e-3#
    paris_exponent =0.1#### 1.03#1.25
  []
  [esm]
   type = CrackTipNodeLayerSubdomainModifier
   crack_front_definition = crackFrontDefinition
   enriched_subdomain_id =2
   base_subdomain_id =1# 0.7
   num_node_layers = ${nl}
  []
  [side_updater]
    type = SidesetAroundSubdomainUpdater
    inner_subdomains = 2
    outer_subdomains = 1
    update_boundary_name = enriched_interface
  []
[]

[Variables]
  [disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [saved_x]
  []
  [saved_y]
  []
  [stress_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_xy]
    order = CONSTANT
    family = MONOMIAL
  []
  [vonmises]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [TensorMechanics]
    use_displaced_mesh = false
    displacements = 'disp_x disp_y'
    add_variables = true
    incremental = true
    generate_output = 'stress_xx stress_yy'
    use_automatic_differentiation =true
  []
[]

[AuxKernels]
  [stress_xx]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  []
  [stress_yy]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  []
  [stress_xy]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
    execute_on = timestep_end
  []
  [vonmises]
    type = ADRankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    scalar_type = vonmisesStress
    execute_on = timestep_end
  []
[]


[BCs]
  [right_x]
    type = DirichletBC
    boundary = 'right upper_node lower_node'
    variable = disp_x
    value = 0
  []
  [right_y]
    type = DirichletBC
    boundary = 'right'
    variable = disp_y
    value = 0
  []
[]

[DiracKernels]
  [point1]
    type = ConstantPointSource
    variable = disp_y
    point = '-32 8.8 0'
    value = 9 # P = 10
  []
  [point2]
    type = ConstantPointSource
    variable = disp_y
    point = '-32 -8.8 0'
    value = -9 # P = 10
  []
[]


[Materials]
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 132000.0
    poissons_ratio = 0.3
  []
  [enrich_strain]
    type = ADComputeCrackTipEnrichmentIncrementalStrain
    crack_front_definition = crackFrontDefinition
    enrichment_displacements = 'enrich1_x enrich2_x enrich3_x enrich4_x enrich1_y enrich2_y enrich3_y enrich4_y'
    block=2
  []
  [strain]
    type = ADComputeIncrementalStrain
    block=1
  []
  [./radial_return_stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'powerlawcrp'
  [../]
  [./powerlawcrp]
    type = ADPowerLawCreepStressUpdate
    coefficient =1.0e-8#2e-23 #
    n_exponent =9#6.25## 5.4
    m_exponent = 0.0
    activation_energy = 0.0
   #relative_tolerance=1e-6
   # absolute_tolerance=1e-4
  [../]
[]


[Executioner]
  type = Transient

  solve_type = 'Newton'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  #automatic_scaling = true
  [Quadrature]
    type = GAUSS
    order = SIXTH
  []

  #l_max_its = 50
  nl_max_its = 40

  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-9
  start_time = 0.0
  #dt = 1
  [TimeStepper]
    type = FunctionDT
    function = '10'
    growth_factor =2
    cutback_factor_at_failure = 0.99
  []
  end_time = 30
  max_xfem_update = 1
[]

[Outputs]
  [xfemcutter]
    type = XFEMCutMeshOutput
    xfem_cutter_uo = cut_mesh2
  []
  file_base =ct_h2_c1em10_e9_dt8_m0p045_n0p1
  exodus = true
  csv=true
  [./console]
    type = Console
    output_linear = true
  [../]
[]
