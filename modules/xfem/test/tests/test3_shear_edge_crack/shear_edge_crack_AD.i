L = 16.0
W = 7
a = 3.5
nx=81
ny=161
ri = 24
ro = 30
nl = 10
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
    points = '-0.01 0 0
              ${a} 0 0'
    loop = false
    num_edges_between_points = 1
    save_with_name = cut_mesh
  []
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx =${nx} # to reproduce literature use 240
    ny = ${ny} # to reproduce literature use 81
    xmin = 0
    xmax = ${W}
    ymin = '-${fparse L/2}'
    ymax = '${fparse L/2}'
    elem_type = QUAD4
  []
  [main_domain]
    type = SubdomainBoundingBoxGenerator
    input = gen
    block_id = 1
    bottom_left = '-0.1 -10.1 0'
    top_right = '10.1 10.1 0'
  []
  [enriched_domain]
    type = SubdomainBoundingBoxGenerator
    input = main_domain
    block_id = 2
    bottom_left = '3 -5 0'
    top_right = '7 5 0'
  []
  [nset_enriched]
    type = SideSetsAroundSubdomainGenerator
    input = enriched_domain
    block = '2'
    new_boundary = 'enriched_interface'
  []
  final_generator = 'nset_enriched'
[]

[DomainIntegral]
  integrals = 'JIntegral InteractionIntegralKI InteractionIntegralKII'
  crack_front_points_provider = cut_mesh2
  2d = true
  number_points_from_provider = 1
  crack_direction_method = CurvedCrackFront
  radius_inner = '${fparse W/nx*ri}'
  radius_outer = '${fparse W/nx*ro}'
  youngs_modulus = 100
  poissons_ratio = 0.3
  output_q = true
  output_vpp = false
  incremental = true
  use_automatic_differentiation = true
  enrichment_displacements = 'enrich1_x enrich2_x enrich3_x enrich4_x enrich1_y enrich2_y enrich3_y enrich4_y'
  enriched_subdomain_id = 2
[]

[UserObjects]
  [cut_mesh2]
    type = MeshCut2DFractureUserObject
    #mesh_file = lined5p71_l0p15.e
    #mesh_file = horline3p94.e
    mesh_generator_name = 'cut_mesh'
    growth_increment = 0.0
    ki_vectorpostprocessor = "II_KI_1"
    kii_vectorpostprocessor = "II_KII_1"
    k_critical = 1e14
   # block=2
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
    boundary = 'bottom'
    variable = disp_x
    value = 0
  []
  [right_y]
    type = DirichletBC
    boundary = 'bottom'
    variable = disp_y
    value = 0
  []
  [bottom_left_elem_y]
    type = NeumannBC
    boundary = 'top'
    variable = disp_x
    value = 1
  []
[]


[Materials]
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 100 #3e7
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
  [stress]
    type = ADComputeFiniteStrainElasticStress
  []
[]


[Executioner]
  type = Transient

  solve_type = 'Newton'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search=none
  automatic_scaling = true
  [Quadrature]
    type = GAUSS
    order = SIXTH
  []

  #l_max_its = 50
  nl_max_its = 800

  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-8
  start_time = 0.0
  dt = 0.1
  end_time = 0.1
  max_xfem_update = 0
[]

[Outputs]
  #[xfemcutter]
  #  type = XFEMCutMeshOutput
  #  xfem_cutter_uo = cut_mesh2
  #[]
  file_base ='shear_edge_crack_nx${nx}_ny${ny}_ri${ri}_ro${ro}_nl${nl}'
  exodus = true
  csv=true
  [./console]
    type = Console
    output_linear = true
  [../]
[]
