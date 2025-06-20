[GlobalParams]
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = true
[]
[XFEM]
  geometric_cut_userobjects = 'cut_mesh2'
  qrule = volfrac
  output_cut_plane = true
  use_crack_tip_enrichment =true
  crack_front_definition = crack_tip
  enrichment_displacements = 'enrich1_x enrich2_x enrich3_x enrich4_x enrich1_y enrich2_y enrich3_y enrich4_y'
  displacements = 'disp_x disp_y'
  cut_off_boundary = all
  cut_off_radius = 0.25
[]
[Mesh]
   file = geo.e
[]



[DomainIntegral]
  integrals = 'Jintegral InteractionIntegralKI InteractionIntegralKII'
  displacements = 'disp_x disp_y'
  crack_front_points_provider = cut_mesh2
  2d = true
  number_points_from_provider = 2
  crack_direction_method = CurvedCrackFront
  radius_inner = '0.1 0.3'
  radius_outer = '0.3 0.5'
  poissons_ratio = 0.3
  youngs_modulus = 207000
  used_by_xfem_to_grow_crack = false
  output_q = false
  output_vpp = false
  incremental = true
[]

[UserObjects]
  [cut_mesh2]
    type = MeshCut2DFractureUserObject
    mesh_file = geo1d.e
    growth_increment = 0.0032
    ki_vectorpostprocessor = "II_KI_1"
    kii_vectorpostprocessor = "II_KII_1"
    k_critical = 0
  []
  [./crack_tip]
    type = CrackFrontDefinition
    crack_direction_method = CrackDirectionVector
    crack_front_points = '0.25 1.0 0'
    crack_direction_vector = '1 0 0'
    2d = true
    axis_2d = 2
  [../]
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
    displacements = 'disp_x disp_y'
        add_variables = true
    incremental = true
    generate_output = 'stress_xx stress_yy stress_zz'
  [../]
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
  [./pull]
    type = PiecewiseLinear
    x='0  50   100'
    y='0  0.02 0.1'
  [../]
  [dt_func]
    type = PiecewiseLinear
    x = '0.    1e-12 1e20'
    y = '1e-12 0.1 0.1'
  [../]
[]

[BCs]
  [./fix_y]
    type = DirichletBC
    boundary = bottom
    variable = disp_y
    value = 0.0
  [../]
  [./fix_x]
    type = DirichletBC
    boundary = bottom
    variable = disp_x
    value = 0.0
  [../]
  [./topx]
    type = DirichletBC
    boundary = top
    variable = disp_x
    value = 0.0
  [../]
  [./topy]
    type = FunctionDirichletBC
    boundary = top
    variable = disp_y
    function = pull
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 207000
    poissons_ratio = 0.3
  [../]
  [./strain]
   # type = ComputeIncrementalStrain
    type = ComputeCrackTipEnrichmentIncrementalStrain
    displacements = 'disp_x disp_y'
    crack_front_definition = crack_tip
    enrichment_displacements = 'enrich1_x enrich2_x enrich3_x enrich4_x enrich1_y enrich2_y enrich3_y enrich4_y'
  [../]
  [./stress]
    type = ComputeFiniteStrainElasticStress
    #type = ComputeStrainIncrementBasedStress
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
   petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  line_search = 'none'

 # solve_type = NEWTON
 # petsc_options_iname = '-ksp_type -pc_type'
 # petsc_options_value = 'preonly lu'

  # Since we do not sub-triangularize the tip element,
  # we need to use higher order quadrature rule to improve
  # integration accuracy.
  # Here second = SECOND is for regression test only.
  # However, order = SIXTH is recommended.
  [./Quadrature]
    type = GAUSS
    order = SIXTH
  [../]

  [./Predictor]
    type = SimplePredictor
    scale = 1.0
  [../]

  # controls for linear iterations
  l_max_its = 100
  l_tol = 1e-8

  # controls for nonlinear iterations
  nl_max_its = 16
  nl_rel_tol = 1e-4 #11
  nl_abs_tol = 1e-4#12

  # time control
  start_time = 0.0
  dt = 0.125
  end_time = 150
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]


[Outputs]
  file_base = edge_crack_2d_enriched_elas_out
  exodus = true
  csv = true
  [./console]
    type = Console
    output_linear = true
  [../]
[]
