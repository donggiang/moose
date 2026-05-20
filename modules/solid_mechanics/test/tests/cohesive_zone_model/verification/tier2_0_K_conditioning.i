# Tier 2.0 — K-conditioning study, single-element mode-I pull (AS4/PEEK)
# See verification_plans/tier2_beam_benchmarks.md
#
# Same single-element setup as Tier 1.1 but with Turon's AS4/PEEK material
# parameters. Sweep `Materials/czm/penalty_stiffness` (K) via cli_args:
#
#     N        = 80    MPa     normal_strength (sigma_3^0)
#     S        = 100   MPa     shear_strength
#     GI_c     = 0.969 N/mm    mode-I fracture energy
#     GII_c    = 1.719 N/mm    mode-II fracture energy
#     eta (BK) = 2.284
#     E (bulk) = 122.7 GPa = 122.7e3 MPa  (AS4/PEEK E1, used isotropically)
#
# Validity constraint on K:
#     delta_n_0 < delta_n_f
#     N/K < 2 GI_c / N
#     K   > N^2 / (2 GI_c) = 80^2 / 1.938 = 3303
#
# Default K = 1e5 for the standalone run; sweep via:
#     -i ... Materials/czm/penalty_stiffness=<K>
#            Outputs/file_base=tier2_0_K_<K>
#
# Loading function spans 0 -> 1.5 * delta_n_f = 0.0363 mm (independent of K)
# with breakpoints chosen so the elastic regime is well-sampled across
# K = 1e4 -> 1e6 (delta_n_0 = N/K spans 8e-3 down to 8e-5).

[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 2
    xmax = 1
    ymax = 2
    nx = 1
    ny = 2
  []
  [block1]
    type = SubdomainBoundingBoxGenerator
    input = msh
    bottom_left = '0 0 0'
    top_right = '1 1 0'
    block_id = 1
    block_name = Block1
  []
  [block2]
    type = SubdomainBoundingBoxGenerator
    input = block1
    bottom_left = '0 1 0'
    top_right = '1 2 0'
    block_id = 2
    block_name = Block2
  []
  [split]
    type = BreakMeshByBlockGenerator
    input = block2
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        add_variables = true
        use_automatic_differentiation = true
      []
    []
    [CohesiveZone]
      [czm]
        strain = SMALL
        boundary = 'Block1_Block2'
        generate_output = 'normal_traction tangent_traction normal_jump tangent_jump'
      []
    []
  []
[]

[Functions]
  [stretch]
    # Three segments to keep elastic regime resolved across all K in the sweep.
    #   t in [0.00, 0.10]:  0     -> 1.0e-4   (~ 0 -> delta_n_0 for K=8e5)
    #   t in [0.10, 0.40]:  1e-4  -> 8.0e-3   (~ delta_n_0 region for K=1e5)
    #   t in [0.40, 1.00]:  8e-3  -> 3.63e-2  (softening + decohesion for all K)
    type = PiecewiseLinear
    x = '0    0.10   0.40   1.00'
    y = '0    1e-4   8e-3   3.63e-2'
  []
[]

[BCs]
  [fix_y_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
    preset = true
  []
  [fix_x_bottom]
    type = DirichletBC
    variable = disp_x
    boundary = bottom
    value = 0
    preset = true
  []
  [fix_x_top]
    type = DirichletBC
    variable = disp_x
    boundary = top
    value = 0
    preset = true
  []
  [pull_y_top]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = stretch
    preset = true
  []
[]

[Materials]
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 122.7e3   # AS4/PEEK E1, MPa
    poissons_ratio = 0.25
  []
  [stress]
    type = ADComputeLinearElasticStress
  []
  [czm]
    type = BiLinearMixedModeTraction
    boundary = 'Block1_Block2'
    penalty_stiffness = 1e5    # default; override via cli_args for the sweep
    GI_c  = 0.969
    GII_c = 1.719
    normal_strength = 80
    shear_strength  = 100
    eta = 2.284
    displacements = 'disp_x disp_y'
    viscosity = 0
  []
[]

[Postprocessors]
  [normal_jump]
    type = SideAverageValue
    variable = normal_jump
    boundary = 'Block1_Block2'
  []
  [normal_traction]
    type = SideAverageValue
    variable = normal_traction
    boundary = 'Block1_Block2'
  []
  [tangent_jump]
    type = SideAverageValue
    variable = tangent_jump
    boundary = 'Block1_Block2'
  []
  [tangent_traction]
    type = SideAverageValue
    variable = tangent_traction
    boundary = 'Block1_Block2'
  []
  [applied_disp]
    type = FunctionValuePostprocessor
    function = stretch
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = none
  automatic_scaling = true

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  nl_max_its = 25
  l_tol = 1e-10

  start_time = 0.0
  end_time   = 1.0
  dt         = 0.005
[]

[Outputs]
  exodus = false
  csv = true
  print_linear_residuals = false
[]
