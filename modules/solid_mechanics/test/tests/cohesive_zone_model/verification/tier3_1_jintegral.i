# Tier 3.1 — J-integral / dissipated-energy conservation
# See verification_plans/tier3_jintegral_checks.md
#
# Geometry, mesh, bulk, and loading are identical to Tier 1.1: a single-element
# mode-I monotonic pull driven past delta_n_f. The verifier checks that the
# cumulative work done by the cohesive traction on the interface equals the
# input mode-I fracture energy G_Ic at full decohesion (Rice contour identity
# in the small-strain limit).
#
# Why this catches things Tier 1.1 doesn't:
#   1.1 only verifies the pointwise t_n(delta_n) curve; it cannot detect a
#   mistaken interface area, a sign error in the traction direction during
#   integration, or a finite-strain reference/current configuration mix-up.
#   This test exercises the *integral* of t * d(delta) on the interface.
#
# Reference law parameters (same as Tier 1.1):
#   K  = 5e3   N = 50   GIc = 0.5
#   delta_n_0 = 1.0e-2   delta_n_f = 2.0e-2
#
# Interface area (per unit out-of-plane thickness in 2D):
#   A_iface = 1   (interface segment length, x in [0,1])
#
# Expected: W_diss(delta_n >= delta_n_f) = G_Ic * A_iface = 0.5

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
        generate_output = 'normal_traction tangent_traction normal_jump tangent_jump traction_y jump_y'
      []
    []
  []
[]

[Functions]
  [stretch]
    # Same piecewise-linear ramp as Tier 1.1 — well-resolved both branches.
    type = PiecewiseLinear
    x = '0    0.10  1.00'
    y = '0    0.01  0.03'
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
    youngs_modulus = 1e9
    poissons_ratio = 0.3
  []
  [stress]
    type = ADComputeLinearElasticStress
  []
  [czm]
    type = BiLinearMixedModeTraction
    boundary = 'Block1_Block2'
    penalty_stiffness = 5e3
    GI_c = 0.5
    GII_c = 0.5
    normal_strength = 50
    shear_strength = 30
    eta = 1.45
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
  [iface_length]
    # Reports the interface measure (length in 2D = "area" in the line integral).
    type = AreaPostprocessor
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

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
  l_tol = 1e-10

  start_time = 0.0
  end_time = 1.0
  dt = 0.005
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]
