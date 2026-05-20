# Tier 1.4 — Single-element mixed-mode loading at fixed mode-mix beta = 1
# See verification_plans/tier1_single_element.md
#
# Same mesh and bulk material as tier1_1. The top face is pulled simultaneously
# in +x and +y by the SAME function so that the cohesive interface sees
# delta_x = delta_y (i.e. delta_shear = delta_n, hence beta = delta_s/delta_n = 1).
#
# Reference law parameters:
#   K = 5e3   N = 50   S = 30   GI_c = GII_c = 0.5   eta = 1.45
#   delta_n_0 = N/K = 1.0e-2     delta_s_0 = S/K = 6.0e-3
#
# Camanho-Davila mixed-mode landmarks at beta = 1 (with S = T):
#   delta_m_0 = delta_n_0 * delta_s_0 * sqrt((1+beta^2)/(delta_s_0^2 + (beta*delta_n_0)^2))
#             = 7.27607e-3
#   delta_m_f = (2/(K*delta_m_0)) * [GI_c + (GII_c-GI_c)*(beta^2/(1+beta^2))^eta]
#             = 2.74876e-2     (since GI_c == GII_c the BK correction vanishes)
#   peak traction in mixed mode  : K * delta_m_0 = 36.380
#
# Per-axis loading endpoints (delta_m = sqrt(2) * u with bulk much stiffer than the interface):
#   u_0  = delta_m_0 / sqrt(2)  = 5.1452e-3
#   u_f  = delta_m_f / sqrt(2)  = 1.9437e-2
#   u_max = 1.5 * u_f           = 2.9155e-2

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
        generate_output = 'normal_traction tangent_traction normal_jump tangent_jump traction_x traction_y jump_x jump_y'
      []
    []
  []
[]

[Functions]
  [pull]
    # Breakpoint at the per-axis mixed-mode onset displacement
    # u_0 = delta_m_0 / sqrt(2) = 5.1452e-3.
    type = PiecewiseLinear
    x = '0    0.1        1.0'
    y = '0    5.1452e-3  2.9155e-2'
  []
[]

[BCs]
  # Bottom face fully clamped.
  [fix_x_bottom]
    type = DirichletBC
    variable = disp_x
    boundary = bottom
    value = 0
    preset = true
  []
  [fix_y_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
    preset = true
  []
  # Top face pulled by the SAME function in both directions to enforce beta = 1.
  [pull_x_top]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = top
    function = pull
    preset = true
  []
  [pull_y_top]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = pull
    preset = true
  []
[]

[Materials]
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    # Bumped from 1e9 to 1e13 so bulk compliance is ~1e-9 of interface compliance
    # in both normal and shear, reducing the residual normal/shear asymmetry that
    # would otherwise break beta = 1 by O(1e-5). Verifier still measures beta.
    youngs_modulus = 1e13
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
  [jump_x]
    type = SideAverageValue
    variable = jump_x
    boundary = 'Block1_Block2'
  []
  [jump_y]
    type = SideAverageValue
    variable = jump_y
    boundary = 'Block1_Block2'
  []
  [traction_x]
    type = SideAverageValue
    variable = traction_x
    boundary = 'Block1_Block2'
  []
  [traction_y]
    type = SideAverageValue
    variable = traction_y
    boundary = 'Block1_Block2'
  []
  [applied_disp]
    type = FunctionValuePostprocessor
    function = pull
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
