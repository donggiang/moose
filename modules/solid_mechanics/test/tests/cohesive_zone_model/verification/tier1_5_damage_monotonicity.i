# Tier 1.5 — Damage monotonicity stress test
# See verification_plans/tier1_single_element.md
#
# Same mesh / material / solver as tier1_3, but a five-segment loading
# function that intentionally cycles the jump down and up several times,
# advancing the running maximum jump (and therefore damage) each loading
# pass:
#
#   t = [0.00, 0.10] : 0.000  -> 0.012   load past delta_0 (d climbs to ~0.33)
#   t = [0.10, 0.20] : 0.012  -> 0.006   partial unload (d frozen)
#   t = [0.20, 0.40] : 0.006  -> 0.018   reload past previous max (d -> ~0.89)
#   t = [0.40, 0.60] : 0.018  -> 0.000   full unload (d frozen at 0.89)
#   t = [0.60, 1.00] : 0.000  -> 0.030   reload past delta_f (d -> 1)
#
# Verifier checks:
#   * traction at every step equals (1 - d(delta_max)) * K * delta_n
#   * d(delta_max) is non-decreasing (irreversibility holds)
#   * final cumulative dissipated energy equals G_Ic
#
# Reference law parameters (same as tier1_1):
#   K  = 5e3   N = 50   GIc = 0.5
#   delta_n_0 = 1.0e-2   delta_n_f = 2.0e-2

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
    type = PiecewiseLinear
    x = '0    0.10   0.20   0.40   0.60   1.00'
    y = '0    0.012  0.006  0.018  0.000  0.030'
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
