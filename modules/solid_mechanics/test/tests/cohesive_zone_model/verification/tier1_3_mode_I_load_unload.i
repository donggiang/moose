# Tier 1.3 — Single-element mode-I load / unload / reload
# See verification_plans/tier1_single_element.md
#
# Same mesh, materials, and solver as tier1_1, but the loading function has
# four phases:
#
#   t = [0.00, 0.10]: 0      -> 0.010   elastic loading on the bilinear envelope
#   t = [0.10, 0.20]: 0.010  -> 0.015   continued loading on the softening branch
#   t = [0.20, 0.40]: 0.015  -> 0.000   unload (secant of slope K*(1-d_max))
#   t = [0.40, 0.60]: 0.000  -> 0.015   reload along the same secant
#   t = [0.60, 1.00]: 0.015  -> 0.030   continue softening past delta_f
#
# Reference law parameters (same as tier1_1):
#   K  = 5e3   N = 50   GIc = 0.5
#   delta_n_0 = N/K       = 1.0e-2
#   delta_n_f = 2 GIc/N   = 2.0e-2
#   delta_n_max held at 0.015 (midpoint of softening branch)
#
# Damage at delta_n_max = 0.015:
#   d_max = delta_f * (delta_max - delta_0) / (delta_max * (delta_f - delta_0))
#         = 0.02 * 0.005 / (0.015 * 0.01) = 2/3
#   secant slope = K (1 - d_max) = 5e3 * 1/3 = 1666.67
#   traction at delta_max:  N * (delta_f - delta_max)/(delta_f - delta_0) = 25

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
    x = '0    0.10  0.20   0.40   0.60   1.00'
    y = '0    0.01  0.015  0.0    0.015  0.03'
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
