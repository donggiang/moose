# Tier 1.1 — Single-element mode-I monotonic pull
# See verification_plans/tier1_single_element.md
#
# Two stacked unit-square blocks split by BreakMeshByBlockGenerator.
# Top face pulled in +y from 0 to 1.5 * delta_n_f.
# Bulk is stiff linear elastic so the cohesive jump dominates the global displacement.
#
# Reference law parameters:
#   K  = 5e3          delta_n_0 = N/K       = 1.0e-2  (softening onset)
#   N  = 50           delta_n_f = 2 GIc/N   = 2.0e-2  (full decohesion)
#   GI_c = 0.5
#
# K was chosen so the elastic-softening transition coincides with the
# breakpoint (t=0.1, u=0.01) of the loading function below — both branches
# of the bilinear law are then well-resolved on the same dt.
#
# Loading is piecewise linear so both the elastic (K*delta) and softening
# branches are well-resolved by the same dt.

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
    # Piecewise linear: 10 % of the time goes to the elastic regime
    # (0 -> 1.5*delta_n_0 = 7.5e-5), the rest covers the full softening
    # branch out to 1.5 * delta_n_f = 0.03.
    type = PiecewiseLinear
    x = '0 0.1 1.0'
    y = '0 0.01 0.03'
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
    # No viscous regularization — we want the analytical bilinear law exactly
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
  dt = 0.01
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]
