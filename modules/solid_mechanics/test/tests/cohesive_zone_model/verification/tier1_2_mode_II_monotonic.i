# Tier 1.2 — Single-element mode-II monotonic shear
# See verification_plans/tier1_single_element.md
#
# Same mesh and bulk material as tier1_1. The top face is sheared in +x
# while disp_y is fixed on both top and bottom faces. Pure isotropic
# linear elasticity in plane strain produces u_y = 0 everywhere, so the
# interface sees pure mode-II loading: jump_y = 0, jump_x = u_*(t).
#
# Reference law parameters (mode II):
#   K  = 5e3   S = 30   GII_c = 0.5
#   delta_s_0 = S/K        = 6.0e-3   (shear softening onset)
#   delta_s_f = 2 GIIc/S   = 3.333e-2 (full decohesion)
#
# Loading function breakpoint at (t=0.1, u=delta_s_0) so the slope change
# coincides with the elastic-softening transition.

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
  [shear]
    type = PiecewiseLinear
    x = '0 0.1 1.0'
    y = '0 0.006 0.05'
  []
[]

[BCs]
  # Bottom face fully clamped.
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
  # Top face: shear in +x, no opening in y.
  [fix_y_top]
    type = DirichletBC
    variable = disp_y
    boundary = top
    value = 0
    preset = true
  []
  [shear_x_top]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = top
    function = shear
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
    function = shear
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
