# Tier 2.1 -- DCB CG-CZM, MOOSE-generated 2-block mesh (upper / lower),
# displacement control with OUTER point-pin (y = +/-1.98).
#
# 4-block experiment failed because BreakMeshByBlock leaves the crack-tip
# corner (95, 0) shared by ALL FOUR blocks, tying upper and lower arms at
# the tip and producing rigid-bond behavior. The Cubit wedge mesh sidesteps
# this with a physical wedge gap.
#
# 2-block fix: split y in [-1.98, 1.98] into 'upper' (y>0) and 'lower' (y<0)
# only. BreakMeshByBlock then splits the WHOLE y=0 interface (x in [0,150])
# cleanly. Cohesive law applied only on the bonded portion (x <= 95) via a
# sub-sideset; the precrack portion (x > 95) is left as a free surface.
#
# Material identical to the wedge variant.

[Mesh]
  [base]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 150
    ymin = -1.98
    ymax =  1.98
    nx = 300
    ny = 8
  []
  [upper]
    type = SubdomainBoundingBoxGenerator
    input = base
    block_id = 1
    block_name = upper
    bottom_left = '0   0    0'
    top_right   = '150 1.98 0'
  []
  [lower]
    type = SubdomainBoundingBoxGenerator
    input = upper
    block_id = 2
    block_name = lower
    bottom_left = '0   -1.98 0'
    top_right   = '150  0    0'
  []
  [break]
    type = BreakMeshByBlockGenerator
    input = lower
    block_pairs = 'upper lower'
    split_interface = true
  []
  # Restrict CZM to the bonded portion (x <= 95) using a sub-sideset
  # of the upper-side split interface 'upper_lower'.
  [bonded_iface]
    type = ParsedGenerateSideset
    input = break
    new_sideset_name = bonded_iface
    combinatorial_geometry = 'abs(y) < 1e-6 & x < 95.0 + 1e-6'
    included_subdomains = 'upper'
  []
  # Sidesets for clamping the left edge of each arm.
  [left_top]
    type = ParsedGenerateSideset
    input = bonded_iface
    new_sideset_name = left_top
    combinatorial_geometry = 'abs(x) < 1e-6 & y > 0.001'
  []
  [left_bot]
    type = ParsedGenerateSideset
    input = left_top
    new_sideset_name = left_bot
    combinatorial_geometry = 'abs(x) < 1e-6 & y < -0.001'
  []
  [pull_upper]
    type = ExtraNodesetGenerator
    input = left_bot
    new_boundary = pull_upper
    coord = '150  1.98 0'
  []
  [pull_lower]
    type = ExtraNodesetGenerator
    input = pull_upper
    new_boundary = pull_lower
    coord = '150 -1.98 0'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
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
  [fy]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    save_in = fy
  []
[]

[Physics]
  [SolidMechanics]
    [CohesiveZone]
      [czm]
        strain = SMALL
        boundary = 'bonded_iface'
        generate_output = 'normal_traction tangent_traction normal_jump tangent_jump'
      []
    []
  []
[]

[Functions]
  [stretch]
    type = PiecewiseLinear
    x = '0 1.0'
    y = '0 8.0'
  []
  [stretch_neg]
    type = PiecewiseLinear
    x = '0 1.0'
    y = '0 -8.0'
  []
[]

[BCs]
  [left_top_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left_top'
    value = 0
    preset = true
  []
  [left_top_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'left_top'
    value = 0
    preset = true
  []
  [left_bot_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left_bot'
    value = 0
    preset = true
  []
  [left_bot_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'left_bot'
    value = 0
    preset = true
  []
  [pull_upper_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = pull_upper
    function = stretch
    preset = true
  []
  [pull_lower_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = pull_lower
    function = stretch_neg
    preset = true
  []
[]

[Materials]
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 150e3
    poissons_ratio = 0.25
  []
  [strain]
    type = ADComputeSmallStrain
  []
  [stress]
    type = ADComputeLinearElasticStress
  []
  [czm]
    type = BiLinearMixedModeTraction
    boundary = 'bonded_iface'
    penalty_stiffness = 1e6
    GI_c   = 0.268
    GII_c  = 1.45
    normal_strength = 30.0
    shear_strength  = 40.0
    eta = 2.284
    displacements = 'disp_x disp_y'
    viscosity = 0
  []
[]

[Postprocessors]
  [P_top_per_width]
    type = NodalSum
    variable = fy
    boundary = pull_upper
  []
  [P_bot_per_width]
    type = NodalSum
    variable = fy
    boundary = pull_lower
  []
  [P_top_N]
    type = ParsedPostprocessor
    expression = 'abs(P_top_per_width) * 20.0'
    pp_names = 'P_top_per_width'
  []
  [P_bot_N]
    type = ParsedPostprocessor
    expression = 'abs(P_bot_per_width) * 20.0'
    pp_names = 'P_bot_per_width'
  []
  [tip_disp_top]
    type = NodalExtremeValue
    variable = disp_y
    boundary = pull_upper
  []
  [tip_disp_bot]
    type = NodalExtremeValue
    variable = disp_y
    boundary = pull_lower
    value_type = min
  []
  [opening]
    type = ParsedPostprocessor
    expression = 'tip_disp_top - tip_disp_bot'
    pp_names   = 'tip_disp_top tip_disp_bot'
  []
  [max_normal_traction]
    type = SideExtremeValue
    variable = normal_traction
    boundary = 'bonded_iface'
  []
  [max_normal_jump]
    type = SideExtremeValue
    variable = normal_jump
    boundary = 'bonded_iface'
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

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  nl_max_its = 30
  l_tol = 1e-10

  start_time = 0.0
  end_time   = 1.0

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.25
    cutback_factor = 0.5
    cutback_factor_at_failure = 0.5
  []
  dtmin = 1e-6
  dtmax = 0.02
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]
