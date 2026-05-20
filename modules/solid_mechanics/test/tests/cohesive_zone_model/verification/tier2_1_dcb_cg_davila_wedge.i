# Tier 2.1 -- DCB CG-CZM under DISPLACEMENT control with POINT-load pin
# Uses Cubit mesh msh_dcb_h0p5_2blocks_rf3.e: 2 blocks ('upper' / 'lower')
# with a PHYSICAL wedge gap representing the precrack (x in [95, 150]).
# The bonded zone (y=0 for x in [0, 95]) shares LAGRANGE nodes between
# the two blocks and is split by BreakMeshByBlockGenerator to host the
# cohesive interface.
#
# Material / loading match dg_notes/dg_phase5_dcb_davila_dgczm.i:
#   E = 150 GPa, nu = 0.25
#   CZM: K = 1e6, GIc = 0.268, GIIc = 1.45, T = 30, S = 40, eta_BK = 2.284
#   Pin pulls top arm at (150, +0.9906), bottom at (150, -0.9906)
#   (closest mesh nodes to y = +/-0.99 used previously).
#
# Reaction: NodalSum of the saved y-momentum kernel residual at each pin.

[Mesh]
  [file]
    type = FileMeshGenerator
    file = dg_notes/msh_dcb_h0p5_2blocks_rf3.e
  []
  # Break the bonded portion of the upper-lower interface to insert the
  # cohesive zone. The precrack region is already a wedge gap (no shared
  # nodes), so nothing happens there.
  [break]
    type = BreakMeshByBlockGenerator
    input = file
    block_pairs = 'upper lower'
    split_interface = true
  []
  # Left-edge sidesets (one per arm, excluding y=0 duplicates from the split)
  [left_top]
    type = ParsedGenerateSideset
    input = break
    new_sideset_name = left_top
    combinatorial_geometry = 'abs(x) < 1e-6 & y > 0.001'
  []
  [left_bot]
    type = ParsedGenerateSideset
    input = left_top
    new_sideset_name = left_bot
    combinatorial_geometry = 'abs(x) < 1e-6 & y < -0.001'
  []
  # Single-NODE pin nodesets at the OUTER right-edge corners of each arm
  # (Davila-style loading pin on the back face, not the arm centroid).
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
        boundary = 'upper_lower'
        generate_output = 'normal_traction tangent_traction normal_jump tangent_jump'
      []
    []
  []
[]

[Functions]
  # Top pin moves +stretch(t), bottom pin moves -stretch(t).
  # End-opening = 2*stretch(t).  Davila Figure 8 reaches ~16 mm.
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
  # Whole left edge fully clamped, both arms.
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
  # Single-node point load: vertical only, x-translation free at the pin.
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
    boundary = 'upper_lower'
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
    boundary = 'upper_lower'
  []
  [max_normal_jump]
    type = SideExtremeValue
    variable = normal_jump
    boundary = 'upper_lower'
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

  # Adaptive time stepping: dt SHRINKS on Newton failure but does NOT
  # grow back. Keeping dt small in the softening branch is critical --
  # large dt lets Newton lock onto a non-physical elastic-like branch
  # (it found peak 73 N at opening 7.6 mm vs the true 58.5 N at 3.84 mm
  # path). growth_factor=1.0 disables growth.
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.25
    cutback_factor = 0.5
    cutback_factor_at_failure = 0.5
  []
  dtmin = 1e-6
  dtmax=0.02
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
  file_base =tier2_1_dcb_cg_davila_wedge_h0p5_out
[]
