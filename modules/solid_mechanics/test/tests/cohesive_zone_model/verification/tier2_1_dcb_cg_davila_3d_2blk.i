# Tier 2.1 -- DCB CG-CZM 3D, MOOSE-generated 2-block mesh (upper/lower)
# extended from tier2_1_dcb_cg_davila_2blk.i to 3D.
#
# Geometry: 150 mm x 3.96 mm x 20 mm (width).  The 2D plane-strain pin
# pull (point load) becomes a 3D LINE pull across the through-thickness
# at the outer right edges (y = +/-1.98, x = 150, z in [0, 20]).
# Same material / CZM as the 2D variant.
#
# Mesh: nx=150 ny=4 nz=2  (coarse 3D smoke test; h_x = 1 mm, h_y ~ 1 mm,
# h_z = 10 mm).  Refine ny/nx later if results look right.

[Mesh]
  [base]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 0
    xmax = 150
    ymin = -1.98
    ymax =  1.98
    zmin = 0
    zmax = 20
    nx = 150
    ny = 4
    nz = 4
  []
  [upper]
    type = SubdomainBoundingBoxGenerator
    input = base
    block_id = 1
    block_name = upper
    bottom_left = '0   0    0'
    top_right   = '150 1.98 20'
  []
  [lower]
    type = SubdomainBoundingBoxGenerator
    input = upper
    block_id = 2
    block_name = lower
    bottom_left = '0   -1.98 0'
    top_right   = '150  0    20'
  []
  [break]
    type = BreakMeshByBlockGenerator
    input = lower
    block_pairs = 'upper lower'
    split_interface = true
  []
  # Cohesive zone applied only on the bonded portion (x <= 95) of the
  # 'upper_lower' split interface (face-set on the upper side).
  [bonded_iface]
    type = ParsedGenerateSideset
    input = break
    new_sideset_name = bonded_iface
    combinatorial_geometry = 'abs(y) < 1e-6 & x < 95.0 + 1e-6'
    included_subdomains = 'upper'
  []
  # Clamp left edges of each arm (full 2D faces in 3D).
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
  # Pull-pin LINE nodesets: in 3D the 2D "point pin" becomes a line of
  # nodes along z at (x=150, y=+/-1.98). All nz+1 nodes are loaded.
  [pull_upper]
    type = ParsedGenerateNodeset
    input = left_bot
    new_nodeset_name = pull_upper
    expression = 'abs(x - 150) < 1e-6 & abs(y - 1.98) < 1e-6'
  []
  [pull_lower]
    type = ParsedGenerateNodeset
    input = pull_upper
    new_nodeset_name = pull_lower
    expression = 'abs(x - 150) < 1e-6 & abs(y + 1.98) < 1e-6'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
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
  [disp_z]
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
  [solid_z]
    type = ADStressDivergenceTensors
    variable = disp_z
    component = 2
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
  [left_top_z]
    type = DirichletBC
    variable = disp_z
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
  [left_bot_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'left_bot'
    value = 0
    preset = true
  []
  # Pull pin: u_y prescribed along the line of nodes; u_x and u_z free.
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
    displacements = 'disp_x disp_y disp_z'
    viscosity = 0
  []
[]

[Postprocessors]
  # Total reaction = NodalSum over the entire line of pull nodes.
  # No per-width factor needed -- in 3D the full out-of-plane width is in the mesh.
  [P_top_N]
    type = NodalSum
    variable = fy
    boundary = pull_upper
  []
  [P_bot_N]
    type = NodalSum
    variable = fy
    boundary = pull_lower
  []
  [P_top_abs]
    type = ParsedPostprocessor
    expression = 'abs(P_top_N)'
    pp_names = 'P_top_N'
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
