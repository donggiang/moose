# Tier 2.1 -- DCB FULL-DG with DG cohesive zone (ADDGCohesiveZone).
#
# Companion to tier2_1_dcb_cg_davila_2blk.i (CG-CZM).
# Key differences from the CG variant:
#   * MONOMIAL p=1 displacements (DG)
#   * NO BreakMeshByBlock -- DG already has jumps at every face
#   * NO [Physics/SolidMechanics/CohesiveZone] block
#   * ADDGElasticity on intra-subdomain interior faces (penalty + consistency)
#   * ADDGCohesiveZone on the bonded interface (Neumann-type t* = T_coh, u* = u)
#   * No kernel at all on the precrack interface -- free surfaces by absence
#   * Abedi star-value Dirichlet BCs on left clamp and right loading pin
#
# Pin is a thin face on the outer edge (y in [1.85, 1.98] for upper arm,
# y in [-1.98, -1.85] for lower) -- ADDGElasticityAbediDirichletBC is a face-
# based weak BC (no nodes in DG), so we use a small face rather than a true
# point node.
#
# Mesh: nx = 300, ny = 8 (matches the CG h=0.5 mm refinement).
# Material identical to tier2_1_dcb_cg_davila_2blk.i.

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
    elem_type = QUAD4
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
  # bonded_iface: y=0 sub-face for x <= 95, both sides tagged (upper element
  # with lower neighbor) so the bulk DG kernel can exclude it bidirectionally.
  [bonded_iface]
    type = ParsedGenerateSideset
    input = lower
    new_sideset_name = bonded_iface
    combinatorial_geometry = 'abs(y) < 1e-6 & x < 95.0 + 1e-6'
    included_subdomains = 'upper'
    included_neighbors = 'lower'
  []
  # precrack_iface: y=0 sub-face for x > 95. No kernel applied here -> the
  # bulk DG penalty must be excluded to let the precrack open freely.
  [precrack_iface]
    type = ParsedGenerateSideset
    input = bonded_iface
    new_sideset_name = precrack_iface
    combinatorial_geometry = 'abs(y) < 1e-6 & x > 95.0 - 1e-6'
    included_subdomains = 'upper'
    included_neighbors = 'lower'
  []
  # Pin BC faces on the outer edges of each arm.
  # With ny=8 the top-element centroid on the right edge is y=1.7325; the
  # threshold 'y > 1.4' picks exactly the top-most element face (and its
  # negative for right_bot), giving a small face analogous to the cantilever
  # NIPG file's right_tab.
  [right_top]
    type = ParsedGenerateSideset
    input = precrack_iface
    new_sideset_name = right_top
    combinatorial_geometry = 'abs(x - 150) < 1e-6 & y > 1.4'
  []
  [right_bot]
    type = ParsedGenerateSideset
    input = right_top
    new_sideset_name = right_bot
    combinatorial_geometry = 'abs(x - 150) < 1e-6 & y < -1.4'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
    family = MONOMIAL
    order = FIRST
  []
  [disp_y]
    family = MONOMIAL
    order = FIRST
  []
[]

[Kernels]
  [stress_div_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
  []
  [stress_div_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
  []
[]

# Bulk DG kernel on intra-subdomain interior faces only. Exclude the bonded
# interface so the DG penalty does NOT compete with the cohesive law there;
# the precrack portion of the upper-lower face is ALSO excluded so the two
# arms separate freely.
[DGKernels]
  [dg_x]
    type = ADDGElasticity
    variable = disp_x
    component = 0
    epsilon = +1
    eta = 5
    penalty_modulus = G_mat
    exclude_boundary = 'bonded_iface precrack_iface'
  []
  [dg_y]
    type = ADDGElasticity
    variable = disp_y
    component = 1
    epsilon = +1
    eta = 5
    penalty_modulus = G_mat
    exclude_boundary = 'bonded_iface precrack_iface'
  []
[]

# Cohesive interface kernel (Neumann-type star traction, t* = T_coh, u* = u).
# ONE InterfaceKernel per displacement component.
[InterfaceKernels]
  [czm_x]
    type = ADDGCohesiveZone
    variable = disp_x
    neighbor_var = disp_x
    boundary = 'bonded_iface'
    component = 0
    displacements = 'disp_x disp_y'
    penalty_stiffness = 1e6
    normal_strength = 30.0
    shear_strength = 40.0
    GI_c = 0.268
    GII_c = 1.45
    eta_BK = 2.284
    viscosity = 1e-3
  []
  [czm_y]
    type = ADDGCohesiveZone
    variable = disp_y
    neighbor_var = disp_y
    boundary = 'bonded_iface'
    component = 1
    displacements = 'disp_x disp_y'
    penalty_stiffness = 1e6
    normal_strength = 30.0
    shear_strength = 40.0
    GI_c = 0.268
    GII_c = 1.45
    eta_BK = 2.284
    viscosity = 1e-3
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
  [zero]
    type = ConstantFunction
    value = 0
  []
[]

[BCs]
  # Whole left edge weakly clamped via Abedi star-value Dirichlet.
  [left_x]
    type = ADDGElasticityAbediDirichletBC
    variable = disp_x
    component = 0
    boundary = 'left'
    function = 'zero zero'
    epsilon = +1
  []
  [left_y]
    type = ADDGElasticityAbediDirichletBC
    variable = disp_y
    component = 1
    boundary = 'left'
    function = 'zero zero'
    epsilon = +1
  []
  # Right-tab pin: u_y prescribed, u_x left as natural (no constraint).
  [pull_top_y]
    type = ADDGElasticityAbediDirichletBC
    variable = disp_y
    component = 1
    boundary = 'right_top'
    function = 'zero stretch'
    epsilon = +1
  []
  [pull_bot_y]
    type = ADDGElasticityAbediDirichletBC
    variable = disp_y
    component = 1
    boundary = 'right_bot'
    function = 'zero stretch_neg'
    epsilon = +1
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
  [shear_modulus]
    # G = E / (2 (1 + nu)) = 150e3 / 2.5 = 60e3 MPa
    type = ADGenericConstantMaterial
    prop_names = 'G_mat'
    prop_values = '60e3'
  []
  # Diagnostic-only interface material: writes TSL state at face QPs so
  # MaterialRealAux can project to the boundary AuxVariables. Mirrors
  # ADDGCohesiveZone's constitutive law (uses the same parameters).
  [czm_output]
    type = ADDGCohesiveZoneOutput
    boundary = 'bonded_iface'
    displacements = 'disp_x disp_y'
    penalty_stiffness = 1e6
    normal_strength = 30.0
    shear_strength = 40.0
    GI_c = 0.268
    GII_c = 1.45
    eta_BK = 2.284
  []
[]

[AuxVariables]
  [sxy]
    family = MONOMIAL
    order = CONSTANT
  []
  # Diagnostic cohesive-state aux variables. MONOMIAL CONSTANT so each
  # interface-touching element gets a value; visible as a colored strip
  # along the bonded interface in ParaView.
  [normal_jump]
    family = MONOMIAL
    order = CONSTANT
  []
  [shear_jump]
    family = MONOMIAL
    order = CONSTANT
  []
  [czm_damage]
    family = MONOMIAL
    order = CONSTANT
  []
  [normal_traction]
    family = MONOMIAL
    order = CONSTANT
  []
  [shear_traction]
    family = MONOMIAL
    order = CONSTANT
  []
  [czm_active_zone]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [sxy_aux]
    type = ADRankTwoAux
    variable = sxy
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    execute_on = 'TIMESTEP_END'
  []
  # Project the InterfaceMaterial output to the aux variables on bonded_iface.
  # check_boundary_restricted=false so the elemental aux var can be touched
  # by a boundary-restricted aux kernel.
  [aux_normal_jump]
    type = MaterialRealAux
    variable = normal_jump
    property = normal_jump
    boundary = 'bonded_iface'
    check_boundary_restricted = false
    execute_on = 'TIMESTEP_END'
  []
  [aux_shear_jump]
    type = MaterialRealAux
    variable = shear_jump
    property = shear_jump
    boundary = 'bonded_iface'
    check_boundary_restricted = false
    execute_on = 'TIMESTEP_END'
  []
  [aux_czm_damage]
    type = MaterialRealAux
    variable = czm_damage
    property = czm_damage
    boundary = 'bonded_iface'
    check_boundary_restricted = false
    execute_on = 'TIMESTEP_END'
  []
  [aux_normal_traction]
    type = MaterialRealAux
    variable = normal_traction
    property = normal_traction
    boundary = 'bonded_iface'
    check_boundary_restricted = false
    execute_on = 'TIMESTEP_END'
  []
  [aux_shear_traction]
    type = MaterialRealAux
    variable = shear_traction
    property = shear_traction
    boundary = 'bonded_iface'
    check_boundary_restricted = false
    execute_on = 'TIMESTEP_END'
  []
  [aux_czm_active_zone]
    type = MaterialRealAux
    variable = czm_active_zone
    property = czm_active_zone
    boundary = 'bonded_iface'
    check_boundary_restricted = false
    execute_on = 'TIMESTEP_END'
  []
[]

[Postprocessors]
  # Reactions at the upper- and lower-arm pin faces -- match dg_cantilever_NIPG.i.
  # ADSidesetReaction integrates n.(sigma.dir) on a face; verified eta-sweet-spot
  # (eta ~ 10-15) gives the correct DG-Abedi reaction within ~3% of CG NodalSum.
  [P_top_per_width]
    type = ADSidesetReaction
    direction = '0 1 0'
    stress_tensor = stress
    boundary = right_top
  []
  [P_bot_per_width]
    type = ADSidesetReaction
    direction = '0 1 0'
    stress_tensor = stress
    boundary = right_bot
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
  # Prescribed (target) and actual tip displacements.  In DG the Abedi BC
  # enforces u_y = stretch(t) weakly, so use SideAverageValue for the actual
  # value rather than the function value.
  [tip_disp]
    type = FunctionValuePostprocessor
    function = stretch
  []
  [tip_disp_top]
    type = SideAverageValue
    variable = disp_y
    boundary = right_top
  []
  [tip_disp_bot]
    type = SideAverageValue
    variable = disp_y
    boundary = right_bot
  []
  [opening]
    type = ParsedPostprocessor
    expression = 'tip_disp_top - tip_disp_bot'
    pp_names   = 'tip_disp_top tip_disp_bot'
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
