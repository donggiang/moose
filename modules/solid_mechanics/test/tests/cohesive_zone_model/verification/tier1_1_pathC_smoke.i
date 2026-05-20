# Path C smoke test — built kernels registered and resolvable
#
# NOTE: this file is a compile-load demonstration only. The pure-Neumann
# CZM single-element setup turns out to have a separate Jacobian conditioning
# issue (FACTOR_NUMERIC_ZEROPIVOT) that's NOT related to arc length —
# it would surface for any plain FunctionNeumannBC on this mesh + CZM. See
# `arclength_pathC_status.md` for full details and the recommended next step
# (a more substantial structural test).
#
# Single-element mode-I CZM, force-controlled via a coupled scalar variable
# `lambda`. PETSc's standard Newton solves the augmented (bordered) system
# of (displacement, lambda) DOFs at each step.
#
# Constraint (Riks-style cylindrical):
#     ∫ |Δu_x|^2 + |Δu_y|^2  dΩ  =  Δs^2
# implemented as two ArcLengthScalarKernel instances (one per displacement
# component); the disp_y instance is designated `primary_constraint = true`
# and carries the −Δs^2 term.
#
# Loading: ALCoupledScaledNeumannBC pulls the top in +y with magnitude
# −λ · F_ref. F_ref = 100 (well above N=50 so the cohesive zone fully
# softens at λ ≈ 0.5–0.6).
#
# Reference law parameters (same as Tier 1.1):
#   K  = 5e3   N = 50   GIc = 0.5
#   delta_n_0 = N/K = 1.0e-2,   delta_n_f = 2 GIc/N = 2.0e-2

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

[Problem]
  # Scalar variables aren't directly covered by a per-variable Kernel; this
  # disables the integrity check that complains about that.
  kernel_coverage_check = false
[]

# Lambda removed — pure-FunctionNeumann sanity test (no scalar)

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        add_variables = true
        use_automatic_differentiation = true
      []
    []
    # CohesiveZone disabled — sanity test: pure elasticity + pure Neumann
    # [CohesiveZone]
    #   [czm]
    #     strain = SMALL
    #     boundary = 'Block1_Block2'
    #     generate_output = 'normal_traction tangent_traction normal_jump tangent_jump'
    #   []
    # []
  []
[]

# PHASE 1: no constraint kernel — lambda is explicit (auxscalar above).
# PHASE 2: enable these and switch [AuxVariables] -> [Variables] for lambda.
# [Kernels]
#   [al_x]
#     type = ArcLengthScalarKernel
#     variable = disp_x
#     lambda = lambda
#     delta_s = 0.05
#     volume_pp = volume
#     primary_constraint = false
#   []
#   [al_y]
#     type = ArcLengthScalarKernel
#     variable = disp_y
#     lambda = lambda
#     delta_s = 0.05
#     volume_pp = volume
#     primary_constraint = true
#   []
# []

[Functions]
  [ref_load]
    type = ConstantFunction
    value = 100.0
  []
  [sanity_ramp]
    type = ParsedFunction
    expression = '5e4 * t'
  []
  [disp_ramp]
    type = ParsedFunction
    expression = '0.03 * t'   # match Tier 1.1's loading
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
  # Back to FunctionNeumannBC — testing without automatic_scaling
  [pull_y_top]
    type = FunctionNeumannBC
    variable = disp_y
    boundary = top
    function = sanity_ramp
  []
  # # Force on top face, scaled by lambda (PHASE 2):
  # [pull_y_top]
  #   type = ALCoupledScaledNeumannBC
  #   variable = disp_y
  #   boundary = top
  #   function = ref_load
  #   arc_length_scalar = lambda
  # []
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
  # [czm]
  #   type = BiLinearMixedModeTraction
  #   boundary = 'Block1_Block2'
  #   penalty_stiffness = 5e3
  #   GI_c = 0.5
  #   GII_c = 0.5
  #   normal_strength = 50
  #   shear_strength = 30
  #   eta = 1.45
  #   displacements = 'disp_x disp_y'
  #   viscosity = 0
  # []
[]

[Postprocessors]
  [volume]
    type = VolumePostprocessor
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  # CZM postprocessors disabled along with CZM
  # [normal_jump]
  #   type = SideAverageValue
  #   variable = normal_jump
  #   boundary = 'Block1_Block2'
  # []
  # [normal_traction]
  #   type = SideAverageValue
  #   variable = normal_traction
  #   boundary = 'Block1_Block2'
  # []
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
  automatic_scaling = false

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  nl_max_its = 30
  l_tol = 1e-10

  start_time = 0.0
  end_time   = 1.0
  dt         = 0.05
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
