# Snap-back CZM benchmark — ARC-LENGTH (Box 4.4) version.
#
# Same geometry / material as snap_back_disp.i; difference is the loading:
# instead of prescribed displacement on the top, we apply a λ-scaled Neumann
# traction (ALCoupledScaledNeumannBC) and let the Box 4.4 arc-length solver
# determine λ at each step. This traces the entire (u, F) curve including
# the snap-back branch where standard displacement-control diverges.

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
  extra_tag_vectors = 'arc_length_load'
[]

[AuxVariables]
  [lambda]
    family = SCALAR
    order = FIRST
    initial_condition = 0.0
  []
[]

# Tiny initial offset on Block2 disp_y so the cohesive interface starts with
# delta_n > 0. This avoids the bilinear law's degenerate zero-tangent corner
# (BiLinearMixedModeTraction.C line 112: ddelta_active_ddelta(0,0) = 0 at
# delta_n=0). With delta_n > 0 from t=0 the normal tangent is (1-d)*K and
# K_T is regular.
[ICs]
  [shift_block2]
    type = ConstantIC
    variable = disp_y
    block = 2
    value = 1.0e-7
  []
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
        generate_output = 'normal_traction normal_jump tangent_traction tangent_jump'
      []
    []
  []
[]

[Functions]
  [ref_load]
    type = ConstantFunction
    value = 100.0
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
    type = ALCoupledScaledNeumannBC
    variable = disp_y
    boundary = top
    function = ref_load
    arc_length_scalar = lambda
    extra_vector_tags = 'arc_length_load'
  []
[]

[Materials]
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1.0e4
    poissons_ratio = 0.3
  []
  [stress]
    type = ADComputeLinearElasticStress
  []
  [czm]
    type = BiLinearMixedModeTraction
    boundary = 'Block1_Block2'
    penalty_stiffness = 1.0e5
    GI_c   = 0.25
    GII_c  = 0.25
    normal_strength = 100.0
    shear_strength  = 100.0
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
  [top_disp]
    type = NodalExtremeValue
    variable = disp_y
    boundary = top
    value_type = max
  []
  [lambda_pp]
    type = ScalarVariable
    variable = lambda
  []
  [reaction_y]
    type = ADSidesetReaction
    direction = '0 1 0'
    stress_tensor = stress
    boundary = top
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = ArcLengthTransient
  solve_type = NEWTON
  scalar_variable = lambda
  delta_s = 0.005           # arc-length step in displacement L2 norm
  arc_length_load_tag = arc_length_load
  al_max_iter = 30
  al_rel_tol = 1e-6

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  start_time = 0.0
  end_time = 60.0          # 60 outer arc-length steps
  dt = 1.0
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
