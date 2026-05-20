# Arc-length kernel + BC test on a 2D bar in tension.
# Simplest possible setup to exercise the new objects without CZM:
#
#                 +----------+
#         clamp ->|          | <- λ·F_ref pulled in +x
#                 |          |
#                 +----------+
#                x=0          x=L
#
# Linear elastic, plane strain. The right-end NeumannBC is supplied by
# `ALCoupledScaledNeumannBC` and scaled by the scalar variable `lambda`.
# In Phase 1 below `lambda` is an AuxScalarVariable driven by t — pure
# force control with explicit λ — to verify the BC works cleanly. Phase 2
# (commented) replaces the AuxScalar with a Variable + the
# `ArcLengthScalarKernel` constraint to test the augmented system.

[Mesh]
  [bar]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 10
    ymin = 0
    ymax = 1
    nx = 20
    ny = 4
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Problem]
  kernel_coverage_check = false   # lambda needs no per-var Kernel here
[]

# Phase 1: lambda is an AuxScalar driven by t — sanity check
#   ALCoupledScaledNeumannBC works (verified: PASSES, see status doc).
# Phase 2 (commented): lambda is a Variable + ArcLengthScalarKernel.
#   Hits first-step degeneracy of cylindrical formulation. Needs custom
#   Predictor (see arclength_pathC_status.md, recommendation #2).
[AuxVariables]
  [lambda]
    family = SCALAR
    order = FIRST
    initial_condition = 0.0
  []
[]

[AuxScalarKernels]
  [lambda_ramp]
    type = FunctionScalarAux
    variable = lambda
    function = 't'
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        add_variables = true
        use_automatic_differentiation = true
        generate_output = 'stress_xx'
      []
    []
  []
[]

[BCs]
  # Left end fully clamped.
  [fix_left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
    preset = true
  []
  [fix_left_y]
    type = DirichletBC
    variable = disp_y
    boundary = left
    value = 0
    preset = true
  []
  # Right end pulled in +x by λ · F_ref.
  [pull_right]
    type = ALCoupledScaledNeumannBC
    variable = disp_x
    boundary = right
    function = ref_load
    arc_length_scalar = lambda
  []
[]

[Functions]
  [ref_load]
    type = ConstantFunction
    value = 1.0e3      # MPa traction reference (corresponds to ~1 % strain at λ=1)
  []
[]

[Materials]
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1.0e5    # MPa
    poissons_ratio = 0.3
  []
  [stress]
    type = ADComputeLinearElasticStress
  []
[]

[Postprocessors]
  [lambda_pp]
    type = ScalarVariable
    variable = lambda
  []
  [tip_disp_x]
    type = PointValue
    variable = disp_x
    point = '10 0.5 0'
  []
  [stress_xx_avg]
    type = ElementAverageValue
    variable = stress_xx
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
  nl_abs_tol = 1e-10
  nl_max_its = 20
  l_tol = 1e-10

  start_time = 0.0
  end_time   = 1.0
  dt         = 0.1
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
