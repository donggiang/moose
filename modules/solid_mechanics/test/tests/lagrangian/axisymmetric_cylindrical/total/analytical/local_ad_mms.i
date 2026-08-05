# Axisymmetric manufactured solution for the local-AD Neo-Hookean material.
# The annulus excludes r = 0, and the diagonal deformation keeps the forcing compact while
# exercising the axisymmetric hoop stretch F_theta_theta = 1 + u_r / r.
# The coefficients a = 0.05 and b = 0.03 produce moderate, spatially varying positive stretches;
# lambda = mu = 400 gives a compressible response without affecting the manufactured solution.

[GlobalParams]
  displacements = 'disp_r disp_z'
  large_kinematics = true
  stabilize_strain = false
[]

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 1
    xmax = 2
    ymin = 0
    ymax = 1
    nx = 8
    ny = 8
  []
  coord_type = RZ
[]

[Variables]
  [disp_r]
  []
  [disp_z]
  []
[]

[Functions]
  [exact_r]
    type = ParsedFunction
    expression = 'a*x^2'
    symbol_names = 'a'
    symbol_values = '0.05'
  []
  [exact_z]
    type = ParsedFunction
    expression = 'b*y^2'
    symbol_names = 'b'
    symbol_values = '0.03'
  []
  [radial_body_force]
    type = ParsedFunction
    expression = '-a*(3*mu + lambda*(2/fr + 1/ft)/fr - 2*q/fr^2 - q/(fr*ft))'
    symbol_names = 'a mu lambda fr ft q'
    symbol_values = '0.05 400 400 radial_stretch hoop_stretch pressure_term'
  []
  [axial_body_force]
    type = ParsedFunction
    expression = '-2*b*(mu + (lambda-q)/fz^2)'
    symbol_names = 'b mu lambda fz q'
    symbol_values = '0.03 400 400 axial_stretch pressure_term'
  []
  [radial_stretch]
    type = ParsedFunction
    expression = '1 + 2*a*x'
    symbol_names = 'a'
    symbol_values = '0.05'
  []
  [axial_stretch]
    type = ParsedFunction
    expression = '1 + 2*b*y'
    symbol_names = 'b'
    symbol_values = '0.03'
  []
  [hoop_stretch]
    type = ParsedFunction
    expression = '1 + a*x'
    symbol_names = 'a'
    symbol_values = '0.05'
  []
  [pressure_term]
    type = ParsedFunction
    expression = 'lambda*log(fr*fz*ft) - mu'
    symbol_names = 'lambda mu fr fz ft'
    symbol_values = '400 400 radial_stretch axial_stretch hoop_stretch'
  []
[]

[Kernels]
  [radial_stress]
    type = TotalLagrangianStressDivergenceAxisymmetricCylindrical
    variable = disp_r
    component = 0
  []
  [axial_stress]
    type = TotalLagrangianStressDivergenceAxisymmetricCylindrical
    variable = disp_z
    component = 1
  []
  [radial_force]
    type = BodyForce
    variable = disp_r
    function = radial_body_force
  []
  [axial_force]
    type = BodyForce
    variable = disp_z
    function = axial_body_force
  []
[]

[BCs]
  [radial_exact]
    type = FunctionDirichletBC
    variable = disp_r
    boundary = 'left right top bottom'
    function = exact_r
  []
  [axial_exact]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'left right top bottom'
    function = exact_z
  []
[]

[Materials]
  [lame_parameters]
    type = GenericConstantMaterial
    prop_names = 'lambda mu'
    prop_values = '400 400'
  []
  [stress]
    type = ComputeLagrangianADNeoHookeanStress
  []
  [strain]
    type = ComputeLagrangianStrainAxisymmetricCylindrical
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  # Resolve the nonlinear system well below the finest-mesh discretization error.
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-12
[]

[Postprocessors]
  [radial_l2_error]
    type = ElementL2Error
    variable = disp_r
    function = exact_r
  []
  [axial_l2_error]
    type = ElementL2Error
    variable = disp_z
    function = exact_z
  []
[]

[Outputs]
  csv = true
[]
