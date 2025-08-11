#This tests the J-Integral evaluation capability.
#This is a 2d plane strain model
#The analytic solution for J1 is 2.434.  This model
#converges to that solution with a refined mesh.
#Reference: National Agency for Finite Element Methods and Standards (U.K.):
#Test 1.1 from NAFEMS publication "Test Cases in Linear Elastic Fracture
#Mechanics" R0020.


[Mesh]
  file = square_h025.e
[]
[GlobalParams]
  order = FIRST
  family = LAGRANGE
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = true
[]

[Physics/SolidMechanics/QuasiStatic]
  [./master]
    strain = FINITE
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz vonmises_stress'
    planar_formulation = PLANE_STRAIN
  [../]
[]


[AuxVariables]
  [./SERD]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./SED]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./SERD]
    type = MaterialRealAux
    variable = SERD
    property = strain_energy_rate_density
    execute_on = timestep_end
  [../]
  [./SED]
    type = MaterialRealAux
    variable = SED
    property = strain_energy_density
    execute_on = timestep_end
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 200000
    poissons_ratio = 0.3
  [../]
  [./radial_return_stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'powerlawcrp'
  [../]
  [./powerlawcrp]
    type = PowerLawCreepStressUpdate
    coefficient = 5.0e-12 #3.125e-21 # 7.04e-17 # what does it means?
    n_exponent =3.0# 5.4
    m_exponent = 0.0
    activation_energy = 0.0
  [../]
[]

[Functions]
  [./rampConstant]
    type = PiecewiseLinear
    x = '0. 1e-12 1.0'
    y = '0. 1.0 1.0'
    scale_factor = -1e2
  [../]
  [dt_func]
    type = PiecewiseLinear
    x = '0.    1e-12 1e20'
    y = '1e-12 0.1 0.1'
  [../]
[]
[BCs]
  [./crack_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.0
  [../]

  [./no_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0.0
  [../]

  [./Pressure]
    [./Side1]
      boundary = 'top'
      function = rampConstant
    [../]
  [../]
[]



[DomainIntegral]
  integrals = CIntegral #'JIntegral InteractionIntegralKI'
  boundary = 'tip'
  crack_direction_method = CrackDirectionVector
  crack_direction_vector = '1 0 0'
  2d = true
  axis_2d = 2
  q_function_rings=true
  q_function_type=Topology
  ring_first=1
  ring_last=6
  incremental = true
  inelastic_models = 'powerlawcrp'
  output_q = false
  symmetry_plane = 1
  output_vpp = false
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  line_search = 'none'

  l_max_its = 50
  nl_max_its = 40

  nl_rel_step_tol= 1e-8
  nl_rel_tol = 1e-12

  start_time = 0.0
  #dt =0.05

  end_time = 300
  num_steps =3000
  [TimeStepper]
    type = FunctionDT
    function = dt_func
    min_dt = 1e-12
  []
[]

[Outputs]
  exodus = true
  csv = true
[]

[Preconditioning]
  [./smp]
    type = SMP
    pc_side = left
    ksp_norm = preconditioned
    full = true
  [../]
[]
