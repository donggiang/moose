
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 8
  ny = 8
  xmax = 1
  ymax = 1
  elem_type = QUAD4
[]

[Variables]
  [T]
    order = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  [exact]
    type = ParsedFunction
    expression = 'sin(pi*x)*sin(pi*y)'
  []
  [source_shape]
    type = ParsedFunction
    expression = 'sin(pi*x)*sin(pi*y)'
  []
[]

[Kernels]
  [heat]
    type = HeatConduction
    variable = T
  []
  [source]
    type = BodyForce
    variable = T
    function = source_shape
    value = 30.06424226692854
  []
[]

[BCs]
  [exact_boundary]
    type = FunctionDirichletBC
    variable = T
    boundary = 'left right bottom top'
    function = exact
  []
[]

[Postprocessors]
  [nodal_error]
    type = NodalL2Error
    variable = T
    function = exact
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [elemental_error]
    type = ElementL2Error
    variable = T
    function = exact
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [center_temperature]
    type = PointValue
    variable = T
    point = '0.5 0.5 0'
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
[]

[VectorPostprocessors]
  [centerline]
    type = LineValueSampler
    variable = T
    start_point = '0 0.5 0'
    end_point = '1 0.5 0'
    num_points = 41
    sort_by = x
    execute_on = FINAL
  []
[]

[Outputs]
  exodus = false
  file_base = 'macro_reference'
  [csv]
    type = CSV
    execute_on = FINAL
  []
[]

[Materials]
  [effective_conductivity]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity specific_heat density'
    prop_values = '1.5230723058976 1 1'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -ksp_gmres_restart'
  petsc_options_value = 'lu       101'
  line_search = none
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-11
  l_max_its = 50
[]
