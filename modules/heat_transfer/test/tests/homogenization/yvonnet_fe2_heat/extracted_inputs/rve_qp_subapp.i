
[Mesh]
  [square_boundary]
    type = PolyLineMeshGenerator
    points = '0 0 0
              1 0 0
              1 1 0
              0 1 0'
    loop = true
    nums_edges_between_points = 24
  []
  [fiber_boundary]
    type = PolyLineMeshGenerator
    points = '0.80946117639254 0.5 0
              0.80681369312908 0.54039278899446 0
              0.79891654251135 0.58009444617022 0
              0.78590484697593 0.61842566516564 0
              0.76800124024096 0.65473058819627 0
              0.74551205798148 0.68838802796344 0
              0.71882209634113 0.71882209634113 0
              0.68838802796344 0.74551205798148 0
              0.65473058819627 0.76800124024096 0
              0.61842566516564 0.78590484697593 0
              0.58009444617022 0.79891654251135 0
              0.54039278899446 0.80681369312908 0
              0.5 0.80946117639254 0
              0.45960721100554 0.80681369312908 0
              0.41990555382978 0.79891654251135 0
              0.38157433483436 0.78590484697593 0
              0.34526941180373 0.76800124024096 0
              0.31161197203656 0.74551205798148 0
              0.28117790365887 0.71882209634113 0
              0.25448794201852 0.68838802796344 0
              0.23199875975904 0.65473058819627 0
              0.21409515302407 0.61842566516564 0
              0.20108345748865 0.58009444617022 0
              0.19318630687092 0.54039278899446 0
              0.19053882360746 0.5 0
              0.19318630687092 0.45960721100554 0
              0.20108345748865 0.41990555382978 0
              0.21409515302407 0.38157433483436 0
              0.23199875975904 0.34526941180373 0
              0.25448794201852 0.31161197203656 0
              0.28117790365887 0.28117790365887 0
              0.31161197203656 0.25448794201852 0
              0.34526941180373 0.23199875975904 0
              0.38157433483436 0.21409515302407 0
              0.41990555382978 0.20108345748865 0
              0.45960721100554 0.19318630687092 0
              0.5 0.19053882360746 0
              0.54039278899446 0.19318630687092 0
              0.58009444617022 0.20108345748865 0
              0.61842566516564 0.21409515302407 0
              0.65473058819627 0.23199875975904 0
              0.68838802796344 0.25448794201852 0
              0.71882209634113 0.28117790365887 0
              0.74551205798148 0.31161197203656 0
              0.76800124024096 0.34526941180373 0
              0.78590484697593 0.38157433483436 0
              0.79891654251135 0.41990555382978 0
              0.80681369312908 0.45960721100554 0'
    loop = true
    nums_edges_between_points = 1
  []
  [fiber]
    type = XYDelaunayGenerator
    boundary = fiber_boundary
    refine_boundary = false
    desired_area = 0.02
    tri_element_type = TRI6
    interior_points = '0.5 0.5 0'
    output_subdomain_id = 2
    output_subdomain_name = fiber
    output_boundary = fiber_boundary
  []
  [cell]
    type = XYDelaunayGenerator
    boundary = square_boundary
    holes = fiber
    stitch_holes = true
    refine_holes = false
    verify_holes = false
    refine_boundary = false
    desired_area = 0.02
    tri_element_type = TRI6
    output_subdomain_id = 1
    output_subdomain_name = matrix
    output_boundary = outer
  []
  [left]
    type = ParsedGenerateSideset
    input = cell
    new_sideset_name = left
    included_boundaries = outer
    combinatorial_geometry = 'x < 1e-10'
  []
  [right]
    type = ParsedGenerateSideset
    input = left
    new_sideset_name = right
    included_boundaries = outer
    combinatorial_geometry = 'x > 1 - 1e-10'
  []
  [bottom]
    type = ParsedGenerateSideset
    input = right
    new_sideset_name = bottom
    included_boundaries = outer
    combinatorial_geometry = 'y < 1e-10'
  []
  [top]
    type = ParsedGenerateSideset
    input = bottom
    new_sideset_name = top
    included_boundaries = outer
    combinatorial_geometry = 'y > 1 - 1e-10'
  []
  [pin]
    type = ExtraNodesetGenerator
    input = top
    new_boundary = pin
    coord = '0.5 0.5 0'
    use_closest_node = true
  []
[]

[Variables]
  [chi_x]
    order = FIRST
    family = LAGRANGE
  []
  [chi_y]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [heat_x]
    type = HeatConduction
    variable = chi_x
  []
  [heat_y]
    type = HeatConduction
    variable = chi_y
  []
  [rhs_x]
    type = HomogenizedHeatConduction
    variable = chi_x
    component = 0
  []
  [rhs_y]
    type = HomogenizedHeatConduction
    variable = chi_y
    component = 1
  []
[]

[BCs]
  [Periodic]
    [left_right]
      primary = left
      secondary = right
      translation = '1 0 0'
    []
    [bottom_top]
      primary = bottom
      secondary = top
      translation = '0 1 0'
    []
  []
  [pin_x]
    type = DirichletBC
    variable = chi_x
    boundary = pin
    value = 0
  []
  [pin_y]
    type = DirichletBC
    variable = chi_y
    boundary = pin
    value = 0
  []
[]

[Materials]
  [matrix]
    type = HeatConductionMaterial
    block = matrix
    specific_heat = 1
    thermal_conductivity = 1
  []
  [fiber]
    type = HeatConductionMaterial
    block = fiber
    specific_heat = 1
    thermal_conductivity = 5
  []
[]


[Executioner]
  type = Transient
  num_steps = 1
  dt = 1
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -ksp_gmres_restart'
  petsc_options_value = 'lu       101'
  line_search = none
  nl_abs_tol = 1e-11
  nl_rel_tol = 1e-10
  l_max_its = 50
[]


[Postprocessors]
  [total_area]
    type = VolumePostprocessor
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [fiber_area]
    type = VolumePostprocessor
    block = fiber
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [actual_fraction]
    type = ParsedPostprocessor
    pp_names = 'fiber_area total_area'
    pp_symbols = 'Af A'
    expression = 'Af/A'
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [k11]
    type = HomogenizedThermalConductivity
    chi = 'chi_x chi_y'
    row = 0
    col = 0
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
[]


[Outputs]
  exodus = false

[]

