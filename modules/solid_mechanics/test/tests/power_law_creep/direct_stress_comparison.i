!include ad_smallstrain.i

[AuxVariables]
  [local_iterations]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [local_iterations]
    type = MaterialRealAux
    variable = local_iterations
    property = radial_return_local_iterations
    execute_on = 'initial timestep_end'
  []
[]

[Postprocessors]
  [creep_strain_yy_average]
    type = ElementAverageValue
    variable = creep_strain_yy
  []
  [local_iterations_average]
    type = ElementAverageValue
    variable = local_iterations
  []
[]

[Outputs]
  csv = true
[]
