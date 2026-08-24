!include ad_aniso_iso_creep_x_3d.i

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
  [local_iterations_average]
    type = ElementAverageValue
    variable = local_iterations
  []
[]
