# ADDirectPowerLawCreepStressUpdate

!syntax description /Materials/ADDirectPowerLawCreepStressUpdate

`ADDirectPowerLawCreepStressUpdate` integrates isotropic von Mises power-law creep with a
backward-Euler local solve whose unknown is the six-component stress tensor. The accumulated
`creep_strain` material property stores the viscoplastic strain history. This provides a direct
stress implementation of the general viscoplastic update in Box 11.4 of de Souza Neto, Peric, and
Owen and is intended for comparison with the scalar radial-return formulation.

!syntax parameters /Materials/ADDirectPowerLawCreepStressUpdate

!syntax inputs /Materials/ADDirectPowerLawCreepStressUpdate

!syntax children /Materials/ADDirectPowerLawCreepStressUpdate

