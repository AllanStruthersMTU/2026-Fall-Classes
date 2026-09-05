# Julia experiment scaffold

1. `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
2. `julia --project=. run_experiments.jl`

The direct-Q baselines are wired to `OwrenZen4()` and `Tsit5()`.
`midpoint_extrapolation.jl` contains the fixed-point Cayley midpoint kernel and the two-half/full extrapolated macrostep.
The next implementation pass should add the adaptive macrostep controller, simultaneous log-R accumulation, the step-local Cayley embedded-RK wrapper, all benchmark systems, CSV logging, and standardized plots.

For fair solver comparison, all Q methods should use the same high-accuracy dense state trajectory y(t), or all should integrate the coupled (y,Q,ell) problem with identical state tolerances. The former is preferable for isolating QR/sensitivity integrator behavior.
