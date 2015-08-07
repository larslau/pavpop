# pavpop

An R package for fitting phase and amplitude varying population pattern (PAVPOP) models on curve data. 

### Version 0.8 release note
The `make_cov_fct` function now allows a new argument `ns` for creating a flexible class of non-stationary covariance functions from a stationary covariance function by adding one or more new parameters describing the change over time. The parameters will be automatically be estimated using maximum likelihood estimation. This is particularly useful in growth curve analysis, where the scale of for example height variation is known to vary with age.

The `pavpop` function now allows the logical argument `warped_amp`, which, if set to true, fits the model where the amplitude effects are modeled in warped time. The likelihood function used for estimation does not yet include a linearization term for the random effect, which may give rise to a slight bias in the likelihood estimates.

### Version 0.7 release note
Parallel computations and clustering are back, and single point samples are now are supported.

### Version 0.6 release note
The package now supports a variety of warping functions (shifts, linear warps, piecewise linear, smooth) and specification of models with fixed warping effects. Parallel computations and clustering are temporarily suspended, but expected to reemerge in the new structure very soon!
