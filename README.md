# pavpop

An R package for fitting phase and amplitude varying population pattern (PAVPOP) models on curve data. 


### Version 0.5 release note
The package now allows parallel computations by setting the `n_cores` argument of `estimate_generic` to the number of cores that should be used.

### Version 0.3 release note
As of version 0.3, the structure of the package has changed. Unfortunately this structure is not compatible with the previous structure, so previous code has to be rewritten to work with newer versions of the package. Furthermore, it is now possible to do probabilistic cluster analysis within the framework.
