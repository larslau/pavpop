<!-- README.md is generated from README.Rmd. Please edit that file -->
pavpop
======

pavpop is an R package for fitting phase and amplitude varying population pattern (pavpop) models on curve data.

So, what does it do?
--------------------

If you have a set of curves that all represent the same underlying object, but are not aligned, pavpop may be able to help. But it can go further than other curve regirstation methods, since it can also register curves that contain systematic differences in shape (i.e. basically any data set of functional data). Consider for example the 12 repetitions of acceleration signals arising from a from a pen writing.

![raw acceleration data](http://i.imgur.com/4HW9lsR.png)

Installation
------------

You can install the latest version of pavpop by running the following code in R.

``` r
if (packageVersion("devtools") < 1.6) {
  install.packages("devtools")
}
devtools::install_github('larslau/pavpop')
```

News
----

### Version 0.9 release note

The method for predicting warps can now be controlled using the argument `warp_optim_method` in `pavpop`. It defaults to conjugate gradient which is markedly slower than before (where Nelder-Mead was used) but does seem to be more robust.

A new option `amp_fct` to `pavpop` gives the possibility of specifying a functional basis in which to express the amplitude variation. As a special case, this opens up the possibility for fitting sitar-like models. For more information, run the vignettes `amplitude_basis.Rmd` and `sitar.Rmd`.

### Version 0.8 release note

The `make_cov_fct` function now allows a new argument `ns` for creating a flexible class of non-stationary covariance functions from a stationary covariance function by adding one or more new parameters describing the change over time. The parameters will be automatically be estimated using maximum likelihood estimation. This is particularly useful in growth curve analysis, where the scale of for example height variation is known to vary with age.

The `pavpop` function now allows the logical argument `warped_amp`, which, if set to true, fits the model where the amplitude effects are modeled in warped time. The likelihood function used for estimation does not yet include a linearization term for the random effect, which may give rise to a slight bias in the likelihood estimates.

### Version 0.7 release note

Parallel computations and clustering are back, and single point samples are now are supported.

### Version 0.6 release note

The package now supports a variety of warping functions (shifts, linear warps, piecewise linear, smooth) and specification of models with fixed warping effects. Parallel computations and clustering are temporarily suspended, but expected to reemerge in the new structure very soon!
