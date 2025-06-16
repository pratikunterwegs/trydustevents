
<!-- README.md is generated from README.Rmd. Please edit that file -->

# trydustevents

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the development version of trydustevents from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("pratikunterwegs/trydustevents")
```

## Modifying ODEs

This is a multi-step process involving a `dust2` object:

1.  Modify the ODE system in `inst/dust/model.cpp` — this is a member
    function of the `sirode` class. This class name can be changed, but
    remember to change it within `R/model.R` as well as it must be
    consistent throughout the package;

2.  Run `dust2::dust_package(".")` from within the package directory;

3.  Run `devtools::load_all()` as usual with any other package.
