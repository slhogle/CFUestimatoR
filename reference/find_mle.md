# Maximum likelihood estimator of CFU

Maximum likelihood estimator that accounts for saturation of a spot or
plate at roughly `N` countable colonies (a most probable number style
estimate).

## Usage

``` r
find_mle(counts, dil, N = 5000, V = 1)
```

## Arguments

- counts:

  Numeric vector of colony counts.

- dil:

  Numeric vector of dilution factors (e.g. `1e-5`), same length as
  `counts`.

- N:

  Estimated maximum number of countable colonies, e.g. the area of the
  plate or spot divided by the average colony size. Defaults to 5000,
  the value recommended for 100 mm petri dishes (about 50 for a 2 ul
  spot).

- V:

  Sampling volume. The default of 1 returns absolute CFU in the amount
  plated or spotted.

## Value

A one-row data frame with columns `r`, `stderr` and `type`. If any count
exceeds `N` the estimate is `Inf`. If all counts are zero the estimate
is 0 with an `NA` standard error.

## Examples

``` r
find_mle(c(2, 4, 3, 25, 32), c(1e-6, 1e-6, 1e-6, 1e-5, 1e-5), N = 100)
#>         r   stderr         type
#> 1 3308623 408870.9 ML_estimator
```
