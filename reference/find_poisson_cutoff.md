# Poisson cutoff estimator of CFU

Pools all observations with counts below `N` (counts at or above `N` are
distrusted) and estimates the rate as total count divided by total
dilution times volume.

## Usage

``` r
find_poisson_cutoff(counts, dil, N = 300, V = 1)
```

## Arguments

- counts:

  Numeric vector of colony counts.

- dil:

  Numeric vector of dilution factors (e.g. `1e-5`), same length as
  `counts`.

- N:

  Counts greater than or equal to `N` are excluded. Defaults to 300, the
  value recommended for 100 mm petri dishes. Use `Inf` to keep all
  counts (naive Poisson estimate).

- V:

  Sampling volume. The default of 1 returns absolute CFU in the amount
  plated or spotted.

## Value

A one-row data frame with columns `r`, `stderr` and `type`. `r` is `NA`
if no observations fall below `N`.

## Details

Setting `N = Inf` excludes nothing and gives the naive Poisson estimate,
which pools every observation.

## Examples

``` r
counts <- c(2, 4, 3, 25, 32)
dil <- c(1e-6, 1e-6, 1e-6, 1e-5, 1e-5)
find_poisson_cutoff(counts, dil, N = 30)
#>         r   stderr           type
#> 1 2615385 448534.8 poisson_cutoff

# No cutoff: naive Poisson estimate
find_poisson_cutoff(counts, dil, N = Inf)
#>         r   stderr           type
#> 1 2869565 353219.1 poisson_cutoff
```
