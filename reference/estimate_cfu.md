# Estimate CFU from a dilution series

Applies either
[`find_mle()`](https://slhogle.github.io/CFUestimatoR/reference/find_mle.md)
(the default) or
[`find_poisson_cutoff()`](https://slhogle.github.io/CFUestimatoR/reference/find_poisson_cutoff.md)
to each group of a data frame of colony counts and converts the result
to CFU/ml.

## Usage

``` r
estimate_cfu(
  data,
  count = count,
  dilution = dilution,
  group = NULL,
  method = c("mle", "poisson"),
  N = NULL,
  sampling_vol = 0.1
)
```

## Arguments

- data:

  A data frame with one row per counted spot or plate. Input data is
  expected to have a header row, e.g. as read by
  [`utils::read.csv()`](https://rdrr.io/r/utils/read.table.html) or
  `readr::read_csv()`.

- count, dilution:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  The single column holding colony counts and the single column holding
  dilution factors, e.g. `count = n_colonies`. Strings such as
  `"n_colonies"` also work. Default to columns named `count` and
  `dilution`. Dilutions must be numeric (e.g. `1e-6`).

- group:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Columns to group by. Any combination of tidyselect expressions can be
  used, e.g. `c(strain, rep)`, `starts_with("treat")`,
  `where(is.character)` or `all_of(vars)`. Leave as `NULL` for no
  grouping. Strings such as `"strain"` also work.

- method:

  Estimator to use: `"mle"` for
  [`find_mle()`](https://slhogle.github.io/CFUestimatoR/reference/find_mle.md)
  (default) or `"poisson"` for
  [`find_poisson_cutoff()`](https://slhogle.github.io/CFUestimatoR/reference/find_poisson_cutoff.md).

- N:

  Threshold passed to the chosen estimator. For `"mle"`, the maximum
  number of countable colonies (see
  [`find_mle()`](https://slhogle.github.io/CFUestimatoR/reference/find_mle.md));
  defaults to 5000, as recommended for 100 mm petri dishes. For
  `"poisson"`, counts at or above `N` are excluded (see
  [`find_poisson_cutoff()`](https://slhogle.github.io/CFUestimatoR/reference/find_poisson_cutoff.md));
  defaults to 300, also as recommended for 100 mm petri dishes.

- sampling_vol:

  Volume plated or spotted, in milliliters (ml), used to convert
  absolute CFU to CFU/ml. Either a single positive number applied to all
  rows (e.g. `0.002` for a 2 ul spot), or the name of a column holding
  the volume for each row (e.g. `vol_plated_ml`). Defaults to 0.1 (100
  ul). The volume must be constant within each group.

## Value

A data frame with one row per group: the grouping columns, `estimator`
(the estimator used), `CFU_ml` (CFU/ml) and `CFU_ml_stderr` (its
standard error, also per ml).

## Details

Absolute CFU in the amount plated or spotted is divided by
`sampling_vol`, which is always in milliliters, so `CFU_ml` is in
CFU/ml. The default of 0.1 (100 ul) is the usual volume plated on a 100
mm petri dish (use 0.002 for a 2 ul spot). The same volume must be used
for all dilutions.

## Examples

``` r
path <- system.file("extdata", "cfu_example.tsv", package = "CFUestimatoR")
cfu <- read.delim(path)

# Plate counts on 100 mm dishes; replicates are pooled within each group
estimate_cfu(
  cfu,
  count = colony_count,
  group = c(measurement, time_hours, growth_media, sample_or_strain),
  N = 5000,
  sampling_vol = vol_plated_ml
)
#>               measurement time_hours growth_media sample_or_strain    estimator
#> 1 SGT calibration plating         48      PPY+Str          control ML_estimator
#> 2          direct plating         48          PPY               E1 ML_estimator
#> 3          direct plating         48          PPY               E2 ML_estimator
#> 4          direct plating         48      PPY+Str               E1 ML_estimator
#> 5          direct plating         48      PPY+Str               E2 ML_estimator
#> 6  starting-density check          0          PPY            1977A ML_estimator
#> 7  starting-density check          0      PPY+Str            1977E ML_estimator
#>       CFU_ml CFU_ml_stderr
#> 1 2117968546     110561892
#> 2 1469416498      82274187
#> 3 2154889790      99937499
#> 4  360327325      40540075
#> 5  589537227      51906132
#> 6 3179163302     121933373
#> 7 3891949271     135366618

# Poisson cutoff estimator, grouping by any tidyselect expression
estimate_cfu(
  cfu,
  count = colony_count,
  group = c(sample_or_strain, growth_media),
  method = "poisson",
  N = 500,
  sampling_vol = vol_plated_ml
)
#>   sample_or_strain growth_media      estimator     CFU_ml CFU_ml_stderr
#> 1            1977A          PPY poisson_cutoff 3090909091     118530953
#> 2            1977E      PPY+Str poisson_cutoff 3759090909     130716399
#> 3               E1          PPY poisson_cutoff 1450000000      81184414
#> 4               E1      PPY+Str poisson_cutoff  359090909      40400884
#> 5               E2          PPY poisson_cutoff 2113636364      98017539
#> 6               E2      PPY+Str poisson_cutoff  586363636      51626440
#> 7          control      PPY+Str poisson_cutoff 2085227273     108847978
```
