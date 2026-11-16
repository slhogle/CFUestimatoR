# CFUestimatoR

`CFUestimatoR` is a package for implementing the Maximum Likelihood and
Poisson estimators from:

Martini KM, Boddu SS, Nemenman I, Vega NM. 2024. Maximum likelihood
estimators for colony-forming units. *Microbiology Spectrum*
12(9):e03946-23. <https://doi.org/10.1128/spectrum.03946-23>

Make sure to cite this paper if you use CFUestimatoR.

## Installation

You can install the development version from GitHub with:

``` r

# install.packages("pak")
pak::pkg_install("slhogle/CFUestimatoR")
```

## Estimating colony forming units

[`estimate_cfu()`](https://slhogle.github.io/CFUestimatoR/reference/estimate_cfu.md)
estimates CFU/ml from a dilution series of colony counts. Input data
must have a header row, with one row per counted spot or plate. The
bundled example, `cfu_example.tsv`, has columns for the sample, media,
`dilution`, `vol_plated_ml` and `colony_count`. Replicates are pooled by
leaving `replicate` out of `group`.

``` r

library(CFUestimatoR)

path <- system.file("extdata", "cfu_example.tsv", package = "CFUestimatoR")
cfu <- read.delim(path)

# Maximum likelihood estimator (default), using the plated volume in each row
estimate_cfu(
  cfu,
  count = colony_count,
  group = c(measurement, time_hours, growth_media, sample_or_strain),
  N = 5000,
  sampling_vol = vol_plated_ml
)

#>               measurement time_hours growth_media sample_or_strain    estimator
#> 1          direct plating         48          PPY               E1 ML_estimator
#> 2          direct plating         48          PPY               E2 ML_estimator
#> 3          direct plating         48      PPY+Str               E1 ML_estimator
#> 4          direct plating         48      PPY+Str               E2 ML_estimator
#> 5 SGT calibration plating         48      PPY+Str          control ML_estimator
#> 6  starting-density check          0          PPY            1977A ML_estimator
#> 7  starting-density check          0      PPY+Str            1977E ML_estimator
#>       CFU_ml CFU_ml_stderr
#> 1 1469416498      82274187
#> 2 2154889790      99937499
#> 3  360327325      40540075
#> 4  589537227      51906132
#> 5 2117968546     110561892
#> 6 3179163302     121933373
#> 7 3891949271     135366618

# Poisson cutoff estimator
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
#> 3          control      PPY+Str poisson_cutoff 2085227273     108847978
#> 4               E1          PPY poisson_cutoff 1450000000      81184414
#> 5               E1      PPY+Str poisson_cutoff  359090909      40400884
#> 6               E2          PPY poisson_cutoff 2113636364      98017539
#> 7               E2      PPY+Str poisson_cutoff  586363636      51626440
```

Notes:

- `sampling_vol` is always in milliliters. Defaults to 0.1 (100 ul), the
  usual volume plated on a 100 mm dish. Give either a single number
  (e.g. `0.002` for a 2 ul spot) or a column holding the volume for each
  row (e.g. `vol_plated_ml`). The volume must be the same for all
  dilutions within a group.
- `N` is the estimator threshold. For `"mle"` it is the maximum
  countable colonies (default 5000, recommended for a 100 mm dish; about
  50 for a 2 ul spot). Counts above `N` give an estimate of `Inf` for
  `"mle"`, so set `N` to suit your assay. For `"poisson"` counts at or
  above `N` are excluded (default 300, recommended for a 100 mm dish).
  Use `N = Inf` to exclude nothing, which gives the naive Poisson
  estimate.
- `group` accepts any tidyselect expression,
  e.g. `group = c(sample_or_strain, growth_media)`,
  `starts_with("treat")` or `where(is.character)`. Omit it for no
  grouping.
- The count and dilution columns default to `count` and `dilution`. Pass
  `count` and `dilution` a column name from your data if yours differ,
  e.g. `count = colony_count`.
