# Getting started with CFUestimatoR

``` r

library(CFUestimatoR)
```

## Overview

`CFUestimatoR` estimates colony forming units (CFU) from colony counts
across a dilution series, using the maximum likelihood estimators from
Martini et al. (2024). Please cite this paper if you use `CFUestimatoR`.

Martini KM, Boddu SS, Nemenman I, Vega NM. 2024. Maximum likelihood
estimators for colony-forming units. *Microbiology Spectrum*
12(9):e03946-23. <https://doi.org/10.1128/spectrum.03946-23>

`CFUestimatoR` offers two estimators:

- `"mle"` (default): a maximum likelihood estimator that accounts for
  saturation of a spot or plate at roughly `N` countable colonies.
- `"poisson"`: a Poisson estimator that pools all observations with
  counts below the cutoff `N`.

## A single dilution series

The estimators work directly on vectors of counts and dilutions. Here
are three replicates at each of three dilutions. Setting `V = 1` returns
CFU in the volume plated (i.e., no back calculation to concentration in
the original sample).

``` r

dilution <- rep(c(1e-7, 1e-6, 1e-5), each = 3)
count <- c(5, 8, 4, 62, 77, 51, 456, 597, 356)

find_mle(count, dilution, N = 5000, V = 1)
#>          r  stderr         type
#> 1 50741067 1262704 ML_estimator
```

The Poisson estimator discards counts at or above `N`, which protects
against crowded plates. Here the three plates at the 1e-5 dilution
(counts of 356 to 597) are excluded, so only the six plates at 1e-7 and
1e-6 are used.

``` r

find_poisson_cutoff(count, dilution, N = 300, V = 1)
#>          r  stderr           type
#> 1 62727273 4359847 poisson_cutoff
```

Setting `N = Inf` includes all observations (even those with very high
colony counts) and gives the naive Poisson estimate from the paper.

``` r

find_poisson_cutoff(count, dilution, N = Inf, V = 1)
#>          r  stderr           type
#> 1 48528529 1207192 poisson_cutoff
```

The naive estimate is lower than the cutoff estimate because colonies on
crowded plates are undercounted. The maximum likelihood estimator uses
all the data and models the saturation, so it needs no cutoff but it
does need a sensible `N`.

### Choosing N

`N` means different things for the two different estimators.

| Method      | Meaning of `N`                      | Default |
|-------------|-------------------------------------|---------|
| `"mle"`     | Maximum countable colonies          | 5000    |
| `"poisson"` | Counts at or above `N` are excluded | 300     |

- The cutoff number in the Poisson approach should be a number less than
  when crowding starts to be important on the plate.
- For `"mle"`, a good starting estimate is the ratio of the area of the
  plate to the average colony size. Any count above `N` gives an
  estimate of `Inf`.

*Note:* sometimes you may use a drop-spotting approach (e.g., [Miles and
Misra](https://en.wikipedia.org/wiki/Miles_and_Misra_method)) if you
have miniaturized your assay. For example if you’ve plated a 2 μl spot
then a value of 50 is a good choice. In general, the drop spotting
approach will usually produce noisier data.

## Many samples at once

[`estimate_cfu()`](https://slhogle.github.io/CFUestimatoR/reference/estimate_cfu.md)
applies one estimator to every group in a data frame. The bundled
example data has one row per plate.

``` r

path <- system.file("extdata", "cfu_example.tsv", package = "CFUestimatoR")
cfu <- read.delim(path)
head(cfu)
#>              measurement time_hours growth_media sample_or_strain replicate
#> 1 starting-density check          0          PPY            1977A         a
#> 2 starting-density check          0          PPY            1977A         b
#> 3 starting-density check          0          PPY            1977A         a
#> 4 starting-density check          0          PPY            1977A         b
#> 5 starting-density check          0      PPY+Str            1977E         a
#> 6 starting-density check          0      PPY+Str            1977E         b
#>   dilution vol_plated_ml colony_count
#> 1    1e-06           0.1          337
#> 2    1e-06           0.1          259
#> 3    1e-07           0.1           40
#> 4    1e-07           0.1           44
#> 5    1e-06           0.1          341
#> 6    1e-06           0.1          391
```

Name the count column with `count` (the dilution column is `dilution`
here, which is the default). Use `group` to say which columns define a
sample. Leaving `replicate` out of `group` pools replicates.

``` r

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
```

`sampling_vol` is always in ml and defaults to 0.1 (100 ul), which is
the standard volume plated on a 100 mm petri dish. It can be a single
number or, as here, a column giving the volume plated for each row
(e.g. `vol_plated_ml`). Note: The volume must be constant within each
group!

Choose the Poisson estimator with `method = "poisson"`, where `N` is now
the cutoff:

``` r

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

## Notes:

- `group` accepts any tidyselect expression,
  e.g. `group = c(sample_or_strain, growth_media)`,
  `starts_with("treat")` or `where(is.character)`. Omit it for no
  grouping.
- The count and dilution columns default to `count` and `dilution`. Pass
  `count` and `dilution` a column name from your data if yours differ,
  e.g. `count = colony_count`.
