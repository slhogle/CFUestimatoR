# input validation errors are informative

    Code
      find_poisson_cutoff(1:2, 1, N = 50)
    Condition
      Error:
      ! `counts` and `dil` must be the same length.

---

    Code
      find_mle(-1, 1e-05, N = 50)
    Condition
      Error:
      ! `counts` must be >= 0 and `dil` must be > 0.

---

    Code
      find_mle("a", 1e-05, N = 50)
    Condition
      Error:
      ! `counts` and `dil` must be numeric.

# estimate_cfu validates arguments

    Code
      estimate_cfu(d, group = nope)
    Condition
      Error in `estimate_cfu()`:
      ! Can't select columns that don't exist.
      x Column `nope` doesn't exist.

---

    Code
      estimate_cfu(d, count = nope)
    Condition
      Error in `estimate_cfu()`:
      ! Can't select columns that don't exist.
      x Column `nope` doesn't exist.

---

    Code
      estimate_cfu(d, count = c(count, dilution))
    Condition
      Error:
      ! `count` must select exactly one column, not 2.

---

    Code
      estimate_cfu(d, group = CFU_ml)
    Condition
      Error:
      ! `group` cannot include columns named `CFU_ml`; rename them first.

---

    Code
      estimate_cfu(d, sampling_vol = -1)
    Condition
      Error:
      ! `sampling_vol` must be a single positive number or the name of a numeric column of positive volumes in ml.

---

    Code
      estimate_cfu(d, sampling_vol = nope)
    Condition
      Error:
      ! object 'nope' not found

---

    Code
      estimate_cfu(varying, sampling_vol = v)
    Condition
      Error:
      ! `sampling_vol` must be constant within each group, but varies in at least one group.

---

    Code
      estimate_cfu(d, method = "bogus")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "mle", "poisson"

