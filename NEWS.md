---
editor_options: 
  markdown: 
    wrap: 72
---

# structr (development version)

-   Initial CRAN submission.
-   Error fix in `RGN_plot()`
-   New CI estimation for RGN method: The original interval was the
    standard error of the bootstrap mean, so it narrowed as you added
    resamples, which is not a measure of uncertainty in Wm. The new
    percentile interval spans the actual spread of the bootstrap
    estimates, and it handles θ error through per-clast classification
    probabilities instead of a hard 3° margin. Measurement error still
    isn’t in the interval itself, so the p_min sweep is what shows how
    sensitive the result is to that error.
