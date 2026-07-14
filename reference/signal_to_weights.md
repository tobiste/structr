# Convert a raw signal vector to weights with NA imputation.

Imputation is done before scaling, in the raw signal space.
Normalisation is deliberately omitted — done by combine_weights() or
internally by wissi().

## Usage

``` r
signal_to_weights(signal, scale_fn = function(x) 1/x^2, pessimistic = FALSE)
```

## Arguments

- signal:

  Numeric vector (e.g. field ranks 1-5, errors in degrees). May contain
  NAs.

- scale_fn:

  Function mapping imputed signal -\> raw weights. Default: 1/x^2
  (inverse variance).

- pessimistic:

  If TRUE, impute NAs with the 75th percentile of the observed signal
  (use when missingness is informative).
