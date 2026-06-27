# Scaling weightings

Scales quantitative errors or qualities of orientation measurements.
Useful for some statistical summaries or fault-slip inversion.

## Usage

``` r
scale_weights(
  e,
  error_type = c("rank", "angle", "rup"),
  scaling = c("lin", "inv_lin", "inv_square", "exp"),
  replace_na = TRUE,
  norm = FALSE
)
```

## Arguments

- e:

  numeric. Weights

- error_type:

  character. One of `"rank"` (a numeric value ranking measurement
  quality), `"angle"` (a reading error expressed as angles in degrees),
  and `"rup"` (RUP values from a previous fault-slip inversion)

- scaling:

  character. Scaling function to use. One of `"lin"` (leaves weights
  \\e\\ as is), `"inv_lin"` (\\1/e\\), `"inv_square"` (\\1/e^2\\), and
  `"exp"` (\\\exp{-(x-1)}\\)

- replace_na:

  logical. Imputation? Whether `NA` should be replaced by the mean of
  the weights? (`TRUE` by default)

- norm:

  logical. Whether the scaled weights should be normalized by their
  mean? (`FALSE` by default)

## Value

numeric.

## See also

[`slip_inversion_angelier()`](https://tobiste.github.io/structr/reference/slip_inversion_angelier.md)

## Examples

``` r
set.seed(20250411)
# Generate some random weights from 1 (poor) to 5 (good)
err <- sample(1:5, size = 10, replace = TRUE)

# Introduce 3 random NAs
err[sample(length(err), 3)] <- NA

scale_weights(err, error_type = 'rank', scaling = 'inv_square', norm = TRUE)
#>  [1] 0.8095728 0.8095728 0.8095728 0.4817784 3.0111151 0.7527788 0.7527788
#>  [8] 0.4817784 1.3382734 0.7527788
```
