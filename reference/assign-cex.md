# Assign plotting size (cex values) to a vector

`assign_cex()` maps the character expansion (size). The cex is most
commonly used for points and text, and humans perceive the area of
points (not their radius), so this provides for optimal perception. The
argument `area` ensures that a value of `0` is mapped to a size of `0`;
`assign_cex_binned()` is a binned version, and `assign_cex_d()` assigns
cex values to discrete values.

## Usage

``` r
assign_cex(x, range = c(0.25, 2), area = FALSE)

assign_cex_binned(x, range = c(0.25, 2), breaks = 5, area = FALSE)

assign_cex_d(x, values = NULL, range = c(0.25, 2))

legend_cex(
  x,
  range = c(0.25, 2),
  breaks = 5,
  values = NULL,
  area = FALSE,
  position = "topright",
  pch = 16,
  ...
)
```

## Arguments

- x:

  vector

- range:

  numeric 2-element vector. Output range of `cex` values. Minimum value
  must be greater than 0.

- area:

  logical. Whether `cex` should be proportional to the area (`TRUE`) or
  the radius (`FALSE`, the default) of the plotting character.

- breaks:

  integer giving the desired number of intervals. Non-integer values are
  rounded down.

- values:

  numeric. `cex` values to manually assign to `x`. Must be at least the
  number of unique values in `x`.

- position:

  Legend position. Either a two-column vector of the x and y
  coordinates, or a keyword from the list `"bottomright"`, `"bottom"`,
  `"bottomleft"`, `"left"`, `"topleft"`, `"top"`, `"topright"`,
  `"right"` and `"center"`.

- pch:

  plotting character to be used in legend.

- ...:

  arguments passed to
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html)

## Value

numeric vector

## Details

The character expansion `cex` is a number indicating the amount by which
plotting text and symbols should be scaled relative to the default.
`1`=default, `1.5` is 50% larger, `0.5` is 50% smaller, etc.

## See also

Other assign:
[`assign-color`](https://tobiste.github.io/structr/reference/assign-color.md),
[`assign-pch`](https://tobiste.github.io/structr/reference/assign-pch.md)

## Examples

``` r
set.seed(20250411)

# example for continuous colors:
x <- rvmf(100, mu = Line(120, 50), k = 5)
key <- runif(100)
plot(x, cex = assign_cex(key), grid.params = list(guides = FALSE))
legend_cex(key, position = 'topright', area = TRUE)
```
