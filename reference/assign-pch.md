# Assigns plotting characters (pch values) to a vector

`assign_pch()` maps discrete variables to six easily discernible shapes.
If you have more than six levels, you will get a warning message, and
the seventh and subsequent levels will not appear on the plot. You can
not map a continuous variable to shape unless `assign_pch_binned()` is
used.

## Usage

``` r
assign_pch(x, solid = TRUE)

assign_pch_binned(x, solid = TRUE, breaks = 6)

legend_pch(x, solid = TRUE, breaks = NULL, position = "topright", ...)
```

## Arguments

- x:

  vector.

- solid:

  Should the plotting character be solid, `TRUE` (the default), or
  hollow, `FALSE`?

- breaks:

  integer giving the desired number of intervals. Non-integer values are
  rounded down.

- position:

  Legend position. Either a two-column vector of the x and y
  coordinates, or a keyword from the list `"bottomright"`, `"bottom"`,
  `"bottomleft"`, `"left"`, `"topleft"`, `"top"`, `"topright"`,
  `"right"` and `"center"`.

- ...:

  arguments passed to
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html)

## Value

named integer vector

## See also

Other assign:
[`assign-cex`](https://tobiste.github.io/structr/reference/assign-cex.md),
[`assign-color`](https://tobiste.github.io/structr/reference/assign-color.md)

## Examples

``` r
set.seed(20250411)

# example for discrete colors
x <- rvmf(5, mu = Line(120, 50), k = 5)
key <- sample(letters, 5, replace = TRUE)
plot(x, pch = assign_pch(key), grid.params = list(guides = FALSE))
legend_pch(key)
```
