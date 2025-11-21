# Helper functions to assign plotting colors to a vector

Helper functions to assign plotting colors to a vector

## Usage

``` r
assign_col_d(x, pal = viridis::viridis, ...)

assign_col(x, n = length(x), pal = viridis::viridis, ...)

assign_col_binned(x, breaks = 5, pal = viridis::viridis, ...)

legend_col(breaks, title = NULL, pal = viridis::viridis, cex = 1, ...)

legend_col_d(fill, legend = names(fill), position = "topright", ...)
```

## Arguments

- x:

  vector to colorize

- pal:

  color function; Default is
  [`viridis::viridis()`](https://sjmgarnier.github.io/viridis/reference/reexports.html)

- ...:

  arguments passed to color function

- n:

  integer. The number of colors (\\\ge\\1) to be in the palette.

- breaks:

  integer giving the desired number of intervals. Non-integer values are
  rounded down.

- title:

  character. Legend title

- cex:

  character expansion factor relative to current par("cex"). Used for
  text in legend.

- fill:

  color vector

- legend:

  character vector. Names of discrete colors. Can be ignored when `cols`
  is a named vector.

- position:

  Legend position. Either a two-column vector of the x and y
  coordinates, or a keyword from the list `"bottomright"`, `"bottom"`,
  `"bottomleft"`, `"left"`, `"topleft"`, `"top"`, `"topright"`,
  `"right"` and `"center"`.

## Value

character vector of colors in hexadecimal code

## See also

Other assign:
[`assign-cex`](https://tobiste.github.io/structr/reference/assign-cex.md),
[`assign-pch`](https://tobiste.github.io/structr/reference/assign-pch.md)

## Examples

``` r
set.seed(20250411)

# example for discrete colors
x <- rvmf(5, mu = Line(120, 50), k = 5)
key <- letters[round(runif(5, 1, 26))]
plot(x, col = assign_col_d(key), grid.params = list(guides = FALSE))
legend_col_d(assign_col_d(key))


# example for continuous colors:
x <- rvmf(100, mu = Line(120, 50), k = 5)
plot(x, col = assign_col(runif(100)), grid.params = list(guides = FALSE))
legend_col(seq(0, 1, .1), title = "test")
```
