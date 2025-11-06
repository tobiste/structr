# The cone or plane best fit of conically or cylindrical disposed s-plane poles

Finding the best fit pole of rotation for a given set of points that are
assumed to lie on a mutual small or great circle circle using Ramsay
1967 algorithm

## Usage

``` r
regression_cone_ramsay(x)

regression_cone_ramsay2(x)
```

## Arguments

- x:

  matrix, where the rows are the observations and the columns are the
  coordinates of points

## Value

numeric vector with

- `x`,`y`,`z`:

  Cartesian coordinates of best fit pole of plane or cone axis,

- `e`:

  residual of the sum of square of the deviations of the observed poles
  to the planes from the best fit pole, and

- `K`:

  (only for cones) half apical angle of best fit cone (in radians).

## References

Ramsay, 1967, p. 18-21

Ramsay, J. G. (1967). Folding and Fracturing of Rocks. McGraw-Hill.

## Examples

``` r
if (FALSE) { # \dontrun{
# example from Ramsay, 1967, p. 20
x <- rbind(
  c(-67, -31, -71),
  c(-62, -53, -50),
  c(-62, -75, -34),
  c(-58, 85, -34),
  c(-79, 40, -52),
  c(90, 14, -75),
  c(80, 10, 90)
) |> acoscartesian_to_cartesian()
regression_cone_ramsay(x) # expect: c(0.856, -0.157, -0.492, NA, 1.56207)
regression_plane_ramsay(x) # expect: c(0.852, -0.154, -0.502, 1-1.002)
} # }
```
