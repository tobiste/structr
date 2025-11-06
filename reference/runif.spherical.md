# Uniformly distributed vectors

Create uniformly distributed vectors using the algorithm *Spherical
Fibonacci Spiral points on a sphere* algorithm (John Burkardt) or
*Golden Section Spiral points on a sphere*.

## Usage

``` r
runif.spherical(
  n = 100,
  class = c("Vec3", "Line", "Plane"),
  method = c("gss", "sfs", "rotasym")
)
```

## Arguments

- n:

  integer. number of random samples to be generated

- class:

  character. Coordinate class of the output vectors.

- method:

  character. The algorithm for generating uniformly distributed vectors.
  Either `"sfs"` for the "Spherical Fibonacci Spiral points on a
  sphere", `"gss"` for "Golden Section Spiral points on a sphere", or
  the algorithm
  [`rotasym::r_unif_sphere()`](https://rdrr.io/pkg/rotasym/man/unif.html)
  from the rotasym package.

## Value

object of class specified by `"class"` argument

## Details

`"sfs"` algorithm is from on John Burkardt
(http://people.sc.fsu.edu/~jburkardt/), `"gss` is from
http://www.softimageblog.com/archives/115

## See also

[`rvmf()`](https://tobiste.github.io/structr/reference/vonmises-fisher.md)
to draw samples from the von Mises Fisher distribution around a
specified mean vector.
[`rkent()`](https://tobiste.github.io/structr/reference/rkent.md) to
draw from a Kent-distribution.
[`rfb()`](https://tobiste.github.io/structr/reference/rfb.md) to draw
from a Fisher-Bingham distribution.

## Examples

``` r
set.seed(20250411)
x1 <- runif.spherical(n = 100, "Line", method = "sfs")
plot(x1)


x2 <- runif.spherical(n = 100, "Line", method = "gss")
plot(x2)

x3 <- runif.spherical(n = 100, "Line", method = "rotasym")
plot(x3)
```
