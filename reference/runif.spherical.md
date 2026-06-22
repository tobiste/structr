# Uniformly Distributed Vectors on the Sphere

Create uniformly distributed vectors on the sphere

## Usage

``` r
runif.spherical(
  n = 100,
  class = c("Vec3", "Ray", "Line", "Plane"),
  method = c("cart", "gss", "sfs", "rotasym")
)
```

## Arguments

- n:

  integer. number of random samples to be generated

- class:

  character. Coordinate class of the output vectors.

- method:

  character. The algorithm for generating uniformly distributed vectors.
  Either `"cart"` (the default) for generating random points in
  Cartesian coordinates (as in the `geologyGeometry` package), `"sfs"`
  for the "Spherical Fibonacci Spiral points on a sphere", `"gss"` for
  "Golden Section Spiral points on a sphere", or the algorithm
  [`rotasym::r_unif_sphere()`](https://rdrr.io/pkg/rotasym/man/unif.html)
  from the rotasym package.

## Value

object of class specified by `"class"` argument

## Details

`"sfs"` algorithm is from on John Burkardt
(http://people.sc.fsu.edu/~jburkardt/), `"gss` is from
http://www.softimageblog.com/archives/115

## See also

Other random:
[`rbing`](https://tobiste.github.io/structr/reference/rbing.md),
[`rfb()`](https://tobiste.github.io/structr/reference/rfb.md),
[`rkent()`](https://tobiste.github.io/structr/reference/rkent.md),
[`rrot()`](https://tobiste.github.io/structr/reference/rrot.md),
[`vonmises-fisher`](https://tobiste.github.io/structr/reference/vonmises-fisher.md)

## Examples

``` r
set.seed(20250411)
x1 <- runif.spherical(n = 100, "Ray", method = "sfs")
contour(x1)


x2 <- runif.spherical(n = 100, "Ray", method = "gss")
contour(x2)


x3 <- runif.spherical(n = 100, "Ray", method = "rotasym")
contour(x3)


x4 <- runif.spherical(n = 100, "Ray", method = "cart")
contour(x4)
```
