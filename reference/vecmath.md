# Vector math operations

Vector math operations

## Usage

``` r
vector_length(x)

vector_norm(x)

# S3 method for class 'Vec3'
crossprod(x, y = NULL, ...)

# S3 method for class 'Line'
crossprod(x, y = NULL, ...)

# S3 method for class 'Ray'
crossprod(x, y = NULL, ...)

# S3 method for class 'Plane'
crossprod(x, y = NULL, ...)

dotprod(x, y)

angle(x, y)

# S3 method for class 'spherical'
angle(x, y)

project(x, y)

# S3 method for class 'spherical'
project(x, y)

reject(x, y)

# Default S3 method
reject(x, y)

transform_linear(x, A, norm = FALSE)
```

## Arguments

- x, y:

  objects of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

- ...:

  arguments passed to function call

- A:

  numeric 3x3 matrix. Transformation matrix.

- norm:

  logical. If `TRUE`, the transformed vectors are normalized to unit
  length.

## Value

objects of same class as `x`, i.e. one of `"Vec3"`, `"Line"`, `"Ray"` or
`"Plane"`. `vector_length()` and `%*%` return a real number. `angle()`
returns a numeric angle (in degrees, unless `x` is class `"Vec3"`).

## Details

- `vector_length()`:

  the length of a vector: \\\|\|x\|\| = \sqrt{x_1^2 + x_2^2 + x_3^2}\\

- `vector_norm()`:

  the normalized vector: \\\hat{x} = \frac{x}{\|\|x\|\|}\\

- [`crossprod()`](https://rdrr.io/r/base/crossprod.html):

  the cross-product of two vectors, i.e. the vector perpendicular to the
  2 vectors. If `y = NULL` is taken to be the sam,e vector as `x`: \$\$x
  \times y = (x_2 y_3 - x_3 y_2, x_3 y_1 - x_1 y_3, x_1 y_2 - x_2
  y_1)\$\$

- `dotprod()`:

  the dot product of two vectors: \\x \cdot y = x_1 y_1 + x_2 y_2 + x_3
  y_3\\

- `angle()`:

  angle between two vectors: \\\theta = \arccos{\frac{x \cdot
  y}{\|\|x\|\| \\ \|\|y\|\|}}\\

- `project()`:

  projection of one vector onto the other (changes the vector length of
  second vector, unless their are unit vectors): \$\$proj_y(x) = \frac{x
  \cdot y}{\|\|y\|\|^2} y\$\$

- `transform_linear()`:

  Linear transformation of a vector by a 3x3 matrix: \\x' = A x\\

## Examples

``` r
vec1 <- Vec3(1, 0, 0)
vec2 <- Vec3(0, 0, 1)

vector_length(vec1) # length of a vector
#> [1] 1
crossprod(vec1, vec2) # cross product
#> Vector (Vec3) object (n = 1):
#>  x  y  z 
#>  0 -1  0 
dotprod(vec1, vec2) # dot product
#> [1] 0
angle(vec1, vec2) # angle between vectors
#> [1] 1.570796
project(vec1, vec2) # projection of a vector
#> Vector (Vec3) object (n = 1):
#> x y z 
#> 0 0 0 
transform_linear(vec1, matrix(runif(9), 3, 3)) # linear transformation
#> Vector (Vec3) object (n = 1):
#>         x         y         z 
#> 0.1714069 0.1418151 0.1022888 
```
