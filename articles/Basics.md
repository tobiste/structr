# Vector Algebra

This tutorial introduces the basic mathematical operations of the
[structr](https://tobiste.github.io/structr/) package for analyzing
orientation data.

``` r
library(structr)
```

## Vector operations

Since the spherical or vector objects are easily convertible, they can
be used for all sort of vector operations, such as the magnitude (or
length), the angle between vectors, dot product, cross product,
projection and rotation.

Define some example vectors:

``` r
line1 <- Line(120, 50)
line2 <- Line(10, 30)
```

The **vector length** (or magnitude):

``` r
vector_length(line1)
#> [1] 1
```

> Orientation vectors are by definition unit vectors, i.e. their length
> is equal to 1.

The **angle** between two vectors

``` r
angle(line1, line2)
#> [1] 78.89371
```

The **dot product** (or scalar product) of two vectors

``` r
dotprod(line1, line2)
#> [1] 0.1926297
```

Intuitively, the dot product tells us how much two vectors point in the
same direction.

The **cross product** of two vectors:

``` r
crossprod(line1, line2)
#> Line object (n = 1):
#>   azimuth    plunge 
#> 258.66786  32.21399
```

This gives the vector that is perpendicular to the plane spanned by the
two vectors.

The **projection** of a vector on another vector:

``` r
project(line1, line2)
#> Line object (n = 1):
#> azimuth  plunge 
#>      10      30
```

> Because the vectors are both unit vectors, the projected vector is
> equal to the second vector.

The **rotation** of a vector about another vector (rotation axis) by a
specified rotation angle:

``` r
rotate(line1, rotaxis = line2, rotangle = 45)
#> Line object (n = 1):
#>   azimuth    plunge 
#> 210.50391  70.01332
```

**Linear transformation** transforms vectors using a transformation
matrix (second-order tensor).

``` r
trans_mat <- matrix(runif(9), 3, 3)
transform_linear(line1, trans_mat)
#> Vector (Vec3) object (n = 1):
#>           x           y           z 
#>  0.44288029 -0.04205583  0.62796680
```
