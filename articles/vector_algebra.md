# Vector Algebra

This tutorial introduces the basic mathematical operations of the
[structr](https://tobiste.github.io/structr/) package for analyzing
orientation data.

``` r

library(structr)
```

## Vector operations

Since the spherical and vector objects are easily convertible, they can
be used for all sorts of vector operations, such as the magnitude (or
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

> **Note**
>
> Orientation vectors are by definition unit vectors, i.e. their length
> is equal to 1.

The **angle** between two vectors:

``` r

angle(line1, line2)
#> [1] 78.89371
```

The **dot product** (or scalar product) of two vectors:

``` r

dotprod(line1, line2)
#> [1] 0.1926297
```

Intuitively, the dot product tells us how much two vectors point in the
same direction. For unit vectors, it equals the cosine of the angle
between them.

The **cross product** of two vectors:

``` r

crossprod(line1, line2)
#> Line object (n = 1):
#>   azimuth    plunge 
#> 258.66786  32.21399
```

This gives the vector that is perpendicular to the plane spanned by the
two vectors. Geologically, the cross product of two lines is the pole to
the plane containing both lines, and the cross product of two poles to
planes is the line of intersection of the two planes.

The **projection** of a vector onto another vector:

``` r

project(line1, line2)
#> Line object (n = 1):
#> azimuth  plunge 
#>      10      30
```

> **Note**
>
> The projection of a vector onto a unit vector is parallel to that unit
> vector, with a length equal to their dot product (the cosine of the
> angle between them). As an orientation, it therefore points along the
> second vector.

The **rotation** of a vector about another vector (rotation axis) by a
specified rotation angle:

``` r

rotate(line1, rotaxis = line2, rotangle = 45)
#> Line object (n = 1):
#>   azimuth    plunge 
#> 210.50391  70.01332
```

> Right-hand-rule convention: Rotation is counterclockwise for positive
> angles.

**Linear transformation** transforms vectors using a transformation
matrix (second-order tensor):

``` r

trans_mat <- matrix(runif(9), 3, 3)
transform_linear(line1, trans_mat)
#> Vector (Vec3) object (n = 1):
#>         x         y         z 
#> 0.4287424 1.0590220 0.4075062
```
