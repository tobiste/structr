# Rate and spin of velocity gradient tensor

The velocity gradient tensor **L** can be decomposed into a symmetric
matrix **S** (the rate or stretching tensor) and the skew-symmetric
matrix **W** (the spin or vorticity tensor).

## Usage

``` r
velgrad_rate(x)

velgrad_spin(x)
```

## Source

`apsg` by O. Lexa

## Arguments

- x:

  3x3 matrix. Velocity gradient tensor.

## Value

3x3 matrix

## Details

The velocity gradient tensor\\\mathbf{L}\\ can be decomposed into the
sum of a symmetric matrix \\\mathbf{\dot{S}}\\ and a skew-symmetric
matrix \\\mathbf{W}\\ \$\$\mathbf{L} = \mathbf{\dot{S}} + \mathbf{W}\$\$

where \\\mathbf{\dot{S}}\\ is the stretching tensor (or strain-rate
tensor) that describes the portion of the deformation that over time
produces strain. \\\mathbf{W}\\ is the vorticity or spin tensor and
describes the internal rotation (rate) during the deformation.

## See also

[`velgrad()`](https://tobiste.github.io/structr/reference/gradient.md)

## Examples

``` r
R <- defgrad_from_comp(xx = 2, xy = 1, zz = 0.5)
L <- velgrad(R, time = 10)
velgrad_rate(L)
#> Velocity gradient tensor
#>            [,1]       [,2]        [,3]
#> [1,] 0.06931472 0.03465736  0.00000000
#> [2,] 0.03465736 0.00000000  0.00000000
#> [3,] 0.00000000 0.00000000 -0.06931472
velgrad_spin(L)
#> Velocity gradient tensor
#>             [,1]       [,2] [,3]
#> [1,]  0.00000000 0.03465736    0
#> [2,] -0.03465736 0.00000000    0
#> [3,]  0.00000000 0.00000000    0
```
