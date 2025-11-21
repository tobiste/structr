# Fault displacement tensor

Creates the fault displacement tensor from displacement components. If
the dip direction is know, the tensor will be rotated into the
geographic reference frame.

## Usage

``` r
displacement_tensor(h, s, v, dip_direction = NULL)

displacement_tensor_decomposition(ftensor, dip_direction = NULL)
```

## Arguments

- h, s, v:

  numeric. the heave, strike-slip, vertical throw displacement

- dip_direction:

  (optional) dip direction in degrees.

- ftensor:

  Fault displacement tensor. A 3x3 matrix. If `NULL`, the fault tensor
  will be given in the fault displacement coordinates. Otherwise, the
  tensor will be in the geographic reference frame.

## Value

`displacement_tensor()` returns a 3x3 matrix of class `"ftensor"`
containing the fault displacement tensor.
`displacement_tensor_decomposition()` returns a list containing the
principal fault displacement tensor and the fault orientation.

## Details

x axis of tensor = heave, y = strike slip, z = vertical throw (positive
for thrusting, negative for normal faulting) is the principal fault
displacement tensor. This can be rotated in the fault plane orientation
to retrieve slip components and rake.

The fault displacement tensor gives the displacements in all directions
and has the following properties:

- The square root of tensor's trace (i.e. the sum of the diagonal
  elements) equals the **net slip** on the fault plane.

- The determinant of the tensor relates to the volumetric strain by:
  det(F) - 1, where

(`displacement_tensor_decomposition()`) retrieves the principal fault
displacement tensor using Singular Value Decomposition of a Matrix and
the fault orientation if the dip direction is known.

The orientation of the net-slip vector is the lineation component of the
fault orientation.

## Examples

``` r
A_princ <- displacement_tensor(s = 2, v = -5, h = 3)
print(A_princ)
#>      [,1] [,2] [,3]
#> [1,]    3    0    0
#> [2,]    0    2    0
#> [3,]    0    0   -5
#> attr(,"class")
#> [1] "matrix"  "array"   "ftensor"
det(A_princ)
#> [1] -30

A_geo <- displacement_tensor(s = 2, v = -5, h = 3, dip_direction = 45)
print(A_geo)
#>         [,1]      [,2] [,3]
#> [1,] 2.12132  1.414214    0
#> [2,] 2.12132 -1.414214    0
#> [3,] 0.00000  0.000000   -5
#> attr(,"class")
#> [1] "matrix"  "array"   "ftensor"
det(A_geo)
#> [1] 30

displacement_tensor_decomposition(A_geo, dip_direction = 45)
#> $displacements
#>           dip    delta     rake verticalthrow horizontalthrow heave  dipslip
#> [1,] 59.03624 56.30993 71.06818            -5        3.605551     3 5.830952
#>      strikeslip  netslip
#> [1,]          2 6.164414
#> 
#> $fault
#> Fault object (n = 1):
#> dip_direction           dip       azimuth        plunge         sense 
#>      45.00000      59.03624      78.69007      54.20424      -1.00000 
#> 
#> $strain_tensor
#>      [,1] [,2] [,3]
#> [1,]  3.0 -1.5 -1.5
#> [2,] -1.0  1.0 -1.0
#> [3,]  2.5  2.5 15.0
#> 
#> $volumetric_strain
#> [1] 19
#> 
#> $shear_strain
#> [1] 5.196152
#> 
```
