# Orientation matrix from fault planes and slip directions

Converts a set of planes and lines into a list of rotation matrices in
SO(3).

## Usage

``` r
pair2rot(p)

is.Rotation(x)

as.Rotation(x)

rot2pair(x, fault = FALSE)
```

## Arguments

- p:

  object of class `"Pair"` or `"Fault"`

- x:

  object of class `"Rotation"`, a 3x3 matrix

- fault:

  logical. Whether to coerces to a fault or a pair object.

## Value

list of rotation matrices

## Examples

``` r
my_fault <- Fault(
  c("a" = 120, "b" = 120, "c" = 100),
  c(60, 60, 50),
  c(110, 25, 30),
  c(58, 9, 23),
  c(1, -1, 1)
)
pair2rot(my_fault)
#> [[1]]
#>            [,1]       [,2]       [,3]
#> [1,]  0.4306074 -0.7432627 0.51199399
#> [2,] -0.1752467  0.4876291 0.85528153
#> [3,] -0.8853620 -0.4580158 0.07972233
#> 
#> [[2]]
#>            [,1]       [,2]       [,3]
#> [1,]  0.3674773 -0.7890351  0.4923252
#> [2,] -0.8695531 -0.4792705 -0.1190680
#> [3,]  0.3299058 -0.3843481 -0.8622289
#> 
#> [[3]]
#>            [,1]       [,2]      [,3]
#> [1,]  0.1290480 -0.7567321 0.6408613
#> [2,]  0.7965464  0.4640353 0.3875372
#> [3,] -0.5906441  0.4604648 0.6626550
#> 
#> attr(,"class")
#> [1] "Rotation" "list"    
```
