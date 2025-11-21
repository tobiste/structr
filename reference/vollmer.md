# Fabric Intensity and Shape of Orientation Tensor

Fabric intensity and shape parameters of the orientation tensor based on
Vollmer (1990)

## Usage

``` r
vollmer(x)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

## Value

numeric vector containing the fabric shape and intensity indices:

- `P`:

  Point (Vollmer 1990). Range: (0, 1)

- `G`:

  Girdle (Vollmer 1990). Range: (0, 1)

- `R`:

  Random (Vollmer 1990). Range: (0, 1)

- `B`:

  Cylindricity (Vollmer 1990). Range: (0, 1)

- `C`:

  Cylindricity or Fabric strength (Woodcock 1977). Range: (0, `Inf`)

- `I`:

  Cylindricity or Fabric intensity (Lisle 1985). Range: (0, 5)

- `D`:

  "Distance" from uniformity, linear from R to P, and R to G (Vollmer
  2020). Range: (0, 1). End members are: uniform D = 0, girdle D = 0.5,
  cluster D = 1. The 99% level for a test against uniformity for a
  sample size of 300 is D = 0.1.

- `U`:

  Uniformity statistic of Mardia (1972)

## References

Lisle, Richard J. (1985): "The use of the orientation tensor for the
description and statistical testing of fabrics." Journal of Structural
Geology 7.1: 115-117.

Mardia, Kantilal Varichand (1975): "Statistics of directional data."
Journal of the Royal Statistical Society Series B: Statistical
Methodology 37.3: 349-371.

Vollmer, Frederick W. (1990): "An application of eigenvalue methods to
structural domain analysis." Geological Society of America Bulletin
102.6: 786-791.

Vollmer, Frederick W. (2020): "Representing Progressive Fabric Paths on
a Triangular Plot Using a Fabric Density Index and Crystal Axes
Eigenvector Barycenters." Geological Society of America Abstracts. Vol.
52.

Woodcock, Nigel H. (1977): "Specification of fabric shapes using an
eigenvalue method." Geological Society of America Bulletin 88.9:
1231-1236.

## See also

[`shape_params()`](https://tobiste.github.io/structr/reference/strain_shape.md),
[`ortensor()`](https://tobiste.github.io/structr/reference/ortensor.md),
[`vollmer_plot()`](https://tobiste.github.io/structr/reference/vollmer-plot.md)

## Examples

``` r
set.seed(20250411)
mu <- Line(120, 50)
x <- rvmf(100, mu = mu, k = 1)
vollmer(x)
#>           P           G           R           B           C           I 
#>  0.13716550  0.06777223  0.79506227  0.20493773  0.49800008  2.45026231 
#>           D           U 
#>  0.15687782 12.30532499 

# Pair objects:
vollmer(simongomez)
#>            P            G            R            B            C            I 
#>   0.06798346   0.93155564   0.20395038   0.99953910   2.18056889   2.42554552 
#>            D            U 
#>   0.52339268 109.57596054 
```
