# Principal Stretches, Strain and Shape Parameters based on the Orientation Tensor.

Principal Stretches, Strain and Shape Parameters based on the
Orientation Tensor.

## Usage

``` r
principal_stretch(x, ...)

# Default S3 method
principal_stretch(x, ...)

# S3 method for class 'Vec3'
principal_stretch(x, ...)

# S3 method for class 'Line'
principal_stretch(x, ...)

# S3 method for class 'Ray'
principal_stretch(x, ...)

# S3 method for class 'Plane'
principal_stretch(x, ...)

# S3 method for class 'ortensor'
principal_stretch(x, ...)

# S3 method for class 'ellipsoid'
principal_stretch(x, ...)

principal_strain(x)

shape_params(x, ...)

# Default S3 method
shape_params(x, ...)

# S3 method for class 'spherical'
shape_params(x, ...)

# S3 method for class 'ortensor'
shape_params(x, ...)

# S3 method for class 'ellipsoid'
shape_params(x, ...)
```

## Arguments

- x:

  object of class `"ortensor"`, `"ellipsoid"`, or `"Vec3"`, `"Line"`,
  `"Ray"`, `"Plane"`, `"Pair"`, or `"Fault"` where the rows are the
  observations and the columns are the coordinates.

- ...:

  optional parameters passed to
  [`ortensor()`](https://tobiste.github.io/structr/reference/ortensor.md)

## Value

list

## Details

- `stretch_ratios`:

  Sqrt of eigenvalue ratios

- `strain_ratios`:

  Log of stretch ratios

- `Ramsay`:

  strain symmetry (Ramsay, 1983)

- `Woodcock`:

  Woodcock shape

- `Flinn`:

  Flinn strain intensity

- `Nadai`:

  natural octahedral unit strain and shear (Nadai, 1963)

- `Lisle_intensity`:

  Intensity index (Lisle, 1985)

- `Waterson_intensity`:

  strain intensity (Watterson, 1968)

- `lode`:

  Lode parameter (Lode, 1926)

- `kind`:

  Descriptive type of ellipsoid: `"O"` - isotrope, `"L"` - L-tectonite,
  `"LLS"` - oblate L-tectonite, `"S"` - S-tectonite, `"SSL"` - prolate
  S-tectonite

- `MAD`:

  maximum angular deviation (Kirschvink, 1980)

## References

Flinn, Derek.(1963): "On the statistical analysis of fabric diagrams."
Geological Journal 3.2: 247-253.

Jelinek, Vit. "Characterization of the magnetic fabric of rocks."
Tectonophysics 79.3-4 (1981): T63-T67.

Kirschvink, J. (1980): The least-squares line and plane and the analysis
of palaeomagnetic data. Geophysical Journal International, 62(3),
699-718.

Lisle, Richard J. (1985): "The use of the orientation tensor for the
description and statistical testing of fabrics." Journal of Structural
Geology 7.1: 115-117.

Lode, Walter (1926): "Versuche über den Einfluß der mittleren
Hauptspannung auf das Fließen der Metalle Eisen, Kupfer und Nickel“
(*"Experiments on the influence of the mean principal stress on the flow
of the metals iron, copper and nickel"*\], Zeitschrift für Physik, vol.
36 (November), pp. 913–939, DOI: 10.1007/BF01400222

Mardia, Kantilal Varichand. (1975): "Statistics of directional data."
Journal of the Royal Statistical Society Series B: Statistical
Methodology 37.3: 349-371.

Nadai, A., and Hodge, P. G., Jr. (1963): "Theory of Flow and Fracture of
Solids, vol. II." ASME. J. Appl. Mech. December 1963; 30(4): 640.
https://doi.org/10.1115/1.3636654

Ramsay, John G. (1967): "Folding and fracturing of rocks." Mc Graw Hill
Book Company 568.

Vollmer, Frederick W. (1990): "An application of eigenvalue methods to
structural domain analysis." Geological Society of America Bulletin
102.6: 786-791.

Vollmer, Frederick W. (2020): "Representing Progressive Fabric Paths on
a Triangular Plot Using a Fabric Density Index and Crystal Axes
Eigenvector Barycenters." Geological Society of America Abstracts. Vol.
52.

Watterson, Juan. (1968): "Homogeneous deformation of the gneisses of
Vesterland, south-west Greenland". No. 78. CA Reitzel.

Woodcock, N. H. (1977): "Specification of fabric shapes using an
eigenvalue method." Geological Society of America Bulletin 88.9:
1231-1236.

## See also

[`vollmer()`](https://tobiste.github.io/structr/reference/vollmer.md)
for Vollmer 1990 shape parameters of the orientation tensor

More details on shape parameters:
[`lode()`](https://tobiste.github.io/structr/reference/ellipsoid-params.md),
[`nadai()`](https://tobiste.github.io/structr/reference/ellipsoid-params.md),
[`jelinek()`](https://tobiste.github.io/structr/reference/ellipsoid-params.md),
[`flinn()`](https://tobiste.github.io/structr/reference/ellipsoid-params.md).
Details on eigenvectors from orientation tensors:
[`ot_eigen()`](https://tobiste.github.io/structr/reference/ot_eigen.md).

Other ortensor:
[`ortensor()`](https://tobiste.github.io/structr/reference/ortensor.md),
[`ot_eigen()`](https://tobiste.github.io/structr/reference/ot_eigen.md)

## Examples

``` r
set.seed(1)
mu <- Line(120, 50)
x <- rvmf(100, mu = mu, k = 20)
principal_stretch(x)
#>        S1        S2        S3 
#> 0.9511749 0.2220799 0.2143520 
principal_strain(x)
#>          e1          e2          e3 
#> -0.05005729 -1.50471812 -1.54013586 
shape_params(x)
#> $stretch_ratios
#>      Rxy      Ryz      Rxz 
#> 4.283031 1.036052 4.437444 
#> 
#> $strain_ratios
#>        e12        e13        e23 
#> 1.45466083 1.49007857 0.03541774 
#> 
#> $Flinn
#> $Flinn$k
#> [1] 91.0627
#> 
#> $Flinn$d
#> [1] 3.283228
#> 
#> 
#> $Ramsay
#> intensity  symmetry 
#>  2.117293 41.071535 
#> 
#> $Woodcock
#>  strength     shape 
#>  1.490079 41.071535 
#> 
#> $Watterson_intensity
#> [1] 4.319083
#> 
#> $Lisle_intensity
#> [1] 3.67315
#> 
#> $Nadai
#>     goct     eoct 
#> 1.388465 1.202446 
#> 
#> $Lode
#> [1] -0.9524619
#> 
#> $kind
#> [1] "L"
#> 
#> $MAD
#> [1] 17.97803
#> 
#> $Jellinek
#> [1] 5.476767
#> 
```
