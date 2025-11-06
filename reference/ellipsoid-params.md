# Ellipsoid shape parameters

Ellipsoid shape parameters

## Usage

``` r
volume(x)

# Default S3 method
volume(x)

# S3 method for class 'ellipsoid'
volume(x)

# S3 method for class 'ortensor'
volume(x)

lode(x)

# Default S3 method
lode(x)

# S3 method for class 'ellipsoid'
lode(x)

# S3 method for class 'ortensor'
lode(x)

nadai(x)

# Default S3 method
nadai(x)

# S3 method for class 'ellipsoid'
nadai(x)

# S3 method for class 'ortensor'
nadai(x)

jelinek(x)

# Default S3 method
jelinek(x)

# S3 method for class 'ellipsoid'
jelinek(x)

# S3 method for class 'ortensor'
jelinek(x)

flinn(x)

# Default S3 method
flinn(x)

# S3 method for class 'ortensor'
flinn(x)

# S3 method for class 'ellipsoid'
flinn(x)

size_invariant(x)

# Default S3 method
size_invariant(x)

# S3 method for class 'ortensor'
size_invariant(x)

# S3 method for class 'ellipsoid'
size_invariant(x)

strain_invariant(x)

# Default S3 method
strain_invariant(x)

# S3 method for class 'ortensor'
strain_invariant(x)

# S3 method for class 'ellipsoid'
strain_invariant(x)

shape_invariant(x)

# Default S3 method
shape_invariant(x)

# S3 method for class 'ortensor'
shape_invariant(x)

# S3 method for class 'ellipsoid'
shape_invariant(x)

kind(x)

# Default S3 method
kind(x)

# S3 method for class 'ortensor'
kind(x)

# S3 method for class 'ellipsoid'
kind(x)
```

## Arguments

- x:

  numeric. either a 3-element vector giving the ellipsoid's semi-axis
  lengths (in any order), an object of class `"ellipsoid"`, or an object
  of class `"ortensor"`.

## Value

positive numeric

## Details

\$\$e_i = \log s_i\$\$ with \\s1 \geq s2 \geq s3\\ the semi-axis lengths
of the ellipsoid.

Lode's shape parameter: \$\$\nu = \frac{2 e_2 - e_1 - e_3}{e_1-e_3}\$\$
with \\e_1 \geq e_2 \geq e_3\\. Note that \\\nu\\ is undefined for
spheres, but we arbitrarily declare \\\nu=0\\ for them. Otherwise \\-1
\geq \nu \geq 1\\. \\\nu=-1\\ for prolate spheroids and \\\nu=1\\ for
oblate spheroids.

Octahedral shear strain \\e_s\\ (Nadai 1963): \$\$e_s =
\sqrt{\frac{(e_1 - e_2)^2 + (e_2 - e_3)^2 + (e_1 - e_3)^2 }{3}}\$\$

Strain symmetry (Flinn 1963): \$\$k = \frac{s_1/s_2 - 1}{s_2/s_3 -
1}\$\$

and strain intensity (Flinn 1963): \$\$d = \sqrt{(s_1/s_2 - 1)^2 +
(s_2/s_3 - 1)^2}\$\$

Jelinek (1981)'s \\P_j\\ parameter: \$\$P_j = e^{\sqrt{2 \vec{v}\cdot
\vec{v}}}\$\$ with \\\vec{v} = e_i - \frac{\sum e_i}{3}\\

## References

Flinn, Derek.(1963): "On the statistical analysis of fabric diagrams."
Geological Journal 3.2: 247-253.

Lode, Walter (1926): "Versuche über den Einfluß der mittleren
Hauptspannung auf das Fließen der Metalle Eisen, Kupfer und Nickel“
(*"Experiments on the influence of the mean principal stress on the flow
of the metals iron, copper and nickel"*\], Zeitschrift für Physik, vol.
36 (November), pp. 913–939,
[doi:10.1007/BF01400222](https://doi.org/10.1007/BF01400222)

Nadai, A., and Hodge, P. G., Jr. (1963): "Theory of Flow and Fracture of
Solids, vol. II." ASME. J. Appl. Mech. December 1963; 30(4): 640.
[doi:10.1115/1.3636654](https://doi.org/10.1115/1.3636654)

Jelinek, Vit. "Characterization of the magnetic fabric of rocks."
Tectonophysics 79.3-4 (1981): T63-T67.

## See also

[`shape_params()`](https://tobiste.github.io/structr/reference/strain_shape.md),
[`ot_eigen()`](https://tobiste.github.io/structr/reference/ot_eigen.md)

## Examples

``` r
# Generate some random data
set.seed(20250411)
dat <- rvmf(100, k = 20)
s <- principal_stretch(dat)

# Volume of ellipsoid
volume(s)
#> [1] 0.171276

#  Size-related tensor invariant of ellipsoids
size_invariant(s)
#> [1] -3.196891

# Strain-related tensor invariant of ellipsoids
strain_invariant(s)
#> [1] 2.60509

# Shape-related tensor invariant of ellipsoids
shape_invariant(s)
#> [1] -0.1155451

# Lode's shape parameter
lode(s)
#> [1] -0.63562

# Nadai's octahedral shear strain
nadai(s)
#> [1] 1.266186

# Jelinek Pj parameter
jelinek(s)
#> [1] 5.993392

# Flinn's intensity and symmetry parameters
flinn(s)
#> $k
#> [1] 8.243876
#> 
#> $d
#> [1] 2.975821
#> 

kind(s)
#> [1] "LLS"
```
