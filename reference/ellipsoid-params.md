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

The natural strain is \$\$\bar{\epsilon}\_i = \log s_i = \log(1 +
\epsilon_i)\$\$ with \\s_1 \geq s_2 \geq s_3\\ the semi-axis lengths of
the ellipsoid, or \\\epsilon_i\\ the strains (elongation) given by
\\\epsilon = \frac{l-l_0}{l_0}\\ (hence \\s = \frac{l}{l_0}\\).

Lode's parameter for **strain symmetry**: \$\$\nu = \frac{2
\bar{\epsilon}\_2 - \bar{\epsilon}\_1 -
\bar{\epsilon}\_3}{\bar{\epsilon}\_1-\bar{\epsilon}\_3}\$\$ with
\\\bar{\epsilon}\_1 \geq \bar{\epsilon}\_2 \geq \bar{\epsilon}\_3\\.
Note that \\\nu\\ is undefined for spheres, but we arbitrarily declare
\\\nu=0\\ for them (plane strain). Otherwise \\-1 \geq \nu \geq 1\\.
\\\nu=-1\\ for prolate spheroids (constriction) and \\\nu=1\\ for oblate
spheroids (flattening).

**Octahedral shear strain** \\\bar{\epsilon}\_s\\ (Nádai 1963):
\$\$\bar{\epsilon}\_s = \sqrt{\frac{(\bar{\epsilon}\_1 -
\bar{\epsilon}\_2)^2 + (\bar{\epsilon}\_2 - \bar{\epsilon}\_3)^2 +
(\bar{\epsilon}\_1 - \bar{\epsilon}\_3)^2 }{3}}\$\$ This is the amount
of strain assuming coaxial deformation (pure-shear).

**Strain symmetry** (Flinn 1963): \$\$k = \frac{s_1/s_2 - 1}{s_2/s_3 -
1}\$\$ The value ranges from 0 to \\\infty\\, and is 0 for oblate
ellipsoids (flattening), 1 for plane strain and \\\infty\\ for prolate
ellipsoids (constriction).

and **strain intensity** (Flinn 1963): \$\$d = \sqrt{(s_1/s_2 - 1)^2 +
(s_2/s_3 - 1)^2}\$\$ This is analogous to Nadai's strain parameter.

Jelinek (1981)'s \\P_j\\ parameter: \$\$P_j = \bar{\epsilon}^{\sqrt{2
\vec{v}\cdot \vec{v}}}\$\$ with \\\vec{v} = \bar{\epsilon}\_i -
\frac{\sum \bar{\epsilon}\_i}{3}\\

## References

Flinn, Derek.(1963): "On the statistical analysis of fabric diagrams."
Geological Journal 3.2: 247-253.

Lode, Walter (1926): "Versuche über den Einfluß der mittleren
Hauptspannung auf das Fließen der Metalle Eisen, Kupfer und Nickel“
(*"Experiments on the influence of the mean principal stress on the flow
of the metals iron, copper and nickel"*\], Zeitschrift für Physik, vol.
36 (November), pp. 913–939,
[doi:10.1007/BF01400222](https://doi.org/10.1007/BF01400222)

Nádai, A., and Hodge, P. G., Jr. (1963): "Theory of Flow and Fracture of
Solids, vol. II." ASME. J. Appl. Mech. December 1963; 30(4): 640.
[doi:10.1115/1.3636654](https://doi.org/10.1115/1.3636654)

Jelinek, Vit. "Characterization of the magnetic fabric of rocks."
Tectonophysics 79.3-4 (1981): T63-T67.

## See also

[`shape_params()`](https://tobiste.github.io/structr/reference/strain_shape.md),
[`ot_eigen()`](https://tobiste.github.io/structr/reference/ot_eigen.md),
[`hsu_plot()`](https://tobiste.github.io/structr/reference/hsu_plot.md),
[`flinn_plot()`](https://tobiste.github.io/structr/reference/flinn_plot.md)

Other ellipsoid:
[`ellipsoid-class`](https://tobiste.github.io/structr/reference/ellipsoid-class.md),
[`ellipsoid_from_stretch()`](https://tobiste.github.io/structr/reference/ellipsoid_from_stretch.md)

## Examples

``` r
# Generate some random orientation data
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

# Lode's shape parameter for the strain symmetry ratio
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

# Ellipsoid data
hossack_ell <- lapply(seq.int(nrow(hossack1968)), function(i) {
  ellipsoid_from_stretch(hossack1968[i, 3], hossack1968[i, 2], hossack1968[i, 1])
})
sapply(hossack_ell, nadai)
#>  [1] 2.197985 2.462677 2.507189 1.986605 2.824592 2.421925 2.580754 2.712021
#>  [9] 2.573880 1.877145 2.556829 2.518613 2.236213 2.643440 2.652371 2.383755
#> [17] 2.361141 1.785130 2.169690 2.349450 2.362748 2.193710 2.382355 2.139940
#> [25] 1.970644 2.192330 1.945885 2.334462 1.978659 2.161851 2.792102 2.188473
#> [33] 1.895198 2.071682 2.346818 2.032335 2.653843 2.943093 2.147554 2.295867
#> [41] 2.638900 2.342928 2.626458 2.172841 2.280597 1.971667 2.457973 2.248606
#> [49] 2.032136 2.358794 2.152572 1.537893 1.142119 1.258516 1.403769 1.312832
#> [57] 1.249068 1.466010 1.599231 2.079253 1.380805 1.420466 1.690162 1.497369
#> [65] 1.443027 1.632871 1.081520 1.313677 1.499537 1.702833 1.466010 1.602093
#> [73] 1.523021 1.542330 1.814373 1.793242 1.847809 1.898982 1.969607 1.987016
#> [81] 1.680192
```
