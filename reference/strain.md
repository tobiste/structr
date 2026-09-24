# Strain quantities

Converts different quantifications of strain.

## Usage

``` r
longitudinal_strain(l, l0)

angular_strain(psi, degree = TRUE)

quadratic_elongation(e = NULL, l = NULL, l0 = NULL)

stretch(e = NULL, l = NULL, l0 = NULL, lambda = NULL)
```

## Arguments

- l, l0:

  numeric. Final and original length, respectively

- psi:

  numeric. Angle

- degree:

  logical. Whether `psi` is given in degree (the default) or radians

- e:

  numeric. Elongation

- lambda:

  numeric. Quadratic elongation

## Details

*Longitudianal strain* (also **elongation**) is the change in length
divided by the original length: \$\$e = (l-l_o)/l_o\$\$

*Angular strain* is the change in angle between two lines that were
initially perpendicular" \$\$\gamma = \tan \psi\$\$

*Stretch* is the ratio of the final length and the original length:
\$\$s = \sqrt{\lambda} = l/l_o = 1 + e\$\$

*Quadratic elongation* is the quadratic stretch: \$\$\lambda = s^2 =
(l/l_o)^2 = (1+e)^2\$\$

## Examples

``` r
longitudinal_strain(l=10,l0=5) 
#> [1] 1
angular_strain(25)
#> [1] 0.4663077
quadratic_elongation(l=10,l0=5)
#> [1] 4
stretch(l=10,l0=5)
#> [1] 2
```
