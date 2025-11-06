# Stress components

Calculates some stress components

## Usage

``` r
diff_stress(sigma1, sigma3)

mean_stress(sigma1, sigma3)

shear_stress(sigma1, sigma3, theta)

normal_stress(sigma1, sigma3, theta)

fracture_angle(mu)

slip_tendency(sigma_s, sigma_n)

dilatation_tendency(sigma1, sigma3, sigma_n)
```

## Arguments

- sigma1, sigma3:

  numeric. Magnitudes of maximum and minimum principal stress
  (\\\sigma_1\\ and \\\sigma_3\\), respectively.

- theta:

  numeric. Angle \\\theta\\ between fracture and \\\sigma_1\\.

- mu:

  numeric. Coefficient of internal friction \\\mu\\.

- sigma_s, sigma_n:

  numeric. Magnitudes of shear and normal stress (\\\sigma_s\\ and
  \\\sigma_n\\), respectively.

## Value

numeric

## Examples

``` r
s1 <- 1025
s3 <- 250

diff_stress(s1, s3)
#> [1] 775
mean_stress(s1, s3)
#> [1] 637.5
ss <- shear_stress(s1, s3, theta = 35)
print(ss)
#> [1] 364.1309
sn <- normal_stress(s1, s3, theta = 35)
print(sn)
#> [1] 504.9672
fracture_angle(mu = 0.6)
#> [1] 60.48188

slip_tendency(ss, sn)
#> [1] 0.7210981
dilatation_tendency(s1, s3, sn)
#> [1] 0.6710101
```
