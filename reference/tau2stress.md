# Principal Stresses from Stress Tensor

The eigenvector and eigenvalues of a stress tensor give the orientations
and relative magnitudes of the principal stress axes.

## Usage

``` r
tau2stress(tau)
```

## Arguments

- tau:

  symmetric 3x3 matrix. The (reduced) stress tensor.

## Value

list with the following components

- `"sigma_vals"`:

  numeric. The relative magnitudes of the principal stress axes.

- `"principal_axes"`:

  The principal stress axes as `Line` objects.

## See also

[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md)

Other stress-tensor:
[`fault_instability_criterion()`](https://tobiste.github.io/structr/reference/fault_instability_criterion.md),
[`reduced_stress()`](https://tobiste.github.io/structr/reference/reduced_stress.md),
[`stress_shape()`](https://tobiste.github.io/structr/reference/stress_shape.md),
[`tau-comp`](https://tobiste.github.io/structr/reference/tau-comp.md),
[`tau2rup()`](https://tobiste.github.io/structr/reference/tau2rup.md)

## Examples

``` r
f <- angelier1990$TYM
tau <- reduced_stress(f)
tau2stress(tau)  
#> $sigma_vals
#>     sigma1     sigma2     sigma3 
#>  1.5186760 -0.6378677 -0.8808083 
#> 
#> $principal_axes
#> Ray object (n = 3):
#>          azimuth    plunge
#> sigma1 225.71104 84.219427
#> sigma2  62.26628  5.542466
#> sigma3 332.10740  1.636857
#> 
```
