# Fault Instability Criterion

Calculates the instability criterion \\I\\ after Vavrycuk (2013, Eq. 3).
Instability ranges from 0 (most stable) to 1 (most unstable) The most
unstable fault is the optimally oriented fault for shear faulting.

## Usage

``` r
fault_instability_criterion(fault, R, friction = 0.6)
```

## Arguments

- fault:

  `"Fault"` object where the rows are the observations, and the columns
  the coordinates.

- R:

  numeric. Stress ratio after Gephart and Forsyth (1984): \\(\sigma_1 -
  \sigma_2)/(\sigma_1 - \sigma_3)\\

- friction:

  numeric. Coefficient of friction (0.6 by default)

## Value

numeric. Instability ranges from 0 (most stable) to 1 (most unstable).
The most unstable fault is the optimally oriented fault for shear
faulting.

## Details

\$\$I = \frac{\tau - \mu(\sigma - \sigma_1)}{\tau_c - \mu(\sigma_c -
\sigma_1)}\$\$

where \\\tau_c\\ and \\\sigma_c\\ are the shear traction and effective
normal traction along the optimally oriented fault, and \\\tau\\ and
\\\sigma\\ are the shear traction and effective normal traction along
the analysed fault plane.

## References

Vavrycuk, V., Bouchaala, F. & Fischer, T., 2013. High-resolution fault
image from accurate locations and focal mechanisms of the 2008 swarm
earthquakes in West Bohemia, Czech Republic, Tectonophysics, 590,
189–195.

## See also

Other stress-tensor:
[`reduced_stress()`](https://tobiste.github.io/structr/reference/reduced_stress.md),
[`stress_shape()`](https://tobiste.github.io/structr/reference/stress_shape.md),
[`tau-comp`](https://tobiste.github.io/structr/reference/tau-comp.md),
[`tau2rup()`](https://tobiste.github.io/structr/reference/tau2rup.md),
[`tau2stress()`](https://tobiste.github.io/structr/reference/tau2stress.md)

## Examples

``` r
f <- angelier1990$TYM
tau <- reduced_stress(f)
s <- stress_shape(tau)

fault_instability_criterion(f, s$R)
#>  [1] 0.80583583 0.96786690 0.05589429 0.57145143 0.76778700 0.15189430
#>  [7] 0.13812261 0.81158280 0.45559237 0.66081621 0.58543179 1.13776326
#> [13] 0.87881003 0.41970424 0.63313399 0.96055789 0.97720506 0.88549935
#> [19] 0.28601284 0.92598051 1.39315476 1.04755339 0.30091177 1.38021623
#> [25] 1.04697429 0.70405761 0.64255950 0.59863018 0.46312371 1.08484041
#> [31] 0.90212159 0.82693097 1.00838327 0.87372833 0.34157203 1.09388247
#> [37] 0.22224539 0.19233285
```
