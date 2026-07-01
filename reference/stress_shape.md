# Shape parameters of the stress ellipsoid

Calculates the stress shape (or stress ratio) after Gephart & Forsyth
(1984), Angelier (1979), and Bott (1959) from a given stress tensor. The
parameter represents the specific shape of the stress ellipsoid, which
characterizes a stress state. The Generalized stress shape ratio
(\\A\_\phi\\) after Simpson (1997) distinguishes between the Andersonian
fault regimes (normal, strike-slip, reverse) based on the number of
principal stresses larger than the vertical stress.

## Usage

``` r
stress_shape(tau)
```

## Arguments

- tau:

  symmetric 3x3 matrix. The (reduced) stress tensor.

## Value

list. Its components are the three stress shape parameters `R` (after
Gephart & Forsyth, 1984), `phi` (Angelier, 1979), `bott` (Bott, 1959),
`A_phi` (Simpson, 1997), and the fault regime `type` (N (normal), S
(strike-slip), T (reverse/thrust)).

## Details

Stress shape ratio (\\\Phi\\) after Angelier (1979):

\$\$\Phi = (\sigma_2 - \sigma_3)/(\sigma_1 - \sigma_3)\$\$ Values range
between 0 (\\\sigma_2 = \sigma_3\\) and 1 (\\\sigma_2 = \sigma_1\\). For
\\\Phi = 0\\, the stress ellipsoids takes a prolate geometry (i.e.
\\\sigma_1 \> \sigma_2 = \sigma_3\\), and for \\\Phi = 1\\, it takes a
oblate geometry (\\\sigma_1 = \sigma_2\\). Intermediate shapes for \\0
\< \Phi \< 1\\, i.e. \\\sigma_1 \> \sigma_2 \> \sigma_3\\.

Stress shape (\\R\\ or \\\phi\\) ratio after Gephart & Forsyth (1984):
\$\$R = (\sigma_1 - \sigma_2)/(\sigma_1 - \sigma_3)\$\$ Values ranging
from 0 to 1, with 0 being \\\sigma_1 = \sigma_2\\ and 1 being \\\sigma_2
= \sigma_3\\.

Stress shape ratio (\\R\\) after Bott (1959): \$\$\R = (\sigma_3 -
\sigma_1)/(\sigma_2 - \sigma_1)\$\$ Values range between \\-\infty\\ and
\\+\infty\\.

Generalized stress shape ratio (\\A\_\phi\\) after Angelier (Simpson
1997): \$\$A\_\phi = (n + 0.5) + (-1)^n \* (\phi - 0.5)\$\$ where \\n\\
equals the number of principal stresses larger than the vertical stress.
Values range from 0 to 1 for normal, 1 to 2 for strike-slip, and 2 to 3
for reverse faults.

## References

Angelier, J., 1979. Determination of the mean principal directions of
stresses for a given fault population. Tectonophysics 56, T17–T26.

Bott, M.H.P., 1959. The mechanics of oblique slip faulting. Geol. Mag.
96, 109–117.

Gephart, J.W., Forsyth, D.W., 1984. An improved method for determining
the regional stress tensor using earthquake focal mechanism data:
application to the San Fernando earthquake sequence. J. Geophys. Res.
Solid Earth 89, 9305–9320.

Simpson, R. W. (1997). Quantifying Anderson’s fault types. Journal of
Geophysical Research: Solid Earth, 102(B8), 17909–17919.
https://doi.org/10.1029/97JB01274

## See also

[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md)

Other stress-tensor:
[`fault_instability_criterion()`](https://tobiste.github.io/structr/reference/fault_instability_criterion.md),
[`reduced_stress()`](https://tobiste.github.io/structr/reference/reduced_stress.md),
[`tau-comp`](https://tobiste.github.io/structr/reference/tau-comp.md),
[`tau2rup()`](https://tobiste.github.io/structr/reference/tau2rup.md),
[`tau2stress()`](https://tobiste.github.io/structr/reference/tau2stress.md)

## Examples

``` r
f <- angelier1990$TYM
tau <- reduced_stress(f)
stress_shape(tau)
#> $R
#> [1] 0.898753
#> 
#> $phi
#> [1] 0.101247
#> 
#> $bott
#> [1] 1.112653
#> 
#> $A_phi
#> [1] 0.101247
#> 
#> $type
#> [1] "N"
#> 
```
