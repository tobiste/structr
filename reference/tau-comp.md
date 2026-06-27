# Resolved Shear and Normal Stress

[`tau2stress()`](https://tobiste.github.io/structr/reference/tau2stress.md)
calculate normal and shear stress components, while `tau2tendency()`
computes the tendency for slip and dilatency for a given set of faults
and a given strerss tensor.

## Usage

``` r
tau2shearnorm(tau, fault, friction = 0.6)

tau2tendency(tau, fault, friction = 0.6)
```

## Arguments

- tau:

  symmetric 3x3 matrix. The (reduced) stress tensor.

- fault:

  `"Fault"` object where the rows are the observations, and the columns
  the coordinates.

- friction:

  numeric. Coefficient of friction (0.6 by default)

## Value

2-column numeric array giving the relative normal and shear stress
components for each fault in `fault`.

## See also

[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md)

Other stress-tensor:
[`fault_instability_criterion()`](https://tobiste.github.io/structr/reference/fault_instability_criterion.md),
[`reduced_stress()`](https://tobiste.github.io/structr/reference/reduced_stress.md),
[`stress_shape()`](https://tobiste.github.io/structr/reference/stress_shape.md),
[`tau2rup()`](https://tobiste.github.io/structr/reference/tau2rup.md),
[`tau2stress()`](https://tobiste.github.io/structr/reference/tau2stress.md)

## Examples

``` r
f <- angelier1990$TYM
tau <- reduced_stress(f)
tau2shearnorm(tau, f)  
#>             normal      shear
#>  [1,] -0.319049467 -1.0074373
#>  [2,] -0.199966082 -1.0595722
#>  [3,] -0.683886779  0.5641852
#>  [4,] -0.516062452  0.8527739
#>  [5,] -0.368199436  0.9627699
#>  [6,] -0.703458404  0.3081223
#>  [7,] -0.712066093  0.5942078
#>  [8,] -0.288765979 -1.0276038
#>  [9,] -0.194934248  1.0740371
#> [10,] -0.218113823 -1.0672641
#> [11,] -0.611589412 -0.7495115
#> [12,]  0.103272184  1.1484876
#> [13,]  0.002922413 -1.1571691
#> [14,] -0.440070773 -0.8983697
#> [15,] -0.311052993  1.0178778
#> [16,] -0.311065410  0.9719439
#> [17,] -0.224411787 -1.0421138
#> [18,] -0.345461213 -0.9740412
#> [19,] -0.702405649 -0.6191082
#> [20,] -0.039455163 -1.1427998
#> [21,]  0.189119716  1.1009707
#> [22,]  0.365529916  1.1962668
#> [23,] -0.652338281 -0.4231349
#> [24,]  0.405893693  1.1501004
#> [25,] -0.132550302  1.0693392
#> [26,] -0.716263098 -0.5471174
#> [27,] -0.368438977 -0.9833086
#> [28,] -0.121031192  1.1152451
#> [29,] -0.567179606 -0.8080235
#> [30,] -0.115269578 -1.0803975
#> [31,] -0.490745306 -0.8277931
#> [32,] -0.366994676 -0.9677083
#> [33,] -0.419461649 -0.8666342
#> [34,]  0.287785173  1.1988409
#> [35,] -0.504018233  0.8705438
#> [36,] -0.506658988  0.7120882
#> [37,] -0.413885186  0.8822828
#> [38,] -0.480009299  0.8488841

tau2tendency(tau, f)  
#>       slip_tendency dilatation_tendency
#>  [1,]     3.1576209           0.7658835
#>  [2,]     5.2987594           0.7162548
#>  [3,]    -0.8249688           0.9179317
#>  [4,]    -1.6524626           0.8479899
#>  [5,]    -2.6148054           0.7863671
#>  [6,]    -0.4380106           0.9260883
#>  [7,]    -0.8344841           0.9296756
#>  [8,]     3.5586043           0.7532627
#>  [9,]    -5.5097402           0.7141577
#> [10,]     4.8931521           0.7238180
#> [11,]     1.2255142           0.8878014
#> [12,]    11.1209776           0.5898783
#> [13,]  -395.9636215           0.6316997
#> [14,]     2.0414208           0.8163199
#> [15,]    -3.2723613           0.7625509
#> [16,]    -3.1245642           0.7625561
#> [17,]     4.6437570           0.7264427
#> [18,]     2.8195385           0.7768908
#> [19,]     0.8814111           0.9256496
#> [20,]    28.9645181           0.6493609
#> [21,]     5.8215543           0.5541008
#> [22,]     3.2726920           0.4805808
#> [23,]     0.6486434           0.9047837
#> [24,]     2.8335014           0.4637589
#> [25,]    -8.0674216           0.6881588
#> [26,]     0.7638497           0.9314248
#> [27,]     2.6688506           0.7864669
#> [28,]    -9.2145260           0.6833582
#> [29,]     1.4246342           0.8692933
#> [30,]     9.3727894           0.6809570
#> [31,]     1.6868081           0.8374388
#> [32,]     2.6368457           0.7858650
#> [33,]     2.0660629           0.8077309
#> [34,]     4.1657495           0.5129814
#> [35,]    -1.7272070           0.8429704
#> [36,]    -1.4054585           0.8440710
#> [37,]    -2.1317091           0.8054069
#> [38,]    -1.7684743           0.8329645
```
