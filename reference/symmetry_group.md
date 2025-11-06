# Symmetry groups

Frequently orientations are subject to some symmetry group
\\\mathbb{G}\\, which is a finite set of rotations (satisfying certain
properties). For any rotation \\Q\\ in \\\mathbb{G}\\, the rotation
matrix \\Q \cdot R\\ represents the same orientation as \\R\\ does.

## Usage

``` r
symmetry_group(
  group = c("triclinic", "ray_in_plane", "line_in_plane", "monoclinic",
    "trigonal-trapezohedral", "hexagonal-trapezohedral", "trivial")
)
```

## Arguments

- group:

  character (symmetry class) or integer number of the enantiomorphic
  point group. See table below for details. Also accepts
  `""ray_in_plane"` (equivalent to triclinic symmetry) and
  `"line_in_plane"` (monoclinic).

## Value

list of rotation parameters

## Details

|                |                            |                                                                             |
|----------------|----------------------------|-----------------------------------------------------------------------------|
| Symmetry       | Enantiomorphic Point Group | Example                                                                     |
| `triclinic`    | 1                          | Ray in plane, plagioclase                                                   |
| `monoclinic`   | 2                          | Line in plane, orthoclase, gypsum, muscovite, clinopyroxene, clinoamphibole |
| `orthorhombic` | 222                        | olivine, aragonite, marcasite, orthopyroxenes                               |
| `tetragonal`   | 4                          | Pyramidal: zircon                                                           |
|                | 422                        | Trapezohedral                                                               |
| `trigonal`     | 3                          | Pyramidal, Rhombohedral                                                     |
|                | 32                         | Trapezohedral: \\\alpha\\-Quartz                                            |
| `hexagonal`    | 6                          | Pyramidal, Rhombohedral                                                     |
|                | 622                        | Trapezohedral: \\\beta\\-Quartz                                             |
| `cubic`        | 23                         | Tetartoidal                                                                 |
|                | 432                        | Hexoctahedral: galena, pyrite, fluorite                                     |

## Examples

``` r
symmetry_group("triclinic")
#> [[1]]
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> 
#> [[2]]
#>      [,1] [,2] [,3]
#> [1,]   -1    0    0
#> [2,]    0    1    0
#> [3,]    0    0   -1
#> 

symmetry_group(2)
#> [[1]]
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> 
#> [[2]]
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0   -1    0
#> [3,]    0    0   -1
#> 
#> [[3]]
#>      [,1] [,2] [,3]
#> [1,]   -1    0    0
#> [2,]    0    1    0
#> [3,]    0    0   -1
#> 
#> [[4]]
#>      [,1] [,2] [,3]
#> [1,]   -1    0    0
#> [2,]    0   -1    0
#> [3,]    0    0    1
#> 
```
