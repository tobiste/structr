# Random Rotation Matrices

Random sample of matrices in SO(3).

## Usage

``` r
rrot(n)
```

## Arguments

- n:

  integer. number of random samples to be generated

## Value

list of rotation matrices

## See also

Other random:
[`rbing`](https://tobiste.github.io/structr/reference/rbing.md),
[`rfb()`](https://tobiste.github.io/structr/reference/rfb.md),
[`rkent()`](https://tobiste.github.io/structr/reference/rkent.md),
[`runif.spherical()`](https://tobiste.github.io/structr/reference/runif.spherical.md),
[`rwatson()`](https://tobiste.github.io/structr/reference/rwatson.md),
[`vonmises-fisher`](https://tobiste.github.io/structr/reference/vonmises-fisher.md)

## Examples

``` r
set.seed(20250411)
# Generate 10 random SO(3) rotation matrices
r <- rrot(10)

# convert SO(3) matrices to "Pair"
rp <- lapply(r, rot2pair) |> lapply(unclass)
rp <- do.call(rbind, args = rp) |>
  as.Pair()
  
# plot pairs
plot(rp)
```
