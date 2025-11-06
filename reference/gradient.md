# Velocity gradient gradient tensors

The velocity gradient tensor describes the velocity of particles at any
instant during the deformation. Velocity gradient tensor from
deformation gradient tensor.

## Usage

``` r
is.velgrad(x)

as.velgrad(object)

velgrad(x, time, ...)

# Default S3 method
velgrad(x, ...)

# S3 method for class 'defgrad'
velgrad(x, time = 1, ...)
```

## Arguments

- x:

  3x3 matrix. Deformation gradient tensor.

- object:

  3x3 `"matrix"`

- time:

  numeric. Total time (default is 1)

- ...:

  parameters passed to function call

## Value

3x3 matrix. If steps is \> 1, then a list of matrices is returned.

## Details

`velgrad()` calculates the velocity gradient tensor as the matrix
logarithm of the deformation gradient tensor divided by given time, and
the deformation gradient tensor accumulated after some time.

## See also

[`defgrad()`](https://tobiste.github.io/structr/reference/defgrad.md),
[`stereo_path()`](https://tobiste.github.io/structr/reference/stereo_path.md)
for plotting

## Examples

``` r
d <- defgrad_from_generalshear(k = 2.5, gamma = 0.9)
l <- velgrad(d, time = 10)
d_steps <- defgrad(l, time = 10, steps = 2)

# apply on orientation data
set.seed(20250411)
v <- rvmf(100, mu = Line(0, 90), k = 100)
v_trans <- lapply(d_steps, function(i) {
  transform_linear(v, i)
})

# plot in stereonet
axes <- Vec3(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
stereo_path(v_trans, type = "l", add = FALSE)
stereo_path(v_trans, type = "p", col = assign_col(seq_along(v_trans)), pch = 16, cex = .4)
points(axes, pch = 15)
text(axes, labels = c("x", "y", "z"), pos = 1)
```
