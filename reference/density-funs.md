# Density estimation

The methods of density estimation in a stereonet.

## Usage

``` r
exponential_kamb(cos_dist, sigma = 3)

linear_inverse_kamb(cos_dist, sigma = 3)

square_inverse_kamb(cos_dist, sigma = 3)

kamb_count(cos_dist, sigma = 3)

schmidt_count(cos_dist, sigma = NULL)
```

## Arguments

- cos_dist:

  cosine distances

- sigma:

  (optional) numeric. The number of standard deviations defining the
  expected number of standard deviations by which a random sample from a
  uniform distribution of points would be expected to vary from being
  evenly distributed across the hemisphere. This controls the size of
  the counting circle, and therefore the degree of smoothing. Higher
  sigmas will lead to more smoothing of the resulting density
  distribution. This parameter only applies to Kamb-based methods.
  Defaults to 3.

## Value

list

## Details

`exponential_kamb()`: Kamb with exponential smoothing A modified Kamb
method using exponential smoothing (ref1). Units are in numbers of
standard deviations by which the density estimate differs from uniform.

`linear_inverse_kamb()`: Kamb with linear smoothing A modified Kamb
method using linear smoothing (ref1). Units are in numbers of standard
deviations by which the density estimate differs from uniform.

`square_inverse_kamb()`: Kamb with squared smoothing A modified Kamb
method using squared smoothing (ref1). Units are in numbers of standard
deviations by which the density estimate differs from uniform.

`kamb_count()`: Kamb with no smoothing Kamb's method (ref2) with no
smoothing. Units are in numbers of standard deviations by which the
density estimate differs from uniform.

`schmidt_count()`: 1% counts. The traditional "Schmidt" (a.k.a. 1%)
method. Counts points within a counting circle comprising 1% of the
total area of the hemisphere. Does not take into account sample size.
Units are in points per 1% area.

## References

Vollmer, 1995. C Program for Automatic Contouring of Spherical
Orientation Data Using a Modified Kamb Method. Computers & Geosciences,
Vol. 21, No. 1, pp. 31–49.

Kamb, 1959. Ice Petrofabric Observations from Blue Glacier, Washington,
in Relation to Theory and Experiment. Journal of Geophysical Research,
Vol. 64, No. 11, pp. 1891–1909.
