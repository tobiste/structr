# Fault displacement components

Calculate fault displacement from (at least three) components of of the
fault displacement system

## Usage

``` r
fault_displacements(
  dip = NULL,
  delta = NULL,
  rake = NULL,
  verticalthrow = NULL,
  dipslip = NULL,
  heave = NULL,
  netslip = NULL,
  horizontalthrow = NULL,
  strikeslip = NULL
)
```

## Arguments

- dip:

  fault's dip angle (in degrees)

- delta:

  angle between horizontal displacement vector and fault's strike (in
  degrees)

- rake:

  (in degrees)

- verticalthrow:

  vertical throw

- dipslip:

  dip slip component

- heave:

  apparent horizontal offset perpendicular to strike

- netslip:

  offset on fault plane parallel to the fault's motion direction

- horizontalthrow:

  apparent horizontal offset parallel to fault's motion direction

- strikeslip:

  strike-slip component

## Value

array

## Details

see vignette for description of fault displacement components

## Examples

``` r
if (FALSE) { # \dontrun{
fault_displacements(strikeslip = 2, verticalthrow = -5, heave = 3)
} # }
```
