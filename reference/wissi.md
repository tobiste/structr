# Weighted Iterative Sigma-Space Inversion (WISSI)

Combines the four classic fault slip inversion algorithms into a single
coherent framework operating in the Yamaji-Sato 5-sphere sigma-space.

## Usage

``` r
wissi(
  normals,
  slips,
  weights = NULL,
  sigma_alpha_deg = 10,
  gamma_max = 10,
  n_anneal = 8L,
  max_iter = 50L,
  tol_deg = 1e-04,
  run_stage4 = TRUE
)
```

## Arguments

- normals:

  N x 3 matrix of unit fault plane normals.

- slips:

  N x 3 matrix of unit slip direction vectors (parallel to shear
  traction; direction of hanging-wall motion).

- weights:

  Optional length-N weight vector. NAs imputed with observed mean.
  Normalised internally (mean = 1). Use signal_to_weights() and
  combine_weights() to build from field ranks, measurement errors, and
  prior RUP values.

- sigma_alpha_deg:

  Estimated slip direction measurement error in degrees. Used in Stage 4
  analytic uncertainty. Default 10.

- gamma_max:

  Maximum sense annealing sharpness parameter. Higher values commit more
  strongly to the predicted slip sense. Default 10. Set to 0 to disable
  sense annealing (fully sense-agnostic, like Hansen 2013).

- n_anneal:

  Number of annealing steps (outer loop of Stage 3). Default 8.

- max_iter:

  Maximum inner iterations per annealing step. Default 50.

- tol_deg:

  Convergence tolerance in angular stress distance (degrees). Default
  1e-4.

- run_stage4:

  Logical. Compute analytic uncertainty (Stage 4). Default TRUE.

## Value

A named list with: sigma : 3x3 reduced stress tensor (Cartesian frame) y
: 6D unit y-vector on S^5 sigma1/2/3 : unit vectors of principal stress
axes (max to min) eigenvalues : eigenvalues of sigma (decreasing) Phi :
shape ratio (sigma1-sigma2)/(sigma1-sigma3) in \\\[0,1\]\\ alpha_deg :
per-fault angular misfit (unsigned, 0-90 deg) alpha_signed_deg :
per-fault signed misfit (0-180 deg) mean_alpha : mean angular misfit
across all faults (deg) suspected_flipped : row indices where
alpha_signed \> 90 deg n_flipped_sense : number of faults whose sense
was corrected in Stage 3 slips_corrected : sense-corrected slip matrix
used in final inversion mu : per-fault magnitude weights from Stage 2/3
phi_sense : per-fault tanh sense confidence from Stage 3 eigenvalue_gap
: lambda_2 - lambda_1 of M5 (condition number proxy) M5_eigvals : all 5
eigenvalues of final M5 unc : Stage 4 uncertainty list (if run_stage4 =
TRUE): Cov5 : 5x5 covariance matrix in sigma-space Cov_y6 : 6x6
covariance matrix (y-space) eigval_gap : eigenvalue gap (same as above)
cov_eigvals : eigenvalues of Cov5 sigma1_unc_deg : approx 1-sigma
uncertainty on sigma1 orientation Phi_unc : approx 1-sigma uncertainty
on Phi n_iter_total : total number of inner iterations
