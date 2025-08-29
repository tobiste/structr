#' Fault displacement components
#'
#' Calculate fault displacement from (at least three) components of of the
#' fault displacement system
#'
#' @param dip fault's dip angle (in degrees)
#' @param delta angle between horizontal displacement vector and fault's strike (in degrees)
#' @param rake (in degrees)
#' @param heave apparent horizontal offset perpendicular to strike
#' @param horizontalthrow apparent horizontal offset parallel to fault's motion direction
#' @param netslip offset on fault plane parallel to the fault's motion direction
#' @param verticalthrow vertical throw
#' @param dipslip dip slip component
#' @param strikeslip  strike-slip component
#'
#' @details
#' see vignette for description of fault displacement components
#'
#' @return array
#'
#' @examples
#' \dontrun{
#' fault_displacements(strikeslip = 2, verticalthrow = -5, heave = 3)
#' }
fault_displacements <- function(
    dip = NULL, delta = NULL, rake = NULL,
    verticalthrow = NULL, dipslip = NULL, heave = NULL,
    netslip = NULL, horizontalthrow = NULL, strikeslip = NULL) {
  # Calculate based on available parameters
  if (!is.null(dip) && !is.null(verticalthrow) && !is.null(delta)) {
    # Slip components in the vertical plane perpendicular to strike
    dipslip <- verticalthrow / sind(dip)
    heave <- sqrt(dipslip^2 - verticalthrow^2)

    # Slip components in the horizontal plane
    horizontalthrow <- heave / sind(delta)
    strikeslip <- sqrt(horizontalthrow^2 - heave^2)

    # Slip in the fault plane
    netslip <- sqrt(strikeslip^2 + dipslip^2)
    rake <- asind(dipslip / netslip)
  } else if (!is.null(delta) && !is.null(horizontalthrow) && !is.null(dip)) {
    # Slip components in the horizontal plane
    strikeslip <- abs(cosd(delta) * horizontalthrow)
    heave <- sqrt(horizontalthrow^2 - strikeslip^2)

    # Slip components in the vertical plane perpendicular to strike
    dipslip <- heave / cosd(dip)
    verticalthrow <- sqrt(dipslip^2 - heave^2)

    # Slip in the fault plane
    netslip <- sqrt(strikeslip^2 + dipslip^2)
    rake <- atan2d(dipslip, strikeslip)
  } else if (!is.null(rake) && !is.null(verticalthrow) && !is.null(dip)) {
    dipslip <- verticalthrow / sind(dip)
    heave <- sqrt(dipslip^2 - verticalthrow^2)
    strikeslip <- dipslip / tand(rake)
    horizontalthrow <- sqrt(strikeslip^2 + heave^2)
    netslip <- sqrt(strikeslip^2 + dipslip^2)
  } else if (!is.null(strikeslip) && !is.null(verticalthrow) && !is.null(heave)) {
    dipslip <- sqrt(heave^2 + verticalthrow^2)
    dip <- atan2d(verticalthrow, heave)
    horizontalthrow <- sqrt(heave^2 + strikeslip^2)
    delta <- atan2d(heave, strikeslip)
    netslip <- sqrt(strikeslip^2 + dipslip^2)
    rake <- atan2d(dipslip, strikeslip)
  } else if (!is.null(netslip) && !is.null(rake) && !is.null(dip)) {
    rake_rad <- deg2rad(rake)
    strikeslip <- netslip * cos(rake_rad)
    dipslip <- netslip * sin(rake_rad)
    dip_rad <- deg2rad(dip)
    heave <- dipslip * cos(dip_rad)
    verticalthrow <- dipslip * sin(dip_rad)
    delta <- atan2d(heave, strikeslip)
    horizontalthrow <- sqrt(heave^2 + strikeslip^2)
  } else {
    stop("Insufficient or incompatible parameter combinations provided for computation.")
  }

  cbind(
    "dip" = abs(dip),
    "delta" = delta,
    "rake" = rake,
    "verticalthrow" = verticalthrow,
    "horizontalthrow" = horizontalthrow,
    "heave" = heave,
    "dipslip" = dipslip,
    "strikeslip" = strikeslip,
    "netslip" = netslip
  )
}



#' Fault displacement tensor
#'
#' Creates the fault displacement tensor from displacement components.
#' If the dip direction is know, the tensor will be rotated into the
#' geographic reference frame.
#'
#' @param h,s,v numeric. the heave, strike-slip, vertical throw displacement
#' @param dip_direction (optional) dip direction in degrees.
#' @param ftensor Fault displacement tensor. A 3x3 matrix.
#' If `NULL`, the fault tensor will be given in the fault displacement
#' coordinates. Otherwise, the tensor will be in the geographic reference frame.
#'
#' @returns `fault_tensor()` returns a 3x3 matrix of class `"ftensor"` containing the fault displacement tensor.
#' `fault_tensor_analysis()` returns a list containing the principal fault displacement tensor and the fault orientation.
#'
#' @details
#' x axis of tensor = heave, y = strike slip, z = vertical throw (positive for
#' thrusting, negative for normal faulting) is the principal fault displacement tensor.
#' This can be rotated in the fault plane orientation to retrieve slip components and rake.
#'
#' The fault displacement tensor gives the displacements in all directions and has the following properties:
#'
#' - The square root of tensor's trace (i.e. the sum of the diagonal elements) equals the **net slip** on the fault plane.
#'
#' - The determinant of the tensor relates to the volumetric strain by: det(F) - 1, where
#'
#' ([fault_tensor_decomposition()]) retrieves the principal fault displacement tensor using Singular Value Decomposition of a Matrix
#' and the fault orientation if the dip direction is known.
#'
#' The orientation of the net-slip vector is the lineation component of the fault orientation.
#'
#' @name fault-tensor
#'
#' @examples
#' A_princ <- fault_tensor(s = 2, v = -5, h = 3)
#' print(A_princ)
#' det(A_princ)
#'
#' A_geo <- fault_tensor(s = 2, v = -5, h = 3, dip_direction = 45)
#' print(A_geo)
#' det(A_geo)
#'
#' fault_tensor_decomposition(A_geo, dip_direction = 45)
NULL

#' @rdname fault-tensor
#' @export
fault_tensor <- function(h, s, v, dip_direction = NULL) {
  A <- diag(c(h, s, v), 3, 3)
  # A_det <- det(A)

  if (!is.null(dip_direction)) {
    dipdir_rad <- deg2rad(dip_direction)
    # coordinates are measured from E!!!
    # c(1, 0, 0) = E
    # c(0, 1, 0) = N
    A[, 1] <- vrotate(t(A[, 1]), t(A[, 3]), -dipdir_rad) |> as.vector()
    A[, 2] <- vrotate(t(A[, 2]), t(A[, 3]), pi - dipdir_rad) |> as.vector()
  }
  # A <- A/det(A) * A_det

  class(A) <- append(class(A), "ftensor")
  A
}

#' @rdname fault-tensor
#' @export
fault_tensor_decomposition <- function(ftensor, dip_direction = NULL) {
  stopifnot(inherits(ftensor, "ftensor") || is.matrix(ftensor))
  A_svd <- svd(ftensor)
  A_d <- A_svd$d
  A_v <- A_svd$v

  A_pr0 <- A_v * A_d
  A_pr <- diag(c(-A_pr0[2, 3], -A_pr0[3, 1], A_pr0[1, 2]), 3, 3)

  strikeslip <- A_pr[2, 2]
  verticalthrow <- A_pr[3, 3]
  heave <- A_pr[1, 1]

  fd <- fault_displacements(strikeslip = strikeslip, verticalthrow = verticalthrow, heave = heave)

  if (!is.null(dip_direction)) {
    fault_p <- Plane(dip_direction, fd[, "dip"])
    fault <- Fault_from_rake(fault_p, abs(fd[, "rake"]) * sign(verticalthrow))
  } else {
    fault <- NULL
  }

  # A_norm = A/det(A)
  # det(A_norm)

  E <- 0.5 * (t(A_pr) %*% A_pr - diag(A_pr))
  E_norm <- E / det(E)
  volumetric_strain <- sum(diag(E))
  # volumetric_strain2 <- det(A_pr) - 1
  shear_strain <- sqrt(2 * sum(E[lower.tri(E)]^2))
  # gamma_xy <- 2 * E[1, 2]
  # gamma_xz <- 2 * E[1, 3]
  # gamma_yz <- 2 * E[2, 3]

  list(
    displacements = fd,
    fault = fault,
    strain_tensor = E,
    volumetric_strain = volumetric_strain,
    # volumetric_strain2 = volumetric_strain2,
    shear_strain = shear_strain
  )
}


#' #' Title
#' #'
#' #' @param n,l,m direction cosines of the plane referred to the stress axes (<U+03C3>x, <U+03C3>y, <U+03C3>z) system and <U+03B8> is the striation (resolved shear stress) rake on the fault plane
#' #' @param R stress ratio R=(<U+03C3>z<U+2212><U+03C3>x)/(<U+03C3>y<U+2013><U+03C3>x)
#' #'
#' #' @returns  striation (resolved shear stress) rake on the fault plane
#' #' @export
#' #'
#' #' @examples
#' bott <- function(n, l, m, R) {
#'   tan_theta <- (n / (l * m)) * (m^2 - (1 - n^2) * R)
#'   theta <- atan(tan_theta) * 180 / pi
#'   return(theta)
#' }
#'
#' #' Bott's stress ratio
#' #'
#' #' @param lambda angle between y (azimuth of <U+03C3>y) and the azimuth (<U+03B1>) of the
#' #' fault, therefore <U+03BB>=<U+03B1><U+00B1>y (it is added or subtracted according to the sense of movement: dextral or sinistral, respectively) (in degree)
#' #' @param theta  striation (resolved shear stress) rake on the fault plane (in degree)
#' #' @param phi fault dip (in degree)
#' #'
#' #' @references Calvin et al. (2014)
#' #'
#' #' @return numeric
#' #' @export
#' #'
#' #' @examples
#' stress_ratio_simongomez <- function(lambda, theta, phi) {
#'   if (theta == 90 | theta == 180) {
#'     theta - 10^-6
#'   }
#'   if (theta == 0) {
#'     theta + 10^-6
#'   }
#'
#'    sind(lambda)^2 - ((tand(theta) * sind(2 * lambda)) / 2 * cosd(phi))
#' }
#'
#' rake <- function(fault) {
#'   ve <- fault[, 5]
#'   p <- Plane(fault[, 1], fault[, 2])
#'   a <- cbind(p[, 1] - 90, rep(0, nrow(p))) |>
#'     as.line() |>
#'     line2vec()
#'   b <- Line(fault[, 3], fault[, 4]) |> line2vec()
#'
#'   cos_theta <- vdot(a, b) / (sqrt(a^2) * sqrt(b^2))
#'   theta <- ve * acos(theta) |> rad2deg()
#'   # theta = vangle(a, b)
#'   theta
#' }
#'
#' Ry_analysis <- function(fault, ny = 1){
#'   theta = rake(fault)
#'   strike = fault[, 1]
#'   dip = fault[, 2]
#'   ve <- fault[, 5]
#'
#'   res <- cbind(fault = integer(), y = numeric(), R = numeric())
#'   for(f in 1:strike){
#'     for(y in seq(0, 360, ny)){
#'     R <- stress_ratio_simongomez(lambda = strike[f] + ve[f] * y, theta = theta[f], phi = dip[f])
#'     Rp = R
#'     if(!is.na(R)){
#'     if(R > 1){
#'       Rp = 1 + ((R-1)/R)
#'     } else if(R < 0){
#'       Rp = -R/(R-1)
#'     }
#'     }
#'
#'     res <- rbind(res, cbind(fault = f, y, R = Rp))
#'     }
#'   }
#'   res
#' }
