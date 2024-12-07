#' Fault displacement components
#'
#' Calculate fault displacement from (at least three) components of of the
#' fault displacement system
#'
#' @param dip fault's dip angle (in degrees)
#' @param delta angle between horizontal displacement vector and fault's strike
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
#' @return data.frame
#'
#' @examples
#' \dontrun{
#' fault_displacements(strikeslip = 2, verticalthrow = -5, heave = 3)
#' }
fault_displacements <- function(dip = NULL, delta = NULL, rake = NULL, verticalthrow = NULL, dipslip = NULL, heave = NULL, netslip = NULL, horizontalthrow = NULL, strikeslip = NULL) {
  horizontal <- NULL
   if (!is.null(dip) & !is.null(verticalthrow) & !is.null(delta)) {
    # Slip components in the vertical plane perpendicular to the strike of the fault

    dipslip <- verticalthrow /sind(dip)
    heave <- sqrt(dipslip^2 - verticalthrow^2)

    # Slip components in the horizontal plane
    horizontalthrow <- heave /sind(delta)
    strikeslip <- sqrt(horizontal^2 - heave^2)

    # Slip components in the fault plane plane
    netslip <- sqrt(strikeslip^2 + dipslip^2)
    rake <-asind(dipslip / netslip)
  } else if (!is.null(delta) & !is.null(horizontalthrow) & !is.null(dip)) {
    # Slip components in the horizontal plane
    strikeslip <- abs(cosd(delta) * horizontalthrow)
    heave <- sqrt(horizontalthrow^2 - strikeslip^2)

    # Slip components in the vertical plane perpendicular to the strike of the fault
    dipslip <- heave /cosd(dip)
    verticalthrow <- sqrt(dipslip^2 - heave^2)

    # Slip components in the fault plane plane
    netslip <- sqrt(strikeslip^2 + dipslip^2)
    rake <-atand(dipslip / strikeslip)
  } else if (!is.null(rake) & !is.null(verticalthrow) & !is.null(dip)) {
    dipslip <- verticalthrow /sind(dip)
    heave <- sqrt(dipslip^2 - verticalthrow^2)
    strikeslip <- dipslip / tand(rake)
    horizontalthrow <- sqrt(strikeslip^2 + heave^2)
    netslip <- sqrt(strikeslip^2 + dipslip^2)
  } else if (!is.null(strikeslip) & !is.null(verticalthrow) & !is.null(heave)) {
    dipslip <- sqrt(heave^2 + verticalthrow^2)
    dip <- atand(verticalthrow / heave)

    horizontalthrow <- sqrt(heave^2 + strikeslip^2)
    delta <-atand(heave / strikeslip)

    # netslip = sqrt(dipslip^2 + strikeslip^2)
    netslip <- sqrt(strikeslip^2 + verticalthrow^2 + heave^2)
    rake <-atand(dipslip / strikeslip)
  }


  data.frame(
    "dip" = dip,
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



#' Fault deformation gradient tensor
#'
#' Creates the fault deformation gradient tensor from displacement components.
#' If the dip direction is know, the tensor will be rotated into the
#' geographic reference frame.
#'
#' @param displacements data.frame containing the strike slip, the
#' heave, and the vertical throw of the fault's displacement components.
#' @param dip_direction (optional) dip direction in degrees.
#' If `NULL`, the fault tensor will be given in the fault displacement
#' coordinates. Otherwise, the tensor will be in the geographic reference frame.
#'
#' @details
#' x axis of tensor = heave, y = strike slip, z = vertical throw (positive for 
#' thrusting, negative for normal faulting)
#'
#' @examples
#' \dontrun{
#' fault_tensor(displacements = 
#'   data.frame(strikeslip = 2, verticalthrow = -5, heave = 3), 
#'   dip_direction = 0)
#' }
fault_tensor <- function(displacements, dip_direction = NULL) {
  A <- diag(c(displacements$heave, displacements$strikeslip, displacements$verticalthrow), 3, 3)

  if (!is.null(dip_direction)) {
    dipdir_rad <- deg2rad(dip_direction)
    # coordinates are measured from E!!!
    # c(1, 0, 0) = E
    # c(0, 1, 0) = N
    A[, 1] <- vrotate(A[, 1], A[, 3], -dipdir_rad) |> as.vector()
    A[, 2] <- vrotate(A[, 2], A[, 3], pi - dipdir_rad) |> as.vector()
  }

  A
}

# t <- fault_tensor(displacements = data.frame(strikeslip = 2, verticalthrow = -5, heave = 3))
# test <- fault_tensor(displacements = data.frame(strikeslip = 2, verticalthrow = -5, heave = 3)) |> svd()
# test$values
#
# test2 <- fault_tensor(displacements = data.frame(strikeslip = 2, verticalthrow = -5, heave = 3), dip_direction = 45)
# test2 |> t() |> vlength()
#
# vec2line(test2[, 2])
#
# test3 <- svd(test2)
# test3$u[, 2] |> vec2line()
# test3$u[, 3] |> vec2line()
#
# vrotate(x, z, (-dipdir_rad)) |> as.vector() |> vec2line()
# vrotate(y, z, (pi-dipdir_rad)) |> as.vector() |> vec2line()
#
#
# eigen(test2)
#
# diag(t)^2 |> sum() |> sqrt()
#
# vec2line(test2[, 1])
# vec2line(test2[, 2])
#
# # netslip vector
# diag(test2) |> vec2line()
#
# # project netslip vector onto horizontal plan = horizontal throw
# vtransform(
#   diag(test2),
#   cbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0))
#   ) |>
#   vec2line()
#
# 326%%180
#
#
#
# 146-45
# c(test2[3, 1], test2[2, 2], test2[1, 3]) |> vec2line()
# c(test2[3, 1], test2[2, 2], test2[1, 3]) |> vec2line()
# c(test2[1, 1], test2[1, 2], test2[1, 3]) |> vec2line()
# c(test2[2, 1], test2[2, 2], test2[2, 3]) |> vec2line()
# c(test2[3, 1], test2[3, 2], test2[3, 3]) |> vec2line()
#
# vtransform((diag(test2)), (cbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0))))



#' #' Title
#' #'
#' #' @param n,l,m direction cosines of the plane referred to the stress axes (σx, σy, σz) system and θ is the striation (resolved shear stress) rake on the fault plane
#' #' @param R stress ratio R=(σz−σx)/(σy–σx)
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
#' #' @param lambda angle between y (azimuth of σy) and the azimuth (α) of the
#' #' fault, therefore λ=α±y (it is added or subtracted according to the sense of movement: dextral or sinistral, respectively) (in degree)
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
