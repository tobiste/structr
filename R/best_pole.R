#' The cone or plane best fit of conically or cylindrical disposed s-plane poles
#'
#' Finding the best fit pole of rotation for a given set of points that are
#' assumed to lie on a mutual small or great circle circle

#' @param x matrix. Cartesian coordinates of points
#' @importFrom tectonicr deg2rad rad2deg
#' @importFrom dplyr mutate summarise
#' @name best_pole
#' @source Ramsay, 1967, p. 18-21
#' @returns numeric vector with 
#' \describe{
#' \item{`x`,`y`,`z`}{Cartesian coordinates of best fit pole of plane or cone axis,}
#' \item{`K`}{(only for cones) half apical angle of best fit cone, and}
#' \item{`e`}{misfit value of the least square of the deviations of the observed poles to the planes from the best fit pole.}
#' }
#' @examples
#' # example from Ramsay, 1967, p. 20
#' x <- rbind(
#'   c(-67, -31, -71),
#'   c(-62, -53, -50),
#'   c(-62, -75, -34),
#'   c(-58, 85, -34),
#'   c(-79, 40, -52),
#'   c(90, 14, -75),
#'   c(80, 10, 90)
#' ) |> acoscartesian_to_cartesian()
#' best_cone(x) # expect: c(acoscartesian_to_cartesian(cbind(31.1, -81, -60.5)), 89.5, NA)
#' best_plane(x) # expect: c(acoscartesian_to_cartesian(cbind(31.6, -81.1, -60)), 1-1.002)
NULL

#' @rdname best_pole
#' @export
best_cone <- function(x) {
  xsum <- data.frame(l = x[, 1], m = x[, 2], n = x[, 3]) |>
    dplyr::mutate(
      l2 = l^2,
      m2 = m^2,
      lm = l * m,
      ln = l * n,
      mn = m * n
    ) |>
    dplyr::summarise(
      l = sum(l),
      m = sum(m),
      n = sum(n),
      l2 = sum(l2),
      m2 = sum(m2),
      lm = sum(lm),
      ln = sum(ln),
      mn = sum(mn)
    )
  N <- nrow(xsum)

  # Rule of Sarrus:
  D <- Da <- Db <- Dc <- matrix(nrow = 3, ncol = 3)
  D[1, 1] <- xsum$l2
  D[1, 2] <- D[2, 1] <- xsum$lm
  D[1, 3] <- D[3, 1] <- xsum$l
  D[2, 2] <- xsum$m2
  D[2, 3] <- D[3, 2] <- xsum$m
  D[3, 3] <- N

  Da[1, 1] <- -xsum$ln
  Da[1, 2] <- D[1, 2]
  Da[1, 3] <- D[1, 3]
  Da[2, 1] <- -xsum$mn
  Da[2, 2] <- D[2, 2]
  Da[2, 3] <- Da[3, 2] <- D[2, 3]
  Da[3, 1] <- -xsum$n
  Da[3, 3] <- N

  Db[1, 1] <- D[1, 1]
  Db[1, 2] <- Da[1, 1]
  Db[1, 3] <- D[1, 3]
  Db[2, 1] <- D[1, 2]
  Db[2, 2] <- Da[2, 1]
  Db[2, 3] <- D[2, 3]
  Db[3, 1] <- D[1, 3]
  Db[3, 2] <- Da[3, 1]
  Db[3, 3] <- N

  Dc[1, 1] <- D[1, 1]
  Dc[1, 2] <- D[1, 2]
  Dc[1, 3] <- Da[1, 1]
  Dc[2, 1] <- D[1, 2]
  Dc[2, 2] <- D[2, 2]
  Dc[2, 3] <- Da[2, 1]
  Dc[3, 1] <- D[1, 3]
  Dc[3, 2] <- D[2, 3]
  Dc[3, 3] <- Da[3, 1]

  A <- det(Da) / det(D)
  B <- det(Db) / det(D)
  C <- det(Dc) / det(D)

  # gamma <- -acos((1 + A^2 + B^2)^(-1 / 2))
  # alpha <- pi - acos(A * (1 + A^2 + B^2)^(-1 / 2))
  # beta <- -acos(B * (1 + A^2 + B^2)^(-1 / 2))
  cos_gamma <- 1/sqrt(1 + A^2 + B^2)
  cos_alpha <- A * cos_gamma
  cos_beta <- B * cos_gamma
  cos_K <- -C * cos_gamma
  
  alpha <- acos(cos_alpha)
  beta <- acos(cos_beta)
  gamma <- acos(cos_gamma)
  
  # half apical angle
  K <- tectonicr:::acosd(cos_K) #|> tectonicr::deviation_norm()
  
  e <- cos(alpha)^2 + cos(beta)^2 + cos(gamma)^2

  # correct for lower hemisphere and convert to Cartesian coordinates
  cart <- cbind(pi-alpha, -beta, -gamma) |>
    tectonicr::rad2deg() |>
    acoscartesian_to_cartesian()
  #names(cart) <- NULL

  return(c(cart[, 1], cart[, 2], cart[, 3], "K" = K, "e" = 1-e))
}

#' @rdname best_pole
#' @export
best_plane <- function(x) {
  xsum <- data.frame(l = x[, 1], m = x[, 2], n = x[, 3]) |>
    dplyr::mutate(
      l2 = l^2,
      m2 = m^2,
      lm = l * m,
      ln = l * n,
      mn = m * n
    ) |>
    dplyr::summarise(
      l = sum(l),
      m = sum(m),
      n = sum(n),
      l2 = sum(l2),
      m2 = sum(m2),
      lm = sum(lm),
      ln = sum(ln),
      mn = sum(mn)
    )

  t <- 1/(xsum$l2 * xsum$m2 - (xsum$lm)^2)
  A <- (xsum$lm * xsum$mn - xsum$ln * xsum$m2) * t
  B <- (xsum$lm * xsum$ln - xsum$mn * xsum$l2) * t

  cos_gamma <- 1 / sqrt(1 + A^2 + B^2)
  cos_alpha <- A * cos_gamma
  cos_beta <- B * cos_gamma

  alpha = acos(cos_alpha) 
  beta = acos(cos_beta)
  gamma = acos(cos_gamma)
  
  cart <- cbind(-alpha, -beta, -gamma) |>
    tectonicr::rad2deg() |>
    acoscartesian_to_cartesian()
  #names(cart) <- NULL
  
  e <- cos(alpha)^2 + cos(beta)^2 + cos(gamma)^2

  return(
    c(cart[, 1], cart[, 2], cart[, 3], e = 1-e)
  )
}

#' Spherical coordinate conversion
#'
#' @inheritParams best_cone
#' @details acoscartesian coordinates are given as
#' \describe{
#' \item{\code{a}}{deviation from x-axis E-W horizontal (angle in degrees)}
#' \item{\code{b}}{deviation from y-axis N-S horizontal (angle in degrees)}
#' \item{\code{c}}{deviation from z-axis vertical (angle in degrees)}
#' }
#' @source Ramsay, 1967, p. 15-16
#' @name ramsay_coords
#' @examples
#' Stereographic coordinates (angle notation):
#' x <- rbind(
#'   c(58, 78, 35),
#'   c(47, 72, 47),
#'   c(57, 68, 41),
#'   c(56, 64, 45),
#'   c(52, 65, 47),
#'   c(61, 62, 51),
#'   c(55, 60, 49),
#'   c(61, 62, 42),
#'   c(68, 52, 46),
#'   c(42, 58, 65)
#' )
#' acoscartesian_to_cartesian(x)
#' acoscartesian_to_cartesian(x) |> to_spherical()
#' acoscartesian_to_cartesian(x) |> cartesian_to_acoscartesian()
NULL

#' @rdname ramsay_coords
cartesian_to_acoscartesian <- function(x) {
  #acos(x) |> tectonicr::rad2deg()
  a <- acos(x[, 1])
  b <- acos(x[, 2])
  c <- acos(x[, 3])
   
  a <- ifelse(x[, 1] < 0, -a, a)
  b <- ifelse(x[, 2] < 0, -b, b)
  c <- ifelse(x[, 3] < 0, -c, c)
  tectonicr::rad2deg(cbind(a, b, c))
}

#' @rdname ramsay_coords
acoscartesian_to_cartesian <- function(x) {
  #tectonicr::deg2rad(x) |> cos()
  x <- tectonicr::deg2rad(x)
  cx <- cos(x[, 1])
  cy <- cos(x[, 2])
  cz <- cos(x[, 3])
  
  cx <- ifelse(x[, 1] < 0, -cx, cx)
  cy <- ifelse(x[, 2] < 0, -cy, cy)
  cz <- ifelse(x[, 3] < 0, -cz, cz)
  cbind(x= cx, y=cy, z=cz)
}
