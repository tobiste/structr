#' Least-square fit of small and great circles to spherically projected data
#'
#' Finds the best small and great circles using the algorithm by Gray et al. (1980)
#'
#' @param x object of class `"Vec3"`, `"Line"`, or `"Plane"`.
#'
#' @returns list
#' @references Gray, N.H., Geiser, P.A., Geiser, J.R. (1980). On the
#' least-square fit of small and great circles to spherically projected data.
#' Mathematical Geology, Vol. 12, No. 3, 1980.
#' @export
#'
#' @examples
#' data("gray_example")
#' gray_example1 <- gray_example |>
#'   dplyr::mutate(
#'     dipdir = Strike + 90,
#'     dip = Dip,
#'     id = seq_along(dip)
#'   )
#' gray_cleavage <- dplyr::filter(gray_example1, Type == "Cleavage")
#' gray_bedding <- dplyr::filter(gray_example1, Type == "Bedding")
#' test_clea <- Plane(gray_cleavage$dipdir, gray_cleavage$dip)
#' test_bedd <- Plane(gray_bedding$dipdir, gray_bedding$dip)
#' best_clea <- best_fit_plane(test_clea)
#' best_bedd <- best_fit_plane(test_bedd)
#'
#' stereoplot()
#' points(test_clea, col = "blue")
#' points(test_bedd, col = "red")
#' lines(best_clea$axis_c, best_clea$cone_angle, col = "lightblue")
#' lines(best_clea$axis_g, 90, lty = 2, col = "blue")
#' lines(best_bedd$axis_c, best_bedd$cone_angle, col = "sienna")
#' lines(best_bedd$axis_g, 90, lty = 2, col = "red")
best_fit_plane <- function(x) {
  xv <- Vec3(x) |> unclass()

  p <- nrow(xv)
  s_res <- gray_algorithm(xv, sm = TRUE)
  g_res <- gray_algorithm(xv, sm = FALSE)
  K_s <- acos(s_res$cos_K) |> as.numeric()
  K_g <- pi / 2 # == acos(g_res$cos_K) |> as.numeric()


  # R <- K - acos(A_eigen3[1] * t(a1) + A_eigen3[2] * t(a2) + A_eigen3[3] * t(a3))
  R_s <- K_s - acos(s_res$eig3 %*% t(xv)) # residual
  r_s <- sum(R_s^2) # sum of squares of residuals

  R_g <- K_g - acos(g_res$eig3 %*% t(xv)) # residual
  r_g <- sum(R_g^2) # sum of squares of residuals

  E_s <- (p - 2) * r_s^2
  E_g <- (p - 3) * r_g^2

  Vr <- (p - 3) * ((r_g - r_s) / r_s)

  axis_s <- as.Vec3(s_res$eig3)
  axis_g <- as.Vec3(g_res$eig3)

  if (is.Line(x) | is.Plane(x)) {
    # if (is.Line(x)) {
    axis_s <- Line(axis_s)
    axis_g <- Line(axis_g)
    # } else if (is.Plane(x)) {
    #  axis_s <- Plane(axis_s)
    #  axis_g <- Plane(axis_g)
    # }
    K_s <- rad2deg(K_s)
  }

  return(
    list(
      axis_c = axis_s,
      axis_g = axis_g,
      cone_angle = K_s,
      r_s = r_s,
      r_g = r_g,
      E_s = E_s,
      E_g = E_g,
      Vr = Vr
    )
  )
}


gray_algorithm <- function(x, sm = TRUE) {
  a1 <- x[, 1]
  a2 <- x[, 2]
  a3 <- x[, 3]

  p <- length(a1)
  if (sm) {
    a1_bar <- sum(a1 / p)
    a2_bar <- sum(a2 / p)
    a3_bar <- sum(a3 / p)
  } else {
    a1_bar <- a2_bar <- a3_bar <- 0
  }

  a1_d <- a1 - a1_bar
  a2_d <- a2 - a2_bar
  a3_d <- a3 - a3_bar

  A1 <- rbind(
    sum(a1 * a1_d),
    sum(a2 * a1_d),
    sum(a3 * a1_d)
  )
  A2 <- rbind(
    sum(a1 * a2_d),
    sum(a2 * a2_d),
    sum(a3 * a2_d)
  )
  A3 <- rbind(
    sum(a1 * a3_d),
    sum(a2 * a3_d),
    sum(a3 * a3_d)
  )
  A <- cbind(A1, A2, A3)

  eig3 <- -eigen(A)$vectors[, 3]

  # cos_K <- eig3[1] * a1_bar + eig3[2] * a2_bar + eig3[3] * a3_bar
  cos_K <- eig3 %*% t(cbind(a1_bar, a2_bar, a3_bar))

  return(list(eig3 = eig3, cos_K = cos_K))
}





#' The cone or plane best fit of conically or cylindrical disposed s-plane poles
#'
#' Finding the best fit pole of rotation for a given set of points that are
#' assumed to lie on a mutual small or great circle circle
#'
#' @param x matrix. Cartesian coordinates of points
#' @importFrom dplyr mutate summarise
#' @references Ramsay, 1967, p. 18-21
#' @returns numeric vector with
#' \describe{
#' \item{`x`,`y`,`z`}{Cartesian coordinates of best fit pole of plane or cone axis,}
#' \item{`e`}{residual of the sum of square of the deviations of the observed poles to the planes from the best fit pole, and}
#' \item{`K`}{(only for cones) half apical angle of best fit cone (in radians).}
#' }
#' @name best_pole
#' @examples
#' \dontrun{
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
#' best_cone_ramsay(x) # expect: c(0.856, -0.157, -0.492, NA, 1.56207)
#' best_plane_ramsay(x) # expect: c(0.852, -0.154, -0.502, 1-1.002)
#' }
NULL

#' @rdname best_pole
# best_cone_ramsay <- function(x) {
#   l <- m <- n <- l2 <- m2 <- lm <- ln <- mn <- numeric()
#   xsum <- data.frame(l = x[, 1], m = x[, 2], n = x[, 3]) |>
#     dplyr::mutate(
#       l2 = l^2,
#       m2 = m^2,
#       lm = l * m,
#       ln = l * n,
#       mn = m * n
#     ) |>
#     dplyr::summarise(
#       l = sum(l),
#       m = sum(m),
#       n = sum(n),
#       l2 = sum(l2),
#       m2 = sum(m2),
#       lm = sum(lm),
#       ln = sum(ln),
#       mn = sum(mn)
#     )
#   N <- nrow(xsum)
# 
#   D <- Da <- Db <- Dc <- matrix(nrow = 3, ncol = 3)
#   D[1, 1] <- xsum$l2
#   D[1, 2] <- D[2, 1] <- xsum$lm
#   D[1, 3] <- D[3, 1] <- xsum$l
#   D[2, 2] <- xsum$m2
#   D[2, 3] <- D[3, 2] <- xsum$m
#   D[3, 3] <- N
# 
#   Da[1, 1] <- -xsum$ln
#   Da[1, 2] <- D[1, 2]
#   Da[1, 3] <- D[1, 3]
#   Da[2, 1] <- -xsum$mn
#   Da[2, 2] <- D[2, 2]
#   Da[2, 3] <- Da[3, 2] <- D[2, 3]
#   Da[3, 1] <- -xsum$n
#   Da[3, 3] <- N
# 
#   Db[1, 1] <- D[1, 1]
#   Db[1, 2] <- Da[1, 1]
#   Db[1, 3] <- D[1, 3]
#   Db[2, 1] <- D[1, 2]
#   Db[2, 2] <- Da[2, 1]
#   Db[2, 3] <- D[2, 3]
#   Db[3, 1] <- D[1, 3]
#   Db[3, 2] <- Da[3, 1]
#   Db[3, 3] <- N
# 
#   Dc[1, 1] <- D[1, 1]
#   Dc[1, 2] <- D[1, 2]
#   Dc[1, 3] <- Da[1, 1]
#   Dc[2, 1] <- D[1, 2]
#   Dc[2, 2] <- D[2, 2]
#   Dc[2, 3] <- Da[2, 1]
#   Dc[3, 1] <- D[1, 3]
#   Dc[3, 2] <- D[2, 3]
#   Dc[3, 3] <- Da[3, 1]
# 
#   A <- det(Da) / det(D)
#   B <- det(Db) / det(D)
#   C <- det(Dc) / det(D)
# 
#   # gamma <- -acos((1 + A^2 + B^2)^(-1 / 2))
#   # alpha <- pi - acos(A * (1 + A^2 + B^2)^(-1 / 2))
#   # beta <- -acos(B * (1 + A^2 + B^2)^(-1 / 2))
#   cos_gamma <- 1 / sqrt(1 + A^2 + B^2)
#   cos_alpha <- A * cos_gamma
#   cos_beta <- B * cos_gamma
#   cos_K <- -C * cos_gamma
# 
#   # half apical angle
#   K <- acos(cos_K)
# 
#   cart <- cbind(x = -cos_alpha, y = -cos_beta, z = -cos_gamma)
#   e <- cos_alpha^2 + cos_beta^2 + cos_gamma^2
# 
#   # alpha <- acos(cos_alpha)
#   # beta <- acos(cos_beta)
#   # gamma <- acos(cos_gamma)
# 
#   # correct for lower hemisphere and convert to Cartesian coordinates
#   # cart <- cbind(pi-alpha, -beta, -gamma) |>
#   #   tectonicr::rad2deg() |>
#   #   acoscartesian_to_cartesian()
#   # #names(cart) <- NULL
#   # e <- cos(alpha)^2 + cos(beta)^2 + cos(gamma)^2
# 
# 
#   return(c(cart[, 1], cart[, 2], cart[, 3], "e" = 1 - e, "K" = K))
# }
best_cone_ramsay <- function(x) {
  # ensure x is a matrix
  x <- as.matrix(x)
  l <- x[, 1]
  m <- x[, 2]
  n <- x[, 3]
  
  # precompute sums directly (avoids mutate/summarise)
  l2 <- sum(l^2)
  m2 <- sum(m^2)
  lm <- sum(l * m)
  ln <- sum(l * n)
  mn <- sum(m * n)
  l_sum <- sum(l)
  m_sum <- sum(m)
  n_sum <- sum(n)
  N <- nrow(x)
  
  # construct matrices directly
  D  <- matrix(c(l2, lm, l_sum,
                 lm, m2, m_sum,
                 l_sum, m_sum, N), nrow = 3, byrow = TRUE)
  
  Da <- matrix(c(-ln,  lm,  l_sum,
                 -mn,  m2,  m_sum,
                 -n_sum, m_sum, N), nrow = 3, byrow = TRUE)
  
  Db <- matrix(c(l2,   -ln, l_sum,
                 lm,   -mn, m_sum,
                 l_sum, -n_sum, N), nrow = 3, byrow = TRUE)
  
  Dc <- matrix(c(l2,  lm,  -ln,
                 lm,  m2,  -mn,
                 l_sum, m_sum, -n_sum), nrow = 3, byrow = TRUE)
  
  # determinants
  detD  <- det(D)
  A <- det(Da) / detD
  B <- det(Db) / detD
  C <- det(Dc) / detD
  
  # direction cosines
  cos_gamma <- 1 / sqrt(1 + A^2 + B^2)
  cos_alpha <- A * cos_gamma
  cos_beta  <- B * cos_gamma
  cos_K     <- -C * cos_gamma
  
  # half apical angle
  K <- acos(cos_K)
  
  cart <- c(x = -cos_alpha, y = -cos_beta, z = -cos_gamma)
  e <- cos_alpha^2 + cos_beta^2 + cos_gamma^2
  
  c(cart, e = 1 - e, K = K)
}


#' @rdname best_pole
best_cone_ramsay2 <- function(x) {
  # ensure matrix
  x <- as.matrix(x)
  l <- x[, 1]
  m <- x[, 2]
  n <- x[, 3]
  
  # compute sums directly
  l2 <- sum(l^2)
  m2 <- sum(m^2)
  lm <- sum(l * m)
  ln <- sum(l * n)
  mn <- sum(m * n)
  
  # precompute
  t <- 1 / (l2 * m2 - lm^2)
  A <- (lm * mn - ln * m2) * t
  B <- (lm * ln - mn * l2) * t
  
  # direction cosines
  cos_gamma <- 1 / sqrt(1 + A^2 + B^2)
  cos_alpha <- A * cos_gamma
  cos_beta  <- B * cos_gamma
  
  # Cartesian coordinates
  cart <- c(x = -cos_alpha, y = -cos_beta, z = -cos_gamma)
  e <- cos_alpha^2 + cos_beta^2 + cos_gamma^2
  
  c(cart, e = 1 - e)
}

