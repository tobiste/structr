#' Least-square fit of small and great circles to spherically projected data
#'
#' Finds the best small and great circles using the algorithm by Gray et al. (1980)
#'
#' @param x numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
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
#' gray_cleavage <- filter(gray_example1, Type == "Cleavage")
#' gray_bedding <- filter(gray_example1, Type == "Bedding")
#' test_clea <- cbind(gray_cleavage$dipdir, gray_cleavage$dip) |> as.plane()
#' test_bedd <- cbind(gray_bedding$dipdir, gray_bedding$dip) |> as.plane()
#' best_clea <- best_fit_plane((test_clea))
#' best_bedd <- best_fit_plane((test_bedd))
#'
#' stereoplot()
#' stereo_point(test_clea, col = "blue")
#' stereo_point(test_bedd, col = "red")
#' stereo_smallcircle(best_clea$axis_c, best_clea$cone_angle, col = "lightblue")
#' stereo_smallcircle(best_clea$axis_g, 90, lty = 2, col = "blue")
#' stereo_smallcircle(best_bedd$axis_c, best_bedd$cone_angle, col = "sienna")
#' stereo_smallcircle(best_bedd$axis_g, 90, lty = 2, col = "red")
best_fit_plane <- function(x) {
  transform <- FALSE
  if (is.spherical(x)) {
    x <- as.line(x) |> to_vec()
    transform <- TRUE
  }

  p <- nrow(x)
  s_res <- gray_algorithm(x, sm = TRUE)
  g_res <- gray_algorithm(x, sm = FALSE)
  K_s <- acos(s_res$cos_K) |> as.numeric()
  K_g <- pi / 2 # == acos(g_res$cos_K) |> as.numeric()


  # R <- K - acos(A_eigen3[1] * t(a1) + A_eigen3[2] * t(a2) + A_eigen3[3] * t(a3))
  R_s <- K_s - acos(s_res$eig3 %*% t(x)) # residual
  r_s <- sum(R_s^2) # sum of squares of residuals

  R_g <- K_g - acos(g_res$eig3 %*% t(x)) # residual
  r_g <- sum(R_g^2) # sum of squares of residuals

  E_s <- (p - 2) * r_s^2
  E_g <- (p - 3) * r_g^2

  Vr <- (p - 3) * ((r_g - r_s) / r_s)

  axis_s <- vec2mat(s_res$eig3)
  axis_g <- vec2mat(g_res$eig3)
  if (transform) {
    axis_s <- to_spherical(axis_s)
    axis_g <- to_spherical(axis_g)
    K_s <- K_s * 180 / pi
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
