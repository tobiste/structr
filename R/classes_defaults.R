vec2mat <- function(x) {
  vm <- if (is.null(dim(x))) {
    as.matrix(t(x))
  } else {
    as.matrix(x)
  }
  rownames(vm) <- rownames(x)
  vm
}


fol2vec <- function(azi, inc) {
  radians <- pi / 180

  azi <- azi * radians
  inc <- inc * radians

  # res <- cbind(
  #   x = -cos(azi) * sin(inc),
  #   y = -sin(azi) * sin(inc),
  #   z = cos(inc)
  # )

  ca <- cos(azi)
  sa <- sin(azi)
  ci <- cos(inc)
  si <- sin(inc)

  matrix(
    c(
      -ca * si,
      -sa * si,
      ci
    ),
    ncol = 3,
    dimnames = list(
      rownames(azi),
      c("x", "y", "z")
    )
  )
  # rownames(res) <- rownames(azi)
  # res
}

lin2vec <- function(azi, inc) {
  radians <- pi / 180
  azi <- azi * radians
  inc <- inc * radians

  ca <- cos(azi)
  sa <- sin(azi)
  ci <- cos(inc)
  si <- sin(inc)

  matrix(
    c(
      ca * ci,
      sa * ci,
      si
    ),
    ncol = 3,
    dimnames = list(
      rownames(azi),
      c("x", "y", "z")
    )
  )
}

vec2lin <- function(x, y, z) {
  xyz <- matrix(
    c(x, y, z),
    ncol = 3
  )
  n <- vnorm(xyz) # normalized vector
  nz <- n[, 3]
  azimuth_rad <- atan2(n[, 2], n[, 1]) %% (2 * pi)
  plunge_rad <- asin(nz)

  degree <- 180 / pi
  azimuth <- azimuth_rad * degree
  plunge <- plunge_rad * degree

  # res <- mapply(correct_inc, azi = azimuth, inc = plunge) |> t()

  i <- plunge < 0
  plunge[i] <- -plunge[i]
  azimuth[i] <- (azimuth[i] + 180) %% 360

  # res <- cbind(
  #   azimuth = azimuth,
  #   plunge = plunge
  # )
  # rownames(res) <- rownames(x)
  # colnames(res) <- c("azimuth", "plunge")
  # res

  matrix(
    c(
      azimuth,
      plunge
    ),
    ncol = 2,
    dimnames = list(
      rownames(x),
      c("azimuth", "plunge")
    )
  )
}

vec2fol <- function(x, y, z) {
  n <- vnorm(cbind(x, y, z)) # normalized vector
  nz <- n[, 3]

  dip_direction_rad <- (atan2(n[, 2], n[, 1]) + pi) %% (2 * pi)
  dip_rad <- pi / 2 - asin(nz)

  degree <- 180 / pi
  dip_direction <- dip_direction_rad * degree
  dip <- dip_rad * degree

  i <- dip < 0
  dip[i] <- -dip[i]
  dip_direction[i] <- (dip_direction[i] + 180) %% 360

  # res <- cbind(
  #   dip_direction = dip_direction,
  #   dip = dip
  # )
  # res <- mapply(correct_inc, azi = dip_direction, inc = dip) |> t()
  # rownames(res) <- rownames(x)
  # colnames(res) <- c("dip_direction", "dip_direction")
  # res

  matrix(
    c(
      dip_direction,
      dip
    ),
    ncol = 2,
    dimnames = list(
      rownames(x),
      c("dip_direction", "dip")
    )
  )
}

correct_inc <- function(azi, inc) {
  if (inc > 90) {
    inc <- 180 - inc
    azi <- (azi + 180) %% 360
  } else if (inc < 0) {
    inc <- abs(inc)
    azi <- (azi + 180) %% 360
  }
  c(azi, inc)
}
