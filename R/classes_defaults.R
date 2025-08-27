vec2mat <- function(x) {
  if (is.null(dim(x))) {
    vm <- as.matrix(t(x))
  } else {
    vm <- as.matrix(x)
  }
  rownames(vm) <- rownames(x)
  vm
}


fol2vec <- function(azi, inc) {
  azi <- azi * pi / 180
  inc <- inc * pi / 180
  res <- cbind(
    x = -cos(azi) * sin(inc),
    y = -sin(azi) * sin(inc),
    z = cos(inc)
  )
  rownames(res) <- rownames(azi)
  res
}

lin2vec <- function(azi, inc) {
  azi <- azi * pi / 180
  inc <- inc * pi / 180
  res <- cbind(
    x = cos(azi) * cos(inc),
    y = sin(azi) * cos(inc),
    z = sin(inc)
  )
  rownames(res) <- rownames(azi)
  res
}

vec2lin <- function(x, y, z) {
  n <- vnorm(cbind(x, y, z)) # normalized vector
  nz <- n[, 3]
  azimuth_rad <- atan2(n[, 2], n[, 1]) %% (2*pi)
  plunge_rad <- asin(nz)
  azimuth <- azimuth_rad * 180 / pi
  plunge <- plunge_rad * 180 / pi

  res <- mapply(correct_inc, azi = azimuth, inc = plunge) |> t()
  rownames(res) <- rownames(x)
  colnames(res) <- c("azimuth", "plunge")
  res
}

vec2fol <- function(x, y, z) {
  n <- vnorm(cbind(x, y, z)) # normalized vector
  nz <- n[, 3]

  dip_direction_rad <- (atan2(n[, 2], n[, 1]) + pi) %% (2*pi)
  dip_rad <- pi/2 - asin(nz)

  dip_direction <- dip_direction_rad * 180 / pi
  dip <- dip_rad * 180 / pi

  res <- mapply(correct_inc, azi = dip_direction, inc = dip) |> t()
  rownames(res) <- rownames(x)
  colnames(res) <- c("dip_direction", "dip")
  res
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
