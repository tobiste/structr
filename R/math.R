vlength <- function(x) {
  sqrt(x[1]^2 + x[2]^2 + x[3]^2)
}

vnorm <- function(x) {
  x / vlength(x)
}

vcross <- function(x, y) {
  vx <- x[2] * y[3] - x[3] * y[2]
  vy <- x[3] * y[1] - x[1] * y[3]
  vz <- x[1] * y[2] - x[2] * y[1]
  cbind(x = vx, y = vy, z = vz)
}


vrotate <- function(x, rotaxis, rotangle) {
  vax <- vcross(rotaxis, x)
  x + vax * sin(rotangle) + vcross(rotaxis, vax) * 2 * (sin(rotangle / 2))^2
}

fol2vec <- function(azi, inc) {
  azi <- azi * DEG2RAD()
  inc <- inc * DEG2RAD()
  c(
    -cos(azi) * sin(inc),
    -sin(azi) * sin(inc),
    cos(inc)
  )
}

lin2vec <- function(azi, inc) {
  azi <- azi * DEG2RAD()
  inc <- inc * DEG2RAD()
  c(
    cos(azi) * cos(inc),
    sin(azi) * cos(inc),
    sin(inc)
  )
}

vec2lin <- function(v) {
  n <- vnorm(v) # normalized vector
  if (n[3] < 0) {
    n <- -n
  }
  c(
    (atan2(n[2], n[1]) / DEG2RAD()) %% 360,
    asin(n[3]) / DEG2RAD()
  )
}

vec2fol <- function(v) {
  n <- vnorm(v) # normalized vector
  if (n[3] < 0) {
    n <- -n
  }
  c(
    ((atan2(n[2], n[1]) / DEG2RAD()) + 180) %% 360,
    90 - (asin(n[3]) / DEG2RAD())
  )
}
