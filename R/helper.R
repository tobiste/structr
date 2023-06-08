rad2deg <- function(rad) {
  rad * 180 / pi
}

deg2rad <- function(deg) {
  deg * pi / 180
}

DEG2RAD <- function() {
  pi / 180
}


sind <- function(x) {
  sinpi(x / 180)
}

cosd <- function(x) {
  cospi(x / 180)
}

tand <- function(x) {
  sinpi(x / 180) / cospi(x / 180)
}


asind <- function(x) {
  asin(x) * 180 / pi
}

acosd <- function(x) {
  acos(x) * 180 / pi
}

atand <- function(x) {
  atan(x) * 180 / pi
}

atan2d <- function(x1, x2) {
  atan2(x1, x2) * 180 / pi
}

cot <- function(x) {
  1 / tan(x)
}

cotd <- function(x) {
  1 / tand(x)
}
