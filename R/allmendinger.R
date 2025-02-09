#' Computes the coordinates of a line in an equal angle or equal area stereonet
#' of unit radius
#'
#' @param trd  trend of line
#' @param plg plunge of line
#' @param sttype character indicating the type of stereonet. `"eangle"` for equal angle
#' and `"earea"` for equal area
#'
#' @returns xp and yp are the coordinates of the line in the stereonet plot
#'
#' @note trend and plunge should be entered in radians
#' @export
#'
StCoordLine <- function(trd, plg, sttype = c("earea", "eangle")) {
  sttype <- match.arg(sttype)

  # Take care of negative plunges
  if (plg < 0) {
    trd <- ZeroTwoPi(trd + pi)
    plg <- -plg
  }

  # Some constants
  piS4 <- pi / 4
  s2 <- sqrt(2)
  plgS2 <- plg / 2

  if (sttype == "eangle") {
    # Equal angle stereonet: From Equation 1.5 above
    # Also see Pollard and Fletcher (2005), eq.2.72
    xp <- tan(piS4 - plgS2) * sin(trd)
    yp <- tan(piS4 - plgS2) * cos(trd)
  } else if (sttype == "earea") {
    # Equal area stereonet: From Equation 1.6 above
    # Also see Pollard and Fletcher (2005), eq.2.90
    xp <- s2 * sin(piS4 - plgS2) * sin(trd)
    yp <- s2 * sin(piS4 - plgS2) * cos(trd)
  }
  cbind(xp, yp)
}


#' constrains azimuth to lie between 0 and 2*pi radians
#'
#' @param a azimuth in radians
#'
#' @returns azimuth from 0 to 2*pi
#' @export
#'
ZeroTwoPi <- function(a) {
  b <- a
  twopi <- 2 * pi
  if (b < 0) {
    b <- b + twopi
  } else if (b >= twopi) {
    b <- b - twopi
  }
  b
}


#' converts from spherical to cartesian coordinates
#'
#' @param trd,plg Angles in radians
#' @param k is an integer to tell whether the trend and plunge of a line
#' (k = 0) or strike and dip of a plane in right hand rule (k = 1) are being
#' sent in the trd and plg slots. In this last case, the direction cosines of
#' the pole to the plane are returned
#'
#' @returns the north (cn), east (ce), and down (cd) direction cosines of a line
#' @export
#'
SphToCart <- function(trd, plg, k) {
  if (k == 0) {
    # If line (see Table 2.1)
    cd <- sin(plg)
    ce <- cos(plg) * sin(trd)
    cn <- cos(plg) * cos(trd)
  } else if (k == 1) {
    # Else pole to plane (see Table 2.1)
    cd <- cos(plg)
    ce <- -sin(plg) * cos(trd)
    cn <- sin(plg) * sin(trd)
  }
  cbind(cn, ce, cd)
}


#' Converts from cartesian to spherical coordinates
#'
#' @param cn,ce,cd north (cn), east (ce), and down (cd) direction cosines
#'
#' @returns trend and plunge in radians
#' @export
#'
CartToSph <- function(cn, ce, cd) {
  # Plunge (see Table 2.1)
  plg <- asin(cd)

  # Trend
  # If north direction cosine is zero, trend is east or west
  # Choose which one by the sign of the east direction cosine
  if (cn == 0) {
    if (ce < 0) {
      trd <- 3 / 2 * pi # trend is west
    } else {
      trd <- pi / 2 # trend is east
    }
  } else {
    # Else use Table 2.1
    trd <- atan(ce / cn)
    if (cn < 0) {
      # Add pi
      trd <- trd + pi
    }
    # Make sure trd is between 0 and 2*pi
    trd <- ZeroTwoPi(trd)
  }
  cbind(plg, trd)
}


#' #' calculates the mean vector for a given series of lines
#' #'
#' #' @param t,p numeric vectors of trends and plunges (in radians)
#' #'
#' #' @returns calculates the trend (trd) and plunge (plg) of the mean vector (in radians),
#' #' its normalized length, and Fisher statistics (concentration factor (conc),
#' #' 99 (d99) and 95 (d95) % uncertainty cones, in radians)
#' #' @export
#' #'
#' #' @importFrom shiny p
#' CalcMV <- function(t, p) {
#'   # Number of lines
#'   nlines <- length(t)
#' 
#'   # Initialize the 3 direction cosines which contain the sums of the
#'   # individual vectors (i.e. the coordinates of the resultant vector)
#'   CNsum <- 0
#'   CEsum <- 0
#'   CDsum <- 0
#' 
#'   # Now add up all the individual vectors
#'   for (i in 1:nlines) {
#'     cs <- SphToCart(t(i), p(i), 0)
#'     cn <- cs[, 1]
#'     ce <- cs[, 2]
#'     cd <- cs[, 3]
#'     CNsum <- CNsum + cn
#'     CEsum <- CEsum + ce
#'     CDsum <- CDsum + cd
#'   }
#' 
#'   # R is the length of the resultant vector and Rave is the length of
#'   # the resultant vector normalized by the number of lines
#'   R <- sqrt(CNsum * CNsum + CEsum * CEsum + CDsum * CDsum)
#'   Rave <- R / nlines
#' 
#'   # If Rave is lower than 0.1, the mean vector is insignificant, return error
#'   if (Rave < 0.1) {
#'     stop("Mean vector is insignificant")
#'   } else {
#'     # Divide the resultant vector by its length to get the average
#'     # unit vector
#'     CNsum <- CNsum / R
#'     CEsum <- CEsum / R
#'     CDsum <- CDsum / R
#' 
#'     # Use the following 'if' statement if you want to convert the
#'     # mean vector to the lower hemisphere
#'     if (CDsum < 0) {
#'       CNsum <- -CNsum
#'       CEsum <- -CEsum
#'       CDsum <- -CDsum
#'     }
#' 
#'     # Convert the mean vector from direction cosines to trend and plunge
#'     vec <- CartToSph(CNsum, CEsum, CDsum)
#'     trd <- vec[, 1]
#'     plg <- vec[, 2]
#' 
#'     # If there are enough measurements calculate the Fisher Statistics
#'     # For more information on these statistics see Fisher et al. (1987)
#'     if (R < nlines) {
#'       if (nlines < 16) {
#'         afact <- 1 - (1 / nlines)
#'         conc <- (nlines / (nlines - R)) * afact^2
#'       } else {
#'         conc <- (nlines - 1) / (nlines - R)
#'       }
#'     }
#'     if (Rave >= 0.65 && Rave < 1) {
#'       afact <- 1 / 0.01
#'       bfact <- 1 / (nlines - 1)
#'       d99 <- acos(1 - ((nlines - R) / R) * (afact^bfact - 1.0))
#'       afact <- 1 / 0.05
#'       d95 <- acos(1 - ((nlines - R) / R) * (afact^bfact - 1.0))
#'     }
#'   }
#'   cbind(trd, plg, Rave, conc, d99, d95)
#' }



#' Angles between vectors
#'
#' Angles calculates the angles between two lines, between two planes,
#' the line which is the intersection of two planes, or the plane
#' containing two apparent dips
#'
#' @param trd1
#' @param plg1
#' @param trd2
#' @param plg2
#' @param ans0 a character that tells the function what to calculate:
#'  `ans0 = 'a'` -> the orientation of a plane given two apparent dips
#'  `ans0 = 'l'` -> the angle between two lines
#'  In the above two cases, the user sends the trend and plunge of two lines
#'  `ans0 = 'i'` -> the intersection of two planes
#'  `ans0 = 'p'` -> the angle between two planes
#'
#' @export
#'
Angles <- function(trd1, plg1, trd2, plg2, ans0 = c("a", "l", "i", "p")) {
  ans0 <- match.arg(ans0)

  # If planes have been entered
  if (ans0 == "i" || ans0 == "p") {
    k <- 1
  } else if (ans0 == "a" || ans0 == "l") {
    # Else if lines have been entered
    k <- 0
  }

  # Calculate the direction cosines of the lines or poles to planes
  c1 <- SphToCart(trd1, plg1, k)
  cn1 <- c1[, 1]
  ce1 <- c1[, 2]
  cd1 <- c1[, 3]

  c2 <- SphToCart(trd2, plg2, k)
  cn2 <- c2[, 1]
  ce2 <- c2[, 2]
  cd2 <- c2[, 3]

  # If angle between 2 lines or between the poles to 2 planes
  if (ans0 == "l" || ans0 == "p") {
    # Use dot product = Sum of the products of the direction cosines
    ans1 <- acos(cn1 * cn2 + ce1 * ce2 + cd1 * cd2)
    ans2 <- pi - ans1
  }

  # If intersection of two planes or pole to a plane containing two apparent dips
  if (ans0 == "a" || ans0 == "i") {
    # If the 2 planes or apparent dips are parallel, return an error
    if (trd1 == trd2 && plg1 == plg2) {
      stop("lines or planes are parallel")
    } else {
      # Else use cross product
      cn <- ce1 * cd2 - cd1 * ce2
      ce <- cd1 * cn2 - cn1 * cd2
      cd <- cn1 * ce2 - ce1 * cn2

      # Make sure the vector points down into the lower hemisphere
      if (cd < 0) {
        cn <- -cn
        ce <- -ce
        cd <- -cd
      }

      # Convert vector to unit vector by dividing it by its length
      r <- sqrt(cn * cn + ce * ce + cd * cd)

      # Calculate line of intersection or pole to plane
      vec <- CartToSph(cn / r, ce / r, cd / r)
      trd <- vec[, 1]
      plg <- vec[, 2]

      # If intersection of two planes
      if (ans0 == "i") {
        ans1 <- trd
        ans2 <- plg
      } else if (ans0 == "a") {
        # Else if plane containing two dips, calculate plane from its pole
        p <- Pole(trd, plg, 0)
        ans1 <- p[, 1]
        ans2 <- p[, 2]
      }
    }
  }
  cbind(ans1, ans2)
}


#' returns the pole to a plane or the plane which correspond to a pole
#'
#' @param trd,plg trend and plunge of a ploe, or strike and dip of a plane
#' @param k k is an integer that tells the program what to calculate.
#'
#' @returns If `k = 0`, returns the strike (trd1) and dip (plg1) of a plane, given the trend (trd)
#' and plunge (plg) of its pole.
#' If `k = 1`, returns the trend (trd1) and plunge (plg1) of a pole, given the
#' strike (trd) and dip (plg) of its plane.
#' @export
Pole <- function(trd, plg, k) {
  # Some constants
  east <- pi / 2

  # Calculate plane given its pole
  if (k == 0) {
    if (plg >= 0) {
      plg1 <- east - plg
      dipaz <- trd - pi
    } else {
      plg1 <- east + plg
      dipaz <- trd
    }

    # Calculate trd1 and make sure it is between 0 and 2*pi
    trd1 <- ZeroTwoPi(dipaz - east)
  } else if (k == 1) {
    # Else calculate pole given its plane
    c1 <- SphToCart(trd, plg, k)
    cn <- c1[, 1]
    ce <- c1[, 2]
    cd <- c1[, 3]
    vec <- CartToSph(cn, ce, cd)
    trd1 <- vec[, 1]
    plg1 <- vec[, 2]
  }
  cbind(trd1, plg1)
}



#' Constructs the down plunge projection of a bed
#'
#' @param bedseg The array `bedseg` is a two-dimensional array of size `npoints` x 3
#' which holds `npoints` on the digitized bed, each point defined by 3
#' coordinates: `X1` = East, `X2` = North, `X3` = Up
#' @param trd,plg trend (`trd`) and plunge (`plg`) of the fold axis (in radians)
#'
#' @returns ?
DownPlunge <- function(bedseg, trd, plg) {
  # Number of points in bed
  nvtex <- length(bedseg[, 1])

  # Allocate some arrays
  a <- matrix(0, nrow = 3, ncol = 3)
  dpbedseg <- matrix(0, nrow = nvtex, ncol = 3)

  # Calculate the transformation matrix a(i,j). The convention is that
  # the first index refers to the new axis and the second to the old axis.
  # The new coordinate system is with X3<U+2019> parallel to the fold axis, X1'
  # perpendicular to the fold axis and in the same vertical plane, and
  # X2' perpendicular to the fold axis and parallel to the horizontal. See
  # equation 3.10
  a[1, 1] <- sin(trd) * sin(plg)
  a[1, 2] <- cos(trd) * sin(plg)
  a[1, 3] <- cos(plg)
  a[2, 1] <- cos(trd)
  a[2, 2] <- -sin(trd)
  # a[2, 3] <- 0
  a[3, 1] <- sin(trd) * cos(plg)
  a[3, 2] <- cos(trd) * cos(plg)
  a[3, 3] <- -sin(plg)

  # The east, north, up coordinates of each point to be rotated already define
  # the coordinates of vectors. Thus we don't need to convert them to
  # direction cosines (and don't want to either because they are not unit vectors)
  # The following nested do-loops perform the coordinate transformation on the
  # bed. The details of this algorithm are described in Chapter 4

  for (nv in 1:nvtex) {
    for (i in 1:3) {
      dpbedseg[nv, i] <- 0
      for (j in 1:3) {
        dpbedseg[nv, i] <- a[i, j] * bedseg[nv, j] + dpbedseg[nv, i]
      }
    }
  }
  dpbedseg
}


#' rotates a line by performing a coordinate transformation on
#'
#' @param raz  trend of rotation axis
#' @param rdip plunge of rotation axis
#' @param rot magnitude of rotation
#' @param trd trend of the vector to be rotated
#' @param plg plunge of the vector to be rotated
#' @param axis logical. whether the line to be rotated is an axis (`TRUE`) a
#' vector (`FALSE`)
#' @note All angles are in radians
#'
#' @returns trd and plg of rotated vector
#' @export
Rotate <- function(raz, rdip, rot, trd, plg, axis = TRUE) {
  ans <- match.arg(ans0)
  # Allocate some arrays
  a <- matrix(nrow = 3, ncol = 3) # Transformation matrix
  pole <- rep(0, 3) # Direction cosines of rotation axis
  plotr <- rep(0, 3) # Direction cosines of rotated vector
  temp <- rep(0, 3) # Direction cosines of unrotated vector

  # Convert rotation axis to direction cosines. Note that the convention here
  # is X1 = North, X2 = East, X3 = Down
  p <- SphToCart(raz, rdip, 0)
  pole[1] <- p[, 1]
  pole[2] <- p[, 2]
  pole[3] <- p[, 3]

  # Calculate the transformation matrix
  x <- 1 - cos(rot)
  sinRot <- sin(rot) # Just reduces the number of calculations
  cosRot <- cos(rot)
  a[1, 1] <- cosRot + pole[1] * pole[1] * x
  a[1, 2] <- -pole[3] * sinRot + pole[1] * pole[2] * x
  a[1, 3] <- pole[2] * sinRot + pole[1] * pole[3] * x
  a[2, 1] <- pole[3] * sinRot + pole[2] * pole[1] * x
  a[2, 2] <- cosRot + pole[2] * pole[2] * x
  a[2, 3] <- -pole[1] * sinRot + pole[2] * pole[3] * x
  a[3, 1] <- -pole[2] * sinRot + pole[3] * pole[1] * x
  a[3, 2] <- pole[1] * sinRot + pole[3] * pole[2] * x
  a[3, 3] <- cosRot + pole[3] * pole[3] * x

  # Convert trend and plunge of vector to be rotated into direction cosines
  t <- SphToCart(trd, plg, 0)
  temp[1] <- t[, 1]
  temp[2] <- t[, 3]
  temp[3] <- t[, 3]

  # The following nested loops perform the coordinate transformation
  for (i in 1:3) {
    plotr[i] <- 0
    for (j in 1:3) {
      plotr[i] <- a[i, j] * temp[j] + plotr[i]
    }
  }

  # Convert to lower hemisphere projection if data are axes (ans0 = 'a')
  if (plotr[3] < 0 && axis) {
    plotr[1] <- -plotr[1]
    plotr[2] <- -plotr[2]
    plotr[3] <- -plotr[3]
  }

  # Convert from direction cosines back to trend and plunge
  r <- CartToSph(plotr[1], plotr[2], plotr[3])
  rtrd <- r[, 1]
  rplg <- r[, 2]
  cbind(rtrd, rplg)
}

#' computes the great circle path of a plane in an equal angle or equal area
#' stereonet of unit radius
#'
#' @param strike strike of plane
#' @param dip dip of plane
#' @param sttype character indicating the type of stereonet. `"eangle"` for equal angle
#' and `"earea"` for equal area
#'
#' @return
#' @export
#'
#' @examples
GreatCircle <- function(strike, dip, sttype = c("earea", "eangle")) {
  # Compute the pole to the plane. This will be the axis of rotation to make
  # the great circle
  p <- Pole(strike, dip, 1)
  trda <- p[, 1]
  plga <- p[, 2]

  # Now pick a line at the intersection of the great circle with the primitive
  # of the stereonet
  trd <- strike
  plg <- 0

  # To make the great circle, rotate the line 180 degrees in increments
  # of 1 degree
  rot <- seq(0, 180, 1) * pi / 180

  path <- matrix(0, nrow = length(rot), ncol = 2)

  for (i in 1:length(rot)) {
    # Avoid joining ends of path
    if (rot[i] == pi) {
      rot[i] <- rot[i] * 0.9999
    }

    # Rotate line
    r <- Rotate(trda, plga, rot[i], trd, plg)
    rtrd <- r[, 1]
    rplg <- r[, 2]

    # Calculate stereonet coordinates of rotated line and add to great
    # circle path

    t <- StCoordLine(rtrd, rplg, sttype)

    path[i, 1] <- t[, 1]
    path[i, 2] <- t[, 2]
    path
  }
}



#' computes the paths of a small circle defined by its axis and cone angle,
#' for an equal angle or equal area stereonet of unit radius
#'
#' @param trda trend of axis
#' @param plga plunge of axis
#' @param coneAngle cone angle
#' @param sttype character indicating the type of stereonet. `"eangle"` for equal angle
#' and `"earea"` for equal area
#'
#' @return
#' @export
#'
#' @examples
SmallCircle <- function(trda, plga, coneAngle, sttype = c("earea", "eangle")) {
  # Find where to start the small circle
  if ((plga - coneAngle) >= 0) {
    trd <- trda
    plg <- plga - coneAngle
  } else {
    if (plga == pi / 2) {
      plga <- plga * 0.9999
    }
    angle <- acos(cos(coneAngle) / cos(plga))
    trd <- ZeroTwoPi(trda + angle)
    plg <- 0
  }

  # To make the small circle, rotate the starting line 360 degrees in
  # increments of 1 degree
  rot <- seq(0, 360, 1) * pi / 180
  path1 <- path2 <- matrix(0, nrow = length(rot), ncol = 3)
  np1 <- np2 <- 0

  for (i in 1:length(rot)) {
    # Rotate line: Notice that here the line is considered as a vector
    r <- Rotate(trda, plga, rot[i], trd, plg, FALSE)
    rtrd <- r[, 1]
    rplg <- r[, 2]

    # Add to the right path
    # If plunge of rotated line is positive add to first path
    if (rplg >= 0) {
      np1 <- np1 + 1

      # Calculate stereonet coordinates and add to path
      t1 <- StCoordLine(rtrd, rplg, sttype)
      path1[np1, 1] <- t1[, 1]
      path1[np1, 2] <- t1[, 2]
    } else {
      # If plunge of rotated line is negative add to second path

      np2 <- np2 + 1
      # Calculate stereonet coordinates and add to path
      t2 <- StCoordLine(rtrd, rplg, sttype)
      path2[np2, 1] <- t2[, 1]
      path2[np2, 2] <- t2[, 2]
    }
  }
  cbind(path1, path2, np1, np2)
}



#' #' multiplies two conformable matrices
#' #'
#' #' @param a,b matrices
#' #'
#' #' @return matrix
#' MultMatrix <- function(a,b){
#'   aRow = dim(a)[1] # Number of rows in a
#'   aCol = dim(a)[2] # Number of columns in a
#'   bRow = dim(b)[1] # Number of rows in b
#'   bCol = dim(b)[2] # Number of columns in b
#'
#' # If the multiplication is conformable
#' if(aCol == bRow){
#'   # Initialize C
#' C = matrix(0, nrow = aRow, ncol = bCol)
#' for(i in 1:aRow){ # note the use of the nested loops
#'   for(j in 1:bCol){ # to do the matrix multiplication
#'     for(k in 1:aCol){
#'       C[i,j] = a[i,k]*b[k,j] + C[i,j]
#'     }
#'   }
#' }
#' } else {
#'   stop('Error: Matrices are not conformable')
#'   }
#' C
#' }

#' calculates all of the cofactor elements for a 3 x 3 matrix
#'
#' @param a matrix
CalcCofac <- function(a) {
  # Number of rows and columns in a
  n <- dim(a)[1]
  m <- dim(a)[2]

  # If matrix is 3 x 3
  if (n == 3 & m == 3) {
    # Initialize cofactor
    cofac <- matrix(0, nrow = 3, ncol = 3)

    # Calculate cofactor. When i+j is odd, the cofactor is negative
    cofac[1, 1] <- a[2, 2] * a[3, 3] - a[2, 3] * a[3, 2]
    cofac[1, 2] <- -(a[2, 1] * a[3, 3] - a[2, 3] * a[3, 1])
    cofac[1, 3] <- a[2, 1] * a[3, 2] - a[2, 2] * a[3, 1]
    cofac[2, 1] <- -(a[1, 2] * a[3, 3] - a[1, 3] * a[3, 2])
    cofac[2, 2] <- a[1, 1] * a[3, 3] - a[1, 3] * a[3, 1]
    cofac[2, 3] <- -(a[1, 1] * a[3, 2] - a[1, 2] * a[3, 1])
    cofac[3, 1] <- a[1, 2] * a[2, 3] - a[1, 3] * a[2, 2]
    cofac[3, 2] <- -(a[1, 1] * a[2, 3] - a[1, 3] * a[2, 1])
    cofac[3, 3] <- a[1, 1] * a[2, 2] - a[1, 2] * a[2, 1]
  } else {
    stop("Matrix is not 3 x 3")
  }
  cofac
}




#' Cauchy
#'
#' Given the stress tensor in a X1,X2,X3 coordinate system of any orientation,
#' Cauchy computes the X1,X2,X3 tractions on an arbitrarily oriented plane
#'
#' @param stress Symmetric 3 x 3 stress tensor
#' @param tX1,pX1 trend and plunge of X1
#' @param tX3 trend of X3
#' @param strike,dip strike and dip of plane
#'
#' @note plane orientation follows the right hand rule. Input/Output angles are
#' in radians
#'
#' @returns list. `t` = 1 x 3 vector with tractions in X1, X2 and X3.
#' `pT` = 1 x 3 vector with direction cosines of pole to plane transformed
#'  to X1,X2,X3 coordinates
Cauchy <- function(stress, tX1, pX1, tX3, strike, dip) {
  # Compute direction cosines of X1,X2,X3
  dC <- DirCosAxes(tX1, pX1, tX3)

  # Calculate direction cosines of pole to plane
  p <- matrix(nrow = 1, ncol = 3)

  cart <- SphToCart(strike, dip, 1)
  p[1] <- cart[, 1]
  p[2] <- cart[, 2]
  p[3] <- cart[, 3]

  # Transform pole to plane to stress coordinates X1,X2,X3
  # The transformation matrix is just the direction cosines of X1,X2,X3
  pT <- matrix(nrow = 1, ncol = 3)
  for (i in 1:3) {
    for (j in 1:3) {
      pT[i] <- dC[i, j] * p[j] + pT[i]
    }
  }

  # Convert transformed pole to unit vector
  r <- sqrt(pT[1] * pT[1] + pT[2] * pT[2] + pT[3] * pT[3])
  # for(i in 1:3){
  #   pT[i] = pT[i]/r
  # }
  pT <- pT / r

  # Calculate the tractions in stress coordinates X1,X2,X3
  t <- matrix(nrow = 1, ncol = 3) # Initialize t
  # Compute tractions using Cauchy's law (Eq. 6.7b)
  for (i in 1:3) {
    for (j in 1:3) {
      t[i] <- stress[i, j] * pT[j] + t[i]
    }
  }

  list(t = t, pT = pT)
}


#' direction cosines
#'
#' calculates the direction cosines of a right handed, orthogonal X1,X2,X3
#' cartesian coordinate system of any orientation with respect to North-East-Down
#'
#' @param tX1,pX1 trend and plunge of X1
#' @param tX3 trend of X3
#'
#' @returns 3 x 3 matrix containing the direction cosines of X1 (row 1),
#' X2 (row 2), and X3 (row 3)
#' @note Input angles should be in radians
DirCosAxes <- function(tX1, pX1, tX3) {
  # Some constants
  east <- pi / 2
  west <- 1.5 * pi

  # Initialize matrix of direction cosines
  dC <- matrix(nrow = 3, ncol = 3)

  # Direction cosines of X1
  cart <- SphToCart(tX1, pX1, 0)
  dC[1, 1] <- cart[, 1]
  dC[1, 2] <- cart[, 2]
  dC[1, 3] <- cart[, 3]

  # Calculate plunge of axis 3
  # If axis 1 is horizontal
  if (pX1 == 0) {
    if (abs(tX1 - tX3) == east | abs(tX1 - tX3) == west) {
      pX3 <- 0
    } else {
      pX3 <- east
    }
  } else {
    # From Equation 2.14 and with theta equal to 90 degrees
    pX3 <- atan(-(dC[1, 1] * cos(tX3) + dC[1, 2] * sin(tX3)) / dC[1, 3])
  }

  # Direction cosines of X3
  crt <- SphToCart(tX3, pX3, 0)
  dC[3, 1] <- cart[, 1]
  dC[3, 2] <- cart[, 2]
  dC[3, 3] <- cart[, 3]

  # Compute direction cosines of X2 by the cross product of X3 and X1
  dC[2, 1] <- dC[3, 2] * dC[1, 3] - dC[3, 3] * dC[1, 2]
  dC[2, 2] <- dC[3, 3] * dC[1, 1] - dC[3, 1] * dC[1, 3]
  dC[2, 3] <- dC[3, 1] * dC[1, 2] - dC[3, 2] * dC[1, 1]

  # Convert X2 to a unit vector
  r <- sqrt(dC[2, 1] * dC[2, 1] + dC[2, 2] * dC[2, 2] + dC[2, 3] * dC[2, 3])

  # for(i in 1:3){
  #   dC[2,i] = dC[2,i]/r
  # }
  dc[2, ] <- dC[2, ] / r

  dC
}
