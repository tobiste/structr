# Add an Ellipse to existing plot

Add an Ellipse to existing plot

## Usage

``` r
ellipse(
  x = 0,
  y = x,
  radius.x = 1,
  radius.y = radius.x,
  rot = 0,
  nv = 512,
  border = par("fg"),
  col = par("bg"),
  lty = par("lty"),
  lwd = par("lwd"),
  plot = TRUE
)
```

## Arguments

- x, y:

  the x and y co-ordinates for the center(s) of the ellipse(s).

- radius.x:

  scalar or a vector giving the semi-major axis of the ellipse.

- radius.y:

  a scalar or a vector giving the semi-minor axis of the ellipse.

- rot:

  angle of rotation in radians.

- nv:

  number of vertices to draw the ellipses.

- border:

  color for borders. The default is par("fg"). Use border = NA to omit
  borders.

- col:

  color(s) to fill or shade the annulus sector with. The default NA (or
  also NULL) means do not fill (say draw transparent).

- lty:

  line type for borders and shading; defaults to "solid".

- lwd:

  line width for borders and shading.

- plot:

  logical. If TRUE the structure will be plotted. If FALSE only the
  points are calculated and returned. Use this if you want to combine
  several geometric structures to a single polygon.

## Value

The function invisibly returns a list of the calculated coordinates for
all shapes.

## Examples

``` r
plot(c(0, 1), c(0, 1), type = "n")
ellipse(.5, .5, radius.x = 0.5, radius.y = .25, col = "darkgreen", border = "red")
```
