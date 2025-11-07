# Import Data

``` r
library(structr)
```

## General

Data can be imported with already in R implemented functions such as
[`read.table()`](https://rdrr.io/r/utils/read.table.html),
[`read.csv()`](https://rdrr.io/r/utils/read.table.html) or from other
packages functions (e.g. [{readr}](https://readr.tidyverse.org/) and
[{data.table}](https://cran.r-project.org/web/packages/data.table/index.html)).
Just make sure that the import produces a matrix or data.frame like
object with each measurement stored in a row, and columns representing
dip directions, dip angles, plunge, etc.

For example your `.txt` file is tab-separated and may look like this:

|     | Dip direction | Dip |
|-----|---------------|-----|
| 1   | 120           | 50  |
| 2   | 60            | 12  |
| 3   | 287           | 82  |

you could import the file like

``` r
imported_data <- read.table(
  "path/to/my/file.xt",
  header = TRUE,
  sep = "",
  row.names = 1
)
```

This imports the tab-separated file as an `"data.frame"` object with the
first column representing the row names, and the two other columns the
headers of the table. Say these measurements represent plane
measurements (e.g. bedding or fault plane orientation), we just need to
coerce that data.frame into a `"Plane"` object:

``` r
my_planes <- as.Plane(imported_data)
```

If you want a `"Line"` object, coerce the data.frame using
[`as.Line()`](https://tobiste.github.io/structr/reference/classes.md).
For `"Pair"` (line-on-plane) and `"Fault"` (line-on-plane with sense of
motion) objects, you’ll need a four and five-column table, respectively,
representing dip directions, dip angle, trend, and plunge and sense
measurements. Then use
[`as.Pair()`](https://tobiste.github.io/structr/reference/classes.md) or
[`as.Fault()`](https://tobiste.github.io/structr/reference/classes.md)
to parse the object into {structr}.

> Note that the **dip direction** is the preferred notation for plane
> measurements in {structr} as it undoubtly indicates the orientation by
> using only 2 parameters.

### Some helpers

There are several other conventions for the notation of plane and line
orientations from compass measurements. Here are some functions for you
to correctly convert your measurements into the notation required for
the use in {structr}.

#### Right-hand rule

To converts strike measurements into dip directions using right-hand
rule:

``` r
strike_measurements <- c(270, 315, 0, 45, 90, 135, 180, 225, 270)
rhr2dd(strike_measurements)
#> [1]   0  45  90 135 180 225 270 315   0

# or from dip direction to strike (using right-hand rule):
dip_directions <- c(0, 45, 90, 135, 180, 225, 270, 315, 360)
dd2rhr(dip_directions)
#> [1] 270 315   0  45  90 135 180 225 270
```

#### Quadrant notation

Sometimes, strike notation doesn’t follow the right hand rule. Instead
the quadrant of the dip direction is indicated by a Cardinal letter,
i.e. “N”, “E”, “S”, and “W”. In that case you can use the function
[`quadrant2dd()`](https://tobiste.github.io/structr/reference/quadrant2dd.md)

``` r
strike_direction <- c(270, 315, 0, 45, 90, 135, 180, 225, 270) # strike in left-hand-rule
dip_quadrtant <- c("N", "E", "E", "S", "S", "W", "W", "N", "N") # dip quadrant
quadrant2dd(strike_direction, dip_quadrtant)
#> [1]   0  45  90 135 180 225 270 315   0
```

If your table contains the strike and dip-quadrant measurement in a
single column, e.g. “270N”, you can conveniently split the column into
two by using
[`split_strike()`](https://tobiste.github.io/structr/reference/split.md)

``` r
split_strike("270N")
#> $measurement
#> 270N 
#>  270 
#> 
#> $direction
#> 270N 
#>  "N"
```

#### Fault notation

Fault (and Pair) measurements are usually a combination of Plane and Ray
or Line measurements. Then they can be defined using the
[`Fault()`](https://tobiste.github.io/structr/reference/classes.md) and
[`Pair()`](https://tobiste.github.io/structr/reference/classes.md)
functions. However, sometimes the Line or Ray component is given by the
rake angle, which is the angle between the fault strike and the
lineation. Unfortunately, there are several different ways how to
indicate the proper orientation of the lineation.

##### Fault plane and rake (or pitch)

This is the standard notation. Here, the rake is the angle between
lineation and the right-handrule strike of the fault plane. The angle is
measured on the fault plane, clockwise from the strike, where
down-plunging is positive. Rake values range between 0 and 360° (or
−180° and 180°).

If your datra follows this convention, use the function
[`Fault_from_rake()`](https://tobiste.github.io/structr/reference/fault_from_rake.md)

``` r
fault_plane <- Plane(c(120, 120, 100, 0), c(60, 60, 50, 40))
fault_pitch  <- c(84.7202, -10, 30, 180)
Fault_from_rake(fault_plane, rake = fault_pitch)
#> Fault object (n = 4):
#>      dip_direction dip   azimuth       plunge sense
#> [1,]           120  60 109.52858 5.958159e+01     1
#> [2,]           120  60 204.96163 8.649165e+00    -1
#> [3,]           100  50  30.36057 2.252101e+01     1
#> [4,]             0  40  90.00000 1.487542e-14     1
```

##### Fault plane, rake angle and plunge quadrant

Here, the rake angle is measured in the fault plane between the strike
given by either right or left-hand rule and the lineation. The angle is
recorded in a clockwise sense (looking down upon the fault plane) and
has a range from 0 to 180%deg;. The quadrant of plunge indicates the
direction of the strike from which the rake angle is measured.

If this is the notation used, call the function
[`Fault_from_rake_quadrant()`](https://tobiste.github.io/structr/reference/Fault_from_rake_quadrant.md)
and set `type="plunge"`

``` r
dip <- c(5, 10, 15, 30, 40, 55, 65, 75, 90)
dip_dir <- c(180, 225, 270, 315, 360, 0, 45, 90, 135)
rake1 <- c(0, 45, 90, 135, 180, 45, 90, 135, 180)
plunge_quadrant <- c("E", "S", "W", "N", "E", "W", "E", "S", "W")
Fault_from_rake_quadrant(Plane(dip_dir, dip), rake1, plunge_quadrant, type = "plunge")
#> Fault object (n = 9):
#>       dip_direction dip    azimuth       plunge sense
#>  [1,]           180   5  90.000000 0.000000e+00     1
#>  [2,]           225  10 179.561451 7.053022e+00     1
#>  [3,]           270  15 270.000000 1.500000e+01     1
#>  [4,]           315  30   4.106605 2.070481e+01     1
#>  [5,]             0  40  90.000000 1.487542e-14     1
#>  [6,]             0  55 299.837566 3.539626e+01     1
#>  [7,]            45  65  45.000000 6.500000e+01     1
#>  [8,]            90  75 165.489181 4.307952e+01     1
#>  [9,]           135  90 225.000000 7.016709e-15     1
```

##### Fault plane, rake angle and rake quadrant

Here, the rake is the **acute** angle measured in the fault plane
between the strike of the fault and the lineation. Starting from the
strike line, the angle is measured in a sense which is down the dip of
the plane. Quadrant of rake indicate the direction of the strike from
which the rake angle is measured, i.e. whether right-hand or left-hand
rule is followed. Angle ranges from 0 to 90 °. Use `sense` argument to
specify the sense of motion.

If this is the notation used, call the function
[`Fault_from_rake_quadrant()`](https://tobiste.github.io/structr/reference/Fault_from_rake_quadrant.md)
and set `type="rake"`

``` r
rake2 <- c(0, 45, 90, 45, 0, 45, 90, 45, 0)
rake_quadrant <- c("E", "S", "S", "E", "E", "W", "N", "S", "W")
Fault_from_rake_quadrant(Plane(dip_dir, dip), rake2, rake_quadrant, type = "rake")
#> Fault object (n = 9):
#>       dip_direction dip    azimuth    plunge sense
#>  [1,]           180   5  90.000000  0.000000     1
#>  [2,]           225  10 179.561451  7.053022     1
#>  [3,]           270  15 270.000000 15.000000     1
#>  [4,]           315  30   4.106605 20.704811     1
#>  [5,]             0  40 270.000000  0.000000     1
#>  [6,]             0  55 299.837566 35.396260     1
#>  [7,]            45  65  45.000000 65.000000     1
#>  [8,]            90  75 165.489181 43.079517     1
#>  [9,]           135  90  45.000000  0.000000     1
```

## StraboSpot

The package {structr} can import all the collected field data from your
[Strabospot](https://strabospot.org) project.

The best way is to import the .json file database of the StraboSpot
project. Go to your field data **My StraboField Data** \> scroll down to
your project \> click on **Options…** \> **Download Project in Strabo
JSON Format**

Now you can import the downloaded file via
[`read_strabo_JSON()`](https://tobiste.github.io/structr/reference/strabo.md):

``` r
strabo_data <- read_strabo_JSON("path/to/my/file.json")
```

The import function produces a `list` object with all the meta data
(`data`), the geographic locations (`spots`), the used tags (`tags`),
and all the plane (`planes`) and line measurements (`lines`) already
converted into {structr} data objects.

Additional to descriptions and comments, the meta data (`data`) also
contains information on the dataset name and time of measurement.

**IMPORTANT**: The meta data and the plane and line measurements all
share the same row indices. Thus, planes and lines with identical row
indices have been measured simultaneously (e.g. as a fault).

> This import allows that the connection of simultaneously measured
> plane and lines (such as faults and their striae) will be preserved.
> Unfortunately, if you **export** your StraboSpot field data into a
> `.csv` or `.xls` file, this connection is lost…

Alternatively, the function
[`read_strabo_xls()`](https://tobiste.github.io/structr/reference/strabo.md)
and
[`read_strabo_mobile()`](https://tobiste.github.io/structr/reference/strabo.md)
provide import of `.xls` and any character-separated table files
(e.g. `.csv` or `.txt`).

Keep in mind that these import options do not properly identify whether
your lineation measurements are ray-like or line-like vectors, nor does
it combine simultaneously lineation-plane measurements to either pairs
or faults. Thus, you’ll need to carefully convert these data dataypes
after the import to move forward.

## Drill core data

Orientations in drill-cores are usually given by α and β angles
(lineations on a plane additionally have a γ angle) which describe
orientations with respect to the drill orientation. To convert these
angles from the “drillcore coordinate reference system” to our
geographical reference system, you may use the function
[`drillcore_transformation()`](https://tobiste.github.io/structr/reference/drillcore.md).
Learn more about it in this
[tutorial](https://tobiste.github.io/structr/articles/Oriented_Drill_Cores.html).
