# Import Data

``` r
library(structr)
```

## General

Data can be imported with already in R implemented functions such as
[`read.table()`](https://rdrr.io/r/utils/read.table.html),
[`read.csv()`](https://rdrr.io/r/utils/read.table.html) or from other
packages functions (e.g. [{readr}](https://readr.tidyverse.org/) and
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
> measurements in {structr}. To convert from right-hand-rule strike
> measurements to dip directions, you can simply use the
> [`rhr2dd()`](https://tobiste.github.io/structr/reference/rhr.md)
> function for that column:

``` r
imported_data[, 1] <- rhr2dd(imported_data[, 1])
as.Plane(imported_data)
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
