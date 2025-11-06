# Import orientation data from *StraboSpot*

Reads the `XLS` format export of field book data, the `JSON` project
file, or the `txt` export of StraboMobile data from
[strabospot.org/my_data](https://tobiste.github.io/structr/reference/strabospot.org/my_data)
and creates a list with the metadata, and the line or plane
orientations.

## Usage

``` r
read_strabo_xls(file, tag_cols = FALSE, sf = TRUE)

read_strabo_mobile(file, sf = TRUE)

read_strabo_JSON(file, sf = TRUE)
```

## Arguments

- file:

  the name of the file which the data are to be read from.

- tag_cols:

  logical. Whether the Tag columns should be summarized in a single
  column (may lead to duplicate rows).

- sf:

  logical. Whether the output should be a spatial `"sf"` object using
  the Longitude and Latitude columns.

## Value

`list` containing the following objects:

- `data`:

  `"tbl_df"` object. Metadata.

- `spots`:

  `"tbl_df"` object or `"sf"` if `sf == TRUE`. Locations and spot
  descriptions.

- `tags`:

  `"tbl_df"` object. Tags and their descriptions.

- `planar`:

  Plane elements. Same row IDs as in `data`.

- `linear`:

  Line elements. Same row IDs as in `data`.

## See also

[`drillcore_transformation()`](https://tobiste.github.io/structr/reference/drillcore.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# import from excel file
read_strabo_xls("path/to/my/file.xlsx")

# import from text file
read_strabo_mobile("path/to/my/file.txt")

# import from .json file
read_strabo_JSON("path/to/my/file.json")
} # }
```
