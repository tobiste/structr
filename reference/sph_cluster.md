# Cluster spherical data

Finds k groups of clusters using the angular distance matrix

## Usage

``` r
sph_cluster(
  x,
  k,
  method = c("hclust", "kmeans", "diana", "agnes", "pam", "clara", "fanny"),
  ...
)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

- k:

  integer. Number of desired clusters.

- method:

  character. Clustering method to be applied. Currently implemented are

  `"hclust"`

  :   Hierarchical Clustering using
      [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html), the
      default)

  `"kmeans"`

  :   K-Means Clustering using
      [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html))

  `"pam"`

  :   Partitioning Around Medoids using
      [`cluster::pam()`](https://rdrr.io/pkg/cluster/man/pam.html)

  `"agnes"`

  :   Agglomerative hierarchical clustering using
      [`cluster::agnes()`](https://rdrr.io/pkg/cluster/man/agnes.html)

  `"diana"`

  :   Divisive hierarchical clustering using
      [`cluster::diana()`](https://rdrr.io/pkg/cluster/man/diana.html)

  `"clara"`

  :   Clustering Large Applications using
      [`cluster::clara()`](https://rdrr.io/pkg/cluster/man/clara.html)

  `"fanny"`

  :   Fuzzy Analysis Clustering using
      [`cluster::fanny()`](https://rdrr.io/pkg/cluster/man/fanny.html)

- ...:

  optional arguments passed to cluster algorithm.

## Value

output of applied cluster function

## See also

[`dist()`](https://rdrr.io/r/stats/dist.html)

## Examples

``` r
set.seed(20250411)
x1 <- rvmf(100, mu = Line(90, 0), k = 20)
x2 <- rvmf(100, mu = Line(0, 0), k = 20)
x3 <- rvmf(100, mu = Line(0, 90), k = 20)
x123 <- rbind(x1, x2, x3)
cl <- sph_cluster(x123, k = 3)
plot(x123, col = cl$cluster)


# sph_cluster(simongomez, k = 3)
```
