% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/math.R
\name{strain_shape}
\alias{strain_shape}
\alias{principal_stretch}
\alias{principal_strain}
\alias{or_shape_params}
\title{Principal Stretches, Strain and Shape Parameters based on the Orientation Tensor.}
\usage{
principal_stretch(x)

principal_strain(x)

or_shape_params(x)
}
\arguments{
\item{x}{numeric. Can be three element vector, three column array, or an
object of class \code{"line"} or \code{"plane"}}
}
\value{
list
}
\description{
Principal Stretches, Strain and Shape Parameters based on the Orientation Tensor.
}
\details{
\describe{
\item{\code{stretch_ratios}}{Sqrt of eigenvalue ratios}
\item{\code{strain_ratios}}{Log of stretch ratios}
\item{\code{Ramsay}}{strain symmetry (Ramsay, 1983)}
\item{\code{Woodcock}}{Woodcock shape}
\item{\code{Flinn}}{Flinn strain intensity}
\item{\code{Vollmer}}{Point, Girdle, Random, Cylindricity (B), and Uniform Distance (D) Indices (Vollmer 1990; 2020). \code{D} is a measure of the "distance" from uniformity, and is linear from R to P, and R to G. End members are: uniform D = 0, girdle D = 0.5, cluster D = 1. The 99\% level for a test against uniformity for a sample size of 300 is D = 0.1.}
\item{\code{Nadai}}{natural octahedral unit strain and shear (Nadai, 1963)}
\item{\code{Lisle_intensity}}{Intensity index (Lisle, 1985)}
\item{\code{Waterson_intensity}}{strain intensity (Watterson, 1968)}
\item{\code{lode}}{Lode parameter (Lode, 1926)}
\item{\code{kind}}{Descriptive type of ellipsoid}
\item{\code{MAD}}{maximum angular deviation (Kirschvink, 1980)}
\item{\code{US}}{Uniformity statistic of Mardia (1972)}
}
}
\examples{
set.seed(1)
mu <- Line(120, 50)
x <- rvmf(100, mu = mu, k = 20)
principal_stretch(x)
principal_strain(x)
or_shape_params(x)
}
\references{
Flinn, Derek.(1963): "On the statistical analysis of fabric diagrams." Geological Journal 3.2: 247-253.

Kirschvink, J. (1980): The least-squares line and plane and the analysis of palaeomagnetic data. Geophysical Journal International, 62(3), 699-718.

Lisle, Richard J.  (1985): "The use of the orientation tensor for the description and statistical testing of fabrics." Journal of Structural Geology 7.1: 115-117.

Lode, Walter (1926): "Versuche über den Einfluß der mittleren Hauptspannung auf das Fließen der Metalle Eisen, Kupfer und Nickel“
(\emph{"Experiments on the influence of the mean principal stress on the flow of the metals iron, copper and nickel"}], Zeitschrift für Physik, vol. 36 (November), pp. 913–939, DOI: 10.1007/BF01400222

Mardia, Kantilal Varichand. (1975): "Statistics of directional data." Journal of the Royal Statistical Society Series B: Statistical Methodology 37.3: 349-371.

Nadai, A., and Hodge, P. G., Jr. (1963): "Theory of Flow and Fracture of Solids, vol. II." ASME. J. Appl. Mech. December 1963; 30(4): 640. https://doi.org/10.1115/1.3636654

Ramsay, John G. (1967): "Folding and fracturing of rocks." Mc Graw Hill Book Company 568.

Vollmer, Frederick W. (1990): "An application of eigenvalue methods to structural domain analysis." Geological Society of America Bulletin 102.6: 786-791.

Vollmer, Frederick W. (2020): "Representing Progressive Fabric Paths on a Triangular Plot Using a Fabric Density Index and Crystal Axes Eigenvector Barycenters." Geological Society of America Abstracts. Vol. 52.

Watterson, Juan. (1968): "Homogeneous deformation of the gneisses of Vesterland, south-west Greenland". No. 78. CA Reitzel.

Woodcock, N. H.  (1977): "Specification of fabric shapes using an eigenvalue method." Geological Society of America Bulletin 88.9: 1231-1236.
}
\seealso{
\code{\link[=ortensor]{ortensor()}}, \code{\link[=or_eigen]{or_eigen()}}, \code{\link[=fabric_indexes]{fabric_indexes()}}
}
