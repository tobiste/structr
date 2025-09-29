if (requireNamespace("rgl", quietly = TRUE)) {
  options(rgl.useNULL = TRUE)
}

# # nocov start
# #' @importFrom vctrs s3_register
# register_all_s3_methods <- function() {
#   s3_register("utils::rbind", "spherical")
#   s3_register("utils::head", "spherical")
#   s3_register("utils::tail", "spherical")
#   s3_register("stats::sd", "spherical")
#   s3_register("stats::var", "spherical")
#   # s3_register("eigen", "spherical")
#   # s3_register("mean", "spherical")
#   # s3_register("[`", "spherical")
#   # s3_register("print", "spherical")
#   s3_register("graphics::points", "spherical")
#   s3_register("graphics::text", "spherical")
#   s3_register("graphics::lines", "spherical")
#   s3_register("graphics::contour", "spherical")
#   s3_register("graphics::filled.contour", "spherical")
#   s3_register("graphics::image", "spherical")
#   s3_register("stats::density", "spherical")
#   # register_vctrs_methods()
# }
#
# .onLoad <- function(...) {
#   register_all_s3_methods()
# }

.onLoad <- function(libname, pkgname) {
  # if (!interactive())
  options(rgl.useNULL = TRUE)
}
