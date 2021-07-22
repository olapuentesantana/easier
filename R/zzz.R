#' @import ExperimentHub
#' @importFrom AnnotationHub query
.onLoad <- function(libname, pkgname) {
  # Retrieve internal data
  eh <- ExperimentHub::ExperimentHub()
  easierdata_eh <- AnnotationHub::query(eh, c("easierData"))
}
