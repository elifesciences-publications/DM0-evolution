# Set up some useful options
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.diauxicgrowthcurver <- list(
    diauxicgrowthcurver.path = "~/diauxicgrowthcurver",
    diauxicgrowthcurver.install.args = "",
    diauxicgrowthcurver.name = "Rohan Maddamsetti",
    diauxicgrowthcurver.desc.author = '"Rohan Maddamsetti <rmaddams@odu.edu> [aut], forked from Kathleen Sprouffske <sprouffske@gmail.com> [cre]"',
    diauxicgrowthcurver.desc.license = "GPL (>=2)",
    diauxicgrowthcurver.desc.suggests = NULL,
    diauxicgrowthcurver.des = list()
  )
  toset <- !(names(op.diauxicgrowthcurver) %in% names(op))
  if(any(toset)) options(op.diauxicgrowthcurver[toset])

  invisible()
}
