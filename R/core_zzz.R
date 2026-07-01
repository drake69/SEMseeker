.pkgglobalenv <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  assign("ssEnv", list(), envir = .pkgglobalenv)

  # Turn on R's partial-match warnings while SEMseeker is loaded. Surfaces a
  # class of silent bugs typical of dispatcher-style functions with `...`
  # forwarding: e.g. passing `areas = "GENE"` to a function with a formal
  # `areas_selection` causes R to partial-match `areas` -> `areas_selection`
  # with no diagnostic, and the runtime value ends up bound to the wrong
  # parameter. With these options on, every partial match emits a warning at
  # the call site.
  pkg_opts <- c("warnPartialMatchArgs",
                "warnPartialMatchAttr",
                "warnPartialMatchDollar")
  saved <- stats::setNames(
    lapply(pkg_opts, function(opt) getOption(opt, default = NA)),
    pkg_opts
  )
  assign("partial_match_opts_saved", saved, envir = .pkgglobalenv)

  options(warnPartialMatchArgs   = TRUE,
          warnPartialMatchAttr   = TRUE,
          warnPartialMatchDollar = TRUE)

  invisible()
}

.onUnload <- function(libpath) {
  # Restore the partial-match options as they were before .onLoad ran.
  # NA in saved means "the option was unset" (R's default is unset = off);
  # convert that to NULL when calling options().
  saved <- .pkgglobalenv$partial_match_opts_saved
  if (is.null(saved)) return(invisible())
  for (opt in names(saved)) {
    val <- saved[[opt]]
    args <- list()
    args[[opt]] <- if (identical(val, NA)) NULL else val
    do.call(options, args)
  }
  invisible()
}
