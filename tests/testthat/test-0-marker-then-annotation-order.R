# Meta-test: genomic annotation must run AFTER derived-marker computation.
#
# Invariant (richiesto esplicitamente): le annotazioni genomiche (Body, TSS,
# Island shores/shelves, ...) prodotte da annotate_position_pivots() devono
# essere SEMPRE successive al calcolo dei marker derivati DELTAQ/DELTARQ/
# DELTAP/DELTARP fatto da deltaX_get(). Inoltre, dove i POSITION pivot di base
# vengono materializzati con create_position_pivots(), quella chiamata deve
# precedere deltaX_get() (deltaX_get legge i pivot DELTAS/DELTAR base).
#
# annotate_position_pivots() NON ricalcola i marker: legge il POSITION/WHOLE
# pivot di ciascun marker (inclusi i derivati) e lo proietta su AREA/SUBAREA.
# Se qualcuno invertisse l'ordine, le annotazioni dei marker derivati
# leggerebbero pivot inesistenti o stantii. Questo test sorgente blocca
# l'inversione a tempo di test, prima della pipeline.

# Walk del parse-tree del corpo di una funzione: restituisce, per ciascun
# nome-di-chiamata richiesto, l'indice di prima comparsa in ordine sorgente
# (NA se assente). L'indice è un contatore globale incrementato a ogni nodo,
# quindi rispetta l'ordinamento testuale degli statement.
.first_call_orders <- function(fun_body, call_names) {
  counter <- 0L
  first <- stats::setNames(rep(NA_integer_, length(call_names)), call_names)
  walk <- function(node) {
    counter <<- counter + 1L
    if (is.call(node)) {
      head <- node[[1L]]
      if (is.symbol(head)) {
        nm <- as.character(head)
        if (nm %in% call_names && is.na(first[[nm]])) {
          first[[nm]] <<- counter
        }
      }
      for (i in seq_along(node)) {
        el <- node[[i]]
        if (!missing(el) && (is.call(el) || is.pairlist(el))) walk(el)
      }
    } else if (is.pairlist(node)) {
      for (i in seq_along(node)) {
        el <- node[[i]]
        if (!missing(el) && (is.call(el) || is.pairlist(el))) walk(el)
      }
    }
    invisible(NULL)
  }
  walk(fun_body)
  first
}

.extract_fun_body <- function(r_file, fun_name) {
  exprs <- parse(r_file)
  for (e in as.list(exprs)) {
    if (length(e) >= 3L && is.symbol(e[[1L]]) &&
        as.character(e[[1L]]) %in% c("<-", "=", "<<-") &&
        is.symbol(e[[2L]]) && as.character(e[[2L]]) == fun_name &&
        is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.symbol("function"))) {
      return(e[[3L]])  # the function() expression (formals + body)
    }
  }
  NULL
}

test_that("deltaX_get() precedes annotate_position_pivots() in every caller", {
  pkg_root <- testthat::test_path("..", "..")
  r_dir <- file.path(pkg_root, "R")

  # caller file -> function name
  callers <- list(
    c("semseeker_core.R",     "semseeker_core"),
    c("recover.R",            "recover"),
    c("association_analysis.R", "association_analysis")
  )

  for (caller in callers) {
    r_file <- file.path(r_dir, caller[[1L]])
    expect_true(file.exists(r_file), info = paste("missing source:", caller[[1L]]))

    fun <- .extract_fun_body(r_file, caller[[2L]])
    expect_false(is.null(fun),
      info = paste("could not locate function", caller[[2L]], "in", caller[[1L]]))

    ord <- .first_call_orders(
      fun,
      c("create_position_pivots", "deltaX_get", "annotate_position_pivots")
    )

    # deltaX_get() and annotate_position_pivots() must BOTH appear and be ordered.
    expect_false(is.na(ord[["deltaX_get"]]),
      info = paste(caller[[2L]], "must call deltaX_get()"))
    expect_false(is.na(ord[["annotate_position_pivots"]]),
      info = paste(caller[[2L]], "must call annotate_position_pivots()"))
    expect_lt(ord[["deltaX_get"]], ord[["annotate_position_pivots"]])

    # Where create_position_pivots() is present, it must precede deltaX_get()
    # (deltaX_get reads the base DELTAS/DELTAR POSITION pivots it writes).
    if (!is.na(ord[["create_position_pivots"]])) {
      expect_lt(ord[["create_position_pivots"]], ord[["deltaX_get"]])
    }
  }
})
