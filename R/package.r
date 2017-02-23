#' DeLorean.
#'
#' @name DeLorean
#' @docType package
#'
#' @import Rcpp
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' @importFrom rstan stan sflist2stanfit vb monitor optimizing log_prob
#'   unconstrain_pars extract
#' @importFrom stats aov as.dist cmdscale cor dgamma
#'   dist dlnorm median prcomp predict runif
#'   sd t.test var weighted.mean
#' @importFrom utils head tail
#'
NULL

# To avoid R CMD check --as-cran NOTES
utils::globalVariables(names=c(
  '.',
  'Mgk',
  'Omega',
  'Psi',
  'S.hat',
  'V',
  'Vgk',
  'adjustment',
  'capture',
  'capture',
  'cell',
  'cell.expr',
  'cell.meta',
  'cell.sizes',
  'cell.type',
  'dynamic.var',
  'expr',
  'expr.l',
  'fit',
  'g',
  'gene',
  'gene.expr',
  'gene.map',
  'gene.map',
  'gene.meta',
  'gene.time.expr',
  'gene.var',
  'held.out.expr',
  'hyper',
  'initialisation',
  'is.held.out',
  'iter',
  'median.tau',
  'mu',
  'name',
  'num.capture',
  'obstime',
  'omega',
  'omega.hat',
  'omega.mean',
  'omega.sd',
  'opts',
  'opts',
  'order.inits',
  'ordering',
  'ordering.label',
  'p.value',
  'pca.df',
  'phi',
  'phi.hat',
  'predictedmean',
  'predictedvar',
  'prop.expr',
  'pseudotime',
  'pseudotime.label',
  'psi',
  'psi.hat',
  'psi.mean',
  'psi.sd',
  'sample.iter',
  'samples.l',
  'size',
  'stan.data',
  'tau',
  'term',
  'time.range',
  'type',
  'x',
  'x.hat',
  'x.hat.sd',
  'x.mean',
  'x.sd',
  'x.var',
  'z',
  'z.hat'))

.onLoad <- function(libname, pkgname) { # nocov start
  # loadRcppModules()
} # nocov end

#' Initialise DeLorean object
#'
#' @param expr Expression array
#' @param gene.meta Data frame of meta data for genes
#' @param cell.meta Data frame of meta data for cells
#'
#' @examples
#' data(WindramDeLorean)
#' dl <- de.lorean(windram.expr, windram.gene.meta, windram.cell.meta)
#'
#' @export
#'
de.lorean <- function(expr, gene.meta, cell.meta) {
    stopifnot(
        "gene" %in% names(gene.meta),
        is.factor(gene.meta$gene),
        "cell" %in% names(cell.meta),
        "capture" %in% names(cell.meta),
        "obstime" %in% names(cell.meta),
        is.factor(cell.meta$cell),
        is.factor(cell.meta$capture),
        is.numeric(cell.meta$obstime),
        is.matrix(expr),
        nrow(expr) == nrow(gene.meta),
        ncol(expr) == nrow(cell.meta),
        ! is.null(rownames(expr)),
        ! is.null(colnames(expr)),
        all(rownames(expr) == gene.meta$gene),
        all(colnames(expr) == cell.meta$cell))
    result <- list(
        expr = expr,
        gene.meta = gene.meta,
        cell.meta = cell.meta,
        opts = list())
    class(result) <- c("de.lorean", class(result))
    result
}

#' Is a DeLorean object?
#'
#' @param x de.lorean object
#'
#' @export
#'
is.de.lorean <- function(x) inherits(x, "de.lorean")

# Is internally consistent?
#
check.de.lorean.consistent <- function(dl) with(dl, {
    if (exists('gene.map')) {
        # Check gene map is in correct order
        stopifnot(all(arrange(gene.map, g)$g == 1:nrow(gene.map)))
        gene.names <- as.character(gene.map$gene)
        stopifnot(all(rownames(expr) == gene.names))
    }
    if (exists('cell.map')) {
        # Check cell map is in correct order
        stopifnot(all(arrange(cell.map, c)$c == 1:nrow(cell.map)))
        cell.names <- as.character(cell.map$cell)
        stopifnot(all(colnames(expr) == cell.names))
    }
})

#' Print details of DeLorean object
#'
#' @param x de.lorean object
#' @param ... Extra arguments
#'
#' @export
#'
print.de.lorean <- function(x, ...) {
    print(sapply(x, head))
}

#' Dimensions of DeLorean object
#'
#' @param x De lorean object
#'
#' @export
#'
dim.de.lorean <- function(x) {
    dim(x$expr)
}

#' Default number of cores to use.
#'
#' @export
#'
default.num.cores <- function() {
  if ('Windows' == Sys.info()['sysname']) {
    return(1)
  } else {
    return(getOption("DL.num.cores", max(parallel::detectCores() - 1, 1)))
  }
}

# Summarise DeLorean object
#
# @param object de.lorean object
# @param ... Extra arguments
#
# @export
#
summary.de.lorean <- function(object, ...) {
    print(sapply(x, summary))
}
