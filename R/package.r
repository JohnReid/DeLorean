#' DeLorean.
#'
#' @name DeLorean
#' @docType package
NULL

#' Initialise DeLorean object
#'
#' @param expr Expression array
#' @param gene.meta Data frame of meta data for genes
#' @param cell.meta Data frame of meta data for cells
#'
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' @import rstan
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

# To avoid R CMD check --as-cran NOTES
globalVariables(c(
  'cell.meta',
  'cell',
  'capture',
  'x',
  'size',
  'tau',
  'predictedmean',
  'predictedvar',
  'adjustment',
  'samples.l',
  'gene.map',
  'g',
  'omega',
  'psi',
  'phi',
  'expr',
  'gene.meta',
  'gene',
  'x.mean',
  'S.hat',
  'x.hat',
  'x.sd',
  'x.hat.sd',
  'capture',
  'x.var',
  'Vgc',
  'Mgc',
  'num.capture',
  'Psi',
  'Omega',
  'mu',
  'stan.data',
  'prop.expr',
  '.',
  'gene.time.expr',
  'psi.hat',
  'omega.hat',
  'cell.expr',
  'gene.expr',
  'gene.var',
  'opts',
  'fit',
  'phi.hat',
  'obstime',
  'time.range',
  'hyper',
  'dynamic.var',
  'opts',
  'name',
  'gene.map',
  'is.held.out',
  'sample.iter',
  'held.out.expr',
  'type',
  'iter',
  'V',
  'median.tau'))
