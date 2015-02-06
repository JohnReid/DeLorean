#' DeLorean.
#'
#' @name DeLorean
#' @docType package
NULL

#' Initialise DeLorean object
#'
#' @export
#'
de.lorean <- function(expr, gene.meta, cell.meta) {
    stopifnot("gene" %in% names(gene.meta))
    stopifnot(is.factor(gene.meta$gene))
    stopifnot("cell" %in% names(cell.meta))
    stopifnot("capture" %in% names(cell.meta))
    stopifnot("obstime" %in% names(cell.meta))
    stopifnot(is.factor(cell.meta$cell))
    stopifnot(is.factor(cell.meta$capture))
    stopifnot(is.numeric(cell.meta$obstime))
    stopifnot(nrow(expr) == nrow(gene.meta))
    stopifnot(ncol(expr) == nrow(cell.meta))
    stopifnot(! is.null(rownames(expr)))
    stopifnot(! is.null(colnames(expr)))
    stopifnot(all(rownames(expr) == gene.meta$gene))
    stopifnot(all(colnames(expr) == cell.meta$cell))
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
#' @export
#'
is.de.lorean <- function(dl) inherits(dl, "de.lorean")

#' Print details of DeLorean object
#'
#' @export
#'
print.de.lorean <- function(dl) {
    print(sapply(dl, head))
}

#' Dimensions of DeLorean object
#'
#' @export
#'
dim.de.lorean <- function(dl) {
    dim(dl$expr)
}

#' Summarise DeLorean object
#'
#' @param dl de.lorean object
#'
#' @export
#'
summary.de.lorean <- function(dl) {
    print(sapply(dl, summary))
}

