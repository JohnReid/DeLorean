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
    stopifnot(
        "gene" %in% names(gene.meta),
        is.factor(gene.meta$gene),
        "cell" %in% names(cell.meta),
        "capture" %in% names(cell.meta),
        "obstime" %in% names(cell.meta),
        is.factor(cell.meta$cell),
        is.factor(cell.meta$capture),
        is.numeric(cell.meta$obstime),
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

