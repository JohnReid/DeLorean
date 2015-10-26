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
#' @export
#'
is.de.lorean <- function(dl) inherits(dl, "de.lorean")

#' Is internally consistent?
#'
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

