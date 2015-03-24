#' Calculate the roughness of the vector. The roughness is the RMS
#' of the differences between consecutive points.
#'
#' @param x Values
#'
#' @export
#'
calc.roughness <- function(x) {
    N <- length(x)
    stopifnot(N > 0)
    if (0 == N) return(0)
    S <- sd(x)
    if (S == 0) return(0)
    sqrt(sum((x[1:(N-1)] - x[2:N])**2)) / S / (N-1)
}


#' Permute a data frame, x. If group.col is given it should name an ordered
#' factor that the order of the permutation should respect.
#'
#' @param .df Data frame
#' @param group.col Name of an ordered factor that the permuation
#'    should respect.
#'
#' @export
#'
permute.df <- function(.df, group.col=NULL) {
    if (is.null(group.col)) {
        return(.df[sample(nrow(.df)),])
    } else {
        stopifnot(is.ordered(.df[,group.col]))
        return(
            .df
            %>% group_by_(as.symbol(group.col))
            %>% do(permute.df(.)))
    }
}

