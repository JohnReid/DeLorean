#' Estimate the cell sizes according to Anders & Huber
#' Differential expression analysis for sequence count data
#'
#' @param x.l Melted expression values.
#'
#' @export
#'
anders.huber.cell.sizes <- function(x.l) {
    (
        x.l
        %>% group_by(gene)
        %>% summarise(mu=mean(x))
        %>% left_join(x.l)
        %>% group_by(cell)
        %>% summarise(size=median(x-mu))
    )
}
