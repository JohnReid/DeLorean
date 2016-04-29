#' Estimate the cell sizes according to Anders & Huber
#' Differential expression analysis for sequence count data
#'
#' @param expr.l Melted expression values.
#'
anders.huber.cell.sizes <- function(expr.l) {
    (
        expr.l
        %>% group_by(gene)
        %>% dplyr::summarise(mu=mean(x))
        %>% left_join(expr.l)
        %>% group_by(cell)
        %>% dplyr::summarise(S.hat=median(x-mu))
    )
}


#' Estimate the cell sizes per capture. Only uses genes that are expressed
#' in more than half the cells.
#'
#' @param expr.l Melted expression values.
#'
#' @export
#'
estimate.capture.cell.sizes <- function(expr.l) (
    expr.l
    %>% group_by(gene)
    %>% dplyr::summarise(prop.expr=mean(x>0))
    %>% filter(prop.expr > .5)
    %>% left_join(expr.l)
    %>% group_by(capture)
    %>% do(anders.huber.cell.sizes(.))
)


#' Estimate the cell sizes. We only consider genes that are expressed in
#' a certain proportion of cells.
#'
#' @param dl de.lorean object.
#' @param cell.prop The proportion of cells a gene must be expressed in to be
#'   considered for cell size estimation
#' @param expr.threshold The threshold at which we consider a gene to be expressed
#'
#' @export
#'
estimate.cell.sizes <- function(dl, cell.prop=.5, expr.threshold=0) within(dl, {
    expr.l <- melt.expr(dl)
    cell.sizes <- anders.huber.cell.sizes(
      expr.l %>%
      group_by(gene) %>%
      dplyr::summarise(prop.expr=mean(x>expr.threshold)) %>%
      filter(prop.expr > cell.prop) %>%
      left_join(expr.l))
})


#' Adjust the expression by the estimated cell sizes.
#'
#' @param dl de.lorean object.
#'
#' @export
#'
adjust.by.cell.sizes <- function(dl) within(dl, {
    expr.before.adj <- expr
    expr <- cast.expr(
        expr.l
        %>% left_join(cell.sizes)
        %>% mutate(x=x-S.hat))
})
