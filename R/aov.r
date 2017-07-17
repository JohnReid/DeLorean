#' Perform an analysis of variance to select genes for the DeLorean model.
#'
#' @param dl DeLorean object.
#' @param expr.l Melted (with melt.expr()) expression values to use
#'
#' @export
#'
aov.dl <- function(dl, expr.l = melt.expr(dl)) {
  dl$aov <-
    expr.l %>%
    dplyr::left_join(dl$cell.meta) %>%
    dplyr::group_by(gene) %>%
    dplyr::do(broom::tidy(aov(x ~ capture, .))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(term == 'capture') %>%
    dplyr::arrange(p.value)
  dl
}
