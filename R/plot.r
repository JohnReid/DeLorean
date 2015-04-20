#' Various DeLorean object plots
#'
#' @param dl de.lorean object
#'
#' @export
#'
plot.de.lorean <- function(dl, type="profiles", ...) {
    result <- switch(type,
        profiles=plot.profiles(dl, ...),
        S.posteriors=plot.S.posteriors(dl, ...),
        pseudotime=plot.pseudotime(dl, ...),
        convergence=plot.convergence(dl, ...),
        expr.data=plot.expr.data(dl, ...),
        roughnesses=plot.roughnesses(dl, ...),
        tau.offsets=plot.tau.offsets(dl, ...)
    )
    if (is.null(result)) {
        stop('Unknown plot type')
    }
    result
}

#' Calculate a suitable value for a rug plot given the
#' number of points
#' @param dl de.lorean object
#'
#' @export
#'
alpha.for.rug <- function(n, scale=100) {
    1 / (max(1, n / scale))
}

#' Plot posterior for marginal log likelihoods of expression profiles
#'
#' @param dl de.lorean object
#'
#' @export
#'
plot.marg.like <- function(dl) {
    with(dl, {
        gp <- (ggplot(samples.l$logmarglike %>% left_join(gene.map),
                      aes(x=gene,
                          y=logmarglike,
                          colour=is.held.out),
                      environment=environment())
            + geom_boxplot()
        )
    })
}

#' Plot pseudotime (tau) against observed capture time
#'
#' @param dl de.lorean object
#'
#' @export
#'
plot.pseudotime <- function(dl) {
    with(dl, {
        gp <- (ggplot(samples.l$tau %>% filter(iter == best.sample),
                      aes(x=tau, y=obstime, color=capture),
                      environment=environment())
            + geom_point()
            + scale_x_continuous(name="Pseudotime (tau)")
            + scale_y_continuous(name="Observed (capture) time")
        )
    })
}


#' Plot the tau offsets.
#'
#' @param dl de.lorean object
#'
#' @export
#'
plot.tau.offsets <- function(dl, rug.alpha=.3) {
    with(dl,
         ggplot(samples.l$tau, aes(x=tau.offset, color=capture))
         + geom_density()
         + geom_rug(alpha=rug.alpha)
         + stat_function(fun=Curry(dnorm, sd=hyper$sigma_tau),
                         linetype='dashed',
                         color="black")
    )
}


#' Plot the Rhat convergence statistics. examine.convergence must be called
#' before this plot can be made.
#'
#' @param dl de.lorean object
#'
#' @export
#'
plot.convergence <- function(dl) {
    with(dl, {
        rhat.df <- data.frame(
            rhat=rhat.sorted,
            param=names(rhat.sorted),
            parameter=str_match(names(rhat.sorted), "^[[:alpha:]]+"))
        gp <- (ggplot(rhat.df,
                      aes(y=rhat, x=parameter),
                      environment=environment())
            + geom_boxplot()
        )
    })
}

#' Plot posterior for cell sizes
#'
#' @param dl de.lorean object
#'
#' @export
#'
plot.S.posteriors <- function(dl) {
    with(dl, {
        gp <- (ggplot(samples.l$S
                      %>% left_join(cell.map)
                      %>% arrange(capture),
                      aes(x=factor(cell,
                                   levels=arrange(cell.map, capture)$cell),
                          y=S,
                          color=capture),
                      environment=environment())
            + geom_boxplot()
            + coord_flip()
        )
    })
}

#' Plot a comparison of the profiles from several de.lorean objects
#'
#' @param ... Named de.lorean objects
#' @param genes Genes to plot (defaults to genes.high.psi of first de.lorean
#'   object)
#' @param sample.iter Which sample to plot
#'
#' @export
#'
plot.cmp.profiles <- function(...,
                              genes = NULL) {
    dls <- list(...)
    dl.levels <- names(dls)
    stopifnot(! is.null(dl.levels))  # Must have names for de.lorean objects
    if (is.null(genes)) {
        genes <- dls[[1]]$genes.high.psi
    }
    get.mean <- function(.name) {
        with(dls[[.name]], (
            predictions
            %>% filter(best.sample == iter)
            %>% left_join(gene.map)
            %>% filter(gene %in% genes)
            %>% mutate(name=factor(.name, levels=dl.levels))
        ))
    }
    means <- do.call(rbind, lapply(dl.levels, get.mean))
    gp <- ggplot(mutate.profile.data(means),
                 aes(x=tau),
                 environment=environment())
    line.alpha <- .8
    ribbon.alpha <- .2
    (
        gp
        + geom_line(aes(x=x, y=mean, color=name),
                    alpha=line.alpha)
        + geom_ribbon(aes(x=x,
                          ymin=mean-2*sqrt(var),
                          ymax=mean+2*sqrt(var),
                          fill=name),
                      alpha=ribbon.alpha)
        + facet_wrap(~ gene)
        + scale_x_continuous(name="Pseudotime",
                             breaks=unique(dls[[1]]$cell.meta$obstime))
        + scale_y_continuous(name="Expression")
    )
}

#' Plot best sample predicted expression
#'
#' @param dl de.lorean object
#' @param genes Genes to plot (defaults to genes.high.psi)
#' @param add.data Add actual expression data to plot
#' @param sample.iter Which sample to plot
#'
#' @export
#'
plot.profiles <- function(dl,
                          genes=dl$genes.high.psi,
                          profile.color='black',
                          add.data=T,
                          sample.iter=dl$best.sample,
                          ...) {
    varargs <- list(...)
    with(dl, {
        if (opts$periodic) {
            modulo.period <- function(t) ( t - floor(t / opts$period)
                                                * opts$period )
        } else {
            modulo.period <- function(t) { t }
        }
        gp <- (ggplot(predictions
                      %>% filter(sample.iter == iter)
                      %>% left_join(gene.map)
                      %>% filter(gene %in% genes),
                      environment=environment()))
        profile.data <- (
            predictions
            %>% filter(sample.iter == iter)
            %>% left_join(dl$gene.map)
            %>% filter(gene %in% genes)
        )
        # stopifnot(! any(is.na(profile.data %>% select(-cbRank, -cbPeaktime))))
        gp <- (
            plot.add.mean.and.variance(
                gp,
                .data=mutate.profile.data(profile.data),
                color=profile.color)
            + facet_wrap(~ gene)
            + scale_x_continuous(name="Pseudotime",
                                 breaks=unique(cell.meta$obstime))
            + scale_y_continuous(name="Expression")
        )
        if (add.data) {
            expr.data <- (
                gene.map
                %>% filter(gene %in% genes)
                %>% left_join(melt(unname(expr),
                                   varnames=c("g", "c"),
                                   value.name="expr"))
                %>% left_join(samples.l$tau
                              %>% filter(sample.iter == iter)
                              %>% mutate(tau=modulo.period(tau))))
            if (! is.null(varargs$cell.size.adj) && varargs$cell.size.adj) {
                expr.data <- (
                    expr.data
                    %>% left_join(samples.l$S)
                    %>% mutate(expr=expr - S))
            }
            gp <- plot.add.expr(gp, .data=expr.data)
        }
        gp
    })
}

#' Mutate the profile data into shape compatible with GP plot function
#'
mutate.profile.data <- function(.data) {
    (
        .data
        %>% mutate(x=tau, mean=predictedmean+phi, var=predictedvar)
        %>% dplyr::select(-tau, -predictedmean, -phi, -predictedvar)
    )
}


#` Adjust the predicted mean with the predictions from the model.
#'
adjust.predictions <- function(.data, adjust.model) {
    # print(names(.data))
    adjustments <- (
        .data
        %>% group_by(t)
        %>% summarise(tau=tau[1]))
    # print(tail(adjustments))
    adjustments$adjustment <- predict(adjust.model,
                                        newdata=adjustments)
    # .T <- nrow(adjustments)
    # print(adjustments$adjustment[1:(.T-1)]-adjustments$adjustment[2:.T])
    (
        .data
        %>% left_join(dplyr::select(adjustments, -tau))
        %>% mutate(predictedmean=predictedmean+adjustment))
}

#' Add expression data to a plot
#'
#' @param gp Plot object
#' @param dl de.lorean object
#' @param genes Genes to plot
#' @param sample.iter Which sample to plot
#' @param cell.size.adj Adjust expression by posterior estimate of cell size
#'
#' @export
#'
plot.add.expr <- function(gp, .data=NULL)
{
    (gp + geom_point(data=.data,
                     aes(x=tau,
                         y=expr,
                         color=capture),
                     size=4,
                     alpha=.7))
}


#' Plot the expression data by the capture points
#'
#' @param dl de.lorean object
#' @param genes Genes to plot. If NULL plots some random varying genes.
#'
#' @export
#'
plot.expr.data <- function(dl, genes=NULL, num.genes=12) {
    with(dl, {
         if (is.null(genes)) {
             num.to.sample <- min(nrow(expr), num.genes * 10)
             sample.genes <- sample(rownames(expr), num.to.sample)
             expr.l <- (
                 expr[sample.genes,]
                 %>% melt(varnames=c("gene", "cell"), value.name="x"))
             variation <- (
                 expr.l
                 %>% group_by(gene)
                 %>% summarise(var=var(x))
                 %>% arrange(-var))
             if (nrow(variation > num.genes)) {
                 variation <- variation %>% head(num.genes)
             }
             genes <- variation$gene
         } else {
            expr.l <- (
                expr[genes,]
                %>% melt(varnames=c("gene", "cell"), value.name="x"))
         }
         stopifnot(all(genes %in% rownames(expr)))
         expr.l <- (
            expr.l
            %>% mutate(cell=factor(cell, levels=levels(cell.meta$cell)))
            %>% left_join(cell.meta))
        (ggplot(expr.l %>% filter(gene %in% genes),
                aes(x=capture, y=x))
            #+ geom_boxplot()
            + geom_violin()
            + stat_summary(fun.y=mean,
                           colour="red",
                           aes(group=gene),
                           geom="line")
            + facet_wrap(~ gene))
    })
}

