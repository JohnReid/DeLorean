# Is a package installed?
#
is.installed <- function(mypkg)
  is.element(mypkg, utils::installed.packages()[,1])

# Configure ggplot2 colours
if (is.installed('ggthemes')) {
  scale_colour_discrete <- function(...) ggthemes::scale_colour_few(...)
  scale_fill_discrete <- function(...) ggthemes::scale_fill_few(...)
}

#' Various DeLorean object plots
#'
#' @param x de.lorean object
#' @param type Type of plot:
#'   \itemize{
#'     \item 'expr.data': The expression data plotted by capture time.
#'       See \code{\link{expr.data.plot}}.
#'     \item 'Rhat': \eqn{hat{R}} convergence statistics
#'       See \code{\link{Rhat.plot}}.
#'     \item 'pseudotime': Pseudotimes in best posterior sample
#'       See \code{\link{pseudotime.plot}}.
#'     \item 'profiles': Gene expression profiles for best posterior sample
#'       See \code{\link{profiles.plot}}.
#'     \item 'tau.offsets': Offsets of pseudotimes to assess the prior
#'       See \code{\link{tau.offsets.plot}}.
#'     \item 'marg.like': Plot the posterior of the marginal likelihoods
#'       for individual genes.
#'       See \code{\link{marg.like.plot}}.
#'     \item 'roughnesses': Roughnesses of the pseudotime posterior
#'       See \code{\link{roughnesses.plot}}.
#'     \item 'init.vs.pseudotimes': Plot the initialisations used against the pseudotimes estimated
#'       See \code{\link{init.orderings.vs.pseudotimes.plot}}.
#'   }
#' @param ... Extra arguments to plot function
#'
#' @method plot de.lorean
#' @export
#'
plot.de.lorean <- function(x, type = "profiles", ...) {
    result <- switch(type,
        profiles = profiles.plot(x, ...),
        pseudotime = pseudotime.plot(x, ...),
        Rhat = Rhat.plot(x, ...),
        expr.data = expr.data.plot(x, ...),
        roughnesses = roughnesses.plot(x, ...),
        marg.like = marg.like.plot(x, ...),
        orderings = orderings.plot(x, ...),
        tau.offsets = tau.offsets.plot(x, ...),
        gene.params = gene.params.plot(x, ...),
        init.vs.pseudotimes = init.orderings.vs.pseudotimes.plot(x, ...)
    )
    if (is.null(result)) {
        stop('Unknown plot type')
    }
    result
}

#' Calculate a suitable value for a rug plot given the
#' number of points
#'
#' @param n Number of points.
#' @param scale Scale the value.
#'
alpha.for.rug <- function(n, scale=100) {
    1 / (max(1, n / scale))
}

#' Plot posterior for marginal log likelihoods of individual gene's
#' expression profiles
#'
#' @param dl de.lorean object
#'
#' @export
#'
marg.like.plot <- function(dl) with(dl,
  ggplot(samples.l$logmarglike %>% left_join(gene.map),
         aes(x = reorder(gene, logmarglike, FUN = median),
             y = logmarglike,
             colour = is.held.out),
         environment=environment()) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)))

#' Plot posterior of psi, the temporal variation parameter for all genes.
#'
#' @param dl de.lorean object
#'
#' @export
#'
post.psi.plot <- function(dl) with(dl,
  ggplot(dl$samples.l$psi, aes(x = log(psi))) +
    geom_density() +
    geom_rug() +
    stat_function(fun = function(x) dnorm(x,
                                          mean = hyper$mu_psi,
                                          sd = hyper$sigma_psi),
                  colour = "blue", alpha = .7, linetype = "dashed"))

#' Boxplot posterior of psi, the temporal variation parameter for each gene.
#'
#' @param dl de.lorean object
#'
#' @export
#'
post.gene.psi.plot <- function(dl) with(dl,
  ggplot(samples.l$psi %>% left_join(gene.map),
         aes(x = reorder(gene, psi, FUN = median),
             y = log(psi))) +
    geom_boxplot() +
    geom_hline(yintercept = stan.data$mu_psi, linetype = "dashed"))

#' Plot posterior of omega, the noise parameter for all genes.
#'
#' @param dl de.lorean object
#'
#' @export
#'
post.omega.plot <- function(dl) with(dl,
  ggplot(samples.l$omega, aes(x = log(omega))) +
    geom_density() +
    geom_rug() +
    stat_function(fun = function(x) dnorm(x,
                                          mean = hyper$mu_omega,
                                          sd = hyper$sigma_omega),
                  colour = "blue", alpha = .7, linetype = "dashed"))

#' Boxplot posterior of omega, the noise parameter for each gene.
#'
#' @param dl de.lorean object
#'
#' @export
#'
post.gene.omega.plot <- function(dl) with(dl,
  ggplot(samples.l$omega %>% left_join(gene.map),
         aes(x = reorder(gene, omega, FUN = median),
             y = log(omega))) +
    geom_boxplot() +
    geom_hline(yintercept = stan.data$mu_omega, linetype = "dashed"))

#' Plot pseudotime (tau) against observed capture time.
#'
#' @param dl de.lorean object
#' @param sample.iter Which sample to take pseudotimes from
#'
#' @export
#'
pseudotime.plot <- function(dl, sample.iter=dl$best.sample) {
    with(dl, {
        gp <- (ggplot(samples.l$tau %>% filter(iter == sample.iter),
                      aes(x=tau, y=obstime, color=capture),
                      environment=environment())
            + geom_point()
            + geom_vline(data=cell.meta %>% group_by(capture),
                         aes(xintercept=obstime, color=capture),
                         linetype=2,
                         alpha=.8)
            + scale_x_continuous(name="Pseudotime (tau)")
            + scale_y_continuous(name="Observed (capture) time")
        )
    })
}


#' Plot the tau offsets, that is how much the pseudotimes (tau) differ
#' from their prior means over the full posterior.
#'
#' @param dl de.lorean object
#' @param rug.alpha Alpha parameter for rug geom
#'
#' @export
#'
tau.offsets.plot <- function(dl, rug.alpha=.3) with(dl, {
  prior.scale <- nrow(samples.l$tau) / length(levels(samples.l$tau$capture))
  ggplot(samples.l$tau, aes(x=tau.offset, color=capture, fill=capture)) +
    geom_histogram(position='dodge') +
    geom_rug(alpha=rug.alpha) +
    stat_function(fun=function(x) prior.scale * dnorm(x, sd=hyper$sigma_tau),
                  linetype='dashed')
})


#' Plot the posterior of the gene parameters.
#'
#' @param dl de.lorean object
#' @param alpha Transparency for points and error bars
#'
#' @export
#'
gene.params.plot <- function(dl, alpha = .4)
  ggplot(dl$var.post.stats,
         aes(x = psi.mean, xmin = psi.mean - psi.sd, xmax = psi.mean + psi.sd,
             y = omega.mean, ymin = omega.mean - omega.sd, ymax = omega.mean + omega.sd,
             label = gene)) +
  geom_point(alpha = alpha) +
  geom_errorbar(alpha = alpha) +
  geom_errorbarh(alpha = alpha) +
  geom_label(alpha = .8)


#' Plot the Rhat convergence statistics. \code{\link{examine.convergence}}
#' must be called before this plot can be made.
#'
#' @param dl de.lorean object
#'
#' @export
#'
Rhat.plot <- function(dl) {
    with(dl, {
        rhat.df <- data.frame(
            rhat=rhat.sorted,
            param=names(rhat.sorted),
            parameter=stringr::str_match(names(rhat.sorted), "^[[:alpha:]]+"))
        gp <- (ggplot(rhat.df,
                      aes(y=rhat, x=parameter),
                      environment=environment())
            + geom_boxplot()
        )
    })
}

#' Plot a comparison of the profiles from several de.lorean objects
#'
#' @param ... Named de.lorean objects
#' @param genes Genes to plot (defaults to genes.high.psi of first de.lorean
#'   object)
#'
#' @export
#'
cmp.profiles.plot <- function(..., genes = NULL) {
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

#' Plot best sample predicted expression.
#'
#' @param dl de.lorean object
#' @param genes Genes to plot (defaults to genes.high.psi)
#' @param profile.color Colour for the profile
#' @param add.data Add actual expression data to plot
#' @param sample.iter Which sample to plot
#' @param ignore.cell.sizes Ignore cell sizes if the model has estimated them
#' @param ... Extra arguments
#'
#' @export
#'
profiles.plot <- function(
  dl,
  genes=dl$genes.high.psi,
  profile.color='black',
  add.data=T,
  sample.iter=dl$best.sample,
  ignore.cell.sizes=FALSE,
  ...)
{
  varargs <- list(...)
  with(dl, {
    if (opts$periodic) {
      modulo.period <- function(t) ( t - floor(t / opts$period)
                                          * opts$period )
    } else {
      modulo.period <- function(t) { t }
    }
    gp <-
      ggplot(
        predictions %>%
          filter(sample.iter == iter) %>%
          left_join(gene.map) %>%
          filter(gene %in% genes),
        environment=environment())
    profile.data <-
      predictions %>%
      filter(sample.iter == iter) %>%
      left_join(dl$samples.l$omega) %>%
      left_join(dl$gene.map) %>%
      left_join(dl$gene.expr) %>%
      filter(gene %in% genes)
    # stopifnot(! any(is.na(profile.data %>% select(-cbRank, -cbPeaktime))))
    gp <-
      plot.add.mean.and.variance(
        gp,
        .data=mutate.profile.data(profile.data),
        color=profile.color,
        add.noise=TRUE) +
      facet_wrap(~ gene) +
      scale_x_continuous(
        name="Pseudotime",
        breaks=unique(cell.meta$obstime)) +
       scale_y_continuous(name="Expression")
    if (add.data) {
        expr.data <-
          gene.map %>%
          filter(gene %in% genes) %>%
          left_join(
            melt(
              unname(expr),
              varnames=c("g", "c"),
              value.name="expr")) %>%
          left_join(
            samples.l$tau %>%
              filter(sample.iter == iter) %>%
              mutate(tau=modulo.period(tau)))
        #
        # Adjust for cell sizes if they are there and we have not
        # been asked to ignore them
        if ('S' %in% names(samples.l) && ! ignore.cell.sizes) {
            expr.data <-
              expr.data %>%
              left_join(samples.l$S) %>%
              mutate(expr=expr - S)
        }
        gp <- plot.add.expr(gp, .data=expr.data)
    }
    gp
  })
}


#' Mutate the profile data into shape compatible with GP plot function
#'
#' @param .data The data
#'
mutate.profile.data <- function(.data) .data %>%
  mutate(x=tau, mean=predictedmean+phi.hat, var=predictedvar, noise2=omega) %>%
  dplyr::select(-tau, -predictedmean, -phi.hat, -predictedvar)


#' Add expression data to a plot
#'
#' @param gp Plot object
#' @param .data Expression data to add
#'
plot.add.expr <- function(gp, .data=NULL) {
    gp + geom_point(data=.data, aes(x=tau, y=expr, color=capture),
                    size=4, alpha=.7)
}


#' Plot the expression data by the capture points
#'
#' @param dl de.lorean object
#' @param genes Genes to plot. If NULL plots some random varying genes
#' @param num.genes Number of genes to plot
#'
#' @export
#'
expr.data.plot <- function(dl, genes=NULL, num.genes=12) {
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
                 %>% dplyr::summarise(var=var(x))
                 %>% arrange(-var))
             if (nrow(variation) > num.genes) {
                 variation <- variation %>% head(num.genes)
             }
             genes <- variation$gene
         } else {
            expr.l <- (
                expr[genes,]
                %>% melt(varnames=c("gene", "cell"), value.name="x"))
         }
         stopifnot(all(genes %in% rownames(expr)))
         expr.l <- expr.l %>%
             mutate(cell=factor(cell, levels=levels(cell.meta$cell))) %>%
             left_join(cell.meta) %>%
             dplyr::filter(gene %in% genes)
         ggplot(expr.l, aes(x=capture, y=x)) +
             # geom_boxplot() +
             geom_violin() +
             stat_summary(fun.y=mean,
                          colour="red",
                          aes(group=gene),
                          geom="line") +
             facet_wrap(~ gene)
    })
}

# Scale and shift x to match the range of the reference.
#
match.range <- function(x, reference) {
  min.x <- min(x)
  max.x <- max(x)
  width.x <- max.x - min.x
  min.r <- min(reference)
  max.r <- max(reference)
  width.r <- max.r - min.r
  min.r + (x - min.x) / width.x * width.r
}

#' Plot the orderings for initialisation against the estimated pseudotime.
#'
#' @param dl The DeLorean object
#' @param sample.iter Which sample to take pseudotimes from
#'
#' @export
#'
init.orderings.vs.pseudotimes.plot <- function(dl, sample.iter = dl$best.sample) {
  pseudotimes <- tau.for.sample(dl, sample.iter = sample.iter)
  ordering.pseudotime <- order(pseudotimes)
  ith.best.ordering <- function(i) {
    init.order <- dl$best.orderings[[i]]
    data.frame(
      idx = 1:length(ordering.pseudotime),
      ordering.rank = i,
      init.method = init.order$method.name,
      LL = -init.order$ll,
      initialisation = match(init.order$ser.order, ordering.pseudotime))
  }
  orderings <-
      do.call(rbind, lapply(1:length(dl$best.orderings), ith.best.ordering)) %>%
      mutate(initialisation = match.range(initialisation, pseudotimes)) %>%
      left_join(data.frame(
          idx = 1:length(ordering.pseudotime),
          capture = dl$cell.map$capture[ordering.pseudotime],
          pseudotime = pseudotimes,
          pseudotime.order = ordering.pseudotime))
  ggplot2::ggplot(orderings, aes(y=pseudotime, yend=initialisation, color=capture)) +
      ggplot2::geom_segment(x = 0, xend = 1) +
      facet_wrap(~ LL + init.method)
}

#' Plot two sets of pseudotimes against each other.
#'
#' @param dl The DeLorean object
#' @param fits Fit indexes
#'
#' @export
#'
pseudotimes.pair.plot <- function(dl, fits=NULL) {
  stopifnot(all(dim(dl$best.orderings) == dim(dl$vb.tau)))
  stopifnot(! is.null(dl$vb.tau))  # Must have estimated tau.
  #
  # Create a data frame with the best orderings
  best.o <- sapply(dl$best.orderings, function(O) O$ser.order)
  best.o.m <- reshape2::melt(best.o, varnames=c('c', 'fit'),
                             value.name='ordering.idx') %>%
              dplyr::left_join(data.frame(
                ordering.idx=1:dl$stan.data$C,
                ordering=even.tau.spread(dl)))
  #
  # Create a data frame with the best pseudotimes
  best.tau.m <- reshape2::melt(dl$vb.tau, varnames=c('c', 'fit'),
                               value.name='pseudotime')
  #
  # If not plotting all fits then filter data frames
  if (! is.null(fits)) {
    best.o.m <- filter(best.o.m, fit %in% fits)
    best.tau.m <- filter(best.tau.m, fit %in% fits)
  }
  #
  # Combine the data frames and join other data
  df. <- best.o.m %>%
    left_join(best.tau.m) %>%
    dplyr::mutate(
      ordering.label=factor('ordering', c('ordering', 'pseudotime')),
      pseudotime.label=factor('pseudotime', c('ordering', 'pseudotime'))) %>%
    dplyr::left_join(dl$cell.map)
  ggplot2::ggplot(
      df.,
      aes(x=ordering.label, xend=pseudotime.label,
          y=ordering, yend=pseudotime,
          color=capture)) +
    ggplot2::geom_segment(alpha=.3) + facet_wrap(~ fit)
}


#' Plot the posterior of the branching model
#'
#' @param dl The DeLorean object
#' @param post The posterior (output of branching.genes.post)
#' @param facets Facets to wrap plots (passed to facet_wrap())
#'
#' @export
#'
branching.post.plot <- function(dl, post, facets = ~ gene) with(dl, {
  ggplot(post, aes(x = tau, y = z, fill = post.mean)) +
    geom_tile(width = tau.step, height = z.step) +
    scale_fill_gradient2(low = scales::muted("blue"), mid = "white", high = scales::muted("red")) +
    facet_wrap(facets)
})


#' Plot the posterior of the temporal variation and show the prior
#'
#' @param dl The DeLorean object
#'
#' @export
#'
psi.post.plot <- function(dl) with(dl,
  ggplot(samples.l$psi, aes(x = log(psi))) +
    geom_density() +
    geom_rug() +
    stat_function(
      fun = function(x) dnorm(x, mean = hyper$mu_psi, sd = hyper$sigma_psi),
      colour = "blue", alpha = .7, linetype = "dashed"))


#' Plot the posterior of the noise levels and show the prior
#'
#' @param dl The DeLorean object
#'
#' @export
#'
omega.post.plot <- function(dl) with(dl,
  ggplot(samples.l$omega, aes(x = log(omega))) +
    geom_density() +
    geom_rug() +
    stat_function(
      fun = function(x) dnorm(x, mean = hyper$mu_omega, sd = hyper$sigma_omega),
      colour = "blue", alpha = .7, linetype = "dashed"))
