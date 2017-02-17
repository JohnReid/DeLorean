#
# Code to initialise latent variables (particularly tau) with sensible values.
#


#' Returns a function that constructs parameter settings with good tau if they are available.
#'
#' @param dl de.lorean object
#'
make.init.fn <- function(dl) {
  function(chain_id) {
    #
    # Create random parameters
    pars <- init.from.tau.prior(dl)
    #
    # Use initialisations for tau if we have them
    if (chain_id <= length(dl$tau.inits)) {
      #
      # Replace tau with good tau
      pars$tau <- dl$tau.inits[[chain_id]]$tau
      #
      # Also replace the tau offsets
      pars$tauoffset <- dl$tau.inits[[chain_id]]$tau - dl$cell.map$obstime
    } else {
      message('Using tau initialisation from prior.')
    }
    pars
  }
}


# Define a function to initialise the chains
#
# @param dl de.lorean object
#
make.init.from.prior.fn <- function(dl) {
  function() {
    with(dl$stan.data, {
      # message("Creating initialisation")
      init <- list(
        S = dl$cell.map$S.hat,
        tauoffset = rnorm(C, mean = 0, sd = sigma_tau),
        z = rnorm(C, mean = 0, sd = 1),
        psi = rlnorm(G, meanlog = mu_psi, sdlog = sigma_psi),
        omega = rlnorm(G, meanlog = mu_omega, sdlog = sigma_omega)
      )
      init$tau <- init$tauoffset + time
      init
    })
  }
}


# Choose an initialisation by sampling tau from the prior. Use estimated
# values for all other parameters.
#
# @param dl de.lorean object
#
init.from.tau.prior <- function(dl) {
  with(dl$stan.data, {
    init <- list(
      S = dl$cell.map$S.hat,
      z = dl$cell.map$z.hat,
      psi = dl$gene.map$psi.hat[1:G],
      omega = dl$gene.map$omega.hat[1:G],
      tauoffset = rnorm(C, 0, sd = sigma_tau)
    )
    init$tau <- init$tauoffset + time
    init
  })
}


# Spread tau values from min(time) - sigma_tau to max(time) + sigma_tau
#
even.tau.spread <- function(dl) {
    with(dl$stan.data,
         seq(min(time) - sigma_tau,
             max(time) + sigma_tau,
             length = C))
}

# From http://stackoverflow.com/questions/11094822/numbers-in-geometric-progression
geom.series <- function(base, max) {
  base^(0:floor(log(max, base)))
}

#' Use seriation package to find good orderings
#'
#' @param dl de.lorean object
#' @param .methods The seriation methods to apply
#' @param scaled Whether to use the scaled and/or unscaled expression data
#' @param dim.red Dimension reduction methods to apply
#' @param dims Number of dimensions to reduce to
#' @param num.cores Number of cores to use in parallel
#' @param num.tau.to.keep How many initialisations to keep
#'
#' @export
#'
seriation.find.orderings <- function(
  dl,
  # .methods = c("ARSA", "TSP", "R2E", "HC", "GW", "OLO"),
  .methods = c("TSP", "R2E", "HC", "GW", "OLO"),
  scaled = c('scaled', 'unscaled'),
  dim.red = c('none', 'pca', 'kfa', 'ica', 'mds'),
  # dim.red = c('mds'),
  dims = geom.series(base=2, max=min(8, nrow(dl$expr)-1)),
  num.cores = default.num.cores(),
  num.tau.to.keep = default.num.cores())
{
  #
  # Calculate all combinations of parameters
  combinations <- expand.grid(method = .methods,
                              scaled = scaled,
                              dim.red = dim.red,
                              dims = dims,
                              stringsAsFactors = FALSE) %>%
    # Only keep one combination with no dimensionality reduction
    dplyr::filter(dim.red != 'none' | 1 == dims) %>%
    # Cannot do KFA with 1 dimension
    dplyr::filter(dim.red != 'kfa' | 1 != dims)
  get.expr <- function(scaled) switch(scaled,
    scaled = scale(t(dl$expr)),
    unscaled = t(dl$expr),
    stop('scaled must be "scaled" or "unscaled"')
  )
  get.expr.mem <- memoise::memoise(get.expr)
  get.red <- function(scaled, dim.red, dims=5) {
    expr <- get.expr(scaled)
    switch(dim.red,
      none = expr,
      pca = prcomp(expr)$x[,1:dims],
      ica = fastICA::fastICA(expr, n.comp = dims)$S,
      kfa = t(kernlab::kfa(t(expr), features = dims)@xmatrix),
      mds = cmdscale(dist(expr), k = dims),
      stop('dim.red must be "none", "ica", "kfa", "mds" or "pca"'))
  }
  get.red.mem <- memoise::memoise(get.red)
  get.dist <- function(scaled, dim.red, dims = 5) {
    dist(get.red.mem(scaled, dim.red, dims))
  }
  get.dist.mem <- memoise::memoise(get.dist)
  result <- parallel::mclapply(
    1:nrow(combinations),
    function(i) with(combinations[i,], within(list(), {
      method.name <- stringr::str_c(method,
                                    ':', scaled,
                                    ':', dim.red,
                                    ':', dims)
      elapsed <- system.time(
        per <- seriation::seriate(get.dist.mem(scaled, dim.red, dims),
                                   method = method))
      ser.order <- seriation::get_order(per)
      ll <- dl$ordering.ll(ser.order)
      message(method.name, '; time=', elapsed[3], 's; LL=', ll)
    })),
    mc.cores = num.cores)
  result
}


#' Run a find good ordering method and append results to existing orderings
#'
#' @param dl de.lorean object
#' @param method Function that runs the method
#' @param ... Any other arguments for the method
#'
#' @export
#'
find.good.ordering <- function(dl, method, ...)
{
  dl <- create.ordering.ll.fn(dl)
  dl$order.inits <- c(
    dl$order.inits,
    lapply(
      method(dl, ...),
      function(i) {
        i$ser.order <- rev.order.if.better(dl, i$ser.order)
        i
      }))
  dl
}


#' Plot likelihoods of orderings against elapsed times taken
#' to generate them
#'
#' @param dl The DeLorean object
#'
#' @export
#'
orderings.plot <- function(dl) with(dl, {
  results.df <- data.frame(
    method = sapply(order.inits, function(r) r$method.name),
    elapsed = sapply(order.inits, function(r) r$elapsed[3]),
    ll = sapply(order.inits, function(r) r$ll))
  ggplot2::ggplot(results.df, aes(x = elapsed, y = ll, label = method)) +
    ggplot2::geom_text()
})

#' Use Magda's code to find good orderings
#'
#' @param dl de.lorean object
#' @param number_paths Number of paths for each starting point
#'
#' @export
#'
magda.find.orderings <- function(
    dl,
    number_paths = 5)
{
  #
  # Determine which cell indexes have either the highest or lowest
  # observation times, we will use these as starting points
  max.obs <- max(dl$cell.map$obstime)
  min.obs <- min(dl$cell.map$obstime)
  starting_points <-
    (1:dl$.C)[which(dl$cell.map$obstime %in% c(max.obs, min.obs))]
  #
  # Call Magda's function to generate paths from these starting points
  elapsed <- system.time(
    magda.paths <- CombfuncPaths(
      dl$expr,
      starting_points = starting_points,
      number_paths = number_paths))
  #
  # Convert Magda's paths into our format
  lapply(
    1:ncol(magda.paths),
    function(c) within(list(), {
      method.name <- stringr::str_c('Magda:', c)
      elapsed <- elapsed / ncol(magda.paths)
      ser.order <- magda.paths[,c]
      ll <- dl$ordering.ll(ser.order)
    }))
}

# Reverse ordering if it is better correlated with observed times
#
rev.order.if.better <- function(dl, ser.order) {
  rev.order <- rev(ser.order)
  if (cor(ser.order, dl$cell.map$obstime) < cor(rev.order, dl$cell.map$obstime)) {
    ser.order <- rev.order
  }
  ser.order
}

# Deduplicate orderings
#
deduplicate.orderings <- function(dl) {
  orderings <- sapply(dl$order.inits, function(i) i$ser.order)
  dupd <- duplicated(t(orderings))
  dl$order.inits <- dl$order.inits[!dupd]
  dl
}

#' Convert best orderings into initialisations
#'
#' @param dl The DeLorean object
#' @param num.to.keep The number to keep (defaults to default.num.cores())
#'
#' @export
#'
pseudotimes.from.orderings <- function(
  dl,
  num.to.keep = default.num.cores())
within(deduplicate.orderings(dl), {
  order.orderings <- order(sapply(order.inits, function(i) i$ll),
                           decreasing = TRUE)
  #
  # Make sure we don't try to keep too many (due to duplicates, etc...)
  actually.keep <- min(length(order.orderings), num.to.keep)
  if (actually.keep < num.to.keep) {
    warning("Don't have enough ", num.to.keep, " pseudotimes to keep,",
            " only have ", actually.keep)
  }
  best.orderings <- order.inits[order.orderings[1:actually.keep]]
  tau.inits <- lapply(
    best.orderings,
    function(O) {
      message('Using ordering ', O$method.name, '; LL = ', O$ll)
      # Create an initialisation using the ordering
      init <- init.from.tau.prior(dl)
      init$tau <- even.tau.spread(dl)[O$ser.order]
      init
    })
})


#' Find best order of the samples assuming some smooth GP prior on the
#' expression profiles over this ordering.
#'
#' @param dl de.lorean object
#' @param psi Temporal variation
#' @param omega Noise
#' @param num.cores Number of cores to run on. Defaults to default.num.cores()
#' @param num.tau.to.try How many initialisations to try
#' @param num.tau.to.keep How many initialisations to keep
#' @param method Method to use "maximise" or "metropolis"
#' @param ... Extra arguments to method
#'
#' @export
#'
find.smooth.tau <- function(
    dl,
    psi = exp(dl$hyper$mu_psi),
    omega = exp(dl$hyper$mu_omega),
    num.cores = default.num.cores(),
    num.tau.to.try = num.cores,
    num.tau.to.keep = num.cores,
    method = "metropolis",
    ...
) {
  dl <- create.ordering.ll.fn(dl)
  log.likelihood <- dl$ordering.ll
  dl$tau.inits <- with(dl, {
    # Maximise the sum of the log marginal likelihoods
    ordering.search <- function(seed) {
      set.seed(seed)
      # Choose a starting point by random projection
      expr.centre <- t(scale(t(expr), center = T, scale = F))
      init.ordering <- order(rnorm(nrow(expr.centre)) %*% expr.centre)
      # init.ordering <- sample(stan.data$C)
      metropolis.fn <- function(ordering, log.likelihood, ...) {
        mh.run <- ordering.metropolis.hastings(
          ordering,
          log.likelihood,
          proposal.fn = ordering.random.block.move,
          ...)
        best.sample <- which.max(mh.run$log.likelihoods)
        #ordering.maximise(mh.run$chain[best.sample,], log.likelihood)
        mh.run$chain[best.sample,]
      }
      method.fn <- switch(method,
                          "maximise" = ordering.maximise,
                          "metropolis" = metropolis.fn,
                          NA)
      ordering <- method.fn(init.ordering, log.likelihood, ...)
      stopifnot(! is.null(ordering))
      # Reverse the ordering if it makes it correlate better with
      # the capture times
      capture.order <- order(stan.data$time)
      if (cor(capture.order, ordering) <
          cor(capture.order, rev(ordering)))
      {
        ordering <- rev(ordering)
      }
      ordering
    }
    # Choose seeds
    seeds <- sample.int(.Machine$integer.max, num.tau.to.try)
    # Run in parallel or not?
    if (num.cores > 1) {
      orderings <- parallel::mclapply(seeds,
                                      mc.cores=num.cores,
                                      ordering.search)
    } else {
      orderings <- lapply(seeds, ordering.search)
    }
    # Order the taus by the best orderings
    lls <- sapply(orderings, function(o) -log.likelihood(o))
    best.order <- order(lls)
    # Make the complete chain initialisation with the tau.
    lapply(
      orderings[best.order[1:num.tau.to.keep]],
      function(ordering) {
        init <- init.from.tau.prior(dl)
        init$tau <- even.tau.spread(dl)[ordering.invert(ordering)]
        init
      })
  })
  dl
}


#' Test ordering Metropolis-Hastings sampler.
#'
#' @param dl de.lorean object
#' @param psi Temporal variation
#' @param omega Noise
#' @param num.cores Number of cores to run on.
#'          Defaults to getOption("DL.num.cores", max(parallel::detectCores()-1, 1))
#' @param iterations Number of iterations
#' @param thin Thin the samples
#'
test.mh <- function(
    dl,
    psi = mean(dl$gene.map$psi.hat),
    omega = mean(dl$gene.map$omega.hat),
    num.cores = getOption("DL.num.cores", max(parallel::detectCores() - 1, 1)),
    iterations = 1000,
    thin = 15
) {
    dl <- create.ordering.ll.fn(dl)
    log.likelihood <- dl$ordering.ll
    dl$tau.inits <- with(dl, {
        # Maximise the sum of the log marginal likelihoods
        ordering.search <- function(seed) {
            set.seed(seed)
            # Choose a random starting point
            init.ordering <- sample(stan.data$C)
            mh.run <- ordering.metropolis.hastings(
                init.ordering,
                log.likelihood,
                proposal.fn = ordering.random.block.move,
                iterations = iterations,
                thin = thin)
        }
        # Choose seeds
        seeds <- sample.int(.Machine$integer.max, num.cores)
        # Run in parallel or not?
        if (num.cores > 1) {
            orderings <- parallel::mclapply(seeds,
                                  mc.cores=num.cores,
                                  ordering.search)
        } else {
            orderings <- lapply(seeds, ordering.search)
        }
        orderings
    })
}


#' Find best tau to initialise chains with by sampling tau from the prior
#' and using empirical Bayes parameter estimates for the other parameters.
#'
#' @param dl de.lorean object
#' @param num.tau.candidates How many candidates to examine. Defaults to 6000.
#' @param num.tau.to.keep How many candidates to keep. Defaults to num.cores.
#' @param num.cores Number of cores to run on. Defaults to default.num.cores()
#'
#' @export
#'
find.best.tau <- function(
  dl,
  num.tau.candidates = 6000,
  num.tau.to.keep = num.cores,
  num.cores = default.num.cores())
within(dl, {
  #
  # Define a function that calculates log probability for
  # random seeded tau
  try.tau.init <- function(i) {
    set.seed(i)
    pars <- init.from.tau.prior(dl)
    lp <- rstan::log_prob(fit, rstan::unconstrain_pars(fit, pars),
                          adjust_transform = FALSE)
    list(lp = lp, tau = pars$tau)
  }
  #
  # Choose tau several times and calculate log probability
  if (num.cores > 1) {
    tau.inits <- parallel::mclapply(1:num.tau.candidates,
                                    mc.cores = num.cores,
                                    try.tau.init)
  } else {
    tau.inits <- lapply(1:num.tau.candidates, try.tau.init)
  }
  #
  # Which tau gave highest log probability?
  tau.inits.order <- order(sapply(tau.inits, function(init) -init$lp))
  #
  # Just keep so many best tau inits
  tau.inits <- tau.inits[tau.inits.order[1:num.tau.to.keep]]
  rm(tau.inits.order, try.tau.init)
})
