# Update a factor based on the levels of another factor.
#
# @param .factor Factor
# @param reference.factor Factor whose levels to use.
# @param levels New levels
#
update.levels <- function(.factor,
                          reference.factor=NULL,
                          .levels=levels(reference.factor),
                          ordered=is.ordered(reference.factor))
{
    factor(.factor, levels=.levels, ordered=ordered)
}


# Retrieve the estimated tau for the given sample.
#
# @param dl de.lorean object
# @param sample.iter Which sample to use, defaults to best sample
#
tau.for.sample <- function(dl, sample.iter = dl$best.sample) {
    (
        dl$samples.l$tau
        %>% filter(sample.iter == iter)  # Filter correct iteration
        %>% arrange(c)  # Sort by cell
    )$tau
}


# The levels of the gene factor
#
# @param dl The de.lorean object.
#
gene.levels <- function(dl) levels=levels(dl$gene.meta$gene)


# The levels of the cell factor
#
# @param dl The de.lorean object.
#
cell.levels <- function(dl) levels=levels(dl$cell.meta$cell)


#' Melt an expression matrix.
#'
#' @param dl The de.lorean object.
#' @param expr Matrix of expression values.
#'
#' @export melt.expr
#'
melt.expr <- function(dl, expr=dl$expr) (
    expr
    %>% melt(varnames=c("gene", "cell"), value.name="x")
    %>% mutate(gene=factor(gene, levels=gene.levels(dl)),
               cell=factor(cell, levels=cell.levels(dl)))
)


# Cast an expression matrix.
#
# @param expr.l Expression values in long format
#
cast.expr <- function(expr.l) expr.l %>% acast(gene ~ cell, value.var="x")


#' Analyse variance of expression between and within capture times.
#'
#' @param dl de.lorean object
#' @param adjust.cell.sizes Choose whether to adjust the expression values by
#'        the cell size estimates
#'
#' @export
#'
analyse.variance <- function(dl, adjust.cell.sizes) {
  #
  # First melt expression data into long format
  expr.adj <-
    melt.expr(dl) %>%
    left_join(dl$cell.sizes) %>%
    # Adjust the expression by the cell size estimates if asked to
    mutate(x.hat=x-adjust.cell.sizes*S.hat)
  stopifnot(! is.na(expr.adj))
  within(dl, {
    #
    # Resummarise the adjusted expression data by gene
    gene.expr <-
      expr.adj %>%
      group_by(gene) %>%
      dplyr::summarise(phi.hat=mean(x.hat))
    stopifnot(! is.na(gene.expr))
    #
    # Examine the expression at each gene/capture time combination
    gene.time.expr <-
      expr.adj %>%
      left_join(cell.meta) %>%
      group_by(gene, capture) %>%
      dplyr::summarise(
        num.capture=n(),
        Mgk=mean(x.hat),
        Vgk=var(x.hat)) %>%
      filter(! is.na(Vgk))  # Remove variances from capture times with one cell
    stopifnot(! is.na(gene.time.expr))
    stopifnot(nrow(gene.time.expr) > 0)  # Must have some rows left
    #
    # Expected variance of samples from zero-mean Gaussian with covariance K.obs
    hyper.capture <- unique(gene.time.expr$capture)
    hyper.obs <- unique(filter(cell.meta, capture %in% hyper.capture)$obstime)
    V.obs <-
      expected.sample.var(
        cov.matern.32(cov.calc.dists(hyper.obs), dl$opts$length.scale))
    #
    # Variance within and between times
    gene.var <-
      gene.time.expr %>%
      group_by(gene) %>%
      dplyr::summarise(omega.hat=mean(Vgk), psi.hat=var(Mgk)/V.obs) %>%
      filter(! is.na(psi.hat),
             ! is.na(omega.hat),
             omega.hat > 0,
             psi.hat > 0)
  })
}


#' Estimate hyperparameters for model using empirical Bayes.
#'
#' @param dl de.lorean object
#' @param sigma.tau Noise s.d. in temporal dimension, that is prior s.d. for tau
#' @param length.scale Length scale for stationary GP covariance function.
#'   Defaults to the range of the observed capture times.
#' @param model.name The model's name:
#'   \itemize{
#'     \item 'exact': The model without a low rank approximation
#'       that does not estimate the cell sizes.
#'     \item 'exact-sizes': The model without a low rank approximation
#'       that does estimate the cell sizes.
#'     \item 'lowrank': Low rank approximation to the 'exact' model.
#'     \item 'lowrank-sizes': Low rank approximation to the 'exact-sizes' model.
#'   }
#' @param adjust.cell.sizes Adjust by the cell sizes for better estimates of the hyperparameters
#'
#' @export
#'
estimate.hyper <- function(
    dl,
    sigma.tau = .5,
    length.scale = NULL,
    model.name = 'exact',
    adjust.cell.sizes = TRUE
) {
  #
  # Estimate the cell sizes
  # We will need these to estimate the hyper-parameters for the cell size prior if
  # the model estimates or cell sizes, or to adjust the data by if it doesn't
  dl <- estimate.cell.sizes(dl)
  dl <- within(dl, {
    #
    # Remember options that depend on the model
    opts$model.name <- model.name
    opts$model.estimates.cell.sizes <- switch(
      opts$model.name,
      "exact-sizes" = TRUE,
      "exact" = FALSE,
      "lowrank-sizes" = TRUE,
      "lowrank" = FALSE,
      "branching" = FALSE,
      "branching-lowrank" = FALSE,
      stop('Unknown model name'))
    opts$model.is.branching <- switch(
      opts$model.name,
      "exact-sizes" = FALSE,
      "exact" = FALSE,
      "lowrank-sizes" = FALSE,
      "lowrank" = FALSE,
      "branching" = TRUE,
      "branching-lowrank" = TRUE,
      stop('Unknown model name'))
    #
    # Set up temporal hyper-parameters
    opts$sigma.tau <- sigma.tau
    time.range <- range(cell.meta$obstime)
    time.width <- time.range[2] - time.range[1]
    if (is.null(length.scale)) {
        opts$length.scale <- time.width / 2
    } else {
        opts$length.scale <- length.scale
    }
  })
  dl <- analyse.variance(dl, adjust.cell.sizes = adjust.cell.sizes)
  within(dl, {
    hyper <- list(
      mu_S = mean(cell.sizes$S.hat),
      sigma_S = sd(cell.sizes$S.hat),
      mu_phi = mean(gene.expr$phi.hat),
      sigma_phi = sd(gene.expr$phi.hat),
      mu_psi = mean(log(gene.var$psi.hat), na.rm = TRUE),
      sigma_psi = sd(log(gene.var$psi.hat), na.rm = TRUE),
      mu_omega = mean(log(gene.var$omega.hat), na.rm = TRUE),
      sigma_omega = sd(log(gene.var$omega.hat), na.rm = TRUE),
      sigma_tau = opts$sigma.tau,
      l = opts$length.scale
    )
    # Check we don't have any NAs
    stopifnot(all(! sapply(hyper, is.na)))
  })
}


#' Filter genes
#'
#' @param dl de.lorean object
#' @param number Number to sample if filter function or genes not supplied.
#' @param genes The genes to keep.
#' @param .filter Function that gakes a list of genes as input and returns
#'     a vector of TRUE/FALSE
#'
#' @export
#'
filter.genes <- function(dl,
                         .filter = function(x) x %in% genes,
                         number = NULL,
                         genes = sample(rownames(dl$expr), number))
{
    within(dl, {
        expr <- expr[.filter(rownames(expr)),]
        message("Have ", nrow(expr), " genes after filtering")
    })
}


#' Filter cells
#'
#' @param dl de.lorean object
#' @param number Number to sample if filter function or cells not supplied.
#' @param cells The cells to keep.
#' @param .filter Function that gakes a list of cells as input and returns
#'     a vector of TRUE/FALSE
#'
#' @export
#'
filter.cells <- function(dl,
                         .filter = function(x) x %in% cells,
                         number = NULL,
                         cells = sample(colnames(dl$expr), number))
{
    within(dl, {
        expr <- expr[,.filter(colnames(expr))]
        message("Have ", ncol(expr), " cells after filtering")
    })
}




# Sample so many cells per capture time.
#
# @param dl de.lorean object
# @param number Number to sample from each capture time
#
sample.per.capture <- function(dl, cells.per.capture) {
    sample.at.most <- function(.df, number) {
        sample_n(.df, min(number, nrow(.df)))
    }
    sampled.cells <- (
        dl$cell.meta
        %>% filter(cell %in% colnames(dl$expr))
        %>% group_by(capture)
        # %>% do(sample_n(., min(cells.per.capture, length(cell))))
        %>% do(sample.at.most(., cells.per.capture))
    )
    filter.cells(dl, cells = sampled.cells$cell)
}


# Sample genes and cells
#
# @param dl de.lorean object
#
sample.genes.and.cells <- function(
    dl,
    max.cells = 0,
    max.genes = 0)
{
    within(dl, {
        opts$max.cells <- max.cells
        opts$max.genes <- max.genes
        if (opts$max.cells && ncol(expr) > opts$max.cells) {
            expr <- expr[,sample(ncol(expr), opts$max.cells)]
            # Remove genes that are not expressed in at least two cells
            num.cells.expr <- rowSums(! is.na(expr))
            expr <- expr[num.cells.expr > 1,]
        }
        if (opts$max.genes && nrow(expr) > opts$max.genes) {
            expr <- expr[sample(nrow(expr), opts$max.genes),]
        }
    })
}


#' Calculate inducing pseudotimes for sparse approximation
#'
#' @param dl de.lorean object
#' @param num.inducing.tau Number of inducing points
#' @param num.inducing.z Number of inducing latent dimension points
#' @param period Period of expression patterns
#' @param num.sd.border The size of the border of the inducing inputs
#'            around the capture times in units of number of standard
#'            deviations
#'
#' @export
#'
calc.inducing.pseudotimes <- function(dl, num.inducing.tau, num.inducing.z, period = 0, num.sd.border = 7) with(dl, {
  if (period > 0) {
    stopifnot(! opts$model.is.branching)  # Cannot have periodic branching processes
    seq(from = 0,
        to = period * (1 - 1/num.inducing.tau),
        length.out = num.inducing.tau)
  } else {
    if (opts$model.is.branching) {
      taupseudo <- seq(
        from = time.range[1] - num.sd.border * hyper$sigma_tau,
        to = time.range[2] + num.sd.border * hyper$sigma_tau,
        length.out = num.inducing.tau)
      zpseudo <- seq(
        from = - num.sd.border,
        to = num.sd.border,
        length.out = num.inducing.z)
      t(expand.grid(taupseudo, zpseudo))
    } else {
      seq(from = time.range[1] - num.sd.border * hyper$sigma_tau,
          to = time.range[2] + num.sd.border * hyper$sigma_tau,
          length.out = num.inducing.tau)
    }
  }
})


#' Prepare for Stan
#'
#' @param dl de.lorean object
#' @param num.test Number of test points to consider
#' @param num.inducing.tau Number of inducing pseudotime points
#' @param num.inducing.z Number of inducing latent dimension points
#' @param period Period of expression patterns
#' @param hold.out Number genes to hold out for generalisation tests
#' @param num.sd.border The size of the border of the inducing inputs
#'            around the capture times in units of number of standard
#'            deviations
#' @param l.z length scale for latent dimension
#'
#' @export prepare.for.stan
#'
prepare.for.stan <- function(
  dl,
  num.test = 101,
  num.inducing.tau = 30,
  num.inducing.z = 5,
  period = 0,
  hold.out = 0,
  num.sd.border = 7,
  l.z = 6)
within(dl, {
  opts$num.test <- num.test
  opts$period <- period
  opts$periodic <- opts$period > 0
  stopifnot(hold.out < nrow(expr))
  .G <- nrow(expr) - hold.out
  .C <- ncol(expr)
  #
  # Permute genes to make held out genes random
  expr <- expr[sample(.G+hold.out),]
  #
  # Calculate the map from gene indices to genes and their meta data
  gene.map <- (
    data.frame(g=1:(.G+hold.out),
               gene=factor(rownames(expr), levels=levels(gene.meta$gene)))
    %>% mutate(is.held.out=g>.G)
    %>% left_join(gene.expr)
    %>% left_join(gene.var)
    %>% left_join(gene.meta))
  stopifnot(! is.na(gene.map[
    c("g", "gene", "phi.hat", "psi.hat", "omega.hat")
  ]))
  #
  # Calculate the map from cell indices to genes and their meta data
  cell.map <-
    data.frame(
      c=1:.C,
      cell=factor(colnames(expr), levels=levels(cell.meta$cell))) %>%
    left_join(cell.meta) %>%
    left_join(cell.sizes)
  stopifnot(! is.na(cell.map %>% dplyr::select(cell, capture, obstime)))
  #
  # Add the z.hat estimates to the cell map if PCA data frame exists
  if (! is.null(dl$pca.df)) {
    cell.map <- dplyr::left_join(cell.map, dplyr::select(pca.df, cell, z.hat))
  }
  #
  # Calculate the time points at which to make predictions
  if (opts$periodic) {
    test.input <- seq(0, opts$period, length.out=num.test)
  } else {
    test.input <- seq(
      time.range[1] - num.sd.border * opts$sigma.tau,
      time.range[2] + num.sd.border * opts$sigma.tau,
      length.out = num.test)
  }
  #
  # Gather all the data into one list
  stan.data <- c(
    # Hyper-parameters
    hyper,
    list(
      # Dimensions
      C = .C,
      G = .G,
      H = hold.out,
      M = num.inducing.tau,
      # Data
      time = cell.map$obstime,
      expr = expr,
      phi = gene.map$phi.hat,
      # Inducing pseudotimes
      u = calc.inducing.pseudotimes(
        dl,
        num.inducing.tau = num.inducing.tau,
        num.inducing.z = num.inducing.z,
        period = period,
        num.sd.border = num.sd.border),
      # Held out parameters
      heldout_psi = filter(gene.map, g > .G)$psi.hat,
      heldout_omega = filter(gene.map, g > .G)$omega.hat,
      # Generated quantities
      numtest = opts$num.test,
      testinput = test.input,
      # Periodic?
      periodic = opts$periodic,
      period = opts$period
    )
  )
  #
  # If we have a branching model fill in all the details we will need
  if (opts$model.is.branching) {
    stan.data$lengthscales = c(stan.data$l, l.z)  # Length scale for pseudotime and z
    stan.data$M <- num.inducing.tau * num.inducing.z
    cell.map <- cell.map %>% left_join(dplyr::select(pca.df, cell, z.hat))
  }
})

#' Compile the model and cache the DSO to avoid unnecessary recompilation.
#'
#' @param dl de.lorean object
#'
#' @export
#'
compile.model <- function(dl) {
  stan.model.file <-
    system.file(file.path('Stan', sprintf('%s.stan', dl$opts$model.name)),
                package = 'DeLorean',
                mustWork = TRUE)
  data.dir <- system.file('extdata', package='DeLorean')
  compiled.model.file <- paste(data.dir,
                               sprintf("%s.rds", dl$opts$model.name),
                               sep = '/')
  within(dl, {
    if (file.exists(compiled.model.file)
        &&
        file.info(compiled.model.file)$mtime
            > file.info(stan.model.file)$mtime)
    {
      message("Loading pre-compiled model from ", compiled.model.file)
      compiled <- readRDS(compiled.model.file)
    } else {
      message("Compiling model")
      compiled <- rstan::stan(file=stan.model.file, chains=0,
                              data=stan.data)
      message("Saving compiled model to ", compiled.model.file)
      saveRDS(compiled, compiled.model.file)
    }
    # Try one iteration to check everything is OK
    # message("Trying iteration")
    fit <- rstan::stan(
      fit = compiled,
      data = stan.data,
      init = list(init.from.prior(dl)),
      warmup = 1,
      iter = 1,
      chains = 1)
  })
}


# The covariance function for the DeLorean object.
cov.fn.for <- function(dl) {
    cov.matern.32
}


#' Perform all the steps necessary to fit the model.
#' - prepare the data
#' - compile the model
#' - find suitable initialisations
#' - fit the model using the specified method (sampling or variational Bayes)
#' - process the posterior.
#'
#' @param dl de.lorean object
#' @param method Fitting method:
#'   \itemize{
#'     \item 'sample': Use a Stan sampler.
#'       See \code{\link{fit.model.sample}}.
#'     \item 'vb': Use Stan ADVI variational Bayes algorithm.
#'       See \code{\link{fit.model.vb}}.
#'   }
#' @param ... Extra arguments for fitting method
#'
#' @export
#'
fit.dl <- function(
    dl,
    method = 'sample',
    ...)
{
  dl <- prepare.for.stan(dl)
  dl <- compile.model(dl)
  dl <- find.good.ordering(dl, seriation.find.orderings)
  dl <- pseudotimes.from.orderings(dl)
  dl <- fit.model(dl, method=method, ...)
  dl <- process.posterior(dl)
  dl <- analyse.noise.levels(dl)
}


#' Fit the model using specified method (sampling or variational Bayes).
#'
#' @param dl de.lorean object
#' @param method Fitting method:
#'   \itemize{
#'     \item 'sample': Use a Stan sampler.
#'       See \code{\link{fit.model.sample}}.
#'     \item 'vb': Use Stan ADVI variational Bayes algorithm.
#'       See \code{\link{fit.model.vb}}.
#'   }
#' @param ... Extra arguments for method
#'
#' @export
#'
fit.model <- function(
    dl,
    method = 'sample',
    ...)
{
    switch(
        method,
        "sample" = fit.model.sample(dl, ...),
        "vb" = fit.model.vb(dl, ...),
        stop('Unknown method'))
}


#' Fit the model using Stan sampler
#'
#' @param dl de.lorean object
#' @param num.cores Number of cores to run on.
#'   Defaults to getOption("DL.num.cores", max(parallel::detectCores()-1, 1))
#' @param chains Number of chains to run on each core
#' @param thin How many samples to generate before retaining one
#' @param ... Extra arguments for rstan::stan() sampling call
#'
#' @export
#'
fit.model.sample <- function(
    dl,
    num.cores = getOption("DL.num.cores", max(parallel::detectCores() - 1, 1)),
    chains = num.cores,
    thin = 50,
    ...)
{
  dl <- create.inits(dl, num.cores)
  dl$fit <-
    rstan::stan(
      fit = dl$fit,
      data = dl$stan.data,
      thin = thin,
      init = get.init.fn(dl),
      chains = chains,
      cores = num.cores,
      refresh = -1,
      ...)
  dl$compiled <- NULL  # Delete large unneeded object
  return(dl)
}


#' Average across a parameters samples.
#'
#' @param s An array of any dimension in which the first dimensions
#'   indexes the samples
#'
avg.par.samples <- function(s) {
  n.dims <- length(dim(s))
  if (n.dims > 1) {
    apply(s, 2:n.dims, mean)
  } else {
    mean(s)
  }
}

#' Get posterior mean of samples
#'
#' @param extract A named list of samples
#'
#' @export
#'
get.posterior.mean <- function(extract) lapply(extract, avg.par.samples)


#' Fit the model using Stan variational Bayes
#'
#' @param dl de.lorean object
#' @param num.cores Number of cores to run on. Defaults to default.num.cores()
#' @param num.inits Number initialisations to try. Defaults to num.cores
#' @param init.idx Which initialisation to use if only using one
#' @param ... Extra arguments for rstan::vb()
#'
#' @export
#'
fit.model.vb <- function(
    dl,
    num.cores = default.num.cores(),
    num.inits = num.cores,
    init.idx = 1,
    ...)
{
  dl <- create.inits(dl, num.cores)
  init.fn <- get.init.fn(dl)
  if (num.cores > 1) {
    #
    # Run variational Bayes in parallel
    sflist <- parallel::mclapply(
      1:num.inits,
      mc.cores = num.cores,
      # mc.cores = 1,
      function(i) within(list(), {
        fit <- rstan::vb(
          rstan::get_stanmodel(dl$fit),
          data = dl$stan.data,
          seed = i,
          init = init.fn(i),
          ...)
        pars <- get.posterior.mean(rstan::extract(fit))
        upars <- rstan::unconstrain_pars(fit, pars)
        lp <- rstan::log_prob(fit, upars, adjust_transform = TRUE)
        lp.unadj <- rstan::log_prob(fit, upars, adjust_transform = FALSE)
      }))
    #
    # Only keep results that worked
    sflist <- Filter(function(x) is.recursive(x) && ! is.null(x$lp), sflist)
    n.worked <- length(sflist)
    if (0 == n.worked) {
      stop('No VB fits worked.')
    } else if (n.worked < num.inits) {
      warning('Only ', n.worked, '/', num.inits, ' VB fits succeeded.')
    }
    #
    # Get the log likelihoods as vectors
    dl$vb.lls <- sapply(sflist, function(sf) sf$lp)
    dl$vb.lls.unadj <- sapply(sflist, function(sf) sf$lp.unadj)
    #
    # Calculate which run had best lp for posterior mean parameters
    dl$best.vb.id <- which.max(dl$vb.lls.unadj)
    #
    # Save the estimated tau for analysis
    dl$vb.tau <- sapply(sflist, function(sf) sf$pars$tau)
    #
    # Use those results
    dl$fit <- sflist[[dl$best.vb.id]]$fit
  } else {
    # Run single variational Bayes
    dl$fit <- rstan::vb(
      rstan::get_stanmodel(dl$fit),
      data = dl$stan.data,
      seed = init.idx,
      init = init.fn(init.idx),
      ...)
    dl$best.vb.id <- init.idx
  }
  return(dl)
}


#' Analyse the samples and gather the convergence statistics. Note this
#' only makes sense if a sampling method was used to fit the model as
#' opposed to variational Bayes.
#'
#' @param dl de.lorean object
#'
#' @export
#'
examine.convergence <- function(dl) {
    within(dl, {
        pars <- c("tau", "psi", "S", "omega")
        summ <- rstan::monitor(fit,
                        print=FALSE,
                        pars=pars)
        ignore.names <- stringr::str_detect(rownames(summ),
                                   "^(predictedvar|predictedmean)")
        rhat.sorted <- sort(summ[! ignore.names, "Rhat"])
        rhat.df <- data.frame(
            rhat=rhat.sorted,
            param=names(rhat.sorted),
            parameter=stringr::str_match(names(rhat.sorted), "^[[:alpha:]]+"))
        rm(summ)
        rm(pars)
    })
}


# The dimensions of the model parameters
#
# @param dl de.lorean object
#
model.parameter.dimensions <- function(dl) {
  sample.dims <- list(
    lp__ = c(),
    S = c("c"),
    tau = c("c"),
    z = c("c"),
    psi = c("g"),
    omega = c("g"),
    predictedmean = c("g", "t"),
    predictedvar = c("g", "t"),
    logmarglike = c("g")
  )
  if (! dl$opts$model.estimates.cell.sizes) {
    sample.dims$S <- NULL
  }
  if (! dl$opts$model.is.branching) {
    sample.dims$z <- NULL
  } else {
    sample.dims$logmarglike <- NULL
    sample.dims$predictedmean <- NULL
    sample.dims$predictedvar <- NULL
  }
  sample.dims
}


# Sample melter
#
# @param dl de.lorean object
#
sample.melter <- function(dl, include.iter = TRUE) {
    function(sample.list, sample.dims) {
        melt.var <- function(param) {
            # message(param)
            if (include.iter) {
                varnames <- c("iter", sample.dims[[param]])
            } else {
                varnames <- sample.dims[[param]]
            }
            melt(sample.list[[param]], varnames, value.name=param)
        }
        sapply(names(sample.dims), melt.var)
    }
}


# Join extra data to tau samples.
#
# @param dl de.lorean object
#
join.tau.samples <- function(dl, tau.samples) {
    with(dl,
         tau.samples
             %>% left_join(cell.map)
             %>% mutate(tau.offset=tau-obstime))
}


#' Process the posterior, that is extract and reformat the samples from
#' Stan. We also determine which sample has the highest likelihood, this
#' is labelled as the 'best' sample.
#'
#' @param dl de.lorean object
#'
#' @export
#'
process.posterior <- function(dl) {
  within(dl, {
    # Define a function to melt samples into a long format
    samples.l <- sample.melter(dl)(rstan::extract(dl$fit, permuted = TRUE),
                                   model.parameter.dimensions(dl))
    best.sample <- which.max(samples.l$lp__$lp__)
    if (TRUE %in% samples.l$logmarglike$is.held.out) {
      mean.held.out.marg.ll <- mean(
        (samples.l$logmarglike
        %>% left_join(gene.map)
        %>% filter(is.held.out))$logmarglike)
      message('Mean held out marginal log likelihood per cell: ',
              mean.held.out.marg.ll / stan.data$C)
    }
    # Include meta data in tau samples
    samples.l$tau <- join.tau.samples(dl, samples.l$tau)
  })
}


# Optimise the log posterior starting at a particular sample or from some
# other set of parameters.
#
# @param dl de.lorean object
# @param sample.iter Sample to optimise (defaults to best sample).
#
optimise.sample <- function(
    dl,
    parameters=sample.parameters(dl, sample.iter=sample.iter),
    sample.iter = dl$best.sample,
    ...)
{
    with(dl, {
        optimised <-
          rstan::optimizing(
            dl$fit@stanmodel,
            data = dl$stan.data,
            init = parameters,
            as_vector=FALSE,
            ...)
        dims <- model.parameter.dimensions(dl)
        # Don't melt lp__ value
        dims$lp__ <- NULL
        samples.l <- sample.melter(dl, include.iter=FALSE)(optimised$par, dims)
        # Include meta data in tau samples
        samples.l$tau <- join.tau.samples(dl, samples.l$tau)
        # Include lp__ value
        samples.l$lp__ <- data.frame(lp__=optimised$value)
        samples.l
    })
}


# Bind a sample.
#
# @param dl de.lorean object
# @param samples Samples to bind to existing samples.
# @param sample.iter Iteration (defaults to -1).
#
bind.sample <- function(dl, samples, sample.iter=-1)
{
    within(
        dl,
        for (param in names(samples.l)) {
            samples.l[[param]] <- rbind(samples[[param]]
                                            %>% mutate(iter=sample.iter),
                                        samples.l[[param]])
        })
}


#' Optimise the best sample and update the best.sample index.
#'
#' @param dl de.lorean object
#' @param sample.to.opt Sample to optimise
#' @param new.best.sample Update to best sample index
#'
#' @export
#'
optimise.best.sample <- function(
    dl,
    sample.to.opt=dl$best.sample,
    new.best.sample=-1)
{
    dl$sample.optimised <- sample.to.opt
    dl <- bind.sample(dl, optimise.sample(dl, sample.iter=sample.to.opt))
    dl$sample.old.best <- dl$best.sample
    dl$best.sample <- new.best.sample
    dl
}


#' Analyse noise levels and assess which genes have the greatest
#' ratio of temporal variance to noise. This are labelled as the
#' 'gene.high.psi' genes.
#'
#' @param dl de.lorean object
#' @param num.high.psi How many genes with high variance to examine
#'
#' @export
#'
analyse.noise.levels <- function(dl, num.high.psi=25) {
  within(dl, {
    noise.levels <- (
      with(samples.l, left_join(psi, omega))
      %>% left_join(gene.map))
    # Summarise by gene
    gene.noise.levels <- (
      noise.levels
      %>% group_by(g)
      %>% dplyr::summarise(omega=mean(omega), psi=mean(psi))
      %>% left_join(gene.map)
      %>% arrange(-psi/omega))
    genes.high.psi <- head(gene.noise.levels$gene, num.high.psi)
    # Calculate some statistics of the posterior of the gene parameters
    var.post.stats <-
      with(samples.l, left_join(psi, omega)) %>%
      group_by(g) %>%
      summarise(
        psi.mean   = mean(psi),
        psi.sd     = sd(psi),
        omega.mean = mean(omega),
        omega.sd   = sd(omega)) %>%
      left_join(gene.map)
  })
}


# Get the sampled parameter for the gene
#
# @param dl de.lorean object
# @param gene.idx Gene index
# @param param Parameter
# @param sample.iter Iteration to use (defaults to best.sample)
#
sampled.gene.param <- function(dl,
                               gene.idx,
                               param,
                               sample.iter = dl$best.sample) {
    with(dl, {
        filter(samples.l[[param]], gene.idx == g, sample.iter == iter)[[param]]
    })
}


#' Make predictions
#'
#' @param dl de.lorean object
#'
#' @export
#'
make.predictions <- function(dl) {
    within(dl, {
        predictions <- with(samples.l,
                            predictedmean
                            %>% left_join(predictedvar)
                            # %>% left_join(S)
                            %>% mutate(tau=test.input[t]))
    })
}


#' Fit held out genes
#'
#' @param dl de.lorean object
#' @param expr.held.out The expression matrix including the held out genes
#' @param sample.iter The sample to use to fit with
#'
fit.held.out <- function(
    dl,
    expr.held.out,
    sample.iter = dl$best.sample)
{
    with(dl, {
        if (opts$model.estimates.cell.sizes) {
            cell.posterior <- samples.l$S %>% filter(sample.iter == iter)
            expr.held.out <- t(t(expr.held.out) - cell.posterior$S)
        }
        #' Calculate covariance over pseudotimes and capture times
        calc.K <- functional::Curry(cov.calc.gene,
                        dl,
                        include.test=F,
                        psi=exp(stan.data$mu_psi),
                        omega=exp(stan.data$mu_omega))
        tau <- tau.for.sample(dl, sample.iter=sample.iter)
        obstime <- cell.map$obstime
        K.tau <- calc.K(tau=tau)
        K.capture <- calc.K(tau=obstime)
        #' Evaluate the held out gene under the GP model using pseudotimes
        #' and a model without.
        #'
        calc.gp.marginals <- function(expr) {
            c(
                gp.log.marg.like(expr, K.tau),
                gp.log.marg.like(expr, K.capture))
        }
        fit.model <- function(expr, model=loess) {
            list(tau=model(expr~s(tau)))
        }
        apply(expr.held.out, 1, calc.gp.marginals)
        # list(
            # gp.marginals=sapply(held.out.genes, calc.gp.marginals),
            # loess=lapply(held.out.genes, functional::Curry(fit.model, model=gam)))
            # gam=lapply(held.out.genes, functional::Curry(fit.model, model=gam)))
    })
}


# Parameter values for sample
#
# @param dl de.lorean object
# @param sample.iter The sample we want the parameters for.
#
sample.parameters <- function(dl,
                              sample.iter = dl$best.sample,
                              param.names = names(dl$samples.l))
{
    parameters <- lapply(param.names,
                         function(param)
                             filter(dl$samples.l[[param]],
                                    iter == sample.iter)[[param]])
    names(parameters) <- param.names
    parameters
}


#' Test fit for log normal and gamma
#'
#' @param vars Data to fit
#'
#' @export
#'
test.fit <- function(vars) {
    fit.gamma <- MASS::fitdistr(vars, 'gamma')
    fit.lognormal <- MASS::fitdistr(vars, 'lognormal')
    gp <- (
        ggplot(data.frame(V=vars), aes(x=V))
        + geom_density()
        + stat_function(fun=functional::Curry(dgamma,
                                shape=fit.gamma$estimate['shape'],
                                rate=fit.gamma$estimate['rate']),
                        linetype='dashed')
        + stat_function(fun=functional::Curry(dlnorm,
                                meanlog=fit.lognormal$estimate['meanlog'],
                                sdlog=fit.lognormal$estimate['sdlog']),
                        linetype='dotted')
    )
    list(gamma=fit.gamma, lognormal=fit.lognormal, gp=gp)
}


# The samples
#
# @param dl de.lorean object
#
sample.iters <- function(dl) dl$samples.l$lp__$iter


#' Make a fit valid by running one iteration of the sampler.
#'
#' @param dl de.lorean object
#'
#' @export
#'
make.fit.valid <- function(dl) rstan::stan(fit=dl$fit, data=dl$stan.data, iter=1, chains=1)
