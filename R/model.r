#' Update a factor based on the levels of another factor.
#'
#' @param .factor Factor
#' @param reference.factor Factor whose levels to use.
#' @param levels New levels
#'
#' @export
#'
update.levels <- function(.factor,
                          reference.factor=NULL,
                          .levels=levels(reference.factor),
                          ordered=is.ordered(reference.factor))
{
    factor(.factor, levels=.levels, ordered=ordered)
}


#' Retrieve the estimated tau for the given sample.
#'
#' @param dl de.lorean object
#' @param sample.iter Which sample to use, defaults to best sample
#'
#' @export
#'
tau.for.sample <- function(dl, sample.iter=dl$best.sample) {
    (
        dl$samples.l$tau
        %>% filter(sample.iter == iter)  # Filter correct iteration
        %>% arrange(c)  # Sort by cell
    )$tau
}


#' The levels of the gene factor
#'
#' @param dl The de.lorean object.
#'
#' @export
#'
gene.levels <- function(dl) levels=levels(dl$gene.meta$gene)


#' The levels of the cell factor
#'
#' @param dl The de.lorean object.
#'
#' @export
#'
cell.levels <- function(dl) levels=levels(dl$cell.meta$cell)


#' Melt an expression matrix.
#'
#' @param dl The de.lorean object.
#' @param expr Matrix of expression values.
#'
#' @export
#'
melt.expr <- function(dl, expr=dl$expr) (
    expr
    %>% melt(varnames=c("gene", "cell"), value.name="x")
    %>% mutate(gene=factor(gene, levels=gene.levels(dl)),
               cell=factor(cell, levels=cell.levels(dl)))
)


#' Cast an expression matrix.
#'
#' @param dl The de.lorean object.
#' @param expr.l Expression values in long format.
#'
#' @export
#'
cast.expr <- function(expr.l) expr.l %>% acast(gene ~ cell, value.var="x")


#' Analyse variance of expression between and within capture times.
#'
#' @param dl de.lorean object
#'
#' @export
#'
analyse.variance <- function(dl, adjust.cell.sizes) {
    within(dl, {
        #
        # First melt expression data into long format
        #
        expr.l <- melt(expr, varnames=c("gene", "cell"), value.name="x")
        expr.l$gene <- factor(expr.l$gene, levels=levels(gene.meta$gene))
        expr.l$cell <- factor(expr.l$cell, levels=levels(cell.meta$cell))
        #
        # Estimate a pseudo reference mean for each gene
        #
        gene.expr <- (expr.l
            %>% group_by(gene)
            %>% summarise(x.mean=mean(x))
        )
        stopifnot(! is.na(gene.expr))
        # Estimate the cell size by the median of the expression
        # adjust by the gene's mean
        cell.expr <- (expr.l
            %>% left_join(gene.expr)
            %>% group_by(cell)
            %>% summarise(S.hat=median(x - x.mean))
        )
        stopifnot(! is.na(cell.expr))
        expr.l <- expr.l %>% left_join(cell.expr)
        if (adjust.cell.sizes) {
            # Adjust the expression by the cell size estimates
            expr.l <- expr.l %>% mutate(x.hat=x - S.hat)
        } else {
            # Use the raw expression values
            expr.l <- expr.l %>% mutate(x.hat=x)
        }
        stopifnot(! is.na(expr.l))
        # Resummarise the adjusted expression data
        gene.expr <- (expr.l
            %>% group_by(gene)
            %>% summarise(x.mean=mean(x),
                          x.sd=sd(x),
                          phi.hat=mean(x.hat),
                          x.hat.sd=sd(x.hat))
            %>% filter(! is.na(x.sd), ! is.na(x.hat.sd))
        )
        stopifnot(! is.na(gene.expr))
        stopifnot(nrow(gene.expr) > 0)  # Must have some rows left
        # Examine the variation within and between times for each gene
        gene.time.expr <- (expr.l
            %>% left_join(cell.meta)
            %>% group_by(gene, capture)
            %>% summarise(x.mean=mean(x.hat),
                          x.var=mean((x.hat-mean(x.hat))**2),
                          num.capture=n())
            %>% filter(! is.na(x.var))
        )
        stopifnot(! is.na(gene.time.expr))
        stopifnot(nrow(gene.time.expr) > 0)  # Must have some rows left
        # Decomposition of variance within and between time.
        gene.var <- (gene.time.expr
            %>% group_by(gene)
            %>% summarise(Omega=weighted.mean(x.var, num.capture, na.rm=TRUE),
                          Psi=weighted.mean((x.mean-mean(x.mean))**2,
                                            num.capture, na.rm=TRUE))
            %>% filter(! is.na(Psi), Omega > 0 | Psi > 0)
        )
        stopifnot(! is.na(gene.var))
        stopifnot(nrow(gene.var) > 0)  # Must have some rows left
        # No longer needed
        rm(expr.l)
    })
}


#' Estimate hyperparameters for model using empirical Bayes
#'
#' @param dl de.lorean object
#' @param sigma.tau Noise s.d. in temporal dimension
#' @param length.scale Length scale for stationary GP covariance function
#'
#' @export
#'
estimate.hyper <- function(
    dl,
    sigma.tau = .5,
    delta = .5,
    length.scale = NULL,
    model.name = 'simplest-model'
) {
    dl <- within(dl, {
        #
        # Remember options that depend on the model
        #
        opts$model.name <- model.name
        opts$estimate.phi <- switch(
            opts$model.name,
            "simple-model" = TRUE,
            "simplest-model" = FALSE,
            NA)
        opts$adjust.cell.sizes <- switch(
            opts$model.name,
            "simple-model" = TRUE,
            "simplest-model" = FALSE,
            NA)
        #
        # Set up temporal hyper-parameters
        #
        opts$delta <- delta
        opts$sigma.tau <- sigma.tau
        time.range <- range(cell.meta$obstime)
        time.width <- time.range[2] - time.range[1]
        if (is.null(length.scale)) {
            opts$length.scale <- time.width
        } else {
            opts$length.scale <- length.scale
        }
        # message("Length scale: ", opts$length.scale)
    })
    dl <- analyse.variance(dl, dl$opts$adjust.cell.sizes)
    within(dl, {
        #
        # Work out the variance that would be expected from
        # a sample from the GP using a sample of tau from the prior
        #
        # First, which cells are we using?
        cells <- (
            cell.meta
            # Select cells in use
            %>% filter(cell %in% colnames(expr))
            # Sample tau from prior
            %>% mutate(tau.sample=rnorm(n(), mean=obstime, sd=sigma.tau)))
        # Calculate the covariance over the pseudotimes
        K.psi <- cov.matern.32(cov.calc.dl.dists(dl,
                                                 tau=cells$tau.sample,
                                                 include.test=FALSE),
                               opts$length.scale)
        colnames(K.psi) <- rownames(K.psi) <- cells$cell
        # The overall expected sample variance from the covariance over
        # the pseudotimes
        V.psi <- expected.sample.var(K.psi)
        # message('V.psi: ', V.psi)
        # The expected sample variance within each capture time
        # var.capture <- (
            # cells
            # %>% group_by(capture)
            # %>% summarise(exp.sample.var=expected.sample.var(K.psi[as.character(cell),as.character(cell)]),
                          # num.capture=n()))
        # V.omega <- with(var.capture, weighted.mean(exp.sample.var,
                                                   # num.capture))
        # message('V.omega: ', V.omega)
        # message('V ratio: ', (V.psi - V.omega) / V.omega)
        gene.var <- (
            gene.var
            %>% mutate(var.ratio=Omega/Psi,
                       lambda=(1+var.ratio)/(V.psi+var.ratio),
                       psi.hat=lambda*Psi,
                       omega.hat=lambda*Omega))
        hyper <- list(
            mu_S=mean(cell.expr$S.hat),
            sigma_S=sd(cell.expr$S.hat),
            mu_phi=mean(gene.expr$phi.hat),
            sigma_phi=sd(gene.expr$phi.hat),
            mu_psi=mean(log(gene.var$psi.hat), na.rm=TRUE),
            sigma_psi=sd(log(gene.var$psi.hat), na.rm=TRUE),
            mu_omega=mean(log(gene.var$omega.hat), na.rm=TRUE),
            sigma_omega=sd(log(gene.var$omega.hat), na.rm=TRUE),
            sigma_tau=opts$sigma.tau,
            l=opts$length.scale
        )
        stopifnot(all(! sapply(hyper, is.na)))
        # No longer needed
        rm(V.psi, K.psi, cells)
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
                         .filter=function(x) x %in% genes,
                         number=NULL,
                         genes=sample(rownames(dl$expr), number))
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
                         .filter=function(x) x %in% cells,
                         number=NULL,
                         cells=sample(colnames(dl$expr), number))
{
    within(dl, {
        expr <- expr[,.filter(colnames(expr))]
        message("Have ", ncol(expr), " cells after filtering")
    })
}


#' Sample genes and cells
#'
#' @param dl de.lorean object
#'
#' @export
#'
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


#' Format for Stan
#'
#' @param dl de.lorean object
#' @param num.test Number of test points to consider
#' @param period Period of expression patterns
#' @param hold.out Number genes to hold out for generalisation tests
#'
format.for.stan <- function(
    dl,
    num.test = 101,
    period = 0,
    hold.out = 0
) {
    within(dl, {
        opts$num.test <- num.test
        opts$period <- period
        opts$periodic <- opts$period > 0
        stopifnot(hold.out < nrow(expr))
        .G <- nrow(expr) - hold.out
        .C <- ncol(expr)
        #
        # Permute genes to make held out genes random
        stan.m <- expr[sample(.G+hold.out),]
        #
        # Calculate the map from gene indices to genes and their meta data
        gene.map <- (data.frame(g=1:(.G+hold.out),
                                gene=factor(rownames(stan.m),
                                            levels=levels(gene.meta$gene)))
                    %>% mutate(is.held.out=g>.G)
                    %>% left_join(gene.expr)
                    %>% left_join(gene.var)
                    %>% left_join(gene.meta))
        stopifnot(! is.na(gene.map[
            c("g", "gene", "x.mean", "x.sd",
            "phi.hat", "x.hat.sd", "Omega", "Psi",
            "omega.hat", "psi.hat")
        ]))
        #
        # Rename phi.hat if we are not estimating phi
        if (! opts$estimate.phi) {
            gene.map <- gene.map %>% rename(phi=phi.hat)
        }
        #
        # Calculate the map from cell indices to genes and their meta data
        cell.map <- (data.frame(c=1:.C,
                                cell=factor(colnames(stan.m),
                                            levels=levels(cell.meta$cell)))
                    %>% left_join(cell.meta)
                    %>% left_join(cell.expr))
        stopifnot(! is.na(cell.map %>% dplyr::select(cell, capture, obstime)))
        #
        # Calculate the time points at which to make predictions
        test.input <- (
            time.range[1] - 2 * opts$sigma.tau
            + (time.width + 4 * opts$sigma.tau)
                * (0:(opts$num.test-1)) / (opts$num.test-1))
        #
        # Gather all the data into one list
        stan.data <- c(
            # Hyper-parameters
            hyper,
            list(
                # Dimensions
                C=.C,
                G=.G,
                H=hold.out,
                # Data
                time=cell.map$obstime,
                expr=stan.m,
                # Held out parameters
                heldout_psi=filter(gene.map, g > .G)$psi.hat,
                heldout_omega=filter(gene.map, g > .G)$omega.hat,
                heldout_phi=filter(gene.map, g > .G)$phi.hat,
                # Generated quantities
                numtest=opts$num.test,
                testinput=test.input,
                # Periodic?
                periodic=opts$periodic,
                period=opts$period
            )
        )
        # If we're not estimating phi, add it to the data
        if (! opts$estimate.phi) {
            stan.data$phi <- gene.map$x.mean
        }
    })
}

#' Define and compile the model.
#'
#' @param dl de.lorean object
#'
#' @export
#'
compile.model <- function(dl) {
    stan.model.file <- system.file(sprintf('inst/Stan/%s.stan',
                                           dl$opts$model.name),
                                   package='DeLorean')
    stopifnot(! is.null(stan.model.file))
    stopifnot("" != stan.model.file)
    data.dir <- system.file('data', package='DeLorean')
    compiled.model.file <- paste(data.dir,
                                 sprintf("%s.rds", dl$opts$model.name),
                                 sep='/')
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
            compiled <- stan(file=stan.model.file, chains=0)
            message("Saving compiled model to ", compiled.model.file)
            saveRDS(compiled, compiled.model.file)
        }
        # Try one iteration to check everything is OK
        # message("Trying iteration")
        fit <- stan(fit=compiled,
                    data=stan.data,
                    init=make.chain.init.fn(dl),
                    iter=1,
                    chains=1)
    })
}


#' Define a function to initialise the chains
#'
#' @param dl de.lorean object
#'
#' @export
#'
make.chain.init.fn <- function(dl) {
    function() {
        with(dl$stan.data, {
            # message("Creating initialisation")
            init <- list(
                delta=runif(1),
                S=dl$cell.map$S.hat,
                tau=rnorm(C, mean=time, sd=sigma_tau),
                psi=rlnorm(G, meanlog=mu_psi, sdlog=sigma_psi),
                omega=rlnorm(G, meanlog=mu_omega, sdlog=sigma_omega)
            )
            # If estimating phi, include it.
            if (dl$opts$estimate.phi) {
                init$phi <- rnorm(G, mean=mu_phi, sd=sigma_phi)
            }
            init
        })
    }
}


even.tau.spread <- function(dl) {
    with(dl$stan.data,
         seq(min(time) - sigma_tau,
             max(time) + sigma_tau,
             length=C))
}


#' Find best order of the samples assuming some smooth GP prior on the
#' expression profiles over this ordering.
#'
#' @param dl de.lorean object
#' @param use.parallel Calculate in parallel
#' @param num.cores Number of cores to run on.
#'          Defaults to getOption("DL.num.cores", max(detectCores()-1, 1))
#' @param num.tau.to.keep How many initialisations to keep.
#' @param method Method to use "maximise" or "metropolis"
#'
#' @export
#'
find.smooth.tau <- function(
    dl,
    psi = exp(dl$hyper$mu_psi),
    omega = exp(dl$hyper$mu_omega),
    use.parallel = TRUE,
    num.cores = getOption("DL.num.cores", max(detectCores() - 1, 1)),
    num.tau.to.try = num.cores,
    num.tau.to.keep = num.cores,
    method = "metropolis",
    ...
) {
    log.likelihood <- ordering.log.likelihood.fn(dl)
    dl$tau.inits <- with(dl, {
        # Maximise the sum of the log marginal likelihoods
        ordering.search <- function(seed) {
            set.seed(seed)
            # Choose a starting point by random projection
            expr.centre <- t(scale(t(expr), center=T, scale=F))
            init.ordering <- order(rnorm(nrow(expr.centre)) %*% expr.centre)
            # init.ordering <- sample(stan.data$C)
            metropolis.fn <- function(ordering, log.likelihood, ...) {
                mh.run <- ordering.metropolis.hastings(
                    ordering,
                    log.likelihood,
                    proposal.fn=ordering.random.block.move,
                    ...)
                best.sample <- which.max(mh.run$log.likelihoods)
                #ordering.maximise(mh.run$chain[best.sample,], log.likelihood)
                mh.run$chain[best.sample,]
            }
            method.fn <- switch(method,
                                "maximise"=ordering.maximise,
                                "metropolis"=metropolis.fn,
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
        if (use.parallel) {
            orderings <- mclapply(seeds,
                                  mc.cores=num.cores,
                                  ordering.search)
        } else {
            orderings <- lapply(seeds, ordering.search)
        }
        # Order the taus by the best orderings
        lls <- sapply(orderings, function(o) -log.likelihood(o))
        best.order <- order(lls)
        # Make the complete chain initialisation with the tau.
        lapply(orderings[best.order[1:num.tau.to.keep]],
               function(ordering) {
                   init <- init.chain.sample.tau(dl)
                   init$tau <- even.tau.spread(dl)[ordering.invert(ordering)]
                   init
               })
    })
    dl
}


#' Test ordering Metropolis-Hastings sampler.
#'
#' @param dl de.lorean object
#' @param use.parallel Calculate in parallel
#' @param num.cores Number of cores to run on.
#'          Defaults to getOption("DL.num.cores", max(detectCores()-1, 1))
#'
#' @export
#'
test.mh <- function(
    dl,
    psi = mean(dl$gene.map$psi.hat),
    omega = mean(dl$gene.map$omega.hat),
    use.parallel = TRUE,
    num.cores = getOption("DL.num.cores", max(detectCores() - 1, 1)),
    iterations = 1000,
    thin = 15
) {
    log.likelihood <- ordering.log.likelihood.fn(dl)
    dl$tau.inits <- with(dl, {
        # Maximise the sum of the log marginal likelihoods
        ordering.search <- function(seed) {
            set.seed(seed)
            # Choose a random starting point
            init.ordering <- sample(stan.data$C)
            mh.run <- ordering.metropolis.hastings(
                init.ordering,
                log.likelihood,
                proposal.fn=ordering.random.block.move,
                iterations=iterations,
                thin=thin)
        }
        # Choose seeds
        seeds <- sample.int(.Machine$integer.max, num.cores)
        # Run in parallel or not?
        if (use.parallel) {
            orderings <- mclapply(seeds,
                                  mc.cores=num.cores,
                                  ordering.search)
        } else {
            orderings <- lapply(seeds, ordering.search)
        }
        orderings
    })
}


#' Create a log likelihood function suitable for evaluating smooth orderings.
#'
#' @param dl de.lorean object
#'
#' @export
#'
ordering.log.likelihood.fn <- function(
    dl,
    psi = mean(dl$gene.map$psi.hat),
    omega = mean(dl$gene.map$omega.hat),
    cov.fn = cov.matern.32)
{
    with(dl, {
        # Evenly spread tau over range of capture times
        even.tau <- even.tau.spread(dl)
        # Calculate the distances
        r <- outer(even.tau, even.tau, FUN="-")
        # Make periodic if necessary
        if (opts$periodic) {
            r <- cov.periodise(r, opts$period)
        }
        # Use the same kernel for each gene
        K <- (
            psi * cov.fn(r, stan.data$l)
            + omega * identity.matrix(nrow(r)))
        # Do Cholesky decomposition once and use in each subsequent smoothing
        U <- chol(K)
        # Make every gene zero mean
        expr.centre <- t(scale(t(dl$expr), scale=FALSE, center=TRUE))
        # Return the function that is the log likelihood of the ordering
        function(o) {
            stopifnot(! is.null(o))
            tryCatch({
                sum(sapply(1:stan.data$G,
                           function(g) {
                               gp.log.marg.like(expr.centre[g,o], U=U)
                           }))
            },
            warning=function(condition) {
                message(
                    "WARNING: log.likelihood\n\t",
                    condition,
                    "\n",
                    o)
                warning(condition)
            }, error = function(condition) {
                message(
                    "ERROR: log.likelihood\n\t",
                    condition,
                    "\n\tOrdering: ",
                    o)
                stop(condition)
            })
        }
    })
}


#' Choose an initialisation by sampling tau from the prior.
#'
#' @param dl de.lorean object
#'
#' @export
#'
init.chain.sample.tau <- function(dl) {
    with(dl$stan.data, {
        init <- list(
            alpha=dl$cell.map$alpha.hat,
            beta=rep(0, G),
            S=dl$cell.map$S.hat,
            tau=rnorm(C, time, sd=sigma_tau),
            phi=dl$gene.map$phi.hat[1:G],
            psi=dl$gene.map$psi.hat[1:G],
            omega=dl$gene.map$omega.hat[1:G],
            delta=.5
        )
        # If not estimating phi, don't include it.
        if (! dl$opts$estimate.phi) {
            init$phi <- NULL
        }
        init
    })
}


#' Find best tau to initialise chains with by using empirical Bayes parameter
#' estimates and sampling tau from its prior.
#'
#' @param dl de.lorean object
#' @param num.tau.candidates How many candidates to examine. Defaults to 6000.
#' @param num.tau.to.keep How many candidates to keep. Defaults to num.cores.
#' @param use.parallel Calculate in parallel
#' @param num.cores Number of cores to run on.
#'          Defaults to getOption("DL.num.cores", max(detectCores()-1, 1))
#'
#' @export
#'
find.best.tau <- function(dl,
                          num.tau.candidates = 6000,
                          num.tau.to.keep = NULL,
                          use.parallel = TRUE,
                          num.cores = getOption("DL.num.cores",
                                                max(detectCores() - 1, 1))
) {
    if (is.null(num.tau.to.keep)) {
        num.tau.to.keep <- num.cores
    }
    within(dl, {
           # Define a function that calculates log probability for random seeded tau
           try.tau.init <- function(i) {
               set.seed(i)
               pars <- init.chain.sample.tau(dl)
               list(lp=log_prob(fit, unconstrain_pars(fit, pars)),
                    tau=pars$tau)
           }
           # Choose tau several times and calculate log probability
           if (use.parallel) {
               tau.inits <- mclapply(1:num.tau.candidates,
                                     mc.cores=num.cores,
                                     try.tau.init)
           } else {
               tau.inits <- lapply(1:num.tau.candidates, try.tau.init)
           }
           # qplot(sapply(tau.inits, function(init) init$lp))
           # Which tau gave highest log probability?
           tau.inits.order <- order(sapply(tau.inits, function(init) -init$lp))
           # Just keep so many best tau inits
           tau.inits <- tau.inits[tau.inits.order[1:num.tau.to.keep]]
           rm(tau.inits.order, try.tau.init)
    })
}


#' Fit the model
#'
#' @param dl de.lorean object
#' @param num.cores Number of cores to run on.
#'          Defaults to getOption("DL.num.cores", max(detectCores()-1, 1))
#' @param chains Number of chains to run on each core
#' @param iter Number of iterations in each chain
#' @param thin How many samples to generate before retaining one
#'
#' @export
#'
fit.model <- function(
    dl,
    num.cores = getOption("DL.num.cores", max(detectCores() - 1, 1)),
    chains = 1,
    iter = 1000,
    thin = 50)
{
    init.chain.good.tau <- function(chain_id) {
        #
        # Create random parameters
        pars <- make.chain.init.fn(dl)()
        #
        # Replace tau with good tau
        pars$tau <- dl$tau.inits[[chain_id]]$tau
        pars
    }
    # Run the chains in parallel
    sflist <- mclapply(1:num.cores,
                       mc.cores=num.cores,
                       function(i)
                           stan(fit=dl$fit, data=dl$stan.data,
                                thin=thin,
                                init=init.chain.good.tau,
                                iter=iter,
                                seed=i, chains=chains,
                                chain_id=i, refresh=-1))
    dl$fit <- sflist2stanfit(sflist)
    dl$compiled <- NULL
    dl$sflist <- NULL
    return(dl)
}


#' Examine convergence
#'
#' @param dl de.lorean object
#'
#' @export
#'
examine.convergence <- function(dl) {
    within(dl, {
        pars <- c("tau", "psi", "S", "omega")
        if (opts$estimate.phi) {
            pars <- c(pars, "phi")
        }
        summ <- monitor(fit,
                        print=FALSE,
                        pars=pars)
        ignore.names <- str_detect(rownames(summ),
                                   "^(predictedvar|predictedmean)")
        rhat.sorted <- sort(summ[! ignore.names, "Rhat"])
        rm(summ)
        rm(pars)
    })
}


#' The dimensions of the model parameters
#'
#' @param dl de.lorean object
#'
#' @export
#'
model.parameter.dimensions <- function(dl) {
    sample.dims <- list(
        lp__=c(),
        delta=c(),
        S=c("c"),
        tau=c("c"),
        phi=c("g"),
        psi=c("g"),
        omega=c("g"),
        predictedmean=c("g", "t"),
        predictedvar=c("g", "t"),
        logmarglike=c("g")
    )
    if (! dl$opts$estimate.phi) {
        sample.dims$phi <- NULL
        sample.dims$S <- NULL
    }
    sample.dims
}


#' Sample melter
#'
#' @param dl de.lorean object
#'
#' @export
#'
sample.melter <- function(dl, include.iter=TRUE) {
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


#' Join extra data to tau samples.
#'
#' @param dl de.lorean object
#'
#' @export
#'
join.tau.samples <- function(dl, tau.samples) {
    with(dl,
         tau.samples
             %>% left_join(cell.map)
             %>% mutate(tau.offset=tau-obstime))
}


#' Process the posterior
#'
#' @param dl de.lorean object
#'
#' @export
#'
process.posterior <- function(dl) {
    within(dl, {
        # Define a function to melt samples into a long format
        samples.l <- sample.melter(dl)(extract(dl$fit, permuted=TRUE),
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


#' Optimise the log posterior starting at a particular sample or from some
#' other set of parameters.
#'
#' @param dl de.lorean object
#' @param sample.iter Sample to optimise (defaults to best sample).
#'
#' @export
#'
optimise.sample <- function(
    dl,
    parameters=sample.parameters(dl, sample.iter=sample.iter),
    sample.iter=dl$best.sample,
    ...)
{
    with(dl, {
        optimised <- optimizing(dl$fit@stanmodel,
                                data=dl$stan.data,
                                init=parameters,
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


#' Bind a sample.
#'
#' @param dl de.lorean object
#' @param samples Samples to bind to existing samples.
#' @param sample.iter Iteration (defaults to -1).
#'
#' @export
#'
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


#' Analyse noise levels
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
            %>% summarise(omega=mean(omega), psi=mean(psi))
            %>% left_join(gene.map)
            %>% arrange(1/psi))
        genes.high.psi <- head(gene.noise.levels$gene, num.high.psi)
    })
}


#' Get the sampled parameter for the gene
#'
#' @param dl de.lorean object
#' @param gene.idx Gene index
#' @param param Parameter
#' @param sample.iter Iteration to use (defaults to best.sample)
#'
sampled.gene.param <- function(dl,
                               gene.idx,
                               param,
                               sample.iter=dl$best.sample) {
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
        if (opts$estimate.phi) {
            predictions <- predictions %>% left_join(samples.l$phi)
        }
    })
}


#' Permute cells and test roughness of expression.
#'
#' @param dl de.lorean object
#' @param expr.held.out The expression matrix of the held out genes
#'
#' @export
#'
permuted.roughness <- function(
    dl,
    expr.held.out)
{
    permuted <- permute.df(dl$cell.map, "capture")
    apply(expr.held.out[,permuted$c], 1, calc.roughness)
}


#' Fit held out genes
#'
#' @param dl de.lorean object
#' @param expr.held.out The expression matrix including the held out genes
#'
#' @export
#'
fit.held.out <- function(
    dl,
    expr.held.out,
    sample.iter=dl$best.sample)
{
    library(mgcv)
    with(dl, {
        if (opts$adjust.cell.sizes) {
            cell.posterior <- samples.l$S %>% filter(sample.iter == iter)
            expr.held.out <- t(t(expr.held.out) + cell.posterior$S)
        }
        #' Calculate covariance over pseudotimes and capture times
        calc.K <- Curry(cov.calc.gene,
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
            # loess=lapply(held.out.genes, Curry(fit.model, model=gam)))
            # gam=lapply(held.out.genes, Curry(fit.model, model=gam)))
    })
}


#' Parameter values for sample
#'
#' @param dl de.lorean object
#' @param sample.iter The sample we want the parameters for.
#'
#' @export
#'
sample.parameters <- function(dl,
                              sample.iter=dl$best.sample,
                              param.names=names(dl$samples.l))
{
    parameters <- lapply(param.names,
                         function(param)
                             filter(dl$samples.l[[param]],
                                    iter == sample.iter)[[param]])
    names(parameters) <- param.names
    parameters
}
