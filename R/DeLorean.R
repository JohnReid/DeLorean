
#' Initialise DeLorean object
#'
#' @export
#'
de.lorean <- function(expr, gene.meta, cell.meta) {
    stopifnot("gene" %in% names(gene.meta))
    stopifnot(is.factor(gene.meta$gene))
    stopifnot("cell" %in% names(cell.meta))
    stopifnot("capture" %in% names(cell.meta))
    stopifnot("obstime" %in% names(cell.meta))
    stopifnot(is.factor(cell.meta$cell))
    stopifnot(is.factor(cell.meta$capture))
    stopifnot(is.numeric(cell.meta$obstime))
    stopifnot(nrow(expr) == nrow(gene.meta))
    stopifnot(ncol(expr) == nrow(cell.meta))
    stopifnot(! is.null(rownames(expr)))
    stopifnot(! is.null(colnames(expr)))
    stopifnot(all(rownames(expr) == gene.meta$gene))
    stopifnot(all(colnames(expr) == cell.meta$cell))
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
        convergence=plot.convergence(dl, ...)
    )
    if (is.null(result)) {
        stop('Unknown plot type')
    }
    result
}


#' The filename of the R markdown report.
#'
#' @param report.name The report name
#'
#' @export
report.file <- function(report.name) {
    system.file("inst", "Rmd", sprintf("%s.Rmd", report.name),
                package="DeLorean")
}


#' The filename of the R markdown stylesheet
#'
#' @export
#'
de.lorean.stylesheet <- function() {
    system.file("inst", "Rmd", "foghorn.css", package="DeLorean")
}


#' Knit a report, the file inst/Rmd/<report.name>.Rmd must exist in
#' the package directory.
#'
#' @param dl de.lorean object
#' @param report.name The name of the report. Used to locate the R
#'  markdown report file in the package.
#'
#' @export
#'
knit.report <- function(dl, report.name) {
    knit2html(report.file(report.name),
              envir=environment(),
              stylesheet=de.lorean.stylesheet())
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

#' Partition de.lorean object by cells
#'
#' @param dl de.lorean object
#' @param pieces How many pieces to partition into
#'
partition.de.lorean <- function(
    dl,
    pieces = 2
) {
    partition <- (1:(dim(dl)[2]) %% pieces) + 1
    cells <- sample(colnames(dl$expr))
    get.piece <- function(p) {
        partition.cells <- cells[partition == p]
        cell.filter <- function(cells) cells %in% partition.cells
        filter.cells(dl, cell.filter)
    }
    lapply(1:pieces, get.piece)
}

#' Test robustness of pseudotime estimation on subsets of de.lorean
#' object
#'
#' @param dl de.lorean object
#' @param pieces How many pieces to partition into
#'
test.robustness.de.lorean <- function(
    dl,
    pieces = 2
) {
    # Partition de.lorean object into several pieces
    partition <- partition.de.lorean(dl, pieces)
    # Define function to fit each piece
    run.model <- function(piece) {
        piece <- format.for.stan(piece)
        piece <- compile.model.simple(piece)
        piece <- find.best.tau(piece)
        piece <- fit.model(piece)
        process.posterior(piece)
    }
    # Fit full de.lorean object
    full.model <- run.model(dl)
    # Get tau posterior
    full.tau <- (full.model$samples.l$tau
                    %>% dplyr::select(-c)
                    %>% mutate(fit="part"))
    # Fit each piece
    piece.models <- lapply(partition, run.model)
    # Get tau posterior from each piece and adjust mean to match full fit
    # mean
    get.tau.posterior <- function(piece.model) {
        posterior <- (piece.model$samples.l$tau
            %>% dplyr::select(-c)
            %>% mutate(fit="full"))
        cells <- unique(posterior$cell)
        full.mean <- mean((full.tau %>% filter(cell %in% cells))$tau)
        piece.mean <- mean(posterior$tau)
        posterior$tau <- posterior$tau - piece.mean + full.mean
        posterior
    }
    # Gather tau posteriors from each piece
    pieces.tau <- do.call(rbind,
                          lapply(piece.models, get.tau.posterior))
    # Combine fit pieces with fit full model
    tau.posterior <- rbind(full.tau, pieces.tau)
    # Sort by median tau
    cells <- (tau.posterior
        %>% group_by(cell)
        %>% summarise(median.tau=median(tau))
        %>% arrange(median.tau)
    )
    # Create result list with plot
    list(tau.posterior = tau.posterior,
         gp = (ggplot(tau.posterior
                      %>% mutate(cell=factor(cell, levels=cells$cell)),
                      aes(x=cell, y=tau, fill=fit))
                  + geom_boxplot()
              ))
}

#' Estimate hyperparameters for model using empirical Bayes
#'
#' @param dl de.lorean object
#' @param sigma.tau Noise s.d. in temporal dimension
#' @param delta Proportion of within time variance to relabel as between time
#' @param min.sd Minimum s.d. used for drop-out effects (to avoid s.d. of 0 when no drop outs)
#' @param length.scale Length scale for stationary GP covariance function
#'
#' @export
#'
estimate.hyper <- function(
    dl,
    sigma.tau = .5,
    delta = .5,
    min.sd = 1e-10,
    length.scale = NULL
) {
    within(dl, {
        #
        # Set up temporal hyper-parameters
        #
        opts$delta <- delta
        opts$sigma.tau <- sigma.tau
        time.range <- range(cell.meta$obstime)
        time.width <- time.range[2] - time.range[1]
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
        # Adjust the expression by the cell size estimates
        expr.l <- (expr.l
            %>% left_join(cell.expr)
            %>% mutate(x.hat=x - S.hat)
        )
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
        # Examine the variation within genes and times
        gene.time.expr <- (expr.l
            %>% left_join(cell.meta)
            %>% group_by(gene, capture)
            %>% summarise(x.mean=mean(x.hat),
                          x.var=var(x.hat))
            %>% filter(! is.na(x.var))
        )
        stopifnot(! is.na(gene.time.expr))
        stopifnot(nrow(gene.time.expr) > 0)  # Must have some rows left
        # Decomposition of variance within and between time.
        gene.var <- (gene.time.expr
            %>% group_by(gene)
            %>% summarise(omega.bar=mean(x.var, na.rm=TRUE),
                          psi.bar=var(x.mean, na.rm=TRUE))
            %>% filter(! is.na(psi.bar))
            %>% mutate(within.time.mislabelled = opts$delta * omega.bar,
                       omega.hat = omega.bar - within.time.mislabelled,
                       psi.hat = psi.bar + within.time.mislabelled)
        )
        stopifnot(! is.na(gene.var))
        stopifnot(nrow(gene.var) > 0)  # Must have some rows left
        if (is.null(length.scale)) {
            length.scale <- time.width / 2
        }
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
            l=length.scale
        )
        rm(expr.l)  # No longer needed
    })
}


#' Filter genes
#'
#' @param dl de.lorean object
#'
#' @export
#'
filter.genes <- function(dl, gene.filter) {
    within(dl, {
        if( ! is.null(gene.filter) ) {
            # .genes.filtered <- rownames(expr)[opts$gene.filter(rownames(expr))]
            expr <- expr[gene.filter(rownames(expr)),]
            # rownames(expr) <- .genes.filtered
            message("Have ", nrow(expr), " genes after filtering")
        }
    })
}


#' Filter cells
#'
#' @param dl de.lorean object
#'
#' @export
#'
filter.cells <- function(dl, cell.filter) {
    within(dl, {
        if( ! is.null(cell.filter) ) {
            expr <- expr[,cell.filter(colnames(expr))]
            message("Have ", ncol(expr), " cells after filtering")
        }
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
            "phi.hat", "x.hat.sd", "omega.bar", "psi.bar",
            "within.time.mislabelled", "omega.hat", "psi.hat")
        ]))
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
                heldout_phi=filter(gene.map, g > .G)$phi.hat,
                heldout_psi=filter(gene.map, g > .G)$psi.hat,
                heldout_omega=filter(gene.map, g > .G)$omega.hat,
                # Generated quantities
                numtest=opts$num.test,
                testinput=test.input,
                # Periodic?
                periodic=opts$periodic,
                period=opts$period
            )
        )
    })
}

#' Define and compile a simple model without dropout
#'
#' @param dl de.lorean object
#'
#' @export
#'
compile.model.simple <- function(dl) {
    stan.model.file <- system.file('inst/Stan/simple-model.stan',
                                   package='DeLorean')
    data.dir <- system.file('data', package='DeLorean')
    compiled.model.file <- paste(data.dir, "simple-model.rds", sep='/')
    within(dl, {
        if (file.exists(compiled.model.file)) {
            # message("Loading pre-compiled model from ", compiled.model.file)
            compiled <- readRDS(compiled.model.file)
        } else {
            # message("Compiling model")
            compiled <- stan(file=stan.model.file, chains=0)
            # message("Saving compiled model to ", compiled.model.file)
            saveRDS(compiled, compiled.model.file)
        }
        # Try one iteration to check everything is OK
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
            list(S=dl$cell.map$S.hat,
                 tau=rnorm(C, mean=time, sd=sigma_tau),
                 phi=rnorm(G, mean=mu_phi, sd=sigma_phi),
                 psi=rlnorm(G, meanlog=mu_psi, sdlog=sigma_psi),
                 omega=rlnorm(G, meanlog=mu_omega, sdlog=sigma_omega)
            )
        })
    }
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
find.best.tau <- function(
    dl,
    num.tau.candidates = 6000,
    num.tau.to.keep = NULL,
    use.parallel = TRUE,
    num.cores = getOption("DL.num.cores", max(detectCores() - 1, 1))
) {
    if (is.null(num.tau.to.keep)) {
        num.tau.to.keep <- num.cores
    }
    within(dl, {
        # Define a function that chooses tau
        init.chain.find.tau <- function() {
            with(stan.data, {
                list(alpha=cell.map$alpha.hat,
                     beta=rep(0, G),
                     S=cell.map$S.hat,
                     tau=rnorm(C, time, sd=sigma_tau),
                     phi=gene.map$phi.hat[1:G],
                     psi=gene.map$psi.hat[1:G],
                     omega=gene.map$omega.hat[1:G]
                )
            })
        }
        # Define a function that calculates log probability for random seeded tau
        try.tau.init <- function(i) {
            set.seed(i)
            pars <- init.chain.find.tau()
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
        rm(tau.inits.order, init.chain.find.tau, try.tau.init)
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
        # print(chain_id)
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
        summ <- monitor(fit,
                        print=FALSE,
                        pars=c("tau", "psi", "S", "phi", "omega"))
        ignore.names <- str_detect(rownames(summ),
                                   "^(predictedvar|predictedmean)")
        rhat.sorted <- sort(summ[! ignore.names, "Rhat"])
        rm(summ)
    })
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
        melt.samples <- function(sample.list, sample.dims) {
            melt.var <- function(var.name) {
                melt(sample.list[[var.name]],
                    c("iter", sample.dims[[var.name]]),
                    value.name=var.name)
            }
            sapply(names(sample.dims), melt.var)
        }
        # The dimensions of each set of samples
        sample.dims <- list(
            lp__=c(),
            S=c("c"),
            tau=c("c"),
            phi=c("g"),
            psi=c("g"),
            omega=c("g"),
            predictedmean=c("g", "t"),
            predictedvar=c("g", "t"),
            logmarglike=c("g")
        )
        samples.l <- melt.samples(extract(dl$fit, permuted=TRUE),
                                  sample.dims)
        best.sample <- which.max(samples.l$lp__$lp__)
        if (TRUE %in% samples.l$logmarglike$is.held.out) {
            mean.held.out.marg.ll <- mean(
                (samples.l$logmarglike
                %>% left_join(gene.map)
                %>% filter(is.held.out))$logmarglike)
            message('Mean held out marginal log likelihood per cell: ',
                    mean.held.out.marg.ll / stan.data$C)
        }
        samples.l$tau <- (samples.l$tau
            %>% left_join(cell.map)
            %>% mutate(tau.offset=tau - obstime)
        )
        rm(melt.samples)
    })
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

#` Makes a distance periodic
#'
#' @param r Distance
#' @param period The period
#'
#' @export
#'
cov.periodise <- function(r, period) {
    period * sin(r * pi / period) / 2;
}

#` Matern 3/2 covariance function
#'
#' @param r Distance
#' @param l Length scale
#'
#' @export
#'
cov.matern.32 <- function(r, l) {
    x <- sqrt(3) * abs(r / l);
    (1 + x) * exp(- x);
}

#' Calculate distances between vectors of time points
#'
#' @param tau.1 First vector of time points
#' @param tau.2 Second vector of time points (defaults to first if not given)
#' @param period Period if periodic
#'
#' @export
#'
cov.calc.dists <- function(tau.1, tau.2=NULL, period=NULL) {
    if (is.null(tau.2)) {
        tau.2 <- tau.1
    }
    # Calculate the distances
    r <- outer(tau.1, tau.2, FUN="-")
    # Make periodic if necessary
    if (! is.null(period)) {
        r <- cov.periodise(r, period)
    }
    r
}

#' Calculate distances over estimated pseudotimes and
#' test inputs.
#'
#' @param dl de.lorean object
#' @param sample.iter Iteration to use (defaults to best.sample)
#'
#' @export
#'
cov.calc.dl.dists <- function(dl, sample.iter=NULL) {
    with(dl, {
        if (is.null(sample.iter)) {
            sample.iter <- best.sample
        }
        tau <- (
            samples.l$tau
            %>% filter(sample.iter == iter)  # Filter correct iteration
            %>% arrange(c)  # Sort by cell
        )$tau
        timepoints <- c(tau, test.input)
        if (opts$periodic) {
            cov.calc.dists(timepoints, opts$period)
        } else {
            cov.calc.dists(timepoints)
        }
    })
}

#' Calculate covariance structure for gene over pseudotimes and test
#' inputs.
#'
#' @param dl de.lorean object
#' @param gene.idx Gene index
#' @param cov.fn Covariance function (defaults to cov.matern.32)
#' @param sample.iter Iteration to use (defaults to best.sample)
#'
#' @export
#'
cov.calc.gene <- function(dl, gene.idx, cov.fn=NULL, sample.iter=NULL) {
    if (is.null(cov.fn)) {
        cov.fn <- cov.matern.32
    }
    with(dl, {
        if (is.null(sample.iter)) {
            sample.iter <- best.sample
        }
        r <- cov.calc.dl.dists(dl, sample.iter)
        psi   <- sampled.gene.param(dl, gene.idx, "psi"  , sample.iter)
        omega <- sampled.gene.param(dl, gene.idx, "omega", sample.iter)
        (
            psi * cov.fn(r, length.scale)  # Structure
            + omega * identity.matrix(nrow(r))  # Noise
        )
    })
}

#' Calculate covariance for gene over test inputs when conditioned on
#' data at estimated pseudotimes.
#'
#' @param dl de.lorean object
#' @param gene.idx Gene index
#' @param cov.fn Covariance function (defaults to cov.matern.32)
#' @param sample.iter Iteration to use (defaults to best.sample)
#'
#' @export
#'
cov.calc.gene.conditioned <- function(dl,
                                      gene.idx,
                                      cov.fn=NULL,
                                      sample.iter=NULL)
{
    if (is.null(cov.fn)) {
        cov.fn <- cov.matern.32
    }
    with(dl, {
        if (is.null(sample.iter)) {
            sample.iter <- best.sample
        }
        Sigma <- cov.calc.gene(dl, gene.idx, cov.fn, sample.iter)
        num.cells <- nrow(dl$cell.map)
        # From Appendix A.2 of Rasmussen/Williams GP book
        slice.obsd <- 1:num.cells
        slice.test <- (num.cells+1):nrow(Sigma)
        .B <- Sigma[slice.obsd,slice.obsd]
        .A <- Sigma[slice.test,slice.test]
        .C <- Sigma[slice.test,slice.obsd]
        .A - .C %*% qr.solve(.B, t(.C))
    })
}

#' Calculate covariances for all genes when conditioned on data at
#' estimated pseudotimes.
#'
#' @param dl de.lorean object
#' @param cov.fn Covariance function (defaults to cov.matern.32)
#' @param sample.iter Iteration to use (defaults to best.sample)
#'
#' @export
#'
cov.all.genes.conditioned <- function(dl,
                                      cov.fn=NULL,
                                      sample.iter=NULL)
{
    vapply(1:nrow(dl$gene.map),
           function(gene.idx) cov.calc.gene.conditioned(dl,
                                                        gene.idx,
                                                        cov.fn,
                                                        sample.iter),
           FUN.VALUE=matrix(0,
                            nrow=length(dl$test.input),
                            ncol=length(dl$test.input)))
}

#' Get the sampled parameter for the gene
#'
#' @param dl de.lorean object
#' @param gene.idx Gene index
#' @param param Parameter
#' @param sample.iter Iteration to use (defaults to best.sample)
#'
sampled.gene.param <- function(dl, gene.idx, param, sample.iter=NULL) {
    with(dl, {
        if (is.null(sample.iter)) {
            sample.iter <- best.sample
        }
        filter(samples.l[[param]], gene.idx == g, sample.iter == iter)[[param]]
    })
}

#' Identity matrix
#'
#' @param N size of matrix
#'
identity.matrix <- function(N) {
    result <- matrix(0, nrow=N, ncol=N)
    diag(result) <- 1
    result
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
                            %>% left_join(phi)
                            %>% mutate(tau=test.input[t]))
    })
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
        gp <- (ggplot(samples.l$S %>% left_join(cell.map),
                      aes(x=factor(cell,
                                   levels=arrange(cell.map, capture)$cell),
                          y=S,
                          color=capture),
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
#' @param sample.iter Which sample to plot
#'
#' @export
#'
plot.cmp.profiles <- function(...,
                              sample.iter = NULL,
                              genes = NULL) {
    if (is.null(sample.iter)) {
        sample.iter <- dl$best.sample
    }
    dls <- list(...)
    dl.levels <- names(dls)
    if (is.null(genes)) {
        genes <- dls[[1]]$genes.high.psi
    }
    stopifnot(! is.null(names(dls)))  # Must have names for de.lorean objects
    get.mean <- function(.name) {
        with(dls[[.name]], (
            predictions
            %>% filter(sample.iter == iter)
            %>% left_join(gene.map)
            %>% filter(gene %in% genes)
            %>% mutate(name=factor(.name, levels=dl.levels))
        ))
    }
    means <- do.call(rbind, lapply(names(dls), get.mean))
    (ggplot(means,
                  aes(x=tau),
                  environment=environment())
        + geom_line(alpha=.8,
                    aes(y=predictedmean + phi,
                        color=name))
        + geom_ribbon(aes(y=predictedmean + phi,
                          ymin=predictedmean+phi-2*sqrt(predictedvar),
                          ymax=predictedmean+phi+2*sqrt(predictedvar),
                          fill=name),
                      alpha=.2)
        + facet_wrap(~ gene)
        + scale_x_continuous(name="Pseudotime",
                             breaks=unique(dls[[1]]$obstime))
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
                          genes=NULL,
                          profile.color='black',
                          add.data=T,
                          sample.iter=NULL,
                          ...) {
    if (is.null(genes)) {
        genes <- dl$genes.high.psi
    }
    if (is.null(sample.iter)) {
        sample.iter <- dl$best.sample
    }
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
        gp <- (
            plot.add.profiles(gp, dl, color=profile.color, genes=genes, ...)
            + facet_wrap(~ gene)
            + scale_x_continuous(name="Pseudotime",
                                 breaks=unique(cell.meta$obstime))
            + scale_y_continuous(name="Expression")
        )
        if (add.data) {
            gp <- plot.add.expr(gp, dl, genes, ...)
        }
        gp
    })
}

#' Add expression profiles to a plot
#'
#' @param gp Plot object
#' @param dl de.lorean object
#' @param genes Genes to plot
#' @param color Color to use
#' @param sample.iter Which sample to plot
#' @param adjust.model A fitted model to adjust expression levels by
#'
#' @export
#'
plot.add.profiles <- function(gp,
                              dl,
                              color='black',
                              genes=NULL,
                              sample.iter=NULL,
                              adjust.model=NULL,
                              ...) {
    if (is.null(genes)) {
        genes <- dl$genes.high.psi
    }
    if (is.null(sample.iter)) {
        sample.iter <- dl$best.sample
    }
    .data <- (
        dl$predictions
        %>% filter(sample.iter == iter)
        %>% left_join(dl$gene.map)
        %>% filter(gene %in% genes)
    )
    if (! is.null(adjust.model)) {
        print(names(.data))
        adjustments <- (
            .data
            %>% group_by(t)
            %>% summarise(tau=tau[1]))
        print(tail(adjustments))
        adjustments$adjustment <- predict(adjust.model,
                                          newdata=adjustments)
        .T <- nrow(adjustments)
        print(adjustments$adjustment[1:(.T-1)]-adjustments$adjustment[2:.T])
        .data <- (
            .data
            %>% left_join(dplyr::select(adjustments, -tau))
            %>% mutate(predictedmean=predictedmean+adjustment))
    }
    (gp
        + geom_line(alpha=.3,
                    color=color,
                    data=.data,
                    aes(x=tau, y=predictedmean + phi))
        + geom_ribbon(aes(x=tau,
                          ymin=predictedmean+phi-2*sqrt(predictedvar),
                          ymax=predictedmean+phi+2*sqrt(predictedvar)),
                      data=.data,
                      fill=color,
                      alpha=.1))
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
plot.add.expr <- function(gp,
                          dl,
                          genes,
                          sample.iter=NULL,
                          cell.size.adj=T,
                          ...)
{
    if (is.null(sample.iter)) {
        sample.iter <- dl$best.sample
    }
    with(dl, {
        .data <- (
            samples.l$tau
            %>% filter(sample.iter == iter)
            %>% left_join(melt(unname(stan.m),
                          varnames=c("g", "c"),
                          value.name="expr"))
            %>% left_join(gene.map)
            %>% filter(gene %in% genes))
        if (cell.size.adj) {
            .data <- (
                .data
                %>% left_join(samples.l$S)
                %>% mutate(expr=expr - S))
        }
        gp + geom_point(data=.data,
                        aes(x=modulo.period(tau),
                            y=expr,
                            color=capture),
                        size=4,
                        alpha=.7)
    })
}

#' Single cell expression data and meta data from Trapnell et al. (2014).
#'
#' @docType data
#' @keywords datasets
#'
#' @name trapnell.expr
#' @aliases trapnell.cell.meta
#' @aliases trapnell.gene.meta
#'
#' @usage data(TrapnellDeLorean)
#'
#' @format There are three objects in this data:
#' \itemize{
#'   \item trapnell.expr A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item trapnell.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item trapnell.cell.meta A data frame containing meta-data
#'     about the cells
#' }
#'
#' @source \url{http://www.nature.com/nbt/journal/v32/n4/full/nbt.2859.html}
#'
NULL

#' Single cell expression data and meta data from McDavid et al. (2014).
#' They investigated differential expression in actively
#' cycling cells: "expression of 333 genes was interrogated in 930
#' cells, across three cell lines: H9 (HTB-176), MDA-MB-231 (HTB-26),
#' and PC3 (CRL-1435)".
#'
#' @docType data
#' @keywords datasets
#'
#' @name mcdavid.expr
#' @aliases mcdavid.cell.meta
#' @aliases mcdavid.gene.meta
#'
#' @usage data(McDavidDeLorean)
#'
#' @format There are three objects in this data:
#' \itemize{
#'   \item mcdavid.expr A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item mcdavid.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item mcdavid.cell.meta A data frame containing meta-data
#'     about the cells
#' }
#'
#' @source \url{http://www.ploscompbiol.org/article/info\%3Adoi\%2F10.1371\%2Fjournal.pcbi.1003696}
#'
NULL

#' Single cell expression data and meta data from Guo et al. (2012).
#' They investigated the expression of 48 genes in 500 mouse embryonic cells.
#'
#' @docType data
#' @keywords datasets
#'
#' @name guo.expr
#' @aliases guo.cell.meta
#' @aliases guo.gene.meta
#'
#' @usage data(GuoDeLorean)
#'
#' @format There are three objects in this data:
#' \itemize{
#'   \item guo.expr A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item guo.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item guo.cell.meta A data frame containing meta-data
#'     about the cells
#' }
#'
#' @source \url{http://www.sciencedirect.com/science/article/pii/S1534580710001103}
#'
NULL

#' Kouno et al. investigated the transcriptional network controlling how
#' THP-1 human myeloid monocytic leukemia cells differentiate into
#' macrophages. They provide expression values for 45 genes in 960 single
#' cells captured across 8 distinct time points.
#'
#' @docType data
#' @keywords datasets
#'
#' @name kouno.expr
#' @aliases kouno.cell.meta
#' @aliases kouno.gene.meta
#'
#' @usage data(KounoDeLorean)
#'
#' @format There are three objects in this data:
#' \itemize{
#'   \item kouno.expr A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item kouno.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item kouno.cell.meta A data frame containing meta-data
#'     about the cells
#' }
#'
#' @source \url{http://genomebiology.com/2013/14/10/R118/abstract}
#'
NULL

#' Windram et al. investigated the defense response in Arabidopsis
#' thaliana to the necrotrophic fungal pathogen Botrytis cinerea.
#' They collected data at 24 time points in two conditions for
#' 30336 genes.
#'
#' @docType data
#' @keywords datasets
#'
#' @name windram.expr
#' @aliases windram.cell.meta
#' @aliases windram.gene.meta
#'
#' @usage data(WindramDeLorean)
#'
#' @format There are three objects in this data:
#' \itemize{
#'   \item windram.expr A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item windram.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item windram.cell.meta A data frame containing meta-data
#'     about the cells
#' }
#'
#' @source \url{http://www.plantcell.org/content/24/9/3530.long}
#'
NULL
