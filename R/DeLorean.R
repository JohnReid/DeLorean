
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
plot.de.lorean <- function(dl, type="best.predictions") {
    result <- switch(type,
        best.predictions=plot.best.predictions(dl),
        S.posteriors=plot.S.posteriors(dl)
    )
    if (is.null(result)) {
        error('Unknown plot type')
    }
    result
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
    report.path <- system.file("inst", "Rmd", sprintf("%s.Rmd", report.name),
                               package="DeLorean")
    stylesheet.path <- system.file("inst", "Rmd", "foghorn.css",
                                   package="DeLorean")
    within(dl, {
        knit2html(report.path,
                  # output=paste(output.dir,
                  #              sprintf('%s.html', report.name),
                  #              sep='/'),
                  envir=environment(),
                  stylesheet=stylesheet.path)
    })
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

#' Estimate hyperparameters for model using empirical Bayes
#'
#' @param sigma.tau Noise s.d. in temporal dimension
#' @param delta Proportion of within time variance to relabel as between time
#' @param min.sd Minimum s.d. used for drop-out effects (to avoid s.d. of 0 when no drop outs)
#'
#' @param dl de.lorean object
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
        )
        stopifnot(! is.na(gene.expr))
        # Examine the variation within genes and times
        gene.time.expr <- (expr.l
            %>% left_join(cell.meta)
            %>% group_by(gene, capture)
            %>% summarise(x.mean=mean(x.hat),
                          x.var=var(x.hat))
        )
        stopifnot(! is.na(gene.time.expr))
        # Decomposition of variance within and between time.
        gene.var <- (gene.time.expr
            %>% group_by(gene)
            %>% summarise(omega.bar=mean(x.var, na.rm=TRUE),
                          psi.bar=var(x.mean, na.rm=TRUE))
            %>% mutate(within.time.mislabelled = opts$delta * omega.bar,
                       omega.hat = omega.bar - within.time.mislabelled,
                       psi.hat = psi.bar + within.time.mislabelled)
        )
        stopifnot(! is.na(gene.var))
        if (is.null(length.scale)) {
            length.scale <- time.width
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
            l_pe=length.scale
        )
        rm(expr.l)  # No longer needed
    })
}


#' Filter genes and cells
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
        stan.m <- expr
        .G <- nrow(stan.m)
        .C <- ncol(stan.m)
        dimnames(stan.m)
        gene.map <- (data.frame(g=1:.G,
                                gene=factor(rownames(stan.m),
                                            levels=levels(gene.meta$gene)))
                    %>% left_join(gene.expr)
                    %>% left_join(gene.var)
                    %>% left_join(gene.meta))
        print(names(gene.map))
        stopifnot(! is.na(gene.map[
            c("g", "gene", "x.mean", "x.sd",
            "phi.hat", "x.hat.sd", "omega.bar", "psi.bar",
            "within.time.mislabelled", "omega.hat", "psi.hat")
        ]))
        cell.map <- (data.frame(c=1:.C,
                                cell=factor(colnames(stan.m),
                                            levels=levels(cell.meta$cell)))
                    %>% left_join(cell.meta)
                    %>% left_join(cell.expr))
        stopifnot(! is.na(cell.map))
        # Have one missing value per gene
        get.missing.for.gene <- function(g) {
            sample(which(! is.na(stan.m[g,])), 1)
        }
        missing.idx <- sapply(1:.G, get.missing.for.gene)
        test.input <- (
            time.range[1] - 2 * opts$sigma.tau
            + (time.width + 4 * opts$sigma.tau)
                * (0:(opts$num.test-1)) / (opts$num.test-1))
        stan.data <- c(
            # Hyper-parameters
            hyper,
            list(
                # Dimensions
                C=ncol(stan.m),
                G=nrow(stan.m),
                # Data
                time=cell.map$obstime,
                expr=stan.m,
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
    within(dl, {
        stan.code <- '
        functions {
            #
            # Periodic function
            #
            real
            periodise(real r, real period) {
                return period * sin(r * pi() / period) / 2;
            }
            #
            # Squared exponential covariance function
            #
            real
            se_cov(real r, real l) {
                return exp(- pow(r / l, 2) / 2);
            }
            #
            # Matern nu=3/2 covariance function
            #
            real
            matern32_cov(real r, real l) {
                real x;
                x <- sqrt(3) * fabs(r / l);
                return (1 + x) * exp(- x);
            }
            #
            # Matern nu=5/2 covariance function
            #
            real
            matern52_cov(real r, real l) {
                real x;
                x <- sqrt(5) * fabs(r / l);
                return (1 + x + pow(x, 2) / 3) * exp(- x);
            }
            #
            # Covariance function
            #
            real
            cov_fn(real r, int periodic, real period, real l) {
                real rprime;
                if (periodic) {
                    rprime <- periodise(r, period);
                } else {
                    rprime <- r;
                }
                return matern32_cov(rprime, l);
                // return se_cov(rprime, l);
                // return se_cov(rprime, l);
                // return matern52_cov(rprime, l);
                // return matern52_cov(rprime, l);
            }
            #
            # Calculate symmetric covariance matrix
            #
            matrix
            cov_symmetric(row_vector tau, int periodic, real period, real l) {
                matrix[cols(tau),cols(tau)] result;
                for (c1 in 1:cols(tau)) {
                    for (c2 in 1:c1) {
                        result[c1,c2] <- cov_fn(tau[c2] - tau[c1], periodic, period, l);
                        if(c1 != c2) {
                            result[c2,c1] <- result[c1,c2];
                        }
                    }
                }
                return result;
            }
            #
            # Calculate covariance matrix
            #
            matrix
            cov(row_vector tau1,
                row_vector tau2,
                int periodic,
                real period,
                real l
            ) {
                matrix[cols(tau1),cols(tau2)] result;
                for (c1 in 1:cols(tau1)) {
                    for (c2 in 1:cols(tau2)) {
                        result[c1,c2] <- cov_fn(tau2[c2] - tau1[c1], periodic, period, l);
                    }
                }
                return result;
            }
            void
            pretty_print_tri_lower(matrix x) {
                if (rows(x) == 0) {
                    print("empty matrix");
                    return;
                }
                print("rows=", rows(x), " cols=", cols(x));
                for (m in 1:rows(x))
                    for (n in 1:m)
                        print("[", m, ",", n, "]=", x[m,n]);
            }
        }
        data {
            #
            # Dimensions
            #
            int<lower=2> C;  // Number of cells
            int<lower=2> G;  // Number of genes
            #
            # Data
            #
            # Time data
            int periodic; // Are the expression patterns periodic?
            real period;  // Cyclic period
            row_vector<lower=0>[C] time;  // Time index for cell c
            # Expression data
            vector[C] expr[G];
            #
            # Hyperparameters
            #
            real mu_S;  // Mean of cell size factor, S
            real<lower=0> sigma_S;  // S.d. of cell size factor, S
            real mu_phi;  // Mean of gene mean, phi
            real<lower=0> sigma_phi;  // S.d. of gene mean, phi
            real mu_psi;  // Mean of log between time variation, psi
            real<lower=0> sigma_psi;  // S.d. of log between time variation, psi
            real mu_omega;  // Mean of log within time variation, omega
            real<lower=0> sigma_omega;  // S.d. of log within time variation, omega
            real l_pe;  // Length scale squared for phi
            real<lower=0> sigma_tau;  // Standard deviation for pseudotime
            #
            # Generated quantities
            #
            # Test inputs for predicted mean
            int numtest;
            row_vector[numtest] testinput;
        }
        transformed data {
            //
            // Unit diagonal covariance matrix
            cov_matrix[C] identity;
            //
            // Transformations of expression
            //
            // Unit diagonal covariance matrix
            identity <- diag_matrix(rep_vector(1, C));
        }
        parameters {
            row_vector[C] S;      // Cell-size factor for expression
            row_vector[C] tau;    // Pseudotime
            row_vector[G] phi;    // Mean positive expression for each gene
            row_vector<lower=0>[G] psi;    // Between time PE variance
            row_vector<lower=0>[G] omega;  // Within time PE variance
        }
        model {
            //
            // Sample cell-specific factors
            S ~ normal(mu_S, sigma_S);  // Cell size factors
            //
            // Sample gene-specific factors
            phi ~ normal(mu_phi, sigma_phi);
            psi ~ lognormal(mu_psi, sigma_psi);
            omega ~ lognormal(mu_omega, sigma_omega);
            //
            // Sample pseudotime
            tau ~ normal(time, sigma_tau);  // Pseudotime
            //
            // For each gene
            for (g in 1:G) {
                expr[g] ~ multi_normal(
                              S + phi[g],
                              psi[g] * cov_symmetric(tau,
                                                     periodic,
                                                     period,
                                                     l_pe)
                                  + omega[g] * identity);
            }
        }
        generated quantities {
            row_vector[numtest] predictedmean[G];
            vector[numtest] predictedvar[G];
            #
            # For each gene
            for (g in 1:G) {
                matrix[C,C] L_g;
                row_vector[C] a;
                vector[C] v;
                matrix[C,numtest] kstartest;
                matrix[C,numtest] vtest;
                //
                // Cholesky decompose the covariance of the inputs
                L_g <- cholesky_decompose(
                            psi[g] * cov_symmetric(tau,
                                                   periodic,
                                                   period,
                                                   l_pe)
                            + omega[g] * identity);
                a <- mdivide_right_tri_low(
                        mdivide_left_tri_low(
                            L_g,
                            expr[g] - S\' - phi[g])\',
                        L_g);
                //
                // Calculate predicted mean on test inputs
                kstartest <- psi[g] * cov(tau,
                                          testinput,
                                          periodic,
                                          period,
                                          l_pe);
                predictedmean[g] <- a * kstartest;
                //
                // Calculate predicted variance on test inputs
                vtest <- mdivide_left_tri_low(L_g, kstartest);
                predictedvar[g] <- (psi[g]
                                    + omega[g]
                                    - diagonal(vtest\' * vtest));
            }
        }
        '
        compiled <- stan(model_code=stan.code, chains=0)
        # Define a function to initialise the chains
        init.chain <- function() {
            with(stan.data, {
                list(S=cell.map$S.hat,
                     tau=rnorm(C, mean=time, sd=sigma_tau),
                     phi=rnorm(G, mean=mu_phi, sd=sigma_phi),
                     psi=rlnorm(G, meanlog=mu_psi, sdlog=sigma_psi),
                     omega=rlnorm(G, meanlog=mu_omega, sdlog=sigma_omega)
                )
            })
        }
        # Try one iteration to check everything is OK
        fit <- stan(fit=compiled,
                    data=stan.data,
                    init=init.chain,
                    iter=1,
                    chains=1)
    })
}


#' Find best tau
#'
#' @param dl de.lorean object
#'
#' @export
#'
find.best.tau <- function(dl, num.tau.candidates = 6000) {
    within(dl, {
        # Define a function that chooses tau
        init.chain.find.tau <- function() {
            with(stan.data, {
                list(alpha=cell.map$alpha.hat,
                     beta=rep(0, G),
                     S=cell.map$S.hat,
                     tau=rnorm(C, time, sd=sigma_tau),
                     phi=gene.map$phi.hat,
                     psi=gene.map$psi.hat,
                     omega=gene.map$omega.hat
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
        tau.inits <- lapply(1:num.tau.candidates, try.tau.init)
        # qplot(sapply(tau.inits, function(init) init$lp))
        # Which tau gave highest log probability?
        tau.inits.order <- order(sapply(tau.inits, function(init) -init$lp))
        sapply(tau.inits[tau.inits.order], function(init) init$lp)[1:10]
        tau.inits[[tau.inits.order[1]]]
    })
}


#' Fit the model
#'
#' @param num.cores Number of cores to run on. Defaults to max(detectCores()-1, 1)
#' @param chain Number of chains to run on each core
#' @param iter Number of iterations in each chain
#' @param thin How many samples to generate before retaining one
#' @param dl de.lorean object
#'
#' @export
#'
fit.model <- function(
    dl,
    num.cores = NULL,
    chains = 1,
    iter = 1000,
    thin = 50)
{
    if (is.null(num.cores)) {
        num.cores <- max(detectCores() - 1, 1)
    }
    init.chain.good.tau <- function(chain_id) {
        # print(chain_id)
        #
        # Create random parameters
        pars <- dl$init.chain()
        #
        # Replace tau with good tau
        pars$tau <- dl$tau.inits[[dl$tau.inits.order[chain_id]]]$tau
        pars
    }
    # Run the chains in parallel
    sflist <- mclapply(1:num.cores,
                       mc.cores=num.cores,
                       function(i)
                           stan(fit=dl$compiled, data=dl$stan.data,
                                thin=thin,
                                init=init.chain.good.tau,
                                # init=init.chain,
                                iter=iter,
                                seed=i, chains=chains,
                                chain_id=i, refresh=-1))
    dl$fit <- sflist2stanfit(sflist)
    dl$la <- extract(dl$fit, permuted=TRUE)
    # mean(as.vector(la$ll))
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
            # ll=c("g"),
            predictedmean=c("g", "t"),
            predictedvar=c("g", "t")
        )
        samples.l <- melt.samples(la, sample.dims)
        best.sample <- which.max(samples.l$lp__$lp__)
        # print(mean(as.vector(la$ll)))
        # print(mean(la$ll[best.sample,]))
        samples.l$tau <- (samples.l$tau
            %>% left_join(cell.map)
            %>% mutate(tau.offset=tau - obstime)
        )
        samples.all <- (
            Reduce(left_join, samples.l[! str_detect(names(samples.l),
                                                     "^predicted")])
            %>% left_join(melt(unname(stan.m),
                            varnames=c("g", "c"),
                            value.name="expr"))
            %>% left_join(gene.map)
        )
        # Just the sample with the best log probability
        sample.best <- samples.all %>% filter(best.sample == iter)
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
        best.mean <- filter(predictions, best.sample == iter)
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
        gp <- (ggplot(samples.all,
                      aes(x=factor(cell,
                                   levels=arrange(cell.map, capture)$cell),
                          y=S,
                          color=capture),
                      environment=environment())
            + geom_boxplot()
        )
    })
}

#' Plot best sample predicted expression
#'
#' @param dl de.lorean object
#' @param genes Genes to plot (defaults to genes.high.psi)
#'
#' @export
#'
plot.best.predictions <- function(dl, genes=NULL) {
    if (is.null(genes)) {
        genes <- dl$genes.high.psi
    }
    with(dl, {
        if (opts$periodic) {
            modulo.period <- function(t) ( t - floor(t / opts$period)
                                                * opts$period )
        } else {
            modulo.period <- function(t) { t }
        }
        gp <- (ggplot(best.mean
                          %>% left_join(gene.map)
                          %>% filter(gene %in% genes),
                      aes(x=modulo.period(tau), y=predictedmean + phi),
                      environment=environment())
            + geom_line(alpha=.3)
            + geom_ribbon(aes(x=modulo.period(tau),
                              y=predictedmean + phi,
                              ymin=predictedmean+phi-2*sqrt(predictedvar),
                              ymax=predictedmean+phi+2*sqrt(predictedvar)),
                        alpha=.1)
            + geom_point(aes(x=modulo.period(tau),
                             y=expr - S,
                             color=capture),
                        data=sample.best %>% filter(gene %in% genes),
                        size=4,
                        alpha=.7)
            + facet_wrap(~ gene)
            + scale_x_continuous(name="Pseudotime",
                                 breaks=unique(cell.meta$obstime))
            + scale_y_continuous(name="Expression (log10)")
        )
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
