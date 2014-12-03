
#' Initialise DeLorean object
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
        cell.meta = cell.meta)
    class(result) <- c("de.lorean", class(result))
    result
}

#' Is a DeLorean object?
#'
is.de.lorean <- function(dl) inherits(dl, "de.lorean")

#' Print details of DeLorean object
#'
print.de.lorean <- function(dl) {
    print(sapply(dl, head))
}

#' Dimensions of DeLorean object
#'
dim.de.lorean <- function(dl) {
    dim(dl$expr)
}

#' Summarise DeLorean object
#'
summary.de.lorean <- function(dl) {
    print(sapply(dl, summary))
}

#' Various DeLorean object plots
#'
plot.de.lorean <- function(dl, type="best.predictions") {
    result <- switch(type,
        best.predictions=plot.best.predictions(dl)
    )
    if (is.null(result)) {
        error('Unknown plot type')
    }
}

#' Knit a report, the file inst/Rmd/<report.name>.Rmd must exist in
#' the package directory.
#'
knit.report <- function(dl, report.name) {
    report.path <- system.file("inst", "Rmd", sprintf("%s.Rmd", report.name),
                               package="DeLorean")
    stylesheet.path <- system.file("inst", "Rmd", "foghorn.css",
                                   package="DeLorean")
    with(dl, {
        knit2html(report.path,
                  # envir=globalenv(),
                  stylesheet=stylesheet.path)
    })
}

#' Calculate a suitable value for a rug plot given the
#' number of points
alpha.for.rug <- function(n, scale=100) {
    1 / (max(1, n / scale))
}

#' Estimate hyperparameters for model using empirical Bayes
#'
#' @param sigma.tau Noise s.d. in temporal dimension
#' @param delta Proportion of within time variance to relabel as between time
#' @param min.sd Minimum s.d. used for drop-out effects (to avoid s.d. of 0 when no drop outs)
#'
estimate.hyper <- function(
    dl,
    sigma.tau = .5,
    delta = .5,
    min.sd = 1e-10,
    length.scale = NULL
) {
    result <- within(dl, {
        #
        # Set up temporal hyper-parameters
        #
        opts <- list(
            delta = delta,
            sigma.tau = sigma.tau
        )
        time.range <- range(cell.meta$obstime)
        time.width <- time.range[2] - time.range[1]
        #
        # First melt expression data into long format
        #
        expr.l <- melt(expr, varnames=c("gene", "cell"), value.name="x")
        expr.l$gene <- factor(expr.l$gene, levels=levels(gene.meta$gene))
        expr.l$cell <- factor(expr.l$cell, levels=levels(cell.meta$cell))
        #
        # Estimate drop out rates: theta
        #
        # Define a function to estimate theta with 1 pseudo-count for each
        # outcome (ON/OFF)
        estimate.theta <- function(x) {
            (sum(! is.na(x)) + 1) / (length(x) + 2)
        }
        # Cell-specific thetas
        expr.cell <- (expr.l
            %>% group_by(cell)
            %>% summarise(theta=estimate.theta(x))
            %>% mutate(alpha.hat=log(theta/(1-theta)))
        ) %>% left_join(cell.meta)
        # Gene-specific thetas
        expr.gene <- (expr.l
            %>% group_by(gene)
            %>% summarise(theta=estimate.theta(x))
            %>% mutate(beta.hat=log(theta/(1-theta)))
        )
        #
        # Look at positive expression
        #
        expr.pos <- expr.l %>% filter(! is.na(x)) %>% left_join(cell.meta)
        # Estimate a pseudo reference mean for each gene
        gene.pos <- (expr.pos
            %>% group_by(gene)
            %>% summarise(x.mean=mean(x))
        )
        # Estimate the cell size by the median of the expression
        # adjust by the gene's mean
        cell.pos <- (expr.pos
            %>% left_join(gene.pos)
            %>% group_by(cell)
            %>% summarise(S.hat=median(x - x.mean))
        )
        stopifnot(! is.na(cell.pos))
        # Adjust the positive expression by the cell size estimates
        expr.pos <- (expr.pos
            %>% left_join(cell.pos)
            %>% mutate(x.hat=x - S.hat)
        )
        stopifnot(! is.na(expr.pos))
        # Resummarise the adjusted positive expression data
        gene.pos <- (expr.pos
            %>% group_by(gene)
            %>% summarise(x.mean=mean(x),
                          x.sd=sd(x),
                          phi.hat=mean(x.hat),
                          x.hat.sd=sd(x.hat),
                          num.pos=length(x))
            %>% filter(num.pos > 1)
        )
        stopifnot(! is.na(gene.pos))
        # Examine the variation
        gene.time.pos <- (expr.pos
            %>% group_by(gene, capture)
            %>% summarise(x.mean=mean(x.hat),
                          x.var=var(x.hat),
                          num.pos=length(x.hat))
            %>% filter(num.pos > 1)
        )
        # Decomposition of variance within and between time.
        gene.var <- (gene.time.pos
            %>% group_by(gene)
            %>% summarise(omega.bar=mean(x.var, na.rm=TRUE),
                          psi.bar=var(x.mean, na.rm=TRUE),
                          num.pos=length(x.var))
            %>% filter(num.pos > 1)
            %>% mutate(within.time.mislabelled = opts$delta * omega.bar,
                       omega.hat = omega.bar - within.time.mislabelled,
                       psi.hat = psi.bar + within.time.mislabelled)
        )
        stopifnot(! is.na(gene.var))
        if (is.null(length.scale)) {
            length.scale <- time.width / 2
        }
        hyper <- list(
            mu_S=mean(cell.pos$S.hat),
            sigma_S=sd(cell.pos$S.hat),
            mu_alpha=mean(expr.cell$alpha.hat),
            sigma_alpha=max(sd(expr.cell$alpha.hat), min.sd),
            sigma_beta=max(sd(expr.gene$beta.hat), min.sd),
            mu_phi=mean(gene.pos$phi.hat),
            sigma_phi=sd(gene.pos$phi.hat),
            mu_psi=mean(log(gene.var$psi.hat), na.rm=TRUE),
            sigma_psi=sd(log(gene.var$psi.hat), na.rm=TRUE),
            mu_omega=mean(log(gene.var$omega.hat), na.rm=TRUE),
            sigma_omega=sd(log(gene.var$omega.hat), na.rm=TRUE),
            sigma_tau=opts$sigma.tau,
            l_pe=length.scale
        )
    })
}


#' Filter genes and cells
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
format.for.stan <- function(
    dl,
    num.test = 101,  # Number of test points to consider
    period = 0  # Period of expression patterns
) {
    within(dl, {
        opts$num.test <- num.test
        opts$period <- period
        opts$periodic <- opts$period > 0
        min.expr <- min(expr, na.rm=TRUE)
        if( min.expr > 0 ) {
            stan.minexpr <- -2 * min.expr
        } else {
            stan.minexpr <-  2 * min.expr
        }
        stan.m <- expr
        .G <- nrow(stan.m)
        .C <- ncol(stan.m)
        dimnames(stan.m)
        gene.map <- (data.frame(g=1:.G,
                                gene=factor(rownames(stan.m),
                                            levels=levels(gene.meta$gene)))
                    %>% left_join(gene.pos %>% select(-num.pos))
                    %>% left_join(gene.var %>% select(-num.pos))
                    %>% left_join(expr.gene)
                    %>% left_join(gene.meta))
        stopifnot(! is.na(gene.map))
        cell.map <- (data.frame(c=1:.C,
                                cell=factor(colnames(stan.m),
                                            levels=levels(cell.meta$cell)))
                    %>% left_join(cell.meta)
                    %>% left_join(cell.pos)
                    %>% left_join(expr.cell))
        stopifnot(! is.na(cell.map))
        expr.l <- (expr.l
                %>% left_join(gene.map)
                %>% left_join(cell.map))
        # Have one missing value per gene
        get.missing.for.gene <- function(g) {
            sample(which(! is.na(stan.m[g,])), 1)
        }
        missing.idx <- sapply(1:.G, get.missing.for.gene)
        # Replace NA with large negative number
        stan.m[is.na(stan.m)] <- 2 * stan.minexpr
        test.input <- (
            time.range[1] - 2 * opts$sigma.tau
            + (time.width + 4 * opts$sigma.tau)
                * (0:(opts$num.test-1)) / (opts$num.test-1))
        expr.f.stan <- with(expr.l, {
            c(
                # Hyper-parameters
                hyper,
                list(
                    # Dimensions
                    C=ncol(stan.m),
                    G=nrow(stan.m),
                    # Data
                    time=cell.map$obstime,
                    expr=stan.m,
                    minexpr=stan.minexpr,
                    heldout=missing.idx,
                    # Generated quantities
                    numtest=opts$num.test,
                    testinput=test.input,
                    periodic=opts$periodic,
                    period=opts$period
                )
            )
        })
    })
}

#' Define and compile a simple model without dropout
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
            real minexpr;
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
                #
                # Evaluate the log likelihood of the heldout data
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
            with(expr.f.stan, {
                list(alpha=rnorm(C, mean=mu_alpha, sd=sigma_alpha),
                     beta=rnorm(G, sd=sigma_beta),
                     S=cell.map$S.hat,
                     tau=rnorm(C, mean=time, sd=sigma_tau),
                     phi=rnorm(G, mean=mu_phi, sd=sigma_phi),
                     psi=rlnorm(G, meanlog=mu_psi, sdlog=sigma_psi),
                     omega=rlnorm(G, meanlog=mu_omega, sdlog=sigma_omega)
                )
            })
        }
        # Try one iteration to check everything is OK
        fit <- stan(fit=compiled,
                    data=expr.f.stan,
                    init=init.chain,
                    iter=1,
                    chains=1)
    })
}

#' Define and compile the model
#'
compile.model <- function(dl) {
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
            #
            # Slice the indexed elements from a vector
            #
            row_vector
            slice_vector(row_vector vec, int[] indexes) {
                row_vector[size(indexes)] result;
                for (w1 in 1:size(indexes)) {
                    result[w1] <- vec[indexes[w1]];
                }
                return result;
            }
            #
            # Slice the indexed rows and columns from a symmetric matrix
            #
            matrix
            slice_symmetric_matrix(matrix mat, int[] indexes) {
                matrix[size(indexes),size(indexes)] result;
                for (w1 in 1:size(indexes)) {
                    int c1;
                    c1 <- indexes[w1];
                    for (w2 in 1:w1) {
                        int c2;
                        c2 <- indexes[w2];
                        result[w1,w2] <- mat[c1,c2];
                        if (c1 != c2) {
                            result[w2,w1] <- mat[c1,c2];
                        }
                    }
                }
                return result;
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
            real minexpr;
            # Held out data, one cell per gene
            int heldout[G];
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
            real mu_alpha;  // Mean of cell drop out rate, alpha
            real<lower=0> sigma_alpha;  // S.d. of alpha
            real<lower=0> sigma_beta;  // S.d. of beta (mean is 0)
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
            vector[G] heldoutexpr;  // Held out expression
            vector[C] geneexpr[G];  // Vector of expression
            int<lower=0,upper=1> isexpr[G,C]; // Is the gene expressed in the cell
            int<lower=1,upper=C> numexpr[G]; // How many cells is the gene expressed in
            int<lower=0,upper=C> whichexpr[G,C]; // The cells the gene is expressed in
            //
            // Unit diagonal covariance matrix
            identity <- diag_matrix(rep_vector(1, C));
            //
            // Process expr values
            for (g in 1:G) {
                int NE;  // Number expressed for this gene
                NE <- 0;
                for(c in 1:C) {
                    # If the gene is expressed
                    if (expr[g,c] >= minexpr) {
                        isexpr[g,c] <- 1;
                        # Is the data heldout?
                        if (heldout[g] == c) {
                            heldoutexpr[g] <- expr[g,c];
                        } else {
                            NE <- NE + 1;
                            geneexpr[g,NE] <- expr[g,c];
                            whichexpr[g,NE] <- c;
                        }
                    } else {
                        isexpr[g,c] <- 0;
                    }
                }
                numexpr[g] <- NE;
            }
        }
        parameters {
            row_vector[C] alpha;  // Dropout cell-specific coefficient
            row_vector[G] beta;   // Dropout gene-specific coefficient
            row_vector[C] S;      // Cell-size factor for expression
            row_vector[C] tau;    // Pseudotime
            row_vector[G] phi;    // Mean positive expression for each gene
            row_vector<lower=0>[G] psi;    // Between time PE variance
            row_vector<lower=0>[G] omega;  // Within time PE variance
        }
        transformed parameters {
            //
            // Covariance matrix
            # cov_matrix[C] Sigma_pe;  // PE covariance function
            matrix[C,C] Sigma_pe;  // PE covariance function
            //
            // Calculate covariances
            Sigma_pe <- cov_symmetric(tau, periodic, period, l_pe);
        }
        model {
            //
            // Sample cell-specific factors
            S ~ normal(mu_S, sigma_S);  // Cell size factors
            alpha ~ normal(mu_alpha, sigma_alpha);  // Cell-specific dropout effect
            //
            // Sample gene-specific factors
            beta ~ normal(0, sigma_beta);
            phi ~ normal(mu_phi, sigma_phi);
            psi ~ lognormal(mu_psi, sigma_psi);
            omega ~ lognormal(mu_omega, sigma_omega);
            //
            // Sample pseudotime
            tau ~ normal(time, sigma_tau);  // Pseudotime
            //
            // For each gene
            for (g in 1:G) {
                matrix[numexpr[g],numexpr[g]] Sigma_pe_g;
                row_vector[numexpr[g]] S_g;
                //
                // Extract covariance structure and mean for expressed cells
                // from overall covariance structure and mean and add
                // noise covariance
                Sigma_pe_g <- (psi[g] * slice_symmetric_matrix(
                                            Sigma_pe,
                                            head(whichexpr[g], numexpr[g]))
                                + omega[g] * diag_matrix(rep_vector(1, numexpr[g])));
                S_g <- slice_vector(S, head(whichexpr[g], numexpr[g]));
                //
                // Sample whether genes are expressed
                isexpr[g] ~ bernoulli_logit(alpha + beta[g]);
                //
                // Sample expressed cells
                head(geneexpr[g], numexpr[g]) ~ multi_normal(S_g + phi[g], Sigma_pe_g);
            }
        }
        generated quantities {
            row_vector[G] ll;
            vector[G] heldoutmu;
            vector[G] heldoutvar;
            row_vector[numtest] predictedmean[G];
            vector[numtest] predictedvar[G];
            #
            # For each gene
            for (g in 1:G) {
                int NE;
                NE <- numexpr[g];
                {
                    #
                    # Evaluate the log likelihood of the heldout data
                    int whichslice[NE];
                    matrix[NE,NE] Sigma_pe_g;
                    row_vector[NE] S_g;
                    row_vector[NE] kstar;
                    matrix[NE,NE] L;
                    row_vector[NE] a;
                    vector[NE] v;
                    matrix[NE,numtest] kstartest;
                    matrix[NE,numtest] vtest;
                    //
                    // Slice which cells are expressed array
                    whichslice <- head(whichexpr[g], NE);
                    //
                    // Extract covariance structure and mean for expressed cells
                    // from overall covariance structure and mean
                    Sigma_pe_g <- (psi[g] * slice_symmetric_matrix(
                                                Sigma_pe,
                                                whichslice)
                                    + omega[g] * diag_matrix(rep_vector(1, NE)));
                    S_g <- slice_vector(S, whichslice);
                    //
                    // Cholesky decompose the covariance of the inputs
                    L <- cholesky_decompose(Sigma_pe_g);
                    a <- mdivide_right_tri_low(
                            mdivide_left_tri_low(
                                L,
                                head(geneexpr[g], NE) - S_g\' - phi[g])\',
                            L);
                    //
                    // Calculate covariance between inputs and test inputs
                    kstar <- slice_vector(psi[g] * Sigma_pe[heldout[g]],
                                        whichslice);
                    v <- mdivide_left_tri_low(L, kstar\');
                    //
                    // Now calculate the mean and variance to use for the log likelihood
                    heldoutmu[g] <- S[heldout[g]] + phi[g] + dot_product(kstar, a);
                    # print(psi[g], ", ", dot_self(v));
                    heldoutvar[g] <- psi[g] + omega[g] - dot_self(v);
                    ll[g] <- normal_log(heldoutexpr[g],
                                        heldoutmu[g],
                                        sqrt(heldoutvar[g]));
                    # print(heldoutexpr[g], ", ", heldoutmu[g], ", ", sqrt(heldoutvar[g]), ", ", ll[g]);
                    //
                    // Calculate predicted mean on test inputs
                    kstartest <- psi[g] * cov(slice_vector(tau, whichslice),
                                            testinput,
                                            periodic,
                                            period,
                                            l_pe);
                    predictedmean[g] <- a * kstartest;
                    //
                    // Calculate predicted variance on test inputs
                    vtest <- mdivide_left_tri_low(L, kstartest);
                    predictedvar[g] <- (psi[g]
                                        + omega[g]
                                        - diagonal(vtest\' * vtest));
                }
            }
        }
        '
        compiled <- stan(model_code=stan.code, chains=0)
        # Define a function to initialise the chains
        init.chain <- function() {
            with(expr.f.stan, {
                list(alpha=rnorm(C, mean=mu_alpha, sd=sigma_alpha),
                     beta=rnorm(G, sd=sigma_beta),
                     S=cell.map$S.hat,
                     tau=rnorm(C, mean=time, sd=sigma_tau),
                     phi=rnorm(G, mean=mu_phi, sd=sigma_phi),
                     psi=rlnorm(G, meanlog=mu_psi, sdlog=sigma_psi),
                     omega=rlnorm(G, meanlog=mu_omega, sdlog=sigma_omega)
                )
            })
        }
        # Try one iteration to check everything is OK
        fit <- stan(fit=compiled,
                    data=expr.f.stan,
                    init=init.chain,
                    iter=1,
                    chains=1)
    })
}


#' Find best tau
#'
find.best.tau <- function(dl, num.tau.candidates = 6000) {
    within(dl, {
        # Define a function that chooses tau
        init.chain.find.tau <- function() {
            with(expr.f.stan, {
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
    within(dl, {
        init.chain.good.tau <- function(chain_id) {
            # print(chain_id)
            #
            # Create random parameters
            pars <- init.chain()
            #
            # Replace tau with good tau
            pars$tau <- tau.inits[[tau.inits.order[chain_id]]]$tau
            pars
        }
        # Run the chains in parallel
        sflist <- mclapply(1:num.cores,
                        mc.cores=num.cores,
                        function(i)
                            stan(fit=compiled, data=expr.f.stan, thin=thin,
                                    init=init.chain.good.tau,
                                    # init=init.chain,
                                    iter=iter,
                                    seed=i, chains=chains,
                                    chain_id=i, refresh=-1))
        fit <- sflist2stanfit(sflist)
        la <- extract(fit, permuted=TRUE)
        # mean(as.vector(la$ll))
    })
}


#' Examine convergence
#'
examine.convergence <- function(dl) {
    within(dl, {
        summ <- monitor(fit,
                        print=FALSE,
                        pars=c("tau", "psi", "S", "phi", "omega"))
        ignore.names <- str_detect(rownames(summ),
                                   "^(predictedvar|predictedmean")
        rhat.sorted <- sort(summ[! ignore.names, "Rhat"])
    })
}


#' Process the posterior
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
            # alpha=c("c"),
            S=c("c"),
            tau=c("c"),
            # beta=c("g"),
            phi=c("g"),
            psi=c("g"),
            omega=c("g"),
            # ll=c("g"),
            # heldoutmu=c("g"),
            # heldoutvar=c("g"),
            predictedmean=c("g", "t"),
            predictedvar=c("g", "t")
        )
        samples.l <- melt.samples(la, sample.dims)
        best.sample <- which.max(samples.l$lp__$lp__)
        best.sample
        # print(mean(as.vector(la$ll)))
        # print(mean(la$ll[best.sample,]))
        samples.l$tau <- (samples.l$tau
            %>% left_join(cell.map %>% select(-theta, -alpha.hat))
            %>% mutate(tau.offset=tau - obstime)
        )
        samples.all <- (
            Reduce(left_join, samples.l[! str_detect(names(samples.l),
                                                     "^predicted")])
            %>% left_join(melt(unname(stan.m),
                            varnames=c("g", "c"),
                            value.name="expr"))
            %>% left_join(gene.map %>% select(-theta, -beta.hat))
        )
        # Just the sample with the best log probability
        sample.best <- samples.all %>% filter(best.sample == iter)
    })
}


#' Analyse noise levels
#'
analyse.noise.levels <- function(dl, num.high.psi=16) {
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
            %>% arrange(omega/psi))
        genes.high.psi <- head(gene.noise.levels$gene, num.high.psi)
    })
}


#' Make predictions
#'
make.predictions <- function(dl) {
    within(dl, {
        predictions <- with(samples.l,
                            predictedmean
                            %>% left_join(predictedvar)
                            %>% left_join(S)
                            %>% mutate(tau=test.input[t]))
        best.mean <- (predictions
                    %>% filter(best.sample == iter)
                    %>% left_join(gene.map))
    })
}

#' Plot best sample predicted expression
#'
plot.best.predictions <- function(dl) {
    with(dl, {
        if (opts$periodic) {
            modulo.period <- function(t) ( t - floor(t / opts$period)
                                                * opts$period )
        } else {
            modulo.period <- function(t) { t }
        }
        gp <- (ggplot(best.mean %>% filter(gene %in% genes.high.psi),
                      aes(x=modulo.period(tau), y=predictedmean),
                      environment=environment())
            + geom_line(alpha=.3)
            + geom_ribbon(aes(x=modulo.period(tau),
                              y=predictedmean,
                              ymin=predictedmean-2*sqrt(predictedvar),
                              ymax=predictedmean+2*sqrt(predictedvar)),
                        alpha=.1)
            + geom_point(aes(x=modulo.period(tau),
                             y=expr - phi - S,
                             color=capture),
                        data=sample.best %>% filter(gene %in% genes.high.psi),
                        #data=sample.best %>% filter(gene %in% genes.high.psi,
                        #                            expr > stan.minexpr),
                        size=4,
                        alpha=.7)
            + facet_wrap(~ gene)
            + scale_x_continuous(name="Pseudotime",
                                 breaks=unique(cell.meta$obstime))
            + scale_y_continuous(name="Expression")
        )
    })
}

