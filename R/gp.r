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

#' Identity matrix
#'
#' @param N size of matrix
#'
identity.matrix <- function(N) {
    result <- matrix(0, nrow=N, ncol=N)
    diag(result) <- 1
    result
}


#' Add posterior representation to a plot.
#'
#' @param gp Plot object
#' @param .data Data frame containing variables to plot (mean, var)
#'        phi, predictedvar)
#' @param color Color to use
#' @param line.alpha Alpha to use for mean line
#' @param ribbon.alpha Alpha to use for variance ribbon
#'
#' @export
#'
plot.add.mean.and.variance <- function(gp,
                                       .data=NULL,
                                       color='black',
                                       line.alpha=.3,
                                       ribbon.alpha=.1) {
    (gp
        + geom_line(data=.data,
                    aes(x=x, y=mean),
                    color=color,
                    alpha=line.alpha)
        + geom_ribbon(data=.data,
                      aes(x=x,
                          ymin=mean-2*sqrt(var),
                          ymax=mean+2*sqrt(var)),
                      fill=color,
                      alpha=ribbon.alpha))
}

#' The log marginal likelihood. See "2.3 Varying the Hyperparameters"
#' on page 19 of Rasmumssen and Williams' book.
#'
#' @param y The targets.
#' @param K The covariance matrix (kernel), not needed if U is provided.
#' @param U Cholesky decomposition of K (chol(K)).
#'
#' @export
#'
gp.log.marg.like <- function(y, K=NULL, U=chol(K)) {
    alpha <- backsolve(U, backsolve(U, y, transpose = TRUE))
    -(
        (t(y) %*% alpha) / 2
        + sum(diag(U))
        + length(y) * log(2 * pi)
    )
}
