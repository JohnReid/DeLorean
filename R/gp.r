#' Makes a distance periodic
#'
#' @param r Distance
#' @param period The period
#'
cov.periodise <- function(r, period) {
    period * sin(r * pi / period) / 2;
}

#' Matern 3/2 covariance function
#'
#' @param r Distance
#' @param l Length scale
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
cov.calc.dists <- function(tau.1, tau.2=tau.1, period=NULL) {
    # Calculate the distances
    r <- outer(tau.1, tau.2, FUN="-")
    # Make periodic if necessary
    if (! is.null(period) && period > 0) {
        r <- cov.periodise(r, period)
    }
    r
}


#' Calculate distances over estimated pseudotimes and
#' test inputs.
#'
#' @param dl de.lorean object
#' @param tau The pseudotimes to use
#' @param include.test Also include the pseudotimes for the test inputs
#'
cov.calc.dl.dists <- function(dl,
                              tau=tau.for.sample(dl),
                              include.test=TRUE) {
    with(dl, {
        if (length(tau) == 1 && "capture" == tau) {
            tau <- dl$cell.map$obstime
        }
        if (include.test) {
            tau <- c(tau, test.input)
        }
        if (exists('periodic', opts) && opts$periodic) {
            cov.calc.dists(tau, opts$period)
        } else {
            cov.calc.dists(tau)
        }
    })
}

#' Calculate covariance structure for gene over pseudotimes and test
#' inputs.
#'
#' @param dl de.lorean object
#' @param gene.idx Gene index
#' @param cov.fn Covariance function (defaults to cov.matern.32)
#' @param tau The pseudotimes to use
#' @param include.test Also include the pseudotimes for the test inputs
#' @param psi Temporal variation
#' @param omega Noise
#'
cov.calc.gene <- function(dl,
                          gene.idx,
                          cov.fn=cov.matern.32,
                          tau=tau.for.sample(dl),
                          include.test=TRUE,
                          psi = sampled.gene.param(dl, gene.idx, "psi"),
                          omega=sampled.gene.param(dl, gene.idx, "omega"))
{
    with(dl, {
        r <- cov.calc.dl.dists(dl,
                               tau=tau,
                               include.test=include.test)
        (
            psi * cov.fn(r, opts$length.scale)  # Structure
            + omega * diag(nrow(r))  # Noise
        )
    })
}

#' Calculate covariance for gene over test inputs when conditioned on
#' data at estimated pseudotimes.
#'
#' @param dl de.lorean object
#' @param gene.idx Gene index
#' @param cov.fn Covariance function (defaults to cov.matern.32)
#' @param tau The pseudotimes to use
#'
cov.calc.gene.conditioned <- function(dl,
                                      gene.idx,
                                      cov.fn=NULL,
                                      tau=tau.for.sample(dl))
{
    if (is.null(cov.fn)) {
        cov.fn <- cov.matern.32
    }
    with(dl, {
        Sigma <- cov.calc.gene(dl, gene.idx, cov.fn, tau=tau)
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
#' @param tau The pseudotimes to use
#'
cov.all.genes.conditioned <- function(dl,
                                      cov.fn=NULL,
                                      tau=tau.for.sample(dl))
{
    vapply(1:nrow(dl$gene.map),
           function(gene.idx) cov.calc.gene.conditioned(dl,
                                                        gene.idx,
                                                        cov.fn,
                                                        tau=tau),
           FUN.VALUE=matrix(0,
                            nrow=length(dl$test.input),
                            ncol=length(dl$test.input)))
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
gp.log.marg.like <- function(y, K=NULL, U=chol(K)) {
    alpha <- backsolve(U, backsolve(U, y, transpose = TRUE))
    -(
        (t(y) %*% alpha) / 2
        + sum(diag(U))
        + length(y) * log(2 * pi)
    )
}

#' Predictive mean, variance and log marginal likelihood of a GP.
#' See "2.3 Varying the Hyperparameters"
#' on page 19 of Rasmumssen and Williams' book.
#'
#' @param y The targets.
#' @param K The covariance matrix (kernel) for input points, not needed if U is provided.
#' @param Kstar The cross covariance matrix (kernel)
#' @param Kstarstar The cross covariance matrix (kernel) for test points
#' @param U Cholesky decomposition of K
#'
gp.predict <- function(y, K=NULL, Kstar, Kstarstar, U=chol(K)) {
    alpha <- backsolve(U, backsolve(U, y, transpose = TRUE))
    v <- backsolve(U, Kstar, transpose = TRUE)
    list(
        mu = t(Kstar) %*% alpha,  # fstar
        Sigma = Kstarstar - t(v) %*% v,  # Vstar
        logpy = -(
            (t(y) %*% alpha) / 2
            + sum(diag(U))
            + length(y) * log(2 * pi)
        )
    )
}

#' Convert the output of gp.predict() into a data.frame.
#'
#' @param predictions The predictions
#'
gp.predictions.df <- function(predictions) {
    with(predictions, data.frame(mu=mu,
                                 Sigma=diag(Sigma),
                                 idx=1:length(mu)))
}

#' Condition a Guassian on another.
#' See Eqn. A.6
#' on page 200 of Rasmumssen and Williams' book.
#'
#' @param y y
#' @param mu.x Mean of x
#' @param mu.y Mean of y
#' @param .A Var(X)
#' @param .B Var(Y)
#' @param .C Cov(X, Y)
#' @param U Cholesky decomposition of .B
#'
gaussian.condition <- function(
    y,
    .A,
    .B,
    .C,
    mu.x=rep(0, nrow(.A)),
    mu.y=rep(0, nrow(.B)),
    U=chol(.B))
{
    alpha <- backsolve(U, t(.C))
    print(dim(alpha))
    print(dim(.C))
    print(dim(alpha %*% .C))
    print(dim(.A))
    list(mu = mu.x + as.vector((y - mu.y) %*% alpha),
         Sigma = .A - .C %*% alpha)
}

#' The expected within sample variance of a Gaussian with the given covariance.
#'
#' @param K Covariance
#'
#' @export
#'
expected.sample.var <- function(K) mean(diag(K)) - mean(K)
