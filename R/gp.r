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


