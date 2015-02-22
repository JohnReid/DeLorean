#' Move one item in an ordering and shift the other items.
#'
#' @param ordering The ordering.
#' @param from The index of the item to move.
#' @param to The index it should be moved to.
#'
#' @export
#'
ordering.move <- function(ordering, from, to) {
    if (from != to) {
        item <- ordering[from]
        if (from < to - 1) {
            ordering[from:(to-1)] <- ordering[(from+1):to]
        } else if (from > to + 1) {
            ordering[(to+1):from] <- ordering[to:(from-1)]
        } else {
            # just swap them
            ordering[from] <- ordering[to]
        }
        ordering[to] <- item
    }
    ordering
}


#' Randomly move one item in an ordering to another location
#'
#' @param ordering The ordering.
#'
#' @export
#'
ordering.random.move <- function(ordering) {
    .sample <- sample(length(ordering), 2)
    ordering.move(ordering, .sample[1], .sample[2])
}


#' Find a good ordering in the sense that some function is
#' locally maximised.
#'
#' @param fn A function to maximise.
#' @param ordering The permutation (ordering) to start from.
#'
#' @export
#'
ordering.maximise <- function(ordering, fn) {
    while (TRUE) {
        ordering.new <- ordering.improve(fn, ordering)
        # If we didn't improve, then break the loop
        if (all(ordering.new == ordering)) {
            break
        } else {
            ordering <- ordering.new
            message('Improved ordering: f(ordering) = ', fn(ordering))
        }
    }
    ordering
}


#' Metropolis-Hastings on orderings.
#'
#' @param init.ordering Initial ordering
#' @param loglikelihood Log likelihood function
#' @param proposal.fn Proposal function
#'
#' @export
#'
ordering.metropolis.hastings <- function(
    init.ordering,
    log.likelihood,
    proposal.fn=ordering.random.move,
    iterations=1000)
{
    chain <- array(as.integer(0), dim=c(iterations+1, length(init.ordering)))
    chain[1,] <- init.ordering
    last.ll <- log.likelihood(chain[1,])
    for (i in 1:iterations) {
        proposal <- proposal.fn(chain[i,])
        this.ll <- log.likelihood(proposal)
        probab <- exp(this.ll - last.ll)
        if (runif(1) < probab) {
            chain[i+1,] <- proposal
            last.ll <- this.ll
        } else {
            chain[i+1,] <- chain[i,]
        }
    }
    mcmc(chain, start=1, end=iterations)
}


#' Improve the ordering in the sense that some function is
#' maximised.
#'
#' @param fn A function to maximise.
#' @param ordering The permutation (ordering) to start from.
#'
#' @export
#'
ordering.improve <- function(fn, ordering) {
    Reduce(
        # Reduction function
        function(ordering, from) {
            scores <- sapply(1:length(ordering),
                            function(to) fn(ordering.move(ordering, from, to)))
            best.to <- which.max(scores)
            ordering.move(ordering, from, best.to)
        },
        # Choose all the elements in a random order
        sample(length(ordering)),
        # Initialise the reduction with the initial ordering
        init=ordering)
}


#' Invert the ordering
#'
#' @param ordering The permutation (ordering) to invert.
#'
#' @export
#'
ordering.invert <- function(ordering) {
    result <- rep(0, length(ordering))
    for (n in 1:length(result)) {
        result[ordering[n]] <- n
    }
    result
}


#' Test ordering score: sum every time consecutive items are in
#' order.
#'
#' @param ordering
#'
#' @export
#'
ordering.test.score <- function(ordering) {
    sum(sapply(1:(length(ordering)-1),
               function(n) as.integer(ordering[n] < ordering[n+1])))
}
