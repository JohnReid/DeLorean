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


#' Find a good ordering in the sense that some function is
#' locally maximised.
#'
#' @param fn A function to maximise.
#' @param ordering The permutation (ordering) to start from.
#'
#' @export
#'
ordering.maximise <- function(fn, ordering) {
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
