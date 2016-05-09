#
# Helper functions for analysis of McDavid data.
#

#
# CycleBase cbPeaktimes of 0 (equiv. 100) correspond to the transition
# between M and G1.
# (http://cyclebase2.jensenlab.org/methods.jsp#peaktime-uncertainty)
# Hence cbPeaktime = 0 corresponds to a pseudotime of .5
#

#
# Make times periodic
periodise <- function(tau, period=100) tau - floor(tau/period) * period
# tau.test <- seq(-.5, 3.5, by=.5)
# tau.test
# periodise(tau.test)

#
# Map from pseudotimes (tau) to cbPeaktimes
tau.to.cbtime <- function(tau) 100*periodise((tau - .5)/3, 1)
#
# and vice versa
cbtime.to.tau <- function(cbtime) periodise(cbtime/100*3+.5, 3)
#
# Test mapping functions
tau.test <- seq(-.5, 3.5, by=.5)
sum(abs(periodise(tau.test) - cbtime.to.tau(tau.to.cbtime(tau.test))))
cbtime.test <- seq(-5, 105, by=5)
sum(abs(periodise(cbtime.test, 100) - tau.to.cbtime(cbtime.to.tau(cbtime.test))))

#
# Distance between two points on interval [0, 100] where 0 is connected to 100
peak.distance <- function(peak.1, peak.2) {
    stopifnot(all(peak.1 >= 0.))
    stopifnot(all(peak.1 <= 100.))
    stopifnot(all(peak.2 >= 0.))
    stopifnot(all(peak.2 <= 100.))
    dist <- abs(peak.1 - peak.2)
    ifelse(dist > 50, 100 - dist, dist)
}
#
# Test peak distances
# peak.distance(c(.1, .2, .1, .9, .5, .5), c(.2, .1, .9, .1, .5, .4))

#
# Calculate the root mean square of d
calc.rms <- function(d) sqrt(sum(d**2)/length(d))

#
# Get the peak time for each gene
get.gene.peaks <- function(expr.l) expr.l %>%
  group_by(gene) %>%
  dplyr::summarise(
    peak.idx=which.max(x),  # Find index of maximal expression
    cell=cell[peak.idx],
    capture=capture[peak.idx],
    peak.tau=obstime[peak.idx])
