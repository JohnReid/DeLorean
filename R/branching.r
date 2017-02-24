#
# Code to implement branching processes.
#

# Convert x,y coordinates to radial
#
to.radial <- function(x, y) {
  r = sqrt(x**2 + y**2)
  list(r = r, theta = atan2(y, x))
}

#' Perform a PCA on the expression data. Calculate which direction in PC1-PC2 space
#' correlates best with capture time.
#'
#' @param dl The de.lorean object.
#'
#' @export
#'
pca.expr <- function(dl) with(dl, {
  pca <- prcomp(t(expr), center=TRUE, scale.=FALSE)
  names(pca)
  summary(pca)
  #
  # Convert to data frame
  pca.l <- reshape2::melt(pca$x, varnames=c("cell", "PC"), value.name="x")
  # sample_n(pca.l, 10)
  pca.df <-
    pca.l %>%
    reshape2::dcast(cell ~ PC, value.var = 'x') %>%
    left_join(cell.meta)
  #
  # Calculate direction that best predicts capture time
  lm.fit <- lm(obstime ~ PC1 + PC2, pca.df)
  a <- lm.fit$coefficients[['PC1']]
  b <- lm.fit$coefficients[['PC2']]
  lambda <- sqrt(a**2 + b**2) / 2
  slope.x <- a / 2 * lambda
  slope.y <- b / 2 * lambda
  capturetime.slope <- slope.y / slope.x
  capturetime.angle <- atan2(slope.y, slope.x)
  dist.to.capturetime.axis <- function(x, y) {
    with(to.radial(x, y), {
      radians.diff <- (capturetime.angle - theta) %% (2 * pi)
      sin(radians.diff) * r
    })
  }
  pca.df <-
    pca.df %>%
    mutate(dist = dist.to.capturetime.axis(PC1, PC2),
           z.hat = (dist - mean(dist)) / sd(dist))
  dl$pca <- pca
  dl$lm.fit <- lm.fit
  dl$capturetime.slope <- capturetime.slope
  dl$capturetime.angle <- capturetime.angle
  dl$pca.df <- pca.df
  dl
})


#' Prepare to calculate posteriors across the latent space in the branching model.
#'
#' @param dl The DeLorean object
#' @param sample.iter The sample (iteration) to use
#' @param tau.grid.size The number of grid entries to use in tau dimension
#' @param z.grid.size The number of grid entries to use in z dimension
#'
#' @export
#'
branching.prepare.posterior <- function(dl, sample.iter = dl$best.sample, tau.grid.size = 51, z.grid.size = 51) {
  #
  # To remove CRAN check problem
  cov_symmetric <- NA
  rm(cov_symmetric)
  #
  dl$best.m <- lapply(dl$samples.l, function(s) filter(s, iter == sample.iter))
  dl$tau.z <- with(dl$best.m, left_join(tau, z))
  #
  # Expose Stan functions
  dl$stan.fns <- rstan::expose_stan_functions(dl$fit)
  #
  # Get range for grid
  tau.range <- range(dl$best.m$tau$tau)
  tau.min <- floor(tau.range[1])
  tau.max <- ceiling(tau.range[2])
  tau.grid <- seq(from = tau.min, to = tau.max, length.out = tau.grid.size)
  dl$tau.step <- tau.grid[2] - tau.grid[1]
  z.range <- range(dl$best.m$z$z)
  z.min <- floor(z.range[1])
  z.max <- ceiling(z.range[2])
  z.grid <- seq(from = z.min, to = z.max, length.out = z.grid.size)
  dl$z.step <- z.grid[2] - z.grid[1]
  #
  # Expand into grid
  dl$grid.df <- expand.grid(tau = tau.grid, z = z.grid)
  grid <- t(as.matrix(dl$grid.df))
  #
  # Covariance across grid points
  dl$K.grid <- cov_symmetric(grid, dl$stan.data$lengthscales)
  #
  # Covariance across latent points from posterior
  cell.post <- with(dl$best.m, left_join(tau, z))
  points.post <- t(as.matrix(select(cell.post, tau, z)))
  dl$K.post <- cov_symmetric(points.post, dl$stan.data$lengthscales)
  #
  # Cross-covariance
  dl$K.cross <- cov(points.post, grid, dl$stan.data$lengthscales)
  #
  # Return DeLorean object
  dl
}


#' Calculate the posterior for several genes in the branching model.
#'
#' @param dl The DeLorean object
#' @param genes The genes to calculate the posterior for
#'
#' @export
#'
branching.genes.post <- function(dl, genes) with(dl, {
  with(dl$best.m, left_join(psi, omega) %>% left_join(dl$gene.map)) %>%
    filter(gene %in% genes) %>%
    group_by(g) %>%
    do(branching.gene.post(dl, .$g, .$psi, .$omega)) %>%
    ungroup() %>%
    left_join(dl$gene.map)
})


#' Calculate the posterior for one gene in the branching model.
#'
#' @param dl The DeLorean object
#' @param g The index of the gene in the expression matrix (gene.map)
#' @param psi The temporal variation of the gene
#' @param omega The noise level of the gene
#'
#' @export
#'
branching.gene.post <- function(dl, g, psi, omega) with(dl, {
  y <- dl$expr[g,]
  gp.post <- gp.predict(
    y - mean(y),
    psi * K.post + diag(omega, nrow = nrow(K.post)),
    psi * K.cross,
    psi * K.grid + diag(omega, nrow = nrow(K.grid)))
  dl$grid.df %>% mutate(post.mean = as.vector(gp.post$mu))
})


#' Calculate the log marginal likelihoods of each cell's expression
#' of each gene.
#'
#' @param dl The DeLorean object.
#'
#' @export
#'
branching.calc.log.marg.lik <- function(dl) {
  #
  # To remove CRAN check problem
  cov_symmetric <- NA
  rm(cov_symmetric)
  #
  cell.post <- with(dl$best.m, left_join(tau, z))
  points.post <- t(as.matrix(select(cell.post, tau, z)))
  K.post <- cov_symmetric(points.post, dl$stan.data$lengthscales)
  gene.log.marg.lik <- function(g1) {
    psi = filter(dl$best.m$psi, g == g1)$psi
    omega = filter(dl$best.m$omega, g == g1)$omega
    # message(psi, ' ', omega)
    y = dl$stan.data$expr[g1,]
    K = psi * K.post + diag(omega, nrow = length(y))
    gp.log.marg.like.individual(y, K = K)
  }
  t(sapply(1:dl$stan.data$G, gene.log.marg.lik))
}
