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
  # print(pca)
  # print(pca)
  # summary(pca)
  #
  # Scree plot
  # plot(pca)
  # pdf('Plots/Guo-PCA-scree.pdf', width = 12, height = 8)
  # plot(pca)
  # dev.off()
  #
  # Bi-plot
  # biplot(pca)
  # pdf('Plots/Guo-PCA-biplot.pdf', width = 12, height = 12)
  # biplot(pca)
  # dev.off()
  #
  # Convert to data frame
  pca.l <- reshape2::melt(pca$x, varnames=c("cell", "PC"), value.name="x")
  # sample_n(pca.l, 10)
  pca.df <-
    pca.l %>%
    reshape2::dcast(cell ~ PC, value.var = 'x') %>%
    left_join(guo.cell.meta)
  #
  # Calculate direction that best predicts capture time
  lm.fit <- lm(obstime ~ PC1 + PC2, pca.df)
  lm.fit
  a <- lm.fit$coefficients[['PC1']]
  b <- lm.fit$coefficients[['PC2']]
  lambda <- sqrt(a**2 + b**2) / 2
  slope.x <- a / 2 * lambda
  slope.y <- b / 2 * lambda
  capturetime.slope <- slope.y / slope.x
  capturetime.angle <- atan2(slope.y, slope.x)
  # with(pca.df, cor(obstime, PC1)) / with(pca.df, cor(obstime, PC2))
  # gp.pca <-
  #   ggplot(pca.df, aes(x = PC1, y = PC2, label = cell, colour = capture, shape = cell.type)) +
  #   #geom_text() +
  #   geom_point(size = 4) +
  #   geom_abline(slope = capturetime.slope)
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
