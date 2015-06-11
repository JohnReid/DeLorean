#
# Install extrafont and import fonts from system.
# We just need to do this once.
#
# install.packages("extrafont")
# library(extrafont)
# font_import()
# fonts()

# Load fonts for plotting
library(extrafont)
# loadfonts()
loadfonts(device="postscript")

# font_import(pattern="[A/a]rial")
# loadfonts(device="win")

#
# Widths for saving figures
#
text.width <- 4.7  # LaTeX width in inches
golden.ratio <- 1.618  # Pleasing ratio
fig.width <- text.width
fig.height <- text.width / golden.ratio
plos.width <- 7.2 # LaTeX width for PLoS paper in inches
plos.height <- plos.width / golden.ratio
html5 <- list(width=1300, height=700)  # in pixels
html5$ratio <- with(html5, width / height)
slide.fig.width <- 7

#
# Theme for PLoS
#
# base.family <- "Times"
base.family <- "Arial"
# base.family <- "Liberation Sans"
plos.theme <- theme_classic(base_size=8, base_family=base.family)
