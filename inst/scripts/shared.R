#
# Code that is shared amongst the vignettes.
#
library(ggplot2)
library(ggthemes)
library(dplyr)
library(reshape2)
library(magrittr)
theme_set(theme_few())

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
loadfonts()
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
bioinf.single.w <- 3.38  # Width of single bioinformatics column in inches
bioinf.double.w <- 7.00  # Width of double bioinformatics column in inches
bioinf.single.h <- bioinf.single.w / golden.ratio
bioinf.double.h <- bioinf.double.w / golden.ratio
single.col.width <- 86  # Bioinformatics width in mm
double.col.width <- 178 # Bioinformatics width in mm

#
# Theme for Bioinformatics
#
bioinf.size <- 8
bioinf.family <- "Helvetica"
#
# We use ggthemes::theme_few() if it is installed.
if ('ggthemes' %in% rownames(installed.packages())) {
  bioinf.config <-
    list(ggthemes::theme_few(base_size=bioinf.size, base_family=bioinf.family),
         ggthemes::scale_colour_few(),
         ggthemes::scale_fill_few())
} else {
  bioinf.config <- ggplot2::theme_grey(base_size=bioinf.size, base_family=bioinf.family)
}
bioinf.sizes <-
  list(width=bioinf.single.w,
       height=bioinf.single.h,
       units="in")
suppl.sizes <-
  list(width=345/72,
       height=345/72/golden.ratio,
       dpi=350,
       units="in")

#
# Theme for PLoS
#
base.family <- "Times"
# base.family <- "Arial"  # Arial may not be on all systems
# base.family <- "Liberation Sans"
plos.theme <- ggplot2::theme_classic(base_size=8, base_family=base.family)

#
# Fit model if not defined
#
if (! exists('fit.model')) fit.model <- TRUE

#
# Create data directory if not there
#
if( ! file.exists('Data')) dir.create('Data')

#
# Retrieve N largest objects in global environment.
largest.objects <- function(N = 10, env = globalenv()) {
    z <- sapply(
                    ls(env),
                        function(x) object.size(get(x, env)))
  as.matrix(rev(sort(z))[1:N]) %>% set_colnames('bytes')
}
