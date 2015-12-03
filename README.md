DeLorean
========

R package to model time series accounting for noise in the temporal dimension.
Specifically designed for single cell transcriptome experiments.



Requirements
------------

To render the vignettes you will need a working version (> 1.12.3) of `pandoc` on your machine with the
`pandoc-citeproc` filter. On Ubuntu do:

    sudo apt-get install pandoc pandoc-citeproc



Installation
------------

If you do not already have [`devtools`](https://github.com/hadley/devtools)
installed, then install it by running:

    install.packages('devtools')

Now you can install the development version of `DeLorean` with:

    devtools::install_github('JohnReid/DeLorean')

