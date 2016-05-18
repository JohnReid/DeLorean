[![Travis-CI Build Status](https://travis-ci.org/JohnReid/DeLorean.svg?branch=master)](https://travis-ci.org/JohnReid/DeLorean)

DeLorean
========

R package to model time series accounting for noise in the temporal dimension.
Specifically designed for single cell transcriptome experiments.



Requirements
------------

To render the vignettes you will need a working version (> 1.12.3) of `pandoc` on your machine with the
`pandoc-citeproc` filter. On Ubuntu do:

    sudo apt-get install pandoc pandoc-citeproc



Installation from CRAN
----------------------

Just run the following in an `R` session:

    install.packages('DeLorean')


Installation from source
------------------------

If you prefer to have the very latest version you can install `DeLorean` from source.
If you do not already have [`devtools`](https://github.com/hadley/devtools)
installed, then install it by running:

    install.packages('devtools')

Now you can install the development version of `DeLorean` with:

    devtools::install_github('JohnReid/DeLorean')


Documentation
-------------

Read the vignette:

    vignette('DeLorean')
