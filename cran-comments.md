## Test environments
* local ubuntu 16.04, R 3.2.5
* travis-CI ubuntu 12.04 (old-rel, release and devel)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or NOTEs.

There was one WARNING:

* checking S3 generic/method consistency ... WARNING
filter:
  function(x, filter, method, sides, circular, init)
filter.cells:
  function(dl, .filter, number, cells)

filter:
  function(x, filter, method, sides, circular, init)
filter.genes:
  function(dl, .filter, number, genes)

See section ‘Generic functions and methods’ in the ‘Writing R
Extensions’ manual.

Found the following apparent S3 methods exported but not registered:
  filter.cells filter.genes
See section ‘Registering S3 methods’ in the ‘Writing R Extensions’
manual.

I am deprecating two non-S3 methods and replacing them with methods that
use underscores rather than dots. Until I can remove I have not been able
to work out how to get around this WARNING.
