## Version 1.4.0
- Make changes to how compiled models are cached to satisfy silly CRAN policy. Now compilation happens more often:
  once per model per R session.

## Version 1.3.0
- Remove filter.genes and filter.cells, (renamed to use underscores) in order to pass CRAN checks.

## Version 1.2.6
- Deprecate filter.genes and filter.cells, (renamed to use underscores).

## Version 1.2.5
- Update rstan requirement to 2.14.1 to avoid bug: http://andrewgelman.com/2017/01/01/stan-2-14-released-r-python-fixes-bug-sampler/

## Version 1.2.4
- Increase border of inducing points in sparse approximations
- Fix further Stan bug (feature?) with square root of integers
- Update rstan requirement to 2.12.1 from bad 2.10 version: http://andrewgelman.com/2016/07/31/stan-2-11-good-stan-2-10-bad/

## Version 1.2.3
- Update for new package versions: rstan (2.10.1) and dplyr (0.5.0).

## Version 1.2.2
- Updated code for revision of manuscript. Changed cell size normalisation
  strategy.

## Version 1.2.1
- Fix for sqrt(integer) bug on Solaris x86.

## Version 1.2.0
- First version submitted to CRAN.
