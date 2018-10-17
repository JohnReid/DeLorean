## Test environments
* local ubuntu 16.04, R 3.5.0
* travis-CI ubuntu 14.04 (old-rel, release and devel)
* win-builder (devel and release)


## R CMD check results
There were no ERRORs.

There was one WARNING:

```
Compilation used the following non-portable flag(s):
  ‘-Wdate-time’ ‘-Werror=format-security’ ‘-Wformat’
```


There were four NOTEs:

```
Maintainer: ‘John Reid <johnbaronreid@gmail.com>’

New maintainer:
  John Reid <johnbaronreid@gmail.com>
Old maintainer(s):
  John Reid <john.reid@mrc-bsu.cam.ac.uk>
```
I will send CRAN an email confirming this change.

```
checking installed package size ... NOTE
  installed size is 79.0Mb
  sub-directories of 1Mb or more:
    libs  76.3Mb
```

```
checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.
```
