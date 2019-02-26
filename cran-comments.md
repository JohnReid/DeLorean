## Test environments
* local ubuntu 18.04, R 3.5.2
* travis-CI ubuntu 14.04 (old-rel, release and devel)
* win-builder (devel and release)


## R CMD check results
There were no ERRORs.

There was one WARNING:

```
Compilation used the following non-portable flag(s):
  ‘-Wdate-time’ ‘-Werror=format-security’ ‘-Wformat’
```


There were two NOTEs:

```
checking installed package size ... NOTE
  installed size is 80.6Mb
  sub-directories of 1Mb or more:
    libs  77.9Mb
```

```
checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.
```
