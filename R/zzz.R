.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  # message('.onLoad start')
  # print(modules)
  # print(stanmodels)
  for (m in modules) loadModule(m, what = TRUE)
  # message('.onLoad end')
}

.onAttach <- function(...) {
  packageStartupMessage("- For execution on a local, multicore CPU with excess RAM we recommend calling")
  packageStartupMessage("options(mc.cores = parallel::detectCores())")
}
