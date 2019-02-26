# common variables

get_sf_file <- function(){
  return("data/sf-allresults.RData")
} 

get_sf_scenarios <- function(){
  # todo move priors out of dgp
  scenarios <- expand.grid(
    snr = c(5, 20),
    sparse = c("high", "lo"),
    corr = c(0, .7),
    p = c("strict", "opt", "loose"),
    n = c(200)
  ) #, 1000))
  scenarios$dgp <- 1:nrow(scenarios)
  return(scenarios)
}