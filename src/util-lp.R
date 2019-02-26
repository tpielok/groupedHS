# common variables

get_lp_file <- function(){
  return("data/lp-allresults.RData")
} 

get_lp_scenarios <- function(){
  # todo move priors out of dgp
  scenarios <- expand.grid(
    snr = c(0.1, 1, 5),
    sparse = c("high", "lo"),
    p = c("strict", "opt", "loose"),
    n = c(100)
  ) 
  return(scenarios)
}
