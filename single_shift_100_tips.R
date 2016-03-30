library(TreePar)
library(TreeSim)
library(NELSI)
library(doParallel)
library(foreach)

set.seed(432134)
nspecies <- 100
time <- c(0, 0.5)
rho <- c(1.0, 0.5) #all of the present day are sampled?
lambda <- c(1.5, 4)
mu <- c(1.5, 0.1)
rateshift_trees <- sim.rateshift.taxa(nspecies, 10, lambda = lambda, mu = mu, frac = rho, times = time, complete = F)

# Fit rate shifts
fit_shifts<- function(tr, rho){
    x_times <- sort(intnode.times(tr), decreasing = T)
    start <- min(x_times)
    end <- max(x_times)
    grid <- diff(range(x_times))/10
    res <- bd.shifts.optim(x_times, sampling = rho, grid, start, end, posdiv = T)
    res[[2]]
}

                                        # Pull out parameters for both numbers of shifts, shift time, and best model.
#Debug HERE
#get_params <- function(fit_model){
    likelihoods <- sapply(fit_model, function(x) x[1])
    test <- pchisq(2 * likelihoods[1] - likelihoods[2], 3)
    turnovers <- c(fit_model[[1]][2], fit_model[[2]][c(2, 4)])
    div_rates <-  c( fit_model[[1]][3], fit_model[[2]][c(3, 5)])
    shift_time <- fit_model[[2]][6]
    lambdas <- div_rates / (1-turnovers)
    mus <- lambdas * div_rates
  #  return(list(likelihoods, turnovers, div_rates, shift_time, lambdas, mus, test))
#}


lt1 <- fit_shifts(rateshift_trees[[1]], rho = rho)
mu/lambda
lambda - mu

 # fit model,
# pull out parameters estimated and save in a matrix,
# use parameter estimates to simulate 100 trees under two models

#loop_replicate <- function(tree, rho){


fit_model <-



