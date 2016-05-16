library(TreePar)
library(TreeSim)
library(NELSI)
library(doParallel)
library(foreach)
##
# Simulate tree, fit parameters using one and two shifts, find best model using LRT,
#get optimal parameters for one and two shifts and simulate 100 trees with each

# Fit rate shifts
fit_shifts<- function(tr, rho){
    x_times <- sort(intnode.times(tr), decreasing = T)
    start <- min(x_times)
    end <- max(x_times)
    grid <- diff(range(x_times))/50
    res <- bd.shifts.optim(x_times, sampling = rho, grid, start, end, posdiv = T)
    res[[2]]
}

# Pull out parameters for both numbers of shifts, shift time, and best model.
get_params <- function(fit_model){
    likelihoods <- sapply(fit_model, function(x) x[1])
    test <- pchisq(2 * likelihoods[1] - likelihoods[2], 3)
    turnovers <- c(fit_model[[1]][2], fit_model[[2]][c(2, 4)])
    div_rates <-  c( fit_model[[1]][3], fit_model[[2]][c(3, 5)])
    shift_time <- fit_model[[2]][6]
    lambdas <- sapply(1:length(div_rates), function(x) div_rates[x] / (1-turnovers[x]))
    mus <- lambdas * div_rates



    best_n_shifts<- as.numeric(pchisq(2*(likelihoods[1] - likelihoods[2]), 3) > 0.95)
    return(list(likelihoods = likelihoods, turnovers = turnovers, div_rates = div_rates, shift_time = shift_time, lambdas= lambdas, mus = mus, best_n_shifts = best_n_shifts))
}


set.seed(432134)
nspecies <- 50
time <- c(0, 0.4)
rho <- c(0.5, 1) #sampling proportion
lambda <- c(2, 5)#spec rate
mu <- c(1.5, 0.5)#extinction rate


rateshift_trees <- sim.rateshift.taxa(nspecies, 10, lambda = lambda, mu = mu, frac = rho, times = time, complete = F)

tr_fit <- fit_shifts(rateshift_trees[[1]], rho = rho)

tr_fit_params <- get_params(tr_fit)

# I am having trouble resimulating data because mu is sometimes greater than lambda, so the trees cannot be simulated


# Simulated parameters
mu/lambda
lambda - mu

plot(rateshift_trees[[1]], show.tip.label = F)
node_times <- intnode.times(rateshift_trees[[1]])
nodelabels(round(node_times, 2), bg = 'white')
lines(x = c(max(node_times)- tr_fit_params$shift_time, max(node_times) - tr_fit_params$shift_time), c(0, 100), col = rgb(1, 0, 0, 0.5), lwd = 6)

sim_trees_one_shift <-sim.rateshift.taxa(nspecies, numbsim = 10, mu = tr_fit_params$mus[1], lambda = tr_fit_params$lambdas[1], frac = rho, times = tr_fit_params$shift_time, complete = F)




 # fit model,
# pull out parameters estimated and save in a matrix with: likelihoods of the two models, model selected, simulated values, estimated values, parametric bootstrap likelihood hpd for one shift, parametric bootstrap likelihood hpd for no shifts.

#loop_replicate <- function(tree, rho){
# Include parallelisation here


