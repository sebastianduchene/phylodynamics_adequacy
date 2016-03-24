library(TreePar)
library(TreeSim)
library(NELSI)
library(doParallel)
library(foreach)

# Simulate 10 trees with one rate shift
set.seed(1234)
nspecies <- 500
time <- c(0, 0.5)
rho <- c(1, 0.5)
lambda <- c(1.5, 4)
mu <- c(1.5, 0)
rateshift_trees <- sim.rateshift.taxa(nspecies, 10, lambda = lambda,
                                     mu = mu, frac = rho, times = time, complete = F)

# Check that the trees have a root node age of more than 0.5 (otherwise there is no rate shift)
plot(rateshift_trees[[1]], show.tip.label = F)
nodelabels(round(intnode.times(rateshift_trees[[1]]), 2))
print(  sapply(rateshift_trees, function(x) max(intnode.times(x)))  )

# Check that rate shifts can be estimated for this trees

# Check that the rateshifts can be estimated
tr <- rateshift_trees[[1]]
x_times <- sort(intnode.times(tr), decreasing = T)
start <- min(x_times)
end <- max(x_times)
grid <- diff(range(x_times))/20
res <- bd.shifts.optim(x_times, sampling = c(1, 0.5), grid, start, end, posdiv = T)
res[[2]]


fit_rate_shifts <- function(tree, rho){ # Rho at present and in the past. This parameter needs to be fixed

    x_times <- sort(intnode.times(tree), decreasing = T)
    start <- min(x_times)
    end <- max(x_times)
    grid <- diff(range(x_times))/20
    res <- bd.shifts.optim(x_times, rho, grid, start, end, posdiv = T)[[2]]

    # Find likelihoods, lambda, mu, and rate-shift times
    likelihoods <- sapply(res, function(x) x[1])

    lambda0 <- res[[1]][3] / (1 - res[[1]][2]) # These are the lambda and mu estimates from turover and net
    mu0 <- lambda0 * res[[1]][2]               # speciation for 0 rate shifts. Please check.

    # The following are also computed, but note that some of them are negative and that this might.
    # I couldn't simulate trees using these parameters, maybe because of the negative values?
    lambda11 <- res[[2]][3] / (1 - res[[2]][2])
    mu11 <- lambda11 * res[[2]][2]
    lambda12 <- res[[2]][5] / (1 - res[[2]][4])
    mu12 <- lambda12 * res[[2]][4]
    time1 <- res[[2]][length(res[[2]])]

    return(list(likelihoods, shifts0= c(lambda0, mu0), shifts1=c(lambda11, lambda12, mu11, mu12, time1)))
}

pvals <- vector()
likelihoods_distros <- list()
likelihoods_empirical <- vector()
empirical_tree_param_estimates <- list()

cl <- makeCluster(8)
registerDoParallel(cl)
for(tr in 1:length(rateshift_trees)){
    reference_estimates <- fit_rate_shifts(rateshift_trees[[tr]], rho)

    sim_trees0 <- sim.bd.taxa(n = nspecies, numbsim = 100, lambda = reference_estimates$shifts0[1],
                              mu = reference_estimates$shifts0[2], frac = 0.5, complete = F)
    liks_sim_trees0 <- foreach(mt = sim_trees0, .packages = c('NELSI', 'TreePar')) %dopar% fit_rate_shifts(mt, rho)[[1]][1]
    likelihoods_distros[[tr]] <- liks_sim_trees0
    likelihoods_empirical[tr] <- reference_estimates[[1]][1]
    empirical_tree_param_estimates[[tr]] <- reference_estimates$shifts0
    pvals[tr] <- sum(reference_estimates[[1]][1] > liks_sim_trees0)

}
stopCluster(cl)

pdf('histograms_500_tips.pdf')
par(mfrow = c(3, 3))
for(i in 1:9){
    hist(as.numeric(likelihoods_distros[[i]]), main = '', ylab = '', xlab = '', col = rgb(0, 0, 0.5, 0.3))
    lines(x = c(likelihoods_empirical[i], likelihoods_empirical[i]), y = c(0, 20), col = 'red', lwd = 2)
}
dev.off()
