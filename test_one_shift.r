
library(TreePar)
library(TreeSim)
library(NELSI)

#set.seed(10)
nspecies <- 200
time <- c(0, 0.5, 1) # At time 1 in the past, there is a rate shift
rho <- c(0.5, 0.5, 1) #half of the present day species are sampled (rho[1] = 0.5)
lambda <- c(1.5, 10, 5)# speciation rates, between t[i] and t[i+1] we have a speciation rate lambda
mu <- c(1.5, 1.5, 1.5)# extinction rate. Similar notation as lambda

#Simulate a tree with a two rate shifts
tree <- sim.rateshift.taxa(nspecies, 1, lambda = lambda, mu = mu, frac = rho, times = time, complete = F)


pdf('output_graphs.pdf', useDingbats = F)
plot(tree[[1]], show.tip.label = F)
nodelabels(round(intnode.times(tree[[1]]), 2))

## Extract speciation times from tree
x <- sort(intnode.times(tree[[1]]), decreasing = T)

# When estimateing the rate shift times, t, based on branching times, x, we allow the shift times to be:
# 0.6, 0.8, 1, 1.2, ... 2.4
start <- 0.4
end <- 2
grid <- 0.2
print(seq(from = start, to = end, by = grid))
res <- bd.shifts.optim(x, c(rho, 1), grid, start, end)[[2]]

#Do two rate shifts explain the data better than one rate shift?
test<-pchisq(2*(res[[2]][1]-res[[2]][1]),3)
test

print(res)
print('Turnovers')
print(res[[2]][2:3])
print('Net diversification')
print(res[[2]][4:5])

# simulation vaules:
# Turnover (extinction/speciation):
mu/lambda

# net diversification
# speciation - extinction
lambda - mu

# Fit 0 rate shift model (i.e very underparameterised)
# Calculate lambda and mu:
turnover_est <- res[[1]][2]
net_div_est <- res[[1]][3]

lambda_est <- net_div_est / (1 - turnover_est)
mu_est <- lambda_est * turnover_est

lambda_est
mu_est

# Simulate 200 trees under a constant rate shift model:
constant_sim_trees <- sim.bd.taxa(n = 200, numbsim = 100, lambda = lambda_est, mu = mu_est, frac = rho[1], complete = F)

start <- 10
end <- 100
grid <- 5

# Fit a single rate to trees (note that I am using the correct sampling):
constant_likelihoods <- vector()
for(i in 1:length(constant_sim_trees)){
    tr_times <- sort(intnode.times(constant_sim_trees[[i]]), decreasing = T)
    lik_temp <- bd.shifts.optim(tr_times, sampling = rho, start = start, end = end, grid = grid)[[2]]
    constant_likelihoods[i] <- lik_temp[[1]][1]
}

hist(constant_likelihoods)
print(res)

dev.off()
