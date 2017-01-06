library(NELSI)
library(apTreeshape)

args <- commandArgs(trailingOnly = T)

print(args)





#f <- args[1]

f <- 'cc_largeNe_391.tree'

file_tag <- gsub('.tree', '', f)

# Empirical__data
tr_empirical <- read.tree(paste0(file_tag, '.tree'))
log_empirical_cc <- read.table(paste0(file_tag, '_cc_1.log'), head = T)
log_empirical_ce <- read.table(paste0(file_tag, '_ce_1.log'), head = T)
log_pps_cc <- read.table(paste0(file_tag, '_cc_pps.log'), head = T)
log_pps_ce <- read.table(paste0(file_tag, '_ce_pps.log'), head = T)

# Check that parameters are properly sampled

pdf(paste0(file_tag, '.pdf'), useDingbats = F, width = 7, height = 10)
par(mfrow = c(5, 3))

plot(density(log_empirical_cc$popSize[-(1:100)]), main = 'popSize')
lines(density(log_pps_cc$popSize[-(1:100)]), col = 'red')
plot(density(log_empirical_ce$ePopSize[-(1:100)]), main = 'ePopSize')
lines(density(log_pps_ce$ePopSize[-(1:100)]), col = 'red')
plot(density(log_empirical_ce$growthRate[-(1:100)]), main = 'growthRate')
lines(density(log_pps_ce$growthRate[-(1:100)]), col = 'red')

# Compute empirical test statistics
emp_df <- get.df(tr_empirical)
emp_rh <- max(intnode.times(tr_empirical))
emp_stemmy <- stemmy(tr_empirical)
emp_colless <- colless(as.treeshape(tr_empirical))
emp_sackin <- sackin(as.treeshape(tr_empirical))
ltt_stats <- get.ltt.summary(tr_empirical)
emp_time_lmax <- ltt_stats[1]
emp_slope_ratio <- ltt_stats[4]
emp_cc_likelihood <- mean(log_empirical_cc$likelihood[-(1:100)])
emp_ce_likelihood <- mean(log_empirical_ce$likelihood[-(1:100)])

##### Load pps and calculate test statistics
pps_trees_ce <- read.nexus(paste0(file_tag, '_ce_pps.trees'))[-(1:200)]
pps_df_ce <- sapply(pps_trees_ce, function(x) get.df(x))
pps_rh_ce <- sapply(pps_trees_ce, function(x) max(intnode.times(x)))
pps_stemmy_ce <- sapply(pps_trees_ce, function(x) stemmy(x))
pps_colless_ce <- sapply(pps_trees_ce, function(x) colless(as.treeshape(x)))
pps_sackin_ce <- sapply(pps_trees_ce, function(x) sackin(as.treeshape(x)))
pps_time_lmax_ce <- vector()
pps_slope_ratio_ce <- vector()
for(i in 1:length(pps_trees_ce)){
    ltt_temp <- get.ltt.summary(pps_trees_ce[[i]])
    pps_time_lmax_ce[i] <- ltt_temp[1]
    pps_slope_ratio_ce[i] <- ltt_temp[4]
}

pps_trees_cc <- read.nexus(paste0(file_tag, '_cc_pps.trees'))[-(1:200)]
pps_df_cc <- sapply(pps_trees_cc, function(x) get.df(x))
pps_rh_cc <- sapply(pps_trees_cc, function(x) max(intnode.times(x)))
pps_stemmy_cc <- sapply(pps_trees_cc, function(x) stemmy(x))
pps_colless_cc <- sapply(pps_trees_cc, function(x) colless(as.treeshape(x)))
pps_sackin_cc <- sapply(pps_trees_cc, function(x) sackin(as.treeshape(x)))
pps_time_lmax_cc <- vector()
pps_slope_ratio_cc <- vector()
for(i in 1:length(pps_trees_cc)){
    ltt_temp <- get.ltt.summary(pps_trees_cc[[i]])
    pps_time_lmax_cc[i] <- ltt_temp[1]
    pps_slope_ratio_cc[i] <- ltt_temp[4]
}

par(mar = c(1, 1, 1, 4))
plot(tr_empirical, show.tip.label = F, type = 'phylogram')
nodelabels(text = round(max(intnode.times(tr_empirical)), 0), node = length(tr_empirical$tip.label)+1, cex = 0.4, bg = 'white', frame = 'circ')
plot(pps_trees_ce[[150]], show.tip.label = F, type = 'phylogram')
nodelabels(text = round(mean(pps_rh_ce), 0), node = length(tr_empirical$tip.label)+1, cex = 0.4, bg = 'white', frame = 'circ')
plot(pps_trees_cc[[150]], show.tip.label = F, type = 'phylogram')
nodelabels(text = round(mean(pps_rh_cc), 0), node = length(tr_empirical$tip.label)+1, cex = 0.4, bg = 'white', frame = 'circ')

par(mar = c(4, 4, 4, 4))
log_ce_files <- dir('.', pattern = paste0(file_tag, '.+ce_pps_.+log'))
pps_likelihood_ce <- sapply(log_ce_files, function(x) mean(read.table(x, head = T)$likelihood[-(1:100)]))
log_cc_files <- dir('.', pattern = paste0(file_tag, '.+cc_pps_.+log'))
pps_likelihood_cc <- sapply(log_cc_files, function(x) mean(read.table(x, head = T)$likelihood[-(1:100)]))

## Pvalues:
pval_df_cc <- round(sum(pps_df_cc>emp_df)/length(pps_df_cc), 2)
pval_df_ce <- round(sum(pps_df_ce>emp_df)/length(pps_df_ce), 2)

pval_rh_cc <- round(sum(pps_rh_cc>emp_rh)/length(pps_rh_cc), 2)
pval_rh_ce <- round(sum(pps_rh_ce>emp_rh)/length(pps_rh_ce), 2)

pval_stemmy_cc <- round(sum(pps_stemmy_cc>emp_stemmy)/length(pps_stemmy_cc), 2)
pval_stemmy_ce <- round(sum(pps_stemmy_ce>emp_stemmy)/length(pps_stemmy_ce), 2)

pval_colless_cc <- round(sum(pps_colless_cc>emp_colless)/length(pps_colless_cc), 2)
pval_colless_ce <- round(sum(pps_colless_ce>emp_colless)/length(pps_colless_ce), 2)

pval_sackin_cc <- round(sum(pps_sackin_cc>emp_sackin)/length(pps_sackin_cc), 2)
pval_sackin_ce <- round(sum(pps_sackin_ce>emp_sackin)/length(pps_sackin_ce), 2)

pval_time_lmax_cc <- round(sum(pps_time_lmax_cc>emp_time_lmax)/length(pps_time_lmax_cc), 2)
pval_time_lmax_ce <- round(sum(pps_time_lmax_ce>emp_time_lmax)/length(pps_time_lmax_ce), 2)

pval_slope_ratio_cc <- round(sum(pps_slope_ratio_cc>emp_slope_ratio)/length(pps_slope_ratio_cc), 2)
pval_slope_ratio_ce <- round(sum(pps_slope_ratio_ce>emp_slope_ratio)/length(pps_slope_ratio_ce), 2)

pval_likelihood_ce <- round(sum(pps_likelihood_ce>emp_ce_likelihood)/length(pps_likelihood_ce), 2)
pval_likelihood_cc <- round(sum(pps_likelihood_cc>emp_cc_likelihood)/length(pps_likelihood_cc), 2)

#par(mfrow = c(3, 3))
hist(pps_df_ce, xlim = range(c(pps_df_ce, pps_df_cc)), col = rgb(0, 1, 0, 0.3),
     main = 'DF', border = rgb(0, 1, 0, 0.3),
     xlab = paste0('Pcc=', pval_df_cc, ' Pce=', pval_df_ce,collapse = ''))
hist(pps_df_cc, add = T, col = rgb(0, 0, 1, 0.3), border = rgb(0, 0, 1, 0.3))
lines(c(emp_df, emp_df), c(-10, 150), col = 'red', lwd = 3)

hist(log10(pps_rh_ce), xlim = range(log10(c(pps_rh_ce, pps_rh_cc))), col = rgb(0, 1, 0, 0.3), main = 'Root Height', border = rgb(0, 1, 0, 0.3), xlab = paste0('Pcc=', pval_rh_cc, ' Pce=', pval_rh_ce, collapse = ''))
hist(log10(pps_rh_cc), add = T, col = rgb(0, 0, 1, 0.3), border = rgb(0, 0, 1, 0.3))
lines(log10(c(emp_rh, emp_rh)), c(-10, 150), col = 'red', lwd = 3)

hist(pps_stemmy_ce, xlim = range(c(pps_stemmy_ce, pps_stemmy_cc)), col = rgb(0, 1, 0, 0.3),
     main = 'Stemminess', border = rgb(0, 1, 0, 0.3), xlab = paste0('Pcc=', pval_stemmy_cc, ' Pce=', pval_stemmy_ce, collapse = ''))
hist(pps_stemmy_cc, add = T, col = rgb(0, 0, 1, 0.3), border = rgb(0, 0, 1, 0.3))
lines(c(emp_stemmy, emp_stemmy), c(-10, 100), col = 'red', lwd = 3)

hist(pps_colless_ce, xlim = range(c(pps_colless_cc, pps_colless_ce)), col = rgb(0, 1, 0, 0.3), main = 'Colless', border = rgb(0, 1, 0, 0.3), xlab = paste0('Pcc=', pval_colless_cc, ' Pce=', pval_colless_ce, collapse = ''))
hist(pps_colless_cc, add = T, col = rgb(0, 0, 1, 0.3), border = rgb(0, 0, 1, 0.3))
lines(c(emp_colless, emp_colless), c(-10, 200), lwd = 3, col = 'red')

hist(pps_sackin_ce, xlim = range(c(pps_sackin_ce, pps_sackin_cc)), col = rgb(0, 1, 0, 0.3),
     main = 'Sackin', border = rgb(0, 1, 0, 0.3), xlab = paste0('Pcc=', pval_sackin_cc, ' Pce=', pval_sackin_ce, collapse = ''))
hist(pps_sackin_cc, add = T, col = rgb(0, 0, 1, 0.3), border = rgb(0, 0, 1, 0.3))
lines(c(emp_sackin, emp_sackin), c(-10, 200), lwd = 3, col = 'red')

hist(log10(pps_time_lmax_ce), xlim = range(log10(c(pps_time_lmax_ce, pps_time_lmax_cc))), col = rgb(0, 1, 0, 0.3), main = 'Time of maximum lineages', border = rgb(0, 1, 0, 0.3), xlab = paste0('Pcc=', pval_time_lmax_cc, ' Pce=', pval_time_lmax_ce, collapse = ''))
hist(log10(pps_time_lmax_cc), add = T, col = rgb(0, 0, 1, 0.3), border = rgb(0, 0, 1, 0.3))
lines(log10(c(emp_time_lmax, emp_time_lmax)), c(-10, 200), lwd = 3, col = 'red')

hist(pps_slope_ratio_ce, xlim = range(pps_slope_ratio_ce, pps_slope_ratio_cc), col = rgb(0, 1, 0, 0.3), main = 'LTT slope ratio', border = rgb(0, 1, 0, 0.3), xlab = paste0('Pcc=', pval_slope_ratio_cc, ' Pce=', pval_slope_ratio_ce))
hist(pps_slope_ratio_cc, col = rgb(0, 0, 1, 0.3), add = T, border = rgb(0, 0, 1, 0.3))
lines(c(emp_slope_ratio, emp_slope_ratio), c(-10, 200), lwd = 3, col = 'red')

hist(pps_likelihood_ce, col = rgb(0, 1, 0, 0.3), main = 'Likelihood Exponential-growth', border = rgb(0, 1, 0, 0.3),
     xlab = paste0('Pce=', pval_likelihood_ce))
lines(c(emp_ce_likelihood, emp_ce_likelihood), c(-10, 25), col = 'red', lwd = 3)

hist(pps_likelihood_cc, col = rgb(0, 0, 1, 0.2), main = 'Likelihood Constant-size', border = rgb(0, 0, 1, 0.2),xlab = paste0('Pcc=', pval_likelihood_cc))
lines(c(emp_cc_likelihood, emp_cc_likelihood), c(-10, 25), col = 'red', lwd = 3)
dev.off()

cat(pval_df_cc, pval_df_ce, pval_rh_cc, pval_rh_ce, pval_stemmy_cc, pval_stemmy_ce, pval_colless_cc, pval_colless_ce, pval_sackin_cc, pval_sackin_ce, pval_time_lmax_cc, pval_time_lmax_ce, pval_slope_ratio_cc, pval_slope_ratio_ce, pval_likelihood_cc, pval_likelihood_ce, '\n', file = paste0(file_tag, '.txt'))
