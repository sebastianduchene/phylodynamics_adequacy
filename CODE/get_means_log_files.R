
library(NELSI)
library(apTreeshape)

get_mean_logfiles <- function(file_names){

    means_list <- list()
    for(i in 1:length(file_names)){
        dat_temp <- read.table(file_names[i], head = T)
        dat_temp <- dat_temp[-(1:ceiling(nrow(dat_temp)*0.1)), ]
        means_list[[i]] <- colMeans(dat_temp)
    }
    rbind.list <- function(vlist){
        if(length(vlist) == 1){
            return(vlist[[1]])
        }else if(length(vlist) == 2){
            return(rbind(vlist[[1]], vlist[[2]]))
        }else{
            return(rbind(vlist[[1]], rbind.list(vlist[-1])))
        }
    }

    return(rbind.list(means_list))
}



summarise_pps <- function(empirical_tree, empirical_log, pps_trees, pps_logs, outfile = 'outfile'){

    cols <- c(rgb(0, 1, 0.4, 0.2), rgb(0, 0, 1, 0.2), rgb(1, 0.5, 0, 0.2))

    pdf(paste0(outfile, '_stats.pdf'), useDingbats = F, width = 6, height = 3)
    par(mfrow = c(1, 3))
    df_pps <- list()
    rh_pps <- list()
    for(i in 1:length(pps_trees)){
        df_pps[[i]] <- sapply(pps_trees[[i]], function(tr) get.df(tr))
        rh_pps[[i]] <- log10(sapply(pps_trees[[i]], function(tr) max(intnode.times(tr))))
    }
    names(df_pps) <- names(pps_trees)
    names(rh_pps) <- names(pps_trees)

    emp_df <- get.df(empirical_tree)
    emp_rh <- log10(max(intnode.times(empirical_tree)))
    plot(x = range(unlist(df_pps)), y = range(unlist(rh_pps)), type = 'n', xlab = 'DF', ylab = 'RH')
    for(i in 1:length(df_pps)){
        points(df_pps[[i]], rh_pps[[i]], col = cols[i], pch = 20)
    }
    points(emp_df, emp_rh, pch = 15, col = 'red', cex = 1)

    emp_ltt <- get.ltt.summary(empirical_tree)
    pps_maxLtt <- list()
    pps_slopeLtt <- list()
    for(i in 1:length(pps_trees)){
        ltti <- t(sapply(pps_trees[[i]], function(tr) get.ltt.summary(tr)))
        pps_maxLtt[[i]] <- ltti[, 'time_max_lineages']
        pps_slopeLtt[[i]] <- ltti[, 'ratio_slopes']
    }

    plot(x = range(unlist(pps_maxLtt)), y = range(unlist(pps_slopeLtt)), type = 'n', xlab = 'max Ltt', ylab = 'Ltt slope ratio')
    for(i in 1:length(pps_maxLtt)){
        points(pps_maxLtt[[i]], pps_slopeLtt[[i]], pch = 20, col = cols[i])
    }
    points(emp_ltt['time_max_lineages'], emp_ltt['ratio_slopes'], pch = 15, col = 'red', cex = 1)

    emp_colless <- colless(as.treeshape(empirical_tree))
    emp_sackin <- sackin(as.treeshape(empirical_tree))

    pps_colless <- list()
    pps_sackin <- list()
    for(i in 1:length(pps_trees)){
        pps_colless[[i]] <- sapply(pps_trees[[i]], function(tr) colless(as.treeshape(tr)))
        pps_sackin[[i]] <- sapply(pps_trees[[i]], function(tr) sackin(as.treeshape(tr)))
    }

    plot(range(unlist(pps_colless)), range(unlist(pps_sackin)), type = 'n', xlab = 'Colless', ylab = 'Sackin')
    for(i in 1:length(pps_trees)){
        points(pps_colless[[i]], pps_sackin[[i]], pch = 20, col = cols[i])
    }
    points(emp_colless, emp_sackin, pch = 15, col = 'red', cex = 1)
    dev.off()

    pdf(paste0(outfile, '_lik.pdf'), useDingbats = F, width = 5, height = 5)
    par(mfcol = c(length(pps_logs), 3))
    for(i in 1:length(pps_logs)){
        hist(pps_logs[[i]][, 'likelihood'], main = '', xlab = names(pps_logs)[i], col = cols[i], border = cols[i])
        lines(rep(mean(empirical_log[[i]][-(1:100), 'likelihood']), 2), c(0, 20), col = 'red', lwd = 2)
    }
    for(i in 1:length(pps_trees)){
        sample_tree <- pps_trees[[i]][sample(x = 100:length(pps_trees[[i]]), size = 1)][[1]]
        root_time <- max(intnode.times(sample_tree))
        plot(sample_tree, show.tip.label = F)
        nodelabels(round(root_time, 1), node = length(sample_tree$tip.label)+1, cex = 0.4, bg = 'white', frame = 'circ')
    }
    plot(empirical_tree, show.tip.label = F)
    nodelabels(round(max(intnode.times(empirical_tree)), 1), node = length(empirical_tree$tip.label)+1, cex = 0.4, bg = 'white', frame = 'circ')
    dev.off()

    make_two_sided <- function(x) 0.5 - abs(x - 0.5)

    pvals_df <- make_two_sided(sapply(1:length(df_pps),
                                      function(i) sum(emp_df > df_pps[[i]]) / length(df_pps[[i]])))

    pvals_rh <- make_two_sided(sapply(1:length(rh_pps),
                                      function(i) sum(emp_rh > rh_pps[[i]]) / length(rh_pps[[i]])))

    pvals_maxLtt <- make_two_sided(sapply(1:length(pps_maxLtt),
                                          function(i) sum(emp_ltt['time_max_lineages'] > pps_maxLtt[[i]]) / length(pps_maxLtt[[i]])))

    pvals_slopeLtt <- make_two_sided(sapply(1:length(pps_slopeLtt),
                                            function(i) sum(emp_ltt['ratio_slopes'] > pps_slopeLtt[[i]]) / length(pps_slopeLtt[[i]])))

    pvals_colless <- make_two_sided(sapply(1:length(pps_colless),
                                           function(i) sum(emp_colless > pps_colless[[i]]) / length(pps_colless[[i]])))

    pvals_sackin <- make_two_sided(sapply(1:length(pps_sackin),
                   function(i) sum(emp_sackin > pps_sackin[[i]]) / length(pps_sackin[[i]])))

    lik_pvals <- vector()
    for(i in 1:length(pps_logs)){
        lik_pvals[i] <- make_two_sided(sum(mean(empirical_log[[i]][-(1:100), 'likelihood']) > pps_logs[[i]][, 'likelihood'])/nrow(pps_logs[[i]]))
    }

    output_table <- cbind(pvals_df, pvals_rh, pvals_maxLtt, pvals_slopeLtt, pvals_colless, pvals_sackin, lik_pvals)
    rownames(output_table) <- names(pps_trees)
    write.table(output_table, file = paste0(outfile, '.txt'))
    return(output_table)
}





# To test:
"
ce_log_files <- paste0('../test_data/', dir('../test_data/', pattern = '50taxa.+_ce_pps_.+log'))
cc_log_files <- paste0('../test_data/', dir('../test_data/', pattern = '50taxa.+_cc_pps_.+log'))
bd_log_files <- paste0('../test_data/', dir('../test_data/', pattern = '50taxa.+_bd_pps_.+log'))


ce_means <- get_mean_logfiles(ce_log_files)
cc_means <- get_mean_logfiles(cc_log_files)
bd_means <- get_mean_logfiles(bd_log_files)

par(mfrow = c(1, 2))
hist(ce_means[, 'likelihood'])
hist(cc_means[, 'likelihood'])
hist(bd_means[, 'likelihood'])

empirical_tree <- read.tree('../test_data/50taxa_CE_1e10_08_rep_81.tree')

empirical_log <- list(cc = read.table('../test_data/50taxa_CE_1e10_08_rep_81_cc_1.log', head = T), ce = read.table('../test_data/50taxa_CE_1e10_08_rep_81_ce_1.log', head = T),
                      bd = read.table('../test_data/50taxa_CE_1e10_08_rep_81_bd_1.log', head = T))

pps_trees <- list(cc = read.nexus('../test_data/50taxa_CE_1e10_08_rep_81_cc_pps.trees')[-(1:100)], ce = read.nexus('../test_data/50taxa_CE_1e10_08_rep_81_ce_pps.trees')[-(1:100)],
                  bd = read.nexus('../test_data/50taxa_CE_1e10_08_rep_81_bd_pps.trees')[-(1:100)])

pps_logs <- list(cc = cc_means, ce = ce_means, bd = bd_means)
"
