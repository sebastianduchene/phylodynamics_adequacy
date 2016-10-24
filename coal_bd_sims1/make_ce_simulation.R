library(ape)
# Input a tree to get tip labels and tip heights. This can correspond to the empirical tree.

# Tags: 
# TAXON_SEQS: <sequence id="seq_t100_2001.37" taxon="t100_2001.37" totalcount="4" value="gc"/>
# TAXON_DATES: t1_2000=2000.0,
# E_POP_SIZE
# GROWTH_RATE

xml_simulation_template <- readLines('ce_sim_template.xml')
#input_tree <- read.tree('../cc_sims_subsample.trees')[[1]]
input_tree <- read.tree('dated.tree')


## 
posterior_params_file <- "dated_ce_1.log"
posterior <- read.table(posterior_params_file, head = T)

# Make taxon_seqs:
taxon_seqs <- paste0("<sequence id=\"seq_", input_tree$tip.label, "\" taxon=\"",input_tree$tip.label, "\" totalcount=\"4\" value=\"gc\"/>", collapse = '\n')

taxon_dates <- vector()
for(ta in 1:length(input_tree$tip.label)){
       date <- gsub('.+_', '', input_tree$tip.label[ta])
       taxon_dates[ta] <- paste0(input_tree$tip.label[ta], '=', date)
}
taxon_dates <- paste0(taxon_dates, collapse = ',\n')

xml_temp <- gsub('TAXON_DATES', taxon_dates, gsub('TAXON_SEQS', taxon_seqs, xml_simulation_template))

xml_temp <- gsub('CE_SIM_TREE_FILE', 'ce_predictive_sims', xml_temp)

# Take some 100 random draws from the posterior. Each will generate an xml file. Kill all operators for parameters. Set starting values sample a bit and get a tree.

posterior_samples <- sample(100:nrow(posterior), 100)

simulated_trees <- list()
for(i in 1:length(posterior_samples)){
      r_xml_temp <- gsub('E_POP_SIZE', posterior$ePopSize[posterior_samples[i]], xml_temp)
      r_xml_temp <- gsub('GROWTH_RATE', posterior$growthRate[posterior_samples[i]], r_xml_temp)
      cat(r_xml_temp, file = paste0('ce_get_pps.xml', i), sep = '\n')
      system(paste0('java -jar ~/Desktop/phylo_programs/beast243/lib/beast.jar -overwrite ce_get_pps.xml', i))
      pps_trees <- read.nexus('ce_predictive_sims.trees') 
      simulated_trees[[i]] <- pps_trees[length(pps_trees)]
}

sim_trees <- lapply(simulated_trees, function(x) x[[1]])
class(sim_trees) <- 'multiPhylo'
write.tree(sim_trees, file = gsub('[.]log', '_ce_pps.trees', posterior_params_file))
