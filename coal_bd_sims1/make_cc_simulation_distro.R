# need to get log file, tree, and template

library(ape)

args = commandArgs(trailingOnly = T)

print('to use: Rscript make_cc_simulation_distro.R xml_template log_file tree_file')

xml_file <- args[1]
log_file <- args[2]
#xml_file <- 'cc_sim_template.xml'
#log_file <- 'cc_replicate_24.log'

xml_simulation_template <- readLines(xml_file)
posterior_params_file <- read.table(log_file, head = T)
input_tree <- read.tree(args[3])

# Make taxon_seqs:
taxon_seqs <- paste0("<sequence id=\"seq_", input_tree$tip.label, "\" taxon=\"",input_tree$tip.label, "\" totalcount=\"4\" value=\"gc\"/>", collapse = '\n')

taxon_dates <- vector()
for(ta in 1:length(input_tree$tip.label)){
       date <- gsub('.+_', '', input_tree$tip.label[ta])
       taxon_dates[ta] <- paste0(input_tree$tip.label[ta], '=', date)
}
taxon_dates <- paste0(taxon_dates, collapse = ',\n')

xml_temp <- gsub('TAXON_DATES', taxon_dates, gsub('TAXON_SEQS', taxon_seqs, xml_simulation_template))
xml_temp <- gsub('CC_SIM_TREE_FILE', gsub('[.]log', '_pps', log_file), xml_temp)

#Get posterior params to set as the simulation prior:
popsize_normal <- round(c(mean(posterior_params_file$popSize), sd(posterior_params_file$popSize)), 2)

xml_temp <- gsub('CONSTANT_POPSIZE_MEAN', popsize_normal[1], xml_temp)
xml_temp <- gsub('CONSTANT_POPSIZE_SD', popsize_normal[2], xml_temp)

cat(xml_temp, file = gsub('[.]log', '_pps.xml', log_file), sep = '\n')

