library(ape)

args = commandArgs(trailingOnly = T)

print('to use: Rscript make_cc_simulation_distro.R xml_template log_file tree_file')

# ENABLE ARGS AND TEST

xml_file <- args[1]
log_file <- args[2]

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

xml_temp <- gsub('BD_SIM_TREE_FILE', gsub('[.]log', '_pps', log_file), xml_temp)

origin <- round(c(mean(posterior_params_file$origin), sd(posterior_params_file$origin)), 2)
sampling <- round(c(mean(posterior_params_file$samplingProportion),
                    sd(posterior_params_file$samplingProportion)), 2)
becomeUninfect <- round(c(mean(posterior_params_file$becomeUninfectiousRate),
                          sd(posterior_params_file$becomeUninfectiousRate)), 2)

r0s <- posterior_params_file[, grep('R0', colnames(posterior_params_file))]
r0s_means <- round(colMeans(r0s), 2)
r0s_sds <- round(sapply(1:ncol(r0s), function(x) sd(r0s[, x])), 2)

xml_temp <- gsub('ORIGIN_MEAN', origin[1], xml_temp)
xml_temp <- gsub('ORIGIN_SD', origin[2], xml_temp)

xml_temp <- gsub('BECOME_UNINFECTIOUS_MEAN', becomeUninfect[1], xml_temp)
xml_temp <- gsub('BECOME_UNINFECTIOUS_SD', becomeUninfect[2], xml_temp)

xml_temp <- gsub('SAMPLING_PROP_MEAN', sampling[1], xml_temp)
xml_temp <- gsub('SAMPLING_PROP_SD', sampling[2], xml_temp)

xml_temp <- gsub('R0_MEAN', paste0(r0s_means, collapse = ' '), xml_temp)
xml_temp <- gsub('R0_SD', paste0(r0s_sds, collapse = ' '), xml_temp)

cat(xml_temp, file = gsub('[.]log', '_pps.xml', log_file), sep = '\n')

