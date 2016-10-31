## Run pipeline

# should have a folder for each set of trees (simulated under cc, ce, bd)

# split trees into each folder

# for each tree:
#  - make cc template
#  - run beast analysis
#  - make pps xml file
#  - run pps and collect trees
#  - make cc template for each pps tree
#  - run beast pps

source('pps_pipeline.R')
library(ape)
cc_trees <- read.tree('cc_veronika_tree_simulated.trees')

for(i in 1:length(cc_trees)){
  write.tree(cc_trees[[i]], paste0('cc_cc/cc_rep_', i, '.tree'))
  write.tree(cc_trees[[i]], paste0('cc_ce/cc_rep_', i, '.tree'))
  write.tree(cc_trees[[i]], paste0('cc_bd/cc_rep_', i, '.tree'))
}

cc_files <- c('cc_cc', 'cc_ce', 'cc_bd')

trees_lines <- readLines('cc_veronika_tree_simulated.trees')

setwd('cc_cc')
xml_files <- dir(pattern = 'xml$')

  tr <- 'cc_rep_1.tree' # INPUT_Tree

  bd_run <- function(tr, beast_command){
    file_name <- gsub('[.]tree', '_bd', tr)
    xml_file <- paste0(file_name, '_1.xml')
    make_bd_template(readLines(tr), file_name)
    log_temp <- run_beast_analyses(beast_command, xml_file)
    make_bd_simulation(log_temp, read.tree('dated.tree'), file_name)
    pps <- run_beast_simulation(beast_command, paste0(file_name, '_pps.xml'))
    pps_names <- paste0(file_name, '_pps')
    make_bd_template(write.tree(pps), pps_names)
    pps_files <- paste0(pps_names, '_', 1:length(pps), '.xml')
    pps_results <- run_beast_pps(beast_command, pps_files)
    write.table(pps_results, file = paste0(file_name, '_pps_result.txt'), row.names = F, quote = F)
  }

setwd('..')
