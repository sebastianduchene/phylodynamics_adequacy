# Input stings are:

# INPUT_TREE_STRING
# INPUT_TAXON_SETS, which corresponds to <taxon id="t100_2001.37" spec="Taxon"/>
# POSTERIOR_OUTPUT_FILE
library(ape)


#bd_template <- readLines("bd_estimate_template.xml")
#trees_lines <- readLines('bd_sim_100.trees')

bd_template <- readLines("../bd_estimate_template.xml")
trees_lines <- readLines('random_trees_bd_1_pps.trees')


# All trees sholud have the same taxa
tr1 <- read.tree(text = trees_lines[[1]])


for(i in 1:100){
      tr_temp <- read.tree(text = trees_lines[i])
#      rep_name <- paste0('bd_replicate_', i)
      taxon_sets <- paste0("<taxon id=\"", tr_temp$tip.label, "\" spec=\"Taxon\"/>", collapse = "\n")
      rep_name <- paste0('random_pps1_', i)
      xml_temp <- gsub('POSTERIOR_OUTPUT_FILE', rep_name, gsub("INPUT_TAXON_SETS", taxon_sets, gsub("INPUT_TREE_STRING", trees_lines[i], bd_template)))
      cat(xml_temp, file = paste0(rep_name, '.xml'), sep = '\n')
}



