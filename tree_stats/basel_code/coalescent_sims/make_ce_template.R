# Input stings are:

# INPUT_TREE_STRING
# INPUT_TAXON_SETS, which corresponds to <taxon id="t100_2001.37" spec="Taxon"/>
# POSTERIOR_OUTPUT_FILE
library(ape)
ce_template <- readLines("ce_estimate_template.xml")
#trees_lines <- readLines('1_n100_lambda0.55_delta0.5_p1_N10000_Imodel_CE.newick')
#trees_lines <- readLines('cc_sims_subsample.trees')
#trees_lines <- readLines('cc_sim_rep_1_pps.trees')
#trees_lines <- readLines('../coal_veronika_subset.trees')
trees_lines <- readLines('ce_veronika_ce_pps.trees')

#i <- 1
for(i in 1:100){
      tr_temp <- read.tree(text = trees_lines[i])
      taxon_sets <- paste0("<taxon id=\"", tr_temp$tip.label, "\" spec=\"Taxon\"/>", collapse = "\n")
      rep_name <- paste0('ce_veronika_pps_', i)
      xml_temp <- gsub('POSTERIOR_OUTPUT_FILE', rep_name, gsub("INPUT_TAXON_SETS", taxon_sets, gsub("INPUT_TREE_STRING", trees_lines[i], ce_template)))
      cat(xml_temp, file = paste0(rep_name, '.xml'), sep = '\n')
}

