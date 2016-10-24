# Input settings are:

args = commandArgs(trailingOnly=TRUE)

cc_template <- readLines(args[1])
trees_lines <- readLines(args[2])

# INPUT_TREE_STRING
# INPUT_TAXON_SETS, which corresponds to <taxon id="t100_2001.37" spec="Taxon"/>
# POSTERIOR_OUTPUT_FILE
library(ape)
#cc_template <- readLines("CC_estimate_template.xml")
#trees_lines <- readLines('cc_veronika_tree_simulated.trees')

for(i in 1:length(trees_lines)){
      tr_temp <- read.tree(text = trees_lines[i])
      taxon_sets <- paste0("<taxon id=\"", tr_temp$tip.label, "\" spec=\"Taxon\"/>", collapse = "\n")
      rep_name <- paste0(gsub('[.]trees?', '', args[2]), '_cc_', i)
      xml_temp <- gsub('POSTERIOR_OUTPUT_FILE', rep_name, gsub("INPUT_TAXON_SETS", taxon_sets, gsub("INPUT_TREE_STRING", trees_lines[i], cc_template)))
      cat(xml_temp, file = paste0(rep_name, '.xml'), sep = '\n')
}

