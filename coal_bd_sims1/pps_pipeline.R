## Make

# - read tree and select an analysis template
#

# Input settings are:

library(ape)

trees_lines <- readLines('dated.tree')
rep_name <- 'rep_test'

make_cc_template <- function(trees_lines, rep_name){

  cc_template <- "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?><beast beautitemplate='Standard' beautistatus='' namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" required=\"\" version=\"2.4\">
  <map name=\"Uniform\" >beast.math.distributions.Uniform</map>
  <map name=\"Exponential\" >beast.math.distributions.Exponential</map>
  <map name=\"LogNormal\" >beast.math.distributions.LogNormalDistributionModel</map>
  <map name=\"Normal\" >beast.math.distributions.Normal</map>
  <map name=\"Beta\" >beast.math.distributions.Beta</map>
  <map name=\"Gamma\" >beast.math.distributions.Gamma</map>
  <map name=\"LaplaceDistribution\" >beast.math.distributions.LaplaceDistribution</map>
  <map name=\"prior\" >beast.math.distributions.Prior</map>
  <map name=\"InverseGamma\" >beast.math.distributions.InverseGamma</map>
  <map name=\"OneOnX\" >beast.math.distributions.OneOnX</map>

  <tree id=\"Tree.t:dummy_aln\" name=\"stateNode\" spec=\"beast.util.TreeParser\"
       IsLabelledNewick=\"true\" adjustTipHeights=\"false\"
       newick=\"INPUT_TREE_STRING\">
           <taxonset id=\"TaxonSet.dummy_aln\" spec=\"TaxonSet\">
	     INPUT_TAXON_SETS
           </taxonset>
  </tree>


  <run id=\"mcmc\" spec=\"MCMC\" chainLength=\"5000000\" sampleFromPrior=\"false\">
    <state id=\"state\" storeEvery=\"5000\">
        <parameter id=\"popSize.t:dummy_aln\" name=\"stateNode\">17.0</parameter>
    </state>

    <distribution id=\"posterior\" spec=\"util.CompoundDistribution\">
        <distribution id=\"prior\" spec=\"util.CompoundDistribution\">
            <prior id=\"PopSizePrior.t:dummy_aln\" name=\"distribution\" x=\"@popSize.t:dummy_aln\">
	      <!--Very uninformative prior on population size -->
	      <Normal id=\"Normal.01\" name=\"distr\">
		<parameter id=\"RealParameter.0\" estimate=\"false\" name=\"mean\">100</parameter>
		<parameter id=\"RealParameter.01\" estimate=\"false\" name=\"sigma\">100</parameter>
	      </Normal>
            </prior>
        </distribution>
        <distribution id=\"likelihood\" spec=\"util.CompoundDistribution\" useThreads=\"true\">
            <distribution id=\"CoalescentConstant.t:dummy_aln\" spec=\"Coalescent\">
                <populationModel id=\"ConstantPopulation.t:dummy_aln\" spec=\"ConstantPopulation\" popSize=\"@popSize.t:dummy_aln\"/>
                <treeIntervals id=\"TreeIntervals.t:dummy_aln\" spec=\"TreeIntervals\" tree=\"@Tree.t:dummy_aln\"/>
            </distribution>

        </distribution>
    </distribution>

    <operator id=\"PopSizeScaler.t:dummy_aln\" spec=\"ScaleOperator\" parameter=\"@popSize.t:dummy_aln\" scaleFactor=\"0.75\" weight=\"3.0\"/>
    <logger id=\"tracelog\" fileName=\"POSTERIOR_OUTPUT_FILE.log\" logEvery=\"10000\" model=\"@posterior\" sanitiseHeaders=\"true\" sort=\"smart\">
        <log idref=\"posterior\"/>
        <log idref=\"likelihood\"/>
        <log idref=\"prior\"/>
        <log id=\"TreeHeight.t:dummy_aln\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@Tree.t:dummy_aln\"/>
        <log idref=\"popSize.t:dummy_aln\"/>
        <log idref=\"CoalescentConstant.t:dummy_aln\"/>
    </logger>

    <logger id=\"screenlog\" fileName=\"\" logEvery=\"10000\">
        <log idref=\"posterior\"/>
        <log id=\"ESS.0\" spec=\"util.ESS\" arg=\"@posterior\"/>
        <log idref=\"likelihood\"/>
        <log idref=\"prior\"/>
    </logger>
    </run>
    </beast>
    "

    for(i in 1:length(trees_lines)){
      tr_temp <- read.tree(text = trees_lines[i])
      taxon_sets <- paste0("<taxon id=\"", tr_temp$tip.label, "\" spec=\"Taxon\"/>", collapse = "\n")
      xml_temp <- gsub('POSTERIOR_OUTPUT_FILE', rep_name, gsub("INPUT_TAXON_SETS", taxon_sets, gsub("INPUT_TREE_STRING", trees_lines[i], cc_template)))
      cat(xml_temp, file = paste0(rep_name, '.xml'), sep = '\n')
    }
}

make_ce_template <- function(trees_lines, rep_name){
  ce_template <- "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?><beast beautitemplate='Standard' beautistatus='' namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" required=\"\" version=\"2.4\">

  <map name=\"Uniform\" >beast.math.distributions.Uniform</map>
  <map name=\"Exponential\" >beast.math.distributions.Exponential</map>
  <map name=\"LogNormal\" >beast.math.distributions.LogNormalDistributionModel</map>
  <map name=\"Normal\" >beast.math.distributions.Normal</map>
  <map name=\"Beta\" >beast.math.distributions.Beta</map>
  <map name=\"Gamma\" >beast.math.distributions.Gamma</map>
  <map name=\"LaplaceDistribution\" >beast.math.distributions.LaplaceDistribution</map>
  <map name=\"prior\" >beast.math.distributions.Prior</map>
  <map name=\"InverseGamma\" >beast.math.distributions.InverseGamma</map>
  <map name=\"OneOnX\" >beast.math.distributions.OneOnX</map>

  <tree id=\"Tree.t:dummy_aln\" name=\"stateNode\" spec=\"beast.util.TreeParser\"
        IsLabelledNewick=\"true\" adjustTipHeights=\"false\"
        newick=\"INPUT_TREE_STRING\">
    <taxonset id=\"TaxonSet.dummy_aln\" spec=\"TaxonSet\">
      INPUT_TAXON_SETS
    </taxonset>
  </tree>


  <run id=\"mcmc\" spec=\"MCMC\" chainLength=\"20000000\" sampleFromPrior=\"false\">
      <state id=\"state\" storeEvery=\"5000\">
        <parameter id=\"ePopSize.t:dummy_aln\" name=\"stateNode\">22.8</parameter>
        <parameter id=\"growthRate.t:dummy_aln\" name=\"stateNode\">0.028</parameter>
      </state>

      <distribution id=\"posterior\" spec=\"util.CompoundDistribution\">
          <distribution id=\"prior\" spec=\"util.CompoundDistribution\">
              <prior id=\"ePopSizePrior.t:dummy_aln\" name=\"distribution\" x=\"@ePopSize.t:dummy_aln\">
                <Normal id=\"Normal.01\" name=\"distr\">
                  <parameter id=\"RealParameter.0\" estimate=\"false\" name=\"mean\">22.8</parameter>
                  <parameter id=\"RealParameter.01\" estimate=\"false\" name=\"sigma\">15.0</parameter>
                </Normal>
              </prior>
              <prior id=\"GrowthRatePrior.t:dummy_aln\" name=\"distribution\" x=\"@growthRate.t:dummy_aln\">
                <Normal id=\"Normal.02\" name=\"distr\">
                  <parameter id=\"RealParameter.02\" estimate=\"false\" name=\"mean\">0.028</parameter>
                  <parameter id=\"RealParameter.03\" estimate=\"false\" name=\"sigma\">5</parameter>
                </Normal>
              </prior>
          </distribution>
          <distribution id=\"likelihood\" spec=\"util.CompoundDistribution\" useThreads=\"true\">
  	  <distribution id=\"CoalescentExponential.t:dummy_aln\" spec=\"Coalescent\">
              <populationModel id=\"ExponentialGrowth.t:dummy_aln\" spec=\"ExponentialGrowth\" growthRate=\"@growthRate.t:dummy_aln\" popSize=\"@ePopSize.t:dummy_aln\"/>
              <treeIntervals id=\"TreeIntervals.t:dummy_aln\" spec=\"TreeIntervals\" tree=\"@Tree.t:dummy_aln\"/>
            </distribution>
          </distribution>
      </distribution>
  <!-- Kill topology operators to for analysis -->
      <operator id=\"ePopSizeScaler.t:dummy_aln\" spec=\"ScaleOperator\" parameter=\"@ePopSize.t:dummy_aln\" scaleFactor=\"0.75\" weight=\"3.0\"/>
      <operator id=\"GrowthRateRandomWalk.t:dummy_aln\" spec=\"RealRandomWalkOperator\" parameter=\"@growthRate.t:dummy_aln\" weight=\"3.0\" windowSize=\"1.0\"/>

      <logger id=\"tracelog\" fileName=\"POSTERIOR_OUTPUT_FILE.log\" logEvery=\"10000\" model=\"@posterior\" sanitiseHeaders=\"true\" sort=\"smart\">
          <log idref=\"posterior\"/>
          <log idref=\"likelihood\"/>
          <log idref=\"prior\"/>
          <log id=\"TreeHeight.t:dummy_aln\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@Tree.t:dummy_aln\"/>
          <log idref=\"CoalescentExponential.t:dummy_aln\"/>
          <log idref=\"ePopSize.t:dummy_aln\"/>
          <log idref=\"growthRate.t:dummy_aln\"/>
      </logger>

      <logger id=\"screenlog\" fileName=\"\" logEvery=\"10000\">
          <log idref=\"posterior\"/>
          <log id=\"ESS.0\" spec=\"util.ESS\" arg=\"@posterior\"/>
          <log idref=\"likelihood\"/>
          <log idref=\"prior\"/>
      </logger>


  </run>

  </beast>"
  for(i in 1:length(trees_lines)){
    tr_temp <- read.tree(text = trees_lines[i])
    taxon_sets <- paste0("<taxon id=\"", tr_temp$tip.label, "\" spec=\"Taxon\"/>", collapse = "\n")
    xml_temp <- gsub('POSTERIOR_OUTPUT_FILE', rep_name, gsub("INPUT_TAXON_SETS", taxon_sets, gsub("INPUT_TREE_STRING", trees_lines[i], ce_template)))
    cat(xml_temp, file = paste0(rep_name, '.xml'), sep = '\n')
  }

}

make_bd_template <- function(trees_lines, rep_name){
  bd_template <- "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?><beast beautitemplate='Standard' beautistatus='' namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" required=\"BDSKY v1.3.2\" version=\"2.4\">

  <map name=\"Uniform\" >beast.math.distributions.Uniform</map>
  <map name=\"Exponential\" >beast.math.distributions.Exponential</map>
  <map name=\"LogNormal\" >beast.math.distributions.LogNormalDistributionModel</map>
  <map name=\"Normal\" >beast.math.distributions.Normal</map>
  <map name=\"Beta\" >beast.math.distributions.Beta</map>
  <map name=\"Gamma\" >beast.math.distributions.Gamma</map>
  <map name=\"LaplaceDistribution\" >beast.math.distributions.LaplaceDistribution</map>
  <map name=\"prior\" >beast.math.distributions.Prior</map>
  <map name=\"InverseGamma\" >beast.math.distributions.InverseGamma</map>
  <map name=\"OneOnX\" >beast.math.distributions.OneOnX</map>

         <tree id=\"Tree.t:dummy_aln\" name=\"stateNode\" spec=\"beast.util.TreeParser\"
         IsLabelledNewick=\"true\" adjustTipHeights=\"false\"
         newick=\"INPUT_TREE_STRING\">
             <taxonset id=\"TaxonSet.dummy_aln\" spec=\"TaxonSet\">
  	     INPUT_TAXON_SETS
             </taxonset>
         </tree>

  <run id=\"mcmc\" spec=\"MCMC\" chainLength=\"10000000\" sampleFromPrior=\"false\">
      <state id=\"state\" storeEvery=\"5000\">
          <parameter id=\"origin.s.t:dummy_aln\" lower=\"0.0\" name=\"stateNode\" upper=\"Infinity\">100.0</parameter>
          <parameter id=\"samplingProportion.s.t:dummy_aln\" lower=\"0.0\" name=\"stateNode\" upper=\"1.0\">0.01</parameter>
          <parameter id=\"becomeUninfectiousRate.s.t:dummy_aln\" lower=\"0.0\" name=\"stateNode\" upper=\"Infinity\">4.0</parameter>
          <parameter id=\"R0.s.t:dummy_aln\" dimension=\"10\" lower=\"0.0\" name=\"stateNode\" upper=\"Infinity\">3</parameter>
      </state>

      <distribution id=\"posterior\" spec=\"util.CompoundDistribution\">
          <distribution id=\"prior\" spec=\"util.CompoundDistribution\">
              <prior id=\"RPrior.s.t:dummy_aln\" name=\"distribution\" x=\"@R0.s.t:dummy_aln\">
                  <Normal id=\"Normal.0\" name=\"distr\">
                      <parameter id=\"RealParameter.0\" estimate=\"false\" name=\"mean\">3</parameter>
                      <parameter id=\"RealParameter.01\" estimate=\"false\" name=\"sigma\">10</parameter>
                  </Normal>
              </prior>
              <prior id=\"becomeUninfectiousRatePrior.s.t:dummy_aln\" name=\"distribution\" x=\"@becomeUninfectiousRate.s.t:dummy_aln\">
                  <Normal id=\"Normal.01\" name=\"distr\">
                      <parameter id=\"RealParameter.02\" estimate=\"false\" name=\"mean\">2</parameter>
                      <parameter id=\"RealParameter.03\" estimate=\"false\" name=\"sigma\">10</parameter>
                  </Normal>
              </prior>
              <prior id=\"originPrior.s.t:dummy_aln\" name=\"distribution\" x=\"@origin.s.t:dummy_aln\">
                <Normal id=\"Normal.03\" name=\"distr\">
                  <parameter id=\"RealParameter.06\" estimate=\"false\" name=\"mean\">50</parameter>
                  <parameter id=\"RealParameter.07\" estimate=\"false\" name=\"sigma\">20</parameter>
                </Normal>
              </prior>
              <prior id=\"samplingProportionPrior.s.t:dummy_aln\" name=\"distribution\" x=\"@samplingProportion.s.t:dummy_aln\">
                  <Normal id=\"Normal.02\" name=\"distr\">
                      <parameter id=\"RealParameter.04\" estimate=\"false\" name=\"mean\">0.01</parameter>
                      <parameter id=\"RealParameter.05\" estimate=\"false\" name=\"sigma\">1.0E-5</parameter>
                  </Normal>
              </prior>
          </distribution>
          <distribution id=\"likelihood\" spec=\"util.CompoundDistribution\" useThreads=\"true\">
            <distribution id=\"BirthDeathSkySerial.t:dummy_aln\" spec=\"beast.evolution.speciation.BirthDeathSkylineModel\" R0=\"@R0.s.t:dummy_aln\" becomeUninfectiousRate=\"@becomeUninfectiousRate.s.t:dummy_aln\" origin=\"@origin.s.t:dummy_aln\" samplingProportion=\"@samplingProportion.s.t:dummy_aln\" tree=\"@Tree.t:dummy_aln\">
            </distribution>
          </distribution>
      </distribution>

      <operator id=\"becomeUninfectiousRateScaler.s.t:dummy_aln\" spec=\"ScaleOperator\" parameter=\"@becomeUninfectiousRate.s.t:dummy_aln\" scaleFactor=\"0.75\" weight=\"2.0\"/>
      <!-- Note that the sampling proportion is fixed -->
  <!--
      <operator id=\"samplingScaler.s.t:dummy_aln\" spec=\"ScaleOperator\" parameter=\"@samplingProportion.s.t:dummy_aln\" scaleFactor=\"0.75\" weight=\"2.0\"/>
  -->
      <operator id=\"RScaler.s.t:dummy_aln\" spec=\"ScaleOperator\" parameter=\"@R0.s.t:dummy_aln\" scaleFactor=\"0.75\" weight=\"10.0\"/>

      <operator id=\"updownBD.s.t:dummy_aln\" spec=\"UpDownOperator\" scaleFactor=\"0.75\" weight=\"2.0\">
          <up idref=\"R0.s.t:dummy_aln\"/>
          <down idref=\"becomeUninfectiousRate.s.t:dummy_aln\"/>
      </operator>

      <operator id=\"origScaler.s.t:dummy_aln\" spec=\"ScaleOperator\" parameter=\"@origin.s.t:dummy_aln\" scaleFactor=\"0.75\" weight=\"1.0\"/>

      <logger id=\"tracelog\" fileName=\"POSTERIOR_OUTPUT_FILE.log\" logEvery=\"20000\" model=\"@posterior\" sanitiseHeaders=\"true\" sort=\"smart\">
          <log idref=\"posterior\"/>
          <log idref=\"likelihood\"/>
          <log idref=\"prior\"/>
          <log idref=\"BirthDeathSkySerial.t:dummy_aln\"/>
          <log idref=\"origin.s.t:dummy_aln\"/>
          <log idref=\"samplingProportion.s.t:dummy_aln\"/>
          <log idref=\"becomeUninfectiousRate.s.t:dummy_aln\"/>
          <log idref=\"R0.s.t:dummy_aln\"/>
          <log id=\"birth.t:dummy_aln\" spec=\"beast.math.statistic.RPNcalculator\" expression=\"R0.s.t:dummy_aln becomeUninfectiousRate.s.t:dummy_aln *\">
              <parameter idref=\"becomeUninfectiousRate.s.t:dummy_aln\"/>
              <parameter idref=\"R0.s.t:dummy_aln\"/>
          </log>
          <log id=\"death.t:dummy_aln\" spec=\"beast.math.statistic.RPNcalculator\" expression=\"becomeUninfectiousRate.s.t:dummy_aln 1 samplingProportion.s.t:dummy_aln - *\">
              <parameter idref=\"becomeUninfectiousRate.s.t:dummy_aln\"/>
              <parameter idref=\"samplingProportion.s.t:dummy_aln\"/>
          </log>
          <log id=\"sampling.t:dummy_aln\" spec=\"beast.math.statistic.RPNcalculator\" expression=\"becomeUninfectiousRate.s.t:dummy_aln samplingProportion.s.t:dummy_aln *\">
              <parameter idref=\"becomeUninfectiousRate.s.t:dummy_aln\"/>
              <parameter idref=\"samplingProportion.s.t:dummy_aln\"/>
          </log>
      </logger>

      <logger id=\"screenlog\" fileName=\"\" logEvery=\"20000\">
          <log idref=\"posterior\"/>
          <log id=\"ESS.0\" spec=\"util.ESS\" arg=\"@posterior\"/>
          <log idref=\"likelihood\"/>
          <log idref=\"prior\"/>
      </logger>
  </run>
  </beast>
  "
  for(i in 1:length(trees_lines)){
    tr_temp <- read.tree(text = trees_lines[i])
    taxon_sets <- paste0("<taxon id=\"", tr_temp$tip.label, "\" spec=\"Taxon\"/>", collapse = "\n")
    xml_temp <- gsub('POSTERIOR_OUTPUT_FILE', rep_name, gsub("INPUT_TAXON_SETS", taxon_sets, gsub("INPUT_TREE_STRING", trees_lines[i], bd_template)))
    cat(xml_temp, file = paste0(rep_name, '.xml'), sep = '\n')
  }
}

run_beast_analyses <- function(beast_command, xml_file_name){
  system(paste(beast_command, xml_file_name))
  log_posterior <- read.table(gsub('[.]xml', '.log', xml_file_name), head = T)
  return(log_posterior)
}


beast_command <- '~/Desktop/phylo_programs/BEAST243/bin/beast'
