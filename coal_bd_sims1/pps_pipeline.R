## Make

# - read tree and select an analysis template
#

# Input settings are:

library(ape)

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
      xml_temp <- gsub('POSTERIOR_OUTPUT_FILE', paste0(rep_name, '_', i), gsub("INPUT_TAXON_SETS", taxon_sets, gsub("INPUT_TREE_STRING", trees_lines[i], cc_template)))
      cat(xml_temp, file = paste0(rep_name, '_', i, '.xml'), sep = '\n')
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


  <run id=\"mcmc\" spec=\"MCMC\" chainLength=\"10000000\" sampleFromPrior=\"false\">
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
    xml_temp <- gsub('POSTERIOR_OUTPUT_FILE', paste0(rep_name, '_', i), gsub("INPUT_TAXON_SETS", taxon_sets, gsub("INPUT_TREE_STRING", trees_lines[i], ce_template)))
    cat(xml_temp, file = paste0(rep_name, '_', i, '.xml'), sep = '\n')
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
    xml_temp <- gsub('POSTERIOR_OUTPUT_FILE', paste0(rep_name, '_', i), gsub("INPUT_TAXON_SETS", taxon_sets, gsub("INPUT_TREE_STRING", trees_lines[i], bd_template)))
    cat(xml_temp, file = paste0(rep_name, '_', i, '.xml'), sep = '\n')
  }
}


run_beast_analyses <- function(beast_command, xml_file_name){
  system(paste(beast_command, xml_file_name))
  log_posterior <- read.table(gsub('[.]xml', '.log', xml_file_name), head = T)
  return(log_posterior)
}




make_cc_simulation <- function(posterior_log_data, input_tree, output_name){
  cc_template <- '<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?><beast beautitemplate=\'Standard\' beautistatus=\'\' namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" required=\"\" version=\"2.4\">

      <data id=\"dummy_aln\" name=\"alignment\">
        TAXON_SEQS
      </data>

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

  <run id=\"mcmc\" spec=\"MCMC\" chainLength=\"4000000\" sampleFromPrior=\"true\">
      <state id=\"state\" storeEvery=\"5000\">
          <tree id=\"Tree.t:dummy_aln\" name=\"stateNode\">
            <trait id=\"dateTrait.t:dummy_aln\" spec=\"beast.evolution.tree.TraitSet\" traitname=\"date\">
              TAXON_DATES
                <taxa id=\"TaxonSet.dummy_aln\" spec=\"TaxonSet\">
                  <alignment idref=\"dummy_aln\"/>
                </taxa>
              </trait>
              <taxonset idref=\"TaxonSet.dummy_aln\"/>
          </tree>
          <parameter id=\"clockRate.c:dummy_aln\" name=\"stateNode\">1.0</parameter>
          <parameter id=\"popSize.t:dummy_aln\" name=\"stateNode\">CONSTANT_POPSIZE_MEAN</parameter>
      </state>

      <init id=\"RandomTree.t:dummy_aln\" spec=\"beast.evolution.tree.RandomTree\" estimate=\"false\" initial=\"@Tree.t:dummy_aln\" taxa=\"@dummy_aln\">
          <populationModel id=\"ConstantPopulation0.t:dummy_aln\" spec=\"ConstantPopulation\">
              <parameter id=\"randomPopSize.t:dummy_aln\" name=\"popSize\">1</parameter>
          </populationModel>
      </init>

      <distribution id=\"posterior\" spec=\"util.CompoundDistribution\">
          <distribution id=\"prior\" spec=\"util.CompoundDistribution\">
              <distribution id=\"CoalescentConstant.t:dummy_aln\" spec=\"Coalescent\">
                  <populationModel id=\"ConstantPopulation.t:dummy_aln\" spec=\"ConstantPopulation\" popSize=\"@popSize.t:dummy_aln\"/>
                  <treeIntervals id=\"TreeIntervals.t:dummy_aln\" spec=\"TreeIntervals\" tree=\"@Tree.t:dummy_aln\"/>
              </distribution>
              <prior id=\"ClockPrior.c:dummy_aln\" name=\"distribution\" x=\"@clockRate.c:dummy_aln\">
                  <Uniform id=\"Uniform.0\" name=\"distr\" upper=\"Infinity\"/>
              </prior>
              <prior id=\"PopSizePrior.t:dummy_aln\" name=\"distribution\" x=\"@popSize.t:dummy_aln\">
  	      <Normal id=\"Normal.01\" name=\"distr\">
  		<parameter id=\"RealParameter.02\" estimate=\"false\" name=\"mean\">CONSTANT_POPSIZE_MEAN</parameter>
                  <parameter id=\"RealParameter.03\" estimate=\"false\" name=\"sigma\">CONSTANT_POPSIZE_SD</parameter>
  	      </Normal>
              </prior>
          </distribution>
          <distribution id=\"likelihood\" spec=\"util.CompoundDistribution\" useThreads=\"true\">
              <distribution id=\"treeLikelihood.dummy_aln\" spec=\"ThreadedTreeLikelihood\" data=\"@dummy_aln\" tree=\"@Tree.t:dummy_aln\">
                  <siteModel id=\"SiteModel.s:dummy_aln\" spec=\"SiteModel\">
                      <parameter id=\"mutationRate.s:dummy_aln\" estimate=\"false\" name=\"mutationRate\">1.0</parameter>
                      <parameter id=\"gammaShape.s:dummy_aln\" estimate=\"false\" name=\"shape\">1.0</parameter>
                      <parameter id=\"proportionInvariant.s:dummy_aln\" estimate=\"false\" lower=\"0.0\" name=\"proportionInvariant\" upper=\"1.0\">0.0</parameter>
                      <substModel id=\"JC69.s:dummy_aln\" spec=\"JukesCantor\"/>
                  </siteModel>
                  <branchRateModel id=\"StrictClock.c:dummy_aln\" spec=\"beast.evolution.branchratemodel.StrictClockModel\" clock.rate=\"@clockRate.c:dummy_aln\"/>
              </distribution>
          </distribution>
      </distribution>

      <operator id=\"StrictClockRateScaler.c:dummy_aln\" spec=\"ScaleOperator\" parameter=\"@clockRate.c:dummy_aln\" scaleFactor=\"0.75\" weight=\"3.0\"/>

      <operator id=\"strictClockUpDownOperator.c:dummy_aln\" spec=\"UpDownOperator\" scaleFactor=\"0.75\" weight=\"3.0\">
          <up idref=\"clockRate.c:dummy_aln\"/>
          <down idref=\"Tree.t:dummy_aln\"/>
      </operator>

      <operator id=\"CoalescentConstantTreeScaler.t:dummy_aln\" spec=\"ScaleOperator\" scaleFactor=\"0.5\" tree=\"@Tree.t:dummy_aln\" weight=\"3.0\"/>
      <operator id=\"CoalescentConstantTreeRootScaler.t:dummy_aln\" spec=\"ScaleOperator\" rootOnly=\"true\" scaleFactor=\"0.5\" tree=\"@Tree.t:dummy_aln\" weight=\"3.0\"/>
      <operator id=\"CoalescentConstantUniformOperator.t:dummy_aln\" spec=\"Uniform\" tree=\"@Tree.t:dummy_aln\" weight=\"30.0\"/>
      <operator id=\"CoalescentConstantSubtreeSlide.t:dummy_aln\" spec=\"SubtreeSlide\" tree=\"@Tree.t:dummy_aln\" weight=\"15.0\"/>
      <operator id=\"CoalescentConstantNarrow.t:dummy_aln\" spec=\"Exchange\" tree=\"@Tree.t:dummy_aln\" weight=\"15.0\"/>
      <operator id=\"CoalescentConstantWide.t:dummy_aln\" spec=\"Exchange\" isNarrow=\"false\" tree=\"@Tree.t:dummy_aln\" weight=\"3.0\"/>
      <operator id=\"CoalescentConstantWilsonBalding.t:dummy_aln\" spec=\"WilsonBalding\" tree=\"@Tree.t:dummy_aln\" weight=\"3.0\"/>
      <operator id=\"PopSizeScaler.t:dummy_aln\" spec=\"ScaleOperator\" parameter=\"@popSize.t:dummy_aln\" scaleFactor=\"0.75\" weight=\"3.0\"/>
      <logger id=\"tracelog\" fileName=\"CC_SIM_TREE_FILE.log\" logEvery=\"5000\" model=\"@posterior\" sanitiseHeaders=\"true\" sort=\"smart\">
          <log idref=\"posterior\"/>
          <log idref=\"likelihood\"/>
          <log idref=\"prior\"/>
          <log idref=\"treeLikelihood.dummy_aln\"/>
          <log id=\"TreeHeight.t:dummy_aln\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@Tree.t:dummy_aln\"/>
          <log idref=\"clockRate.c:dummy_aln\"/>
          <log idref=\"popSize.t:dummy_aln\"/>
          <log idref=\"CoalescentConstant.t:dummy_aln\"/>
      </logger>

      <logger id=\"screenlog\" fileName=\"\" logEvery=\"5000\">
          <log idref=\"posterior\"/>
  	<log idref=\"TreeHeight.t:dummy_aln\"/>
          <log id=\"ESS.0\" spec=\"util.ESS\" arg=\"@posterior\"/>
          <log idref=\"likelihood\"/>
          <log idref=\"prior\"/>
      </logger>

      <logger id=\"treelog.t:dummy_aln\" fileName=\"CC_SIM_TREE_FILE.trees\" logEvery=\"5000\" mode=\"tree\">
          <log id=\"TreeWithMetaDataLogger.t:dummy_aln\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@Tree.t:dummy_aln\"/>
      </logger>

  </run>

  </beast>
  '
  taxon_seqs <- paste0("<sequence id=\"seq_", input_tree$tip.label, "\" taxon=\"",input_tree$tip.label, "\" totalcount=\"4\" value=\"gc\"/>", collapse = '\n')

  taxon_dates <- vector()
  for(ta in 1:length(input_tree$tip.label)){
       date <- gsub('.+_', '', input_tree$tip.label[ta])
       taxon_dates[ta] <- paste0(input_tree$tip.label[ta], '=', date)
     }
  taxon_dates <- paste0(taxon_dates, collapse = ',\n')
  popsize_normal <- round(c(mean(posterior_log_data$popSize), sd(posterior_log_data$popSize)), 2)
  xml_temp <- gsub('TAXON_DATES', taxon_dates, gsub('TAXON_SEQS', taxon_seqs, cc_template))
  xml_temp <- gsub('CC_SIM_TREE_FILE', paste0(output_name, '_pps'), xml_temp)
  xml_temp <- gsub('CONSTANT_POPSIZE_MEAN', popsize_normal[1], xml_temp)
  xml_temp <- gsub('CONSTANT_POPSIZE_SD', popsize_normal[2], xml_temp)
  cat(xml_temp, file = paste0(output_name, '_pps.xml'), sep = '\n')
}

make_ce_simulation <- function(posterior_log_data, input_tree, output_name){
  ce_template <- '<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?><beast beautitemplate=\'Standard\' beautistatus=\'\' namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" required=\"\" version=\"2.4\">

      <data id=\"dummy_aln\" name=\"alignment\">
        TAXON_SEQS
      </data>

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

  <run id=\"mcmc\" spec=\"MCMC\" chainLength=\"10000000\" sampleFromPrior=\"true\">
      <state id=\"state\" storeEvery=\"5000\">
        <tree id=\"Tree.t:dummy_aln\" name=\"stateNode\">
          <trait id=\"dateTrait.t:dummy_aln\" spec=\"beast.evolution.tree.TraitSet\" traitname=\"date\">
            TAXON_DATES
            <taxa id=\"TaxonSet.dummy_aln\" spec=\"TaxonSet\">
              <alignment idref=\"dummy_aln\"/>
            </taxa>
          </trait>
          <taxonset idref=\"TaxonSet.dummy_aln\"/>
        </tree>
        <parameter id=\"ePopSize.t:dummy_aln\" name=\"stateNode\">E_POP_SIZE_MEAN</parameter>
        <parameter id=\"growthRate.t:dummy_aln\" name=\"stateNode\">GROWTH_RATE_MEAN</parameter>
        <parameter id=\"clockRate.c:dummy_aln\" name=\"stateNode\">1.0</parameter>
      </state>

      <init id=\"RandomTree.t:dummy_aln\" spec=\"beast.evolution.tree.RandomTree\" estimate=\"false\" initial=\"@Tree.t:dummy_aln\" taxa=\"@dummy_aln\">
          <populationModel id=\"ConstantPopulation0.t:dummy_aln\" spec=\"ConstantPopulation\">
              <parameter id=\"randomPopSize.t:dummy_aln\" name=\"popSize\">1</parameter>
          </populationModel>
      </init>

      <distribution id=\"posterior\" spec=\"util.CompoundDistribution\">
        <distribution id=\"prior\" spec=\"util.CompoundDistribution\">
          <distribution id=\"CoalescentExponential.t:dummy_aln\" spec=\"Coalescent\">
            <populationModel id=\"ExponentialGrowth.t:dummy_aln\" spec=\"ExponentialGrowth\" growthRate=\"@growthRate.t:dummy_aln\" popSize=\"@ePopSize.t:dummy_aln\"/>
            <treeIntervals id=\"TreeIntervals.t:dummy_aln\" spec=\"TreeIntervals\" tree=\"@Tree.t:dummy_aln\"/>
          </distribution>
          <prior id=\"ePopSizePrior.t:dummy_aln\" name=\"distribution\" x=\"@ePopSize.t:dummy_aln\">
            <Normal id=\"Normal.01\" name=\"distr\">
              <parameter id=\"RealParameter.0\" estimate=\"false\" name=\"mean\">E_POP_SIZE_MEAN</parameter>
              <parameter id=\"RealParameter.01\" estimate=\"false\" name=\"sigma\">E_POP_SIZE_SD</parameter>
            </Normal>
          </prior>
          <prior id=\"GrowthRatePrior.t:dummy_aln\" name=\"distribution\" x=\"@growthRate.t:dummy_aln\">
            <Normal id=\"Normal.02\" name=\"distr\">
              <parameter id=\"RealParameter.02\" estimate=\"false\" name=\"mean\">GROWTH_RATE_MEAN</parameter>
              <parameter id=\"RealParameter.03\" estimate=\"false\" name=\"sigma\">GROWTH_RATE_SD</parameter>
            </Normal>
          </prior>
          <prior id=\"ClockPrior.c:dummy_aln\" name=\"distribution\" x=\"@clockRate.c:dummy_aln\">
            <Uniform id=\"Uniform.0\" name=\"distr\" upper=\"Infinity\"/>
          </prior>
        </distribution>

        <distribution id=\"likelihood\" spec=\"util.CompoundDistribution\" useThreads=\"true\">
          <distribution id=\"treeLikelihood.dummy_aln\" spec=\"ThreadedTreeLikelihood\" data=\"@dummy_aln\" tree=\"@Tree.t:dummy_aln\">
            <siteModel id=\"SiteModel.s:dummy_aln\" spec=\"SiteModel\">
              <parameter id=\"mutationRate.s:dummy_aln\" estimate=\"false\" name=\"mutationRate\">1.0</parameter>
              <parameter id=\"gammaShape.s:dummy_aln\" estimate=\"false\" name=\"shape\">1.0</parameter>
              <parameter id=\"proportionInvariant.s:dummy_aln\" estimate=\"false\" lower=\"0.0\" name=\"proportionInvariant\" upper=\"1.0\">0.0</parameter>
              <substModel id=\"JC69.s:dummy_aln\" spec=\"JukesCantor\"/>
            </siteModel>
          <branchRateModel id=\"StrictClock.c:dummy_aln\" spec=\"beast.evolution.branchratemodel.StrictClockModel\" clock.rate=\"@clockRate.c:dummy_aln\"/>
  	</distribution>
        </distribution>
      </distribution>

      <operator id=\"CoalescentExponentialTreeScaler.t:dummy_aln\" spec=\"ScaleOperator\" scaleFactor=\"0.5\" tree=\"@Tree.t:dummy_aln\" weight=\"3.0\"/>
      <operator id=\"CoalescentExponentialTreeRootScaler.t:dummy_aln\" spec=\"ScaleOperator\" rootOnly=\"true\" scaleFactor=\"0.5\" tree=\"@Tree.t:dummy_aln\" weight=\"3.0\"/>
      <operator id=\"CoalescentExponentialUniformOperator.t:dummy_aln\" spec=\"Uniform\" tree=\"@Tree.t:dummy_aln\" weight=\"30.0\"/>
      <operator id=\"CoalescentExponentialSubtreeSlide.t:dummy_aln\" spec=\"SubtreeSlide\" tree=\"@Tree.t:dummy_aln\" weight=\"15.0\"/>
      <operator id=\"CoalescentExponentialNarrow.t:dummy_aln\" spec=\"Exchange\" tree=\"@Tree.t:dummy_aln\" weight=\"15.0\"/>
      <operator id=\"CoalescentExponentialWide.t:dummy_aln\" spec=\"Exchange\" isNarrow=\"false\" tree=\"@Tree.t:dummy_aln\" weight=\"3.0\"/>
      <operator id=\"CoalescentExponentialWilsonBalding.t:dummy_aln\" spec=\"WilsonBalding\" tree=\"@Tree.t:dummy_aln\" weight=\"3.0\"/>
      <operator id=\"ePopSizeScaler.t:dummy_aln\" spec=\"ScaleOperator\" parameter=\"@ePopSize.t:dummy_aln\" scaleFactor=\"0.75\" weight=\"3.0\"/>
      <operator id=\"GrowthRateRandomWalk.t:dummy_aln\" spec=\"RealRandomWalkOperator\" parameter=\"@growthRate.t:dummy_aln\" weight=\"3.0\" windowSize=\"1.0\"/>

      <logger id=\"tracelog\" fileName=\"CE_SIM_TREE_FILE.log\" logEvery=\"1000\" model=\"@posterior\" sanitiseHeaders=\"true\" sort=\"smart\">
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
  <!-- sample trees a bit less often -->
      <logger id=\"treelog.t:dummy_aln\" fileName=\"CE_SIM_TREE_FILE.trees\" logEvery=\"10000\" mode=\"tree\">
          <log id=\"TreeWithMetaDataLogger.t:dummy_aln\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@Tree.t:dummy_aln\"/>
      </logger>
  </run>
  </beast>'
  taxon_seqs <- paste0("<sequence id=\"seq_", input_tree$tip.label, "\" taxon=\"",input_tree$tip.label, "\" totalcount=\"4\" value=\"gc\"/>", collapse = '\n')

  taxon_dates <- vector()
  for(ta in 1:length(input_tree$tip.label)){
       date <- gsub('.+_', '', input_tree$tip.label[ta])
       taxon_dates[ta] <- paste0(input_tree$tip.label[ta], '=', date)
     }
  taxon_dates <- paste0(taxon_dates, collapse = ',\n')
  epopsize <- round(c(mean(posterior_log_data$ePopSize), sd(posterior_log_data$ePopSize)), 2)
  growthrate <- round(c(mean(posterior_log_data$growthRate.), sd(posterior_log_data$growthRate.)), 2)
  xml_temp <- gsub('TAXON_DATES', taxon_dates, gsub('TAXON_SEQS', taxon_seqs, ce_template))
  xml_temp <- gsub('CE_SIM_TREE_FILE', paste0(output_name, '_pps'), xml_temp)
  xml_temp <- gsub('E_POP_SIZE_MEAN', epopsize[1], xml_temp)
  xml_temp <- gsub('E_POP_SIZE_SD', epopsize[2], xml_temp)
  xml_temp <- gsub('GROWTH_RATE_MEAN', growthrate[1], xml_temp)
  xml_temp <- gsub('GROWTH_RATE_SD', growthrate[2], xml_temp)

  cat(xml_temp, file = paste0(output_name, '_pps.xml'), sep = '\n')

}

make_bd_simulation <- function(posterior_log_data, input_tree, output_name){
  bd_template <- '<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?><beast beautitemplate=\'Standard\' beautistatus=\'\' namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" required=\"BDSKY v1.3.2\" version=\"2.4\">

      <data id=\"dummy_aln\" name=\"alignment\">
        TAXON_SEQS
      </data>

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


  <run id=\"mcmc\" spec=\"MCMC\" chainLength=\"5000000\" sampleFromPrior=\"true\">
      <state id=\"state\" storeEvery=\"1000\">
          <tree id=\"Tree.t:dummy_aln\" name=\"stateNode\">
              <trait id=\"dateTrait.t:dummy_aln\" spec=\"beast.evolution.tree.TraitSet\" traitname=\"date\">
  	      TAXON_DATES
  	      <taxa id=\"TaxonSet.dummy_aln\" spec=\"TaxonSet\">
                  <alignment idref=\"dummy_aln\"/>
                </taxa>
              </trait>
              <taxonset idref=\"TaxonSet.dummy_aln\"/>
          </tree>
          <parameter id=\"clockRate.c:dummy_aln\" name=\"stateNode\">1.0</parameter>
          <parameter id=\"origin.t:dummy_aln\" lower=\"0.0\" name=\"stateNode\" upper=\"Infinity\">ORIGIN_MEAN</parameter>
          <parameter id=\"samplingProportion.t:dummy_aln\" lower=\"0.0\" name=\"stateNode\" upper=\"1.0\">SAMPLING_PROP_MEAN</parameter>
          <parameter id=\"becomeUninfectiousRate.t:dummy_aln\" lower=\"0.0\" name=\"stateNode\" upper=\"Infinity\">BECOME_UNINFECTIOUS_MEAN</parameter>
          <parameter id=\"R0.t:dummy_aln\" dimension=\"10\" lower=\"0.0\" name=\"stateNode\" upper=\"Infinity\">R0_MEAN</parameter>
      </state>

      <init id=\"RandomTree.t:dummy_aln\" spec=\"beast.evolution.tree.RandomTree\" estimate=\"false\" initial=\"@Tree.t:dummy_aln\" taxa=\"@dummy_aln\">
          <populationModel id=\"ConstantPopulation0.t:dummy_aln\" spec=\"ConstantPopulation\">
              <parameter id=\"randomPopSize.t:dummy_aln\" name=\"popSize\">1.0</parameter>
          </populationModel>
      </init>

      <distribution id=\"posterior\" spec=\"util.CompoundDistribution\">
          <distribution id=\"prior\" spec=\"util.CompoundDistribution\">
              <distribution id=\"BirthDeathSkySerial.t:dummy_aln\" spec=\"beast.evolution.speciation.BirthDeathSkylineModel\" R0=\"@R0.t:dummy_aln\" becomeUninfectiousRate=\"@becomeUninfectiousRate.t:dummy_aln\" origin=\"@origin.t:dummy_aln\" samplingProportion=\"@samplingProportion.t:dummy_aln\" tree=\"@Tree.t:dummy_aln\"/>
              <prior id=\"RPrior.s.t:dummy_aln\" name=\"distribution\" x=\"@R0.t:dummy_aln\">
                  <Normal id=\"Normal.0\" name=\"distr\">
                      <parameter id=\"RealParameter.0\" estimate=\"false\" name=\"mean\">R0_MEAN</parameter>
                      <parameter id=\"RealParameter.01\" estimate=\"false\" name=\"sigma\">R0_SD</parameter>
                  </Normal>
              </prior>
              <prior id=\"becomeUninfectiousRatePrior.t:dummy_aln\" name=\"distribution\" x=\"@becomeUninfectiousRate.t:dummy_aln\">
                  <Normal id=\"Normal.01\" name=\"distr\">
                      <parameter id=\"RealParameter.02\" estimate=\"false\" name=\"mean\">BECOME_UNINFECTIOUS_MEAN</parameter>
                      <parameter id=\"RealParameter.03\" estimate=\"false\" name=\"sigma\">BECOME_UNINFECTIOUS_SD</parameter>
                  </Normal>
              </prior>
              <prior id=\"ClockPrior.c:dummy_aln\" name=\"distribution\" x=\"@clockRate.c:dummy_aln\">
                  <Uniform id=\"Uniform.0\" name=\"distr\" upper=\"Infinity\"/>
              </prior>
              <prior id=\"originPrior.t:dummy_aln\" name=\"distribution\" x=\"@origin.t:dummy_aln\">

                <Normal id=\"Normal.03\" name=\"distr\">
                  <parameter id=\"RealParameter.06\" estimate=\"false\" name=\"mean\">ORIGIN_MEAN</parameter>
                  <parameter id=\"RealParameter.07\" estimate=\"false\" name=\"sigma\">ORIGIN_SD</parameter>
                </Normal>

              </prior>
              <prior id=\"samplingProportionPrior.t:dummy_aln\" name=\"distribution\" x=\"@samplingProportion.t:dummy_aln\">
                  <Normal id=\"Normal.02\" name=\"distr\">
                      <parameter id=\"RealParameter.04\" estimate=\"false\" name=\"mean\">SAMPLING_PROP_MEAN</parameter>
  		    <parameter id=\"RealParameter.05\" estimate=\"false\" name=\"sigma\">0.001</parameter>
  		    <!--
  			   <parameter id=\"RealParameter.05\" estimate=\"false\" name=\"sigma\">SAMPLING_PROP_SD</parameter>
  			   -->
                  </Normal>
              </prior>
          </distribution>
          <distribution id=\"likelihood\" spec=\"util.CompoundDistribution\" useThreads=\"true\">
              <distribution id=\"treeLikelihood.dummy_aln\" spec=\"ThreadedTreeLikelihood\" data=\"@dummy_aln\" tree=\"@Tree.t:dummy_aln\">
                  <siteModel id=\"SiteModel.s:dummy_aln\" spec=\"SiteModel\">
                      <parameter id=\"mutationRate.s:dummy_aln\" estimate=\"false\" name=\"mutationRate\">1.0</parameter>
                      <parameter id=\"gammaShape.s:dummy_aln\" estimate=\"false\" name=\"shape\">1.0</parameter>
                      <parameter id=\"proportionInvariant.s:dummy_aln\" estimate=\"false\" lower=\"0.0\" name=\"proportionInvariant\" upper=\"1.0\">0.0</parameter>
                      <substModel id=\"JC69.s:dummy_aln\" spec=\"JukesCantor\"/>
                  </siteModel>
                  <branchRateModel id=\"StrictClock.c:dummy_aln\" spec=\"beast.evolution.branchratemodel.StrictClockModel\" clock.rate=\"@clockRate.c:dummy_aln\"/>
              </distribution>
          </distribution>
      </distribution>

      <operator id=\"StrictClockRateScaler.c:dummy_aln\" spec=\"ScaleOperator\" parameter=\"@clockRate.c:dummy_aln\" scaleFactor=\"0.75\" weight=\"3.0\"/>

      <operator id=\"strictClockUpDownOperator.c:dummy_aln\" spec=\"UpDownOperator\" scaleFactor=\"0.75\" weight=\"3.0\">
          <up idref=\"clockRate.c:dummy_aln\"/>
          <down idref=\"Tree.t:dummy_aln\"/>
      </operator>

      <operator id=\"BDSKY_serialtreeScaler.t:dummy_aln\" spec=\"ScaleOperator\" scaleFactor=\"0.5\" tree=\"@Tree.t:dummy_aln\" weight=\"3.0\"/>

      <operator id=\"BDSKY_serialtreeRootScaler.t:dummy_aln\" spec=\"ScaleOperator\" rootOnly=\"true\" scaleFactor=\"0.5\" tree=\"@Tree.t:dummy_aln\" weight=\"3.0\"/>

      <operator id=\"BDSKY_serialUniformOperator.t:dummy_aln\" spec=\"Uniform\" tree=\"@Tree.t:dummy_aln\" weight=\"30.0\"/>

      <operator id=\"BDSKY_serialSubtreeSlide.t:dummy_aln\" spec=\"SubtreeSlide\" tree=\"@Tree.t:dummy_aln\" weight=\"15.0\"/>

      <operator id=\"BDSKY_serialnarrow.t:dummy_aln\" spec=\"Exchange\" tree=\"@Tree.t:dummy_aln\" weight=\"15.0\"/>

      <operator id=\"BDSKY_serialwide.t:dummy_aln\" spec=\"Exchange\" isNarrow=\"false\" tree=\"@Tree.t:dummy_aln\" weight=\"3.0\"/>

      <operator id=\"BDSKY_serialWilsonBalding.t:dummy_aln\" spec=\"WilsonBalding\" tree=\"@Tree.t:dummy_aln\" weight=\"3.0\"/>

      <operator id=\"becomeUninfectiousRateScaler.t:dummy_aln\" spec=\"ScaleOperator\" parameter=\"@becomeUninfectiousRate.t:dummy_aln\" scaleFactor=\"0.75\" weight=\"2.0\"/>

      <operator id=\"samplingScaler.s.t:dummy_aln\" spec=\"ScaleOperator\" parameter=\"@samplingProportion.t:dummy_aln\" scaleFactor=\"0.75\" weight=\"2.0\"/>


      <operator id=\"RScaler.s.t:dummy_aln\" spec=\"ScaleOperator\" parameter=\"@R0.t:dummy_aln\" scaleFactor=\"0.75\" weight=\"10.0\"/>

      <operator id=\"updownBD.s.t:dummy_aln\" spec=\"UpDownOperator\" scaleFactor=\"0.75\" weight=\"2.0\">
          <up idref=\"R0.t:dummy_aln\"/>
          <down idref=\"becomeUninfectiousRate.t:dummy_aln\"/>
      </operator>

      <operator id=\"origScaler.t:dummy_aln\" spec=\"ScaleOperator\" parameter=\"@origin.t:dummy_aln\" scaleFactor=\"0.75\" weight=\"1.0\"/>

      <logger id=\"tracelog\" fileName=\"BD_SIM_TREE_FILE.log\" logEvery=\"1000\" model=\"@posterior\" sanitiseHeaders=\"true\" sort=\"smart\">
          <log idref=\"posterior\"/>
          <log idref=\"likelihood\"/>
          <log idref=\"prior\"/>
          <log idref=\"treeLikelihood.dummy_aln\"/>
          <log id=\"TreeHeight.t:dummy_aln\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@Tree.t:dummy_aln\"/>
          <log idref=\"clockRate.c:dummy_aln\"/>
          <log idref=\"BirthDeathSkySerial.t:dummy_aln\"/>
          <log idref=\"origin.t:dummy_aln\"/>
          <log idref=\"samplingProportion.t:dummy_aln\"/>
          <log idref=\"becomeUninfectiousRate.t:dummy_aln\"/>
          <log idref=\"R0.t:dummy_aln\"/>
          <log id=\"birth.t:dummy_aln\" spec=\"beast.math.statistic.RPNcalculator\" expression=\"R0.t:dummy_aln becomeUninfectiousRate.t:dummy_aln *\">
              <parameter idref=\"becomeUninfectiousRate.t:dummy_aln\"/>
              <parameter idref=\"R0.t:dummy_aln\"/>
          </log>
          <log id=\"death.t:dummy_aln\" spec=\"beast.math.statistic.RPNcalculator\" expression=\"becomeUninfectiousRate.t:dummy_aln 1 samplingProportion.t:dummy_aln - *\">
              <parameter idref=\"becomeUninfectiousRate.t:dummy_aln\"/>
              <parameter idref=\"samplingProportion.t:dummy_aln\"/>
          </log>
          <log id=\"sampling.t:dummy_aln\" spec=\"beast.math.statistic.RPNcalculator\" expression=\"becomeUninfectiousRate.t:dummy_aln samplingProportion.t:dummy_aln *\">
              <parameter idref=\"becomeUninfectiousRate.t:dummy_aln\"/>
              <parameter idref=\"samplingProportion.t:dummy_aln\"/>
          </log>
      </logger>

      <logger id=\"screenlog\" fileName=\"\" logEvery=\"1000\">
          <log idref=\"posterior\"/>
          <log id=\"ESS.0\" spec=\"util.ESS\" arg=\"@posterior\"/>
          <log idref=\"likelihood\"/>
          <log idref=\"prior\"/>
      </logger>

      <logger id=\"treelog.t:dummy_aln\" fileName=\"BD_SIM_TREE_FILE.trees\" logEvery=\"5000\" mode=\"tree\">
          <log id=\"TreeWithMetaDataLogger.t:dummy_aln\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@Tree.t:dummy_aln\"/>
      </logger>

  </run>

  </beast>'

    taxon_seqs <- paste0("<sequence id=\"seq_", input_tree$tip.label, "\" taxon=\"",input_tree$tip.label, "\" totalcount=\"4\" value=\"gc\"/>", collapse = '\n')
    taxon_dates <- vector()
    for(ta in 1:length(input_tree$tip.label)){
         date <- gsub('.+_', '', input_tree$tip.label[ta])
         taxon_dates[ta] <- paste0(input_tree$tip.label[ta], '=', date)
       }
    taxon_dates <- paste0(taxon_dates, collapse = ',\n')

    xml_temp <- gsub('TAXON_DATES', taxon_dates, gsub('TAXON_SEQS', taxon_seqs, bd_template))

    origin <- round(c(mean(posterior_log_data$origin), sd(posterior_log_data$origin)), 2)
    sampling <- round(c(mean(posterior_log_data$samplingProportion),
                        sd(posterior_log_data$samplingProportion)), 2)
    becomeUninfect <- round(c(mean(posterior_log_data$becomeUninfectiousRate),
                              sd(posterior_log_data$becomeUninfectiousRate)), 2)
    r0s <- posterior_log_data[, grep('R0', colnames(posterior_log_data))]
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

    xml_temp <- gsub('BD_SIM_TREE_FILE', paste0(output_name, '_pps'), xml_temp)

    cat(xml_temp, file = paste0(output_name, '_pps.xml'), sep = '\n')


}

# Run beast simulation; run beast and collect trees only

run_beast_simulation <- function(beast_command, xml_file_name){
    system(paste(beast_command, xml_file_name))
    trees_sampled <- read.nexus(gsub('[.]xml', '.trees', xml_file_name))
    return(trees_sampled[sample(1:length(trees_sampled), 100)])# Note that we only sample 100 trees
}

run_beast_pps <- function(beast_command, xml_files){
  library(foreach)
  library(doParallel)
  cl <- makeCluster(5)
  registerDoParallel(cl)

  run_rep <- function(x){
    log_temp <- run_beast_analyses(beast_command, x)
    return(colMeans(log_temp))
   }
  output <- foreach(x = xml_files, .export = 'run_beast_analyses', .combine = rbind) %dopar% run_rep(x)
  stopCluster(cl)
  return(output)
}

## read tree
## make template for one of the model
## Run beast and collect log file
## make simulations
## run simulations
## run beast pps

# make_cc_template
# run_beast_analyses
# make_cc_simulation
# run_beast_simulation
# run_beast_pps -> save logs in an pps file

trees_lines <- readLines('dated.tree')
rep_name <- 'rep_test'
beast_command <- '~/Desktop/phylo_programs/BEAST243/bin/beast'
