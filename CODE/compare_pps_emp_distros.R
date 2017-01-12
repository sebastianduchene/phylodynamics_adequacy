
compare_pps_emp_distros <- function(emp_logs, pps_logs, outfile = 'pps_emp_comp'){
    pdf(paste0(outfile, '.pdf'), useDingbats = F, width = 5, height = 4.7)
    par (mfrow = c(2, 3))
    if('bd' %in% names(pps_logs)){
  # CC
        plot(density(emp_logs$cc$popSize[-(1:100)]), main = '',
             xlab = 'popSize')
        lines(density(pps_logs$cc$popSize[-(1:100)]), col = 'red')
    }

    if('ce' %in% names(pps_logs)){
  #CE
        plot(density(emp_logs$ce$growthRate.[-(1:100)]), main = '',
             xlab = 'growthRate')
        lines(density(pps_logs$ce$growthRate.[-(1:100)]), col = 'red')

        plot(density(emp_logs$ce$ePopSize[-(1:100)]), main = '',
             xlab = 'ePopSize')
        lines(density(pps_logs$ce$ePopSize[-(1:100)]), col = 'red')
    }

    if('bd' %in% names(pps_logs)){
   # BD
        plot(density(emp_logs$bd$becomeUninfectiousRate[-(1:100)]), main = '',
             xlab = 'becomeUninfectiousRate')
        lines(density(pps_logs$bd$becomeUninfectiousRate[-(1:100)]), col = 'red')

        plot(density(emp_logs$bd$R0.1[-(1:100)]), main = '',
             xlab = 'R0.1')
        lines(density(pps_logs$bd$R0.1[-(1:100)]), col = 'red')
    }
    dev.off()
}




#emp_logs <- list(cc = read.table('../test_data/50taxa_CE_1e10_08_rep_81_cc_1.log', head = T), ce = read.table('../test_data/50taxa_CE_1e10_08_rep_81_ce_1.log', head = T), bd = read.table('../test_data/50taxa_CE_1e10_08_rep_81_bd_1.log', head = T))

#pps_logs <- list(cc = read.table('../test_data/50taxa_CE_1e10_08_rep_81_cc_pps.log', head = T), ce = read.table('../test_data/50taxa_CE_1e10_08_rep_81_ce_pps.log', head = T), bd = read.table('../test_data/50taxa_CE_1e10_08_rep_81_bd_pps.log', head = T))


#compare_pps_emp_distros(emp_logs, pps_logs)
