
get_mean_logfiles <- function(file_names){

    means_list <- list()
    for(i in 1:length(file_names)){
        dat_temp <- read.table(file_names[i], head = T)
        dat_temp <- dat_temp[-(1:ceiling(nrow(dat_temp)*0.1)), ]
        means_list[[i]] <- colMeans(dat_temp)
    }
    rbind.list <- function(vlist){
        if(length(vlist) == 1){
            return(vlist[[1]])
        }else if(length(vlist) == 2){
            return(rbind(vlist[[1]], vlist[[2]]))
        }else{
            return(rbind(vlist[[1]], rbind.list(vlist[-1])))
        }
    }

    return(rbind.list(means_list))
}


# To test:
#ce_log_files <- paste0('../test_data/', dir('../test_data/', pattern = '_ce_pps_.+log'))
#cc_log_files <- paste0('../test_data/', dir('../test_data/', pattern = '_cc_.+log'))

#ce_means <- get_mean_logfiles(ce_log_files)

#hist(ce_means[, 'TreeHeight'])
