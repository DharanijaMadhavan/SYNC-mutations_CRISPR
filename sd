#Author: Benedict Rausher
require('mixtools')

pb <- progress_estimated(length(unique(depmap_ceres$symbol)))

#Fitting distributions and generating mean and standard deviation for each gene
mm_set <- depmap_ceres %>% group_by(symbol) %>%
  group_map(~ {
    pb$tick()$print()
    ## fit gaussian mixture
    mm <- tryCatch(normalmixEM(.x$cscore),
                   error=function(cond) return(NULL))
    if(!is.null(mm)){
      ## select largest component
      compl <- as_tibble(mm$posterior) %>%
        mutate(comp = ifelse(comp.1 > comp.2, 1, 2)) %>%
        count(comp) %>% arrange(desc(n)) %>% dplyr::slice(1) %>%
        pull(comp)
      ## extract paramters
      tibble(
        mu_null = mm$mu[compl],
        sigma_null = mm$sigma[compl]
      )
    } else {
      tibble(
        mu_null = NA,
        sigma_null = NA
      )
    }
  }) %>% ungroup()
 
 #Test to compare CERES of gene in cell line of interest vs rest of cell lines
 2*pnorm(x, mean= mm_set2$mu_null, sd= mm_set2$sigma_null ) #x= CERES of gene in cell line analysing
