#' correlation with bootstrap errors
corBoot <- function(data, formula, nboot = 1000, ncpus = 1, ...){
  
  # Function to return correlation estimate
  getCor <- function(data, i, formula, ...){
    # New data
    new_data <- data[i, ]
    
    # Return the correlation estimate 
    cor.test(formula = as.formula(formula), 
             data = new_data, ...)$estimate
  }
  
  # Bootstram them
  cor_results <- boot::boot(data = data, 
                            statistic = getCor, 
                            R = nboot,
                            formula = formula, 
                            ncpus = ncpus, 
                            ...)
  
  # Output tidy results
  broom::tidy(cor_results, conf.int = TRUE, conf.method = "bca") %>% 
    select(statistic, conf.low, conf.high) %>% 
    rename(conflo = conf.low, confhi = conf.high)
}
