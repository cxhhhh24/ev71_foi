
compute_infected_matrix <- function(ages, years, foi, age_risk = NULL, rho = NULL) {

  I_matrix <- matrix(NA, nrow = length(years), ncol = length(ages))
  colnames(I_matrix) <- paste0(ages, "yo")
  rownames(I_matrix) <- years 

  for (year_idx in 1:length(years)) {
    year <- years[year_idx]  
    
    for (age_idx in 1:length(ages)) {
      age <- ages[age_idx]  
      # birth_year <- year - age + 1  
      
      # Model 1/9: independent or two-outbreak model 
      if (is.null(rho) & is.null(age_risk)) {
        
        lambda_values <- foi$lambda[foi$Year %in% (year - (age - seq(1, age)))]
        ar <- numeric(age)
        s_neg <- 1 
 
        for (i in 1:age) {
          
          lambda_i <- lambda_values[i]
          ar[i] <- s_neg * (1 - exp(-lambda_i)) 
          s_neg <- s_neg * exp(-lambda_i) 
          I_matrix[year_idx, age_idx] <- ar[i]
          
        }
      } 
      
      # Model 2/10: independent or two-outbreak model accounting for age-dependent FOI
      if (is.null(rho) & !is.null(age_risk)) {
        
        lambda_values <- foi$lambda[foi$Year %in% (year - (age - seq(1, age)))]
        ar <- numeric(age)
        s_neg <- 1 
        
        for (i in 1:age) {
          
          lambda_age_i <- lambda_values[i] * exp(age_risk * (i - 1))  
          ar[i] <- s_neg * (1 - exp(-lambda_age_i)) 
          s_neg <- s_neg * exp(-lambda_age_i) 
          I_matrix[year_idx, age_idx] <- ar[i]
          
        }
      } 
      
      # Model 3/11: independent or two-outbreak model accounting for seroreversion
      if((!is.null(rho)) && is.null(age_risk))  {
        
        lambda_values <- foi$lambda[foi$Year %in% (year - (age - seq(1, age)))]
        ar <- numeric(age)
        s_pos <- 0
        
        for (i in 1:age) {
          lambda_i <- lambda_values[i]
          s_neg <- 1 - s_pos
          ar[i] <- s_neg * (1 - exp(-lambda_i)) 
          s_pos <- (lambda_i / (lambda_i + rho)) +  (s_pos - (lambda_i / (lambda_i + rho))) * exp(-(lambda_i + rho))
          I_matrix[year_idx, age_idx] <- ar[i]
          
        }
      } 
      
      # Model 4/12: independent or two-outbreak model accounting for age-dependent FOI and seroreversion
      if(!is.null(rho) && !is.null(age_risk))  {
        
        lambda_values <- foi$lambda[foi$Year %in% (year - (age - seq(1, age)))]
        ar <- numeric(age)
        s_pos <- 0
        
        for (i in 1:age) {
          
          lambda_age_i <- lambda_values[i] * exp(age_risk * (i - 1))  
          s_neg <- 1 - s_pos
          ar[i] <- s_neg * (1 - exp(-lambda_age_i)) 
          s_pos <- (lambda_age_i / (lambda_age_i + rho)) +  (s_pos - (lambda_age_i / (lambda_age_i + rho))) * exp(-(lambda_age_i + rho))
          I_matrix[year_idx, age_idx] <- ar[i]
          
        }
      } 
      
    }
  }
  return(I_matrix)
}      




      