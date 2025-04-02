
compute_infected_matrix <- function(ages, years, foi, age_risk = NULL, rho = NULL) {
  
  I_matrix <- matrix(NA, nrow = length(years), ncol = length(ages))
  colnames(I_matrix) <- paste0(ages, "yo")
  rownames(I_matrix) <- years 
  
  for (year_idx in 1:length(years)) {
    year <- years[year_idx]  
    
    for (age_idx in 1:length(ages)) {
      age <- ages[age_idx]  
      birth_year <- year - age + 1  
      
      # Model 1/9: independent or two-outbreak model 
      if (is.null(rho) & is.null(age_risk)) {
        
        lambda_values <- foi$lambda[foi$Year %in% (year - (age - seq(1, age)))]
        P_infected_yearly <- numeric(age)
        P_not_infected_before <- 1 
        
        for (i in 1:age) {
          lambda_i <- lambda_values[i]
          P_infected_yearly[i] <- (1 - exp(-lambda_i)) * P_not_infected_before
          P_not_infected_before <- P_not_infected_before * exp(-lambda_i)  
          I_matrix[year_idx, age_idx] <- P_infected_yearly[i]
          
        }
      } 
      
      # Model 2/10: independent or two-outbreak model accounting for age-dependent FOI
      if (is.null(rho) & !is.null(age_risk)) {
        
        lambda_values <- foi$lambda[foi$Year %in% (year - (age - seq(1, age)))]
        P_infected_yearly <- numeric(age)
        P_not_infected_before <- 1 
        
        for (i in 1:age) {
          
          lambda_age <- lambda_values[i] * exp(age_risk * (i - 1))  
          P_infected_yearly[i] <- (1 - exp(-lambda_age)) * P_not_infected_before
          P_not_infected_before <- P_not_infected_before * exp(-lambda_age)  
          I_matrix[year_idx, age_idx] <- P_infected_yearly[i]
          
        }
      } 
      
      # Model 3/11: independent or two-outbreak model accounting for seroreversion
      if((!is.null(rho)) && is.null(age_risk))  {
        
        lambda_values <- foi$lambda[foi$Year %in% (year - (age - seq(1, age)))]
        P_infected_yearly <- numeric(age)
        X_a <- 0 # fully susceptible
        
        for (i in 1:age) {
          
          lambda_i <- lambda_values[i] 
          P_infected_yearly[i] <- (1 - X_a) * (1 - exp(-lambda_i))
          X_a <- X_a * exp(-(lambda_i + rho)) + (lambda_i / (lambda_i + rho)) * (1 - exp(-(lambda_i + rho)))
          I_matrix[year_idx, age_idx] <- P_infected_yearly[i]
          
        }
      } 
      
      # Model 4/12: independent or two-outbreak model accounting for age-dependent FOI and seroreversion
      if(!is.null(rho) && !is.null(age_risk))  {
        
        lambda_values <- foi$lambda[foi$Year %in% (year - (age - seq(1, age)))]
        P_infected_yearly <- numeric(age)
        X_a <- 0  # fully susceptible
        
        for (i in 1:age) {
          
          lambda_age <- lambda_values[i] * exp(age_risk * (i - 1))  
          P_infected_yearly[i] <- (1 - X_a) * (1 - exp(-lambda_age))
          X_a <- X_a * exp(-(lambda_age + rho)) + (lambda_age / (lambda_age + rho)) * (1 - exp(-(lambda_age + rho)))
          I_matrix[year_idx, age_idx] <- P_infected_yearly[i]
          
        }
      } 
      
    }
  }
  return(I_matrix)
}      




      