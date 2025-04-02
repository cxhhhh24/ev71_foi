
compute_susceptible_matrix <- function(ages, years, foi, age_risk = NULL, rho = NULL) {
  
  
  #- creat empty matrix
  S_matrix <- matrix(NA, nrow = length(years), ncol = length(ages))
  colnames(S_matrix) <- paste0(ages, "yo")
  rownames(S_matrix) <- years 
  
  # loop
  for (year_idx in 1:length(years)) {
    
    year <- years[year_idx]  
    
    for (age_idx in 1:length(ages)) {
      age <- ages[age_idx]  
      birth_year <- year - age + 1  
      
      # Model 1/9: independent or two-outbreak model 
      if (is.null(rho) & is.null(age_risk)) {
        
        if (birth_year < min(years)) {
          lambda_values <- NA 
        } else {
          lambda_values <- foi$lambda[foi$Year %in% (year - (age - seq(1,age)))]
        } 
        lambda_exp <- sum(lambda_values)  
        P_infected <- 1 - exp(-lambda_exp)  
        S_matrix[year_idx, age_idx] <- 1 - P_infected
      } 
      
      # Model 2/10: independent or two-outbreak model accounting for age-dependent FOI 
      if (is.null(rho) & !is.null(age_risk)) {
        
        if (birth_year < min(years)) {
          lambda_values <- NA 
        } else {
          lambda_values <- foi$lambda[foi$Year %in% (year - (age - seq(1,age)))]
        } 
        lambda_exp <- sum(lambda_values * exp(age_risk * (seq(age,1) - 1)))  
        P_infected <- 1 - exp(-lambda_exp)
        S_matrix[year_idx, age_idx] <- 1 - P_infected
        
      }
      
      # Model 3/11: independent or two-outbreak model accounting for seroreversion 
      if ((!is.null(rho)) && is.null(age_risk)) {
        
        if (birth_year < min(years)) {
          next
        }
        X <- 0  
        for (i in birth_year:year) {  
          foi_i <- foi[which(years == i), "lambda"]  
          X <- X * exp(-(foi_i + rho)) + (foi_i / (foi_i + rho)) * (1 - exp(-(foi_i + rho)))
        }
        S_matrix[year_idx, which(ages == age)] <- 1 - X      
        
      }
      
      # Model 4/12: independent or two-outbreak model accounting for age-dependent FOI and seroreversion  
      if (!is.null(rho) && !is.null(age_risk)) {
        
        if (birth_year < min(years)) {
          next
        }
        X <- 0  
        for (y in birth_year:year) {  
          foi_y <- foi[which(years == y), "lambda"]  
          age_risk_y <- exp(age_risk * (age - 1))  
          lambda_y_age <- foi_y * age_risk_y 
          X <- X * exp(-(lambda_y_age + rho)) + (lambda_y_age / (lambda_y_age + rho)) * (1 - exp(-(lambda_y_age + rho)))
        }
        
        S_matrix[year_idx, which(ages == age)] <- 1 - X
      }
      
     }
   }
  return(S_matrix)
}
