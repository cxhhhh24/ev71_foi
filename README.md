##### 'compute_infected_matrix.R' and 'compute_susceptible_matrix.R' are two functions to calcualte susceptible proportion and annual incidence across different models

##### 'fig4_susceptible_inected_matrix_independent_model.R' is the code ,for independent model, to used aboved function to run the calculation and also to calculate the CIs by sampling from posterior distributions, as well as ploting  
- model 1: basic independent model
- model 2: independent model plus age-dependent FOI
- model 3: independent model plus seroreversion
- model 4: independent model plus age-dependent FOI and seroreversion
- line 28-132: a loop to calculate susceptible matrix, infected matrix, and mean age of infection from model 1 to model 4
- line 136-176: data formating, for the next-setp plotting
- line 231-430: fig4 plotting
- ps: 'compute_infected_matrix.R' and 'compute_susceptible_matrix.R' were used in line 82-83
##### 'sup_fig_susceptible_inected_matrix_outbreak_model.R' is the code ,for two-outbreak model, to used aboved function to run the calculation and also to calculate the CIs by sampling from posterior distributions, as well as ploting.  
- model 9: basic two-outbreak model
- model 10: two-outbreak model plus age-dependent FOI
- model 11: two-outbreak model plus seroreversion
- model 12: two-outbreak model plus age-dependent FOI and seroreversion
- line 28-132: a loop to calculate susceptible matrix, infected matrix, and mean age of infection from model 9 to model 12
- line 136-176: data formating, for the next-setp plotting
- line 231-430: Sup_fig plotting
- ps: 'compute_infected_matrix.R' and 'compute_susceptible_matrix.R' were used in line 82-83
