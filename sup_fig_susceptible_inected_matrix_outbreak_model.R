rm(list = ls())

library(Rsero)
library(rstan)
library(loo)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

model_fit <- readRDS("model_fit_12model_0326.rds")
model_plot <- readRDS("model_plot_12model_0326.rds")
source("compute_susceptible_matrix.R") # function to calculate susceptible matrix
source("compute_infected_matrix.R") # function to calculate infected matrix



####- model selection -####
models <- 9:12 # two-outbreak model

####- data input -####
years <- 1958:1973    
ages <- 1:9          
n_years <- length(years)
age_labels <- paste0(ages, "yo")  

####- loop for computation -####
data <- list()
mean_age <- data.frame()
for (model in models) {
  
  posterior_sample <- rstan::extract(model_fit[[model]][[1]], permuted = FALSE)
  
  n_samples <- dim(posterior_sample)[1]  # 2500 after removing warm-up per chain
  n_chains <- dim(posterior_sample)[2]   
  n_iter <- n_samples * n_chains  
  
  # bootstrap
  n_sim <- 1000
  
  # creat sample sets
  foi_samples <- matrix(NA, nrow = n_iter, ncol = n_years)
  for (i in 1:n_years) {

    foi_samples[, i] <- as.vector(posterior_sample[, , paste0("lambda[", 16 - (i - 1), "]")])
  }
  
  if (model %in% c(2, 4)) {
    age_risk_samples <- as.vector(posterior_sample[, , "age_risk"])
  }
  if (model %in% c(3, 4)) {
    rho_samples <- as.vector(posterior_sample[, , "rho"])
  }
  
  # sampling
  S_matrix <- array(NA, dim = c(n_years, length(ages), n_sim))
  I_matrix <- array(NA, dim = c(n_years, length(ages), n_sim))
  
  
  for (i in 1:n_sim) {
    
    foi_sampled_values <- foi_samples[sample(1:n_iter, 1, replace = TRUE), ]
    foi_sampled <- data.frame(
      Year = years,
      lambda = foi_sampled_values
    )
    
    if (model %in% c(2, 4)) {
      age_risk_sampled <- sample(age_risk_samples, 1, replace = TRUE)
    } else {
      age_risk_sampled <- NULL
    }
    
    if (model %in% c(3, 4)) {
      rho_sampled <- sample(rho_samples, 1, replace = TRUE)
    } else {
      rho_sampled <- NULL
    }
    
    # compute susceptible matrix
    S_matrix[, , i] <- compute_susceptible_matrix(ages, years, foi_sampled, age_risk_sampled, rho_sampled)
    I_matrix[, , i] <- compute_infected_matrix(ages, years, foi_sampled, age_risk_sampled, rho_sampled)
    
    
    # extract age
    age_data <- as.data.frame(I_matrix[, , i])
    age_data$year <- years
    age_data <- gather(age_data, key = "age", value = "value", -year)
    age_data$model <- model
    age_data$sim <- i
    age_data$age <- as.numeric(gsub("V", "", age_data$age)) 
    mean_age <- rbind(mean_age, age_data)
    
  }
  
  
  #aggregated S matrix
  median_S_matrix <- apply(S_matrix, c(1, 2), median)
  lower_S_matrix <- apply(S_matrix, c(1, 2), function(x) quantile(x, probs = 0.025, na.rm = TRUE))
  upper_S_matrix <- apply(S_matrix, c(1, 2), function(x) quantile(x, probs = 0.975, na.rm = TRUE))
  
  rownames(median_S_matrix) <- years
  colnames(median_S_matrix) <- age_labels
  rownames(lower_S_matrix) <- years
  colnames(lower_S_matrix) <- age_labels
  rownames(upper_S_matrix) <- years
  colnames(upper_S_matrix) <- age_labels
  
  #aggregated I matrix
  median_I_matrix <- apply(I_matrix, c(1, 2), median)
  lower_I_matrix <- apply(I_matrix, c(1, 2), function(x) quantile(x, probs = 0.025, na.rm = TRUE))
  upper_I_matrix <- apply(I_matrix, c(1, 2), function(x) quantile(x, probs = 0.975, na.rm = TRUE))
  
  rownames(median_I_matrix) <- years
  colnames(median_I_matrix) <- age_labels
  rownames(lower_I_matrix) <- years
  colnames(lower_I_matrix) <- age_labels
  rownames(upper_I_matrix) <- years
  colnames(upper_I_matrix) <- age_labels
  
  
  # output
  data[[model]] <- list(
    median_S = median_S_matrix,
    lower_S = lower_S_matrix,
    upper_S_ = upper_S_matrix,
    
    median_I = median_I_matrix,
    lower_I = lower_I_matrix,
    upper_I = upper_I_matrix
  )
}


####- data transformation and combination -####
long_format <- function(matrix, matrix_lower = NULL, matrix_upper = NULL, model_name) {
  
  data <- as.data.frame(matrix)
  data$year <- years
  data <- gather(data, key = "age", value = "value", -year)
  data$model <- model_name
  data$age <- as.numeric(gsub("yo", "", data$age)) 
  
  if (!is.null(matrix_lower)) {
    data_lower <- as.data.frame(matrix_lower)
    data_lower$year <- years
    data_lower <- gather(data_lower, key = "age", value = "value", -year)
    data$value_lower <- as.data.frame(data_lower)$value
  }
  
  if (!is.null(matrix_upper)) {
    data_upper <- as.data.frame(matrix_upper)
    data_upper$year <- years
    data_upper <- gather(data_upper, key = "age", value = "value", -year)
    data$value_upper <- as.data.frame(data_upper)$value  
  }
  
  return(data)
}

S_matrix_m1 <- long_format(data[[9]]$median_S, data[[9]]$lower_S, data[[9]]$upper_S, "Two-outbreak model")
S_matrix_m2 <- long_format(data[[10]]$median_S, data[[10]]$lower_S, data[[10]]$upper_S, "Two-outbreak model + age risk")
S_matrix_m3 <- long_format(data[[11]]$median_S, data[[11]]$lower_S, data[[11]]$upper_S, "Two-outbreak model + rho")
S_matrix_m4 <- long_format(data[[12]]$median_S, data[[12]]$lower_S, data[[12]]$upper_S, "Two-outbreak model + age risk + rho")
S_data <- bind_rows(S_matrix_m1, S_matrix_m2, S_matrix_m3, S_matrix_m4) %>% 
  mutate(model = factor(model, levels = c("Two-outbreak model", "Two-outbreak model + age risk",
                                          "Two-outbreak model + rho", "Two-outbreak model + age risk + rho")))

I_matrix_m1 <- long_format(data[[9]]$median_I, data[[9]]$lower_I, data[[9]]$upper_I, "Two-outbreak model")
I_matrix_m2 <- long_format(data[[10]]$median_I, data[[10]]$lower_I, data[[10]]$upper_I, "Two-outbreak model + age risk")
I_matrix_m3 <- long_format(data[[11]]$median_I, data[[11]]$lower_I, data[[11]]$upper_I, "Two-outbreak model + rho")
I_matrix_m4 <- long_format(data[[12]]$median_I, data[[12]]$lower_I, data[[12]]$upper_I, "Two-outbreak model + age risk + rho")
I_data <- bind_rows(I_matrix_m1, I_matrix_m2, I_matrix_m3, I_matrix_m4) %>% 
  mutate(model = factor(model, levels = c("Two-outbreak model", "Two-outbreak model + age risk",
                                          "Two-outbreak model + rho", "Two-outbreak model + age risk + rho")))


####- plotting 1: original susceptible matrix -####
# ggplot(S_data %>% filter(year %in% 1966:1973), 
#        aes(x = age, y = value, color = model)) +
#   geom_line(linewidth = 1) +
#   geom_ribbon(aes(ymin = value_lower, ymax = value_upper, fill = model), 
#               alpha = 0.2, linetype = 0) +
#   facet_wrap(~ year, scales = "fixed", ncol = 4) + 
#   scale_x_continuous(breaks = c(1:9), labels = c("≤1", as.character(2:9))) +
#   theme_minimal() +
#   labs(x = "Age (year)", y = "Proportion of susceptible 1-9 yrs children", title = "") +
#   theme_bw() +
#   theme(plot.margin = margin(1,1,1,1, unit = "cm"),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 12),
#         panel.background = element_rect(fill = "transparent"),
#         legend.position = "top",
#         axis.text = element_text(size = 12, color = "black"),
#         axis.title = element_text(size = 14, face = "bold"),
#         strip.text = element_text(size = 12),
#         panel.grid = element_blank(),
#         strip.background = element_rect(fill = "transparent"))   


####- plotting 2: Total proportion susceptible matrix -####
# ggplot(S_data %>% filter(year %in% 1966:1973) %>% 
#          group_by(year, model) %>% 
#          summarise(value_avg = sum(value, na.rm = TRUE) / 9,
#                    value_lower = sum(value_lower, na.rm = TRUE) / 9,
#                    value_upper = sum(value_upper, na.rm = TRUE) / 9), 
#        aes(x = factor(year), y = value_avg, fill = model, colour = model)) +
#   geom_point(position = position_dodge(width = 0.75), size = 2) +
#   geom_errorbar(aes(ymin = value_lower, ymax = value_upper), 
#                 position = position_dodge(width = 0.75), width = 0.5, linewidth = 0.8) +  
#   geom_text(aes(y = value_upper, label = sprintf("%.2f",round(value_avg, 2))), colour = "black",
#             position = position_dodge(width = 0.7), vjust = -1, size = 3.5) + 
#   scale_x_discrete(breaks = seq(1966, 1973, by = 1)) +
#   scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
#   coord_cartesian(ylim = c(0, 1.1), clip = "off") +
#   scale_color_manual(values = brewer.pal(6, "Set3")[c(1,3:5)]) +
#   labs(x = "Year", y = "Proportion of Susceptible 1-9 yrs children") +
#   theme(plot.margin = margin(1,0.5,1,0.5, unit = "cm"),
#         legend.title = element_blank(),
#         legend.position = "top",
#         legend.text =  element_text(size = 12),
#         panel.background = element_rect(fill = "transparent"),
#         panel.border = element_blank(),  
#         axis.text = element_text(size = 12, color = "black"),
#         axis.title = element_text(size = 14, face = "bold"),
#         panel.grid = element_blank(),
#         axis.line = element_line(color = "black")  
#   ) 

####- Sup_FigX_A: : Age-specific proportion susceptible matrix -####
S_data %>% filter(year %in% 1966:1973) %>% 
  group_by(year, model, age) %>%  
  summarise(value_avg = sum(value, na.rm = TRUE) / 9) %>% 
  left_join(S_data %>% filter(year %in% 1966:1973) %>% 
              group_by(year, model) %>% 
              summarise(value_total = sum(value, na.rm = TRUE) / 9)) %>%
  mutate(idx = paste0(year,"-", model)) %>%
  ungroup() %>% 
  mutate(x = c(rep(1:4,each = 9),
               rep(6:9,each = 9),
               rep(11:14,each = 9),
               rep(16:19,each = 9),
               rep(21:24,each = 9),
               rep(26:29,each = 9),
               rep(31:34,each = 9),
               rep(36:39,each = 9))) -> S_data1

ggplot(S_data1, aes(x = x, y = value_avg, fill = factor(age, levels = seq(9,1)))) +  
  geom_bar(stat = "identity", position = "stack", width = 0.95, color = "black", alpha = 0.7) + 
  geom_text(aes(label = sprintf("%.1f",value_avg * 9 * 100)),  
            position = position_stack(vjust = 0.5), size = 3.5) +
  geom_text(aes(y = value_total+0.05, label = sprintf("%.1f",value_total* 100)),  
            size = 4, fontface = "bold") +  
  scale_fill_manual(
    values = brewer.pal(9, "Set3"),
    labels = paste0(rev(unique(S_data$age)), " yo"),
    breaks = seq(9, 1)
  ) +
  annotate("text", x = seq(2.5, 38.5, by = 5), y = -0.08, label = seq(1966, 1973), size = 5.5) +
  annotate("text", x = unique(S_data1$x), y = -0.03, label = rep(c("m1", "m2"," m3", "m4"), 8), size = 5.5) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.5,1), label = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(label = NULL) + 
  coord_cartesian(ylim = c(0, 1.02), clip = "off") +  
  labs(x = " ", title = "A",
       y = "Susceptible Proportion (%)", 
       fill = "Age Group")+
  # caption ="m1: Original   m2: Accouting for age indedpendent FOI    m3: Accouting for seroreversion    m4: Accouting for age indedpendent FOI and seroreversion") +
  guides(fill = guide_legend(nrow = 1)) + 
  theme(
    plot.margin = margin(0.5, 1, 1, 1, unit = "cm"),
    plot.title = element_text(size = 20),
    legend.title = element_blank(),
    legend.position = c(0.65,0.95),
    legend.text =  element_text(size = 16),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_blank(),  
    axis.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), 
    axis.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(vjust = -5),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.caption = element_text(size = 16, hjust = 0, margin =  margin(t = 20, b = -5)) 
  ) -> figX_A



####- Sup_FigX_B: annual incidence -####
ggplot(I_data %>% filter(year %in% 1966:1973) %>% 
         group_by(year, model) %>% 
         summarise(value_avg = sum(value, na.rm = TRUE) / 9,
                   value_lower = sum(value_lower, na.rm = TRUE) / 9,
                   value_upper = sum(value_upper, na.rm = TRUE) / 9), 
       aes(x = factor(year), y = value_avg, fill = model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), 
           width = 0.75, size = 1, colour = "black", alpha = 0.7) +
  geom_errorbar(aes(ymin = value_lower, ymax = value_upper), 
                position = position_dodge(width = 0.75), size = 0.8, width = 0.3) +  
  geom_text(aes(y = value_upper, label = sprintf("%.2f",value_avg)), colour = "black",
            position = position_dodge(width = 0.7), vjust = -1, size = 4.5) + 
  scale_x_discrete(breaks = seq(1966, 1973, by = 1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  coord_cartesian(ylim = c(0, 0.4), clip = "off") +
  scale_fill_manual(values = brewer.pal(6, "Set3")[c(1,3:5)],
                    labels = c("Two-outbreak model",  "Plus age-dedpendent FOI", 
                               "Plus seroreversion", "Plus age-dedpendent FOI and seroreversion")) +
  labs(x = "Year", y = "Annual incidence", title = "B") +
  guides(fill = guide_legend(nrow = 4),
         color = guide_legend(nrow = 4)) + 
  theme(plot.margin = margin(0.5,1,0.5,1, unit = "cm"),
        plot.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = c(0.8,0.85),
        legend.text =  element_text(size = 16),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),  
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black")  
  ) -> figX_B




####- Sup_FigX_C : annual incidence by age group in specific year -####
I_data %>% filter(year %in% c(1968,1969)) %>% 
  group_by(year, model, age) %>%  
  summarise(value_avg = sum(value, na.rm = TRUE) / 9) %>% 
  left_join(I_data %>% filter(year %in% c(1968,1969)) %>% 
              group_by(year, model) %>% 
              summarise(value_total = sum(value, na.rm = TRUE) / 9)) %>%
  mutate(idx = paste0(year,"-", model)) %>%
  ungroup() %>% 
  mutate(x = c(rep(1:4,each = 9),
               rep(6:9,each = 9))) -> I_data1
I_data1$label <- I_data1$value_avg
I_data1$label <- as.character(sprintf("%.2f",I_data1$label * 9))
I_data1$label[I_data1$model == "Two-outbreak model" & I_data1$year == 1968] <- " "


ggplot(I_data1, aes(x = x, y = value_avg, fill = factor(age, levels = seq(9,1)))) +  
  geom_bar(stat = "identity", position = "stack", width = 0.95, color = "black", alpha = 0.7) + 
  geom_text(aes(label = label),  
            position = position_stack(vjust = 0.5), size = 5) +
  geom_text(aes(y = value_total+0.02, label = sprintf("%.2f",value_total)),  
            size = 5, fontface = "bold") +  
  scale_fill_manual(
    values = brewer.pal(9, "Set3"),
    labels = paste0(rev(unique(S_data$age)), " yo"),
    breaks = seq(9, 1)
  ) +
  annotate("text", x = c(2.5, 7.5), y = 0.3, label = seq(1968, 1969), size = 6.5) +
  annotate("text", x = unique(I_data1$x), y = -0.01, 
           label = rep(c("m1", "m2"," m3", "m4"), 2),
           size = 5.5) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.5,0.4), ) +
  scale_x_continuous(label = NULL) + 
  coord_cartesian(ylim = c(0, 0.402), clip = "off") +  
  labs(title = "C",
       x = " ", 
       y = "Annual incidence", 
       fill = "Age Group")+
  # caption ="m1: Original   m2: Accouting for age indedpendent FOI    m3: Accouting for seroreversion    m4: Accouting for age indedpendent FOI and seroreversion") +
  guides(fill = guide_legend(nrow = 2)) + 
  theme(
    plot.margin = margin(0.5, 1, 0.5, 1, unit = "cm"),
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.position = c(0.5,0.95),
    legend.text = element_text(size = 14),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_blank(),  
    axis.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(vjust = -5),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.caption = element_text(size = 12, hjust = 0, margin =  margin(t = 20, b = -5)) 
  )  -> figX_C




####- Sup_FigX_D: mean age of infection -####
infect_age <- mean_age %>%
  filter(year %in% 1966:1973) %>%
  group_by(year, model, sim) %>%
  summarise(
    mean_age = sum(age * value, na.rm = TRUE) / sum(value, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(year, model) %>%
  summarise(meanage = median(mean_age),
            meanage_lower = quantile(mean_age, probs = 0.025),
            meanage_upper = quantile(mean_age, probs = 0.975))


ggplot(infect_age, aes(x = factor(year), y = meanage, 
                       group = as.character(model),
                       color = as.character(model), 
                       fill = as.character(model))) +
  geom_ribbon(aes(ymin = meanage_lower, ymax = meanage_upper), alpha = 0.3, color = NA) +
  geom_line(size = 1.25) +
  geom_point(size = 2.5) +
  scale_color_manual(values = brewer.pal(9, "Set1")[c(1,3,4,5)],
                     label = c("Two-outbreak model",  "Plus age-dedpendent FOI", 
                               "Plus seroreversion", "Plus age-dedpendent FOI and seroreversion")) +
  scale_fill_manual(values = brewer.pal(9, "Set1")[c(1,3,4,5)],
                    label = c("Two-outbreak model",  "Plus age-dedpendent FOI", 
                              "Plus seroreversion", "Plus age-dedpendent FOI and seroreversion")) +
  scale_x_discrete(breaks = seq(1966, 1973, by = 1), expand = c(0,0.1)) +
  scale_y_continuous(expand = c(0,0), limits = c(4, 6)) +
  labs(x = "Year", 
       y = "Mean age of Infection",
       title = "D") +
  guides(fill = guide_legend(nrow = 4),
         color = guide_legend(nrow = 4)) + 
  theme(plot.margin = margin(0.5,1,0.5,1, unit = "cm"),
        plot.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = c(0.6,0.95),
        legend.text = element_text(size = 16),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),  
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black")  
  ) -> figX_D




tiff("../result/sup_figX_model_comparison.tiff", width = 40, height = 42, units = "cm", compression = "lzw", res = 300)

figX_A /figX_B / (figX_C + figX_D)

dev.off()

