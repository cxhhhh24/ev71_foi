rm(list = ls())

library(Rsero)
library(ggplot2)
library(dplyr)
library(readxl)
library(patchwork)
library(ggpubr)
library(rcartocolor)
library(RColorBrewer)
library(wesanderson)
library(rstan)
library(loo)
library(posterior)




model_fit <- readRDS("model_fit_12model_0326.rds")
model_plot <- readRDS("model_plot_12model_0326.rds")

my_colors <- c(
  "#377EB8",  
  "#E41A1C",  
  "#4DAF4A",  
  "#984EA3",  
  "#FF7F00",  
  "#FFFF33",  
  "#A65628",  
  "#F781BF",  
  "#999999"   
)


#### - selecting models -####
m1 <- 2 #independent with age
m2 <- 3 #independent with rho
m3 <- 10 #two-outbreak with age
m4 <- 11 #two-outbreak with rho


#### - fitted seroprevalence -####
models <- c(m1, m2)
foi1 <- data.frame()

for (i in models) {
  
  foi_tmp <- data.frame(
    year = 1958:1973,
    median = rev(parameters_credible_intervals(model_fit[[i]])[,2][1:16]),
    upper = rev(parameters_credible_intervals(model_fit[[i]])[,3][1:16]),
    lower = rev(parameters_credible_intervals(model_fit[[i]])[,1][1:16]),
    model = i
  )
  foi1 <- rbind(foi1, foi_tmp)
  
}

models <- c(m3, m4)
foi2 <- data.frame()

for (i in models) {
  
  posterior_sample <- rstan::extract(model_fit[[i]][[1]], permuted = FALSE)
  foi_tmp <- matrix(NA, nrow = 10000, ncol = 16)
  
  for (j in 1:16) {
    
    foi_tmp[, j] <- as.vector(posterior_sample[, , paste0("lambda[", 16 - (j - 1), "]")])
    
  }
  
  foi_tmp <- data.frame(
    year = 1958:1973,
    median = apply(foi_tmp, 2, function(x) quantile(x, prob = 0.5)),
    upper = apply(foi_tmp, 2, function(x) quantile(x, prob = 0.975)),
    lower = apply(foi_tmp, 2, function(x) quantile(x, prob = 0.025)),
    model = i)
  
  foi2 <- rbind(foi2, foi_tmp)
}

foi <- rbind(foi1, foi2)



# foi_output <- foi %>% mutate(foi = paste0(sprintf("%.2f", median), " (", sprintf("%.2f", lower), "-", sprintf("%.2f", upper), ")"))
# write.csv(foi_output, "../result/foi_output.csv")

####- FOI plotting -####
ggplot(foi %>% filter(model %in% c(2,3)),
       aes(x = year, y = median, color = factor(model))) +
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.3, color = NA) + 
  geom_line(size = 0.75)  +
  scale_fill_manual(name = "", values = my_colors[1:2], 
                    labels = c("Accounting for age-dependent FOI", "Accounting for seroreversion")) + 
  scale_color_manual(name = "", values = my_colors[1:2], 
                     labels = c("Accounting for age-dependent FOI", "Accounting for seroreversion")) +
  scale_y_continuous("Estimated force of infection", limits = c(0,0.6), expand = c(0,0.005)) +
  scale_x_continuous("Year",  breaks = 1958:1973, labels = 1958:1973) + 
  labs(title = "A. Independent model") +
  theme(plot.margin = margin(0.5,1,0.1,1, unit = "cm"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),  
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12)) -> fig2_A



ggplot(foi %>% filter(model %in% c(10,11)),
       aes(x = year, y = median, color = factor(model))) +
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.3, color = NA) + 
  geom_line(size = 1, alpha = .75)  +
  scale_fill_manual(name = "", values = my_colors[3:4], 
                    labels = c("Accounting for age-dependent FOI", "Accounting for seroreversion")) + 
  scale_color_manual(name = "", values = my_colors[3:4], 
                     labels = c("Accounting for age-dependent FOI", "Accounting for seroreversion")) +
  scale_y_continuous("Estimated force of infection",limits = c(0,0.6), expand = c(0,0.005)) +
  scale_x_continuous("Year",  breaks = 1958:1973, labels = 1958:1973) + 
  labs(title = "B. Two-outbreak model") +
  theme(plot.margin = margin(0.5,1,0.1,1, unit = "cm"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),  
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12)) -> fig2_B



# tiff("../result/fig2.tiff", width = 30, height = 24, units = "cm", compression = "lzw", res = 300)
# 
# fig2_A / fig2_B
# 
# dev.off()

