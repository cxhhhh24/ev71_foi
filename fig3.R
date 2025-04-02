
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

#### - selecting models -####
m1 <- 2 #independent with age
m2 <- 3 #independent with rho
m3 <- 10 #two-outbreak with age
m4 <- 11 #two-outbreak with rho


#### - fitted seroprevalence -####
models <- c(m1, m2, m3, m4)
foi <- data.frame()

for (i in models) {
  
    foi_tmp <- data.frame(
      year = plot(model_fit[[i]])[[1]]$layers[[1]]$data$x[1:16],
      median = rev(plot(model_fit[[i]])[[1]]$layers[[2]]$data$y[1:16]),
      upper = plot(model_fit[[i]])[[1]]$layers[[1]]$data$y[1:16],
      lower = rev(plot(model_fit[[i]])[[1]]$layers[[1]]$data$y[17:32]),
      model = i
    )
    foi <- rbind(foi, foi_tmp)
   
}

####- FOI plotting -####
ggplot(foi %>% filter(model %in% c(2,3)),
       aes(x = year, y = median, color = factor(model))) +
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.3, color = NA) + 
  geom_line(size = 1, alpha = 0.5)  +
  scale_fill_manual(name = "", values = c("#377EB8", "#E41A1C"), 
                    labels = c("Accounting for age-dependent FOI", "Accounting for seroreversion")) + 
  scale_color_manual(name = "", values = c("#377EB8", "#E41A1C"), 
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
        axis.text.y = element_text(size = 12)) -> fig3_A



ggplot(foi %>% filter(model %in% c(10,11)),
       aes(x = year, y = median, color = factor(model))) +
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.3, color = NA) + 
  geom_line(size = 1, alpha = 0.5)  +
  scale_fill_manual(name = "", values = c("#377EB8", "#E41A1C"), 
                    labels = c("Accounting for age-dependent FOI", "Accounting for seroreversion")) + 
  scale_color_manual(name = "", values = c("#377EB8", "#E41A1C"), 
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
        axis.text.y = element_text(size = 12)) -> fig3_B



# tiff("../result/fig3.tiff", width = 30, height = 24, units = "cm", compression = "lzw", res = 300)
# 
# fig3_A / fig3_B
# 
# dev.off()





#####---- FOI with individual trajectory  ---#####
plot(model_fit[[m1]], individual_samples = 20)[[1]] +       
  scale_y_continuous("Estimated force of infection",limits = c(0,0.6), expand = c(0,0.005)) +
  scale_x_continuous(breaks = seq(1958,1973,2), labels = seq(1958,1973,2)) + 
  labs(title = "A. Independent model accounting for age-dependent FOI") -> foi_trajectory_m1


plot(model_fit[[m2]], individual_samples = 20)[[1]] +       
  scale_y_continuous("Estimated force of infection",limits = c(0,0.6), expand = c(0,0.005)) +
  scale_x_continuous(breaks = seq(1958,1973,2), labels = seq(1958,1973,2)) + 
  labs(title = "B. Independent model accounting for seroreversion") -> foi_trajectory_m2

plot(model_fit[[m3]], individual_samples = 20)[[1]] +       
  scale_y_continuous("Estimated force of infection",limits = c(0,0.6), expand = c(0,0.005)) +
  scale_x_continuous(breaks = seq(1958,1973,2), labels = seq(1958,1973,2)) + 
  labs(title = "C. Two-outbreak model accounting for age-dependent FOI") -> foi_trajectory_m3

plot(model_fit[[m4]], individual_samples = 20)[[1]] +       
  scale_y_continuous("Estimated force of infection",limits = c(0,0.6), expand = c(0,0.005)) +
  scale_x_continuous(breaks = seq(1958,1973,2), labels = seq(1958,1973,2)) + 
  labs(title = "D. Two-outbreak model accounting for seroreversion") -> foi_trajectory_m4



# tiff("../result/foi_trajectory.tiff", width = 50, height = 24, units = "cm", compression = "lzw", res = 300)
# 
# (foi_trajectory_m1 + foi_trajectory_m2) / (foi_trajectory_m3 + foi_trajectory_m4)
# 
# dev.off()

