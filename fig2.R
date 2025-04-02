

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
years <- 1:8
models <- c(m1, m2, m3, m4)
fitted_seropre <- data.frame()
observed_seropre <- data.frame()
for (i in models) {
  for (j in years) {
    #fitted data
    fitted_seropre_tmp <- data.frame(
      age = model_plot[[i]][[j]]$layers[[1]]$data$x[1:9],
      median = model_plot[[i]][[j]]$layers[[2]]$data$y,
      upper = model_plot[[i]][[j]]$layers[[1]]$data$y[1:9],
      lower = rev(model_plot[[i]][[j]]$layers[[1]]$data$y[10:18]),
      year = j,
      model = i
    )
    fitted_seropre <- rbind(fitted_seropre, fitted_seropre_tmp)
    # observed data
    observed_seropre_tmp <- data.frame(model_plot[[i]][[j]]$layers[[3]]$data)
    observed_seropre_tmp$year <- j
    observed_seropre_tmp$model <- i
    observed_seropre <- rbind(observed_seropre, observed_seropre_tmp)
  }
}
observed_seropre$age[observed_seropre$age == 9] <- 8.5

  
  
# color = brewer.pal(6, "Set3")

ggplot() +
  facet_wrap(~year, ncol = 4, scales = "free_x",
             labeller = labeller(year = c("1" = "1966", "2" = "1967", "3" = "1968", "4" = "1969", 
                                          "5" = "1970", "6" = "1971", "7" = "1972", "8" = "1973"))) +
  geom_point(data = observed_seropre %>% filter(model %in% c(2,3)),
             aes(x = age-0.5, y = mean, fill = factor(model)), size= 1)+
  geom_linerange(data = observed_seropre %>% filter(model %in% c(2,3)),
                 aes(x = age-0.5, ymin = lower, ymax = upper), color = "black",size = 0.3) +
  geom_ribbon(data = fitted_seropre %>% filter(model %in% c(2,3)),
              aes(x = age, ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.3)+
  geom_line(data = fitted_seropre %>% filter(model %in% c(2,3)),
            aes(x = age, y = median, color = factor(model)), size = 1) + 

  scale_fill_manual(name = "", values = c("#377EB8", "#E41A1C"), 
                    labels = c("Accounting for age-dependent FOI", "Accounting for seroreversion")) + 
  scale_color_manual(name = "", values = c("#377EB8", "#E41A1C"), 
                     labels = c("Accounting for age-dependent FOI", "Accounting for seroreversion")) +
  scale_y_continuous("Seroprevalence",limits = c(0,1)) +
  scale_x_continuous("Age(year)", breaks = c(1:8), labels = c("≤1", 2:7, "8-10")) + 
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
        axis.text.y = element_text(size = 12)
  ) -> fig2_A
              
  

ggplot() +
  facet_wrap(~year, ncol = 4, scales = "free_x",
             labeller = labeller(year = c("1" = "1966", "2" = "1967", "3" = "1968", "4" = "1969", 
                                          "5" = "1970", "6" = "1971", "7" = "1972", "8" = "1973"))) +
  geom_point(data = observed_seropre %>% filter(model %in% c(10,11)),
             aes(x = age-0.5, y = mean, fill = factor(model)), size= 1)+
  geom_linerange(data = observed_seropre %>% filter(model %in% c(10,11)),
                 aes(x = age-0.5, ymin = lower, ymax = upper), color = "black",size = 0.3) +
  geom_ribbon(data = fitted_seropre %>% filter(model %in% c(10,11)),
              aes(x = age, ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.3)+
  geom_line(data = fitted_seropre %>% filter(model %in% c(10,11)),
            aes(x = age, y = median, color = factor(model)), size = 1) + 
  scale_fill_manual(name = "", values = c("#377EB8", "#E41A1C"), 
                    labels = c("Accounting for age-dependent FOI", "Accounting for seroreversion")) + 
  scale_color_manual(name = "", values = c("#377EB8", "#E41A1C"), 
                     labels = c("Accounting for age-dependent FOI", "Accounting for seroreversion")) +
  scale_y_continuous("Seroprevalence",limits = c(0,1)) +
  scale_x_continuous("Age(year)", breaks = c(1:8), labels = c("≤1", 2:7, "8-10")) + 
  labs(title = "B. Two-outbreak model") +
  theme(plot.margin = margin(0.1,1,0.5,1, unit = "cm"),
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
        axis.text.y = element_text(size = 12)
  ) -> fig2_B




# tiff("../result/fig2.tiff", width = 30, height = 30, units = "cm", compression = "lzw", res = 300)
# 
# fig2_A / fig2_B
# 
# dev.off()

