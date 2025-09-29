rm(list = ls())

library(Rsero)
library(rstan)
library(loo)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

matrix_ind <- read.csv("I_data_independent_model.csv")
matrix_out <- read.csv("I_data_outbreak_model.csv")
dt <- rbind(
  matrix_ind %>% filter(model %in% c("Independent model + age risk", 
                                     "Independent model + rho")) %>%
    filter(year %in% c(1967:1971)) %>%
    filter(age %in% c(2,9)),
  
  matrix_out %>% filter(model %in% c("Two-outbreak model + age risk", 
                                     "Two-outbreak model + rho")) %>%
    filter(year %in% c(1967:1971)) %>%
    filter(age %in% c(2,9))
)

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
dt$age <- factor(dt$age, levels = c(2, 9), labels = c("2 yo", "9 yo"))
dt$model <- factor(dt$model, levels = c("Independent model + age risk", "Independent model + rho",
                                        "Two-outbreak model + age risk","Two-outbreak model + rho"), 
                   labels = c("Model2: Independent model with age-dependent FOI",
                              "Model3: Independent model with seroreversion",
                              "Model10: Two-outbreak model with age-dependent FOI", 
                              "Model11: Two-outbreak model with seroreversion"))

ggplot(dt, aes(x = year, y = value, color = age, fill = age)) +
  geom_ribbon(aes(ymin = value_lower, ymax = value_upper), 
              alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 1) +
  geom_point(size = 2, show.legend = F) +
  facet_wrap(~ model, ncol = 2) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(x = "Year",
       y = "Infection probability",
       color = "Age Group") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_blank(),  
    axis.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  ) -> p1

# ggsave("../result/infection_probability_by_model.tiff", plot = p1,
#        width = 8, height = 6, units = "in", dpi = 300, compression = "lzw", bg = "white")

