rm(list=ls())

library(Rsero)
library(ggplot2)
library(dplyr)
library(readxl)
library(RColorBrewer)
library(patchwork)



read.csv(dat_indi, "../data/dat_indi_clean.csv") -> dat
read.csv("dat_aggragate.csv") -> dat

#####----- seroprevalence plot -----#######

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


dat$pre <- dat$p / dat$n
dat$Year <- factor(dat$Year)

# purple_palette <- brewer.pal(n = 9, name = "Purples")  


ggplot(data = dat, aes(x = Age, y = pre, group = Year, color = Year)) + 
  geom_point(size = 2) +
  geom_line(size = 1) +
  # geom_linerange(aes(x = age_jitter, ymin = lower, ymax = upper, group = year))
  scale_color_manual(name = "Sampling year", values = my_colors) +
  guides(color = guide_legend(ncol = 3))+
  scale_y_continuous(expand = c(0,0), limits = c(-0.02, 1), name = "Seroprevalence")+
  scale_x_continuous(name = "Age (years)", breaks = c(1:8), labels = c("≤1",2:7,"8-10")) +
  theme(
    panel.background = element_blank(),
    plot.margin = unit(c(1,1,1,1), "cm"),
    axis.line = element_line(colour = "black", size = 0.75),
    axis.ticks = element_line(colour = "black", size = 0.75),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 16),
    legend.position = c(0.7,0.9),
    legend.text = element_text(colour = "black", size = 8),
    legend.title = element_text(colour = "black", size = 8)
  ) -> fig1

quantile(dat$pre[dat$Age == 8], probs = c(0.025, 0.5, 0.975)) 
quantile(dat$pre[dat$Age %in% c(1,2,3)], probs = c(0.025, 0.5, 0.975)) 


dat %>% group_by(Age) %>% summarise(n = sum(n), p = sum(p)) %>% 
  mutate(pre = p/n) %>%
  rowwise() %>%
  mutate(lci = binom.test(p,n)$conf.int[1],
         uci = binom.test(p,n)$conf.int[2]) -> dat_age

ggplot(data = dat_age) + 
  geom_point(aes(x = Age, y = pre), color = my_colors[1], size = 2)+
  geom_linerange(aes(x = Age, ymin = lci, ymax = uci), color = my_colors[1],size = 1) +
  scale_x_continuous(name = "Age (years)", breaks = c(1:8), labels = c("≤1",2:7,"8-10"))+
  scale_y_continuous(expand = c(0,0), limits = c(-0.02, 0.5), name = "Seroprevalence") +
  theme(
    panel.background = element_blank(),
    plot.margin = unit(c(1,1,1,1), "cm"),
    axis.line = element_line(colour = "black", size = 0.75),
    axis.ticks = element_line(colour = "black", size = 0.75),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 16),
    legend.position = c(0.8,0.85),
    legend.text = element_text(colour = "black", size = 12))  -> fig2




# 
# tiff("../result/fig1_cd.tiff", width = 26, height = 10, units = "cm", compression = "lzw", res = 300)
# fig1 + fig2
# dev.off()

