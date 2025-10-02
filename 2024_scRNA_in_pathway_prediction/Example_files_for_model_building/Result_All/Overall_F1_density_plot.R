filepath <- getwd()
setwd(filepath)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)
args = commandArgs(TRUE)
SimulationCSV <- args[1]
value <- args[2]
module <- args[3]

Refervalue <- read.csv(SimulationCSV, header=TRUE, sep=",")
value <- read.csv(value, header = T)

Mean <- mean(value$Global_F1_Score)
SD <- sd(value$Global_F1_Score)
if(mean(Refervalue[[1]]) == 0)  { zscore <- 0}else{zscore <- (Mean-mean(Refervalue[[1]]))/sd(Refervalue[[1]])}
if(zscore > 0)  { pvalue <- pnorm(q = zscore, lower.tail = FALSE) }else{
  pvalue <- pnorm(q = zscore, lower.tail = TRUE)}
zscore_value <- round(zscore,3)
pval <- format(pvalue, scientific = TRUE, digits = 3)

den <- density(Refervalue[[1]])
Max<-max(den$y)*1.1
x_lim = max(value$Global_F1_Score*1.2)
x_start = Mean - SD
x_end = Mean + SD
data_change <- ggplot(Refervalue, aes(x = Refervalue[[1]])) + 
  geom_density(color = "#69b3a2", lwd = 0.6, linetype = 1, fill="#69b3a2", adjust=1.75, alpha=0.5) +
  labs(x = 'F1_score', y = "Density") +
  geom_rect(xmin=x_start, xmax=x_end, ymin=-Inf, ymax=Inf, fill="#E9DCDB", alpha=0.02, linewidth = 0) +
  geom_vline(xintercept = Mean, linetype = 1, color="red", linewidth = 0.6) +
  geom_vline(xintercept = mean(Refervalue[[1]]), linetype = 2, color="grey", linewidth = 0.6) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlim(0, x_lim) +
  ylim(0, Max)+
  annotate("text", x = x_lim*0.9, y = Max, label = format(paste0("z=", zscore_value)), size = 3, col = "black") +
  annotate("text", x = x_lim*0.9, y = Max*0.8, label = format(paste0("p=", pval)), size = 3, col = "black")

p1 <- data_change
ggsave(p1, filename = paste("Overall_F1_densityplot_", paste(module, collapse="_"), ".pdf", sep=""), width = 6, height = 4)
