filepath <- getwd()
setwd(filepath)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)
args = commandArgs(TRUE)
GeneType <- args[1]
value <- args[2]
module <- args[3]

Refervalue <- read.csv(GeneType, header=TRUE, sep=",")
value <- read.csv(value, header = T)

cols <- colnames(Refervalue)[1:14]
pvalue_list <- c()
num <- length(Refervalue)
for (i in 1:length(cols)) {
  count_gt <- sum(Refervalue[[cols[i]]] >= value$Refervalue[i], na.rm = TRUE)
  p <- (count_gt + 1) / (num + 1)
  pvalue_list <- c(pvalue_list, format(p, scientific = TRUE, digits = 3))
  p <- paste0('p', cols[i])
  den <- density(Refervalue[[cols[i]]])
  Max<-max(den$y)*1.1
  x_lim = max(value$Refervalue*1.3)
  data_change <- ggplot(Refervalue, aes(x = !!sym(cols[i]))) + 
    geom_density(color = "#69b3a2", lwd = 0.6, linetype = 1, fill="#69b3a2", adjust=1.75, alpha=0.5) +
    labs(x = cols[i], y = "Density") +
    geom_vline(xintercept = value$Refervalue[i], linetype = 1, color="red", linewidth = 0.6) +
    geom_vline(xintercept = mean(Refervalue[[cols[i]]]), linetype = 2, color="grey", linewidth = 0.6) +
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    xlim(0, x_lim) +
    annotate("text", x = x_lim*0.9, y = Max*0.8, label = format(paste0("p=", pvalue_list[i])), size = 3, col = "black")
  assign(p, data_change)
}

p1 <- pAmines.and.Polyamines+pAmino.Acids+pCarbohydrates+pCofactors+pDetoxification+pEnergy.Metabolism+pFatty.Acids.and.Lipids+pHormones+pInorganic.Nutrients+pIntermediate.Metabolism+pNucleotides+pRedox+pSpecialized.Metabolism+pOther+plot_layout(nrow = 7)
print(p1)
ggsave(p1, filename = paste("Simulation_C_densityplot_", paste(module, collapse="_"), ".pdf", sep=""), width = 9, height = 12)
