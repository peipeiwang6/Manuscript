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
value <- value[, c(1, 12, 13)]

cols <- colnames(Refervalue)[1:14]
zscore_list <- c()
pvalue_list <- c()
for (i in 1:length(cols)) {
  if(mean(Refervalue[[cols[i]]]) == 0)  { zscore <- data.frame(0) }else{
    zscore <- data.frame(value$Mean[i]-mean(Refervalue[[cols[i]]]))/sd(Refervalue[[cols[i]]])}
  zscore_list <- c(zscore_list, round(zscore[[1]],3))
  if(zscore_list[i] > 0)  { pvalue <- pnorm(q = zscore_list[i],lower.tail = FALSE) }else{
    pvalue <- pnorm(q = zscore[[1]],lower.tail = TRUE)}
  pvalue_list <- c(pvalue_list, format(pvalue, scientific = TRUE, digits = 3))
  p <- paste0('p', cols[i])
  den <- density(Refervalue[[cols[i]]])
  Max<-max(den$y)*0.75
  x_lim = max(value$Mean*1.3)
  x_start = value$Mean[i] - value$SD[i]
  x_end = value$Mean[i] + value$SD[i]
  data_change <- ggplot(Refervalue, aes(x = !!sym(cols[i]))) + 
    geom_density(color = "#69b3a2", lwd = 0.6, linetype = 1, fill="#69b3a2", adjust=1.75, alpha=0.5) +
    labs(x = cols[i], y = "Density") +
    geom_rect(xmin=x_start, xmax=x_end, ymin=-Inf, ymax=Inf, fill="#E9DCDB", alpha=0.02, linewidth = 0) +
    geom_vline(xintercept = value$Mean[i], linetype = 1, color="red", linewidth = 0.6) +
    geom_vline(xintercept = mean(Refervalue[[cols[i]]]), linetype = 2, color="grey", linewidth = 0.6) +
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    xlim(0, x_lim) +
    ylim(0, max(den$y))+
    annotate("text", x = x_lim*0.9, y = Max, label = format(paste0("z=", zscore_list[i])), size = 3, col = "black") +
    annotate("text", x = x_lim*0.9, y = Max*0.8, label = format(paste0("p=", pvalue_list[i])), size = 3, col = "black")
  assign(p, data_change)
}

p1 <- pAmines.and.Polyamines+pAmino.Acids+pCarbohydrates+pCofactors+pDetoxification+pEnergy.Metabolism+pFatty.Acids.and.Lipids+pHormones+pInorganic.Nutrients+pIntermediate.Metabolism+pNucleotides+pRedox+pSpecialized.Metabolism+pOther+plot_layout(nrow = 7)
print(p1)
ggsave(p1, filename = paste("Simulation_F1_densityplot_", paste(module, collapse="_"), ".pdf", sep=""), width = 9, height = 12)
dev.off()
