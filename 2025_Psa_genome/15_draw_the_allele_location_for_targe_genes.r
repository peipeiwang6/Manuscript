#install.packages("stringr")
library(stringr)
args = commandArgs(TRUE)
#setwd('YourWorkDirectory')

gene <- args[1]

alle <- read.table('Psa_syri_allelic_fillout_using_McScanX_and_inversion_20250328_combine_tandem_duplicates.txt',head=T,sep='\t')
alle <- cbind(alle, 'Allele_numbers' = 0)
for(i in 1:nrow(alle)) alle[i,5] <- sum(alle[i,1:4]!= '')

gff <- read.table('Psa_GFF_only_for_genes.gff3',head=F,sep='\t')

Find_loc_for_gene_in_alle <- function(gene){
	for(i in 1:4){
		loci <- alle[str_detect(alle[,i], gene),]
		if(nrow(loci) > 0) break
		}
	return(loci)
	}
Alle_location <- function(gene){
	chr = substr(gene, 5,6)
	if(chr %in% c('02','06','08','09','13','18')) number <- 3 else number <- 4
	loci <- Find_loc_for_gene_in_alle(gene)
	left <- alle[1:(as.numeric(rownames(loci))-1),]
	left <- left[left[,5]==number,]
	left_loci <- as.numeric(rownames(left)[nrow(left)])
	right <- alle[(as.numeric(rownames(loci))+1):nrow(alle),]
	right <- right[right[,5]==number,]
	right_loci <- as.numeric(rownames(right)[1])
	suballe <- alle[left_loci:right_loci,]
	gff1 <- gff[gff[,9] >= min(c(strsplit(suballe[1,1],',')[[1]],strsplit(suballe[nrow(suballe),1],',')[[1]])) & gff[,9] <= max(c(strsplit(suballe[1,1],',')[[1]],strsplit(suballe[nrow(suballe),1],',')[[1]])),]
	gff2 <- gff[gff[,9] >= min(c(strsplit(suballe[1,2],',')[[1]],strsplit(suballe[nrow(suballe),2],',')[[1]])) & gff[,9] <= max(c(strsplit(suballe[1,2],',')[[1]],strsplit(suballe[nrow(suballe),2],',')[[1]])),]
	gff3 <- gff[gff[,9] >= min(c(strsplit(suballe[1,3],',')[[1]],strsplit(suballe[nrow(suballe),3],',')[[1]])) & gff[,9] <= max(c(strsplit(suballe[1,3],',')[[1]],strsplit(suballe[nrow(suballe),3],',')[[1]])),]
	if(number == 4) gff4 <- gff[gff[,9] >= min(c(strsplit(suballe[1,4],',')[[1]],strsplit(suballe[nrow(suballe),4],',')[[1]])) & gff[,9] <= max(c(strsplit(suballe[1,4],',')[[1]],strsplit(suballe[nrow(suballe),4],',')[[1]])),]

	if(number == 4) {
		max_x <- max(gff1[nrow(gff1),5]-gff1[1,4], gff2[nrow(gff2),5]-gff2[1,4], gff3[nrow(gff3),5]-gff3[1,4], gff4[nrow(gff4),5]-gff4[1,4])
		}	else {
		max_x <- max(gff1[nrow(gff1),5]-gff1[1,4], gff2[nrow(gff2),5]-gff2[1,4], gff3[nrow(gff3),5]-gff3[1,4])
		}
	h_distance <- max_x%/%12.5
	pdf(paste('Allele_location/Psa_alle_location_',gene,'.pdf',sep=''),width=10, height=4)
	plot.new()
	plot.window(xlim = c(-1000, max_x + 1000), ylim = c(-h_distance, 5*h_distance))
	segments(0, 0, 10000, 0) # draw the bar
	text(5000, -h_distance/3, '10 Kb', cex=0.3)
	segments(0, 4*h_distance, gff1[nrow(gff1),5]-gff1[1,4], 4*h_distance)
	for(i in 1:nrow(gff1)) {
		if(gff1[i,7] == '+') color = 'red' else color = 'blue'
		rect(gff1[i,4] - gff1[1,4], 4*h_distance-h_distance/6, gff1[i,5] - gff1[1,4], 4*h_distance+h_distance/6, col = color)
		text(mean(c(gff1[i,4] - gff1[1,4], gff1[i,5] - gff1[1,4])), 4*h_distance-h_distance/3, gff1[i,9], cex=0.3)
		}
	segments(0, 3*h_distance, gff2[nrow(gff2),5]-gff2[1,4], 3*h_distance)
	for(i in 1:nrow(gff2)) {
		if(gff2[i,7] == '+') color = 'red' else color = 'blue'
		rect(gff2[i,4] - gff2[1,4], 3*h_distance-h_distance/6, gff2[i,5] - gff2[1,4], 3*h_distance+h_distance/6, col = color)
		text(mean(c(gff2[i,4] - gff2[1,4], gff2[i,5] - gff2[1,4])), 3*h_distance-h_distance/3, gff2[i,9], cex=0.3)
		}
	segments(0, 2*h_distance, gff3[nrow(gff3),5]-gff3[1,4], 2*h_distance)
	for(i in 1:nrow(gff3)) {
		if(gff3[i,7] == '+') color = 'red' else color = 'blue'
		rect(gff3[i,4] - gff3[1,4], 2*h_distance-h_distance/6, gff3[i,5] - gff3[1,4], 2*h_distance+h_distance/6, col = color)
		text(mean(c(gff3[i,4] - gff3[1,4], gff3[i,5] - gff3[1,4])), 2*h_distance-h_distance/3, gff3[i,9], cex=0.3)
		}
	if(number == 4)  {
		segments(0, h_distance, gff4[nrow(gff4),5]-gff4[1,4], h_distance)
		for(i in 1:nrow(gff4)) {
			if(gff4[i,7] == '+') color = 'red' else color = 'blue'
			rect(gff4[i,4] - gff4[1,4], h_distance-h_distance/6, gff4[i,5] - gff4[1,4], h_distance+h_distance/6, col = color)
			text(mean(c(gff4[i,4] - gff4[1,4], gff4[i,5] - gff4[1,4])), h_distance-h_distance/3, gff4[i,9], cex=0.3)
			}
		}
	# connection between syntenic genes between Chra and Chrb
	for(i in 1:nrow(gff1)){
		loc1 <- Find_loc_for_gene_in_alle(gff1[i,9])
		for(j in 1:nrow(gff2)){
			loc2 <- Find_loc_for_gene_in_alle(gff2[j,9])
			if(rownames(loc1)[1] == rownames(loc2)[1]){
				segments(mean(c(gff1[i,4] - gff1[1,4], gff1[i,5] - gff1[1,4])), 4*h_distance-h_distance/6, mean(c(gff2[j,4] - gff2[1,4], gff2[j,5] - gff2[1,4])), 3*h_distance+h_distance/6, col = 'gray')
				}
			}
		}
	# connection between syntenic genes between Chrb and Chrc
	for(i in 1:nrow(gff2)){
		loc2 <- Find_loc_for_gene_in_alle(gff2[i,9])
		for(j in 1:nrow(gff3)){
			loc3 <- Find_loc_for_gene_in_alle(gff3[j,9])
			if(rownames(loc2)[1] == rownames(loc3)[1]){
				segments(mean(c(gff2[i,4] - gff2[1,4], gff2[i,5] - gff2[1,4])), 3*h_distance-h_distance/6, mean(c(gff3[j,4] - gff3[1,4], gff3[j,5] - gff3[1,4])), 2*h_distance+h_distance/6, col = 'gray')
				}
			}
		}
	# connection between syntenic genes between Chrc and Chrd
	if(number == 4) {
		for(i in 1:nrow(gff3)){
			loc3 <- Find_loc_for_gene_in_alle(gff3[i,9])
			for(j in 1:nrow(gff4)){
				loc4 <- Find_loc_for_gene_in_alle(gff4[j,9])
				if(rownames(loc3)[1] == rownames(loc4)[1]){
					segments(mean(c(gff3[i,4] - gff3[1,4], gff3[i,5] - gff3[1,4])), 2*h_distance-h_distance/6, mean(c(gff4[j,4] - gff4[1,4], gff4[j,5] - gff4[1,4])), h_distance+h_distance/6, col = 'gray')
					}
				}
			}
		}
	dev.off()
}

Alle_location(gene)















