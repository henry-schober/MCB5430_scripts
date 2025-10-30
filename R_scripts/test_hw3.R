
closest_summits <- read.table("ATF1_summits_closestTSS.txt", stringsAsFactors = F)
closest_motifs <- read.table("ATF1_motif_closestTSS.txt", stringsAsFactors = F)


boxplot(closest_motifs[,13], closest_summits[,12], outline=F, border='black', col=c('green', 'blue'), lwd=2, names=c('Motifs', 'Summits'), ylab= 'Distance', main='Peak and motif distances to the TSS')

# CDF plot

motif_cdf <- ecdf(closest_motifs[,13])
summits_cdf <- ecdf(closest_summits[,12])

plot(motif_cdf, xlim=c(-100000, 100000), col='green', main= "CDF plot of peaks and motif distances to TSS", lwd =3, ylab="Fraction", xlab="Distances")
lines(summits_cdf, col='blue', lwd=2)
legend("topleft", legend = c("Motif", "Summits"), col = c("green", "blue"))


#pie 

summits_overlaps <- read.table("ATF1_summits_overlap.txt", stringsAsFactors = F)

overlaps_sum_df <- data.frame(region = c("promoter", "gene", "intergenic"), count = c(sum(summits_overlaps[,6]),sum(summits_overlaps[,7]), sum(summits_overlaps[,8])))

pct <- round((overlaps_sum_df[,2]/sum(overlaps_sum_df[,2]))*100)

lbls <- c("promoter", "gene", "intergenic")
lbls2 <- paste(lbls, " ", pct, "%", sep = "")

pie(overlaps_sum_df$count, labels=lbls2, main="Simple pie chart w/ percentages", col = cbPalette)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")