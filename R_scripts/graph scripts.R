




#plotting dist. of ATF1 peaks on 11/4/24

closest <- read.table("ATF1_summits_closestTSS.txt", stringsAsFactors = F)

#boxplots

boxplot(closest[,12], outline =F, border = 'blue', col = 'white', lwd=3, xlab="ATF1_peaks", ylab="distance to TSS (bp)", main="Distance of peaks to gene TSSs")

#histogram

#bars
closest_hist <- hist(closest[,12], breaks=500, xlim=c(-25000, 25000), main="Distance of peaks to gene TSSs", xlab="Distance to TSS (bp)")

#points
closest_hist <- hist(closest[,12], breaks=500, xlim=c(-25000, 25000), main="Distance of peaks to gene TSSs", xlab="Distance to TSS (bp)")


plot(closest_hist$mids, closest_hist$counts, xlim=c(-25000, 25000))

#lines
plot(closest_hist$mids, closest_hist$counts, xlim=c(-25000, 25000), type="l", lwd=3, col='blue')

# Bars and Lines
closest_hist <- hist(closest[,12], breaks=500, xlim=c(-25000, 25000), main="Distance of peaks to gene TSSs", xlab="Distance to TSS (bp)")
lines(closest_hist$mids, closest_hist$counts, xlim=c(-25000, 25000), type="l", lwd=3, col='blue')
points(closest_hist$mids, closest_hist$counts, xlim=c(-25000, 25000), pch=19, cex=1.5, col='blue')

# CDF plots
closest_cdf <- ecdf(closest[,12])

plot(closest_cdf, xlim=c(-5000, 5000))
abline(h=0.5, lty=2)
abline(v=(median(closest[,12])), lty=2, col='red')
legend("topleft", c("ATF1"), pch=19, inset=0.1, col='red')
dev.copy(jpeg, filename="ATF1_cdf.jpeg")
dev.off()




