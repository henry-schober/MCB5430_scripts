# Working with mtcars

cyls <- levels(as.factor(mtcars$cyl))

counts <- summary(as.factor(mtcars$cyl))

df <- data.frame(cyls=cyls, count=counts)

cp <- ggplot( df, aes(x="", count, color=cyls, fill=cyls))

cp + geom_bar(stat= "identity")





## Pie charts

slices <- c(10, 12, 4, 8, 16)

lbls <- c("promoter", "enhancers", "exons", "introns", "intergenic")

pct <- round((slices/sum(slices))*100)
lbls2 <- paste(lbls, " ", pct, "%", sep = "")

pie(slices, labels=lbls, main="Simple pie chart")

pie(slices, labels=lbls2, main="Simple pie chart w/ percentages", col=cbPalette)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#plotrix

pie3D(slices, labels = lbls2, explode = 0.1, main="3D pie chart")

#scatter plots

gp <- ggplot(iris, aes(Sepal.Width, Sepal.Length, color=Species, fill=Species)) +geom_point()

gp +facet_wrap(~Species, ncol=1) + geom_smooth(method = "lm", se=F, color="black", lwd=0.5) + stat_cor(method="spearman")


#mtcars -lm

fit <- lm(formula=mpg~cyl, data=mtcars)

corr <- cor(mtcars$cyl, mtcars$mpg, method="spearman")

plot(mpg ~ cyl, data =mtcars)
abline(fit[1], fit[2], lwd=2, col='red')
legend("topright", paste("rho=", corr, sep=""), inset=0.1)

#ATF1 scatterplots

ATF1_rep1 <- read.table("ATF1_rep1_overlap_readsInPeaks_sorted.txt", header = F, stringsAsFactors = F)

ATF1_rep2 <- read.table("ATF1_rep2_overlap_readsInPeaks_sorted.txt", header = F, stringsAsFactors = F)

qplot(log2(ATF1_rep1[, 6]), log2(ATF1_rep2[, 6]))

data = data.frame(ATF1_rep1[, 6], ATF1_rep2[, 6])

ggplot(data, aes(x = log2(data[, 1]), log2(data[, 2]))) + geom_point() + geom_bin2d()

heatscatter(log2(ATF1_rep1[, 6]), log2(ATF1_rep2[,6]))

heatscatter(log2(ATF1_rep1[, 6]), log2(ATF1_rep2[,6]), pch = 20, cex = 0.5, cor = T, greyscale = T, add.contour = T, main = "ATF1 peak intensities")
