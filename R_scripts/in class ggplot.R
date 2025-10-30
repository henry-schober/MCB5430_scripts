ATF1 <- read.table("ATF1_regions_overlap.txt", header =T, stringsAsFactors = F)

overlaps <- data.frame(region = c("promoter", "gene","intergenic"), count = c(sum(ATF1$promoter), sum(ATF1$gene), sum(ATF1$intergenic)))

#bar plots
overlap_graph <- ggplot(overlaps, aes(region, counts, color=region, fill=region))
overlap_graph + geom_bar(stat = "identity")


#IN class exercise 2

corr <- cor(mtcars$mpg, mtcars$hp, method = "spearman")
plot(hp ~ mpg, data = mtcars)
legend("topright", paste("rho=", corr, sep=""), inset=0.1)

ggplot(mtcars, aes(mpg, hp)) + geom_point() + stat_cor(method = "spearmans")
