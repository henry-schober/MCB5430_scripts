# ggplot iris data

gg <- ggplot(iris, aes(Species, Sepal.Width, color=Species, fill=Species)) 

#boxplot

gg + geom_boxplot()
gg + geom_boxplot(notch = T)

#jitter plot

gg + geom_jitter(width=0.1, alpha=0.5)

# violin plot

gg + geom_violin(width=0.5, alpha=0.5, lwd=1, draw_quantiles = c(0.25, .5, .75)) + geom_jitter(width=0.1)

#in class

dn <- ggplot(DNase, aes(conc, density, color=Run))
dn + geom_point() +scale_color_brewer(palette = 'Spectral')

#Dnase density point plots

dn +geom_point(color=DNase$Run) + facet_wrap(~Run)

# histogram
ghist <- ggplot(iris, aes(Sepal.Width, color=Species, fill=Species))

ghist + geom_histogram(binwidth=0.2, alpha=0.2) + facet_wrap(~Species, ncol=1)

ghist + geom_density(alpha=0.3) + facet_wrap(~Species, ncol=1)


#CDF plots
gdist <- ggplot(iris, aes(Sepal.Width, color=Species, fill=Species))

gdist + stat_ecdf(geom="step") + geom_hline(yintercept = 0.5, linetype="dashed") +
  geom_vline(xintercept = c(aggregate(iris$Sepal.Width ~ iris$Species, FUN=median)[,2]), linetype="dashed")

medians <- c(aggregate(iris$Sepal.Width ~ iris$Species, FUN=median)[,2])

gdist + stat_ecdf(geom="step") + geom_hline(yintercept = 0.5, linetype="dashed") +
  geom_vline(xintercept = medians, linetype="dashed")




