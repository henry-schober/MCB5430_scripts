## Functions practice

dat <- c(23, 5, 67, 94, 102)

meanFunc <- function(data)
{
  # enter commands
  m = mean(data)
  
  #print the result to screen to variable return()
  
  return(m)
}


meanFunc(dat)

m <- meanFunc(dat)

printFunc <- function(data, add=T, num=5)
{
  for (ii in data)
  {
    if (add)
    {
      print(ii+num)
    }
    else
    {
      print(ii)
    }
  }
}

printFunc(dat, num=7)


## Heatmap correlations

cormat <- cor(mtcars)

library(RColorBrewer)
heatmap(cormat, Colv = NA, Rowv = NA)

heatmap(cormat, keep.dendro = F, col=colorRampPalette(brewer.pal(9, "Blues"))(9))

legend("left", legend = quantile(cormat, probs = seq(0.5, 0.95, 0.1)), cex=.8, fill = colorRampPalette(brewer.pal(9, "Blues"))(9))


library(corrplot)

rquery.cormat(cormat, type="upper")


ggheatmap <- function(data, low="blue", mid="red", high="green")
{
  library(corrplot)
  library(ggplot2)
  library(reshape2)
  
  cormat <- cor(data)
  
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  # Reorder the correlation matrix
  cormat <- reorder_cormat(cormat)
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = low, high = high, mid = mid, 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  # Print the heatmap
  print(ggheatmap)
}

ggheatmap(mtcars, low="blue", mid="yellow", high="red")




