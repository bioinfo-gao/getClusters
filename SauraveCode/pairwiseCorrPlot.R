
######### Function pairwiseCorrPlot -----------------------------------------
# Input:  
#   dataframe of dmrs merged with cpgs, with variables sample.start, sample.end, CHR, MAPINFO
# Output: 
#   correlation heatmap, a vector of pairwise correlations

pairwiseCorrPlot <- function (cpg.set, sample.start, sample.end, plot) {
  
  # cpg.set <- t
  # sample.start <- 3
  # sample.end <- 16
  
  one <- cpg.set
  
  one.ordered <- one[order(one$CHR, one$MAPINFO) ,]
  
  one.t <-t(one.ordered[, sample.start:sample.end])
  
  correlation<-cor(one.t)
  
  if (plot == TRUE) {
    require(corrplot)
    corrplot(correlation, method="number")
  }
  return (as.data.frame (correlation[upper.tri(correlation)]))
}

############# Example Usage ----------------------------------------------------

combined.dmrcate.cpgs <- readRDS ("combined.dmrcate.cpgs.RDS") # dataset produced by file "combineDMRcpgs.R

##### loop for in all pairwise correlations DMRs
dmrs.cpgs <- combined.dmrcate.cpgs [, c(1, 18:32, 33:34)]  #col 1 = dmr.order, col 18:32 = methylation values, col33:34 = CHR, MAPINFO
temp <- dmrs.cpgs
row.names(temp)<-temp$X.y


corr.all.pairwise <-as.data.frame(matrix(ncol=2,nrow=0)) #initialize

for (i in 50:50){
  
  #i<-50
  
  t <- temp[which(temp$dmr.order==i) ,]

  correlation.vector <- pairwiseCorrPlot (cpg.set = t, sample.start = 3, sample.end = 16, plot=TRUE)
  
  corr.all.pairwise <- rbind( corr.all.pairwise,as.data.frame(cbind(i, correlation.vector)))
  
}

colnames (corr.all.pairwise) <- c("dmr.order", "pairwise.corr")


######### side by side comparison of multiple methods - for later use
#ggplot(cor,aes(x=model,y=pairwise_correlation)) + geom_boxplot(aes(color=model), width=0.5)+ theme_bw() + 
scale_color_manual(values=c("blue", "red")) + theme(legend.position = c(0.3, 0.2)) +
  theme(legend.text=element_text(size=20, face="bold"), legend.key.size = unit(3, "line"),
        axis.title.y=element_text(size=20, face="bold"), axis.text.y = element_text(size=20, face="bold") )