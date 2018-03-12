library("corrplot")

################ correlation heatmap  
# Input:  
#   dataframe of methylation values, 
#     sample.start = methylation value of first sample
#     sample.end = methylation value of last sample
#     CHR, MAPINFO
#   titleinfo = text for title
# Output: 
#   correlation heatmap


corrHeatMap <- function(methylation.data, sample.start, sample.end, titleinfo){

  one <- methylation.data
  
  rownames (one) <- one$ProbeID
  one.ordered <- one[order(one$CHR, one$MAPINFO) ,]
  
  one.t <-t(one.ordered[, sample.start:sample.end])
  
  correlation<-cor(one.t)  

  require (corrplot)
  corrplot(correlation, method="number", title = titleinfo, mar=c(5,5,5,5)) # from corrplt
}


############### Usage 
# combined.dmrcate.cpgs <- readRDS ("combined.dmrcate.cpgs.RDS") # dataset produced by file "combineDMRcpgs.R
# 
# ##### loop for in all pairwise correlations DMRs
# dmrs.cpgs <- combined.dmrcate.cpgs [, c(1, 18:32, 33:34)]  #col 1 = dmr.order, col 18:32 = methylation values, col33:34 = CHR, MAPINFO
# temp <- dmrs.cpgs [which(dmrs.cpgs$dmr.order==50) ,]
# row.names(temp)<-temp$X.y
# 
# corrHeatMap (methylation.data = temp, sample.start = 3, sample.end = 16)
# 
