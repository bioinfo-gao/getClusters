data.dir <- "C:/Users/lxw391/Dropbox (BBSR)/Saurav (Chu) second/DMR-comparions/LWtest/DMRcate_9-25-17/After_LWcorr/"
aclust <- read.csv (paste0(data.dir, "A-clust-results.csv"))

setwd("C:/Users/lxw391/Dropbox (BBSR)/Saurav (Chu) second/DMR-comparions/LWtest/Test_10-3-2017")

source ("C:/Users/lxw391/Dropbox (BBSR)/LW/MethylationFunctions/clusterLinePlot.R")
source ("C:/Users/lxw391/Dropbox (BBSR)/LW/MethylationFunctions/corrHeatMap.R")


pdf ("aclust-line-plots-corrHeatmaps.pdf")

for (i in 1:20)
{
  
  cluster <- aclust[which(aclust$Clusternumber==i) ,] 
  names(cluster)[names(cluster) == 'cpg'] <- 'ProbeID'
  
  clusterLinePlot (methyl.wide = cluster, methyl.value="beta.value", sample.start=9, sample.end=22, 
                  center=TRUE, center.by="sample.no", group.info=NULL, methyl.value.analyze="beta.value",
                  titleinfo = paste0("cluster = ", i))
  
  corrHeatMap (methylation.data = cluster, sample.start=9, sample.end=22, titleinfo=paste0("cluster = ", i) )
  
}

dev.off ()