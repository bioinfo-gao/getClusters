
################### Plot cpg methylation values in a line plot -------------------
# Input: 
#   wide dataset for methylation values, with variables ProbeID, CHR, MAPINFO 
#   methyl.value = type of methylation values provided, can be "beta.value" or "M.value"
#   methyl.value.analyze = type of methylation value for plot
#   sample.start
#   sample.end
#   n.sample = number of samples 
#   center = whether plot needs to be centered, so all cpgs have mean methylation values of 0, 
#             useful to put all lines on the same plot
#   center.by = unit for centering methylation values, can be "sample.no" or a group, e.g. "disease.stage"
#   group.info = dataframe for group, with variables "sample.id", and "group"
#   titleinfo = text about the title, e.g. cluster number  
# Output: line plot of cpgs by samples or groups

clusterLinePlot <- function (methyl.wide, methyl.value, sample.start, sample.end, center, center.by, group.info, methyl.value.analyze, titleinfo)
{
  # methyl.wide <- cluster
  # methyl.value <- "beta.value"
  # sample.start <- 9
  # sample.end <- 22
  # 
  # #sample.start <- 1 
  # #sample.end <- 14  
  # center = TRUE
  # center.by = "sample.no"
  # 
  # group.info <- groups2
  # methyl.value.analyze <- "beta.value"
  
  
  
one <- methyl.wide [, c(grep("ProbeID", colnames(methyl.wide)),
                        grep("CHR", colnames(methyl.wide)), 
                        grep("MAPINFO", colnames(methyl.wide)), 
                        sample.start:sample.end)]
  
# convert to tall format
require (reshape2)
tall<- as.data.frame( melt(one, id=c("ProbeID", "CHR", "MAPINFO")) )

if (methyl.value == "beta.value") {
  colnames (tall) <- c("ProbeID", "CHR", "MAPINFO", "sample", "beta.value")
  tall$M.value <- log2 ( tall$beta.value/ (1- tall$beta.value) )
}

if (methyl.value == "M.value") {
  colnames (tall) <- c("ProbeID", "CHR", "MAPINFO", "sample", "M.value")
  tall$beta.value <- (2^tall$M.value)/ (1 + (2^ tall$M.value) )
}  
  
tall.ordered <- tall[order(tall$CHR, tall$MAPINFO) ,]

tall.ordered$value <- tall.ordered[ ,methyl.value.analyze]

tall.ordered$sample.no <- c(rep(1: (sample.end - sample.start + 1), nrow(one)))

tall.ordered$ProbeID <- as.factor (tall.ordered$ProbeID)

require(ggplot2)

if (center == FALSE){ 
  plot <- ggplot(tall.ordered, aes(x=sample.no, y=value, color=as.factor(ProbeID) ))+ 
    geom_point(size=1, shape=21) + geom_line() + theme_bw() + 
    scale_x_continuous(name="sample number", breaks = unique(tall.ordered$sample.no)) +
    ggtitle(paste0("plotted are ", methyl.value.analyze, ", ", titleinfo))
  print (plot)
}

if (center == TRUE & center.by == "sample.no") {
  means.cpg <-aggregate(value ~ ProbeID, data=tall.ordered, mean)
  names(means.cpg)[names(means.cpg) == 'value'] <- 'av.cpg.value'
  
  means.cpg.center <- merge (means.cpg, tall.ordered, by="ProbeID")
  means.cpg.center$centered.value <- means.cpg.center$value - means.cpg.center$av.cpg.value
  
  plot <- ggplot(means.cpg.center, aes(x=sample.no, y=centered.value, color=ProbeID) ) + geom_point(size=3, shape=21) + 
       geom_line() + theme_bw() + scale_x_continuous(name="sample number", breaks = unique(means.cpg.center$sample.no) ) + 
       ggtitle(paste0("plotted are ", methyl.value.analyze, ", ", titleinfo))
  
  print (plot)
}  

if (!is.null (group.info) )
{
  # merge with group information
  
  tall.group <- merge (tall.ordered, group.info, by.x ="sample") 
  
  tall.group$grp <- tall.group[, center.by]
  
  means.cpg <-aggregate(value ~ grp + ProbeID, data=tall.group, mean)
  
  av <-aggregate (value ~ ProbeID, data=means.cpg, mean) 
  
  names(av)[names(av) == 'value'] <- 'mean.cpg.allsamples'
  
  means.cpg.center <- merge (means.cpg, av, by="ProbeID")
  means.cpg.center$value.centered <-means.cpg.center$value - means.cpg.center$mean.cpg.allsamples
  
  plot <- ggplot(means.cpg.center, aes(x=grp, y=value.centered, color=ProbeID) ) + 
    geom_point(size=3, shape=21) + geom_line() +
    scale_x_continuous(name=center.by, breaks = unique(means.cpg.center$grp) ) + 
    theme_bw() + ggtitle(paste0("plotted are ", methyl.value.analyze, ", ", titleinfo))
  
  print (plot)
}


}

######################### Usage Examples


# 
# setwd("C:/Users/lxw391/Dropbox (BBSR)/Zhen_Gao/MethylationFunctions")
# 
# # an exmaple of 14 samples
# combined.dmrcate.cpgs <- readRDS ("combined.dmrcate.cpgs.RDS") # dataset produced by file "combineDMRcpgs.R
# combined.dmrcate.cpgs$ProbeID <- combined.dmrcate.cpgs$X.y
# one <- combined.dmrcate.cpgs[which(combined.dmrcate.cpgs$dmr.order==50), c(19:32, 33:34, 37)]
# clusterLinePlot (methyl.wide = one, methyl.value="beta.value", sample.start=1, sample.end=14, center=FALSE, center.by="sample.no",
#                  group.info=NULL, methyl.value.analyze="M.value", titleinfo = "test")
# 
# # an exmaple of 110 brain samples
# oneisland <- read.csv ("oneisland.csv")
# groups <- read.csv ("sampleInfo.csv")
# groups$sample <- groups$Sample
# groups2 <- subset (groups, select=-c(Sample,X))
# one<- subset(oneisland, select=-X)
# clusterLinePlot (methyl.wide = one, methyl.value="beta.value", sample.start=2, sample.end=111, center=FALSE, center.by="Stage",
#                  group.info=groups2, methyl.value.analyze="beta.value", titleinfo="test")


# test
# test
# test - don't want to merge to master


# test

# test - don't want to merge to master

