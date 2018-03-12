setwd("C:/Users/lxw391/Box Sync/METHODS-SHARED/MethylationFunctions")

######### 1. Function to merge cpgs with DMRs -----------------------------------------
# Input: 
#   dmrs - a dataframe with variables start, end, seqnames, sig (1,0)
#   cpgs - a dataframe with variables MAPINFO, seqnames
# Output: 
#   a dataframe with dmr info, cpg info, ncpgs

combine.dmrs.cpgs <- function(dmrs, cpgs) { 
  require(GenomicRanges)

  # 1. make ranges of dmr info  
  dmrs <- dmrs[which(dmrs$sig==1) ,]
  sig.ranges <- IRanges(dmrs$start, dmrs$end)
  dmr.ranges <- GRanges (seqname=dmrs$seqnames, ranges=sig.ranges)
  
  # 2. make ranges of cpgs info
  temp.cpg.ranges <- IRanges(cpgs$MAPINFO, cpgs$MAPINFO)
  cpg.ranges <-GRanges (seqname=cpgs$seqname, ranges=temp.cpg.ranges) 
  
  #3. find overlaps
  dmr.cpgs.overlaps <- as.data.frame(findOverlaps(query=dmr.ranges, subject=cpg.ranges, maxgap=0L, minoverlap=1L, type="any"))
  
  dmr.cpgs.overlaps$dmr.order <- as.numeric(dmr.cpgs.overlaps$queryHits)
  dmr.cpgs.overlaps$cpg.order <- as.numeric(dmr.cpgs.overlaps$subjectHits)
  
  # 4. merge with dmr info
  
  dmrs$dmr.row <- as.numeric(seq.int(nrow(dmrs)))
  
  dmrs.info <- merge (x=dmr.cpgs.overlaps, y=dmrs, by.x="dmr.order", by.y="dmr.row")
  
  # 5. merge with cpgs info
  
  cpgs$row <- as.numeric(row.names(cpgs))
  
  dmr.cpg.overlap.info <- merge (dmrs.info, cpgs, by.x="cpg.order", by.y="row")
  
  dmr.cpg.overlap.info <- dmr.cpg.overlap.info[order(dmr.cpg.overlap.info$dmr.order) ,]
  
  # 6. add number of cpgs
  ncpgs <- as.data.frame(table (dmr.cpg.overlap.info$dmr.order))
  
  dmr.cpg.overlap.ncpgs <-merge (dmr.cpg.overlap.info, ncpgs, by.x="dmr.order", by.y="Var1")
  
  
  return (dmr.cpg.overlap.ncpgs)
} 


############## Example Usage ------------------------------------------------------------------------

# Prepare files 
data.dir <- "C:/Users/lxw391/Dropbox (BBSR)/Saurav (Chu) second/DMR-comparions/LWtest/DMRcate_9-25-17/"

require(DMRcate)
dmrs.dmrcate <- read.csv (paste0(data.dir, "temp.ranges.csv"))
betas <- read.csv(paste0(data.dir, "GSE41169meth_betaValue.csv"))

dmrs.dmrcate$sig <- ifelse (dmrs.dmrcate$Stouffer<0.05, 1, 0)

load("fullannotInd.rda")   # data downloaded from https://rforge.net/IMA/ , under "4 Annotation file"
annot<-as.data.frame(fullannot[, c("ILMNID", "CHR", "MAPINFO")])

betas.location <- merge (betas, annot, by.x="X", by.y="ILMNID")

betas.location$MAPINFO <- as.numeric(as.character(betas.location$MAPINFO))

betas.location$seqname <- paste0("chr", betas.location$CHR)

# call function 

combined.dmrcate.cpgs <- combine.dmrs.cpgs (dmrs = dmrs.dmrcate, cpgs = betas.location)

saveRDS(combined.dmrcate.cpgs, "combined.dmrcate.cpgs.RDS")
