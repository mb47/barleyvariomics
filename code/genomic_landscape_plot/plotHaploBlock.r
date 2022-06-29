args<-commandArgs(TRUE)

## USAGE:
## Rscript plotHaploBlock.r <formatted.det.file> <output prefix> <block length threshold in Kb>
##
#read the BED file with the SNP coordinates
blocks <- read.table(args[1],sep="\t",header=T)
prefix <- args[2]

#order the chromosomes
goodChrOrder <- c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H", "chrUn")
#apply this as a factor
blocks$CHR <- factor(blocks$CHR,levels=goodChrOrder)

#subset data to remove chrUn
blocks <- blocks[grep ("chrUn", blocks$CHR, invert=TRUE), ]

library(ggplot2)

#======================================================
# Segment plot of haplotype blocks
#======================================================
blockPlot <- ggplot(blocks, colour=CHR) + 
geom_segment(aes(x=blocks$BP1, y=CHR, xend=blocks$BP2, yend=CHR),size=15) + 
ggtitle("") + 
xlab("Position (bp)") + 
ylab("Chromosomes") + 
theme_void(base_size = 20)

#export to a PNG file
outName <- paste(prefix,".segmentPlot.png",sep="")
png(outName,2000,1000)
print(blockPlot)
dev.off()
#======================================================

#======================================================
# Histogram of haplotype block length distribution
#======================================================
blockPlot <- ggplot(blocks, aes(KB)) + 
geom_histogram(binwidth=100) + 
facet_wrap(~ blocks$CHR,ncol=1) + 
ggtitle("") + 
xlab("Block size (100 kb)") + 
ylab("Chromosomes") + 
theme_bw(base_size = 20)

#export to a PNG file
outName <- paste(prefix,".lengthHistogram.png",sep="")
png(outName,1000,3000)
print(blockPlot)
dev.off()
#======================================================

#
# Subsetting data using a given block length threshold
# 
threshold <- as.numeric(args[3])
blocksSubset <- blocks[ which(blocks$KB > threshold),]

#======================================================
# Segment plot of haplotype blocks
#======================================================
blockPlot <- ggplot(blocksSubset, colour=CHR) + 
geom_segment(aes(x=blocksSubset$BP1, y=CHR, xend=blocksSubset$BP2, yend=CHR),size=15) + 
ggtitle("") + 
xlab("Position (bp)") + 
ylab("Chromosomes") + 
theme_void(base_size = 20)

#export to a PNG file
outName <- paste(prefix,".subsetInKB",threshold,".segmentPlot.png",sep="")
png(outName,2000,1000)
print(blockPlot)
dev.off()
#======================================================

#======================================================
# Histogram of haplotype block length distribution
#======================================================
blockPlot <- ggplot(blocksSubset, aes(KB)) + 
geom_histogram(binwidth=100) + 
facet_wrap(~ blocksSubset$CHR,ncol=1) + 
ggtitle("") + 
xlab("Block size (100 kb)") + 
ylab("Chromosomes") + 
theme_bw(base_size = 20)

#export to a PNG file
outName <- paste(prefix,".subsetInKB",threshold,".lengthHistogram.png",sep="")
png(outName,1000,3000)
print(blockPlot)
dev.off()
#======================================================
