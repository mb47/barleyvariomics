args<-commandArgs(TRUE)

write("read the input file", stdout())
data <- read.table(args[1], header=T, sep="\t")

write("order the chromosomes", stdout())
chromosomeOrder <- c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H", "chrUn")

write("apply this as a factor", stdout())
data$chr <- factor(data$CHROM,levels=chromosomeOrder)

write("subset the data to exclude chrUn", stdout())
dataSubset <- data[grep ("chrUn", data$chr, invert=TRUE), ]

write("plot data", stdout())
library(ggplot2)
library(zoo)

rollAvgWindowSize <- as.integer(args[3])
lineWidth <- as.integer(args[4])
title <- args[5]
yLabel <- args[6]

piPlot<- ggplot(dataSubset, aes(POS/1000000, PI)) + 
geom_line(size = lineWidth, colour="#b2182b",aes(x=dataSubset$POS/1000000, y=rollmean(dataSubset$spontaneums, rollAvgWindowSize, na.pad=TRUE))) + 
geom_line(size = lineWidth, colour="#ef8a62",aes(x=dataSubset$POS/1000000, y=rollmean(dataSubset$landraces, rollAvgWindowSize, na.pad=TRUE))) + 
geom_line(size = lineWidth, colour="#2166ac",aes(x=dataSubset$POS/1000000, y=rollmean(dataSubset$cultivars, rollAvgWindowSize, na.pad=TRUE))) + 
facet_wrap(~ dataSubset$chr,ncol=1) + 
ggtitle(title) + 
xlab("Position (Mbp)") + 
ylab(yLabel) + 
theme_bw(base_size = 40) +
theme(strip.background = element_rect(fill="white"))

write("export to a PNG file", stdout())
png(args[2],1000,3000)
print(piPlot)
dev.off()

