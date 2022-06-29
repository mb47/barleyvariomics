args<-commandArgs(TRUE)

write("read the Fst file", stdout())
data <- read.table(args[1], header=T)

write("order the chromosomes", stdout())
chromosomeOrder <- c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H", "chrUn")

write("apply this as a factor", stdout())
data$chr <- factor(data$CHROM,levels=chromosomeOrder)

write("subset the data to exclude chrUn", stdout())
dataSubset <- data[grep ("chrUn", data$chr, invert=TRUE), ]


write("plot Fst", stdout())
library(ggplot2)
library(zoo)

fstPlot<- ggplot(dataSubset, aes(POS/1000000, WEIR_AND_COCKERHAM_FST)) +
geom_point(size=0.5) +
#geom_line(colour="red",aes(y=rollmean(WEIR_AND_COCKERHAM_FST, 100, na.pad=TRUE))) +
facet_wrap(~ chr,ncol=1) +
ylim(0,1) +
ggtitle("") +
xlab("Position (Mbp)") +
ylab(expression("F"["ST"])) +
theme_bw(base_size = 60) +
theme(strip.background = element_rect(fill="white"))


write("export to a PNG file", stdout())
png(args[2],1500,3000)
print(fstPlot)
dev.off()

