
args<-commandArgs(TRUE)

## I/O
lineWidth <- 1
baseSize <- 45

##==================================================================
## Load libraries
##==================================================================
write("--------Loading libraries",stdout())
library(ggplot2)
library(zoo)
library(gtable)
library(gridExtra)
# library(ggpubr)
library(gridGraphics)
library(cowplot)

##==================================================================
## Chromosome ordering
##==================================================================
chromosomeOrder <- c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H", "chrUn")

##==================================================================
## Coordinates
##==================================================================
write("--------Read and plot coordinates", stdout())
coor <- read.table(args[1], header=T, sep="\t")
coor$chr <- factor(coor$CHROM,levels=chromosomeOrder)

xlab <- c(0,200,400,600)

coorPlot <- ggplot(coor) +
geom_segment(aes(x=coor$START/1000, y=1, xend=coor$END/1000, yend=1),size=12,color=coor$COLOR) +
scale_x_continuous(position = "top",labels=xlab,breaks=c(0,200000,400000,600000)) +
facet_wrap(~ chr,nrow=1) +
theme_bw(base_size = baseSize) +
theme(strip.background=element_rect(fill="white"), panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
theme(strip.text.x=element_blank())+
theme(legend.position = "none")

##==================================================================
## Pi plot
##==================================================================
write("--------Read and plot pi file", stdout())
data <- read.table(args[2], header=T, sep="\t")
data$chr <- factor(data$CHROM,levels=chromosomeOrder)
dataPiSubset <- data[grep ("chrUn", data$chr, invert=TRUE), ]

rollAvgWindowSize <- as.numeric(100)

piPlot<- ggplot(dataPiSubset, aes(POS, PI)) + 
geom_line(size = lineWidth, colour="#b2182b",aes(x=dataPiSubset$POS, y=rollmean(dataPiSubset$spontaneums, rollAvgWindowSize, na.pad=TRUE))) + 
geom_line(size = lineWidth, colour="#ef8a62",aes(x=dataPiSubset$POS, y=rollmean(dataPiSubset$landraces, rollAvgWindowSize, na.pad=TRUE))) + 
geom_line(size = lineWidth, colour="#2166ac",aes(x=dataPiSubset$POS, y=rollmean(dataPiSubset$cultivars, rollAvgWindowSize, na.pad=TRUE))) +
facet_wrap(~ dataPiSubset$chr,nrow=1) + 
scale_y_continuous(labels = function(x) format(x, scientific = FALSE), name="a")+
theme_classic(base_size = baseSize) +
theme(strip.background=element_rect(fill="white"), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank())+
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(strip.text.x=element_blank())

##==================================================================
## Fst plot
##==================================================================

## Cult vs Spot
write("--------Read and plot Fst file", stdout())
data <- read.table(args[3], header=T, sep="\t")
data$chr <- factor(data$CHROM,levels=chromosomeOrder)
dataFst1Subset <- data[grep ("chrUn", data$chr, invert=TRUE), ]

fstPlot1<- ggplot(dataFst1Subset, aes(POS/1000, WEIR_AND_COCKERHAM_FST)) +
geom_point(size=1,color=ifelse(dataFst1Subset$WEIR_AND_COCKERHAM_FST>=0.8,"red","black")) +
facet_wrap(~ chr,nrow=1) +
scale_y_continuous(breaks=c(0,0.5,1),name="b",limits=c(0,1))+
theme_classic(base_size = baseSize) +
theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(strip.text.x=element_blank())

## Land vs Spot
write("--------Read and plot Fst file", stdout())
data <- read.table(args[4], header=T, sep="\t")
data$chr <- factor(data$CHROM,levels=chromosomeOrder)
dataFst2Subset <- data[grep ("chrUn", data$chr, invert=TRUE), ]

fstPlot2<- ggplot(dataFst2Subset, aes(POS/1000, WEIR_AND_COCKERHAM_FST)) +
geom_point(size=1,color=ifelse(dataFst2Subset$WEIR_AND_COCKERHAM_FST>=0.8,"red","black")) +
facet_wrap(~ chr,nrow=1) +
scale_y_continuous(breaks=c(0,0.5,1),name="c",limits=c(0,1))+
theme_classic(base_size = baseSize) +
theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(strip.text.x=element_blank())

## Cult vs Land
write("--------Read and plot Fst file", stdout())
data <- read.table(args[5], header=T, sep="\t")
data$chr <- factor(data$CHROM,levels=chromosomeOrder)
dataFst3Subset <- data[grep ("chrUn", data$chr, invert=TRUE), ]

fstPlot3<- ggplot(dataFst3Subset, aes(POS/1000, WEIR_AND_COCKERHAM_FST)) +
geom_point(size=1,color=ifelse(dataFst3Subset$WEIR_AND_COCKERHAM_FST>=0.8,"red","black")) +
facet_wrap(~ chr,nrow=1) +
scale_y_continuous(breaks=c(0,0.5,1),name="d",limits=c(0,1))+
theme_classic(base_size = baseSize) +
theme(strip.background=element_blank(), panel.border=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank()) +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_text(angle = 0,face="bold"))+
theme(strip.text.x=element_blank())



##==================================================================
## Plotting the panel
##==================================================================
write("--------Running ggplotGrob coordinates",stdout())
g0 <- ggplotGrob(coorPlot)
write("--------Running ggplotGrob pi plot",stdout())
g1 <- ggplotGrob(piPlot)
write("--------Running ggplotGrob fst plot 1",stdout())
g2.1 <- ggplotGrob(fstPlot1)
write("--------Running ggplotGrob fst plot 2",stdout())
g2.2 <- ggplotGrob(fstPlot2)
write("--------Running ggplotGrob fst plot 3",stdout())
g2.3 <- ggplotGrob(fstPlot3)


write("--------Creating panel",stdout())
## cowplot method. easier to specify row height
#plotPanel <- plot_grid(g0,g1,g2.1,g2.2,g2.3,g3,g4,g5.1,g5.2,g5.3,g6,align="v",nrow=11,rel_heights=c(0.5,0.8,1,1,1,0.65,0.8,0.5,0.5,0.5,0.5))
plotPanel <- plot_grid(g0,g1,g2.1,g2.2,g2.3,align="v",nrow=5,rel_heights=c(0.8,1,1,1,1))

write("--------Export to a PNG file", stdout())
png(args[8],3000,1200)
grid.draw(plotPanel)
dev.off()
write("--------Done",stdout())
